import argparse
import logging
import warnings
from pathlib import Path
from typing import Dict
from typing import List
from typing import Literal
from typing import Optional
from typing import Tuple
from typing import Union

import allel
import multiprocessing as mp
import numpy as np
from allel import GenotypeArray
from allel import GenotypeVector
from allel import SortedIndex
from multiprocessing import Queue
from multiprocessing import queues
from multiprocessing.managers import BaseProxy
from numpy.typing import NDArray

from pixy.calc import calc_tajima_d
from pixy.calc import calc_watterson_theta
from pixy.calc import serialize_tajima_d_variant_counts
from pixy.enums import FSTEstimator
from pixy.enums import PixyStat
from pixy.models import PixyTempResult
from pixy.models import TajimaDResult
from pixy.models import WattersonThetaResult
from pixy.stats.summary import compute_summary_dxy
from pixy.stats.summary import compute_summary_fst
from pixy.stats.summary import compute_summary_pi
from pixy.stats.summary import precompute_filtered_variant_array

# Note: the previous pandas-based aggregation pipeline (`aggregate_output()` plus the helpers
# `_group_aggregated_output`, `_calculate_aggregated_stat`, `_coerce_aggregated_output_types`,
# `_aggregate_tajima_d_variant_counts`) has been replaced by the streaming pure-Python
# implementation in `pixy.agg.aggregate_rows()` / `pixy.agg.write_stat_file()`. Those used
# `pandas.read_csv` + `groupby` + `to_csv` to round-trip the tmp file through a DataFrame,
# which scaled memory linearly with the window count; the new path keeps a dict accumulator
# keyed by (pop1, pop2?, bin_idx) and emits output line-by-line.


# function for breaking down large windows into chunks
def assign_subwindows_to_windows(
    window_pre_list: List[List[int]], chunk_size: int
) -> List[List[int]]:
    """
    Splits each window into smaller subwindows of size `chunk_size`.

    Args:
        window_pre_list: a list of windows, where each window is represented
            as a list with two positions (start and stop).
        chunk_size: the size of each subwindow. Each subwindow will have a width of
            `chunk_size`.

    Returns:
        window_lst: A list of non-overlapping subwindows
    """
    # build list of subwindows
    window_lst = []
    for i in range(len(window_pre_list)):
        original_start = window_pre_list[i][0]
        original_stop = window_pre_list[i][1]
        subwindow_pos_1_list = [*range(original_start, original_stop, chunk_size)]
        subwindow_pos_2_list = [(item + chunk_size - 1) for item in subwindow_pos_1_list]

        # end of last window is always the end of the original window
        subwindow_pos_2_list[-1] = original_stop

        sub_windows = [
            list(a) for a in zip(subwindow_pos_1_list, subwindow_pos_2_list, strict=True)
        ]
        window_lst.extend(sub_windows)

    return window_lst


# function for assinging windows to larger chunks
# goal here to reduce I/O bottlenecks by only doing a single read of
# the VCF for a group of small windows
def assign_windows_to_chunks(window_pre_list: List[List[int]], chunk_size: int) -> List[List[int]]:
    """
    Assigns windows to larger chunks.

    Windows that overlap multiple chunks (based on start and stop positions) will be corraled
    into the first chunk.

    Args:
        window_pre_list: a list of windows, where each window is represented
            as a list with two positions (start and stop).
        chunk_size: the size of each subwindow. Each subwindow will have a width of
            `chunk_size`.

    Returns:
        window_lst: contains the start position, end position, and the chunk index for each window
    """
    # the start and end positions of each window
    window_pos_1_list = [item[0] for item in window_pre_list]
    window_pos_2_list = [item[1] for item in window_pre_list]

    # assign starts and ends of windows to chunks using np.floor
    chunk_list_1 = [np.floor(x / chunk_size) for x in window_pos_1_list]
    chunk_list_2 = [np.floor(x / chunk_size) for x in window_pos_2_list]

    # bind the lists back together
    window_lst = [
        list(a)
        for a in zip(window_pos_1_list, window_pos_2_list, chunk_list_1, chunk_list_2, strict=True)
    ]

    # nudge windows that overlap two chunks into the first chunk
    for i in range(len(window_lst)):
        if not window_lst[i][2] == window_lst[i][3]:
            window_lst[i][3] = window_lst[i][2]

    # remove chunk placeholders
    for i in range(len(window_lst)):
        del window_lst[i][3]

    return window_lst


# function for assinging sites to larger chunks
def assign_sites_to_chunks(sites_pre_list: List[int], chunk_size: int) -> List[List[int]]:
    """
    Assigns each site in a list to a chunk based on its position and given chunk size.

    The chunk index for each site is determined by dividing the site position by
    `(chunk_size + 1)` and using `np.floor` to calculate the chunk index. This
    function returns a list of sites, where each site is paired with the chunk
    index it belongs to.

    Args:
        sites_pre_list: the sites (positions) of interest
        chunk_size: the size of each chunk

    Returns:
        sites_list: A list where each element is a pair of a site and its corresponding
            chunk index
    """
    # assign sites to chunks using np.floor
    chunk_list = [np.floor(x / (chunk_size + 1)) for x in sites_pre_list]

    # bind the lists back together
    sites_list = [list(a) for a in zip(sites_pre_list, chunk_list, strict=True)]

    return sites_list


# function for masking non-target sites in a genotype array
# used in conjuctions with the --sites_list command line option
# NB: `SortedIndex` is of type `int`, but is not generic (at least in 3.8)
def mask_non_target_sites(
    gt_array: Union[GenotypeArray, allel.HaplotypeArray],
    pos_array: SortedIndex,
    sites_list_chunk: List[int],
) -> Union[GenotypeArray, allel.HaplotypeArray]:
    """
    Masks non-target sites in a genotype or haplotype array.

    Masked sites are set to missing data (`-1`). Haploid data is represented as a
    ``HaplotypeArray`` (2-D, no ``ploidy``/``n_samples`` attributes), so the missing-row shape is
    chosen per array type.

    Args:
        gt_array: the `GenotypeArray` (diploid+) or `HaplotypeArray` (haploid) of interest
        pos_array: positions corresponding to the sites in `gt_array`
        sites_list_chunk: target site positions that should remain unmasked

    Returns:
        gt_array: a modified array of the same type as the input
    """
    # get the indexes of sites that are NOT the target sites
    # (these will be masked with missing rows to remove them from the calculations)
    masked_sites: List[int] = sorted(set(pos_array) - set(sites_list_chunk))

    gt_mask_indexes: NDArray[np.intp] = np.flatnonzero(pos_array.locate_keys(masked_sites))

    # Build the missing-row mask as a numpy array (rather than a nested Python list) — the
    # previous `[[-1] * ploidy] * n_samples` form aliased the inner list across every sample row,
    # which is a latent bug if downstream code ever mutates row-wise. np.full is also faster.
    missing_row: NDArray[np.int8]
    if isinstance(gt_array, allel.HaplotypeArray):
        missing_row = np.full(gt_array.n_haplotypes, -1, dtype=np.int8)
    else:
        missing_row = np.full((gt_array.n_samples, gt_array.ploidy), -1, dtype=np.int8)

    # apply the mask to all non-target sites
    for pos in gt_mask_indexes:
        gt_array[pos, :] = missing_row

    return gt_array


# function for reading in a genotype matrix from a VCF file
# also filters out all but biallelic SNPs and invariant sites
# returns a genotype matrix, and array of genomic coordinates and
# a logical indicating whether the array(s) are empty
def read_and_filter_genotypes(
    args: argparse.Namespace,
    chromosome: str,
    window_pos_1: int,
    window_pos_2: int,
    sites_list_chunk: Optional[List[int]],
    ploidy: int,
) -> Tuple[bool, Optional[GenotypeArray], Optional[SortedIndex]]:
    """
    Ingests genotypes from a VCF file, retains biallelic SNPs or invariant sites.

    Filters out non-SNPs, multi-allelic SNPs, and non-variant sites. Optionally masks out
    non-target sites based on a provided list (`sites_list_chunk`).

    Args:
        args (Namespace): Command-line arguments, should include the path to the VCF file as
            `args.vcf`
        chromosome (str): The chromosome for which to read the genotype data (e.g., 'chr1')
        window_pos_1 (int): The start position of the genomic window
        window_pos_2 (int): The end position of the genomic window
        sites_list_chunk (list): An optional list of positions in which to mask non-target sites
            If `None`, no masking is applied.
        ploidy (int): The ploidy of the given chromosome, used to correctly shape the genotype
            array read from the VCF. Typically obtained from
            ``args.ploidy_map[chromosome]``.

    Returns:
        Tuple[bool, Optional[GenotypeArray], Optional[SortedIndex]]:
            - A boolean flag (`True` if the VCF callset is empty for the region, `False` otherwise)
            - A GenotypeArray containing the filtered genotype data (or `None` if no valid genotypes
              remain)
            - A SortedIndex containing the positions corresponding to the valid genotypes (or `None`
              if no valid positions remain)

    """
    # a string representation of the target region of the current window
    window_region = chromosome + ":" + str(window_pos_1) + "-" + str(window_pos_2)

    include_multiallelic_snps: bool = args.include_multiallelic_snps

    # read in data from the source VCF for the current window
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, module="allel")
        callset = allel.read_vcf(
            args.vcf,
            region=window_region,
            fields=[
                "CHROM",
                "POS",
                "calldata/GT",
                "variants/is_snp",
                "variants/numalt",
            ],
            numbers={"GT": ploidy},
        )

    # Pre-declare the array vars so mypy widens to the union type rather than narrowing to
    # whichever subclass appears first below (HaplotypeArray vs GenotypeArray).
    gt_array: Optional[GenotypeArray]
    pos_array: Optional[SortedIndex]

    # keep track of whether the callset was empty (no sites for this range in the VCF)
    # used by compute_summary_stats to add info about completely missing sites
    if callset is None:
        callset_is_none = True
        gt_array = None
        pos_array = None

    else:
        # if the callset is NOT empty (None), continue with pipeline
        callset_is_none = False

        # convert to a genotype array object
        if ploidy == 1:
            gt_array = allel.HaplotypeArray(allel.HaplotypeDaskArray(callset["calldata/GT"]))
        else:
            gt_array = allel.GenotypeArray(allel.GenotypeDaskArray(callset["calldata/GT"]))

        # build an array of positions for the region
        pos_array = allel.SortedIndex(callset["variants/POS"])

        # create a mask for biallelic snps and invariant sites
        is_biallelic_snp = np.logical_and(
            callset["variants/is_snp"][:] == 1,
            callset["variants/numalt"][:] == 1,
        )

        is_invariant_site = callset["variants/numalt"][:] == 0

        # build the mask for multiallelic or biallelic snps + invariant sites
        # NB: np.logical_or takes a maximum of TWO arrays

        snp_invar_mask = np.logical_or(is_biallelic_snp, is_invariant_site)

        if include_multiallelic_snps:
            is_multiallelic_snp = np.logical_and(
                callset["variants/is_snp"][:] == 1,
                callset["variants/numalt"][:] > 1,
            )
            snp_invar_mask = np.logical_or(snp_invar_mask, is_multiallelic_snp)

        # remove rows that are NOT snps or invariant sites from the genotype array
        gt_ndarray = np.delete(gt_array, np.where(np.invert(snp_invar_mask)), axis=0)

        # use a haplotype array instead of a genotype array for haploid data
        if ploidy == 1:
            gt_array = allel.HaplotypeArray(gt_ndarray)
        else:
            gt_array = allel.GenotypeArray(gt_ndarray)

        # select rows that ARE snps or invariant sites in the position array
        pos_array = pos_array[snp_invar_mask]

        # TODO: cannot index value of type None
        # if a list of target sites was specified, mask out all non-target sites
        if sites_list_chunk is not None:
            assert pos_array is not None, "pos_array should not be None"
            gt_array = mask_non_target_sites(gt_array, pos_array, sites_list_chunk)

        # extra 'none' check to catch cases where every site was removed by the mask
        if len(gt_array) == 0:
            callset_is_none = True
            gt_array = None
            pos_array = None

    return callset_is_none, gt_array, pos_array


# main pixy function for computing summary stats over a list of windows (for one chunk)


def compute_summary_stats(  # noqa: C901
    args: argparse.Namespace,
    popnames: NDArray,
    popindices: Dict[str, NDArray],  # NB: this is an array[int]  # TODO: add type param in 3.9
    temp_file: Path,
    chromosome: str,
    chunk_pos_1: int,
    chunk_pos_2: int,
    window_list_chunk: List[List[int]],
    q: Union[Queue, Literal["NULL"]],
    sites_list_chunk: Optional[List[int]],
    aggregate: bool,
    window_size: int,
) -> None:
    """
    Calculates summary statistics and writes them out to a temp file.

    Args:
        args: user-specified command-line arguments
        popnames: unique population names derived from `--populations` file
        popindices: the indices of the population-specific samples within the provided VCF
        temp_file: a temporary file in which to hold in-flight `pixy` calculations
        chromosome: the chromosome of interest
        chunk_pos_1: the start position of the genomic window
        chunk_pos_2: the end position of the genomic window
        window_list_chunk: the list of window start:stop that correspond to this chunk
        q: either "NULL" in single-core mode or a `Queue` object in multicore mode
        sites_list_chunk: list of positions in which to mask non-target sites
        aggregate: True if `window_size` > `chunk_size` or the chromosome is longer than the cutoff
        window_size: window size over which to calculate stats (in base pairs)

    Returns:
        None
    """
    pixy_output: List[PixyTempResult] = []

    # look up the ploidy of the current chromosome from the per-contig ploidy map
    # (built once during argument validation).
    ploidy_map: Dict[str, int] = getattr(args, "ploidy_map", {}) or {}
    if chromosome not in ploidy_map:
        raise KeyError(
            f"No ploidy entry for chromosome {chromosome!r} in args.ploidy_map. "
            "This indicates that the per-contig ploidy map was not initialized correctly."
        )
    chrom_ploidy: int = ploidy_map[chromosome]

    # Weir-Cockerham FST requires diploid data; skip FST for non-diploid contigs when WC is
    # requested. Hudson FST works for any ploidy.
    skip_fst_for_chrom: bool = (
        "fst" in args.stats and str(args.fst_type).upper() == "WC" and chrom_ploidy != 2
    )

    # read in the genotype data for the chunk
    callset_is_none, gt_array, pos_array = read_and_filter_genotypes(
        args, chromosome, chunk_pos_1, chunk_pos_2, sites_list_chunk, chrom_ploidy
    )

    # if computing FST, pre-compute a filtered array of variants (only)
    if "fst" in args.stats and not skip_fst_for_chrom:
        # These should only be returned by `read_and_filter_genotypes` as None if `fst` is not a
        # requested stat. The asserts are to narrow the types.
        assert gt_array is not None, "genotype array is None"
        assert pos_array is not None, "position array is None"

        per_site_fst_results, gt_array_fst, pos_array_fst = precompute_filtered_variant_array(
            args=args,
            gt_array=gt_array,
            pos_array=pos_array,
            callset_is_none=callset_is_none,
            window_size=window_size,
            popindices=popindices,
            chromosome=chromosome,
        )

        if per_site_fst_results is not None:
            pixy_output.extend(per_site_fst_results)

    # loop over the windows within the chunk and compute summary stats
    for window_index in range(0, len(window_list_chunk)):
        window_pos_1 = window_list_chunk[window_index][0]
        window_pos_2 = window_list_chunk[window_index][1]

        window_is_empty: bool
        gt_region: Union[GenotypeVector, None]
        loc_region: Union[slice, None]

        if pos_array is None:
            window_is_empty = True
            gt_region = None
        elif len(pos_array[(pos_array >= window_pos_1) & (pos_array <= window_pos_2)]) == 0:
            window_is_empty = True
            gt_region = None
        else:
            window_is_empty = False
            # pull out the genotypes for the window
            loc_region = pos_array.locate_range(window_pos_1, window_pos_2)

            # the assertion here is required to narrow the type on `gt_array` from `Optional`
            assert gt_array is not None, "genotype array is None"
            gt_region = gt_array[loc_region]

            # An empty slice after subsetting means the window has no usable sites.
            if len(gt_region) == 0:
                window_is_empty = True

        # PI:
        # AVERAGE NUCLEOTIDE DIFFERENCES WITHIN POPULATIONS

        if (args.populations is not None) and ("pi" in args.stats):
            pi_results: List[PixyTempResult] = compute_summary_pi(
                popnames=popnames,
                window_is_empty=window_is_empty,
                gt_region=gt_region,
                popindices=popindices,
                chromosome=chromosome,
                window_pos_1=window_pos_1,
                window_pos_2=window_pos_2,
            )

            pixy_output.extend(pi_results)

        # DXY:
        # AVERAGE NUCLEOTIDE DIFFERENCES BETWEEN POPULATIONS

        if (args.populations is not None) and ("dxy" in args.stats):
            # create a list of all pairwise comparisons between populations in the popfile
            dxy_results: List[PixyTempResult] = compute_summary_dxy(
                popnames=popnames,
                window_is_empty=window_is_empty,
                gt_region=gt_region,
                popindices=popindices,
                chromosome=chromosome,
                window_pos_1=window_pos_1,
                window_pos_2=window_pos_2,
            )

            pixy_output.extend(dxy_results)

        # FST:
        # WEIR AND COCKERHAMS FST
        # This is just a loose wrapper around the scikit-allel fst function
        # TBD: explicit fst when data is completely missing
        if (
            (args.populations is not None)
            and ("fst" in args.stats)
            and window_size != 1
            and not skip_fst_for_chrom
        ):
            # disabling these assertions for now, because they are handled by "window_is_empty"
            # i.e. genotype arrays CAN be empty.
            # assert gt_array_fst is not None, "genotype array is None"

            fst_type: FSTEstimator = FSTEstimator(args.fst_type)
            fst_results: List[PixyTempResult] = compute_summary_fst(
                fst_type=fst_type,
                gt_array_fst=gt_array_fst,
                pos_array_fst=pos_array_fst,
                window_is_empty=window_is_empty,
                callset_is_none=callset_is_none,
                aggregate=aggregate,
                popindices=popindices,
                chromosome=chromosome,
                window_pos_1=window_pos_1,
                window_pos_2=window_pos_2,
            )

            pixy_output.extend(fst_results)
        # TODO: fix so that passed args should be a `PixyArgs` object
        if (args.populations is not None) and ("tajima_d" in args.stats):
            for pop in popnames:
                # if the window has no sites in the VCF, assign all NAs,
                # otherwise calculate Tajima's D
                if window_is_empty:
                    tajima_result: TajimaDResult = TajimaDResult.empty()
                else:
                    # subset the window for the individuals in each population
                    assert gt_region is not None, "gt_region is None"
                    gt_pop = gt_region.take(popindices[pop], axis=1)

                    # if the population specific window for this region is empty, report it as such
                    if len(gt_pop) == 0:
                        tajima_result = TajimaDResult.empty()

                    # otherwise compute Tajima's D as normal
                    else:
                        tajima_result = calc_tajima_d(gt_pop)
                # consult the docstring of `PixyTempResult` for more details on overloaded fields
                pixy_results: PixyTempResult = PixyTempResult(
                    pixy_stat=PixyStat.TAJIMA_D,
                    population_1=pop,
                    population_2="NA",
                    chromosome=chromosome,
                    window_pos_1=window_pos_1,
                    window_pos_2=window_pos_2,
                    calculated_stat=tajima_result.tajima_d,
                    shared_sites_with_alleles=tajima_result.num_sites,
                    total_differences=tajima_result.raw_pi,
                    total_comparisons=tajima_result.watterson_theta,
                    total_missing=tajima_result.d_stdev,
                    tajima_d_variant_counts=serialize_tajima_d_variant_counts(
                        tajima_result.variant_gt_counts
                    ),
                )
                pixy_output.append(pixy_results)

        # WATTERSON'S THETA:
        # GENETIC DIVERSITY CALCULATED FROM NUMBER OF SEGREGATING (VARIANT) SITES

        if (args.populations is not None) and ("watterson_theta" in args.stats):
            for pop in popnames:
                # if the window has no sites in the VCF, assign all NAs,
                # otherwise calculate Watterson's theta
                if window_is_empty:
                    watterson_result = WattersonThetaResult.empty()
                else:
                    assert gt_region is not None, "genotype array is None"
                    # subset the window for the individuals in each population
                    gt_pop = gt_region.take(popindices[pop], axis=1)

                    # if the population specific window for this region is empty, report it as such
                    if len(gt_pop) == 0:
                        watterson_result = WattersonThetaResult.empty()
                    # otherwise compute Watterson's theta as normal
                    else:
                        watterson_result = calc_watterson_theta(gt_pop)
                    # consult the docstring of `PixyTempResult`
                    # for more details on overloaded fields
                pixy_results = PixyTempResult(
                    pixy_stat=PixyStat.WATTERSON_THETA,
                    population_1=pop,
                    population_2="NA",
                    chromosome=chromosome,
                    window_pos_1=window_pos_1,
                    window_pos_2=window_pos_2,
                    calculated_stat=watterson_result.avg_theta,
                    shared_sites_with_alleles=watterson_result.num_sites,
                    total_differences=watterson_result.raw_theta,
                    total_comparisons=watterson_result.num_var_sites,
                    total_missing=watterson_result.num_weighted_sites,
                )
                pixy_output.append(pixy_results)

        # OUTPUT
        # if in mc mode, put the results in the writing queue
        # otherwise just write to file

    # ensure the output variable exists in some form
    if pixy_output:
        temp_pixy_content: str = "\n".join(f"{result}" for result in pixy_output)
        if args.n_cores > 1:
            assert isinstance(q, (queues.Queue, BaseProxy)), "q should be a Queue BaseProxy"
            q.put(temp_pixy_content)
        elif args.n_cores == 1:
            outfile = open(temp_file, "a")
            outfile.write(temp_pixy_content + "\n")
            outfile.close()

