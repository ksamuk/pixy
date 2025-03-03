import argparse
import logging
import os
import shutil
import subprocess
import uuid
from pathlib import Path
from typing import Dict
from typing import List
from typing import Literal
from typing import Optional
from typing import Set
from typing import Tuple
from typing import Union

import allel
import multiprocess as mp
import numpy as np
import pandas
from allel import GenotypeArray
from allel import GenotypeVector
from allel import SortedIndex
from multiprocess import Queue
from numpy.typing import NDArray

from pixy.args_validation import get_chrom_list
from pixy.args_validation import validate_bed_path
from pixy.args_validation import validate_output_path
from pixy.args_validation import validate_populations_path
from pixy.args_validation import validate_sites_path
from pixy.args_validation import validate_vcf_path
from pixy.args_validation import validate_window_and_interval_args
from pixy.calc import calc_tajima_d
from pixy.calc import calc_watterson_theta
from pixy.enums import FSTEstimator
from pixy.enums import PixyStat
from pixy.models import PixyTempResult
from pixy.models import TajimaDResult
from pixy.models import WattersonThetaResult
from pixy.stats.summary import compute_summary_dxy
from pixy.stats.summary import compute_summary_fst
from pixy.stats.summary import compute_summary_pi
from pixy.stats.summary import precompute_filtered_variant_array


def aggregate_output(
    df_for_stat: pandas.DataFrame,
    stat: str,
    chromosome: str,
    window_size: int,
    fst_type: str,
) -> pandas.DataFrame:
    """
    Aggregates genomic data into windows and computes summary statistics over each window.

    Summary statistics could be one or more of pi (genetic diversity within a population), dxy
    (genetic diversity between populations), or fst (proportion of the total genetic variance across
    a subpopulation).

    Assumes that column 4 is the position of the variant; columns 7-10 are the alleles
    counts/related metrics.

    Empty statistics are marked as 'NA'. The computed statistics are grouped by population and
    window position in the resultant Dataframe.

    Args:
        df_for_stat: Contains genomic data with columns for position,
            allele counts, and population information.
        stat: The statistic to compute.
        chromosome: The name of the chromosome for the current window.
        window_size: The size of the genomic windows (in base pairs) for aggregation.
        fst_type: The method for calculating fst. One of:
            - 'wc' for Weir and Cockerham's fst
            - 'hudson' for Hudson's fst

    Returns:
        outsorted: a pandas.DataFrame that stores aggregated statistics for each window
    """
    outsorted = df_for_stat.sort_values([4])  # sort by position
    interval_start = df_for_stat[4].min()
    interval_end = df_for_stat[4].max()
    # create a subset of the window list specific to this chunk
    bins = np.arange(interval_start - 1, interval_end + window_size, window_size)
    # bin values into discrete intervals, both left and right inclusive
    assignments, edges = pandas.cut(
        outsorted[4], bins=bins, labels=False, retbins=True, include_lowest=True
    )
    outsorted["label"] = assignments
    outsorted["window_pos_1"] = edges[assignments] + 1
    outsorted["window_pos_2"] = edges[assignments + 1]

    # group by population, window
    if stat == "pi" or stat == "tajima_d":  # pi and tajima_d only have one population field
        outsorted = (
            outsorted.groupby([1, "window_pos_1", "window_pos_2"], as_index=False, dropna=False)
            .agg({7: "sum", 8: "sum", 9: "sum", 10: "sum"})
            .reset_index()
        )
    elif stat == "dxy" or stat == "fst":  # dxy and fst have 2 population fields
        outsorted = (
            outsorted.groupby([1, 2, "window_pos_1", "window_pos_2"], as_index=False, dropna=False)
            .agg({7: "sum", 8: "sum", 9: "sum", 10: "sum"})
            .reset_index()
        )

    if stat == "pi" or stat == "dxy" or stat == "watterson_theta" or stat == "tajima_d":
        outsorted[stat] = outsorted[8] / outsorted[9]
    elif stat == "fst":
        if fst_type == "wc":
            outsorted[stat] = outsorted[8] / (outsorted[8] + outsorted[9] + outsorted[10])
        elif fst_type == "hudson":
            # 'a' is the numerator of hudson and 'b' is the denominator
            # (there is no 'c')
            outsorted[stat] = outsorted[8] / (outsorted[9])

    outsorted[stat].fillna("NA", inplace=True)
    outsorted["chromosome"] = chromosome

    # reorder columns
    if (
        stat == "pi" or stat == "watterson_theta" or stat == "tajima_d"
    ):  # pi, Watterson's theta and Tajima's D only have one population field
        outsorted = outsorted[[1, "chromosome", "window_pos_1", "window_pos_2", stat, 7, 8, 9, 10]]
    else:  # dxy and fst have 2 population fields
        outsorted = outsorted[
            [1, 2, "chromosome", "window_pos_1", "window_pos_2", stat, 7, 8, 9, 10]
        ]

    # make sure sites, comparisons, missing get written as floats
    if stat == "pi" or stat == "dxy" or stat == "watterson_theta" or stat == "tajima_d":
        cols = [7, 8, 9, 10]
    elif stat == "fst":
        cols = [7]
    outsorted[cols] = outsorted[cols].astype("int32")
    return outsorted


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

        sub_windows = [list(a) for a in zip(subwindow_pos_1_list, subwindow_pos_2_list)]
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
        list(a) for a in zip(window_pos_1_list, window_pos_2_list, chunk_list_1, chunk_list_2)
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
    sites_list = [list(a) for a in zip(sites_pre_list, chunk_list)]

    return sites_list


# function for masking non-target sites in a genotype array
# used in conjuctions with the --sites_list command line option
# NB: `SortedIndex` is of type `int`, but is not generic (at least in 3.8)
def mask_non_target_sites(
    gt_array: GenotypeArray, pos_array: SortedIndex, sites_list_chunk: List[int]
) -> GenotypeArray:
    """
    Masks non-target sites in a genotype array.

    Masked sites are set to missing data (`-1`).

    Args:
        gt_array: the `GenotypeArray` of interest
        pos_array: positions corresponding to the sites in `gt_array`
        sites_list_chunk: target site positions that should remain unmasked

    Returns:
        gt_array: a modified `GenotypeArray`
    """
    # get the indexes of sites that are NOT the target sites
    # (these will be masked with missing rows to remove them from the calculations)
    masked_sites_set: Set[int] = set(pos_array) ^ set(sites_list_chunk)
    masked_sites: List[int] = sorted(list(set(pos_array).intersection(masked_sites_set)))

    gt_mask_indexes: List[int] = list(np.flatnonzero(pos_array.locate_keys(masked_sites)))

    # a missing row of data to use as a mask
    missing_row: List[List[int]] = [[-1] * gt_array.ploidy] * gt_array.n_samples

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
) -> Tuple[bool, Optional[GenotypeArray], Optional[SortedIndex]]:
    """
    Ingests genotypes from a VCF file, retains biallelic SNPs or invariant sites.

    Filters out non-SNPs, multi-allelic SNPs, and non-variant sites. Optionally masks out
    non-target sites based on a provided list (`sites_list_chunk`).

    Variants for which the depth of coverage (`DP`) is less than 1 are considered to be missing
    and replaced with `-1`.

    If the VCF contains no variants over the specified genomic region, sets `callset_is_none` to
    `True`.

    Args:
        args (Namespace): Command-line arguments, should include the path to the VCF file as
            `args.vcf`
        chromosome (str): The chromosome for which to read the genotype data (e.g., 'chr1')
        window_pos_1 (int): The start position of the genomic window
        window_pos_2 (int): The end position of the genomic window
        sites_list_chunk (list): An optional list of positions in which to mask non-target sites
            If `None`, no masking is applied.

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

    include_multiallelic_snps: bool = args.include_multiallelic_snps == "yes"

    # read in data from the source VCF for the current window
    callset = allel.read_vcf(
        args.vcf,
        region=window_region,
        fields=[
            "CHROM",
            "POS",
            "calldata/GT",
            "calldata/DP",
            "variants/is_snp",
            "variants/numalt",
        ],
    )

    # keep track of whether the callset was empty (no sites for this range in the VCF)
    # used by compute_summary_stats to add info about completely missing sites
    if callset is None:
        callset_is_none = True
        gt_array = None
        pos_array = None

    else:
        # if the callset is NOT empty (None), continue with pipeline
        callset_is_none = False

        # fix for cursed GATK 4.0 missing data representation
        # forces DP<1 (zero) to be missing data (-1 in scikit-allel)
        # BROKEN HERE <- add if statement to check if DP info is present!
        callset["calldata/GT"][callset["calldata/DP"] < 1, :] = -1

        # convert to a genotype array object
        gt_array = allel.GenotypeArray(allel.GenotypeDaskArray(callset["calldata/GT"]))

        # build an array of positions for the region
        pos_array = allel.SortedIndex(callset["variants/POS"])

        # create a mask for biallelic snps and invariant sites
        is_biallelic_snp = np.logical_and(
            callset["variants/is_snp"][:] == 1,
            callset["variants/numalt"][:] == 1,
        )
        is_multiallelic_snp = np.logical_and(
            callset["variants/is_snp"][:] == 1,
            callset["variants/numalt"][:] > 1,
        )
        is_invariant_site = callset["variants/numalt"][:] == 0

        snp_invar_mask = np.logical_or(
            is_biallelic_snp,
            is_invariant_site,
            np.logical_and(include_multiallelic_snps, is_multiallelic_snp),
        )

        # remove rows that are NOT snps or invariant sites from the genotype array
        gt_ndarray = np.delete(gt_array, np.where(np.invert(snp_invar_mask)), axis=0)
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

    # read in the genotype data for the chunk
    callset_is_none, gt_array, pos_array = read_and_filter_genotypes(
        args, chromosome, chunk_pos_1, chunk_pos_2, sites_list_chunk
    )

    # if computing FST, pre-compute a filtered array of variants (only)
    if "fst" in args.stats:
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

            # double check that the region is not empty after subsetting
            try:
                loc_region  # noqa: B018
            except Exception:
                loc_region = None

            try:
                gt_region  # noqa: B018
            except Exception:
                gt_region = None

            if (gt_region is None) or len(gt_region) == 0 or (loc_region is None):
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
        if (args.populations is not None) and ("fst" in args.stats) and window_size != 1:
            assert gt_array_fst is not None, "genotype array is None"

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
            assert gt_region is not None, "gt_region is None"
            for pop in popnames:
                # if the window has no sites in the VCF, assign all NAs,
                # otherwise calculate Tajima's D
                if window_is_empty:
                    tajima_result: TajimaDResult = TajimaDResult.empty()
                else:
                    # subset the window for the individuals in each population
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
                )
                pixy_output.append(pixy_results)

        # WATTERSON'S THETA:
        # GENETIC DIVERSITY CALCULATED FROM NUMBER OF SEGREGATING (VARIANT) SITES

        if (args.populations is not None) and ("watterson_theta" in args.stats):
            assert gt_region is not None, "genotype array is None"
            assert callset_is_none is not None, "callset is empty"
            for pop in popnames:
                # if the window has no sites in the VCF, assign all NAs,
                # otherwise calculate Watterson's theta
                if window_is_empty:
                    watterson_result = WattersonThetaResult.empty()
                else:
                    # subset the window for the individuals in each population
                    gt_pop = GenotypeArray(gt_region.take(popindices[pop], axis=1))
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
            assert isinstance(q, Queue), "q should be a Queue object when multiprocessing"
            q.put(temp_pixy_content)
        elif args.n_cores == 1:
            outfile = open(temp_file, "a")
            outfile.write(temp_pixy_content + "\n")
            outfile.close()


# function for checking & validating command line arguments
# when problems are detected, throws an error/warning + message


def check_and_validate_args(  # noqa: C901
    args: argparse.Namespace,
) -> Tuple[
    NDArray,
    Dict[str, NDArray],
    List[str],
    List[str],
    Path,
    str,
    str,
    pandas.DataFrame,
    pandas.DataFrame,
]:
    """
    Checks whether user-specific arguments are valid.

    Args:
        args: parsed CLI args specified by the user

    Raises:
        Exception: if the output_folder is not writeable
        Exception: if any of the required input files are missing
        Exception: if the output_prefix contains either forward or backward slashes
        Exception: if any of the provided files do not exist
        Exception: if the VCF file is not compressed with `bgzip` or indexed with `tabix`
        Exception: if the VCF file does not contain variant sites
        Exception: if the provided `--chromosomes` do not occur in the specified VCF file
        Exception: if neither a `--bed-file` nor a `--window_size` is provided
        Exception: if only one of `--interval_start` and `--interval_end` is given
        Exception: if multiple `--chromosomes` and an interval are provided
        Exception: if any rows in either the `--bed-file` or `--sites-path` are missing data
        Exception: if any of the samples provided in the `populations_file` do not exist in the VCF
        Exception: if the `populations_file` does not contain at least 2 populations

    Returns:
        popnames: unique population names derived from `--populations` file
        popindices: the indices of the population-specific samples within the provided VCF
        chrom_list: list of chromosomes that are provided and also found in the provided VCF
        IDs: list of each individual in the `--populations` file
        temp_file: a temporary file in which to hold in-flight `pixy` calculations
        output_folder: the directory to which to write any `pixy` results
        output_prefix: the combination of a given `output_folder` and `output_prefix`
        bed_df: a pandas.DataFrame that represents a given `--bed-file` (the dataframe will be empty
            if no BED file was given)
        sites_df: a pandas.DataFrame that represents a given `--sites-file` (the dataframe will be
            empty if no sites file was given)

    """
    # CHECK FOR TABIX
    logger = logging.getLogger(__name__)
    tabix_path = shutil.which("tabix")

    if tabix_path is None:
        raise Exception(
            "[pixy] ERROR: tabix is not installed (or cannot be located in the path). "
            'Install tabix with "conda install -c bioconda htslib".'
        )

    if args.vcf is None:
        raise Exception("[pixy] ERROR: The --vcf argument is missing or incorrectly specified.")

    if args.populations is None:
        raise Exception(
            "[pixy] ERROR: The --populations argument is missing or incorrectly specified."
        )

    # reformat file paths for compatibility

    populations_path: Path = Path(os.path.expanduser(args.populations))
    populations_df: pandas.DataFrame = validate_populations_path(populations_path)

    vcf_path: str = os.path.expanduser(args.vcf)  # we don't want a Path object just yet because
    # most of the downstream operations require a string
    validate_vcf_path(vcf_path)

    if args.output_folder != "":
        output_folder = args.output_folder + "/"
    else:
        output_folder = os.path.expanduser(os.getcwd() + "/")

    output_prefix = output_folder + args.output_prefix

    # get vcf header info
    vcf_headers = allel.read_vcf_headers(vcf_path)

    logger.info("\n[pixy] Validating VCF and input parameters...")

    # CHECK OUTPUT FOLDER
    logger.info("[pixy] Checking write access...")

    output_folder, output_prefix = validate_output_path(
        output_folder=args.output_folder, output_prefix=args.output_prefix
    )
    # generate a name for a unique temp file for collecting output
    temp_file: Path = Path(output_folder) / f"pixy_tmpfile_{uuid.uuid4().hex}.tmp"
    # check if temp file is writable
    with open(temp_file, "w"):
        pass  # file is created and then closed
    assert os.access(temp_file, os.W_OK), "temp file is not writable"

    # CHECK CPU CONFIGURATION
    logger.info("[pixy] Checking CPU configuration...")

    if args.n_cores > mp.cpu_count():
        logger.warning(
            f"[pixy] WARNING: {args.n_cores} CPU cores requested but only "
            f"{mp.cpu_count()} are available. Using {mp.cpu_count()}."
        )
        args.n_cores = mp.cpu_count()

    # CHECK FOR EXISTENCE OF INPUT FILES

    if args.bed_file is not None:
        bed_path: Path = Path(os.path.expanduser(args.bed_file))
        bed_df: pandas.DataFrame = validate_bed_path(bed_path)

    else:
        bed_df = []

    # VALIDATE THE VCF

    # check if the vcf contains any invariant sites
    # a very basic check: just looks for at least one invariant site in the alt field
    logger.info("[pixy] Checking for invariant sites...")

    if args.bypass_invariant_check == "no":
        alt_list = (
            subprocess.check_output(
                "gunzip -c "
                + vcf_path
                + " | grep -v '#' | head -n 100000 | awk '{print $5}' | sort | uniq",
                shell=True,
            )
            .decode("utf-8")
            .split()
        )
        if "." not in alt_list:
            raise Exception(
                "[pixy] ERROR: the provided VCF appears to contain no invariant sites "
                '(ALT = "."). '
                "This check can be bypassed via --bypass_invariant_check 'yes'."
            )
        if "." in alt_list and len(alt_list) == 1:
            logger.warning(
                "[pixy] WARNING: the provided VCF appears to contain no variable sites in the "
                "first 100 000 sites. It may have been filtered incorrectly, or genetic diversity "
                "may be extremely low. "
                "This warning can be suppressed via --bypass_invariant_check 'yes'.'"
            )
    else:
        if not (len(args.stats) == 1 and (args.stats[0] == "fst")):
            logger.warning(
                "[pixy] EXTREME WARNING: --bypass_invariant_check is set to 'yes'. Note that a "
                "lack of invariant sites will result in incorrect estimates."
            )

    # check if requested chromosomes exist in vcf
    # parses the whole CHROM column (!)

    logger.info("[pixy] Checking chromosome data...")

    chrom_list: List[str] = get_chrom_list(args)

    # INTERVALS
    # check if intervals are correctly specified
    # validate the BED file (if present)

    logger.info("[pixy] Checking intervals/sites...")

    if args.bed_file is None:
        check_message = validate_window_and_interval_args(args)
        logger.info(f"{check_message}")
    else:
        if (
            args.interval_start is not None
            or args.interval_end is not None
            or args.window_size is not None
        ):
            raise Exception(
                "[pixy] ERROR: --interval_start, --interval_end, and --window_size are not valid "
                "when a BED file of windows is provided."
            )

        bed_chrom: List[str] = list(bed_df["chrom"])
        missing = list(set(bed_chrom) - set(chrom_list))
        chrom_list = list(set(chrom_list) & set(bed_chrom))

        if len(missing) > 0:
            logger.warning(
                "[pixy] WARNING: the following chromosomes in the BED file do not occur in the VCF "
                f"and will be ignored: {missing}"
            )
    if args.sites_file is None:
        sites_df: pandas.DataFrame = pandas.DataFrame([])
        chrom_sites: List[str] = []
        missing_sites: List[str] = []

    else:
        sites_path: Path = Path(os.path.expanduser(args.sites_file))
        sites_df = validate_sites_path(sites_path=sites_path)

        # all the chromosomes in the sites file
        chrom_sites = list(sites_df["CHROM"])

        # the difference between the chromosomes in the sites file and the VCF
        missing_sites = list(set(chrom_sites) - set(chrom_list))

        if len(missing_sites) > 0:
            logger.warning(
                "[pixy] WARNING: the following chromosomes in the sites file do not occur in the "
                f"VCF and will be ignored: {missing_sites}"
            )

    # SAMPLES
    # check if requested samples exist in vcf

    logger.info("[pixy] Checking sample data...")

    # - parse + validate the population file
    # - format is IND POP (tab separated)
    # - throws an error if individuals are missing from VCF

    # get a list of samples from the callset
    samples_list = vcf_headers.samples
    # make sure every indiv in the pop file is in the VCF callset
    sample_ids: List[str] = list(populations_df["ID"])
    missing = list(set(sample_ids) - set(samples_list))

    # find the samples in the callset index by matching up the order of samples between the
    # population file and the callset
    # also check if there are invalid samples in the popfile
    try:
        samples_callset_index = [samples_list.index(s) for s in populations_df["ID"]]
    except ValueError as e:
        raise Exception(
            "[pixy] ERROR: the following samples are listed in the population file "
            f"but not in the VCF: {missing}"
        ) from e
    else:
        populations_df["callset_index"] = samples_callset_index

        # use the popindices dictionary to keep track of the indices for each population
        popindices = {}
        popnames = populations_df.Population.unique()
        for name in popnames:
            popindices[name] = populations_df[
                populations_df.Population == name
            ].callset_index.values

    if len(popnames) == 1 and ("fst" in args.stats or "dxy" in args.stats):
        raise Exception(
            "[pixy] ERROR: calculation of fst and/or dxy requires at least two populations to be "
            "defined in the population file."
        )

    logger.info("[pixy] All initial checks passed!")

    return (
        popnames,
        popindices,
        chrom_list,
        sample_ids,
        temp_file,
        output_folder,
        output_prefix,
        bed_df,
        sites_df,
    )
