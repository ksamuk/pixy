import argparse
import warnings
from multiprocessing import Queue
from multiprocessing import queues
from multiprocessing.managers import BaseProxy
from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Optional
from typing import Tuple
from typing import Union

import allel
import numpy as np
from allel import GenotypeArray
from allel import GenotypeVector
from allel import SortedIndex
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
from pixy.wisp import WindowInvariantContributions
from pixy.wisp import WispMask
from pixy.wisp import compute_window_invariant_contributions

# Note: the previous pandas-based aggregation pipeline (`aggregate_output()` plus the helpers
# `_group_aggregated_output`, `_calculate_aggregated_stat`, `_coerce_aggregated_output_types`,
# `_aggregate_tajima_d_variant_counts`) has been replaced by the streaming pure-Python
# implementation in `pixy.agg.aggregate_rows()` / `pixy.agg.write_stat_file()`. Those used
# `pandas.read_csv` + `groupby` + `to_csv` to round-trip the tmp file through a DataFrame,
# which scaled memory linearly with the window count; the new path keeps a dict accumulator
# keyed by (pop1, pop2?, bin_idx) and emits output line-by-line.


def temp_file_listener(q: "Queue", temp_file: str) -> None:
    """
    Drain worker-produced rows from a queue and append them to the temp file.

    Lives in `pixy.core` rather than `pixy.__main__` because workers launched under
    `forkserver` or `spawn` start a fresh interpreter where `__main__` is the worker
    bootstrap, not pixy's CLI entrypoint — so `pickle` cannot resolve attributes of the
    parent's `__main__` module by name. Functions targeted by `pool.apply_async` must live
    in an importable module (which `pixy.core` is) for cross-process pickling to round-trip.
    """
    with open(temp_file, "a") as f:
        while True:
            m = q.get()
            if m == "kill":
                break
            f.write(str(m) + "\n")
            f.flush()  # immediately write data, do not buffer


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


# Default sub-chunk size for the streaming variant-only read path. Matches scikit-allel's
# own DEFAULT_CHUNK_LENGTH so per-iteration parse overhead is unchanged from a single
# `read_vcf` call; the win is that invariant rows are dropped per-sub-chunk, so peak RSS
# scales with `sub_chunk_length + variant_density * full_chunk_length` rather than
# `full_chunk_length`.
_STREAMING_SUB_CHUNK_LENGTH: int = 65536


def _read_filtered_variants_streaming(  # noqa: C901
    vcf_path: str,
    window_region: str,
    ploidy: int,
    include_multiallelic_snps: bool,
) -> Tuple[bool, Optional[GenotypeArray], Optional[SortedIndex]]:
    """
    Stream a VCF region in sub-chunks, keeping only variant rows in memory.

    Used by `read_and_filter_genotypes` when invariant sites are not needed (FST-only
    runs, or any path where the VCF is variants-only by design — e.g. the wisp-mask
    path). For each sub-chunk yielded by `allel.iter_vcf_chunks`, rows that are not
    biallelic SNPs (or multi-allelic SNPs when `include_multiallelic_snps`) are
    dropped before being appended to the accumulator. Peak RSS therefore scales with
    the variant density of the chunk rather than its total length.
    """
    pos_pieces: List[NDArray] = []
    gt_pieces: List[NDArray] = []

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, module="allel")
        _fields, _samples, _headers, it = allel.iter_vcf_chunks(  # type: ignore[attr-defined]
            vcf_path,
            fields=[
                "CHROM",
                "POS",
                "calldata/GT",
                "variants/is_snp",
                "variants/numalt",
            ],
            numbers={"GT": ploidy},
            region=window_region,
            chunk_length=_STREAMING_SUB_CHUNK_LENGTH,
        )

        for chunk_dict, _chunk_length, _chrom, _last_pos in it:
            is_biallelic_snp = np.logical_and(
                chunk_dict["variants/is_snp"][:] == 1,
                chunk_dict["variants/numalt"][:] == 1,
            )
            snp_mask = is_biallelic_snp
            if include_multiallelic_snps:
                is_multiallelic_snp = np.logical_and(
                    chunk_dict["variants/is_snp"][:] == 1,
                    chunk_dict["variants/numalt"][:] > 1,
                )
                snp_mask = np.logical_or(snp_mask, is_multiallelic_snp)

            if not snp_mask.any():
                continue

            pos_pieces.append(chunk_dict["variants/POS"][snp_mask])
            gt_pieces.append(chunk_dict["calldata/GT"][snp_mask])

    if not pos_pieces:
        return True, None, None

    pos_concat = np.concatenate(pos_pieces)
    gt_concat = np.concatenate(gt_pieces, axis=0)

    gt_array: GenotypeArray
    if ploidy == 1:
        gt_array = allel.HaplotypeArray(gt_concat)
    else:
        gt_array = allel.GenotypeArray(gt_concat)
    pos_array = allel.SortedIndex(pos_concat)
    return False, gt_array, pos_array


# function for reading in a genotype matrix from a VCF file
# also filters out all but biallelic SNPs and invariant sites
# returns a genotype matrix, and array of genomic coordinates and
# a logical indicating whether the array(s) are empty
def read_and_filter_genotypes(  # noqa: C901
    args: argparse.Namespace,
    chromosome: str,
    window_pos_1: int,
    window_pos_2: int,
    sites_list_chunk: Optional[List[int]],
    ploidy: int,
    needs_invariants: bool = True,
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
        needs_invariants (bool): Whether to include invariant sites in the returned arrays.
            Should be ``True`` when computing pi, dxy, tajima_d, or watterson_theta (all of
            which use callable-site counts in their denominators), and ``False`` when computing
            FST only (which only needs variant sites).

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

    # Pre-declare the array vars so mypy widens to the union type rather than narrowing to
    # whichever subclass appears first below (HaplotypeArray vs GenotypeArray).
    gt_array: Optional[GenotypeArray]
    pos_array: Optional[SortedIndex]

    # When invariant sites are not needed (FST-only, or any path where the VCF is
    # variants-only by design), stream the region in sub-chunks and drop invariant
    # rows per-sub-chunk so peak RSS scales with variant density rather than chunk
    # length. `sites_list_chunk` masking, if any, is applied below by the shared
    # post-read path; the streaming branch returns early with the assembled arrays.
    if not needs_invariants:
        callset_is_none, gt_array, pos_array = _read_filtered_variants_streaming(
            vcf_path=args.vcf,
            window_region=window_region,
            ploidy=ploidy,
            include_multiallelic_snps=include_multiallelic_snps,
        )
        if not callset_is_none and sites_list_chunk is not None:
            assert pos_array is not None and gt_array is not None
            gt_array = mask_non_target_sites(gt_array, pos_array, sites_list_chunk)
            if len(gt_array) == 0:
                return True, None, None
        return callset_is_none, gt_array, pos_array

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

        # create a mask for biallelic SNPs; invariant sites are added only when needed
        # (pi, dxy, tajima_d, watterson_theta all require them for their denominators;
        # FST only needs variant sites so can skip this step for a meaningful speedup)
        is_biallelic_snp = np.logical_and(
            callset["variants/is_snp"][:] == 1,
            callset["variants/numalt"][:] == 1,
        )

        snp_invar_mask = is_biallelic_snp
        if needs_invariants:
            is_invariant_site = callset["variants/numalt"][:] == 0
            snp_invar_mask = np.logical_or(snp_invar_mask, is_invariant_site)

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
    # `q` is the manager-provided proxy (`multiprocessing.managers.SyncManager.Queue()`) in
    # multicore mode, or the literal string "NULL" in single-core mode. The proxy's type is
    # `queue.Queue` (NOT `multiprocessing.queues.Queue`), and the two type stubs are
    # incompatible — typing this as `Any` accepts both forms; the function does an explicit
    # `isinstance` check on use.
    q: Any,
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

    # Invariant sites are only needed for pi, dxy, tajima_d, and watterson_theta, which use
    # callable-site counts in their denominators. FST only needs variant sites, so skip
    # invariant ingestion when FST is the sole requested statistic.
    #
    # When a wisp mask is supplied, the invariant denominator comes from the mask
    # instead of from VCF invariants — the input VCF is variants-only by design, so
    # `needs_invariants` should be False even for pi/dxy/tajima_d/watterson_theta.
    wisp_mask: Optional[WispMask] = getattr(args, "wisp_mask", None)
    stats_needing_invariants = {"pi", "dxy", "tajima_d", "watterson_theta"}
    needs_invariants: bool = (
        bool(stats_needing_invariants.intersection(args.stats)) and wisp_mask is None
    )

    # read in the genotype data for the chunk
    callset_is_none, gt_array, pos_array = read_and_filter_genotypes(
        args,
        chromosome,
        chunk_pos_1,
        chunk_pos_2,
        sites_list_chunk,
        chrom_ploidy,
        needs_invariants,
    )

    # if computing FST, pre-compute a filtered array of variants (only)
    # Initialize to None so the window loop below can reference them even when the
    # callset is empty (no VCF records in this chunk — common in sparse RAD-seq data).
    gt_array_fst: Optional[GenotypeArray] = None
    pos_array_fst: Optional[SortedIndex] = None
    if "fst" in args.stats and not skip_fst_for_chrom and not callset_is_none:
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

    # Whether the wisp path is active (drives per-window invariant-contribution lookups
    # and the post-hoc Watterson/Tajima merges further down).
    wisp_active: bool = wisp_mask is not None
    pop_names_list: List[str] = [str(p) for p in popnames]

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
            assert gt_region is not None  # narrow for the len() call below

            # An empty slice after subsetting means the window has no usable sites.
            if len(gt_region) == 0:
                window_is_empty = True

        # Wisp invariant contribution for this window (computed once and reused
        # across pi/dxy/tajima/watterson). When the wisp mask is not active, this
        # stays None and the consumers behave exactly as before.
        invariant_contribution: Optional[WindowInvariantContributions] = None
        if wisp_active:
            assert wisp_mask is not None  # narrowed by wisp_active
            # Variant positions within the current window (1-based, pixy convention),
            # used to subtract variant sites from each wisp range when computing the
            # invariant count. `pos_array` is the entire chunk; subset to the window.
            if pos_array is None:
                window_var_positions: List[int] = []
            else:
                in_window = (pos_array >= window_pos_1) & (pos_array <= window_pos_2)
                window_var_positions = [int(p) for p in np.asarray(pos_array)[in_window]]
            invariant_contribution = compute_window_invariant_contributions(
                wisp_mask=wisp_mask,
                chromosome=chromosome,
                window_pos_1=window_pos_1,
                window_pos_2=window_pos_2,
                variant_positions=window_var_positions,
                pop_names=pop_names_list,
                ploidy=chrom_ploidy,
            )
            # If the window had no VCF variants AND we got at least one invariant site
            # back from the mask, treat the window as non-empty: the wisp contribution
            # supplies a real denominator. Without this flip the consumers would NA-out
            # everything and discard the invariant-only data.
            if window_is_empty:
                any_inv = any(p.num_sites > 0 for p in invariant_contribution.pi.values()) or any(
                    d.num_sites > 0 for d in invariant_contribution.dxy.values()
                )
                if any_inv:
                    window_is_empty = False
                    # Build a zero-row genotype region so the per-stat code paths don't
                    # crash on `gt_region.take(...)`. Each compute_summary_X function
                    # handles the empty-population case by emitting NA, after which the
                    # invariant_contribution merge promotes the totals.
                    if gt_array is not None and len(gt_array) > 0:
                        gt_region = gt_array[0:0]
                    else:
                        # No variant rows in this chunk to slice from (sparse VCF, or
                        # the whole chunk is variant-free). Construct a zero-row
                        # GenotypeArray of (0, n_samples, ploidy) so .take(popindices[pop])
                        # still returns an empty pop-specific array. We route the
                        # construction through [0:0] so it picks up the same
                        # type-narrowing as the slice-from-existing branch above.
                        max_idx = -1
                        for idx in popindices.values():
                            if len(idx) > 0:
                                max_idx = max(max_idx, int(np.max(idx)))
                        n_samples = max_idx + 1 if max_idx >= 0 else 0
                        empty_full = GenotypeArray(
                            np.empty((0, n_samples, chrom_ploidy), dtype=np.int8)
                        )
                        gt_region = empty_full[0:0]

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
                invariant_contribution=invariant_contribution,
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
                invariant_contribution=invariant_contribution,
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
                # Add the wisp invariant-site contribution to num_sites only.
                # Invariants don't contribute to raw_pi, Watterson's theta, the per-site
                # variant-class counts, or the d_stdev — those are all variant-only quantities.
                num_sites = tajima_result.num_sites
                if invariant_contribution is not None:
                    inv_t = invariant_contribution.tajima.get(str(pop))
                    if inv_t is not None:
                        num_sites = num_sites + inv_t.num_sites
                # consult the docstring of `PixyTempResult` for more details on overloaded fields
                pixy_results: PixyTempResult = PixyTempResult(
                    pixy_stat=PixyStat.TAJIMA_D,
                    population_1=pop,
                    population_2="NA",
                    chromosome=chromosome,
                    window_pos_1=window_pos_1,
                    window_pos_2=window_pos_2,
                    calculated_stat=tajima_result.tajima_d,
                    shared_sites_with_alleles=num_sites,
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
                # Merge in the wisp invariant contribution. Invariants don't add to
                # raw_theta or to num_var_sites (no variants by definition), but they DO
                # add to num_sites and reshape `num_weighted_sites` because that quantity
                # depends on the window-wide `max(k_haps)` across both variants and
                # invariants. We recompute weighted_sites here from the union distribution.
                num_sites_combined = watterson_result.num_sites
                raw_theta = watterson_result.raw_theta
                num_var_sites = watterson_result.num_var_sites
                num_weighted = watterson_result.num_weighted_sites
                avg_theta = watterson_result.avg_theta
                if invariant_contribution is not None:
                    inv_w = invariant_contribution.watterson.get(str(pop))
                    if inv_w is not None and inv_w.num_sites > 0:
                        # Combine k-distributions. The variant side is available from the
                        # genotype array, but recomputing it costs little and keeps the
                        # logic localized; reach for it only when we actually have an
                        # invariant contribution to merge in.
                        from collections import Counter

                        k_dist: "Counter[int]" = Counter(inv_w.k_haps_counter)
                        if not window_is_empty and gt_region is not None:
                            gt_pop = gt_region.take(popindices[pop], axis=1)
                            if len(gt_pop) > 0:
                                allele_counts = gt_pop.count_alleles()
                                for k_haps in np.asarray(allele_counts).sum(axis=1):
                                    if k_haps > 0:
                                        k_dist[int(k_haps)] += 1
                        num_sites_combined = num_sites_combined + inv_w.num_sites
                        if k_dist:
                            max_k = max(k_dist)
                            num_weighted = float(
                                sum(count * (k / max_k) for k, count in k_dist.items())
                            )
                        # Recompute avg_theta with the new denominator. `raw_theta` from
                        # `calc_watterson_theta` may already be a numeric value or "NA".
                        if num_sites_combined > 0 and not isinstance(raw_theta, str):
                            avg_theta = float(raw_theta) / num_sites_combined
                pixy_results = PixyTempResult(
                    pixy_stat=PixyStat.WATTERSON_THETA,
                    population_1=pop,
                    population_2="NA",
                    chromosome=chromosome,
                    window_pos_1=window_pos_1,
                    window_pos_2=window_pos_2,
                    calculated_stat=avg_theta,
                    shared_sites_with_alleles=num_sites_combined,
                    total_differences=raw_theta,
                    total_comparisons=num_var_sites,
                    total_missing=num_weighted,
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
            # `put` exists on both Queue and the AutoProxy returned by manager.Queue(); the
            # BaseProxy stub doesn't advertise it, so narrow for mypy.
            q.put(temp_pixy_content)  # type: ignore[union-attr]
        elif args.n_cores == 1:
            outfile = open(temp_file, "a")
            outfile.write(temp_pixy_content + "\n")
            outfile.close()
