"""Compute summary statistics."""

import argparse
from itertools import combinations
from typing import Dict
from typing import List
from typing import Tuple
from typing import Union
from typing import cast

import numpy as np
from allel import AlleleCountsArray
from allel import GenotypeArray
from allel import GenotypeVector
from allel import SortedIndex
from allel import windowed_hudson_fst
from allel import windowed_weir_cockerham_fst
from numpy.typing import NDArray

from pixy.args_validation import PixyStat
from pixy.calc import calc_dxy
from pixy.calc import calc_fst
from pixy.calc import calc_fst_persite
from pixy.calc import calc_pi
from pixy.enums import FSTEstimator
from pixy.models import DxyResult
from pixy.models import FstResult
from pixy.models import PiResult
from pixy.models import PixyTempResult


def precompute_filtered_variant_array(
    args: argparse.Namespace,
    gt_array: GenotypeArray,
    pos_array: SortedIndex,
    callset_is_none: bool,
    window_size: int,
    popindices: Dict[str, NDArray[np.int64]],
    chromosome: str,
) -> Tuple[List[PixyTempResult], Union[GenotypeArray, None], Union[SortedIndex, None]]:
    """
    Pre-compute a filtered array of variants for FST calculation.

    Args:
        args: The pixy command-line arguments.
        gt_array: Filtered genotypes ingested from the input VCF.
        pos_array: Positions of each genotype in `gt_array`.
        callset_is_none: True if the input VCF contained no variants in the region of interest.
        window_size: The size (in bp) of the window over which to calculate pixy statistics.
        popindices: Indices of the population-specific samples within the input VCF.
        chromosome: The chromosome of interest.

    Returns:
        A tuple (per_site_fst_results, gt_array_fst, pos_array_fst), where `gt_array_fst` and
        `pos_array_fst` contain the subset of variants in `gt_array` and `pos_array` which are
        biallelic. These objects are `None` if the input VCF contained no variants in the region of
        interest (i.e. `callset_is_none` is True or `gt_array` is empty), or when a `populations`
        file was not specified in the command-line arguments.

        **Per-site FST results are only computed by this function when `window_size` is set to 1.**
        When `window_size=1`, `per_site_fst_results` is a list of per-site FST values for all of the
        filtered sites in `gt_array_fst`. Otherwise, it is empty.
    """
    gt_array_fst: Union[GenotypeArray, None]
    pos_array_fst: Union[SortedIndex, None]

    if (not callset_is_none) and (args.populations is not None) and (len(gt_array) != 0):
        # compute allel freqs
        allele_counts: AlleleCountsArray = gt_array.count_alleles()
        allele_freqs: NDArray[np.float64] = allele_counts.to_frequencies()

        # remove invariant/polyallelic sites
        variants_array: List[bool] = [len(x) == 2 and x[0] < 1 for x in allele_freqs]

        # filter gt and position arrays for biallelic variable sites
        # NB: we are masking `gt_array` and `pos_array` with the filtered variants. Indexing either
        # `GenotypeArray` or `SortedIndex` with a boolean array returns a subset of the initial
        # object, but scikit-allel doesn't type this appropriately, so the cast is needed.
        gt_array_fst = cast(GenotypeArray, gt_array[variants_array])
        pos_array_fst = cast(SortedIndex, pos_array[variants_array])

    else:
        gt_array_fst = None
        pos_array_fst = None

    # if obtaining per-site estimates,
    # compute the FST values for the whole chunk
    # instead of looping over subwindows (below)
    pixy_results: List[PixyTempResult] = []

    # The original code only computed the per-site FST if the window size were 1
    # TODO extract and conslidate the FST-calculations
    # https://github.com/fulcrumgenomics/pixy-dev/issues/47
    if window_size != 1:
        return pixy_results, gt_array_fst, pos_array_fst

    # determine all the possible population pairings
    pop_names: List[str] = list(popindices.keys())
    fst_pop_list: List[Tuple[str, str]] = list(combinations(pop_names, 2))

    # for each pair, compute fst using the filtered gt_array
    for pop_pair in fst_pop_list:
        # the indices for the individuals in each population
        fst_pop_indicies: List[List[int]] = [
            popindices[pop_pair[0]].tolist(),
            popindices[pop_pair[1]].tolist(),
        ]

        # compute FST
        # windowed_weir_cockerham_fst seems to generate (spurious?) warnings about div/0, so
        # suppressing warnings (this assumes that the scikit-allel function is working as
        # intended)
        # TODO: verify these are indeed spurious
        # https://github.com/fulcrumgenomics/pixy-dev/issues/48
        np.seterr(divide="ignore", invalid="ignore")

        # if the genotype matrix is not empty, compute FST
        # otherwise return NA
        if not callset_is_none and gt_array_fst is not None and len(gt_array_fst) > 0:
            assert isinstance(pos_array_fst, SortedIndex), (
                "if gt_array_fst is not None, pos_array_fst should be not None as well"
            )

            per_site_fsts: NDArray[np.float64] = calc_fst_persite(
                gt_array_fst, fst_pop_indicies, args.fst_type
            )
            assert gt_array_fst.shape[0] == pos_array_fst.shape[0] == per_site_fsts.shape[0], (
                "the genotype, position, and FST arrays should have the same length"
            )

            snps: int = 1

            for i, fst in enumerate(per_site_fsts):
                # append trailing NAs so that pi/dxy/fst have the same # of columns
                pixy_result: PixyTempResult = PixyTempResult(
                    pixy_stat=PixyStat.FST,
                    population_1=pop_pair[0],
                    population_2=pop_pair[1],
                    chromosome=chromosome,
                    window_pos_1=pos_array_fst[i],
                    window_pos_2=pos_array_fst[i],
                    calculated_stat=fst,
                    shared_sites_with_alleles=snps,
                    total_differences="NA",
                    total_missing="NA",
                    total_comparisons="NA",
                )

                pixy_results.append(pixy_result)

    return pixy_results, gt_array_fst, pos_array_fst


def compute_summary_pi(
    popnames: NDArray[np.string_],
    window_is_empty: bool,
    gt_region: Union[GenotypeVector, None],
    popindices: Dict[str, NDArray[np.int64]],
    chromosome: str,
    window_pos_1: int,
    window_pos_2: int,
) -> List[PixyTempResult]:
    """Compute pi for all populations in the specified window."""
    pixy_results: List[PixyTempResult] = []

    for pop in popnames:
        pi_result: PiResult
        no_sites: int
        # if the window has no sites in the VCF, assign all NAs,
        # otherwise calculate pi
        if window_is_empty:
            pi_result = PiResult.empty()
            no_sites = 0
        else:
            # NB: in the original implementation this branch appears to be unreachable when
            # `gt_region` is None, the value error is included as a type narrower
            if gt_region is None:
                raise ValueError("gt_region is None")

            # subset the window for the individuals in each population
            gt_pop = gt_region.take(popindices[pop], axis=1)

            # if the population specific window for this region is empty, report it as such
            if len(gt_pop) == 0:
                pi_result = PiResult.empty()
                no_sites = 0

            # otherwise compute pi as normal
            else:
                # number of sites genotyped in the population
                # not directly used in the calculation
                no_sites = np.count_nonzero(np.sum(gt_pop.count_alleles(max_allele=1), 1))
                pi_result = calc_pi(gt_pop)

        # create a `PixyTempResult` composed of pi results to write to file
        pixy_result: PixyTempResult = PixyTempResult(
            pixy_stat=PixyStat.PI,
            population_1=pop,
            population_2="NA",
            chromosome=chromosome,
            window_pos_1=window_pos_1,
            window_pos_2=window_pos_2,
            calculated_stat=pi_result.avg_pi,
            shared_sites_with_alleles=no_sites,
            total_differences=pi_result.total_diffs,
            total_comparisons=pi_result.total_comps,
            total_missing=pi_result.total_missing,
        )

        # append the result to list of `PixyTempResult` objects
        pixy_results.append(pixy_result)

    return pixy_results


def compute_summary_dxy(
    popnames: NDArray[np.string_],
    window_is_empty: bool,
    gt_region: Union[GenotypeVector, None],
    popindices: Dict[str, NDArray[np.int64]],
    chromosome: str,
    window_pos_1: int,
    window_pos_2: int,
) -> List[PixyTempResult]:
    """Compute dxy for all pairwise combinations of populations in the specified window."""
    pixy_results: List[PixyTempResult] = []

    dxy_pop_list = list(combinations(popnames, 2))

    # iterate over all population pairs and compute dxy
    for pop_pair in dxy_pop_list:
        pop1 = pop_pair[0]
        pop2 = pop_pair[1]

        dxy_result: DxyResult
        no_sites: int
        if window_is_empty:
            dxy_result = DxyResult.empty()
            no_sites = 0

        else:
            # NB: in the original implementation this branch appears to be unreachable when
            # `gt_region` is None, the value error is included as a type narrower
            if gt_region is None:
                raise ValueError("gt_region is None")

            # use the pop_gts dictionary to keep track of this region's GTs for each population
            pop_gts = {}
            for name in pop_pair:
                gt_pop = gt_region.take(popindices[name], axis=1)
                pop_gts[name] = gt_pop

            pop1_gt_region = pop_gts[pop1]
            pop2_gt_region = pop_gts[pop2]

            # if either of the two population-specific windows for this region are empty,
            # report it missing
            if len(pop1_gt_region) == 0 or len(pop2_gt_region) == 0:
                dxy_result = DxyResult.empty()
                no_sites = 0

            # otherwise compute dxy as normal
            else:
                # for number of sites (not used in calculation), report the
                # number of sites that have at least one genotype in BOTH populations
                pop1_sites = np.sum(pop1_gt_region.count_alleles(max_allele=1), 1) > 0
                pop2_sites = np.sum(pop2_gt_region.count_alleles(max_allele=1), 1) > 0
                no_sites = np.sum(np.logical_and(pop1_sites, pop2_sites))
                dxy_result = calc_dxy(pop1_gt_array=pop1_gt_region, pop2_gt_array=pop2_gt_region)

        # create a string of for the dxy results
        pixy_result: PixyTempResult = PixyTempResult(
            pixy_stat=PixyStat.DXY,
            population_1=pop1,
            population_2=pop2,
            chromosome=chromosome,
            window_pos_1=window_pos_1,
            window_pos_2=window_pos_2,
            calculated_stat=dxy_result.avg_dxy,
            shared_sites_with_alleles=no_sites,
            total_differences=dxy_result.total_diffs,
            total_comparisons=dxy_result.total_comps,
            total_missing=dxy_result.total_missing,
        )

        # append the result to list of `PixyTempResult` objects
        pixy_results.append(pixy_result)

    return pixy_results


def compute_summary_fst(
    fst_type: FSTEstimator,
    gt_array_fst: GenotypeArray,
    pos_array_fst: Union[SortedIndex, None],
    window_is_empty: bool,
    callset_is_none: bool,
    aggregate: bool,
    popindices: Dict[str, NDArray[np.int64]],
    chromosome: str,
    window_pos_1: int,
    window_pos_2: int,
) -> List[PixyTempResult]:
    """Compute fst for all filtered sites."""
    pixy_results: List[PixyTempResult] = []

    # If there are no valid sites, exit early
    if pos_array_fst is None:
        return pixy_results

    # If there are no valid sites in the current window, exit early
    if not np.logical_and(pos_array_fst >= window_pos_1, pos_array_fst <= window_pos_2).any():
        return pixy_results

    # if there are valid sites, determine all the possible population pairings
    pop_names = list(popindices.keys())
    fst_pop_list = list(combinations(pop_names, 2))

    # for each pair, compute fst using the filtered gt_array
    for pop_pair in fst_pop_list:
        # the indices for the individuals in each population
        fst_pop_indicies = [
            popindices[pop_pair[0]].tolist(),
            popindices[pop_pair[1]].tolist(),
        ]

        # compute FST
        # windowed_weir_cockerham_fst seems to generate (spurious?) warnings about div/0, so
        # suppressing warnings (this assumes that the scikit-allel function is working as
        # intended)
        # TODO: verify these are indeed spurious
        # https://github.com/fulcrumgenomics/pixy-dev/issues/48
        np.seterr(divide="ignore", invalid="ignore")

        # if the genotype matrix is not empty, compute FST
        # other wise return NA
        gt_matrix_is_empty: bool = not (
            not callset_is_none
            and gt_array_fst is not None
            and len(gt_array_fst) > 0
            and not window_is_empty
        )

        fst_results: List[PixyTempResult]
        if aggregate:
            fst_results = _compute_aggregate_fst_for_pair(
                gt_matrix_is_empty=gt_matrix_is_empty,
                fst_type=fst_type,
                gt_array_fst=gt_array_fst,
                fst_pop_indicies=fst_pop_indicies,
                window_pos_1=window_pos_1,
                window_pos_2=window_pos_2,
                pop_pair=pop_pair,
                chromosome=chromosome,
            )
        else:
            fst_results = _compute_individual_fst_for_pair(
                gt_matrix_is_empty=gt_matrix_is_empty,
                fst_type=fst_type,
                pos_array_fst=pos_array_fst,
                gt_array_fst=gt_array_fst,
                fst_pop_indicies=fst_pop_indicies,
                window_pos_1=window_pos_1,
                window_pos_2=window_pos_2,
                pop_pair=pop_pair,
                chromosome=chromosome,
            )

        pixy_results.extend(fst_results)

    return pixy_results


def _compute_aggregate_fst_for_pair(
    gt_matrix_is_empty: bool,
    fst_type: FSTEstimator,
    gt_array_fst: GenotypeArray,
    fst_pop_indicies: List[List[int]],
    window_pos_1: int,
    window_pos_2: int,
    pop_pair: Tuple[str, str],
    chromosome: str,
) -> List[PixyTempResult]:
    """Compute aggregate FST for a pair of populations."""
    result: FstResult
    if gt_matrix_is_empty:
        result = FstResult.empty()
    else:
        result = calc_fst(gt_array_fst, fst_pop_indicies, fst_type)

    pixy_result = PixyTempResult(
        pixy_stat=PixyStat.FST,
        population_1=pop_pair[0],
        population_2=pop_pair[1],
        chromosome=chromosome,
        window_pos_1=window_pos_1,
        window_pos_2=window_pos_2,
        calculated_stat=result.fst,
        shared_sites_with_alleles=result.n_sites,
        total_differences=result.a,
        total_comparisons=result.b,
        total_missing=result.c,
    )

    pixy_results = [pixy_result]

    return pixy_results


def _compute_individual_fst_for_pair(
    gt_matrix_is_empty: bool,
    fst_type: FSTEstimator,
    pos_array_fst: SortedIndex,
    gt_array_fst: GenotypeArray,
    fst_pop_indicies: List[List[int]],
    window_pos_1: int,
    window_pos_2: int,
    pop_pair: Tuple[str, str],
    chromosome: str,
) -> List[PixyTempResult]:
    """Compte individual FST for a pair of populations."""
    if gt_matrix_is_empty:
        pixy_result = PixyTempResult(
            pixy_stat=PixyStat.FST,
            population_1=pop_pair[0],
            population_2=pop_pair[1],
            chromosome=chromosome,
            window_pos_1=window_pos_1,
            window_pos_2=window_pos_2,
            calculated_stat="NA",
            shared_sites_with_alleles=0,
            total_differences="NA",
            total_comparisons="NA",
            total_missing="NA",
        )

        return [pixy_result]

    pixy_results: List[PixyTempResult] = []

    # compute an ad-hoc window size
    fst_window_size = window_pos_2 - window_pos_1

    if fst_type is FSTEstimator.WC:
        fst, window_positions, n_snps = windowed_weir_cockerham_fst(
            pos_array_fst,
            gt_array_fst,
            subpops=fst_pop_indicies,
            size=fst_window_size,
            start=window_pos_1,
            stop=window_pos_2,
        )

    elif fst_type is FSTEstimator.HUDSON:
        ac1 = gt_array_fst.count_alleles(subpop=fst_pop_indicies[0])
        ac2 = gt_array_fst.count_alleles(subpop=fst_pop_indicies[1])
        fst, window_positions, n_snps = windowed_hudson_fst(
            pos_array_fst,
            ac1,
            ac2,
            size=fst_window_size,
            start=window_pos_1,
            stop=window_pos_2,
        )

    else:
        raise ValueError("unreachable")

    for fst_stat, wind, snps in zip(fst, window_positions, n_snps):
        # append trailing NAs so that pi/dxy/fst have the same # of columns
        pixy_result = PixyTempResult(
            pixy_stat=PixyStat.FST,
            population_1=pop_pair[0],
            population_2=pop_pair[1],
            chromosome=chromosome,
            window_pos_1=wind[0],
            window_pos_2=wind[1],
            calculated_stat=fst_stat,
            shared_sites_with_alleles=snps,
            total_differences="NA",
            total_comparisons="NA",
            total_missing="NA",
        )

        pixy_results.append(pixy_result)

    return pixy_results
