from collections import Counter
from functools import lru_cache
from typing import Any
from typing import Counter as CounterType
from typing import List
from typing import Mapping
from typing import Tuple
from typing import Union

import allel
import numpy as np
from allel import AlleleCountsArray
from allel import GenotypeArray
from numpy.typing import NDArray
from typing_extensions import TypeAlias

from pixy.enums import FSTEstimator
from pixy.models import NA
from pixy.models import DxyResult
from pixy.models import FstResult
from pixy.models import PiResult
from pixy.models import TajimaDResult
from pixy.models import WattersonThetaResult

# vectorized functions for calculating pi and dxy
# these are reimplementations of the original functions

# adding 2 `TypeAlias`s for additional type safety and defensiveness
VariantCount: TypeAlias = int
SiteCount: TypeAlias = int
ObservedAlleleCount: TypeAlias = int


@lru_cache(maxsize=1024)
def _harmonic_sum(n: int) -> np.float64:
    """Return sum_{k=1}^{n-1} 1/k. Returns np.float64(0.0) when n <= 1.

    Used by Watterson's theta and Tajima's D, which evaluate this for each variant-site class
    keyed by the per-site observed-allele count. The result depends only on `n`, so caching by
    `n` removes a quadratic-in-window-size cost when many sites share the same observed count.

    The return dtype is `np.float64` (not Python `float`) so that callers can rely on numpy's
    `inf` semantics for the n == 1 case under `np.errstate(divide="ignore")` — Python `float`
    would raise `ZeroDivisionError` instead, which is a behavior change from the previous
    `np.sum(1 / np.arange(1, n))` formulation.
    """
    if n <= 1:
        return np.float64(0.0)
    return np.float64(np.sum(1.0 / np.arange(1, n)))


@lru_cache(maxsize=1024)
def _harmonic_sum_sq(n: int) -> np.float64:
    """Return sum_{k=1}^{n-1} 1/k**2. Returns np.float64(0.0) when n <= 1."""
    if n <= 1:
        return np.float64(0.0)
    return np.float64(np.sum(1.0 / (np.arange(1, n) ** 2)))


@lru_cache(maxsize=1024)
def _tajima_constants(n: int) -> Tuple[float, float]:
    """Return (e1, e2) — the Tajima 1989 coefficients used in the stdev calculation.

    Returns (0.0, 0.0) for `n < 2` (no variance contribution possible).
    """
    if n < 2:
        return 0.0, 0.0
    a1 = _harmonic_sum(n)
    a2 = _harmonic_sum_sq(n)
    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n * n + n + 3) / (9 * n * (n - 1))
    c1 = b1 - (1 / a1)
    c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / (a1**2))
    e1 = c1 / a1
    e2 = c2 / (a1**2 + a2)
    return e1, e2


def _n_haps(gt_array: GenotypeArray) -> int:
    """Return the number of haploid samples represented by a (Haplotype|Genotype)Array.

    `HaplotypeArray` exposes this directly as `n_haplotypes`; for the general case the
    haploid sample count is `n_samples * ploidy`. This branch appeared inline in several
    callers below; factoring it out removes the duplication.
    """
    if isinstance(gt_array, allel.HaplotypeArray):
        return int(gt_array.n_haplotypes)
    return int(gt_array.n_samples * gt_array.ploidy)


def count_diff_comp_missing(row: NDArray[Any], n_haps: int) -> Tuple[int, int, int]:
    """
    Helper function for the calculation of pi.

    For the given site (row of the count table), count the number of differences, number of
    comparisons, and number of missing. The function uses the number of haploid samples (n_haps) to
    compute the number of expected genotype comparisons and determine the count of missing.

    Args:
        row: counts of each allele at a given site
        n_haps: number of haploid samples in the population

    Returns:
        A tuple `(diffs, comps, missing)`, where `diffs` is the number of differences within the
        population, `comps` is the number of comparisons made within the population, and `missing`
        is the difference between the actual number of comparisons and the total possible (based on
        the number of haploid samples).
    """
    # Algebraic simplification of the previous implementation:
    #   diffs = sum_{i<j} row[i] * row[j]  = (n_gts**2 - sum(row**2)) / 2
    #   comps = C(n_gts, 2)                 = n_gts * (n_gts - 1) / 2
    # Both forms produce identical integer results (n_gts and the row entries are nonneg ints,
    # and n_gts**2 - sum(row**2) is even by construction). The old form used an O(n_alleles**2)
    # Python loop plus scipy.special.comb; this version uses pure arithmetic.
    n_gts = int(np.sum(row))
    n_possible_comps = n_haps * (n_haps - 1) // 2

    if n_gts == 0:
        return 0, 0, n_possible_comps

    comps = n_gts * (n_gts - 1) // 2
    diffs = (n_gts * n_gts - int(np.sum(np.asarray(row, dtype=np.int64) ** 2))) // 2
    missing = n_possible_comps - comps
    return diffs, comps, missing


# Largest `n_haps` whose square still fits in a signed int32. Beyond this we must promote
# `allele_counts` to int64 in the vectorized sum-of-squares; below it we can stay in the
# native int32 dtype that scikit-allel's `count_alleles()` returns, halving the intermediate
# memory footprint. Threshold = floor(sqrt(2**31 - 1)).
_INT32_SAFE_NHAPS_MAX = 46340


def _count_diff_comp_missing_vectorized(
    allele_counts: NDArray[Any], n_haps: int
) -> Tuple[NDArray[np.int64], NDArray[np.int64], NDArray[np.int64]]:
    """Vectorized form of `count_diff_comp_missing` over all sites at once.

    Inputs:
        allele_counts: (n_sites, n_alleles) int array (typically int32 from scikit-allel).
        n_haps: number of haploid samples in the population.

    Returns:
        Three (n_sites,) int64 arrays: (diffs, comps, missing).
    """
    # Keep `allele_counts` in its native int32 dtype unless `n_haps` is large enough that
    # `n_gts**2` could overflow int32 (>= 46341). The bulk of the memory in the intermediate
    # `ac * ac` lives here — staying in int32 halves it for typical pixy datasets.
    ac = np.asarray(allele_counts)
    if n_haps > _INT32_SAFE_NHAPS_MAX or ac.dtype.itemsize < 4:
        ac = ac.astype(np.int64)
    # `sum(axis=1)` on int32 returns int64 by default on 64-bit platforms, which gives us the
    # headroom we need for `n_gts * n_gts` even when `ac` itself is still int32.
    n_gts = ac.sum(axis=1, dtype=np.int64)
    n_possible_comps = n_haps * (n_haps - 1) // 2
    comps = n_gts * (n_gts - 1) // 2
    # sum_{i<j} a_i*a_j = (n_gts**2 - sum(a_i**2)) / 2 — see the scalar form above. We use
    # np.einsum so numpy never materializes the full (n_sites, n_alleles) product matrix;
    # it computes the row sums directly.
    sq_sum = np.einsum("ij,ij->i", ac, ac, dtype=np.int64)
    diffs = (n_gts * n_gts - sq_sum) // 2
    missing = n_possible_comps - comps
    return diffs, comps, missing


def calc_pi(gt_array: GenotypeArray) -> PiResult:
    """
    Given a filtered genotype matrix, calculate `pi`.

    This function implements vectorized calculation of pi from a pre-filtered scikit-allel genotype
    matrix. This function does not support filtering of the input by population - it simply
    calculates pi on all of the provided samples.

    Args:
        gt_array: a GenotypeArray representing the calls of each variant at each filtered site in a
            given population. This array must be pre-filtered to the population of interest.

    Returns:
        The average `pi` and total difference, comparison, and missing counts over all sites in the
        input array.
    """
    # counts of each of the two alleles at each site
    allele_counts: AlleleCountsArray = gt_array.count_alleles()

    n_haps = _n_haps(gt_array)

    # Vectorized over all sites at once; previously this used np.apply_along_axis with a Python
    # callback per site, which is the same as a for-loop in numpy clothing.
    diffs, comps, missing = _count_diff_comp_missing_vectorized(np.asarray(allele_counts), n_haps)

    total_diffs = int(diffs.sum())
    total_comps = int(comps.sum())
    total_missing = int(missing.sum())

    # if there are valid data (comparisons between genotypes) at the site, compute average pi
    # otherwise return NA
    avg_pi: Union[float, NA] = total_diffs / total_comps if total_comps > 0 else "NA"

    return PiResult(
        avg_pi=avg_pi,
        total_diffs=total_diffs,
        total_comps=total_comps,
        total_missing=total_missing,
    )


# function for vectorized calculation of dxy from a pre-filtered scikit-allel genotype matrix
def calc_dxy(pop1_gt_array: GenotypeArray, pop2_gt_array: GenotypeArray) -> DxyResult:
    """
    Given a filtered genotype matrix, calculate `dxy`.

    Args:
        pop1_gt_array: the GenotypeArray representing population-specific allele counts
        pop2_gt_array: the GenotypeArray representing population-specific allele counts


    Returns:
        avg_dxy: proportion of total differences across total comparisons. "NA" if no valid data.
        total_diffs: sum of the number of differences between the populations
        total_comps: sum of the number of comparisons between the populations
        total_missing: sum of the number of missing between the populations
    """
    if pop1_gt_array.n_variants != pop2_gt_array.n_variants:
        raise ValueError("Input genotype matrices must have the same number of variants")

    n_sites: int = pop1_gt_array.n_variants

    # the counts of each of the two alleles in each population at each site
    pop1_allele_counts: AlleleCountsArray = pop1_gt_array.count_alleles()
    pop2_allele_counts: AlleleCountsArray = pop2_gt_array.count_alleles()

    pop1_n_haps = _n_haps(pop1_gt_array)
    pop2_n_haps = _n_haps(pop2_gt_array)

    # Per-site differences (= pairwise comparisons between pops that differ in allele).
    # Identity used:
    #     sum_{i != j} a_i * b_j  =  (sum a) * (sum b)  -  sum_i a_i * b_i
    # This replaces an O(n_alleles**2) Python loop with a handful of vectorized numpy ops,
    # and lets us reuse `persite_comps` (== pop1_sum * pop2_sum) instead of recomputing it.
    #
    # Memory: we stay in the native int32 dtype that scikit-allel's `count_alleles()` returns
    # unless `pop1_n_haps * pop2_n_haps` would overflow int32. Sums and products are computed
    # into explicit int64 accumulators (numpy promotes during `sum`/`einsum`), so the only
    # large intermediate that scales with `n_sites * n_alleles` stays in int32.
    pop1_arr = np.asarray(pop1_allele_counts)
    pop2_arr = np.asarray(pop2_allele_counts)
    max_product = pop1_n_haps * pop2_n_haps
    if max_product > 2**31 - 1 or pop1_arr.dtype.itemsize < 4 or pop2_arr.dtype.itemsize < 4:
        pop1_arr = pop1_arr.astype(np.int64)
        pop2_arr = pop2_arr.astype(np.int64)
    pop1_sum = pop1_arr.sum(axis=1, dtype=np.int64)
    pop2_sum = pop2_arr.sum(axis=1, dtype=np.int64)
    persite_comps: NDArray = pop1_sum * pop2_sum
    # einsum 'ij,ij->i' is the row-wise dot product without materializing the elementwise
    # product matrix `pop1_arr * pop2_arr` (which would be (n_sites, n_alleles) of int).
    persite_diffs: NDArray = persite_comps - np.einsum(
        "ij,ij->i", pop1_arr, pop2_arr, dtype=np.int64
    )
    assert persite_comps.shape == (n_sites,)

    total_diffs: int = int(persite_diffs.sum())
    total_comps: int = int(persite_comps.sum())

    # the total count of possible pairwise comparisons at all sites
    total_possible: int = (pop1_n_haps * pop2_n_haps) * n_sites

    # the amount of missing is possible comps - actual ('total') comps
    total_missing: int = total_possible - total_comps

    # if there are valid data (comparisons between genotypes) at the site, compute average dxy
    # otherwise return NA
    avg_dxy: Union[float, NA] = total_diffs / total_comps if total_comps > 0 else "NA"

    return DxyResult(
        avg_dxy=avg_dxy,
        total_diffs=total_diffs,
        total_comps=total_comps,
        total_missing=total_missing,
    )


# function for obtaining fst AND variance components via scikit allel function
# (need variance components for proper aggregation)
# for single sites, this is the final FST calculation
# in aggregation mode, we just want a,b,c and n_sites for aggregating and fst
def calc_fst(
    gt_array_fst: GenotypeArray, fst_pop_indicies: List[List[int]], fst_type: FSTEstimator
) -> FstResult:
    # TODO: update the return type here after refactoring (2 -> 1 return statements)
    """
    Calculates FST according to either Weir and Cockerham (1984) or Hudson (1992).

    FST is a measure of the total genetic variance within a subpopulation relative to total genetic
    variance.

    Args:
        gt_array_fst: allele counts to use for computation of variance
        fst_pop_indicies: sample indices for each subpopulation
        fst_type: one of either WC or Hudson, corresponding to the method of calculation

    Returns:
        If `fst_type` is `wc`, the following will be returned:
            fst: "NA" if no valid data
            a: variance between populations
            b: variance between individuals within populations
            c: variance between gametes within individuals
            n_sites: the number of sites over which variance was measured

        The shape of `a`, `b`, and `c` correspond to [n_variants * n_alleles].

        If `fst_type` is `hudson`, the following will be returned:
            fst: "NA" if no valid data.
            num: divergence between the two populations minus average of diversity within
              each population
            den: divergence between the two populations.
            c: a placeholder of 0 to maintain the shape of the return Tuple
            n_sites: the number of sites over which variance was measured

    """
    # compute basic (multisite) FST via scikit allel
    fst: Union[float, NA]
    result: FstResult

    n_sites: int = len(gt_array_fst)

    # WC 84
    if fst_type is FSTEstimator.WC:
        if gt_array_fst.ploidy != 2:
            raise NotImplementedError(
                "`pixy` does not currently support calculation "
                "of Weir-Cockerham FST in non-diploid genomes"
            )
        a: NDArray
        b: NDArray
        c: NDArray
        a, b, c = allel.weir_cockerham_fst(gt_array_fst, subpops=fst_pop_indicies)

        # compute variance component sums
        a_sum: float = np.nansum(a)
        b_sum: float = np.nansum(b)
        c_sum: float = np.nansum(c)

        # compute fst
        if (a_sum + b_sum + c_sum) > 0:
            fst = a_sum / (a_sum + b_sum + c_sum)
        else:
            fst = "NA"

        result = FstResult(fst=fst, a=a_sum, b=b_sum, c=c_sum, n_sites=n_sites)

    # Hudson 92
    if fst_type is FSTEstimator.HUDSON:
        # following scikit allel docs
        # allel counts for each population
        ac1: AlleleCountsArray = gt_array_fst.count_alleles(subpop=fst_pop_indicies[0])
        ac2: AlleleCountsArray = gt_array_fst.count_alleles(subpop=fst_pop_indicies[1])

        # hudson fst has two components (numerator & denominator)
        num: NDArray
        den: NDArray
        num, den = allel.hudson_fst(ac1, ac2)

        # compute variance component sums
        num_sum: float = np.nansum(num)
        den_sum: float = np.nansum(den)

        # compute fst
        if (num_sum + den_sum) > 0:
            fst = num_sum / den_sum
        else:
            fst = "NA"

        # same abc format as WC84, where 'a' is the numerator and
        # 'b' is the demoninator, and 'c' is a zero placeholder
        result = FstResult(fst=fst, a=num_sum, b=den_sum, c=0, n_sites=n_sites)

    return result


# simplified version of above to handle the case
# of per-site estimates of FST over whole chunks


def calc_fst_persite(
    gt_array_fst: GenotypeArray,
    fst_pop_indicies: List[List[int]],
    fst_type: str,
) -> Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """
    Calculates site-specific FST according to Weir and Cockerham (1984) or Hudson (1992).

    Args:
        gt_array_fst: allele counts to use for computation of variance
        fst_pop_indicies: sample indices for each subpopulation
        fst_type: one of either WC or Hudson, corresponding to the method of calculation

    Returns:
        A tuple of arrays containing site-specific FST and the summed components used to compute it.
        For Weir-Cockerham FST these components are ``a``, ``b``, and ``c``. For Hudson FST they
        are numerator, denominator, and a zero-valued placeholder.

    """
    # compute basic (multisite) FST via scikit allel
    fst: NDArray[np.float64]
    a_sum: NDArray[np.float64]
    b_sum: NDArray[np.float64]
    c_sum: NDArray[np.float64]

    # WC 84
    if fst_type == "wc":
        a, b, c = allel.weir_cockerham_fst(gt_array_fst, subpops=fst_pop_indicies)

        a_raw = np.sum(a, axis=1)
        b_raw = np.sum(b, axis=1)
        c_raw = np.sum(c, axis=1)
        a_sum = np.nansum(a, axis=1)
        b_sum = np.nansum(b, axis=1)
        c_sum = np.nansum(c, axis=1)
        fst = a_raw / (a_raw + b_raw + c_raw)

    # Hudson 92
    elif fst_type == "hudson":
        # following scikit allel docs
        # allel counts for each population
        ac1 = gt_array_fst.count_alleles(subpop=fst_pop_indicies[0])
        ac2 = gt_array_fst.count_alleles(subpop=fst_pop_indicies[1])

        # hudson fst has two components (numerator & denominator)
        num, den = allel.hudson_fst(ac1, ac2)

        a_sum = np.nan_to_num(num, nan=0.0)
        b_sum = np.nan_to_num(den, nan=0.0)
        c_sum = np.zeros_like(num)
        fst = num / den

    else:
        raise ValueError("fst_type must be either 'wc' or 'hudson'")

    return fst, a_sum, b_sum, c_sum


def calc_watterson_theta(gt_array: GenotypeArray) -> WattersonThetaResult:
    """
    Calculates Watterson's Theta for a provided genotype array based on Watterson 1975.

    This equation is also given as Equation 3 in Bailey et al., 2025.

    `num_sites` represents the number of sites in the input array with an observed genotype.
    `num_var_sites` represents the number of sites in the input array with an observed alternate
    genotype (variant count != 0). The `weighted_site_count` represents the number of sites
    weighted by how many genotypes are missing in each site.

    Args:
        gt_array: the genotype array for which to calculate Watterson's Theta

    Returns:
        a WattersonThetaResult that encapsulates the number of sites with an observed genotype,
    the number of sites with an observe alternatve genotype, the averaged Watterson's Theta,
    raw Watterson's Theta, and weighted site count
    """
    # alleles are counted, variant alleles are extracted
    # number of sites is count of sites with more than zero alleles

    # counts of each of the two alleles at each site
    # with max_allele=1    -> this returns an (n x 2) array
    # with max_allele=None -> this returns an (n x a) array, where a is the # of alternate alleles
    allele_counts: AlleleCountsArray = gt_array.count_alleles()

    # counts of only variant sites by excluding sites with variant count 0
    # [ ref alt ]  => ref + alt == # haploid samples (assuming none are missing)
    variant_counts: AlleleCountsArray = allele_counts[allele_counts[:, 1:].sum(axis=1) != 0]

    num_sites = np.count_nonzero(np.sum(allele_counts, 1))
    num_var_sites = np.count_nonzero(np.sum(variant_counts, 1))
    # for variant sites only use Counter to generate dictionary
    # where the key is the number of genotypes and value is number of sites with that many genotypes
    variant_sites_counter: CounterType[VariantCount] = Counter(variant_counts.sum(axis=1))
    all_sites_counter: CounterType[SiteCount] = Counter(allele_counts.sum(axis=1))

    allele_freq_counts: NDArray[np.int64] = np.array(tuple(all_sites_counter.items()))

    # calculate Watterson's theta as sum of equations for differing numbers of genotypes
    # this is calculating Watterson's theta incorporating missing genotypes
    #
    # NB: when only a single (haploid) genotype is observed at a site (num_genotypes == 1),
    # `np.arange(1, 1)` is empty so `reciprocal_sum == 0`, and `site_count / 0` evaluates to
    # `inf`. This is the mathematically correct sentinel — Watterson's theta is undefined when
    # you can't observe variation across samples — and is exactly what the
    # `test_calc_watterson_theta_haploid_singleton` test asserts. `np.errstate` suppresses the
    # accompanying RuntimeWarning so it doesn't pollute pytest output for this known case.
    watterson_theta: float = 0.0
    with np.errstate(divide="ignore"):
        for num_genotypes, site_count in variant_sites_counter.items():
            # `_harmonic_sum(n)` returns 0.0 for n <= 1, preserving the documented "inf"
            # sentinel for the singleton-haploid case (see comment above).
            reciprocal_sum: float = _harmonic_sum(int(num_genotypes))
            watterson_theta += site_count / reciprocal_sum

    # Calculate an auxiliary effective-site count that reflects within-site missingness.
    # The Watterson's theta denominator is `num_sites`; this value is emitted as a diagnostic.

    weighted_sites: float
    if max(all_sites_counter) == 0:
        weighted_sites = 0.0
    else:
        # NB: this is a fractional weighted count — the inner division produces floats.
        weighted_sites = float(
            np.sum(
                np.multiply(
                    allele_freq_counts[:, 1:].sum(axis=1),
                    (allele_freq_counts[:, 0] / max(all_sites_counter)),
                )
            )
        )

    avg_theta: Union[float, NA]
    if num_sites == 0:
        avg_theta = "NA"
    else:
        avg_theta = watterson_theta / num_sites

    return WattersonThetaResult(
        num_sites=int(num_sites),
        num_var_sites=int(num_var_sites),
        avg_theta=avg_theta,
        raw_theta=watterson_theta,
        num_weighted_sites=weighted_sites,
    )


def calc_tajima_d_stdev(
    variant_gt_counts: Mapping[ObservedAlleleCount, SiteCount],
) -> float:
    """
    Calculates the standard deviation term used as Tajima's D denominator.

    The missing-data correction from Bailey et al. 2025 sums the denominator over variant-site
    classes with the same number of observed alleles. `variant_gt_counts` maps each observed allele
    count `n` to the number of segregating sites observed with that count.
    """
    d_stdev = 0.0
    for n, s in variant_gt_counts.items():
        if n < 2 or s <= 0:
            continue
        # (e1, e2) depend only on `n` and are cached so that repeated `n` values across many
        # sites in a window don't re-derive the same a1/a2/b1/b2/c1/c2/e1/e2 chain.
        e1, e2 = _tajima_constants(int(n))
        # np.sqrt (not math.sqrt) to preserve the prior behavior of returning nan on negative
        # variance terms rather than raising ValueError.
        d_stdev += float(np.sqrt((e1 * s) + (e2 * s * (s - 1))))

    return float(d_stdev)


def serialize_tajima_d_variant_counts(
    variant_gt_counts: Mapping[ObservedAlleleCount, SiteCount],
) -> str:
    """Serializes Tajima's D observed-allele-count classes for the temp file."""
    if not variant_gt_counts:
        return "NA"

    return ",".join(f"{n}:{s}" for n, s in sorted(variant_gt_counts.items()))


def deserialize_tajima_d_variant_counts(value: object) -> CounterType[ObservedAlleleCount]:
    """Deserializes Tajima's D observed-allele-count classes from the temp file."""
    if not isinstance(value, str) or value == "NA" or value == "":
        return Counter()

    counts: CounterType[ObservedAlleleCount] = Counter()
    for item in value.split(","):
        n, s = item.split(":")
        counts[int(n)] += int(s)

    return counts


def calc_tajima_d(gt_array: GenotypeArray) -> TajimaDResult:
    """
    Calculates Tajima's D over a provided genotype array based on Tajima 1989.

    This equation is also given as Equation 4 in Bailey et al., 2025. The number of genotyped sites
     is not directly used in the calculation but returned as it's required for the emitted results.

    Args:
        gt_array: the genotype array for which to calculate Tajima's D

    Returns:
        a TajimaDResult that encapsulates the calculated Tajima's D, number of sites with an
        observed genotype, the calculated raw values for pi and Watterson's theta, and the standard
        deviation (Tajima's D denominator)
    """
    # The number of total genotypes observed at a variant site.
    NumGenotypes: TypeAlias = int

    # counts of each of the two alleles at each site
    allele_counts: AlleleCountsArray = gt_array.count_alleles()
    num_sites = np.count_nonzero(np.sum(allele_counts, 1))

    # calculate mean pairwise differences and sum these together
    mpd: NDArray[np.floating] = allel.mean_pairwise_difference(ac=allele_counts, fill=0)
    raw_pi: float = np.sum(mpd)

    # counts of only variant sites by excluding sites with variant count 0
    variant_allele_counts: AlleleCountsArray = allele_counts[allele_counts[:, 1:].sum(axis=1) != 0]

    # for variant sites only, use Counter to generate dictionary
    # where the key is the number of genotypes and value is number of sites with that many genotypes
    variant_gt_counts: CounterType[NumGenotypes] = Counter(variant_allele_counts.sum(axis=1))

    # calculate watterson's theta as sum of equations for differing numbers of genotypes
    # this is calculating Watterson's theta incorporating missing genotypes
    #
    # NB: same singleton-haploid corner as in `calc_watterson_theta` — when `n == 1` the
    # denominator `a1` is 0 and `s / a1` evaluates to `inf`. The
    # `test_calc_tajima_d_haploid_singleton` test asserts exactly that; `np.errstate` suppresses
    # the accompanying RuntimeWarning.
    watterson_theta: float = 0.0
    with np.errstate(divide="ignore"):
        for n, s in variant_gt_counts.items():
            # See comment in calc_watterson_theta — `_harmonic_sum(1) == 0.0` reproduces the
            # documented inf sentinel for the singleton case.
            a1: float = _harmonic_sum(int(n))
            watterson_theta += s / a1

    d_stdev = calc_tajima_d_stdev(variant_gt_counts)

    tajima_d: Union[float, NA]
    if d_stdev > 0 and not any(np.isnan(x) for x in [raw_pi, watterson_theta, d_stdev]):
        tajima_d = (raw_pi - watterson_theta) / d_stdev
    else:
        tajima_d = "NA"

    # return Tajima's D calculation using raw pi and Watterson's theta calculations above
    # also return the raw pi calculation, raw Watterson's theta, and standard deviation of their
    # covariance individually
    # note that the "raw" values of pi and Watterson's theta are needed for Tajima's D, not the
    # ones incorporating sites
    return TajimaDResult(
        tajima_d=tajima_d,
        num_sites=int(num_sites),
        raw_pi=raw_pi,
        watterson_theta=watterson_theta,
        d_stdev=d_stdev,
        variant_gt_counts=dict(variant_gt_counts),
    )
