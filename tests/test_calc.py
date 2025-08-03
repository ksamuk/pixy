import allel
import numpy as np
import pytest
from allel import AlleleCountsArray
from allel import GenotypeArray
from allel import hudson_fst
from allel import mean_pairwise_difference
from allel import mean_pairwise_difference_between
from allel import watterson_theta
from allel import weir_cockerham_fst

from pixy.calc import calc_dxy
from pixy.calc import calc_fst
from pixy.calc import calc_pi
from pixy.calc import calc_tajima_d
from pixy.calc import calc_watterson_theta
from pixy.enums import FSTEstimator
from pixy.models import DxyResult
from pixy.models import FstResult
from pixy.models import PiResult
from pixy.models import TajimaDResult
from pixy.models import WattersonThetaResult

##################
# Helper functions
##################


def formulaic_wattersons_theta(gt_array: GenotypeArray) -> float:
    """
    Calculate Watterson's Theta per Equation 2 in Bailey et al., 2025.

    Args:
        gt_array: the genotype array over which to calculate Watterson's Theta

    Returns:
        the per-site (raw) theta calculation
    """
    # count the number of segregating sites (S)
    ac: AlleleCountsArray = gt_array.count_alleles()
    num_segregating_sites: int = ac.count_segregating()
    n: int = ac.sum(axis=1).max()
    # calculate a_n (the reciprocal sums from 1 to n-1)
    a_n = np.sum([1 / i for i in range(1, n)])

    # calculate Watterson's theta
    raw_theta: float = num_segregating_sites / a_n
    return raw_theta


def test_calc_pi() -> None:
    """Test that our Pi calculation matches the expected values based on the paper's formula."""
    array = GenotypeArray([
        [[0, 0], [0, 1], [1, 1]],  # site 1
        [[0, 1], [0, 1], [1, 1]],  # site 2
    ])

    # sanity
    assert array.n_samples == 3
    assert array.n_variants == 2
    assert array.ploidy == 2

    result: PiResult = calc_pi(array)

    # Section 1.1 of Korunes & Samuk
    # https://pmc.ncbi.nlm.nih.gov/articles/PMC8044049/
    expected_diffs = 9 + 8  # 9 diffs at site 1, 8 at site 2
    expected_comps = 15 * 2  # (6 choose 2) = 15 comparisons at each of two sites

    assert result.total_diffs == expected_diffs
    assert result.total_comps == expected_comps
    assert result.total_missing == 0  # all samples have calls at all sites
    assert result.avg_pi == pytest.approx(expected_diffs / expected_comps)


def test_calc_pi_single_locus() -> None:
    """Test that our Pi calculation matches the expected values based on the paper's formula."""
    array = GenotypeArray([
        [[0, 0], [0, 1], [1, 1]],  # site 1
    ])

    # sanity
    assert array.n_samples == 3
    assert array.n_variants == 1
    assert array.ploidy == 2

    result: PiResult = calc_pi(array)

    # Section 1.1 of Korunes & Samuk
    # https://pmc.ncbi.nlm.nih.gov/articles/PMC8044049/
    expected_diffs = 9  # 3 zeros * 3 ones
    expected_comps = 15  # (6 choose 2) = 15 comparisons

    assert result.total_diffs == expected_diffs
    assert result.total_comps == expected_comps
    assert result.total_missing == 0
    assert result.avg_pi == pytest.approx(expected_diffs / expected_comps)


def test_calc_dxy() -> None:
    """Test that our Dxy calculation matches the expected values based on the paper's formula."""
    pop1_gt_array = GenotypeArray([
        [[0, 0], [0, 1], [1, 1]],  # population 1/site 1
        [[0, 1], [0, 0], [1, 1]],  # population 1/site 2
    ])

    pop2_gt_array = GenotypeArray([
        [[0, 1], [0, 1], [1, 1]],  # population 2/site 1
        [[0, 1], [1, 1], [1, 1]],  # population 2/site 2
    ])

    # sanity
    assert pop1_gt_array.n_samples == pop2_gt_array.n_samples == 3
    assert pop1_gt_array.n_variants == pop2_gt_array.n_variants == 2
    assert pop1_gt_array.ploidy == pop2_gt_array.ploidy == 2  # diploid

    result: DxyResult = calc_dxy(pop1_gt_array=pop1_gt_array, pop2_gt_array=pop2_gt_array)

    # Section 1.1 of Korunes & Samuk
    # https://pmc.ncbi.nlm.nih.gov/articles/PMC8044049/
    expected_diffs_site_1 = 3 * 4 + 3 * 2  # allele_counts are [3,3] for pop1 and [2,4] for pop2
    expected_comps_site_1 = 6 * 6  # number of pairwise comparisons among alleles (6 per pop)

    expected_diffs_site_2 = 3 * 5 + 3 * 1  # allele_counts are [3,3] for pop1 and [1,5] for pop2
    expected_comps_site_2 = 6 * 6  # number of pairwise comparisons among alleles (6 per pop)

    expected_diffs = expected_diffs_site_1 + expected_diffs_site_2
    expected_comps = expected_comps_site_1 + expected_comps_site_2

    assert result.total_diffs == expected_diffs
    assert result.total_comps == expected_comps
    assert result.total_missing == 0
    assert result.avg_dxy == pytest.approx(expected_diffs / expected_comps)


def test_calc_dxy_single_locus() -> None:
    """Test that our Dxy calculation matches the expected values based on the paper's formula."""
    pop1_gt_array = GenotypeArray([
        [[0, 0], [0, 1], [1, 1]],  # population 1/site 1
    ])

    pop2_gt_array = GenotypeArray(
        [[[0, 1], [0, 1], [1, 1]]]  # population 2/site 1
    )

    # sanity
    assert pop1_gt_array.n_samples == pop2_gt_array.n_samples == 3
    assert pop1_gt_array.n_variants == pop2_gt_array.n_variants == 1
    assert pop1_gt_array.ploidy == pop2_gt_array.ploidy == 2  # diploid

    result: DxyResult = calc_dxy(pop1_gt_array=pop1_gt_array, pop2_gt_array=pop2_gt_array)

    # Section 1.1 of Korunes & Samuk
    # https://pmc.ncbi.nlm.nih.gov/articles/PMC8044049/
    expected_diffs = 3 * 4 + 3 * 2  # allele_counts are [3,3] for pop1 and [2,4] for pop2
    expected_comps = 6 * 6  # number of pairwise comparisons among alleles (6 per pop)

    assert result.total_diffs == expected_diffs
    assert result.total_comps == expected_comps
    assert result.total_missing == 0
    assert result.avg_dxy == pytest.approx(expected_diffs / expected_comps)


def test_calc_pi_arbitrary_ploidy() -> None:
    """
    Test that our Pi calculation is valid in the context of arbitrary ploidy.

    This test evaluates   at one site in a tetraploid context.
    Importantly, there are NO missing genotypes here.
    This test asserts the following:
      - the `pixy` calculation of pi matches the scikit-allel implementation
          (which is sufficient with no missing sites)
      - the pairwise comparison short-cut is valid for a given gt_array of arbitrary ploidy
      - the `pixy` calculation matches the publication
    """
    gt_array = GenotypeArray([
        [[0, 0, 0, 0], [0, 1, 0, 0], [1, 1, 0, 0]],
    ])  # tetraploid

    # sanity
    assert gt_array.n_samples == 3
    assert gt_array.n_variants == 1
    assert gt_array.ploidy == 4
    result: PiResult = calc_pi(gt_array)

    # Section 1.1 of Korunes & Samuk
    # https://pmc.ncbi.nlm.nih.gov/articles/PMC8044049/
    expected_diffs = 9 * 3  # total 9 zeros * 3 ones
    expected_comps = 66  # 12 choose 2
    assert result.total_diffs == expected_diffs
    assert result.total_comps == expected_comps
    assert result.total_missing == 0
    assert result.avg_pi == pytest.approx(expected_diffs / expected_comps)

    # assert that `pixy` matches `scikit-allel` in the context of no missing sites
    allele_counts: AlleleCountsArray = gt_array.count_alleles()  # [9,3]
    assert result.avg_pi == pytest.approx(
        mean_pairwise_difference(ac=allele_counts, an=np.sum(allele_counts, axis=1))
    )


def test_calc_dxy_arbitrary_ploidy() -> None:
    """
    Test that our Dxy calculation is valid in the context of arbitrary ploidy.

    Currently this test looks at one site in a tetraploid context only.
    In the future we can parameterize `GenotypeArray` inputs.
    Importantly, there are NO missing genotypes here.
    This test asserts the following:
      - the `pixy` calculation matches a brute-force approach for arbitrary ploidy
      - the pairwise comparison short-cut is valid for given gt_arrays of arbitrary ploidy
      - the `pixy` calculation  matches the publication
    """
    pop1_gt_array = GenotypeArray([
        [[0, 0, 0, 0], [0, 1, 0, 0], [1, 1, 0, 0]],
    ])  # population 1/site 1

    pop2_gt_array = GenotypeArray([
        [[0, 1, 0, 0], [0, 1, 0, 0], [1, 1, 1, 1]],  # population 2/site 1
    ])

    # sanity
    assert pop1_gt_array.n_samples == pop2_gt_array.n_samples == 3
    assert pop1_gt_array.n_variants == pop2_gt_array.n_variants == 1
    assert pop1_gt_array.ploidy == pop2_gt_array.ploidy == 4
    result: DxyResult = calc_dxy(pop1_gt_array=pop1_gt_array, pop2_gt_array=pop2_gt_array)
    expected_diffs = (9 * 6) + (3 * 6)
    expected_comps = 144  # 12 alleles in each population
    assert result.total_diffs == expected_diffs
    assert result.total_comps == expected_comps
    assert result.total_missing == 0
    assert result.avg_dxy == pytest.approx(expected_diffs / expected_comps)

    # scikit-allel provides support for dxy calculation
    # (but unlike pixy does not support missing sites)
    # double-check that we match their implementation
    allel_dxy = mean_pairwise_difference_between(
        ac1=pop1_gt_array.count_alleles(),
        ac2=pop2_gt_array.count_alleles(),
    )
    assert result.avg_dxy == pytest.approx(allel_dxy)


def test_calc_fst_hudson_arbitrary_ploidy() -> None:
    """
    Compare FST calculation (Hudson) to the `scikit-allel` implementation.

    The WC FST calculation for arbitrary ploidy is not currently supported within `scikit-allel`
    and thus is not tested here.

    """
    # scikit-allel's implementation accepts one array for each population
    pop1_gt_array = GenotypeArray([
        [[0, 0, 0, 0], [0, 1, 0, 0], [1, 1, 0, 0]],
        [[0, 1, 0, 0], [0, 1, 0, 0], [1, 1, 1, 1]],  # population 2/site 1
    ])  # population 1/site 1

    pop2_gt_array = GenotypeArray([
        [[0, 1, 0, 0], [0, 1, 0, 0], [1, 1, 1, 1]],  # population 2/site 1
        [[0, 0, 0, 0], [0, 1, 0, 0], [1, 1, 0, 0]],
    ])

    num, den = hudson_fst(
        ac1=pop1_gt_array.count_alleles(),
        ac2=pop2_gt_array.count_alleles(),
    )

    expected_fst = num.sum() / den.sum()

    # `calc_fst` accepts an aggregate array with all populations combined
    combined_gt_array = pop1_gt_array.concatenate(pop2_gt_array, axis=1)
    assert combined_gt_array.n_samples == pop1_gt_array.n_samples + pop2_gt_array.n_samples
    assert combined_gt_array.n_variants == pop1_gt_array.n_variants == pop2_gt_array.n_variants

    result: FstResult = calc_fst(
        gt_array_fst=combined_gt_array,
        fst_pop_indicies=[
            [0, 1, 2],  # population 1
            [3, 4, 5],  # population 2
        ],
        fst_type=FSTEstimator.HUDSON,
    )

    assert result.fst == pytest.approx(expected_fst)
    assert result.a == pytest.approx(num.sum())
    assert result.b == pytest.approx(den.sum())
    assert result.c == 0
    assert result.n_sites == combined_gt_array.n_variants


def test_calc_fst_wc_arbitrary_ploidy_raises() -> None:
    """Ensure we raise a `NotImplementedError` if FST_WC is requested in non-diploid genomes."""
    pop1_gt_array = GenotypeArray([
        [[0, 0, 0, 0], [0, 1, 0, 0], [1, 1, 0, 0]],
        [[0, 1, 0, 0], [0, 1, 0, 0], [1, 1, 1, 1]],  # population 2/site 1
    ])  # population 1/site 1

    pop2_gt_array = GenotypeArray([
        [[0, 1, 0, 0], [0, 1, 0, 0], [1, 1, 1, 1]],  # population 2/site 1
        [[0, 0, 0, 0], [0, 1, 0, 0], [1, 1, 0, 0]],
    ])
    combined_gt_array = pop1_gt_array.concatenate(pop2_gt_array, axis=1)

    with pytest.raises(NotImplementedError):
        calc_fst(
            gt_array_fst=combined_gt_array,
            fst_pop_indicies=[
                [0, 1, 2],  # population 1
                [3, 4, 5],  # population 2
            ],
            fst_type=FSTEstimator.WC,
        )


def test_calc_fst_hudson() -> None:
    """Compare FST calculation to the `scikit-allel` implementation."""
    # scikit-allel's implementation accepts one array for each population
    pop1_gt_array = GenotypeArray([
        [[0, 0], [0, 1], [1, 1]],  # population 1/site 1
        [[0, 1], [0, 0], [1, 1]],  # population 1/site 2
    ])

    pop2_gt_array = GenotypeArray([
        [[0, 1], [0, 1], [1, 1]],  # population 2/site 1
        [[0, 1], [1, 1], [1, 1]],  # population 2/site 2
    ])

    num, den = hudson_fst(
        ac1=pop1_gt_array.count_alleles(),
        ac2=pop2_gt_array.count_alleles(),
    )

    expected_fst = num.sum() / den.sum()

    # `calc_fst` accepts an aggregate array with all populations combined
    combined_gt_array = pop1_gt_array.concatenate(pop2_gt_array, axis=1)
    assert combined_gt_array.n_samples == pop1_gt_array.n_samples + pop2_gt_array.n_samples
    assert combined_gt_array.n_variants == pop1_gt_array.n_variants == pop2_gt_array.n_variants

    result: FstResult = calc_fst(
        gt_array_fst=combined_gt_array,
        fst_pop_indicies=[
            [0, 1, 2],  # population 1
            [3, 4, 5],  # population 2
        ],
        fst_type=FSTEstimator.HUDSON,
    )

    assert result.fst == pytest.approx(expected_fst)
    assert result.a == pytest.approx(num.sum())
    assert result.b == pytest.approx(den.sum())
    assert result.c == 0
    assert result.n_sites == combined_gt_array.n_variants


def test_calc_fst_hudson_single_locus() -> None:
    """Compare FST calculation to the `scikit-allel` implementation."""
    # scikit-allel's implementation accepts one array for each population
    pop1_gt_array = GenotypeArray([
        [[0, 0], [0, 1], [1, 1]],  # population 1/site 1
    ])

    pop2_gt_array = GenotypeArray(
        [[[0, 1], [0, 1], [1, 1]]]  # population 2/site 1
    )

    num, den = hudson_fst(
        ac1=pop1_gt_array.count_alleles(),
        ac2=pop2_gt_array.count_alleles(),
    )

    expected_fst = num.sum() / den.sum()

    # `calc_fst` accepts an aggregate array with all populations combined
    combined_gt_array = pop1_gt_array.concatenate(pop2_gt_array, axis=1)
    assert combined_gt_array.n_samples == pop1_gt_array.n_samples + pop2_gt_array.n_samples
    assert combined_gt_array.n_variants == pop1_gt_array.n_variants == pop2_gt_array.n_variants

    result: FstResult = calc_fst(
        gt_array_fst=combined_gt_array,
        fst_pop_indicies=[
            [0, 1, 2],  # population 1
            [3, 4, 5],  # population 2
        ],
        fst_type=FSTEstimator.HUDSON,
    )

    assert result.fst == pytest.approx(expected_fst)
    assert result.a == pytest.approx(num.sum())
    assert result.b == pytest.approx(den.sum())
    assert result.c == 0
    assert result.n_sites == combined_gt_array.n_variants


def test_calc_fst_wc() -> None:
    """Compare FST calculation to the `scikit-allel` implementation."""
    # scikit-allel's implementation accepts one array for each population
    pop1_gt_array = GenotypeArray([
        [[0, 0], [0, 1], [1, 1]],  # population 1/site 1
        [[0, 1], [0, 0], [1, 1]],  # population 1/site 2
    ])

    pop2_gt_array = GenotypeArray([
        [[0, 1], [0, 1], [1, 1]],  # population 2/site 1
        [[0, 1], [1, 1], [1, 1]],  # population 2/site 2
    ])

    # both the scikit-allel implementation and `calc_fst` accept an aggregate array with all
    # populations combined
    combined_gt_array = pop1_gt_array.concatenate(pop2_gt_array, axis=1)
    assert combined_gt_array.n_samples == pop1_gt_array.n_samples + pop2_gt_array.n_samples
    assert combined_gt_array.n_variants == pop1_gt_array.n_variants == pop2_gt_array.n_variants

    subpops = [
        [0, 1, 2],  # population 1
        [3, 4, 5],  # population 2
    ]

    a, b, c = weir_cockerham_fst(
        g=combined_gt_array,
        subpops=subpops,
        max_allele=1,  # TODO add test coverage over multiallelic sites
    )

    expected_fst = a.sum() / (a.sum() + b.sum() + c.sum())

    result: FstResult = calc_fst(
        gt_array_fst=combined_gt_array,
        fst_pop_indicies=subpops,
        fst_type=FSTEstimator.WC,
    )

    assert result.fst == pytest.approx(expected_fst)
    assert result.a == pytest.approx(a.sum())
    assert result.b == pytest.approx(b.sum())
    assert result.c == pytest.approx(c.sum())
    assert result.n_sites == combined_gt_array.n_variants


def test_calc_fst_wc_single_locus() -> None:
    """Compare FST calculation to the `scikit-allel` implementation."""
    # scikit-allel's implementation accepts one array for each population
    pop1_gt_array = GenotypeArray([
        [[0, 0], [0, 1], [1, 1]],  # population 1/site 1
    ])

    pop2_gt_array = GenotypeArray(
        [[[0, 1], [0, 1], [1, 1]]]  # population 2/site 1
    )

    # both the scikit-allel implementation and `calc_fst` accept an aggregate array with all
    # populations combined
    combined_gt_array = pop1_gt_array.concatenate(pop2_gt_array, axis=1)
    assert combined_gt_array.n_samples == pop1_gt_array.n_samples + pop2_gt_array.n_samples
    assert combined_gt_array.n_variants == pop1_gt_array.n_variants == pop2_gt_array.n_variants

    subpops = [
        [0, 1, 2],  # population 1
        [3, 4, 5],  # population 2
    ]

    a, b, c = weir_cockerham_fst(
        g=combined_gt_array,
        subpops=subpops,
        max_allele=1,  # TODO add test coverage over multiallelic sites
    )

    expected_fst = a.sum() / (a.sum() + b.sum() + c.sum())

    result: FstResult = calc_fst(
        gt_array_fst=combined_gt_array,
        fst_pop_indicies=subpops,
        fst_type=FSTEstimator.WC,
    )

    assert result.fst == pytest.approx(expected_fst)
    assert result.a == pytest.approx(a.sum())
    assert result.b == pytest.approx(b.sum())
    assert result.c == pytest.approx(c.sum())
    assert result.n_sites == combined_gt_array.n_variants


def test_calc_pi_single_locus_multiallelic() -> None:
    """Test that our Pi calculation matches the expected values based on the paper's formula."""
    array = GenotypeArray([
        [[0, 0], [0, 2], [1, 1]],  # site 1
    ])

    # sanity
    assert array.n_samples == 3
    assert array.n_variants == 1
    assert array.ploidy == 2

    result: PiResult = calc_pi(array)

    # Section 1.1 of Korunes & Samuk
    # https://pmc.ncbi.nlm.nih.gov/articles/PMC8044049/
    expected_diffs = 11  # 3 zeros * 2 ones + 3 zeros * 1 two + 2 ones * 1 two
    expected_comps = 15  # (6 choose 2) = 15 comparisons

    assert result.total_diffs == expected_diffs
    assert result.total_comps == expected_comps
    assert result.total_missing == 0
    assert result.avg_pi == pytest.approx(expected_diffs / expected_comps)


def test_calc_dxy_single_locus_multiallelic() -> None:
    """Test that our Dxy calculation matches the expected values based on the paper's formula."""
    pop1_gt_array = GenotypeArray([
        [[0, 0], [0, 1], [1, 2]],  # population 1/site 1
    ])

    pop2_gt_array = GenotypeArray([
        [[0, 1], [0, 2], [1, 1]],  # population 2/site 1
    ])

    # sanity
    assert pop1_gt_array.n_samples == pop2_gt_array.n_samples == 3
    assert pop1_gt_array.n_variants == pop2_gt_array.n_variants == 1
    assert pop1_gt_array.ploidy == pop2_gt_array.ploidy == 2  # diploid

    result: DxyResult = calc_dxy(pop1_gt_array=pop1_gt_array, pop2_gt_array=pop2_gt_array)

    # Section 1.1 of Korunes & Samuk
    # https://pmc.ncbi.nlm.nih.gov/articles/PMC8044049/
    expected_diffs = 23  # 3 * 3 + 3 * 1 + 2 * 2 + 2 * 1 + 1 * 2 + 1 * 3
    expected_comps = 6 * 6  # number of pairwise comparisons among alleles (6 per pop)

    assert result.total_diffs == expected_diffs
    assert result.total_comps == expected_comps
    assert result.total_missing == 0
    assert result.avg_dxy == pytest.approx(expected_diffs / expected_comps)


def test_calc_dxy_multiallelic() -> None:
    """Test that our Dxy calculation matches the expected values based on the paper's formula."""
    pop1_gt_array = GenotypeArray([
        [[0, 0], [0, 1], [1, 2]],  # population 1/site 1
        [[0, 0], [0, 2], [2, 2]],  # population 1/site 2
    ])

    pop2_gt_array = GenotypeArray([
        [[0, 1], [0, 2], [1, 1]],  # population 2/site 1
        [[0, 1], [1, 2], [0, 1]],  # population 2/site 2
    ])

    # sanity
    assert pop1_gt_array.n_samples == pop2_gt_array.n_samples == 3
    assert pop1_gt_array.n_variants == pop2_gt_array.n_variants == 2
    assert pop1_gt_array.ploidy == pop2_gt_array.ploidy == 2  # diploid

    result: DxyResult = calc_dxy(pop1_gt_array=pop1_gt_array, pop2_gt_array=pop2_gt_array)

    # Section 1.1 of Korunes & Samuk
    # https://pmc.ncbi.nlm.nih.gov/articles/PMC8044049/
    expected_diffs_site_1 = 3 * 3 + 3 * 1 + 2 * 2 + 2 * 1 + 1 * 2 + 1 * 3
    expected_diffs_site_2 = 3 * 3 + 3 * 1 + 0 * 2 + 0 * 1 + 3 * 2 + 3 * 3

    expected_diffs = expected_diffs_site_1 + expected_diffs_site_2
    expected_comps = 6 * 6 * 2  # 6 haploids * 6 haploids * 2 sites

    assert result.total_diffs == expected_diffs
    assert result.total_comps == expected_comps
    assert result.total_missing == 0
    assert result.avg_dxy == pytest.approx(expected_diffs / expected_comps)


def test_calc_fst_hudson_multiallelic() -> None:
    """Compare FST calculation to the `scikit-allel` implementation."""
    # scikit-allel's implementation accepts one array for each population
    pop1_gt_array = GenotypeArray([
        [[0, 0], [0, 1], [1, 2]],  # population 1/site 1
        [[0, 0], [0, 2], [2, 2]],  # population 1/site 2
    ])

    pop2_gt_array = GenotypeArray([
        [[0, 1], [0, 2], [1, 1]],  # population 2/site 1
        [[0, 1], [1, 2], [0, 1]],  # population 2/site 2
    ])

    num, den = hudson_fst(
        ac1=pop1_gt_array.count_alleles(),
        ac2=pop2_gt_array.count_alleles(),
    )

    expected_fst = num.sum() / den.sum()

    # `calc_fst` accepts an aggregate array with all populations combined
    combined_gt_array = pop1_gt_array.concatenate(pop2_gt_array, axis=1)
    assert combined_gt_array.n_samples == pop1_gt_array.n_samples + pop2_gt_array.n_samples
    assert combined_gt_array.n_variants == pop1_gt_array.n_variants == pop2_gt_array.n_variants

    result: FstResult = calc_fst(
        gt_array_fst=combined_gt_array,
        fst_pop_indicies=[
            [0, 1, 2],  # population 1
            [3, 4, 5],  # population 2
        ],
        fst_type=FSTEstimator.HUDSON,
    )

    assert result.fst == pytest.approx(expected_fst)
    assert result.a == pytest.approx(num.sum())
    assert result.b == pytest.approx(den.sum())
    assert result.c == 0
    assert result.n_sites == combined_gt_array.n_variants


def test_calc_fst_wc_multiallelic() -> None:
    """Compare FST calculation to the `scikit-allel` implementation."""
    # scikit-allel's implementation accepts one array for each population
    pop1_gt_array = GenotypeArray([
        [[0, 0], [0, 1], [1, 2]],  # population 1/site 1
        [[0, 0], [0, 2], [2, 2]],  # population 1/site 2
    ])

    pop2_gt_array = GenotypeArray([
        [[0, 1], [0, 2], [1, 1]],  # population 2/site 1
        [[0, 1], [1, 2], [0, 1]],  # population 2/site 2
    ])

    # both the scikit-allel implementation and `calc_fst` accept an aggregate array with all
    # populations combined
    combined_gt_array = pop1_gt_array.concatenate(pop2_gt_array, axis=1)
    assert combined_gt_array.n_samples == pop1_gt_array.n_samples + pop2_gt_array.n_samples
    assert combined_gt_array.n_variants == pop1_gt_array.n_variants == pop2_gt_array.n_variants

    subpops = [
        [0, 1, 2],  # population 1
        [3, 4, 5],  # population 2
    ]

    a, b, c = weir_cockerham_fst(
        g=combined_gt_array,
        subpops=subpops,
        max_allele=None,
    )

    expected_fst = a.sum() / (a.sum() + b.sum() + c.sum())

    result: FstResult = calc_fst(
        gt_array_fst=combined_gt_array,
        fst_pop_indicies=subpops,
        fst_type=FSTEstimator.WC,
    )

    assert result.fst == pytest.approx(expected_fst)
    assert result.a == pytest.approx(a.sum())
    assert result.b == pytest.approx(b.sum())
    assert result.c == pytest.approx(c.sum())
    assert result.n_sites == combined_gt_array.n_variants


##############################
# Watterson's Theta Unit-tests

# Expected values for these tests were generated based on math contained
# in https://github.com/ksamuk/pixy/pull/117. Some tests check the formulaic result as well.
##############################


def test_calc_watterson_theta_single_locus() -> None:
    """
    Assert that Watterson's theta calculation produces known outputs with known inputs.

    In this case we test a single locus. Expected values were generated based on math contained
    in https://github.com/ksamuk/pixy/pull/117.

    Additionally, because there is no missing data, we compare
    to the `scikit-allel` implementation and assert that the results match.
    """
    array = GenotypeArray([
        [[0, 0], [0, 1], [1, 1]],  # site 1, 3 samples
    ])
    # formulaic Watterson's theta
    ac: AlleleCountsArray = array.count_alleles()  # [3,3]
    num_segregating_sites: int = ac.count_segregating()  # just 1 site is polymorphic
    n: int = ac.sum(axis=1).max()  # [3,3] -> 6
    # calculate a_n (the reciprocal sums from 1 to n-1)
    a_n = np.sum([1 / i for i in range(1, n)])  # 2.283333333333333
    formulaic_raw_theta: float = num_segregating_sites / a_n  # 1 / 2.283
    watterson_result: WattersonThetaResult = calc_watterson_theta(array)
    # compare `pixy` result to both scikit-allel and expected output
    assert watterson_result.avg_theta == pytest.approx(
        watterson_theta(ac=array.count_alleles(), pos=[1])
    )
    assert watterson_result.avg_theta == pytest.approx(0.43795620437956206)
    # compare `pixy` result to formulaic result and expected output
    assert watterson_result.raw_theta == pytest.approx(formulaic_raw_theta)
    assert watterson_result.raw_theta == pytest.approx(0.43795620437956206)
    assert watterson_result.num_weighted_sites == 1.0


def test_calc_watterson_theta_haploid_singleton() -> None:
    """Test calculation of Watterson's Theta for haploid genomes with singletons."""
    array = GenotypeArray([
        [[0], [0]],
        [[1], [1]],
        [[0], [1]],
        [[0], [0]],
        [[1], [1]],
        [[0], [-1]],
        [[1], [-1]],
        [[-1], [-1]],
    ])
    watterson_result = calc_watterson_theta(array)

    # expected values
    assert np.isinf(watterson_result.avg_theta)
    assert np.isinf(watterson_result.raw_theta)
    assert watterson_result.num_weighted_sites == 6.0


def test_calc_watterson_theta_diploid_biallelic() -> None:
    """
    Assert that Watterson's Theta calculation produces known outputs with known inputs.

    In this case we test diploid genomes with biallelic sites only.
    """
    array = GenotypeArray([
        [[0, 0], [0, 1], [1, 1]],  # site 1
        [[0, 1], [0, 1], [1, 1]],  # site 2
    ])
    # formulaic Watterson's theta
    ac: AlleleCountsArray = array.count_alleles()  # [3,3] at site 1, [2,4] at site 2
    num_segregating_sites: int = ac.count_segregating()  # 2 sites are polymorphic
    n: int = ac.sum(axis=1).max()  # 6
    # calculate a_n (the reciprocal sums from 1 to n-1)
    a_n = np.sum([1 / i for i in range(1, n)])  # 2.283333333333333
    formulaic_raw_theta: float = num_segregating_sites / a_n  # 2 / 2.283

    watterson_result = calc_watterson_theta(array)
    assert watterson_result.avg_theta == pytest.approx(
        watterson_theta(ac=array.count_alleles(), pos=[1, 2])
    )
    assert watterson_result.avg_theta == pytest.approx(0.43795620437956206)

    assert watterson_result.raw_theta == pytest.approx(formulaic_raw_theta)
    assert watterson_result.raw_theta == pytest.approx(0.8759124087591241)

    assert watterson_result.num_weighted_sites == 2.0


def test_calc_watterson_theta_diploid_multiallelic() -> None:
    """
    Assert that Watterson's Theta calculation produces known outputs with known inputs.

    In this case we test diploid genomes with four sites, some of which are multiallelic.
    Because there is no missing data here, we compare to the `scikit-allel` implementation.
    """
    array = GenotypeArray([
        [[0, 0], [0, 1], [1, 1]],  # site 1
        [[0, 2], [0, 1], [1, 1]],  # site 2
        [[0, 0], [0, 0], [0, 0]],  # site 3
        [[0, 2], [0, 1], [1, 1]],  # site 4
    ])

    result = calc_watterson_theta(array)

    assert result.avg_theta == pytest.approx(
        watterson_theta(ac=array.count_alleles(), pos=[1, 2, 3, 4])
    )
    assert result.raw_theta == pytest.approx(1.3138686131386863)
    assert result.num_weighted_sites == 4


def test_calc_watterson_theta_diploid_multiallelic_missing_data() -> None:
    """
    Assert that Watterson's Theta calculation produces known outputs with known inputs.

    In this case we test diploid genomes with two multiallelic sites. Importantly, this genotype
    array has missing data, and so the results from `pixy` are only compared to the formulaic
    Watterson's Theta and not the `scikit-allel` implementation (which cannot handle missing data).
    """
    array = GenotypeArray([
        [[0, 0], [0, 1], [-1, 1]],  # site 1 (missing data)
        [[0, 2], [0, 1], [1, 1]],  # site 2 (multiallelic)
        [[0, 0], [0, 0], [0, 0]],  # site 3 (all ref, non-segregating)
        [[0, 2], [0, 1], [1, 1]],  # site 4 (multiallelic)
    ])

    formulaic_weighted_raw_theta: float = (1 / 2.0833) + (2 / 2.283333)
    formulaic_weighted_avg_theta: float = formulaic_weighted_raw_theta / 4  # num_sites

    watterson_result = calc_watterson_theta(array)
    # avg theta comparisons
    assert watterson_result.avg_theta == pytest.approx(0.3389781021897811)
    assert watterson_result.avg_theta == pytest.approx(formulaic_weighted_avg_theta, rel=1e-05)
    # raw theta comparisons
    assert watterson_result.raw_theta == pytest.approx(1.3559124087591243)
    assert watterson_result.raw_theta == pytest.approx(formulaic_weighted_raw_theta, rel=1e-05)
    assert (
        watterson_result.num_weighted_sites == 3.8333333333333335
    )  # np.sum(array([1, 3]) * array([0.83333333, 3.]))


def test_calc_watterson_theta_tetraploidy() -> None:
    """
    Assert that Watterson's Theta calculation produces known outputs with known inputs.

    In this case we test tetraploid genomes.
    """
    array = GenotypeArray([
        [[0, 0, 0, 0], [0, 1, 0, 0], [1, 1, 1, 1]],  # site 1
        [[0, 1, 0, 0], [0, 1, 1, 1], [1, 1, 0, 0]],  # site 2
    ])

    # formulaic Watterson's theta
    ac: AlleleCountsArray = array.count_alleles()  # [7,5], [6,6]
    num_segregating_sites: int = ac.count_segregating()  # 2 sites are polymorphic
    n: int = ac.sum(axis=1).max()  # 12!
    # calculate a_n (the reciprocal sums from 1 to n-1)
    a_n = np.sum([1 / i for i in range(1, n)])  # 3448773446
    formulaic_raw_theta: float = num_segregating_sites / a_n  # 2 / 3.01987
    result = calc_watterson_theta(array)
    assert result.avg_theta == pytest.approx(0.3311392767975535)
    assert result.raw_theta == pytest.approx(0.662278553595107) == formulaic_raw_theta
    assert result.num_weighted_sites == 2.0


##############################
# Tajima's D Unit-tests
#
# Expected values for these tests were generated based on math contained
# in https://github.com/ksamuk/pixy/pull/117.
##############################


def test_calc_tajima_d_single_locus() -> None:
    """
    Assert that Tajima's D calculation produces known outputs with known inputs.

    In this case we are testing a single locus.
    """
    array = GenotypeArray([
        [[0, 0], [0, 1], [1, 1]],  # site 1
    ])

    result: TajimaDResult = calc_tajima_d(array)

    # Without missing data, the scikit-allel implementation should match pixy's
    ac = array.count_alleles()
    assert result.tajima_d == pytest.approx(allel.tajima_d(ac=ac, min_sites=0))
    assert result.raw_pi == pytest.approx(allel.mean_pairwise_difference(ac=ac).sum())
    assert result.watterson_theta == pytest.approx(allel.watterson_theta(pos=[1], ac=ac))

    # There is no standalone function or helper for the denominator - this was manually calculated
    # from scikit-allel's code
    assert result.d_stdev == pytest.approx(0.1121335)

@pytest.mark.xfail(reason="Updating assertion logic")
def test_calc_tajima_d_haploid_singleton() -> None:
    """
    Assert that Tajima's D calculation produces known outputs with known inputs.

    In this case we are testing haploid genomes with some singleton values. This test case was
    mentioned in the original PR to pixy.
    """
    array = GenotypeArray([
        [[0], [0]],
        [[1], [1]],
        [[0], [1]],
        [[0], [0]],
        [[1], [1]],
        [[0], [-1]],
        [[1], [-1]],
        [[-1], [-1]],
    ])

    result: TajimaDResult = calc_tajima_d(array)

    assert result.tajima_d == "NA"
    assert result.raw_pi == pytest.approx(1.0)
    assert np.isinf(result.watterson_theta)
    assert result.d_stdev == pytest.approx(0.0)


def test_calc_tajima_d_diploid_biallelic() -> None:
    """
    Assert that Tajima's D calculation produces known outputs with known inputs.

    In this case we test diploid genomes with biallelic sites only.
    """
    array = GenotypeArray([
        [[0, 0], [0, 1], [1, 1]],  # site 1
        [[0, 1], [0, 1], [1, 1]],  # site 2
    ])

    result: TajimaDResult = calc_tajima_d(array)

    # Without missing data, the scikit-allel implementation should match pixy's
    ac = array.count_alleles()
    assert result.tajima_d == pytest.approx(allel.tajima_d(ac=ac, min_sites=0))
    assert result.raw_pi == pytest.approx(allel.mean_pairwise_difference(ac=ac).sum())

    # scikit-allel's function returns the average theta over the number of bases - to obtain the
    # raw theta, multiply it back out
    assert result.watterson_theta == pytest.approx(allel.watterson_theta(pos=[1, 2], ac=ac) * 2)

    # There is no standalone function or helper for the denominator - this was manually calculated
    # from scikit-allel's code
    assert result.d_stdev == pytest.approx(0.18485059)


def test_calc_tajima_d_diploid_multiallelic() -> None:
    """
    Assert that Tajima's D calculation produces known outputs with known inputs.

    Here we test diploid genomes with a multiallelic spike-in.
    """
    array = GenotypeArray([
        [[0, 0], [0, 1], [1, 1]],  # site 1
        [[0, 2], [0, 1], [1, 1]],  # site 2
    ])

    result: TajimaDResult = calc_tajima_d(array)

    # Without missing data, the scikit-allel implementation should match pixy's
    ac = array.count_alleles()
    assert result.tajima_d == pytest.approx(allel.tajima_d(ac=ac, min_sites=0))
    assert result.raw_pi == pytest.approx(allel.mean_pairwise_difference(ac=ac).sum())

    # scikit-allel's function returns the average theta over the number of bases - to obtain the
    # raw theta, multiply it back out
    assert result.watterson_theta == pytest.approx(allel.watterson_theta(pos=[1, 2], ac=ac) * 2)

    # There is no standalone function or helper for the denominator - this was manually calculated
    # from scikit-allel's code
    assert result.d_stdev == pytest.approx(0.1848505)


def test_calc_tajima_d_tetraploidy() -> None:
    """Assert that Tajima's D calculation produces known outputs with known inputs (tetraploidy)."""
    array = GenotypeArray([
        [[0, 0, 0, 0], [0, 1, 0, 0], [1, 1, 1, 1]],  # site 1
        [[0, 1, 0, 0], [0, 1, 1, 1], [1, 1, 0, 0]],  # site 2
    ])

    result: TajimaDResult = calc_tajima_d(array)

    # Without missing data, the scikit-allel implementation should match pixy's
    ac = array.count_alleles()
    assert result.tajima_d == pytest.approx(allel.tajima_d(ac=ac, min_sites=0))
    assert result.raw_pi == pytest.approx(allel.mean_pairwise_difference(ac=ac).sum())

    # scikit-allel's function returns the average theta over the number of bases - to obtain the
    # raw theta, multiply it back out
    assert result.watterson_theta == pytest.approx(allel.watterson_theta(pos=[1, 2], ac=ac) * 2)

    # There is no standalone function or helper for the denominator - this was manually calculated
    # from scikit-allel's code
    assert result.d_stdev == pytest.approx(0.2266425)
