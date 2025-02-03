import pytest
from allel import GenotypeArray
from allel import hudson_fst
from allel import weir_cockerham_fst

from pixy.calc import calc_dxy
from pixy.calc import calc_fst
from pixy.calc import calc_pi
from pixy.enums import FSTEstimator
from pixy.models import DxyResult
from pixy.models import FstResult
from pixy.models import PiResult


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
