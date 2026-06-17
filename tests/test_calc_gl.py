"""Unit tests for the GL-based pi estimator in `pixy.calc_gl`."""

import numpy as np
import pytest
from allel import GenotypeArray

from pixy.calc import calc_pi
from pixy.calc_gl import calc_pi_gl
from pixy.calc_gl import likelihoods_to_posteriors
from pixy.calc_gl import reduces_to_hardcall
from pixy.models import PiResult


def _gt_to_confident_pl(gt: np.ndarray) -> np.ndarray:
    """
    Synthesize confident PL triplets from a hard-call genotype array.

    `gt` has shape (n_sites, n_samples, 2) with values in {0, 1} (no missing). Returns a
    PL array of shape (n_sites, n_samples, 3) where each sample's triplet places PL=0 on
    the actual genotype class and PL=255 on the other two — effectively a one-hot posterior
    after conversion via `likelihoods_to_posteriors`.
    """
    n_alt = gt.sum(axis=-1)  # (n_sites, n_samples), values in {0, 1, 2}
    pl = np.full(n_alt.shape + (3,), 255, dtype=np.int16)
    np.put_along_axis(pl, n_alt[..., None], 0, axis=-1)
    return pl


def test_likelihoods_to_posteriors_normalizes_pl() -> None:
    """Random PL → posteriors that sum to 1 with no rows flagged missing."""
    rng = np.random.default_rng(0)
    pl = rng.integers(0, 50, size=(5, 4, 3), dtype=np.int16)
    posteriors, missing = likelihoods_to_posteriors(pl, "PL")
    assert posteriors.shape == (5, 4, 3)
    assert missing.shape == (5, 4)
    assert not missing.any()
    np.testing.assert_allclose(posteriors.sum(axis=-1), 1.0, atol=1e-12)


def test_likelihoods_to_posteriors_normalizes_gl() -> None:
    """Random GL → posteriors that sum to 1 with no rows flagged missing."""
    rng = np.random.default_rng(1)
    gl = -rng.random(size=(5, 4, 3)) * 5.0
    posteriors, missing = likelihoods_to_posteriors(gl, "GL")
    assert posteriors.shape == (5, 4, 3)
    assert not missing.any()
    np.testing.assert_allclose(posteriors.sum(axis=-1), 1.0, atol=1e-12)


def test_likelihoods_to_posteriors_missing_pl_flagged() -> None:
    """Rows of all-(-1) and all-zero PL trip the missing-row sentinel."""
    pl = np.array(
        [
            [[0, 30, 60], [-1, -1, -1], [0, 0, 0]],
        ],
        dtype=np.int16,
    )
    posteriors, missing = likelihoods_to_posteriors(pl, "PL")
    assert missing.shape == (1, 3)
    assert missing[0, 0] == False  # noqa: E712 — explicit elementwise check
    assert missing[0, 1] == True  # noqa: E712
    assert missing[0, 2] == True  # noqa: E712
    np.testing.assert_array_equal(posteriors[0, 1], np.zeros(3))
    np.testing.assert_array_equal(posteriors[0, 2], np.zeros(3))


def test_likelihoods_to_posteriors_missing_gl_flagged() -> None:
    """Rows of NaN and all-equal GL trip the missing-row sentinel."""
    gl = np.array(
        [
            [[0.0, -1.0, -3.0], [np.nan, np.nan, np.nan], [0.0, 0.0, 0.0]],
        ],
        dtype=np.float64,
    )
    _posteriors, missing = likelihoods_to_posteriors(gl, "GL")
    assert missing[0, 0] == False  # noqa: E712
    assert missing[0, 1] == True  # noqa: E712
    assert missing[0, 2] == True  # noqa: E712


def test_likelihoods_to_posteriors_rejects_bad_shape() -> None:
    """Input that isn't (n_sites, n_samples, 3) is a usage error."""
    with pytest.raises(ValueError, match="shape"):
        likelihoods_to_posteriors(np.zeros((5, 4, 2), dtype=np.int16), "PL")


def test_likelihoods_to_posteriors_rejects_bad_kind() -> None:
    """Anything other than PL or GL is rejected at the boundary."""
    with pytest.raises(ValueError, match="lik_kind"):
        likelihoods_to_posteriors(np.zeros((1, 1, 3), dtype=np.int16), "BEAGLE")


def test_confident_pl_reduces_to_hardcall_pi_paper_example() -> None:
    """
    Confident PLs reduce GL pi to hard-call pi on the Korunes & Samuk 2021 paper example.

    Reuses the 3-sample 2-site case from `tests/test_calc.py::test_calc_pi` so the expected
    integer totals (17 diffs, 30 comps) are already verified against the paper formula.
    """
    gt = np.array([
        [[0, 0], [0, 1], [1, 1]],  # site 1
        [[0, 1], [0, 1], [1, 1]],  # site 2
    ])
    array = GenotypeArray(gt)
    hard_result: PiResult = calc_pi(array)

    pl = _gt_to_confident_pl(gt)
    posteriors, missing = likelihoods_to_posteriors(pl, "PL")
    # Sanity: the reducer says we're one-hot at every site.
    assert reduces_to_hardcall(posteriors).all()

    gl_result: PiResult = calc_pi_gl(posteriors, missing, n_haps=array.n_samples * array.ploidy)

    assert gl_result.total_diffs == pytest.approx(hard_result.total_diffs, abs=1e-9)
    assert gl_result.total_comps == pytest.approx(hard_result.total_comps, abs=1e-9)
    assert gl_result.total_missing == pytest.approx(hard_result.total_missing, abs=1e-9)
    assert gl_result.avg_pi == pytest.approx(hard_result.avg_pi)


@pytest.mark.parametrize(
    "gt",
    [
        np.array([[[0, 0], [0, 1], [1, 1]]]),  # single site, balanced
        np.array([[[0, 0], [0, 0], [0, 0]]]),  # all hom-ref (zero diffs)
        np.array([[[1, 1], [1, 1], [1, 1]]]),  # all hom-alt (zero diffs)
        np.array([[[0, 1], [0, 1], [0, 1]]]),  # all het
        np.array(  # multi-site mixed
            [
                [[0, 0], [1, 1], [0, 1]],
                [[0, 1], [0, 0], [1, 1]],
                [[0, 0], [0, 0], [0, 1]],
            ]
        ),
    ],
)
def test_confident_pl_matches_hardcall_param(gt: np.ndarray) -> None:
    """Property check across several configurations: confident PL → identical to hard call."""
    array = GenotypeArray(gt)
    hard_result = calc_pi(array)

    pl = _gt_to_confident_pl(gt)
    posteriors, missing = likelihoods_to_posteriors(pl, "PL")
    gl_result = calc_pi_gl(posteriors, missing, n_haps=array.n_samples * array.ploidy)

    assert gl_result.total_diffs == pytest.approx(hard_result.total_diffs, abs=1e-9)
    assert gl_result.total_comps == pytest.approx(hard_result.total_comps, abs=1e-9)
    assert gl_result.total_missing == pytest.approx(hard_result.total_missing, abs=1e-9)


def test_hand_crafted_posterior_pi_two_individuals() -> None:
    """
    Hand-derive pi for two diploid individuals at one site against the formula in calc_gl.

    Two diploid individuals at one site, hand-crafted posteriors:

      indiv 1: (0.5, 0.5, 0.0)   -> mu_1 = 0.5, h_1 = 0.5
      indiv 2: (0.0, 0.5, 0.5)   -> mu_2 = 1.5, h_2 = 0.5

    Then S = 2.0, H = 1.0, Q = 0.25 + 2.25 = 2.5, n = 2.

    diffs_site = H + 2*(n-1)*S - S**2 + Q
               = 1.0 + 2*1*2.0 - 4.0 + 2.5
               = 1.0 + 4.0 - 4.0 + 2.5
               = 3.5

    comps_site = 2*n*(2*n - 1) / 2 = 2*2*3/2 = 6.0
    """
    posteriors = np.array(
        [[[0.5, 0.5, 0.0], [0.0, 0.5, 0.5]]],
        dtype=np.float64,
    )
    missing = np.zeros((1, 2), dtype=bool)
    result = calc_pi_gl(posteriors, missing, n_haps=4)

    assert result.total_diffs == pytest.approx(3.5)
    assert result.total_comps == pytest.approx(6.0)
    assert result.total_missing == pytest.approx(0.0)
    assert result.avg_pi == pytest.approx(3.5 / 6.0)


def test_uniform_posterior_pi_closed_form() -> None:
    """
    All-uniform posterior (1/3, 1/3, 1/3) for every sample at a single site.

    For one sample: mu_i = 1/3 + 2/3 = 1.0, h_i = 1/3, mu_i**2 = 1.0.
    With n samples: S = n, H = n/3, Q = n.

    diffs_site = n/3 + 2*(n-1)*n - n**2 + n
               = n/3 + 2*n**2 - 2*n - n**2 + n
               = n**2 - n + n/3
               = n**2 - 2n/3

    For n=3: 9 - 2 = 7.
    comps_site = 2n*(2n-1)/2 = 6*5/2 = 15.
    """
    n = 3
    posteriors = np.full((1, n, 3), 1.0 / 3.0, dtype=np.float64)
    missing = np.zeros((1, n), dtype=bool)
    result = calc_pi_gl(posteriors, missing, n_haps=2 * n)

    assert result.total_diffs == pytest.approx(n * n - 2 * n / 3.0)
    assert result.total_comps == pytest.approx(2 * n * (2 * n - 1) / 2)


def test_missing_likelihoods_reduce_effective_n() -> None:
    """
    Missing PL on one individual drops it from S/H/Q sums and from per-site comps.

    Two individuals at one site; first is fully missing (PL = -1). The site behaves
    as if it has only one informative individual, contributing 0 to diffs and 1 to comps
    (since a het individual yields 1 within-pair difference and 1 within-pair comparison).
    """
    # indiv 1: missing PL; indiv 2: confident het
    pl = np.array([[[-1, -1, -1], [255, 0, 255]]], dtype=np.int16)
    posteriors, missing = likelihoods_to_posteriors(pl, "PL")
    assert missing[0, 0] == True  # noqa: E712
    assert missing[0, 1] == False  # noqa: E712

    result = calc_pi_gl(posteriors, missing, n_haps=4)
    # Only the het individual contributes: 1 within-individual diff, 1 within-individual comp.
    assert result.total_diffs == pytest.approx(1.0)
    assert result.total_comps == pytest.approx(1.0)
    # Possible comps at this site = 4*3/2 = 6; missing = 6 - 1 = 5.
    assert result.total_missing == pytest.approx(5.0)


def test_all_missing_returns_na_pi() -> None:
    """All-missing PL yields zero diffs/comps and avg_pi = NA."""
    pl = np.full((2, 3, 3), -1, dtype=np.int16)
    posteriors, missing = likelihoods_to_posteriors(pl, "PL")
    assert missing.all()
    result = calc_pi_gl(posteriors, missing, n_haps=6)
    assert result.total_diffs == pytest.approx(0.0)
    assert result.total_comps == pytest.approx(0.0)
    assert result.avg_pi == "NA"


def test_calc_pi_gl_rejects_bad_shape() -> None:
    """Posterior array that isn't (n_sites, n_samples, 3) is a usage error."""
    with pytest.raises(ValueError, match="shape"):
        calc_pi_gl(np.zeros((1, 1, 2), dtype=np.float64), np.zeros((1, 1), dtype=bool), 2)


def test_calc_pi_gl_rejects_mismatched_mask() -> None:
    """Missing-mask whose shape doesn't match the posterior array is a usage error."""
    with pytest.raises(ValueError, match="missing_mask"):
        calc_pi_gl(
            np.zeros((2, 3, 3), dtype=np.float64),
            np.zeros((2, 2), dtype=bool),
            6,
        )


def test_reduces_to_hardcall_detects_one_hot() -> None:
    """`reduces_to_hardcall` distinguishes one-hot rows from mixed posteriors per site."""
    one_hot = np.array([
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
    ])
    not_one_hot = np.array([
        [[0.5, 0.5, 0.0], [1.0, 0.0, 0.0]],
    ])
    np.testing.assert_array_equal(reduces_to_hardcall(one_hot), np.array([True, True]))
    np.testing.assert_array_equal(reduces_to_hardcall(not_one_hot), np.array([False]))
