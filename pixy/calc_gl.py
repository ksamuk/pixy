"""
Genotype-likelihood-based estimators (v1: pi for diploid biallelic sites).

The hard-call estimator in `pixy.calc.calc_pi` works on per-site allele counts produced by
scikit-allel's `count_alleles()`. At low or uneven sequencing coverage, hard calls bias pi
downward because heterozygotes are preferentially missed (the variant caller falls back to
hom-ref when one of the two reads happens to support the reference). The genotype-likelihood
(GL) approach used by ANGSD avoids that bias by keeping the per-genotype uncertainty and
integrating over it.

This module implements the GL analog of `calc_pi`, restricted to diploid biallelic sites
(matching ANGSD's scope for the same estimator). Multi-allelic and non-diploid sites are
rejected upstream in `pixy.args_validation`.

# Math

Let `n` = number of non-missing individuals at a site, and let `P(g_i = k)` for `k in {0,1,2}`
be the posterior probability that individual `i` has genotype `k` alt alleles (derived from
PL or GL via `likelihoods_to_posteriors`). Define:

    mu_i = P(g_i = 1) + 2 * P(g_i = 2)       (expected alt-allele dose)
    h_i  = P(g_i = 1)                        (heterozygote posterior)
    S    = sum_i mu_i
    H    = sum_i h_i
    Q    = sum_i mu_i**2

Per-site numerator (expected number of pairs of distinct chromosomes that differ):

    diffs_site = H + 2 * (n - 1) * S - S**2 + Q

Derivation: within-individual contribution sums to `H` (one het-pair contribution per
heterozygote). Between-individual contribution over (i, j) with i < j is
4 * [(mu_i/2)(1 - mu_j/2) + (1 - mu_i/2)(mu_j/2)], which expands to
2*(n-1)*S - S**2 + Q after collapsing the double sum via the identity
sum_{i<j} mu_i*mu_j = (S**2 - Q) / 2.

Per-site denominator (number of pairs of distinct chromosomes):

    comps_site = 2*n * (2*n - 1) / 2

Same form as the hard-call path. Reduces to `n_gts*(n_gts-1)/2` from
`pixy.calc._count_diff_comp_missing_vectorized` when the posterior collapses to one-hot
(confident calls), in which case `mu_i in {0,1,2}` and `h_i = 1` iff het, and the formula
reduces algebraically to `(N**2 - sum ac**2) / 2`, matching the hard-call numerator exactly.
This reduction is enforced by tests in `tests/test_calc_gl.py`.

Missing handling: a sample with missing PL/GL at a site contributes nothing to `S`, `H`, `Q`,
and reduces `n` by 1 at that site. `total_missing` mirrors pixy's existing semantics
(possible_comps - actual_comps) using the population-wide `n_haps` as `possible`.

# Prior

v1 uses a flat prior over the three genotypes. ANGSD's `-doPost 2` analog (empirical
allele-frequency-informed prior) is a planned future addition; it would change the conversion
in `likelihoods_to_posteriors` but not the downstream math.
"""

from typing import Literal
from typing import Tuple
from typing import cast

import numpy as np
from numpy.typing import NDArray

from pixy.models import PiResult


def likelihoods_to_posteriors(
    lik: NDArray, lik_kind: str
) -> Tuple[NDArray[np.float64], NDArray[np.bool_]]:
    """
    Convert per-sample genotype likelihoods into normalized posteriors under a flat prior.

    Args:
        lik: array of shape (n_sites, n_samples, 3) of phred-scaled PL (int) or log10 GL (float).
        lik_kind: "PL" or "GL".

    Returns:
        posteriors: float64 array of shape (n_sites, n_samples, 3) with each non-missing row
            summing to 1.0. Missing rows are filled with zeros and flagged in `missing_mask`.
        missing_mask: bool array of shape (n_sites, n_samples) where True marks a missing row.

    Missing convention:
      * PL: row of all `-1` (scikit-allel's missing sentinel) or all-zero is treated as missing.
      * GL: row containing any NaN or all-equal is treated as missing.
    """
    if lik_kind not in ("PL", "GL"):
        raise ValueError(f"lik_kind must be 'PL' or 'GL', got {lik_kind!r}")

    arr = np.asarray(lik)
    if arr.ndim != 3 or arr.shape[-1] != 3:
        raise ValueError(f"lik must have shape (n_sites, n_samples, 3); got {arr.shape}")

    if lik_kind == "PL":
        # scikit-allel encodes missing PL as -1. All-zero triplets are also commonly emitted
        # at sites where the caller refused to compute likelihoods.
        missing_mask = np.logical_or(
            (arr == -1).all(axis=-1),
            (arr == 0).all(axis=-1),
        )
        # PL = -10 * log10(L); L = 10**(-PL/10). Cap at PL=255 (the int8 ceiling some callers
        # use) to keep the exponent bounded — values that large are effectively zero anyway.
        pl = np.clip(arr.astype(np.float64), 0.0, 255.0)
        lik_lin = np.power(10.0, -pl / 10.0)
    else:
        # GL = log10(L); L = 10**GL. NaN means missing for floating-point GL emissions.
        gl = arr.astype(np.float64)
        nan_row = np.isnan(gl).any(axis=-1)
        equal_row = np.logical_and(gl[..., 0] == gl[..., 1], gl[..., 1] == gl[..., 2])
        missing_mask = np.logical_or(nan_row, equal_row)
        # Numerically stable normalization: subtract the row max before exponentiating,
        # then divide by the row sum. This is mathematically equivalent to 10**GL / sum(10**GL)
        # but avoids underflow when all GLs are very negative.
        gl_safe = np.where(np.isnan(gl), 0.0, gl)
        gl_shifted = gl_safe - gl_safe.max(axis=-1, keepdims=True)
        lik_lin = np.power(10.0, gl_shifted)

    row_sum = lik_lin.sum(axis=-1, keepdims=True)
    # Guard against zero rows (extremely low-likelihood numerics) — treat as missing.
    zero_sum = row_sum.squeeze(-1) == 0.0
    missing_mask = np.logical_or(missing_mask, zero_sum)

    # Normalize with a safe divisor; rows flagged as missing get zeroed below.
    safe_sum = np.where(row_sum == 0.0, 1.0, row_sum)
    posteriors = lik_lin / safe_sum

    # Zero out missing rows so callers can sum without masking inside the hot loop.
    posteriors = np.where(missing_mask[..., None], 0.0, posteriors)

    return posteriors.astype(np.float64), missing_mask


def calc_pi_gl(
    posteriors: NDArray[np.float64],
    missing_mask: NDArray[np.bool_],
    n_haps: int,
) -> PiResult:
    """
    Per-site vectorized GL-based pi.

    Args:
        posteriors: float64 array of shape (n_sites, n_samples, 3) with rows summing to 1 for
            non-missing samples and zeros for missing ones (output of
            `likelihoods_to_posteriors`).
        missing_mask: bool array of shape (n_sites, n_samples) — True where the sample is
            missing at the site.
        n_haps: total number of haploid samples in the population
            (= n_samples * 2 for diploid; matches the convention in `pixy.calc._n_haps`).

    Returns:
        A PiResult with float `total_diffs` and `total_comps` (the GL estimator produces
        expected counts, which are real-valued).
    """
    if posteriors.ndim != 3 or posteriors.shape[-1] != 3:
        raise ValueError(
            f"posteriors must have shape (n_sites, n_samples, 3); got {posteriors.shape}"
        )
    if missing_mask.shape != posteriors.shape[:2]:
        raise ValueError(
            "missing_mask shape must match posteriors first two dims; "
            f"got {missing_mask.shape} vs {posteriors.shape[:2]}"
        )

    # Per-individual expected dose mu and heterozygote prob h.
    # Missing rows have posteriors = 0, so they contribute 0 to all sums below.
    mu = posteriors[..., 1] + 2.0 * posteriors[..., 2]  # (n_sites, n_samples)
    h = posteriors[..., 1]

    # Per-site reductions across samples. Names mirror the symbols in the module docstring
    # (uppercase S / H / Q) so the derivation lines up character-for-character; ruff's
    # N806 (uppercase locals) is suppressed for those three.
    s_sum = mu.sum(axis=-1)  # sum of mu_i  (symbol: S)
    h_sum = h.sum(axis=-1)  # sum of h_i   (symbol: H)
    q_sum = np.einsum("ij,ij->i", mu, mu)  # sum of mu_i**2 (symbol: Q)

    # Non-missing individual count per site.
    n = (~missing_mask).sum(axis=-1).astype(np.float64)

    # Per-site numerator and denominator (see module docstring for derivation).
    # diffs_site = H + 2*(n - 1)*S - S**2 + Q
    diffs_site = h_sum + 2.0 * (n - 1.0) * s_sum - s_sum * s_sum + q_sum
    # comps_site = 2n * (2n - 1) / 2
    comps_site = n * (2.0 * n - 1.0)
    # Floating-point arithmetic can leave tiny negative residuals near zero when every
    # posterior at a site is one-hot hom-ref; clip to zero for numerical hygiene.
    diffs_site = np.maximum(diffs_site, 0.0)
    comps_site = np.maximum(comps_site, 0.0)

    total_diffs = float(diffs_site.sum())
    total_comps = float(comps_site.sum())

    # `total_missing` mirrors pixy's hard-call convention: possible_comps - actual_comps,
    # where `possible_comps` is n_haps * (n_haps - 1) / 2 per site.
    n_sites = posteriors.shape[0]
    possible_per_site = n_haps * (n_haps - 1) / 2
    total_possible = float(possible_per_site * n_sites)
    total_missing = total_possible - total_comps

    avg_pi: float | Literal["NA"]
    if total_comps > 0:
        avg_pi = total_diffs / total_comps
    else:
        avg_pi = "NA"

    return PiResult(
        avg_pi=avg_pi,
        total_diffs=total_diffs,
        total_comps=total_comps,
        total_missing=total_missing,
    )


def reduces_to_hardcall(posteriors: NDArray[np.float64], tol: float = 1e-12) -> NDArray[np.bool_]:
    """
    Per-site mask: True where every non-zero row in `posteriors` is one-hot within `tol`.

    Used by the reduction unit test. A row of all zeros (missing) is also considered to
    satisfy the one-hot condition.
    """
    if posteriors.ndim != 3 or posteriors.shape[-1] != 3:
        raise ValueError(
            f"posteriors must have shape (n_sites, n_samples, 3); got {posteriors.shape}"
        )
    row_max = posteriors.max(axis=-1)
    row_sum = posteriors.sum(axis=-1)
    # One-hot: max == sum (the other two entries are zero). Empty rows (all zeros): also fine.
    one_hot = np.logical_or(np.abs(row_max - row_sum) <= tol, row_sum == 0.0)
    return cast(NDArray[np.bool_], one_hot.all(axis=-1))
