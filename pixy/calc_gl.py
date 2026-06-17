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

By default the estimator now uses an empirical Hardy-Weinberg prior with the alt-allele
frequency `p` estimated per site via EM directly on the genotype likelihoods (analogous to
ANGSD's `-doMaf 1 -doPost 2`). `estimate_maf_em` runs the EM; `calc_pi_gl_em` is the
top-level entry point that ties EM + posterior + per-site math together. The flat prior
remains available via `calc_pi_gl` for the reduction property tests (under a flat prior
confident PLs collapse exactly to the hard-call estimator; under HWE that reduction only
holds at sites where `p_hat` saturates to 0 or 1).
"""

from typing import Literal
from typing import Optional
from typing import Tuple
from typing import cast

import numpy as np
from numpy.typing import NDArray

from pixy.models import PiResult


def _pl_or_gl_to_linear(
    lik: NDArray, lik_kind: str
) -> Tuple[NDArray[np.float64], NDArray[np.bool_]]:
    """
    Convert PL or GL to linear-scale likelihoods plus a missing-row mask.

    Shared by `likelihoods_to_posteriors` and `estimate_maf_em` so the numerical conventions
    (PL cap, GL max-subtraction, missing sentinels) only live in one place.

    Args:
        lik: array of shape (n_sites, n_samples, 3) of phred-scaled PL (int) or log10 GL.
        lik_kind: "PL" or "GL".

    Returns:
        lik_lin: float64 array of shape (n_sites, n_samples, 3) of linear-scale likelihoods.
            For GL the result is up to a positive per-row constant (the max-subtracted
            normalization); that constant drops out of both posterior normalization and EM.
        missing_mask: bool array of shape (n_sites, n_samples) where True marks a missing row.

    Missing convention:
      * PL: row of all `-1` (scikit-allel's missing sentinel) or all-zero is treated as missing.
      * GL: row containing any NaN or all-equal is treated as missing.
      * Either: a row whose linear likelihoods sum to exactly 0 (degenerate numerics).
    """
    if lik_kind not in ("PL", "GL"):
        raise ValueError(f"lik_kind must be 'PL' or 'GL', got {lik_kind!r}")

    arr = np.asarray(lik)
    if arr.ndim != 3 or arr.shape[-1] != 3:
        raise ValueError(f"lik must have shape (n_sites, n_samples, 3); got {arr.shape}")

    if lik_kind == "PL":
        missing_mask = np.logical_or(
            (arr == -1).all(axis=-1),
            (arr == 0).all(axis=-1),
        )
        # PL = -10 * log10(L); L = 10**(-PL/10). Cap at PL=255 — values that large are
        # effectively zero anyway.
        pl = np.clip(arr.astype(np.float64), 0.0, 255.0)
        lik_lin = np.power(10.0, -pl / 10.0)
    else:
        gl = arr.astype(np.float64)
        nan_row = np.isnan(gl).any(axis=-1)
        equal_row = np.logical_and(gl[..., 0] == gl[..., 1], gl[..., 1] == gl[..., 2])
        missing_mask = np.logical_or(nan_row, equal_row)
        # Max-subtract for numerical stability — the dropped per-row constant is positive
        # and cancels everywhere downstream (posterior normalization, EM ratio).
        gl_safe = np.where(np.isnan(gl), 0.0, gl)
        gl_shifted = gl_safe - gl_safe.max(axis=-1, keepdims=True)
        lik_lin = np.power(10.0, gl_shifted)

    row_sum = lik_lin.sum(axis=-1)
    missing_mask = np.logical_or(missing_mask, row_sum == 0.0)
    # Zero out missing rows so callers can sum without masking in the hot loop.
    lik_lin = np.where(missing_mask[..., None], 0.0, lik_lin)
    return lik_lin.astype(np.float64), missing_mask


def likelihoods_to_posteriors(
    lik: NDArray,
    lik_kind: str,
    prior: Optional[NDArray[np.float64]] = None,
) -> Tuple[NDArray[np.float64], NDArray[np.bool_]]:
    """
    Convert per-sample genotype likelihoods into normalized posteriors under a given prior.

    Args:
        lik: array of shape (n_sites, n_samples, 3) of phred-scaled PL (int) or log10 GL (float).
        lik_kind: "PL" or "GL".
        prior: optional float64 array of shape (n_sites, 3) giving a per-site prior over the
            three diploid biallelic genotypes. If None, a flat (1/3, 1/3, 1/3) prior is used —
            equivalent to the v1 behavior, since a flat prior cancels in the normalization.

    Returns:
        posteriors: float64 array of shape (n_sites, n_samples, 3) with each non-missing row
            summing to 1.0. Missing rows are filled with zeros and flagged in `missing_mask`.
        missing_mask: bool array of shape (n_sites, n_samples) where True marks a missing row.

    Missing convention: see `_pl_or_gl_to_linear`.
    """
    lik_lin, missing_mask = _pl_or_gl_to_linear(lik, lik_kind)

    if prior is not None:
        prior = np.asarray(prior, dtype=np.float64)
        if prior.shape != (lik_lin.shape[0], 3):
            raise ValueError(
                f"prior must have shape (n_sites, 3) = ({lik_lin.shape[0]}, 3); got {prior.shape}"
            )
        # Per-site, per-genotype multiplier broadcast over the sample axis.
        weighted = lik_lin * prior[:, None, :]
    else:
        weighted = lik_lin

    row_sum = weighted.sum(axis=-1, keepdims=True)
    # A prior with a zero genotype entry could send a previously-informative row to zero;
    # flag those as missing too so downstream sums stay consistent.
    new_zero_rows = np.logical_and(~missing_mask, row_sum.squeeze(-1) == 0.0)
    missing_mask = np.logical_or(missing_mask, new_zero_rows)

    safe_sum = np.where(row_sum == 0.0, 1.0, row_sum)
    posteriors = weighted / safe_sum
    posteriors = np.where(missing_mask[..., None], 0.0, posteriors)

    return posteriors.astype(np.float64), missing_mask


def estimate_maf_em(
    lik_lin: NDArray[np.float64],
    missing_mask: NDArray[np.bool_],
    *,
    max_iter: int = 50,
    tol: float = 1e-6,
) -> NDArray[np.float64]:
    """
    Per-site EM estimate of the alt-allele frequency from genotype likelihoods.

    Analog of ANGSD's `-doMaf 1`: maximizes the marginal likelihood
    `prod_i sum_k L_ik * C(2,k) * p^k * (1-p)^(2-k)` per site under HWE. Vectorized over
    sites — each iteration is two einsums plus a normalize.

    Args:
        lik_lin: float64 array of shape (n_sites, n_samples, 3) of linear-scale likelihoods
            (output of `_pl_or_gl_to_linear`). Missing rows must already be zeroed.
        missing_mask: bool array of shape (n_sites, n_samples) where True marks a missing row.
        max_iter: max EM iterations (default 50).
        tol: convergence tolerance on max |p_new - p_old| across sites (default 1e-6).

    Returns:
        p_hat: float64 array of shape (n_sites,) clipped to `[1/(4n), 1 - 1/(4n)]` per site,
            where `n` is the non-missing sample count at that site (ANGSD-style MAF floor).
            Sites with zero non-missing samples get `p_hat = 0.5` (no information, but the
            returned array still has a defined value so downstream broadcasts work).
    """
    if lik_lin.ndim != 3 or lik_lin.shape[-1] != 3:
        raise ValueError(f"lik_lin must have shape (n_sites, n_samples, 3); got {lik_lin.shape}")
    if missing_mask.shape != lik_lin.shape[:2]:
        raise ValueError(
            "missing_mask shape must match lik_lin first two dims; "
            f"got {missing_mask.shape} vs {lik_lin.shape[:2]}"
        )

    n_per_site = (~missing_mask).sum(axis=-1).astype(np.float64)  # (n_sites,)
    has_data = n_per_site > 0
    # Safe divisor (never divide by zero); for sites with no data we fall back to p=0.5
    # at the end. The intermediate computation may produce nonsense for those sites but is
    # overwritten below.
    safe_n = np.where(has_data, n_per_site, 1.0)

    # --- Initialize p from the flat-prior expected dose ---
    # Under a flat prior the row-normalized posterior is just lik_lin / row_sum.
    row_sum = lik_lin.sum(axis=-1, keepdims=True)
    safe_row_sum = np.where(row_sum == 0.0, 1.0, row_sum)
    flat_post = lik_lin / safe_row_sum  # missing rows stay zero (lik_lin was zeroed)
    mu_flat = flat_post[..., 1] + 2.0 * flat_post[..., 2]
    p = (mu_flat.sum(axis=-1) / (2.0 * safe_n)).astype(np.float64)
    # Clip well inside (0, 1) so the first iteration's prior is non-degenerate.
    p = np.clip(p, 1e-6, 1.0 - 1e-6)

    # --- EM ---
    for _ in range(max_iter):
        prior_0 = (1.0 - p) ** 2
        prior_1 = 2.0 * p * (1.0 - p)
        prior_2 = p**2
        # Weighted likelihood per sample per genotype: L_ik * prior_k (broadcast on sample axis).
        w = lik_lin * np.stack([prior_0, prior_1, prior_2], axis=-1)[:, None, :]
        w_sum = w.sum(axis=-1, keepdims=True)
        safe_w_sum = np.where(w_sum == 0.0, 1.0, w_sum)
        post = w / safe_w_sum  # (n_sites, n_samples, 3); zero rows stay zero
        # Expected alt dose summed over non-missing samples; per-site mean over (2n).
        mu = post[..., 1] + 2.0 * post[..., 2]
        p_new = mu.sum(axis=-1) / (2.0 * safe_n)
        p_new = np.clip(p_new, 1e-6, 1.0 - 1e-6)

        if np.max(np.abs(p_new - p)) < tol:
            p = p_new
            break
        p = p_new

    # ANGSD-style per-site floor: 1/(4n). Sites with no data get 0.5 (uninformative).
    floor = np.where(has_data, 1.0 / (4.0 * safe_n), 0.5)
    p = np.clip(p, floor, 1.0 - floor)
    p = np.where(has_data, p, 0.5)
    return cast(NDArray[np.float64], p.astype(np.float64))


def _pi_from_posteriors(
    posteriors: NDArray[np.float64],
    missing_mask: NDArray[np.bool_],
    n_haps: int,
) -> PiResult:
    """
    Per-site vectorized pi from posteriors. Shared by `calc_pi_gl` and `calc_pi_gl_em`.

    See the module docstring for the derivation of `diffs_site` and `comps_site`.
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

    mu = posteriors[..., 1] + 2.0 * posteriors[..., 2]
    h = posteriors[..., 1]

    s_sum = mu.sum(axis=-1)
    h_sum = h.sum(axis=-1)
    q_sum = np.einsum("ij,ij->i", mu, mu)

    n = (~missing_mask).sum(axis=-1).astype(np.float64)

    diffs_site = h_sum + 2.0 * (n - 1.0) * s_sum - s_sum * s_sum + q_sum
    comps_site = n * (2.0 * n - 1.0)
    diffs_site = np.maximum(diffs_site, 0.0)
    comps_site = np.maximum(comps_site, 0.0)

    total_diffs = float(diffs_site.sum())
    total_comps = float(comps_site.sum())

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


def calc_pi_gl(
    posteriors: NDArray[np.float64],
    missing_mask: NDArray[np.bool_],
    n_haps: int,
) -> PiResult:
    """
    Per-site vectorized GL-based pi from pre-computed posteriors.

    Thin wrapper around `_pi_from_posteriors`. The flat-prior path remains available via
    this entry point so the reduction-to-hard-call property tests keep working; production
    code calls `calc_pi_gl_em`.

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
    return _pi_from_posteriors(posteriors, missing_mask, n_haps)


def calc_pi_gl_em(
    lik: NDArray,
    lik_kind: str,
    n_haps: int,
    *,
    max_iter: int = 50,
    tol: float = 1e-6,
) -> PiResult:
    """
    GL-based pi using a per-site EM-estimated empirical Hardy-Weinberg prior.

    Equivalent to ANGSD `-doMaf 1 -doPost 2` plumbed into pixy's per-window pi reduction.

    Args:
        lik: array of shape (n_sites, n_samples, 3) of PL (int) or GL (float) values.
        lik_kind: "PL" or "GL".
        n_haps: total number of haploid samples in the population (= n_samples * 2 for diploid).
        max_iter: max EM iterations (default 50).
        tol: EM convergence tolerance on |Δp| (default 1e-6).

    Returns:
        PiResult with float `total_diffs`, `total_comps`, `total_missing` for the window.

    Notes:
        At invariant sites the synthesized confident hom-ref PLs drive `p_hat` to the
        per-site floor `1/(4n)`. The HWE prior puts essentially all mass on hom-ref so the
        posterior collapses to one-hot and the site contributes zero diffs / full comps —
        identical to the hard-call denominator semantics pixy already implements.
    """
    lik_lin, missing_mask = _pl_or_gl_to_linear(lik, lik_kind)
    p_hat = estimate_maf_em(lik_lin, missing_mask, max_iter=max_iter, tol=tol)

    # Build the per-site HWE prior triplet: ((1-p)^2, 2p(1-p), p^2).
    prior_0 = (1.0 - p_hat) ** 2
    prior_1 = 2.0 * p_hat * (1.0 - p_hat)
    prior_2 = p_hat**2
    hwe_prior = np.stack([prior_0, prior_1, prior_2], axis=-1)  # (n_sites, 3)

    # E-step on the final p_hat: posteriors = L * prior, row-normalized per sample.
    weighted = lik_lin * hwe_prior[:, None, :]
    row_sum = weighted.sum(axis=-1, keepdims=True)
    new_zero_rows = np.logical_and(~missing_mask, row_sum.squeeze(-1) == 0.0)
    missing_mask = np.logical_or(missing_mask, new_zero_rows)
    safe_sum = np.where(row_sum == 0.0, 1.0, row_sum)
    posteriors = weighted / safe_sum
    posteriors = np.where(missing_mask[..., None], 0.0, posteriors).astype(np.float64)

    return _pi_from_posteriors(posteriors, missing_mask, n_haps)


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
