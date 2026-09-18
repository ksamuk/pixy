"""
GVCF block expansion for pixy.

Pixy's stat pipeline assumes one VCF record per genomic position. A GVCF
collapses runs of consecutive invariant positions into a single block record
with ``INFO/END=<pos>`` and an ``ALT`` of ``<NON_REF>`` (GATK) or ``<*>``
(Strelka). To reuse the existing pipeline unchanged, we expand each block
back into per-site rows immediately after ``allel.read_vcf`` returns, marking
the synthesised rows as invariant (``numalt=0``, ``is_snp=0``) so the
existing biallelic/invariant filter in :mod:`pixy.core` picks them up.

Memory footprint after expansion matches a true all-sites VCF over the same
region, so this approach has the same peak RSS as pixy's existing
all-sites path. A lower-memory analytic strategy (analogous to the wisp
mask) is plausible as a follow-up but is out of scope for v1.
"""

from __future__ import annotations

from typing import Dict
from typing import Optional

import numpy as np
from numpy.typing import NDArray

# Symbolic ALT alleles that stand for "any other allele" rather than an observed one.
SYMBOLIC_ALTS = ("<NON_REF>", "<*>")

# Number of ALT columns to request from scikit-allel for GVCF input: up to three real SNP
# alleles plus the trailing symbolic allele. (scikit-allel's default of 3 would truncate the
# symbolic allele off a triallelic-plus-`<NON_REF>` record and hide it from us.)
GVCF_ALT_NUMBER = 4


def strip_symbolic_alts(callset: Dict[str, NDArray]) -> Dict[str, NDArray]:
    """
    Discount symbolic ALT alleles (``<NON_REF>``, ``<*>``) from ``numalt`` and ``is_snp``.

    GATK writes ``<NON_REF>`` as the last ALT of *every* GVCF record, not just multi-site
    blocks. scikit-allel counts it as a real allele, so a reference-only record reports
    ``numalt=1, is_snp=False`` and a plain SNP (``A -> G,<NON_REF>``) reports
    ``numalt=2, is_snp=False``; pixy's site filter would discard both. This rewrites those
    two fields, for rows carrying a symbolic allele only, to what they would be had the
    symbolic allele not been listed. Genotype calls pointing at the symbolic allele (an
    unobserved, unknown allele) are set to missing.

    Args:
        callset: a dict as returned by :func:`allel.read_vcf` containing ``variants/REF``,
            ``variants/ALT``, ``variants/numalt``, ``variants/is_snp`` and ``calldata/GT``.

    Returns:
        The same dict, modified in place.
    """
    alt: NDArray = np.asarray(callset["variants/ALT"]).astype(str)
    if alt.ndim == 1:
        alt = alt[:, None]
    is_symbolic: NDArray = np.isin(alt, SYMBOLIC_ALTS)
    rows: NDArray = np.any(is_symbolic, axis=1)
    if not rows.any():
        return callset

    is_real: NDArray = (alt != "") & ~is_symbolic
    numalt_real: NDArray = is_real.sum(axis=1)
    ref_len: NDArray = np.char.str_len(np.asarray(callset["variants/REF"]).astype(str))
    # a SNP has a single-base REF and at least one real ALT, all of them single-base
    # (`*`, the spanning-deletion allele, is not a base)
    alt_is_base: NDArray = (np.char.str_len(alt) == 1) & (alt != "*")
    is_snp_real: NDArray = (ref_len == 1) & (numalt_real > 0) & (alt_is_base | ~is_real).all(axis=1)

    numalt: NDArray = callset["variants/numalt"]
    is_snp: NDArray = callset["variants/is_snp"]
    numalt[rows] = numalt_real[rows]
    is_snp[rows] = is_snp_real[rows]

    # allele index of the symbolic ALT (ALT column j is allele j + 1)
    gt: NDArray = callset["calldata/GT"]
    symbolic_allele: NDArray = np.where(rows, is_symbolic.argmax(axis=1) + 1, -2)
    gt[gt == symbolic_allele.reshape((-1,) + (1,) * (gt.ndim - 1))] = -1

    return callset


def expand_blocks(
    callset: Dict[str, NDArray],
    region_start: int,
    region_end: int,
) -> Dict[str, NDArray]:
    """
    Expand GVCF block rows in ``callset`` into per-site invariant rows.

    A block row is one whose ``variants/END`` is greater than ``variants/POS``.
    Each block of length ``L = END - POS + 1`` becomes ``L`` per-site rows whose
    genotype is the block's genotype replicated ``L`` times. The synthesised
    rows carry ``numalt=0`` and ``is_snp=0`` so the existing pipeline treats
    them as invariant sites. Non-block rows are passed through unchanged.

    The output is sorted by position and clipped to
    ``[region_start, region_end]`` (1-based inclusive). Clipping is needed
    because the caller widens the tabix region on the left by
    ``--gvcf_max_block_size`` so blocks that start before the requested
    window are not missed.

    If ``variants/END`` is absent from ``callset`` or no row has
    ``END > POS``, the input dict is returned unchanged (apart from the
    region clip, which the caller still wants).

    Args:
        callset: a dict as returned by :func:`allel.read_vcf` containing at
            least ``variants/POS``, ``calldata/GT``, ``variants/numalt`` and
            ``variants/is_snp``. Should also contain ``variants/END``; if
            absent, no expansion happens.
        region_start: 1-based inclusive start of the window of interest. Sites
            outside ``[region_start, region_end]`` are dropped.
        region_end: 1-based inclusive end of the window of interest.

    Returns:
        A dict with the same key set as the input (minus ``variants/END``),
        sorted by ``POS`` and clipped to the requested region.
    """
    pos: NDArray = callset["variants/POS"]
    end: Optional[NDArray] = callset.get("variants/END")

    # scikit-allel returns -1 for missing INFO ints; treat any END <= POS as
    # "not a block row" (including the all-missing case).
    if end is None:
        return _clip_to_region(callset, pos, region_start, region_end)

    block_mask = end > pos
    if not block_mask.any():
        return _clip_to_region(callset, pos, region_start, region_end)

    return _expand_and_assemble(callset, pos, end, block_mask, region_start, region_end)


def _expand_and_assemble(  # noqa: C901
    callset: Dict[str, NDArray],
    pos: NDArray,
    end: NDArray,
    block_mask: NDArray,
    region_start: int,
    region_end: int,
) -> Dict[str, NDArray]:
    """Vectorised block-expansion + clip. See :func:`expand_blocks` for semantics."""
    gt: NDArray = callset["calldata/GT"]
    numalt: NDArray = callset["variants/numalt"]
    is_snp: NDArray = callset["variants/is_snp"]
    chrom: Optional[NDArray] = callset.get("variants/CHROM")

    nonblock_mask = ~block_mask

    # --- non-block rows: kept as-is
    keep_pos = pos[nonblock_mask]
    keep_gt = gt[nonblock_mask]
    keep_numalt = numalt[nonblock_mask]
    keep_is_snp = is_snp[nonblock_mask]
    keep_chrom = chrom[nonblock_mask] if chrom is not None else None

    # --- block rows: vectorised expansion
    block_pos = pos[block_mask].astype(np.int64, copy=False)
    block_end = end[block_mask].astype(np.int64, copy=False)
    block_gt = gt[block_mask]
    lengths = (block_end - block_pos + 1).astype(np.int64, copy=False)
    total_len = int(lengths.sum())

    # Positions: for each expanded slot k (0..total_len-1) compute
    #   start = block_pos[ block(k) ]
    #   offset_within_block = k - cumulative_start[ block(k) ]
    # Both are derived via np.repeat across `lengths`, fully vectorised.
    starts_expanded = np.repeat(block_pos, lengths)
    cum_starts = np.cumsum(lengths) - lengths
    local_off = np.arange(total_len, dtype=np.int64) - np.repeat(cum_starts, lengths)
    expanded_pos = (starts_expanded + local_off).astype(pos.dtype, copy=False)

    # Genotypes: each block's GT row is replicated `lengths[i]` times.
    expanded_gt = np.repeat(block_gt, lengths, axis=0)

    expanded_numalt = np.zeros(total_len, dtype=numalt.dtype)
    expanded_is_snp = np.zeros(total_len, dtype=is_snp.dtype)

    # --- concatenate non-block + expanded and sort by POS
    new_pos = np.concatenate([keep_pos, expanded_pos])
    new_gt = np.concatenate([keep_gt, expanded_gt], axis=0)
    new_numalt = np.concatenate([keep_numalt, expanded_numalt])
    new_is_snp = np.concatenate([keep_is_snp, expanded_is_snp])

    order = np.argsort(new_pos, kind="stable")
    new_pos = new_pos[order]
    new_gt = new_gt[order]
    new_numalt = new_numalt[order]
    new_is_snp = new_is_snp[order]

    in_region = (new_pos >= region_start) & (new_pos <= region_end)
    new_pos = new_pos[in_region]
    new_gt = new_gt[in_region]
    new_numalt = new_numalt[in_region]
    new_is_snp = new_is_snp[in_region]

    out: Dict[str, NDArray] = {
        "variants/POS": new_pos,
        "calldata/GT": new_gt,
        "variants/numalt": new_numalt,
        "variants/is_snp": new_is_snp,
    }
    if chrom is not None:
        block_chrom = chrom[block_mask]
        expanded_chrom = np.repeat(block_chrom, lengths)
        assert keep_chrom is not None
        new_chrom = np.concatenate([keep_chrom, expanded_chrom])[order][in_region]
        out["variants/CHROM"] = new_chrom
    return out


def _clip_to_region(
    callset: Dict[str, NDArray],
    pos: NDArray,
    region_start: int,
    region_end: int,
) -> Dict[str, NDArray]:
    """
    Drop rows outside ``[region_start, region_end]`` from a callset.

    Used in the no-blocks-found branch of :func:`expand_blocks`. The callset
    may still contain rows from the left-extension that the caller added to
    catch overlapping blocks; even if no blocks were present in this region,
    those extra rows still need to be removed.
    """
    in_region = (pos >= region_start) & (pos <= region_end)
    if in_region.all():
        return callset
    return {
        key: arr[in_region] if isinstance(arr, np.ndarray) and arr.shape[:1] == pos.shape else arr
        for key, arr in callset.items()
    }
