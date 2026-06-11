r"""
Wisp mask support for pixy: variants-only VCF + a quantized callable-sites BED.

Background
----------
Computing pi, dxy, Watterson's theta, or Tajima's D from a VCF requires knowing,
at every site in each window, how many samples per population had an observable
genotype — the "callable" denominator. The canonical way to supply this is an
all-sites VCF (one record per reference position, with `./.` GTs where coverage
fails). All-sites VCFs are large and expensive to call.

The wisp companion tool (https://github.com/samuk-lab/wisp) produces a much smaller
representation: a BED with one column per population whose value is the number of
samples in that population with an "intact" genotype across every site in the
interval. Adjacent intervals with identical counts across all populations are
collapsed into a single row.

When a wisp BED is supplied alongside a variants-only VCF, pixy reads variant
sites from the VCF as before and analytically adds the contribution of invariant
sites to each window's denominator, deriving the per-site sample counts from the
wisp mask. The variant-site arithmetic is unchanged.

Format
------
The first line of the BED carries a JSON header preceded by
``#wisp_mask_metadata\\t``::

    #wisp_mask_metadata\\t{"populations": ["GBR", "YRI"],
                            "population_sample_counts": {"GBR": 10, "YRI": 10},
                            "population_columns": [
                              {"column_number": 4, "name": "GBR"},
                              {"column_number": 5, "name": "YRI"}],
                            "threshold": 30, ...}

Subsequent rows are tab-delimited ``chrom\\tstart\\tend\\tk1\\tk2\\t...`` in
BED 0-based half-open coordinates.
"""

from __future__ import annotations

import gzip
import json
import subprocess
from dataclasses import dataclass
from dataclasses import field
from pathlib import Path
from typing import Dict
from typing import Iterator
from typing import List
from typing import Optional
from typing import Tuple


@dataclass(frozen=True)
class WispMetadata:
    """
    Metadata parsed from the ``#wisp_mask_metadata`` JSON header.

    Attributes:
        populations: per-population labels in their BED-column order.
        population_sample_counts: total sample count (``N``) per population — the maximum
            possible value any data cell can take.
        threshold: depth (or whatever criterion wisp used) below which a sample is
            considered to lack an intact genotype. Stored for diagnostics only.
        column_index_for: mapping from population name to the 0-based index into the data
            row's split fields where that population's count lives. (Wisp's
            ``column_number`` is 1-based; we convert.)
    """

    populations: Tuple[str, ...]
    population_sample_counts: Dict[str, int]
    threshold: int
    column_index_for: Dict[str, int]

    @classmethod
    def from_header_line(cls, line: str) -> "WispMetadata":
        r"""
        Parse the first ``#wisp_mask_metadata`` line of a wisp BED.

        The expected shape is ``#wisp_mask_metadata\t{JSON}\n``.
        """
        if not line.startswith("#wisp_mask_metadata\t"):
            raise ValueError("Wisp BED is missing the required `#wisp_mask_metadata` header line")
        _, json_text = line.rstrip("\n").split("\t", 1)
        try:
            meta = json.loads(json_text)
        except json.JSONDecodeError as e:
            raise ValueError(f"Wisp metadata header is not valid JSON: {e.msg}") from e

        for key in ("populations", "population_sample_counts", "population_columns"):
            if key not in meta:
                raise ValueError(f"Wisp metadata is missing required key {key!r}")

        populations = tuple(str(p) for p in meta["populations"])
        sample_counts = {str(k): int(v) for k, v in meta["population_sample_counts"].items()}
        # `population_columns` is `[{"column_number": int, "name": str}, ...]` — convert to
        # a population -> 0-based field index.
        col_for: Dict[str, int] = {}
        for entry in meta["population_columns"]:
            name = str(entry["name"])
            col = int(entry["column_number"]) - 1
            if col < 3:
                raise ValueError(
                    f"Wisp metadata has population {name!r} mapped to column {col + 1}; "
                    "population columns must come after chrom/start/end (column >= 4)"
                )
            col_for[name] = col
        # Cross-check: every name in `populations` should have a column mapping.
        missing = [p for p in populations if p not in col_for]
        if missing:
            raise ValueError(
                f"Wisp metadata declares populations {missing!r} but no column mapping"
            )
        # And every name should have a sample count.
        missing_counts = [p for p in populations if p not in sample_counts]
        if missing_counts:
            raise ValueError(
                f"Wisp metadata declares populations {missing_counts!r} but no sample counts"
            )
        threshold = int(meta.get("threshold", 0))
        return cls(
            populations=populations,
            population_sample_counts=sample_counts,
            threshold=threshold,
            column_index_for=col_for,
        )


@dataclass(frozen=True)
class WispRange:
    """
    One quantized callable-sites interval from a wisp BED.

    Coordinates are 0-based half-open (BED convention). ``counts`` is parallel to
    ``WispMetadata.populations``.
    """

    start: int
    end: int
    counts: Tuple[int, ...]

    @property
    def length(self) -> int:
        """Number of sites in this interval."""
        return self.end - self.start


@dataclass
class WispMask:
    """
    Tabix-backed reader over a wisp BED.

    Reading is done one window at a time via ``iter_range``; the file itself is never
    fully loaded into memory.
    """

    path: Path
    metadata: WispMetadata
    # In-memory cache of populations we know how to look up — populated lazily.
    _pop_order_cache: Dict[Tuple[str, ...], Tuple[int, ...]] = field(default_factory=dict)

    @classmethod
    def from_path(cls, path: Path) -> "WispMask":
        """
        Open a wisp BED and parse its metadata header.

        The companion ``.tbi`` index is required for window queries; we don't verify
        that here so that the calling validator can produce a pixy-style error message.
        """
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as fh:
            first = fh.readline()
        metadata = WispMetadata.from_header_line(first)
        return cls(path=path, metadata=metadata)

    def column_indices_for(self, populations: Tuple[str, ...]) -> Tuple[int, ...]:
        """Return the 0-based field indices for ``populations``, in the given order."""
        if populations in self._pop_order_cache:
            return self._pop_order_cache[populations]
        try:
            indices = tuple(self.metadata.column_index_for[p] for p in populations)
        except KeyError as e:
            missing = [p for p in populations if p not in self.metadata.column_index_for]
            raise ValueError(
                f"Population(s) {missing!r} are not present in the wisp metadata "
                f"(known: {list(self.metadata.column_index_for)})"
            ) from e
        self._pop_order_cache[populations] = indices
        return indices

    def iter_range(
        self,
        chrom: str,
        start: int,
        end: int,
        column_indices: Tuple[int, ...],
    ) -> Iterator[WispRange]:
        """
        Yield wisp ranges overlapping the half-open interval ``[start, end)``.

        ``column_indices`` is the field-index sequence returned by
        ``column_indices_for`` for the populations you want to pull.

        Implementation note: tabix uses 1-based inclusive coordinates. We translate
        once here and let tabix do the spatial filtering. Returned ranges are NOT
        clipped to the query window — callers must clip if they need exact overlap
        arithmetic.
        """
        if end <= start:
            return
        region = f"{chrom}:{start + 1}-{end}"
        proc = subprocess.run(
            ["tabix", str(self.path), region],
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            # tabix prints to stderr on a missing/unindexed file. Surface that verbatim.
            raise RuntimeError(
                f"tabix lookup failed for wisp BED {self.path}: {proc.stderr.strip()}"
            )
        for raw_line in proc.stdout.splitlines():
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            # Defensive: a row with too few fields can't be parsed safely.
            max_col = max(column_indices) if column_indices else 2
            if len(fields) <= max_col:
                row_pos = fields[1] if len(fields) > 1 else "?"
                raise ValueError(
                    f"Wisp BED {self.path}: row at {chrom}:{row_pos} "
                    f"has {len(fields)} fields but column index {max_col} was requested"
                )
            counts = tuple(int(fields[c]) for c in column_indices)
            yield WispRange(
                start=int(fields[1]),
                end=int(fields[2]),
                counts=counts,
            )


# ---------------------------------------------------------------------------
# Per-window invariant contributions
# ---------------------------------------------------------------------------


@dataclass
class PiInvariantContribution:
    """
    Per-population invariant-site additions to a pi result.

    All sums are over invariant sites (= sites in the window that are NOT in the
    variants-only VCF) in a single window for a single population.

    Attributes:
        num_sites: count of invariant sites where at least one haploid was observed
            (k > 0). Goes into the ``no_sites`` output column.
        comps: sum of pairwise comparisons k_haps * (k_haps - 1) / 2 across invariant
            sites. Added to the variant-derived ``total_comparisons``.
        missing: sum of (N_haps * (N_haps - 1) / 2 - k_haps * (k_haps - 1) / 2). Added
            to the variant-derived ``total_missing``.
    """

    num_sites: int = 0
    comps: int = 0
    missing: int = 0


@dataclass
class DxyInvariantContribution:
    """
    Per-population-pair invariant-site additions to a dxy result.

    Attributes:
        num_sites: count of invariant sites where both populations had at least one
            observed haploid (k1 > 0 and k2 > 0).
        comps: sum of k1_haps * k2_haps across invariant sites.
        missing: sum of (N1_haps * N2_haps - k1_haps * k2_haps).
    """

    num_sites: int = 0
    comps: int = 0
    missing: int = 0


@dataclass
class WattersonInvariantContribution:
    """
    Per-population invariant-site additions to a Watterson's theta result.

    Invariants don't contribute to the numerator (``raw_theta``) or to the variant-site
    classes — they're all reference. They DO contribute to ``num_sites`` (the
    denominator) and to ``weighted_sites`` (an emitted diagnostic that weights each
    site by its observed-haploid count divided by the window-wide maximum).

    Attributes:
        num_sites: count of invariant sites with k_haps > 0.
        k_haps_counter: distribution of observed-haploid counts across invariant sites
            (key = k_haps, value = number of invariant sites with that count). Combined
            with the variant-site distribution to recompute ``weighted_sites`` once the
            window-wide ``max(k_haps)`` is known.
    """

    num_sites: int = 0
    k_haps_counter: Dict[int, int] = field(default_factory=dict)

    def add(self, k_haps: int, n_sites: int) -> None:
        """Record ``n_sites`` invariant sites with ``k_haps`` observed haploids."""
        if k_haps > 0 and n_sites > 0:
            self.num_sites += n_sites
            self.k_haps_counter[k_haps] = self.k_haps_counter.get(k_haps, 0) + n_sites


@dataclass
class TajimaInvariantContribution:
    """
    Per-population invariant-site additions to a Tajima's D result.

    Invariants only nudge ``num_sites`` upward — the numerator (``raw_pi``),
    Watterson's theta term, and the per-site standard-deviation accumulator all derive
    from variant sites only.
    """

    num_sites: int = 0


@dataclass
class WindowInvariantContributions:
    """
    Container for all per-window invariant contributions.

    Built once per window by ``compute_window_invariant_contributions`` and then
    threaded through ``compute_summary_pi`` / ``compute_summary_dxy`` and the inline
    Watterson / Tajima loops in ``compute_summary_stats``.
    """

    pi: Dict[str, PiInvariantContribution] = field(default_factory=dict)
    dxy: Dict[Tuple[str, str], DxyInvariantContribution] = field(default_factory=dict)
    watterson: Dict[str, WattersonInvariantContribution] = field(default_factory=dict)
    tajima: Dict[str, TajimaInvariantContribution] = field(default_factory=dict)


def compute_window_invariant_contributions(  # noqa: C901
    wisp_mask: WispMask,
    chromosome: str,
    window_pos_1: int,
    window_pos_2: int,
    variant_positions: List[int],
    pop_names: List[str],
    ploidy: int,
) -> WindowInvariantContributions:
    """
    Compute analytical invariant-site contributions for one window from a wisp mask.

    Args:
        wisp_mask: the open wisp mask.
        chromosome: chromosome name.
        window_pos_1: 1-based inclusive window start (pixy's convention).
        window_pos_2: 1-based inclusive window end.
        variant_positions: 1-based positions of variant sites already loaded from the
            VCF for this window. Used to compute ``n_invariant_in_range = range_length -
            n_variant_in_range`` per wisp interval. Must contain only positions inside
            the window.
        pop_names: populations to compute contributions for. Must all be present in
            the wisp metadata.
        ploidy: per-contig ploidy (typically 2; pixy supports per-contig variation).

    Returns:
        A ``WindowInvariantContributions`` ready to merge into per-population results.
    """
    # Window in 0-based half-open coordinates for clipping wisp ranges.
    win_start = window_pos_1 - 1
    win_end = window_pos_2  # window_pos_2 is 1-based inclusive ⇒ half-open end == win_pos_2

    contributions = WindowInvariantContributions()
    pop_names_t = tuple(pop_names)
    col_indices = wisp_mask.column_indices_for(pop_names_t)
    sample_counts = [wisp_mask.metadata.population_sample_counts[p] for p in pop_names]

    # Pre-init per-pop / per-pair accumulators so missing populations still get an
    # empty contribution (rather than a KeyError later in the consumer).
    for pop in pop_names:
        contributions.pi[pop] = PiInvariantContribution()
        contributions.watterson[pop] = WattersonInvariantContribution()
        contributions.tajima[pop] = TajimaInvariantContribution()
    for i in range(len(pop_names)):
        for j in range(i + 1, len(pop_names)):
            contributions.dxy[(pop_names[i], pop_names[j])] = DxyInvariantContribution()

    # Sort the variant positions once; we use binary search per wisp range to count
    # how many variants fall inside it (those don't count as invariant sites).
    sorted_var_pos = sorted(int(p) for p in variant_positions)

    for rng in wisp_mask.iter_range(chromosome, win_start, win_end, col_indices):
        # Clip the wisp range to the window in 0-based half-open coordinates.
        eff_start_0 = max(rng.start, win_start)
        eff_end_0 = min(rng.end, win_end)
        if eff_end_0 <= eff_start_0:
            continue
        range_length = eff_end_0 - eff_start_0

        # Count variants in this clipped range — pixy uses 1-based positions, so
        # the half-open 0-based interval [eff_start_0, eff_end_0) becomes the
        # closed 1-based interval [eff_start_0 + 1, eff_end_0].
        var_lo = eff_start_0 + 1
        var_hi = eff_end_0
        # bisect: number of var positions in [var_lo, var_hi]
        lo_idx = _bisect_left(sorted_var_pos, var_lo)
        hi_idx = _bisect_right(sorted_var_pos, var_hi)
        n_variant_in_range = hi_idx - lo_idx
        n_invariant_in_range = range_length - n_variant_in_range
        if n_invariant_in_range <= 0:
            # Either the range contains only variant sites, or (defensively) the VCF
            # claims more variants than the range has sites. Skip either way; the
            # variant contribution is handled by the VCF code path.
            continue

        # Per-population contributions.
        per_pop_k_haps: List[int] = []
        per_pop_total_haps: List[int] = []
        for pop, k_pop, n_pop_samples in zip(pop_names, rng.counts, sample_counts, strict=True):
            k_haps = int(k_pop) * ploidy
            n_total_haps = int(n_pop_samples) * ploidy
            per_pop_k_haps.append(k_haps)
            per_pop_total_haps.append(n_total_haps)

            site_comps = k_haps * (k_haps - 1) // 2
            site_total = n_total_haps * (n_total_haps - 1) // 2
            pi_inv = contributions.pi[pop]
            pi_inv.comps += n_invariant_in_range * site_comps
            pi_inv.missing += n_invariant_in_range * (site_total - site_comps)
            if k_haps > 0:
                pi_inv.num_sites += n_invariant_in_range

            contributions.watterson[pop].add(k_haps=k_haps, n_sites=n_invariant_in_range)
            if k_haps > 0:
                contributions.tajima[pop].num_sites += n_invariant_in_range

        # Per-pair dxy contributions.
        for i in range(len(pop_names)):
            for j in range(i + 1, len(pop_names)):
                k1 = per_pop_k_haps[i]
                k2 = per_pop_k_haps[j]
                n1_total = per_pop_total_haps[i]
                n2_total = per_pop_total_haps[j]
                site_comps_pair = k1 * k2
                site_total_pair = n1_total * n2_total
                dxy_inv = contributions.dxy[(pop_names[i], pop_names[j])]
                dxy_inv.comps += n_invariant_in_range * site_comps_pair
                dxy_inv.missing += n_invariant_in_range * (site_total_pair - site_comps_pair)
                if k1 > 0 and k2 > 0:
                    dxy_inv.num_sites += n_invariant_in_range

    return contributions


# Tiny local bisect helpers — pixy avoids stdlib `bisect` import here for no good
# reason other than keeping the import surface tight. They're identical in behavior.


def _bisect_left(a: List[int], x: int) -> int:
    lo, hi = 0, len(a)
    while lo < hi:
        mid = (lo + hi) // 2
        if a[mid] < x:
            lo = mid + 1
        else:
            hi = mid
    return lo


def _bisect_right(a: List[int], x: int) -> int:
    lo, hi = 0, len(a)
    while lo < hi:
        mid = (lo + hi) // 2
        if a[mid] <= x:
            lo = mid + 1
        else:
            hi = mid
    return lo


def empty_contributions(pop_names: List[str]) -> WindowInvariantContributions:
    """
    Return a zeroed ``WindowInvariantContributions`` for the given populations.

    Used when the wisp path is enabled but a particular window has no wisp ranges
    overlapping it (the contributions are zero but the structure must be present so
    the consumer can look up populations without KeyError).
    """
    contributions = WindowInvariantContributions()
    for pop in pop_names:
        contributions.pi[pop] = PiInvariantContribution()
        contributions.watterson[pop] = WattersonInvariantContribution()
        contributions.tajima[pop] = TajimaInvariantContribution()
    for i in range(len(pop_names)):
        for j in range(i + 1, len(pop_names)):
            contributions.dxy[(pop_names[i], pop_names[j])] = DxyInvariantContribution()
    return contributions


def validate_wisp_data_rows(  # noqa: C901
    wisp_mask: WispMask, peek_rows: int = 5
) -> Optional[str]:
    """
    Sanity-check the wisp BED's data rows.

    Reads up to ``peek_rows`` non-header rows after the metadata header and confirms
    each parses as a well-formed wisp entry:

      * the row has at least ``max(column_number) + 1`` tab-delimited fields,
      * the start/end coordinates are integers with ``end > start``,
      * each per-population count is an integer in ``[0, N_pop]``.

    Returns ``None`` on success, or a human-readable error message on the first
    problem encountered. Errors include a header-only file (no data rows). The
    populations-vs-popfile cross-check is handled separately by
    ``validate_wisp_against_populations``.

    Cheap by design — only the first ``peek_rows`` data rows are inspected, so the
    cost is independent of mask size.
    """
    path = wisp_mask.path
    populations = wisp_mask.metadata.populations
    if populations:
        max_col = max(wisp_mask.metadata.column_index_for[p] for p in populations)
    else:
        max_col = 2

    data_rows_seen = 0
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line_idx, line in enumerate(fh, start=1):
            stripped = line.rstrip("\n")
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split("\t")
            if len(fields) < max_col + 1:
                return (
                    f"Wisp mask {path}: row at line {line_idx} has {len(fields)} "
                    f"tab-delimited fields; the metadata header declared population columns "
                    f"up to index {max_col} (expected at least {max_col + 1} fields)"
                )
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                return (
                    f"Wisp mask {path}: row at line {line_idx} has non-integer "
                    f"coordinates start={fields[1]!r}, end={fields[2]!r}"
                )
            if end <= start:
                return (
                    f"Wisp mask {path}: row at line {line_idx} has end <= start "
                    f"({fields[0]}:{start}-{end})"
                )
            for pop in populations:
                idx = wisp_mask.metadata.column_index_for[pop]
                try:
                    count = int(fields[idx])
                except ValueError:
                    return (
                        f"Wisp mask {path}: row at line {line_idx} has non-integer "
                        f"count for population {pop!r}: {fields[idx]!r}"
                    )
                n_total = wisp_mask.metadata.population_sample_counts[pop]
                if count < 0 or count > n_total:
                    return (
                        f"Wisp mask {path}: row at line {line_idx} has count={count} "
                        f"for population {pop!r} (expected 0..{n_total})"
                    )
            data_rows_seen += 1
            if data_rows_seen >= peek_rows:
                return None

    if data_rows_seen == 0:
        return f"Wisp mask {path}: no data rows found (only the metadata header)"
    return None


def validate_wisp_against_populations(
    wisp_mask: WispMask,
    pop_to_sample_count: Dict[str, int],
) -> Optional[str]:
    """
    Confirm that the wisp mask matches the user's populations file.

    Returns a human-readable error message if there's a mismatch, or ``None`` if
    everything checks out. Mismatches detected:
      * a population in the wisp metadata that's absent from --populations
      * a population in --populations that's absent from the wisp metadata
      * a per-population sample count that doesn't match between the two

    Args:
        wisp_mask: parsed wisp mask.
        pop_to_sample_count: mapping from population name to the number of samples
            in that population, as derived from the user's --populations file.
    """
    wisp_pops = set(wisp_mask.metadata.populations)
    user_pops = set(pop_to_sample_count)

    only_in_wisp = wisp_pops - user_pops
    only_in_users = user_pops - wisp_pops
    if only_in_wisp or only_in_users:
        problems = []
        if only_in_users:
            problems.append(f"only in populations file: {sorted(only_in_users)}")
        if only_in_wisp:
            problems.append(f"only in wisp mask: {sorted(only_in_wisp)}")
        return "Wisp mask populations do not match --populations file. " + "; ".join(problems)

    count_mismatches = []
    for pop in sorted(user_pops):
        wisp_n = wisp_mask.metadata.population_sample_counts.get(pop)
        user_n = pop_to_sample_count[pop]
        if wisp_n != user_n:
            count_mismatches.append(f"{pop}: wisp={wisp_n}, populations_file={user_n}")
    if count_mismatches:
        return (
            "Wisp mask per-population sample counts differ from the populations file: "
            + ", ".join(count_mismatches)
        )
    return None
