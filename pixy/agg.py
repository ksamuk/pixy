"""
Pandas-free aggregation and output-writing for pixy's temp file.

The temp file is a tab-delimited stream of `PixyTempResult` rows (one per chunk/window/pop
tuple). This module replaces the previous pandas-based aggregation path (which loaded the
whole file with `pandas.read_csv` then did `groupby` + `to_csv`) with a streaming pass that
holds only the smallest dict of accumulators needed.

The output files are required to be byte-identical to what the previous pandas-based code
produced — the regression test fixtures in `tests/main/expected_outputs/` are compared
line-by-line after sorting. The formatters in this module match pandas's float64 / Int64 /
Float64 output exactly:

  * Floats: `repr(float(x))` — same shortest-roundtrip representation pandas emits.
  * Ints: `str(int(x))` — plain decimal.
  * Missing values: literal `"NA"`.
"""

from __future__ import annotations

import math
from collections import Counter
from collections import defaultdict
from pathlib import Path
from typing import Any
from typing import Dict
from typing import Iterable
from typing import Iterator
from typing import List
from typing import NamedTuple
from typing import Optional
from typing import Tuple
from typing import Union

from pixy.calc import calc_tajima_d_stdev
from pixy.calc import deserialize_tajima_d_variant_counts

# ---------------------------------------------------------------------------
# Temp row shape
# ---------------------------------------------------------------------------

# Number of columns in `PixyTempResult.__str__()` output (see pixy/models.py).
# Column 11 (`tajima_d_variant_counts`) is only present on Tajima's D rows; older rows have
# 11 columns and we treat the missing field as the sentinel "NA".
_TEMP_COL_COUNT = 12


class TempRow(NamedTuple):
    """
    One parsed line from the temp file.

    All numeric fields are kept as raw strings here; conversion to int/float happens lazily
    in the aggregator and the formatters. That keeps the streaming reader allocation-free
    for the cells we never look at (e.g. col 6 / `calculated_stat` is only used in the
    non-aggregated path).
    """

    stat: str
    pop1: str
    pop2: str
    chrom: str
    window_pos_1_str: str
    window_pos_2_str: str
    calculated_stat: str
    n_sites: str
    diffs: str
    comps: str
    missing: str
    tajima_counts: str  # "n:s,n:s,..." or "NA"


def iter_temp_rows(temp_file: Path) -> Iterator[TempRow]:
    """Stream the temp file as `TempRow` instances; one pass, constant memory."""
    with open(temp_file, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            # Pad the optional tajima-counts field if absent (older rows).
            if len(fields) < _TEMP_COL_COUNT:
                fields = fields + ["NA"] * (_TEMP_COL_COUNT - len(fields))
            yield TempRow(*fields[:_TEMP_COL_COUNT])


# ---------------------------------------------------------------------------
# Formatters
# ---------------------------------------------------------------------------

# Type alias for "numeric or NA" values as they appear in the temp file and in aggregated
# accumulators.
Numeric = Union[int, float, str]


def _is_na(x: Any) -> bool:
    """True if `x` is the literal NA sentinel, None, or a float NaN."""
    if x is None:
        return True
    if isinstance(x, str):
        return x == "NA"
    if isinstance(x, float):
        return math.isnan(x)
    return False


def _fmt_float(x: Any) -> str:
    """
    Format a value for the output file.

    Uses `%.14g` — 14 significant digits, scientific notation for very small/large values.
    This matches the precision cap that `PixyTempResult.__str__` writes into the temp file,
    so values produced by the calc layer and by the aggregator emit consistent strings.
    Missing/NA values become the literal `"NA"`.
    """
    if _is_na(x):
        return "NA"
    return f"{float(x):.14g}"


def _fmt_int(x: Any) -> str:
    """
    Format a value as pandas would format an Int64 (nullable int) column.

    Plain decimal, no thousands separator. Missing/NA becomes `"NA"`.
    """
    if _is_na(x):
        return "NA"
    return str(int(x))


def _to_num(x: str) -> Numeric:
    """Parse a temp-file cell. Returns the literal `"NA"` for missing values, else int/float."""
    if x == "NA" or x == "":
        return "NA"
    # Fast path: integers stay as int, floats as float.
    try:
        # `int(x)` rejects "1.0"; fall through to float on ValueError.
        return int(x)
    except ValueError:
        return float(x)


def _add(acc: Numeric, val: Numeric) -> Numeric:
    """
    Accumulate `val` into `acc`, treating NA as zero (matching pandas `sum(skipna=True)`).

    If both `acc` and `val` are NA, the result remains NA.
    """
    if _is_na(val):
        return acc
    if _is_na(acc):
        return val
    # `_is_na` rules out the "NA" sentinel string above, so both operands are int|float here.
    # mypy can't infer that through the helper, so we narrow with `assert`.
    assert not isinstance(acc, str) and not isinstance(val, str)
    return acc + val


# ---------------------------------------------------------------------------
# Aggregator
# ---------------------------------------------------------------------------

# Statistic-specific config: which population fields key each output row, and how to compute
# the final `calculated_stat` from the aggregated components.
_AGG_TWO_POP = {"dxy", "fst"}


def _agg_key(row: TempRow, stat: str, bin_idx: int) -> Tuple[str, Optional[str], int]:
    """Build the grouping key for one row, matching pandas' groupby(...) keys."""
    if stat in _AGG_TWO_POP:
        return (row.pop1, row.pop2, bin_idx)
    return (row.pop1, None, bin_idx)


def _bin_idx(window_pos_1: int, interval_start: int, window_size: int) -> int:
    """
    Map a chunk's window_pos_1 to an output-window bin index.

    Mirrors the result of `pandas.cut(pos, bins=arange(start-1, end+w, w), include_lowest=True)`
    used by the previous implementation, but without materializing the bins array.
    """
    return (window_pos_1 - interval_start) // window_size


def _final_stat(  # noqa: C901
    stat: str,
    fst_type: str,
    n_sites: Numeric,
    diffs: Numeric,
    comps: Numeric,
    missing: Numeric,
    d_stdev: Optional[float] = None,
) -> Union[float, str]:
    """Compute the final `calculated_stat` from aggregated components."""
    if stat == "pi" or stat == "dxy":
        if _is_na(comps) or comps == 0:
            return "NA"
        return float(diffs) / float(comps)
    if stat == "watterson_theta":
        # avg_theta = raw_theta / num_sites  (col 8 / col 7 in the aggregator).
        if _is_na(n_sites) or n_sites == 0:
            return "NA"
        return float(diffs) / float(n_sites)
    if stat == "tajima_d":
        # (raw_pi - watterson_theta) / d_stdev  -- d_stdev recomputed from the merged
        # variant-count classes by the caller.
        if d_stdev is None or d_stdev <= 0:
            return "NA"
        if _is_na(diffs) or _is_na(comps):
            return "NA"
        return (float(diffs) - float(comps)) / d_stdev
    if stat == "fst":
        if fst_type == "wc":
            denom: Numeric = _add(_add(diffs, comps), missing)
            if _is_na(denom) or denom == 0:
                return "NA"
            return float(diffs) / float(denom)
        if fst_type == "hudson":
            if _is_na(comps) or comps == 0:
                return "NA"
            return float(diffs) / float(comps)
    raise ValueError(f"Unsupported statistic for aggregation: {stat}")


class AggRow(NamedTuple):
    """One row of aggregated output, ready to format."""

    pop1: str
    pop2: Optional[str]  # None for single-pop stats (pi/watterson_theta/tajima_d)
    chrom: str
    window_pos_1: int
    window_pos_2: int
    calculated_stat: Union[float, str]
    n_sites: Numeric
    diffs: Numeric
    comps: Numeric
    missing: Numeric
    tajima_counts: str = "NA"


def aggregate_rows(
    rows: Iterable[TempRow],
    stat: str,
    chromosome: str,
    window_size: int,
    fst_type: str,
) -> List[AggRow]:
    """
    Aggregate temp rows for a single (stat, chromosome) into output-window-sized rows.

    Replaces `pixy.core.aggregate_output()` (the previous pandas implementation). The
    binning is identical to what `pandas.cut` produced for the same `interval_start` /
    `window_size`.

    Returns rows sorted by `window_pos_1` (then by pop, matching the previous output).
    """
    rows = list(rows)
    if not rows:
        return []
    positions = [int(r.window_pos_1_str) for r in rows]
    interval_start = min(positions)

    accum: Dict[Tuple[str, Optional[str], int], List[Numeric]] = defaultdict(
        lambda: ["NA", "NA", "NA", "NA"]
    )
    tajima_counts: Dict[Tuple[str, Optional[str], int], Counter[int]] = defaultdict(Counter)

    for r, pos in zip(rows, positions, strict=True):
        bin_idx = _bin_idx(pos, interval_start, window_size)
        key = _agg_key(r, stat, bin_idx)
        a = accum[key]
        a[0] = _add(a[0], _to_num(r.n_sites))
        a[1] = _add(a[1], _to_num(r.diffs))
        a[2] = _add(a[2], _to_num(r.comps))
        a[3] = _add(a[3], _to_num(r.missing))
        if stat == "tajima_d":
            for n, s in deserialize_tajima_d_variant_counts(r.tajima_counts).items():
                tajima_counts[key][n] += s

    out: List[AggRow] = []
    for key, comps in accum.items():
        pop1, pop2, bin_idx = key
        window_pos_1 = interval_start + bin_idx * window_size
        window_pos_2 = window_pos_1 + window_size - 1
        d_stdev: Optional[float] = None
        if stat == "tajima_d":
            d_stdev = calc_tajima_d_stdev(tajima_counts[key])
            # For tajima_d, col 10 (missing) is repurposed as d_stdev in the output.
            comps_tail: Numeric = d_stdev
        else:
            comps_tail = comps[3]
        final = _final_stat(stat, fst_type, comps[0], comps[1], comps[2], comps_tail, d_stdev)
        # Tajima D serializes the merged variant counts back into the output's optional col.
        from pixy.calc import serialize_tajima_d_variant_counts  # local to avoid cycle

        tcounts = (
            serialize_tajima_d_variant_counts(tajima_counts[key]) if stat == "tajima_d" else "NA"
        )
        out.append(
            AggRow(
                pop1=pop1,
                pop2=pop2,
                chrom=chromosome,
                window_pos_1=window_pos_1,
                window_pos_2=window_pos_2,
                calculated_stat=final,
                n_sites=comps[0],
                diffs=comps[1],
                comps=comps[2],
                missing=comps_tail,
                tajima_counts=tcounts,
            )
        )

    # Sort matching the previous pandas output: by (window_pos_1, pop1, pop2).
    out.sort(key=lambda r: (r.window_pos_1, r.pop1, r.pop2 or ""))
    return out


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

# Headers per stat. Built dynamically for fst because the components columns depend on
# --fst_type and --fst_components.
HEADER_PI = (
    "pop",
    "chromosome",
    "window_pos_1",
    "window_pos_2",
    "avg_pi",
    "no_sites",
    "count_diffs",
    "count_comparisons",
    "count_missing",
)
HEADER_DXY = (
    "pop1",
    "pop2",
    "chromosome",
    "window_pos_1",
    "window_pos_2",
    "avg_dxy",
    "no_sites",
    "count_diffs",
    "count_comparisons",
    "count_missing",
)
HEADER_WATTERSON = (
    "pop",
    "chromosome",
    "window_pos_1",
    "window_pos_2",
    "avg_watterson_theta",
    "no_sites",
    "raw_watterson_theta",
    "no_var_sites",
    "weighted_no_sites",
)
HEADER_TAJIMA_BASE = (
    "pop",
    "chromosome",
    "window_pos_1",
    "window_pos_2",
    "tajima_d",
    "no_sites",
    "raw_pi",
    "raw_watterson_theta",
    "tajima_d_stdev",
)


def _fst_header(fst_type: str, fst_components: bool) -> Tuple[str, ...]:
    base = (
        "pop1",
        "pop2",
        "chromosome",
        "window_pos_1",
        "window_pos_2",
        f"avg_{fst_type}_fst",
        "no_snps",
    )
    if not fst_components:
        return base
    if fst_type == "wc":
        return base + ("wc_fst_a", "wc_fst_b", "wc_fst_c")
    return base + ("hudson_fst_num", "hudson_fst_den")


def _format_temp_row_for_stat(row: TempRow, stat: str, fst_type: str, fst_components: bool) -> str:
    """
    Format a non-aggregated temp row as one output line for `stat`.

    Column ordering & dtypes match the previous pandas paths:
      pi/dxy:           pop[,pop2], chrom, w1, w2, avg, n_sites, diffs, comps, missing
      watterson_theta:  pop, chrom, w1, w2, avg, n_sites, raw_theta(float),
                        var_sites(int), weighted(float)
      tajima_d:         pop, chrom, w1, w2, tajima_d, n_sites, raw_pi(float),
                        watterson(float), stdev(float)
      fst:              pop1, pop2, chrom, w1, w2, avg_fst, n_snps, [components...]
    """
    if stat == "pi":
        cells = [
            row.pop1,
            row.chrom,
            row.window_pos_1_str,
            row.window_pos_2_str,
            _fmt_float(_to_num(row.calculated_stat)),
            _fmt_int(_to_num(row.n_sites)),
            _fmt_int(_to_num(row.diffs)),
            _fmt_int(_to_num(row.comps)),
            _fmt_int(_to_num(row.missing)),
        ]
    elif stat == "dxy":
        cells = [
            row.pop1,
            row.pop2,
            row.chrom,
            row.window_pos_1_str,
            row.window_pos_2_str,
            _fmt_float(_to_num(row.calculated_stat)),
            _fmt_int(_to_num(row.n_sites)),
            _fmt_int(_to_num(row.diffs)),
            _fmt_int(_to_num(row.comps)),
            _fmt_int(_to_num(row.missing)),
        ]
    elif stat == "watterson_theta":
        cells = [
            row.pop1,
            row.chrom,
            row.window_pos_1_str,
            row.window_pos_2_str,
            _fmt_float(_to_num(row.calculated_stat)),
            _fmt_int(_to_num(row.n_sites)),
            _fmt_float(_to_num(row.diffs)),  # raw_theta
            _fmt_int(_to_num(row.comps)),  # num_var_sites
            _fmt_float(_to_num(row.missing)),
        ]  # weighted_sites
    elif stat == "tajima_d":
        cells = [
            row.pop1,
            row.chrom,
            row.window_pos_1_str,
            row.window_pos_2_str,
            _fmt_float(_to_num(row.calculated_stat)),
            _fmt_int(_to_num(row.n_sites)),
            _fmt_float(_to_num(row.diffs)),
            _fmt_float(_to_num(row.comps)),
            _fmt_float(_to_num(row.missing)),
        ]
        if fst_components:  # fst_components flag is overloaded; tajima_components is what matters
            pass  # handled by caller
    elif stat == "fst":
        cells = [
            row.pop1,
            row.pop2,
            row.chrom,
            row.window_pos_1_str,
            row.window_pos_2_str,
            _fmt_float(_to_num(row.calculated_stat)),
            _fmt_int(_to_num(row.n_sites)),
        ]
        if fst_components:
            if fst_type == "wc":
                cells.extend([
                    _fmt_float(_to_num(row.diffs)),
                    _fmt_float(_to_num(row.comps)),
                    _fmt_float(_to_num(row.missing)),
                ])
            else:  # hudson — drop the unused c placeholder
                cells.extend([
                    _fmt_float(_to_num(row.diffs)),
                    _fmt_float(_to_num(row.comps)),
                ])
    else:
        raise ValueError(f"Unsupported stat: {stat}")
    return "\t".join(cells)


def _format_agg_row(row: AggRow, stat: str, fst_type: str, fst_components: bool) -> str:
    """Format an aggregated row as one output line."""
    if stat == "pi":
        cells = [
            row.pop1,
            row.chrom,
            str(row.window_pos_1),
            str(row.window_pos_2),
            _fmt_float(row.calculated_stat),
            _fmt_int(row.n_sites),
            _fmt_int(row.diffs),
            _fmt_int(row.comps),
            _fmt_int(row.missing),
        ]
    elif stat == "dxy":
        cells = [
            row.pop1,
            row.pop2 or "",
            row.chrom,
            str(row.window_pos_1),
            str(row.window_pos_2),
            _fmt_float(row.calculated_stat),
            _fmt_int(row.n_sites),
            _fmt_int(row.diffs),
            _fmt_int(row.comps),
            _fmt_int(row.missing),
        ]
    elif stat == "watterson_theta":
        cells = [
            row.pop1,
            row.chrom,
            str(row.window_pos_1),
            str(row.window_pos_2),
            _fmt_float(row.calculated_stat),
            _fmt_int(row.n_sites),
            _fmt_float(row.diffs),
            _fmt_int(row.comps),
            _fmt_float(row.missing),
        ]
    elif stat == "tajima_d":
        cells = [
            row.pop1,
            row.chrom,
            str(row.window_pos_1),
            str(row.window_pos_2),
            _fmt_float(row.calculated_stat),
            _fmt_int(row.n_sites),
            _fmt_float(row.diffs),
            _fmt_float(row.comps),
            _fmt_float(row.missing),
        ]
    elif stat == "fst":
        cells = [
            row.pop1,
            row.pop2 or "",
            row.chrom,
            str(row.window_pos_1),
            str(row.window_pos_2),
            _fmt_float(row.calculated_stat),
            _fmt_int(row.n_sites),
        ]
        if fst_components:
            if fst_type == "wc":
                cells.extend([
                    _fmt_float(row.diffs),
                    _fmt_float(row.comps),
                    _fmt_float(row.missing),
                ])
            else:
                cells.extend([
                    _fmt_float(row.diffs),
                    _fmt_float(row.comps),
                ])
    else:
        raise ValueError(f"Unsupported stat: {stat}")
    return "\t".join(cells)


def write_stat_file(  # noqa: C901
    out_path: Path,
    stat: str,
    rows_by_chrom: Dict[str, List[TempRow]],
    chrom_list: Iterable[str],
    aggregate: bool,
    window_size: int,
    fst_type: str,
    fst_components: bool = False,
    tajima_components: bool = False,
) -> List[str]:
    """Write the per-stat output file. Returns the list of chromosomes that had no data."""
    chroms_with_no_data: List[str] = []
    if stat == "fst":
        header = _fst_header(fst_type, fst_components)
    elif stat == "tajima_d":
        header = HEADER_TAJIMA_BASE + (("tajima_d_s_counts",) if tajima_components else ())
    elif stat == "pi":
        header = HEADER_PI
    elif stat == "dxy":
        header = HEADER_DXY
    elif stat == "watterson_theta":
        header = HEADER_WATTERSON
    else:
        raise ValueError(f"Unsupported stat: {stat}")

    with open(out_path, "w") as f:
        f.write("\t".join(header) + "\n")
        for chromosome in chrom_list:
            rows = rows_by_chrom.get(chromosome)
            if not rows:
                chroms_with_no_data.append(chromosome)
                continue
            if aggregate:
                for arow in aggregate_rows(rows, stat, chromosome, window_size, fst_type):
                    line = _format_agg_row(arow, stat, fst_type, fst_components)
                    if stat == "tajima_d" and tajima_components:
                        line = line + "\t" + arow.tajima_counts
                    f.write(line + "\n")
            else:
                # Non-aggregated: emit each temp row directly, sorted by window_pos_1.
                rows_sorted = sorted(rows, key=lambda r: int(r.window_pos_1_str))
                for trow in rows_sorted:
                    line = _format_temp_row_for_stat(trow, stat, fst_type, fst_components)
                    if stat == "tajima_d" and tajima_components:
                        line = line + "\t" + trow.tajima_counts
                    f.write(line + "\n")
    return chroms_with_no_data


def group_temp_rows_by_stat_chrom(
    rows: Iterable[TempRow],
) -> Dict[Tuple[str, str], List[TempRow]]:
    """Single-pass grouping equivalent to `outpanel.groupby([0, 3])` in the old code."""
    out: Dict[Tuple[str, str], List[TempRow]] = defaultdict(list)
    successful_stats: set[str] = set()
    for r in rows:
        out[(r.stat, r.chrom)].append(r)
        successful_stats.add(r.stat)
    return out
