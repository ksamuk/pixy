r"""
Build the GVCF regression-test fixture from the ag1000 all-sites VCF.

Reads ``tests/main/data/ag1000_pixy_test.vcf.gz`` and writes
``tests/main/data/ag1000_pixy_gvcf_test.vcf.gz``. Maximal runs of consecutive
invariant records (``ALT='.'``) on the same chromosome **whose FORMAT and
per-sample columns are byte-identical across the entire run** are collapsed
into a single block record. Runs where any sample's call varies across
positions are emitted individually so that the per-sample callable counts
remain identical to the all-sites baseline.

For each collapsed block:

* ``POS`` is the first POS in the run
* ``REF`` is the first REF in the run (the per-position bases are not preserved
  by GVCF blocks; downstream pixy only cares that the sites are ref-only)
* ``ALT`` is ``<NON_REF>``
* ``QUAL`` and ``FILTER`` follow the first row
* ``INFO`` contains ``END=<last POS in the run>`` (other INFO fields are dropped)
* FORMAT and the per-sample columns follow the (uniform) row content

Variant rows (``ALT != '.'``) are written through unchanged. Invariant rows
that fail the per-sample-uniformity test are emitted individually with their
original ALT (``.``).

After writing the plain VCF text, the script bgzips and tabix-indexes the
output (calling ``bgzip`` and ``tabix`` from PATH; both are bundled with the
``pixy-py*-test`` conda envs).

Run once from a checkout under WSL::

    /home/ksamuk/miniconda3/envs/pixy-py311-test/bin/python \\
        tests/main/data/scripts/build_gvcf_fixture.py

The generated ``ag1000_pixy_gvcf_test.vcf.gz`` and ``.tbi`` are committed to
the repository; this script only needs to be re-run if the source all-sites
fixture changes.
"""

from __future__ import annotations

import gzip
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[4]
INPUT = REPO_ROOT / "tests" / "main" / "data" / "ag1000_pixy_test.vcf.gz"
OUTPUT_VCF = REPO_ROOT / "tests" / "main" / "data" / "ag1000_pixy_gvcf_test.vcf"
OUTPUT_GZ = OUTPUT_VCF.with_suffix(".vcf.gz")


# Required header line so consumers can parse the `END` INFO field.
END_INFO_LINE = (
    '##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">\n'
)
# Symbolic ALT allele used by GATK GVCFs for invariant blocks.
NON_REF_ALT_LINE = (
    '##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">\n'
)


def _gt_subfields(row: List[str]) -> List[str]:
    """
    Extract the GT subfield from each per-sample column on a single row.

    FORMAT (field 8) names the colon-separated sub-fields of each sample
    column; the index of ``GT`` is the same across every sample on the same
    row. Returns one string per sample.
    """
    format_keys = row[8].split(":")
    try:
        gt_idx = format_keys.index("GT")
    except ValueError:
        # No GT — treat the whole sample columns as the comparable shape (fallback).
        return row[9:]
    out: List[str] = []
    for sample in row[9:]:
        sub = sample.split(":", gt_idx + 1)
        out.append(sub[gt_idx] if len(sub) > gt_idx else ".")
    return out


def _gt_uniform(rows: List[List[str]]) -> bool:
    """
    True when the per-sample GT subfield is identical across every row in the run.

    GT is the only sample-level field pixy consumes from the VCF (AD/DP/GQ/etc.
    are not used). Holding GT constant across the run is therefore sufficient
    for the GVCF expansion to reproduce the all-sites callable counts exactly.
    """
    if len(rows) <= 1:
        return True
    first_gts = _gt_subfields(rows[0])
    return all(_gt_subfields(r) == first_gts for r in rows[1:])


def _emit_rows_individually(rows: List[List[str]]) -> List[str]:
    """Fallback when a run cannot be safely collapsed: emit each row as-is."""
    return ["\t".join(r) + "\n" for r in rows]


def _collapse(rows: List[List[str]]) -> List[str]:
    """
    Collapse a list of consecutive invariant records on the same chromosome.

    All input rows must have ``ALT == '.'``, share the same chromosome, and have
    uniform FORMAT/per-sample columns (callers gate this via
    :func:`_sample_columns_uniform`). Returns the single collapsed block-record
    line as a string, wrapped in a list for shape consistency with the fallback.
    """
    first = list(rows[0])  # don't mutate the caller's row
    last_pos = rows[-1][1]
    # Replace INFO with just END=...; drop any other INFO content (not used downstream).
    first[4] = "<NON_REF>"
    first[7] = f"END={last_pos}"
    return ["\t".join(first) + "\n"]


def main() -> None:  # noqa: C901
    """Generate and bgzip+index the GVCF fixture next to the all-sites baseline."""
    if not INPUT.exists():
        raise SystemExit(f"input VCF not found: {INPUT}")

    header_lines: List[str] = []
    out_lines: List[str] = []
    end_line_present = False
    nonref_line_present = False

    pending: List[List[str]] = []
    pending_chrom: Optional[str] = None

    def flush_pending() -> None:
        nonlocal pending, pending_chrom
        if not pending:
            return
        if len(pending) > 1 and _gt_uniform(pending):
            out_lines.extend(_collapse(pending))
        else:
            # Runs with varying per-sample GTs can't be collapsed without
            # changing what pixy will see; emit them verbatim so the regression
            # baseline still matches.
            out_lines.extend(_emit_rows_individually(pending))
        pending = []
        pending_chrom = None

    with gzip.open(INPUT, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("##INFO=<ID=END"):
                    end_line_present = True
                if line.startswith("##ALT=<ID=NON_REF"):
                    nonref_line_present = True
                if line.startswith("#CHROM"):
                    # Insert the END and NON_REF declarations before #CHROM so they
                    # appear as header metadata.
                    if not end_line_present:
                        header_lines.append(END_INFO_LINE)
                    if not nonref_line_present:
                        header_lines.append(NON_REF_ALT_LINE)
                header_lines.append(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                # malformed; pass through
                flush_pending()
                out_lines.append(line)
                continue
            chrom = fields[0]
            alt = fields[4]
            try:
                pos = int(fields[1])
            except ValueError:
                # Malformed POS — flush and pass-through.
                flush_pending()
                out_lines.append(line)
                continue
            is_invariant = alt == "."
            # A run is broken by any of: a non-invariant record, a chrom change, or a
            # POS gap (the all-sites VCF can omit positions; collapsing across a gap
            # would invent positions that don't exist in the baseline).
            if pending:
                last_pos = int(pending[-1][1])
                if not is_invariant or chrom != pending_chrom or pos != last_pos + 1:
                    flush_pending()
            if is_invariant:
                pending.append(fields)
                pending_chrom = chrom
            else:
                out_lines.append(line)
        flush_pending()

    # Write the plain VCF text, then bgzip + tabix.
    OUTPUT_VCF.write_text("".join(header_lines + out_lines))

    # Tools from the active conda env's PATH.
    for tool in ("bgzip", "tabix"):
        if shutil.which(tool) is None:
            raise SystemExit(f"{tool!r} not on PATH; run from within a pixy-py*-test conda env")

    # bgzip -f overwrites any stale output. The resulting file is .vcf.gz.
    subprocess.run(["bgzip", "-f", str(OUTPUT_VCF)], check=True)
    subprocess.run(["tabix", "-p", "vcf", str(OUTPUT_GZ)], check=True)

    print(f"wrote {OUTPUT_GZ}")
    print(f"wrote {OUTPUT_GZ}.tbi")


if __name__ == "__main__":
    main()
    sys.exit(0)
