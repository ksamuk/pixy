"""
Tests for the wisp-mask code path.

Covers:
  * `WispMetadata.from_header_line` parsing — happy path + malformed inputs.
  * `validate_wisp_data_rows` — well-formed BEDs pass; truncated, malformed,
    or out-of-range BEDs return a descriptive error.
  * End-to-end pixy run with `--wisp_bed` on a tiny synthetic dataset that
    reproduces the crash scenario from `pixy/core.py:531` (a chromosome where
    the variants-only VCF has zero records but the wisp mask still supplies
    invariant-site coverage).
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from pixy.wisp import WispMask
from pixy.wisp import WispMetadata
from pixy.wisp import validate_wisp_against_populations
from pixy.wisp import validate_wisp_data_rows
from tests.conftest import run_pixy_helper

HEADER_OK = (
    "#wisp_mask_metadata\t"
    '{"populations": ["A", "B"], '
    '"population_sample_counts": {"A": 2, "B": 2}, '
    '"population_columns": [{"column_number": 4, "name": "A"}, '
    '{"column_number": 5, "name": "B"}], '
    '"threshold": 10}'
)


def _write_wisp_bed(tmp_path: Path, body: str, name: str = "test.bed") -> Path:
    """Write a wisp BED (header + body lines), gzip + tabix-index it, return the .gz path."""
    raw = tmp_path / name
    raw.write_text(HEADER_OK + "\n" + body)
    subprocess.run(["bgzip", "-f", str(raw)], check=True)
    bgz = raw.with_suffix(raw.suffix + ".gz")
    subprocess.run(["tabix", "-p", "bed", str(bgz)], check=True)
    return bgz


def _write_plain_wisp(tmp_path: Path, body: str, name: str = "test.bed") -> Path:
    """Write a plain (non-bgzipped) wisp BED. Skips bgzip/tabix; parser-only fixtures."""
    raw = tmp_path / name
    raw.write_text(HEADER_OK + "\n" + body)
    return raw


def _load_wisp_mask(path: Path) -> WispMask:
    """Load a WispMask from `path` (transparent .gz vs plain)."""
    return WispMask.from_path(path)


# ---------------------------------------------------------------------------
# WispMetadata header parsing
# ---------------------------------------------------------------------------


def test_metadata_from_header_line_basic() -> None:
    """A well-formed header line parses into the expected populations / counts / columns."""
    meta = WispMetadata.from_header_line(HEADER_OK + "\n")
    assert meta.populations == ("A", "B")
    assert meta.population_sample_counts == {"A": 2, "B": 2}
    assert meta.threshold == 10
    # `column_number` is 1-based in the JSON; we expose 0-based field indices.
    assert meta.column_index_for == {"A": 3, "B": 4}


def test_metadata_from_header_line_missing_prefix() -> None:
    """A header line lacking the `#wisp_mask_metadata` prefix is rejected."""
    with pytest.raises(ValueError, match="missing the required"):
        WispMetadata.from_header_line('{"populations": ["A"]}\n')


def test_metadata_from_header_line_bad_json() -> None:
    """A malformed JSON payload after the prefix is rejected."""
    with pytest.raises(ValueError, match="not valid JSON"):
        WispMetadata.from_header_line("#wisp_mask_metadata\t{not-json}\n")


def test_metadata_from_header_line_missing_required_key() -> None:
    """A header missing one of the required JSON keys is rejected."""
    line = '#wisp_mask_metadata\t{"populations": ["A"]}\n'
    with pytest.raises(ValueError, match="missing required key"):
        WispMetadata.from_header_line(line)


def test_metadata_from_header_line_column_out_of_range() -> None:
    """`column_number` <= 3 collides with chrom/start/end and is rejected."""
    # column_number=2 → 0-based index 1, before chrom/start/end (must be >= 3).
    bad = (
        "#wisp_mask_metadata\t"
        '{"populations": ["A"], '
        '"population_sample_counts": {"A": 2}, '
        '"population_columns": [{"column_number": 2, "name": "A"}]}\n'
    )
    with pytest.raises(ValueError, match="must come after chrom/start/end"):
        WispMetadata.from_header_line(bad)


# ---------------------------------------------------------------------------
# validate_wisp_data_rows
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    shutil.which("bgzip") is None or shutil.which("tabix") is None,
    reason="bgzip and tabix are required to build wisp BED fixtures",
)
def test_validate_data_rows_ok(tmp_path: Path) -> None:
    """A well-formed (bgzipped, tabix-indexed) wisp BED passes validation."""
    bgz = _write_wisp_bed(
        tmp_path,
        body="chr1\t0\t100\t2\t2\nchr1\t100\t200\t2\t1\nchr2\t0\t500\t1\t2\n",
    )
    assert validate_wisp_data_rows(_load_wisp_mask(bgz)) is None


def test_validate_data_rows_header_only(tmp_path: Path) -> None:
    """A header-only wisp BED (no data rows) is rejected with a clear message."""
    raw = _write_plain_wisp(tmp_path, body="")
    err = validate_wisp_data_rows(_load_wisp_mask(raw))
    assert err is not None and "no data rows found" in err


def test_validate_data_rows_truncated_fields(tmp_path: Path) -> None:
    """A row with fewer columns than the metadata declares is rejected."""
    # Header declares pop column index 4 (1-based 5) but row only has 4 fields.
    raw = _write_plain_wisp(tmp_path, body="chr1\t0\t100\t2\n")
    err = validate_wisp_data_rows(_load_wisp_mask(raw))
    assert err is not None
    assert "tab-delimited fields" in err
    assert "index 4" in err


def test_validate_data_rows_non_integer_coords(tmp_path: Path) -> None:
    """A row whose start/end is not an integer is rejected."""
    raw = _write_plain_wisp(tmp_path, body="chr1\tnope\t100\t2\t2\n")
    err = validate_wisp_data_rows(_load_wisp_mask(raw))
    assert err is not None and "non-integer coordinates" in err


def test_validate_data_rows_inverted_coords(tmp_path: Path) -> None:
    """A row with end <= start is rejected."""
    raw = _write_plain_wisp(tmp_path, body="chr1\t100\t100\t2\t2\n")
    err = validate_wisp_data_rows(_load_wisp_mask(raw))
    assert err is not None and "end <= start" in err


def test_validate_data_rows_count_exceeds_n(tmp_path: Path) -> None:
    """A per-population count larger than the declared sample count is rejected."""
    # Header says N_A=2, but row reports 3 callable samples for A.
    raw = _write_plain_wisp(tmp_path, body="chr1\t0\t100\t3\t2\n")
    err = validate_wisp_data_rows(_load_wisp_mask(raw))
    assert err is not None and "expected 0..2" in err


def test_validate_data_rows_negative_count(tmp_path: Path) -> None:
    """A negative per-population count is rejected."""
    raw = _write_plain_wisp(tmp_path, body="chr1\t0\t100\t-1\t2\n")
    err = validate_wisp_data_rows(_load_wisp_mask(raw))
    assert err is not None and "expected 0..2" in err


def test_validate_data_rows_non_integer_count(tmp_path: Path) -> None:
    """A non-integer per-population count is rejected."""
    raw = _write_plain_wisp(tmp_path, body="chr1\t0\t100\tx\t2\n")
    err = validate_wisp_data_rows(_load_wisp_mask(raw))
    assert err is not None and "non-integer count" in err


# ---------------------------------------------------------------------------
# validate_wisp_against_populations
# ---------------------------------------------------------------------------


def test_validate_against_populations_ok(tmp_path: Path) -> None:
    """Matching populations + sample counts pass the cross-check."""
    raw = _write_plain_wisp(tmp_path, body="chr1\t0\t100\t2\t2\n")
    mask = _load_wisp_mask(raw)
    assert validate_wisp_against_populations(mask, {"A": 2, "B": 2}) is None


def test_validate_against_populations_mismatched_names(tmp_path: Path) -> None:
    """Population names that disagree between the popfile and the mask are surfaced."""
    raw = _write_plain_wisp(tmp_path, body="chr1\t0\t100\t2\t2\n")
    mask = _load_wisp_mask(raw)
    err = validate_wisp_against_populations(mask, {"A": 2, "C": 2})
    assert err is not None
    assert "only in populations file" in err
    assert "only in wisp mask" in err


def test_validate_against_populations_mismatched_counts(tmp_path: Path) -> None:
    """Population sample counts that disagree between the popfile and the mask are surfaced."""
    raw = _write_plain_wisp(tmp_path, body="chr1\t0\t100\t2\t2\n")
    mask = _load_wisp_mask(raw)
    err = validate_wisp_against_populations(mask, {"A": 2, "B": 3})
    assert err is not None and "B: wisp=2, populations_file=3" in err


# ---------------------------------------------------------------------------
# End-to-end pixy regression test
# ---------------------------------------------------------------------------


def _build_wisp_dataset(tmp_path: Path) -> dict[str, Path]:
    """
    Build a 4-sample 2-population synthetic dataset that exercises the wisp path.

    chr1 has variants but chr2 has zero usable VCF records — exactly the pattern that
    crashed pixy's wisp path before the fix in ``pixy/core.py``. The wisp mask
    covers both chroms.
    """
    # chr1 carries usable variants. chr2 carries only an indel — ploidy inference
    # sees a record (so the per-contig ploidy check passes), but pixy's biallelic
    # SNP mask drops it, leaving the chr2 chunk with `gt_array=None`. That's the
    # state the `pixy/core.py` fix has to handle.
    vcf_text = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        "##contig=<ID=chr2,length=1000>\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
        "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/0\n"
        "chr1\t250\t.\tG\tC\t.\tPASS\t.\tGT\t0/1\t0/0\t0/1\t1/1\n"
        "chr2\t500\t.\tAC\tA\t.\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/1\n"
    )
    vcf = tmp_path / "data.vcf"
    vcf.write_text(vcf_text)
    subprocess.run(["bgzip", "-f", str(vcf)], check=True)
    vcf_gz = vcf.with_suffix(".vcf.gz")
    subprocess.run(["tabix", "-p", "vcf", str(vcf_gz)], check=True)

    pop_file = tmp_path / "pops.txt"
    pop_file.write_text("S1\tA\nS2\tA\nS3\tB\nS4\tB\n")

    # Wisp mask: full callability on both chroms (both pops at N).
    bed_raw = tmp_path / "wisp.bed"
    bed_raw.write_text(HEADER_OK + "\n" + "chr1\t0\t1000\t2\t2\n" + "chr2\t0\t1000\t2\t2\n")
    subprocess.run(["bgzip", "-f", str(bed_raw)], check=True)
    bed_gz = bed_raw.with_suffix(".bed.gz")
    subprocess.run(["tabix", "-p", "bed", str(bed_gz)], check=True)

    return {"vcf": vcf_gz, "pops": pop_file, "wisp": bed_gz}


@pytest.mark.skipif(
    shutil.which("bgzip") is None or shutil.which("tabix") is None,
    reason="bgzip and tabix are required to build the synthetic wisp dataset",
)
def test_pixy_with_wisp_mask_runs(tmp_path: Path) -> None:
    """
    Regression: wisp-only chromosomes no longer crash with ``gt_region is None``.

    pixy used to fail when a chromosome had zero usable VCF records but the wisp
    mask still supplied invariant coverage. With the ``pixy/core.py`` fix, a
    synthesized zero-row GenotypeArray keeps the per-stat code paths alive and chr2
    ends up reported with wisp-derived denominators.
    """
    inputs = _build_wisp_dataset(tmp_path)
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    run_pixy_helper(
        pixy_out_dir=out_dir,
        stats=["pi", "dxy", "watterson_theta", "tajima_d"],
        vcf_path=inputs["vcf"],
        populations_path=inputs["pops"],
        window_size=500,
        wisp_bed_path=inputs["wisp"],
        bypass_invariant_check=False,  # wisp mask implies this on
    )

    pi_path = out_dir / "pixy_pi.txt"
    dxy_path = out_dir / "pixy_dxy.txt"
    assert pi_path.exists(), "pixy did not emit pi output"
    assert dxy_path.exists(), "pixy did not emit dxy output"

    pi_lines = pi_path.read_text().splitlines()
    # Header + chr2 rows for both pops (would have been absent / crashed before fix).
    chr2_rows = [ln for ln in pi_lines[1:] if "\tchr2\t" in ln]
    assert len(chr2_rows) >= 2, (
        f"Expected wisp-derived pi rows for chr2 across both populations, got:\n{chr2_rows!r}"
    )
    # And every chr2 row should report a real (non-NA) site count, since the
    # wisp mask supplies invariant denominators even when the VCF is empty.
    for row in chr2_rows:
        cols = row.split("\t")
        no_sites = cols[5]
        assert no_sites != "NA" and int(no_sites) > 0, (
            f"chr2 pi row has NA / zero no_sites with wisp mask active: {row!r}"
        )


@pytest.mark.skipif(
    shutil.which("bgzip") is None or shutil.which("tabix") is None,
    reason="bgzip and tabix are required to build the synthetic wisp dataset",
)
def test_pixy_with_wisp_mask_rejects_header_only_bed(tmp_path: Path) -> None:
    """A wisp BED with no data rows is rejected at args-validation time."""
    inputs = _build_wisp_dataset(tmp_path)
    # Overwrite the wisp BED with a header-only file.
    raw = tmp_path / "header_only.bed"
    raw.write_text(HEADER_OK + "\n")
    subprocess.run(["bgzip", "-f", str(raw)], check=True)
    bed_gz = raw.with_suffix(".bed.gz")
    subprocess.run(["tabix", "-p", "bed", str(bed_gz)], check=True)

    out_dir = tmp_path / "out2"
    out_dir.mkdir()

    with pytest.raises(ValueError, match="no data rows found"):
        run_pixy_helper(
            pixy_out_dir=out_dir,
            stats=["pi"],
            vcf_path=inputs["vcf"],
            populations_path=inputs["pops"],
            window_size=500,
            wisp_bed_path=bed_gz,
        )


@pytest.mark.skipif(
    shutil.which("bgzip") is None or shutil.which("tabix") is None,
    reason="bgzip and tabix are required to build the synthetic wisp dataset",
)
def test_pixy_with_wisp_mask_rejects_population_mismatch(tmp_path: Path) -> None:
    """Wisp metadata that doesn't match --populations is rejected up front."""
    inputs = _build_wisp_dataset(tmp_path)
    # Populations file that names different pops than the wisp mask carries.
    bad_pops = tmp_path / "bad_pops.txt"
    bad_pops.write_text("S1\tX\nS2\tX\nS3\tY\nS4\tY\n")

    out_dir = tmp_path / "out3"
    out_dir.mkdir()

    with pytest.raises(ValueError, match="do not match --populations file"):
        run_pixy_helper(
            pixy_out_dir=out_dir,
            stats=["pi"],
            vcf_path=inputs["vcf"],
            populations_path=bad_pops,
            window_size=500,
            wisp_bed_path=inputs["wisp"],
        )


@pytest.mark.skipif(
    shutil.which("bgzip") is None or shutil.which("tabix") is None,
    reason="bgzip and tabix are required to build the synthetic wisp dataset",
)
def test_pixy_with_wisp_mask_missing_index(tmp_path: Path) -> None:
    """A wisp BED without a .tbi/.csi index is rejected with a clear message."""
    inputs = _build_wisp_dataset(tmp_path)
    # Strip the tabix index for the existing wisp BED.
    tbi = inputs["wisp"].with_suffix(inputs["wisp"].suffix + ".tbi")
    if tbi.exists():
        tbi.unlink()

    out_dir = tmp_path / "out4"
    out_dir.mkdir()

    with pytest.raises(ValueError, match="not indexed"):
        run_pixy_helper(
            pixy_out_dir=out_dir,
            stats=["pi"],
            vcf_path=inputs["vcf"],
            populations_path=inputs["pops"],
            window_size=500,
            wisp_bed_path=inputs["wisp"],
        )
