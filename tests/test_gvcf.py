"""
Unit tests for the GVCF expansion module and its validation hooks.

The expand-block tests build synthetic callset dicts in the shape returned by
:func:`allel.read_vcf` and exercise :func:`pixy.gvcf.expand_blocks` directly,
keeping the suite fast and independent of any VCF fixture. The validation
tests drive :func:`pixy.__main__.main` through the standard
``run_pixy_helper`` harness on real fixtures.
"""

import logging
import shutil
import subprocess
from pathlib import Path
from typing import Any
from typing import Dict

import numpy as np
import pytest
from numpy.typing import NDArray

from pixy.gvcf import expand_blocks
from pixy.gvcf import strip_symbolic_alts
from tests.conftest import assert_files_are_consistent
from tests.conftest import run_pixy_helper


def _make_callset(
    pos: list[int],
    end: list[int],
    numalt: list[int],
    is_snp: list[int],
    gt: NDArray,
    chrom: list[str] | None = None,
) -> Dict[str, NDArray]:
    """Build a callset dict mirroring scikit-allel's return shape."""
    cs: Dict[str, NDArray] = {
        "variants/POS": np.asarray(pos, dtype=np.int32),
        "variants/END": np.asarray(end, dtype=np.int32),
        "variants/numalt": np.asarray(numalt, dtype=np.int32),
        "variants/is_snp": np.asarray(is_snp, dtype=np.int8),
        "calldata/GT": gt,
    }
    if chrom is not None:
        cs["variants/CHROM"] = np.asarray(chrom, dtype=object)
    return cs


def test_no_blocks_passthrough_when_end_missing() -> None:
    """A callset whose END is uniformly -1 (scikit-allel's missing sentinel) is unchanged."""
    gt = np.zeros((3, 2, 2), dtype=np.int8)
    cs = _make_callset(
        pos=[10, 20, 30],
        end=[-1, -1, -1],
        numalt=[1, 0, 1],
        is_snp=[1, 0, 1],
        gt=gt,
        chrom=["chr1"] * 3,
    )
    out = expand_blocks(cs, region_start=1, region_end=100)
    np.testing.assert_array_equal(out["variants/POS"], [10, 20, 30])
    np.testing.assert_array_equal(out["variants/numalt"], [1, 0, 1])
    np.testing.assert_array_equal(out["variants/is_snp"], [1, 0, 1])
    np.testing.assert_array_equal(out["calldata/GT"], gt)


def test_block_expands_to_per_site_rows() -> None:
    """A single block [100..103] yields 4 invariant rows with replicated genotype."""
    block_gt = np.array([[[0, 0], [0, 0]]], dtype=np.int8)  # 1 row, 2 samples, ploidy=2
    cs = _make_callset(
        pos=[100],
        end=[103],
        numalt=[1],  # blocks carry numalt=1 from <NON_REF>
        is_snp=[0],
        gt=block_gt,
        chrom=["chr1"],
    )
    out = expand_blocks(cs, region_start=1, region_end=200)
    np.testing.assert_array_equal(out["variants/POS"], [100, 101, 102, 103])
    # All synthesised rows are now invariant per the post-expansion contract.
    np.testing.assert_array_equal(out["variants/numalt"], [0, 0, 0, 0])
    np.testing.assert_array_equal(out["variants/is_snp"], [0, 0, 0, 0])
    np.testing.assert_array_equal(out["calldata/GT"].shape, (4, 2, 2))
    np.testing.assert_array_equal(out["calldata/GT"][0], block_gt[0])
    np.testing.assert_array_equal(out["calldata/GT"][3], block_gt[0])
    np.testing.assert_array_equal(out["variants/CHROM"], ["chr1"] * 4)


def test_mixed_blocks_and_variants_sorted_by_pos() -> None:
    """Mixed records — variants and blocks — interleave correctly after expansion."""
    # Layout: variant at 95, block 100..102, variant at 105
    gt = np.array(
        [
            [[0, 1], [0, 0]],  # variant @95
            [[0, 0], [0, 0]],  # block @100..102
            [[1, 1], [0, 1]],  # variant @105
        ],
        dtype=np.int8,
    )
    cs = _make_callset(
        pos=[95, 100, 105],
        end=[-1, 102, -1],
        numalt=[1, 1, 1],
        is_snp=[1, 0, 1],
        gt=gt,
    )
    out = expand_blocks(cs, region_start=1, region_end=200)
    np.testing.assert_array_equal(out["variants/POS"], [95, 100, 101, 102, 105])
    np.testing.assert_array_equal(out["variants/is_snp"], [1, 0, 0, 0, 1])
    np.testing.assert_array_equal(out["variants/numalt"], [1, 0, 0, 0, 1])
    # Variant rows preserved verbatim
    np.testing.assert_array_equal(out["calldata/GT"][0], gt[0])
    np.testing.assert_array_equal(out["calldata/GT"][4], gt[2])
    # Block rows carry the replicated GT
    for i in (1, 2, 3):
        np.testing.assert_array_equal(out["calldata/GT"][i], gt[1])


def test_left_overlap_block_clipped_to_region() -> None:
    """
    Blocks starting before ``region_start`` are clipped on the left.

    The caller-side tabix widening is what makes such a block visible in the
    first place; expand_blocks must then drop the overshoot.
    """
    block_gt = np.array([[[0, 0]]], dtype=np.int8)  # 1 row, 1 sample, ploidy=2
    cs = _make_callset(
        pos=[100],
        end=[200],
        numalt=[1],
        is_snp=[0],
        gt=block_gt,
    )
    out = expand_blocks(cs, region_start=150, region_end=175)
    expected = np.arange(150, 176, dtype=np.int32)
    np.testing.assert_array_equal(out["variants/POS"], expected)
    assert out["calldata/GT"].shape == (26, 1, 2)


def test_right_overshoot_block_clipped_to_region() -> None:
    """A block extending past `region_end` is clipped at the right edge."""
    block_gt = np.array([[[0, 0]]], dtype=np.int8)
    cs = _make_callset(
        pos=[50],
        end=[60],
        numalt=[1],
        is_snp=[0],
        gt=block_gt,
    )
    out = expand_blocks(cs, region_start=1, region_end=55)
    np.testing.assert_array_equal(out["variants/POS"], [50, 51, 52, 53, 54, 55])


def test_haploid_block_expansion() -> None:
    """Haploid GTs (no ploidy axis on the inner dim) expand cleanly."""
    block_gt = np.array([[0, 0]], dtype=np.int8)  # 1 row, 2 samples, no ploidy dim
    cs = _make_callset(
        pos=[10],
        end=[12],
        numalt=[1],
        is_snp=[0],
        gt=block_gt,
    )
    out = expand_blocks(cs, region_start=1, region_end=100)
    np.testing.assert_array_equal(out["variants/POS"], [10, 11, 12])
    assert out["calldata/GT"].shape == (3, 2)


def test_empty_after_clip_returns_empty_arrays() -> None:
    """When clipping removes every row, the output has length-0 arrays (not None)."""
    block_gt = np.array([[[0, 0]]], dtype=np.int8)
    cs = _make_callset(
        pos=[100],
        end=[110],
        numalt=[1],
        is_snp=[0],
        gt=block_gt,
    )
    out = expand_blocks(cs, region_start=500, region_end=600)
    assert out["variants/POS"].shape == (0,)
    assert out["calldata/GT"].shape[0] == 0


def test_no_end_field_at_all_is_passthrough() -> None:
    """Missing 'variants/END' key (older callers) leaves the callset alone except for clipping."""
    gt = np.zeros((2, 1, 2), dtype=np.int8)
    cs: Dict[str, NDArray] = {
        "variants/POS": np.array([10, 20], dtype=np.int32),
        "variants/numalt": np.array([1, 0], dtype=np.int32),
        "variants/is_snp": np.array([1, 0], dtype=np.int8),
        "calldata/GT": gt,
    }
    out = expand_blocks(cs, region_start=1, region_end=100)
    np.testing.assert_array_equal(out["variants/POS"], [10, 20])


@pytest.mark.parametrize(
    "region_start,region_end,expected_positions",
    [
        (1, 100, list(range(10, 21))),
        (12, 18, list(range(12, 19))),
        (10, 10, [10]),
    ],
)
def test_parametrized_clipping_boundaries(
    region_start: int, region_end: int, expected_positions: list[int]
) -> None:
    """Inclusive-inclusive clipping holds at the edges."""
    block_gt = np.array([[[0, 0]]], dtype=np.int8)
    cs = _make_callset(pos=[10], end=[20], numalt=[1], is_snp=[0], gt=block_gt)
    out = expand_blocks(cs, region_start=region_start, region_end=region_end)
    np.testing.assert_array_equal(out["variants/POS"], expected_positions)


################################################################################
# Validation paths driven through pixy.__main__.main on real fixtures
################################################################################


def test_gvcf_flag_on_all_sites_vcf_raises(
    pixy_out_dir: Path,
    ag1000_pop_path: Path,
    ag1000_vcf_path: Path,
) -> None:
    """``--gvcf`` against an all-sites VCF must abort with an actionable error."""
    with pytest.raises(ValueError, match="does not appear to be a GVCF"):
        run_pixy_helper(
            pixy_out_dir=pixy_out_dir,
            stats=["pi"],
            window_size=10000,
            vcf_path=ag1000_vcf_path,
            populations_path=ag1000_pop_path,
            gvcf=True,
        )


def test_gvcf_input_without_flag_raises(
    pixy_out_dir: Path,
    ag1000_pop_path: Path,
    ag1000_gvcf_path: Path,
) -> None:
    """A GVCF without ``--gvcf`` must abort and point the user at the flag."""
    with pytest.raises(ValueError, match="--gvcf"):
        run_pixy_helper(
            pixy_out_dir=pixy_out_dir,
            stats=["pi"],
            window_size=10000,
            vcf_path=ag1000_gvcf_path,
            populations_path=ag1000_pop_path,
        )


def test_gvcf_and_wisp_bed_combination_rejected(
    pixy_out_dir: Path,
    ag1000_pop_path: Path,
    ag1000_gvcf_path: Path,
    tmp_path: Path,
) -> None:
    """``--gvcf`` and ``--wisp_bed`` are mutually exclusive; the combination errors out."""
    # The wisp BED is parsed before --gvcf compatibility is checked, but the --gvcf +
    # --wisp_bed rejection fires first in check_and_validate_args, so we just need any
    # path here — the validator never opens it.
    fake_wisp = tmp_path / "fake.bed.gz"
    fake_wisp.write_bytes(b"")
    with pytest.raises(ValueError, match="incompatible with --wisp_bed"):
        run_pixy_helper(
            pixy_out_dir=pixy_out_dir,
            stats=["pi"],
            window_size=10000,
            vcf_path=ag1000_gvcf_path,
            populations_path=ag1000_pop_path,
            gvcf=True,
            wisp_bed_path=fake_wisp,
        )


def test_gvcf_recognised_log_message(
    pixy_out_dir: Path,
    ag1000_pop_path: Path,
    ag1000_gvcf_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """When --gvcf is set on a real GVCF, pixy logs that the input was recognised."""
    caplog.set_level(logging.INFO)
    run_pixy_helper(
        pixy_out_dir=pixy_out_dir,
        stats=["pi"],
        window_size=10000,
        vcf_path=ag1000_gvcf_path,
        populations_path=ag1000_pop_path,
        gvcf=True,
        chromosomes="X",
    )
    assert any("Input recognised as a GVCF" in rec.getMessage() for rec in caplog.records)


################################################################################
# Symbolic ALT alleles (<NON_REF>, <*>)
################################################################################


def test_strip_symbolic_alts_recomputes_numalt_and_is_snp() -> None:
    """Symbolic ALTs are not real alleles: they must not count toward numalt or block is_snp."""
    gt = np.array(
        [
            [[0, 0], [0, 0]],  # single-site block: A -> <NON_REF>
            [[0, 1], [0, 2]],  # SNP: A -> G,<NON_REF>; one call to the symbolic allele
            [[0, 1], [0, 0]],  # plain SNP, untouched
            [[0, 1], [0, 0]],  # indel with a symbolic allele stays a non-SNP
            [[0, 1], [1, 2]],  # multiallelic SNP: A -> G,T,<*>
        ],
        dtype=np.int8,
    )
    cs: Dict[str, NDArray] = {
        "variants/REF": np.array(["A", "A", "A", "AC", "A"], dtype=object),
        "variants/ALT": np.array(
            [
                ["<NON_REF>", "", "", ""],
                ["G", "<NON_REF>", "", ""],
                ["G", "", "", ""],
                ["A", "<NON_REF>", "", ""],
                ["G", "T", "<*>", ""],
            ],
            dtype=object,
        ),
        "variants/numalt": np.array([1, 2, 1, 2, 3], dtype=np.int32),
        "variants/is_snp": np.array([False, False, True, False, False]),
        "calldata/GT": gt,
    }

    out = strip_symbolic_alts(cs)

    np.testing.assert_array_equal(out["variants/numalt"], [0, 1, 1, 1, 2])
    np.testing.assert_array_equal(out["variants/is_snp"], [False, True, True, False, True])
    # the call to the symbolic allele becomes missing; real alleles are untouched
    np.testing.assert_array_equal(out["calldata/GT"][1], [[0, 1], [0, -1]])
    np.testing.assert_array_equal(out["calldata/GT"][4], [[0, 1], [1, 2]])


def _write_indexed_vcf(path: Path, body: str) -> Path:
    header = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        '##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
    )
    path.write_text(header + body)
    subprocess.run(["bgzip", "-f", str(path)], check=True)
    gz = path.with_name(path.name + ".gz")
    subprocess.run(["tabix", "-p", "vcf", str(gz)], check=True)
    return gz


@pytest.mark.skipif(
    shutil.which("bgzip") is None or shutil.which("tabix") is None,
    reason="bgzip and tabix are required to build the synthetic GVCF",
)
def test_pixy_gvcf_gatk_style_records_match_all_sites(tmp_path: Path) -> None:
    """
    GATK-style GVCF records must give the same output as the equivalent all-sites VCF.

    Covers the record shapes real GATK GVCFs contain but the collapsed fixture does not:
    single-site blocks (``END == POS``), ``<NON_REF>`` records without ``END``, and variant
    records that carry a trailing ``<NON_REF>`` allele.
    """
    ref = "0/0\t0/0\t0/0\t0/0"
    all_sites = "".join(f"chr1\t{pos}\t.\tA\t.\t.\tPASS\t.\tGT\t{ref}\n" for pos in range(1, 7))
    all_sites += "chr1\t7\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\t1/1\t0/1\n"
    all_sites += "chr1\t8\t.\tA\tG,T\t.\tPASS\t.\tGT\t0/1\t0/2\t1/1\t0/0\n"
    all_sites += f"chr1\t9\t.\tA\t.\t.\tPASS\t.\tGT\t{ref}\n"

    gvcf = f"chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=5\tGT\t{ref}\n"
    gvcf += f"chr1\t6\t.\tA\t<NON_REF>\t.\tPASS\tEND=6\tGT\t{ref}\n"
    gvcf += "chr1\t7\t.\tA\tG,<NON_REF>\t.\tPASS\t.\tGT\t0/1\t0/0\t1/1\t0/1\n"
    gvcf += "chr1\t8\t.\tA\tG,T,<NON_REF>\t.\tPASS\t.\tGT\t0/1\t0/2\t1/1\t0/0\n"
    gvcf += f"chr1\t9\t.\tA\t<NON_REF>\t.\tPASS\t.\tGT\t{ref}\n"

    all_sites_gz = _write_indexed_vcf(tmp_path / "all_sites.vcf", all_sites)
    gvcf_gz = _write_indexed_vcf(tmp_path / "gvcf.vcf", gvcf)
    pop_file = tmp_path / "pops.txt"
    pop_file.write_text("S1\tA\nS2\tA\nS3\tB\nS4\tB\n")

    all_sites_out = tmp_path / "all_sites_out"
    gvcf_out = tmp_path / "gvcf_out"
    all_sites_out.mkdir()
    gvcf_out.mkdir()

    common: Dict[str, Any] = dict(
        stats=["pi", "dxy", "fst", "watterson_theta", "tajima_d"],
        window_size=100,
        populations_path=pop_file,
        fst_type="hudson",
        include_multiallelic_snps=True,
    )
    run_pixy_helper(pixy_out_dir=all_sites_out, vcf_path=all_sites_gz, **common)
    run_pixy_helper(pixy_out_dir=gvcf_out, vcf_path=gvcf_gz, gvcf=True, **common)

    for name in (
        "pixy_pi.txt",
        "pixy_dxy.txt",
        "pixy_fst.txt",
        "pixy_watterson_theta.txt",
        "pixy_tajima_d.txt",
    ):
        assert_files_are_consistent(gvcf_out / name, all_sites_out / name)

    # guard against both runs being trivially empty: 9 sites, 2 of them variable
    pi_rows = (gvcf_out / "pixy_pi.txt").read_text().splitlines()[1:]
    assert all(row.split("\t")[5] == "9" for row in pi_rows)
    assert all(float(row.split("\t")[4]) > 0 for row in pi_rows)


@pytest.mark.skipif(
    shutil.which("bgzip") is None or shutil.which("tabix") is None,
    reason="bgzip and tabix are required to build the synthetic GVCF",
)
def test_pixy_gvcf_fst_only_keeps_variants_with_symbolic_allele(tmp_path: Path) -> None:
    """The FST-only streaming read path must also see through a trailing ``<NON_REF>``."""
    gvcf = "chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=6\tGT\t0/0\t0/0\t0/0\t0/0\n"
    gvcf += "chr1\t7\t.\tA\tG,<NON_REF>\t.\tPASS\t.\tGT\t0/1\t0/0\t1/1\t0/1\n"
    gvcf_gz = _write_indexed_vcf(tmp_path / "gvcf.vcf", gvcf)
    pop_file = tmp_path / "pops.txt"
    pop_file.write_text("S1\tA\nS2\tA\nS3\tB\nS4\tB\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    run_pixy_helper(
        pixy_out_dir=out_dir,
        stats=["fst"],
        window_size=100,
        vcf_path=gvcf_gz,
        populations_path=pop_file,
        fst_type="hudson",
        gvcf=True,
    )

    fst_rows = (out_dir / "pixy_fst.txt").read_text().splitlines()[1:]
    assert len(fst_rows) == 1
    assert fst_rows[0].split("\t")[6] == "1"  # no_snps
