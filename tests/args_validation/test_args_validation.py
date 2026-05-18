import argparse
import shutil
import subprocess
from pathlib import Path

import pytest

from pixy.args_validation import PixyArgs
from pixy.args_validation import _ploidy_from_gt
from pixy.args_validation import check_and_validate_args
from pixy.args_validation import infer_ploidy_per_contig


@pytest.mark.parametrize("bypass_variant_check", [True, False])
def test_check_and_validate_args(
    ag1000_vcf_path: Path,
    ag1000_pop_path: Path,
    bypass_variant_check: str,
    tmp_path: Path,
) -> None:
    """Assert that provided CLI args are transformed as expected by `check_and_validate_args`."""
    args = argparse.Namespace()
    args.vcf = str(ag1000_vcf_path)
    args.stats = ["dxy", "fst", "pi"]
    args.populations = ag1000_pop_path
    args.output_folder = str(tmp_path / "output")
    args.temp_file = str(tmp_path / "temp_file.txt")
    args.chromosomes = "X"
    args.n_cores = 1
    args.bypass_invariant_check = bypass_variant_check
    args.include_multiallelic_snps = False
    args.fst_type = "wc"
    args.output_prefix = "test"
    args.chunk_size = 100000
    args.bed_file = None
    args.window_size = 1
    args.interval_start = None
    args.interval_end = None
    args.sites_file = None

    generated_pixy_args: PixyArgs = check_and_validate_args(args)
    if bypass_variant_check:
        assert generated_pixy_args.bypass_invariant_check
    else:
        assert not generated_pixy_args.bypass_invariant_check
    # ag1000 test data is uniformly diploid on chromosome X
    assert generated_pixy_args.ploidy_map == {"X": 2}


################################################################################
# Tests for per-contig ploidy inference
################################################################################


def test_ploidy_from_gt_diploid() -> None:
    """Diploid genotypes with `/` or `|` separators return ploidy 2."""
    assert _ploidy_from_gt("0/0") == 2
    assert _ploidy_from_gt("0|1") == 2
    assert _ploidy_from_gt("./.") == 2


def test_ploidy_from_gt_haploid() -> None:
    """A single allele (with or without missing partner) returns ploidy 1."""
    assert _ploidy_from_gt("0") == 1
    assert _ploidy_from_gt("1") == 1
    assert _ploidy_from_gt(".") == 1


def test_ploidy_from_gt_polyploid() -> None:
    """Higher-ploidy genotypes return the correct count."""
    assert _ploidy_from_gt("0/0/0") == 3
    assert _ploidy_from_gt("0|1|0|1") == 4


def _bgzip_and_tabix(tmp_path: Path, vcf_text: str) -> Path:
    """Write `vcf_text` to a `.vcf` file, bgzip it, and tabix-index it. Returns the .vcf.gz path."""
    vcf_path = tmp_path / "mixed_ploidy.vcf"
    vcf_path.write_text(vcf_text)
    subprocess.run(["bgzip", "-f", str(vcf_path)], check=True)
    bgz_path = vcf_path.with_suffix(".vcf.gz")
    subprocess.run(["tabix", "-p", "vcf", str(bgz_path)], check=True)
    return bgz_path


@pytest.mark.skipif(
    shutil.which("bgzip") is None or shutil.which("tabix") is None,
    reason="bgzip and tabix are required to build the mixed-ploidy fixture VCF",
)
def test_infer_ploidy_per_contig_mixed(tmp_path: Path) -> None:
    """`infer_ploidy_per_contig` returns per-contig ploidy for VCFs with variable ploidy."""
    vcf_text = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1>\n"
        "##contig=<ID=chrX>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
        "chr1\t1\t.\tA\t.\t.\tPASS\t.\tGT\t0/0\t0/1\n"
        "chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\t1/1\n"
        "chrX\t1\t.\tA\t.\t.\tPASS\t.\tGT\t0\t1\n"
        "chrX\t2\t.\tG\tA\t.\tPASS\t.\tGT\t1\t0\n"
    )
    bgz_path = _bgzip_and_tabix(tmp_path, vcf_text)

    result = infer_ploidy_per_contig(str(bgz_path), ["chr1", "chrX"])
    assert result == {"chr1": 2, "chrX": 1}


@pytest.mark.skipif(
    shutil.which("bgzip") is None or shutil.which("tabix") is None,
    reason="bgzip and tabix are required to build the mixed-ploidy fixture VCF",
)
def test_infer_ploidy_per_contig_missing_contig_raises(tmp_path: Path) -> None:
    """A contig that has no records in the VCF raises a clear error."""
    vcf_text = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr1\t1\t.\tA\t.\t.\tPASS\t.\tGT\t0/0\n"
    )
    bgz_path = _bgzip_and_tabix(tmp_path, vcf_text)

    with pytest.raises(RuntimeError, match="No genotype records found"):
        infer_ploidy_per_contig(str(bgz_path), ["chr1", "chrMissing"])
