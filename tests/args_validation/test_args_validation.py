import argparse
from pathlib import Path

import pytest

from pixy.args_validation import PixyArgs
from pixy.args_validation import check_and_validate_args


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
