import sys
from pathlib import Path
from typing import List
from typing import Optional
from unittest.mock import patch

import pandas as pd
import pytest

from pixy.__main__ import main

################################################################################
# Fixtures for testing input + output locations
################################################################################


@pytest.fixture
def datadir() -> Path:
    """Path to the input data directory for testing. `pixy` consumes files from this directory."""
    return Path(__file__).parent / "main" / "data"


@pytest.fixture
def expected_outputs() -> Path:
    """
    Path to a directory containing expected outputs for testing comparisons.

    `pixy` compares testing output to the ground-truth files in this directory.
    """
    return Path(__file__).parent / "main" / "expected_outputs"


################################################################################
# Fixtures for invalid input testing data
################################################################################


@pytest.fixture()
def ag1000_invariant_vcf_path(datadir: Path) -> Path:
    """
    Path to ag1000 invariant `vcf_path`.

    This VCF does not contain any variation and is considered a suboptimal input to `pixy`.
    """
    return datadir / Path("ag1000_pixy_invariant_test.vcf.gz")


################################################################################
# Fixtures for valid input testing data
################################################################################


#################
# VCF FILES ###
#################
@pytest.fixture()
def ag1000_vcf_path(datadir: Path) -> Path:
    """Path to ag1000 VCF."""
    return datadir / "ag1000_pixy_test.vcf.gz"


@pytest.fixture()
def ag1000_csi_path(datadir: Path) -> Path:
    """Path to ag1000 VCF CSI index."""
    return datadir / "ag1000_pixy_test.vcf.gz.csi"


@pytest.fixture()
def missing50_vcf_path(datadir: Path) -> Path:
    """Path to a simulated VCF that is known to be missing a small number of genotypes."""
    return datadir / "simulated_data_missing50p_genos.vcf.gz"


@pytest.fixture()
def missing5000_vcf_path(datadir: Path) -> Path:
    """Path to a simulated VCF that is known to be missing a larger number of genotypes."""
    return datadir / "simulated_data_missing5000_sites.vcf.gz"


########################
# POPULATION FILES ###
########################
@pytest.fixture()
def ag1000_pop_path(datadir: Path) -> Path:
    """Path to a valid test `populations` file (derived from Anopheles gambiae 1000 Genomes)."""
    return datadir / "ag1000_populations_file.txt"


@pytest.fixture()
def simulated_pop_pi_path(datadir: Path) -> Path:
    """Path to a simulated test `populations` file representing 1 population."""
    return datadir / "simulated_populations_pi.txt"


@pytest.fixture()
def simulated_pop_dxy_path(datadir: Path) -> Path:
    """Path to a simulated test `populations` file with 2 populations."""
    return datadir / "simulated_populations_dxy.txt"


#################
# BED FILES ###
#################
@pytest.fixture()
def test_bed_path(datadir: Path) -> Path:
    """Path to a valid test BED file containing 2 chromosomes (X, 1)."""
    return datadir / "test_regions_1.bed"


@pytest.fixture()
def test_three_chrom_bed_path(datadir: Path) -> Path:
    """Path to a valid test BED file containing 3 chromosomes (X, 1, 5)."""
    return datadir / "test_regions_2.bed"


###################
# SITES FILES ###
###################
@pytest.fixture()
def test_sites_path(datadir: Path) -> Path:
    """Path to a valid test `sites` file (derived from Anopheles gambiae 1000 Genomes)."""
    return datadir / "test_sites_1.txt"


@pytest.fixture()
def test_small_sites_path(datadir: Path) -> Path:
    """Path to a smaller valid test `sites` file."""
    return datadir / "test_sites_2.txt"


################################################################################
# Shared fixtures: output directory paths, input files, etc.
################################################################################
@pytest.fixture()
def pixy_out_dir(tmp_path: Path) -> Path:
    """Temporary directory for `pixy` output."""
    return tmp_path / "output"


# NB: The complexity level here is fine, we're just unpacking a bunch of Optionals
def run_pixy_helper(  # noqa: C901
    pixy_out_dir: Path,
    stats: List[str],
    vcf_path: Path,
    populations_path: Optional[Path] = None,
    include_multiallelic_snps: bool = False,
    bypass_invariant_check: bool = False,
    window_size: Optional[int] = None,
    interval_start: Optional[int] = None,
    interval_end: Optional[int] = None,
    bed_path: Optional[Path] = None,
    chromosomes: Optional[str] = None,
    sites_path: Optional[Path] = None,
    output_prefix: Optional[str] = None,
    debug: bool = False,
    cores: Optional[int] = None,
    fst_type: Optional[str] = None,
) -> None:
    """
    Run `pixy` with the specified arguments.

    We don't do any error checking here because we want `pixy` to do it. E.g., if `--bed_file` is
    not given to `pixy`, it checks for `--window_size`. If that is not given, then it requires
    `--interval_start` and `--interval_end`. We only check for what is given and map to `pixy` args
     accordingly (we do not check if args are valid).
    """
    if chromosomes is None:
        chromosomes = "all"  # default pixy value

    test_args = [
        "pixy",
        "--stats",
        *stats,
        "--vcf",
        f"{vcf_path}",
        "--populations",
        f"{populations_path}",
        "--output_folder",
        f"{pixy_out_dir}",
        "--chromosomes",
        f"{chromosomes}",
    ]
    if window_size is not None:
        test_args.extend((["--window_size", f"{window_size}"]))

    if sites_path is not None:
        test_args.extend(["--sites_file", f"{sites_path}"])
    # NB: we deliberately permit this helper to run with any type of args for the next 3 flags
    # in order to test that `pixy` correctly raises an error in this case.

    if bed_path is not None:
        test_args.extend(["--bed_file", f"{bed_path}"])
    else:
        test_args.extend(["--bed_file"])

    if interval_start is not None:
        test_args.extend((["--interval_start", f"{interval_start}"]))
    else:
        test_args.extend((["--interval_start"]))

    if interval_end is not None:
        test_args.extend((["--interval_end", f"{interval_end}"]))
    else:
        test_args.extend((["--interval_end"]))

    if output_prefix is not None:
        test_args.extend((["--output_prefix", f"{output_prefix}"]))

    if debug is True:
        test_args.extend((["--debug"]))

    if cores is not None:
        test_args.extend((["--n_cores", f"{cores}"]))

    if fst_type is not None:
        test_args.extend((["--fst_type", f"{fst_type}"]))

    if bypass_invariant_check:
        test_args.extend(["--bypass_invariant_check"])

    if include_multiallelic_snps:
        test_args.extend(["--include_multiallelic_snps"])
    print(f"test_args: {test_args}")
    with patch.object(sys, "argv", test_args):
        main()


def assert_files_are_consistent(gen_file_path: Path, exp_file_path: Path) -> None:
    """
    Helper to display diff if two files are found to be different.

    The basename of the generated file is shown (e.g., `pixy_dxy`). The full path of the expected
    file is given to highlight which test case is failing.

    Args:
         exp_file_path: the path of the expected file (i.e., ground-truth)
         gen_file_path: the path of the generated file

    Raises:
        AssertionError, if the two files are not the same
    """
    if not files_are_consistent(gen_file_path, exp_file_path):
        diff_string: str = show_diff(gen_file_path, exp_file_path)
        raise AssertionError(f"Files differ: {gen_file_path.stem}, {exp_file_path}\n{diff_string}")


def show_diff(expected_file: Path, generated_file: Path) -> str:
    """
    Show the diff between an expected and generated file.

    Useful to examine why a regression test fails.

    Args:
         expected_file: the path of the expected file (i.e., ground-truth)
         generated_file: the path of the generated file

    """
    exp = pd.read_csv(expected_file, sep="\t")
    gen = pd.read_csv(generated_file, sep="\t")
    diff = exp.compare(gen, result_names=("Expected", "Generated"), align_axis=0)
    pretty_display: str = diff.to_string()
    return pretty_display


def files_are_consistent(gen_file_path: Path, exp_file_path: Path) -> bool:
    """
    Helper function to compare non-deterministic files generated by `pixy`.

    Checks that the headers are the same and the length of the files matches.
    Sorts the data to be deterministic before comparing to ground-truth data.

    Used in regression testing to compare specific rows in generated files for reproducibility.

    Raises:
        FileNotFoundError: if one or both files do not exist
        ValueError: if one of the provided paths is not a file
    Returns:
        True if lines in file match each other; False if there is a discrepancy
    """
    if not gen_file_path.exists() or not exp_file_path.exists():
        raise FileNotFoundError("One or both files do not exist.")

    if not gen_file_path.is_file() or not exp_file_path.is_file():
        raise ValueError("Both paths should be files.")

    with open(gen_file_path, "r") as generated_file, open(exp_file_path, "r") as expected_file:
        generated_data: List[str] = generated_file.readlines()
        expected_data: List[str] = expected_file.readlines()

        # lines of data should be the same
        if len(generated_data) != len(expected_data):
            return False

        generated_data.sort()
        expected_data.sort()

        for line1, line2 in zip(generated_data, expected_data):
            if line1 != line2:
                return False
    return True
