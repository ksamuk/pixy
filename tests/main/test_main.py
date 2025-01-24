import filecmp
import logging
import os
import shutil
from pathlib import Path
from typing import List
from typing import Optional
from unittest.mock import patch

import pytest

from tests.conftest import files_are_consistent
from tests.conftest import run_pixy_helper

################################################################################
# Tests for pixy.main(): missing, invalid, or conflicting arguments
################################################################################


@pytest.mark.regression
@pytest.mark.parametrize(
    "bed_str, window_size, chromosomes, interval_start, interval_end, expected_error_msg",
    [
        (
            "test_bed_path",
            10000,
            "all",
            None,
            None,
            "--interval_start, --interval_end, and --window_size",
        ),  # too many args: bed_path and window_size
        (
            "test_bed_path",
            10000,
            "all",
            1,
            20,
            "--interval_start, --interval_end, and --window_size",
        ),  # too many args: bed_path, window_size, interval_start, and interval_stop
        (
            None,
            None,
            "all",
            None,
            None,
            "In the absence of a BED file",
        ),  # no bed file and no window_size
        (
            None,
            1,
            "all",
            1,
            1,
            "--interval_start and --interval_end are not valid when calculating over multiple",
        ),  # can't specify an interval over multiple chromosomes
        (
            None,
            1,
            "X",
            1,
            None,
            "When specifying an interval,",
        ),  # too few args: only `interval_start` provided
        (
            None,
            1,
            "X",
            None,
            1,
            "When specifying an interval,",
        ),  # too few args: only `interval_end` provided
    ],
)
def test_missing_or_conflicting_args(
    bed_str: Optional[str],
    window_size: Optional[int],
    chromosomes: str,
    interval_start: Optional[int],
    interval_end: Optional[int],
    expected_error_msg: str,
    pixy_out_dir: Path,
    ag1000_vcf_path: Path,
    request: pytest.FixtureRequest,
    ag1000_pop_path: Path,
) -> None:
    """
    Assert that we raise an exception when args are missing or if there are conflicting args.

    `vcf_path` and `populations_path` stay the same here and are tested separately elsewhere.
    """
    if bed_str is not None:
        bed_path = request.getfixturevalue(bed_str)
    else:
        bed_path = None
    with pytest.raises(ValueError, match=expected_error_msg):
        run_pixy_helper(
            pixy_out_dir=pixy_out_dir,
            stats=["pi", "fst", "dxy"],
            bed_path=bed_path,
            chromosomes=chromosomes,
            interval_start=interval_start,
            interval_end=interval_end,
            window_size=window_size,
            vcf_path=ag1000_vcf_path,
            populations_path=ag1000_pop_path,
        )


@pytest.mark.regression
def test_vcf_missing_index(
    pixy_out_dir: Path,
    tmp_path: Path,
    ag1000_pop_path: Path,
    ag1000_vcf_path: Path,
) -> None:
    """Assert that we raise an exception when missing .tbi index."""
    missing_index_vcf_path: Path = tmp_path / "ag1000_pixy_test.vcf.gz"
    shutil.copy(ag1000_vcf_path, missing_index_vcf_path)
    with pytest.raises(ValueError, match="The vcf is not indexed."):
        run_pixy_helper(
            pixy_out_dir=pixy_out_dir,
            stats=["pi", "fst", "dxy"],
            window_size=10000,
            vcf_path=missing_index_vcf_path,
            populations_path=ag1000_pop_path,
        )


@pytest.mark.regression
def test_missing_tabix_path_raises(
    pixy_out_dir: Path, ag1000_pop_path: Path, ag1000_vcf_path: Path
) -> None:
    """Assert that we raise an exception when `tabix` is not on the path."""
    with patch("pixy.args_validation.shutil.which", return_value=None):
        with pytest.raises(ValueError, match="`tabix` is not installed"):
            run_pixy_helper(
                pixy_out_dir=pixy_out_dir,
                stats=["pi", "fst", "dxy"],
                window_size=10000,
                vcf_path=ag1000_vcf_path,
                populations_path=ag1000_pop_path,
            )


@pytest.mark.regression
def test_missing_vcf_file_raises(tmp_path: Path, pixy_out_dir: Path, ag1000_pop_path: Path) -> None:
    """Assert that we raise an exception with an uncompressed VCF file."""
    vcf_path: Path = tmp_path / "uncompressed_vcf.vcf"

    with open(vcf_path, "w") as f:
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t\\INFO\tFORMAT")

    with pytest.raises(ValueError, match="The vcf is not compressed"):
        run_pixy_helper(
            pixy_out_dir=pixy_out_dir,
            stats=["pi", "fst", "dxy"],
            window_size=10000,
            vcf_path=vcf_path,
            populations_path=ag1000_pop_path,
        )


@pytest.mark.regression
def test_missing_pop_file_raises(pixy_out_dir: Path, ag1000_vcf_path: Path) -> None:
    """Assert that we raise an exception with a missing `populations_path`."""
    with pytest.raises(FileNotFoundError, match="The specified populations file"):
        run_pixy_helper(
            pixy_out_dir=pixy_out_dir,
            stats=["pi", "fst", "dxy"],
            window_size=10000,
            vcf_path=ag1000_vcf_path,
            populations_path=None,
        )


@pytest.mark.regression
def test_missing_chroms_in_vcf_raises(
    pixy_out_dir: Path, ag1000_vcf_path: Path, ag1000_pop_path: Path
) -> None:
    """Assert that we raise an exception when a given chromosome subset is not found in the VCF."""
    with pytest.raises(ValueError, match="The following chromosomes"):
        run_pixy_helper(
            pixy_out_dir=pixy_out_dir,
            stats=["pi", "fst", "dxy"],
            window_size=10000,
            vcf_path=ag1000_vcf_path,
            populations_path=ag1000_pop_path,
            chromosomes="23",  # 23 is not in the VCF
        )


@pytest.mark.regression
def test_vcf_no_invariant_sites_raises(
    ag1000_pop_path: Path,
    test_bed_path: Path,
    pixy_out_dir: Path,
    ag1000_invariant_vcf_path: Path,
) -> None:
    """Assert that we raise an error when a VCF contains no variable sites."""
    with pytest.raises(ValueError, match="The provided VCF appears to contain no invariant"):
        run_pixy_helper(
            pixy_out_dir=pixy_out_dir,
            stats=["pi", "fst", "dxy"],
            vcf_path=ag1000_invariant_vcf_path,
            populations_path=ag1000_pop_path,
            bed_path=test_bed_path,
        )


################################################################################
# Tests for pixy.main(): invalid or malformed arguments
################################################################################


@pytest.mark.regression
@pytest.mark.parametrize(
    "malformed_populations_input, expected_error_msg",
    [
        (
            """ERS224670\tBFS\nERS224248\tBFS\nERS224089\t""",
            "The specified populations file contains missing data",
        ),  # row missing 2nd column
        (
            """ERS224670\tBFS\nERS224248\tBFS\n""",
            "Calculation of fst and/or dxy requires at least two",
        ),  # only 1 population here; `pixy` requires 2 populations for calculation
        (
            """INVALID_SAMPLE_ID_001\tXFS\nINVALID_SAMPLE_ID_002\tXFS\n""",
            "The following samples are listed",
        ),  # samples in populations.txt are not in VCF
    ],
)
def test_malformed_populations_files_raises(
    malformed_populations_input: str,
    expected_error_msg: str,
    ag1000_vcf_path: Path,
    pixy_out_dir: Path,
    tmp_path: Path,
) -> None:
    """Raise an exception when given a malformed `population.txt` file."""
    malformed_pop_path: Path = tmp_path / "malformed_pop_file.txt"
    with open(malformed_pop_path, "w") as f:
        f.write(malformed_populations_input)
    with pytest.raises(ValueError, match=expected_error_msg):
        run_pixy_helper(
            pixy_out_dir=pixy_out_dir,
            stats=["pi", "fst", "dxy"],
            window_size=10000,
            vcf_path=ag1000_vcf_path,
            populations_path=malformed_pop_path,
        )


@pytest.mark.regression
@pytest.mark.parametrize(
    "malformed_sites_input, expected_error_msg",
    [
        (
            """X\t1\nX\t2\nX\t""",
            "The specified sites file contains missing data",
        ),  # row missing position
        (
            """X\n1\n2\n3\n""",
            "Too many columns specified: expected 2 and found 1",
        ),  # row missing position
    ],
)
def test_malformed_sites_file(
    malformed_sites_input: str,
    expected_error_msg: str,
    ag1000_vcf_path: Path,
    pixy_out_dir: Path,
    tmp_path: Path,
    ag1000_pop_path: Path,
) -> None:
    """
    Raise an exception when given a malformed `sites.txt` file.

    `@pytest.mark.parametrize` is used to support addition of different errors in the future.
    """
    malformed_site_path: Path = tmp_path / "malformed_sites_file.txt"
    with open(malformed_site_path, "w") as f:
        f.write(malformed_sites_input)
    with pytest.raises(ValueError, match=expected_error_msg):
        run_pixy_helper(
            pixy_out_dir=pixy_out_dir,
            stats=["pi", "fst", "dxy"],
            window_size=10000,
            vcf_path=ag1000_vcf_path,
            populations_path=ag1000_pop_path,
            sites_path=malformed_site_path,
        )


@pytest.mark.regression
@pytest.mark.parametrize(
    "malformed_bed_input, expected_error_msg",
    [
        (
            """X\t1\t20\nX\t2\n""",
            "The specified BED file contains missing data",
        ),  # row missing position
    ],
)
def test_malformed_bed_file(
    malformed_bed_input: str,
    expected_error_msg: str,
    pixy_out_dir: Path,
    tmp_path: Path,
    ag1000_vcf_path: Path,
    ag1000_pop_path: Path,
) -> None:
    """
    Raise an exception when given a malformed `intervals.bed` file.

    `@pytest.mark.parametrize` is used to support addition of different errors in the future.
    """
    bed_path: Path = tmp_path / "malformed.bed"

    with open(bed_path, "w") as f:
        f.write(malformed_bed_input)
    with pytest.raises(ValueError, match=expected_error_msg):
        run_pixy_helper(
            pixy_out_dir=pixy_out_dir,
            stats=["pi", "fst", "dxy"],
            vcf_path=ag1000_vcf_path,
            populations_path=ag1000_pop_path,
            bed_path=bed_path,
        )


################################################################################
# Tests for pixy.main(): warnings
################################################################################


@pytest.mark.regression
def test_vcf_bed_chrom_difference_warns(
    tmp_path_factory: pytest.TempPathFactory,
    pixy_out_dir: Path,
    ag1000_pop_path: Path,
    ag1000_vcf_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Assert that chromosomes in .bed but not .vcf yields a warning."""
    temp_dir = tmp_path_factory.mktemp("warn_files")
    bed_path: Path = temp_dir / "not_in_vcf.bed"
    caplog.set_level(logging.WARNING)
    with open(bed_path, "w") as f:
        f.write("""X\t1\t20\nX\t2\t30\n""")  # X is in the VCF
        f.write("""23\t1\t20\n23\t2\t30\n""")  # 23 is not in the VCF

    run_pixy_helper(
        pixy_out_dir=pixy_out_dir,
        stats=["pi", "fst", "dxy"],
        bed_path=bed_path,
        vcf_path=ag1000_vcf_path,
        populations_path=ag1000_pop_path,
    )
    assert (
        "The following chromosomes are in the BED file but do not occur in the VCF "
        "and will be ignored: ['23']" in caplog.messages
    )


@pytest.mark.regression
def test_bypass_invariant_check_warns(
    pixy_out_dir: Path,
    ag1000_pop_path: Path,
    ag1000_vcf_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Assert that `--bypass_invariant_check flag yields a warning."""
    run_pixy_helper(
        pixy_out_dir=pixy_out_dir,
        stats=["pi", "fst", "dxy"],
        window_size=10000,
        vcf_path=ag1000_vcf_path,
        populations_path=ag1000_pop_path,
        bypass_invariant_check="yes",
    )
    assert "EXTREME WARNING: --bypass_invariant_check is set to 'yes'" in caplog.text


################################################################################
# Tests for pixy.main(): valid inputs and expected results
################################################################################

#############################
# Tests for VCF file indexes
#############################


@pytest.mark.regression
def test_pixy_csi_index(
    tmp_path: Path,
    ag1000_pop_path: Path,
    ag1000_vcf_path: Path,
    ag1000_csi_path: Path,
    pixy_out_dir: Path,
    expected_outputs: Path,
) -> None:
    """
    Assert that a VCF can have either a `.tbi` or a `.csi` index with valid inputs.

    The outputs with a `.csi` index should match the outputs of the `.tbi` index.

    NB, we copy `ag1000_pixy_test.vcf.gz` and `ag1000_pixy_test.vcf.gz.csi` into a clean directory
    so that we can be confident there is no interference from the pre-existing `.tbi` file.
    """
    vcf_path: Path = tmp_path / "ag1000_pixy_test.vcf.gz"
    csi_path: Path = tmp_path / "ag1000_pixy_test.vcf.gz.csi"
    shutil.copy(ag1000_vcf_path, vcf_path)
    shutil.copy(ag1000_csi_path, csi_path)

    run_pixy_helper(
        pixy_out_dir=pixy_out_dir,
        window_size=10000,
        vcf_path=vcf_path,
        populations_path=ag1000_pop_path,
        stats=["pi", "dxy", "fst"],
        output_prefix="pixy",
    )

    expected_out_files: List[Path] = [
        Path("pixy_dxy.txt"),
        Path("pixy_fst.txt"),
        Path("pixy_pi.txt"),
    ]
    # this run of `pixy` should match the run using the same inputs and a `.tbi` index
    for file in expected_out_files:
        generated_data_path: Path = pixy_out_dir / file
        exp_data_path: Path = expected_outputs / "baseline" / file
        assert generated_data_path.exists()

        assert filecmp.cmp(generated_data_path, exp_data_path)


#######################################
# Tests for output formatting/creation
#######################################


@pytest.mark.regression
@pytest.mark.parametrize(
    "output_prefix, stats_requested, expected_files",
    [
        (
            None,
            ["pi", "fst", "dxy"],
            ["pixy_pi.txt", "pixy_fst.txt", "pixy_dxy.txt"],
        ),  # default prefix
        (
            "test",
            ["pi", "fst", "dxy"],
            ["test_pi.txt", "test_fst.txt", "test_dxy.txt"],
        ),  # non-default prefix
        ("test", ["pi"], ["test_pi.txt"]),  # only pi
        ("test", ["pi", "fst"], ["test_pi.txt", "test_fst.txt"]),  # both pi and fst
    ],
)
def test_pixy_output_creation(
    output_prefix: Optional[str],
    stats_requested: List[str],
    expected_files: List[str],
    pixy_out_dir: Path,
    ag1000_pop_path: Path,
    ag1000_vcf_path: Path,
) -> None:
    """Given different stats and file prefixes, assert that output creation is as expected."""
    run_pixy_helper(
        pixy_out_dir=pixy_out_dir,
        window_size=10000,
        vcf_path=ag1000_vcf_path,
        populations_path=ag1000_pop_path,
        stats=stats_requested,
        output_prefix=output_prefix,
    )

    for file in expected_files:
        full_path: Path = pixy_out_dir / Path(file)
        assert full_path.exists()

    all_possible_files: List[str] = ["pi", "fst", "dxy"]
    # make sure we do not produce files we did not ask for
    unexpected_files: List[str] = list(set(all_possible_files) - set(stats_requested))
    for file in unexpected_files:
        full_path = pixy_out_dir / Path(f"{output_prefix}_{file}.txt")
        assert not os.path.exists(full_path)


################################################################################
# Tests for pixy.main(): baseline input
################################################################################


@pytest.mark.regression
def test_pixy_main_valid_inputs(
    pixy_out_dir: Path,
    expected_outputs: Path,
    ag1000_pop_path: Path,
    ag1000_vcf_path: Path,
    capsys: pytest.CaptureFixture,
) -> None:
    """
    Given specific input data, assert that outputs do not change.

    Uses `filecmp` library to compare 2 files without opening them and reading line-by-line.
    `filecmp.cmp` returns True if 2 files are equal.
    """
    run_pixy_helper(
        pixy_out_dir=pixy_out_dir,
        stats=["pi", "fst", "dxy"],
        window_size=10000,
        vcf_path=ag1000_vcf_path,
        populations_path=ag1000_pop_path,
    )
    captured = capsys.readouterr()
    assert "Data set contains 2 population(s), 2 chromosome(s), and 36 sample(s)" in captured.out

    expected_out_files: List[Path] = [
        Path("pixy_dxy.txt"),
        Path("pixy_fst.txt"),
        Path("pixy_pi.txt"),
    ]
    for file in expected_out_files:
        generated_data_path: Path = pixy_out_dir / file
        exp_data_path: Path = expected_outputs / "baseline" / file
        assert generated_data_path.exists()

        assert filecmp.cmp(generated_data_path, exp_data_path)


################################################################################
# Tests for pixy.main(): limited/single sites
################################################################################


@pytest.mark.regression
@pytest.mark.parametrize(
    "window_size, sites_file, stats, output_prefix",
    [
        (
            10000,
            "test_small_sites_path",
            ["pi", "dxy"],
            "limited_sites",
        ),
        (1, None, ["fst"], "single_sites"),
    ],
)
def test_pixy_limited_sites(
    pixy_out_dir: Path,
    expected_outputs: Path,
    ag1000_pop_path: Path,
    ag1000_vcf_path: Path,
    request: pytest.FixtureRequest,
    window_size: int,
    sites_file: str,
    stats: List[str],
    output_prefix: str,
) -> None:
    """Run `pixy` with either limited or single sites and compare to ground-truth data."""
    sites_path: Optional[Path] = (
        request.getfixturevalue(sites_file) if sites_file is not None else None
    )
    run_pixy_helper(
        pixy_out_dir=pixy_out_dir,
        stats=stats,
        sites_path=sites_path if sites_file is not None else None,
        window_size=window_size,
        vcf_path=ag1000_vcf_path,
        populations_path=ag1000_pop_path,
        output_prefix=output_prefix,
        debug=True,
        cores=1,
    )

    expected_out_files: List[Path] = [Path(f"{output_prefix}_{stat}.txt") for stat in stats]

    for file in expected_out_files:
        generated_data_path: Path = pixy_out_dir / file
        exp_data_path: Path = expected_outputs / output_prefix / file

        assert generated_data_path.exists()
        assert files_are_consistent(generated_data_path, exp_data_path)


################################################################################
# Tests for pixy.main(): limited BED file
################################################################################
@pytest.mark.regression
def test_pixy_limited_bed_file(
    pixy_out_dir: Path,
    test_three_chrom_bed_path: Path,
    expected_outputs: Path,
    ag1000_pop_path: Path,
    ag1000_vcf_path: Path,
) -> None:
    """
    Test that `pixy` produces deterministic stats given static input data.

    Testing a small BED file in isolation (e.g., separate from `sites` file).
    """
    run_pixy_helper(
        pixy_out_dir=pixy_out_dir,
        bed_path=test_three_chrom_bed_path,
        stats=["dxy", "pi", "fst"],
        vcf_path=ag1000_vcf_path,
        populations_path=ag1000_pop_path,
        output_prefix="limited_bed",
        debug=True,
        cores=1,
    )

    expected_out_files: List[Path] = [
        Path("limited_bed_dxy.txt"),
        Path("limited_bed_pi.txt"),
        Path("limited_bed_fst.txt"),
    ]
    for file in expected_out_files:
        generated_data_path: Path = pixy_out_dir / file
        exp_data_path: Path = expected_outputs / "limited_bed" / file

        assert generated_data_path.exists()
        assert files_are_consistent(generated_data_path, exp_data_path)


################################################################################
# Tests for pixy.main(): limited sites and BED file
################################################################################


@pytest.mark.regression
def test_pixy_limited_sites_bed(
    pixy_out_dir: Path,
    test_three_chrom_bed_path: Path,
    test_small_sites_path: Path,
    expected_outputs: Path,
    ag1000_pop_path: Path,
    ag1000_vcf_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """
    Test that `pixy` produces deterministic stats given static input data.

    Input data here is a small sites file and a small BED file.
    """
    run_pixy_helper(
        pixy_out_dir=pixy_out_dir,
        bed_path=test_three_chrom_bed_path,
        sites_path=test_small_sites_path,
        stats=["dxy", "pi"],
        vcf_path=ag1000_vcf_path,
        populations_path=ag1000_pop_path,
        output_prefix="limited_sites_and_bed",
        debug=True,
        cores=1,
    )
    expected_warnings: List[str] = [
        "The following chromosomes are in the BED file but do not occur in the VCF "
        "and will be ignored: ['5']",
        "The following chromosomes occur in the sites file but do not occur in the VCF "
        "and will be ignored: ['B']",
    ]

    assert set(expected_warnings) == set(caplog.messages)

    expected_out_files: List[Path] = [
        Path("limited_sites_and_bed_dxy.txt"),
        Path("limited_sites_and_bed_pi.txt"),
    ]
    for file in expected_out_files:
        generated_data_path: Path = pixy_out_dir / file
        exp_data_path: Path = expected_outputs / "limited_sites_and_bed" / file
        assert generated_data_path.exists()
        assert files_are_consistent(generated_data_path, exp_data_path)


###############################################################################
# Tests for pixy.main(): Hudson's FST
################################################################################


@pytest.mark.regression
def test_pixy_hudson_fst(
    pixy_out_dir: Path,
    expected_outputs: Path,
    ag1000_pop_path: Path,
    ag1000_vcf_path: Path,
) -> None:
    """Test that pixy produces deterministic Hudson FST stats with known input data."""
    run_pixy_helper(
        pixy_out_dir=pixy_out_dir,
        stats=["fst"],
        window_size=10000,
        vcf_path=ag1000_vcf_path,
        populations_path=ag1000_pop_path,
        output_prefix="hudson",
        debug=True,
        cores=1,
        fst_type="hudson",
    )

    expected_out_files: List[Path] = [
        Path("hudson_fst.txt"),
    ]
    for file in expected_out_files:
        generated_data_path: Path = pixy_out_dir / file
        exp_data_path: Path = expected_outputs / "hudson_fst" / file
        assert generated_data_path.exists()
        assert files_are_consistent(generated_data_path, exp_data_path)
