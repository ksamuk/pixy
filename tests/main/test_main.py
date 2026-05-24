import logging
import os
import shutil
from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Optional
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from pixy.calc import calc_tajima_d_stdev
from pixy.calc import deserialize_tajima_d_variant_counts
from tests.conftest import assert_files_are_consistent
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
        bypass_invariant_check=True,
    )
    assert "EXTREME WARNING: --bypass_invariant_check is set to True" in caplog.text


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
        # shutil.copy(generated_data_path, exp_data_path)
        assert_files_are_consistent(generated_data_path, exp_data_path)


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
        ("test", ["pi", "tajima_d"], ["test_pi.txt", "test_tajima_d.txt"]),  # include tajima_d
        (
            "test",
            ["pi", "watterson_theta"],
            ["test_pi.txt", "test_watterson_theta.txt"],
        ),  # include watterson_theta
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

    all_possible_files: List[str] = ["pi", "fst", "dxy", "watterson_theta", "tajima_d"]
    # make sure we do not produce files we did not ask for
    unexpected_files: List[str] = list(set(all_possible_files) - set(stats_requested))
    for file in unexpected_files:
        full_path = pixy_out_dir / Path(f"{output_prefix}_{file}.txt")
        assert not os.path.exists(full_path)


@pytest.mark.regression
def test_pixy_watterson_theta_aggregation_matches_direct_calculation(
    pixy_out_dir: Path,
    ag1000_pop_path: Path,
    ag1000_vcf_path: Path,
) -> None:
    """Chunk-aggregated Watterson's theta should match direct window calculation."""
    shared_args: Dict[str, Any] = {
        "pixy_out_dir": pixy_out_dir,
        "stats": ["watterson_theta"],
        "window_size": 10000,
        "vcf_path": ag1000_vcf_path,
        "populations_path": ag1000_pop_path,
        "chromosomes": "X",
        "cores": 1,
    }
    run_pixy_helper(
        **shared_args,
        output_prefix="watterson_direct",
        chunk_size=100000,
    )
    run_pixy_helper(
        **shared_args,
        output_prefix="watterson_aggregated",
        chunk_size=5000,
    )

    direct = pd.read_csv(pixy_out_dir / "watterson_direct_watterson_theta.txt", sep="\t")
    aggregated = pd.read_csv(pixy_out_dir / "watterson_aggregated_watterson_theta.txt", sep="\t")
    sort_columns = ["pop", "chromosome", "window_pos_1", "window_pos_2"]
    direct = direct.sort_values(sort_columns).reset_index(drop=True)
    aggregated = aggregated.sort_values(sort_columns).reset_index(drop=True)

    exact_columns = sort_columns + ["no_sites", "no_var_sites"]
    pd.testing.assert_frame_equal(direct[exact_columns], aggregated[exact_columns])

    float_columns = ["avg_watterson_theta", "raw_watterson_theta", "weighted_no_sites"]
    assert np.allclose(direct[float_columns].to_numpy(), aggregated[float_columns].to_numpy())
    assert np.allclose(
        aggregated["avg_watterson_theta"].to_numpy(),
        (aggregated["raw_watterson_theta"] / aggregated["no_sites"]).to_numpy(),
    )


@pytest.mark.regression
def test_pixy_tajima_d_aggregation_matches_direct_calculation(
    pixy_out_dir: Path,
    ag1000_pop_path: Path,
    ag1000_vcf_path: Path,
) -> None:
    """Chunk-aggregated Tajima's D should match direct window calculation."""
    shared_args: Dict[str, Any] = {
        "pixy_out_dir": pixy_out_dir,
        "stats": ["tajima_d"],
        "window_size": 10000,
        "vcf_path": ag1000_vcf_path,
        "populations_path": ag1000_pop_path,
        "chromosomes": "X",
        "cores": 1,
    }
    run_pixy_helper(
        **shared_args,
        output_prefix="tajima_direct",
        chunk_size=100000,
    )
    run_pixy_helper(
        **shared_args,
        output_prefix="tajima_aggregated",
        chunk_size=5000,
    )

    direct = pd.read_csv(pixy_out_dir / "tajima_direct_tajima_d.txt", sep="\t")
    aggregated = pd.read_csv(pixy_out_dir / "tajima_aggregated_tajima_d.txt", sep="\t")
    assert "tajima_d_s_counts" not in direct.columns
    assert "tajima_d_s_counts" not in aggregated.columns
    sort_columns = ["pop", "chromosome", "window_pos_1", "window_pos_2"]
    direct = direct.sort_values(sort_columns).reset_index(drop=True)
    aggregated = aggregated.sort_values(sort_columns).reset_index(drop=True)

    exact_columns = sort_columns + ["no_sites"]
    pd.testing.assert_frame_equal(direct[exact_columns], aggregated[exact_columns])

    float_columns = ["tajima_d", "raw_pi", "raw_watterson_theta", "tajima_d_stdev"]
    assert np.allclose(
        direct[float_columns].to_numpy(),
        aggregated[float_columns].to_numpy(),
        equal_nan=True,
    )

    valid_denominators = aggregated["tajima_d_stdev"] > 0
    assert np.allclose(
        aggregated.loc[valid_denominators, "tajima_d"].to_numpy(),
        (
            (
                aggregated.loc[valid_denominators, "raw_pi"]
                - aggregated.loc[valid_denominators, "raw_watterson_theta"]
            )
            / aggregated.loc[valid_denominators, "tajima_d_stdev"]
        ).to_numpy(),
    )


@pytest.mark.regression
def test_pixy_tajima_d_components_enable_posthoc_aggregation(
    pixy_out_dir: Path,
    ag1000_pop_path: Path,
    ag1000_vcf_path: Path,
) -> None:
    """The --tajima_components flag should emit enough data to aggregate Tajima's D."""
    shared_args: Dict[str, Any] = {
        "pixy_out_dir": pixy_out_dir,
        "stats": ["tajima_d"],
        "vcf_path": ag1000_vcf_path,
        "populations_path": ag1000_pop_path,
        "chromosomes": "X",
        "cores": 1,
        "chunk_size": 100000,
    }
    run_pixy_helper(
        **shared_args,
        window_size=10000,
        output_prefix="tajima_components_10kb",
        tajima_components=True,
    )
    run_pixy_helper(
        **shared_args,
        window_size=20000,
        output_prefix="tajima_direct_20kb",
    )

    components = pd.read_csv(pixy_out_dir / "tajima_components_10kb_tajima_d.txt", sep="\t")
    direct = pd.read_csv(pixy_out_dir / "tajima_direct_20kb_tajima_d.txt", sep="\t")
    assert "tajima_d_s_counts" in components.columns
    assert "tajima_d_s_counts" not in direct.columns

    components["posthoc_window_pos_1"] = ((components["window_pos_1"] - 1) // 20000) * 20000 + 1
    components["posthoc_window_pos_2"] = components["posthoc_window_pos_1"] + 19999
    group_columns = ["pop", "chromosome", "posthoc_window_pos_1", "posthoc_window_pos_2"]

    posthoc_rows = []
    for group_key, group in components.groupby(group_columns, sort=False):
        variant_counts: Dict[int, int] = {}
        for value in group["tajima_d_s_counts"]:
            for n, s in deserialize_tajima_d_variant_counts(value).items():
                variant_counts[n] = variant_counts.get(n, 0) + s

        raw_pi = group["raw_pi"].sum()
        raw_watterson_theta = group["raw_watterson_theta"].sum()
        tajima_d_stdev = calc_tajima_d_stdev(variant_counts)
        tajima_d = (raw_pi - raw_watterson_theta) / tajima_d_stdev if tajima_d_stdev > 0 else np.nan
        pop, chromosome, window_pos_1, window_pos_2 = group_key
        posthoc_rows.append({
            "pop": pop,
            "chromosome": chromosome,
            "window_pos_1": window_pos_1,
            "window_pos_2": window_pos_2,
            "tajima_d": tajima_d,
            "no_sites": group["no_sites"].sum(),
            "raw_pi": raw_pi,
            "raw_watterson_theta": raw_watterson_theta,
            "tajima_d_stdev": tajima_d_stdev,
        })

    posthoc = pd.DataFrame(posthoc_rows)
    sort_columns = ["pop", "chromosome", "window_pos_1", "window_pos_2"]
    posthoc = posthoc.sort_values(sort_columns).reset_index(drop=True)
    direct = direct.sort_values(sort_columns).reset_index(drop=True)

    exact_columns = sort_columns + ["no_sites"]
    pd.testing.assert_frame_equal(direct[exact_columns], posthoc[exact_columns])
    float_columns = ["tajima_d", "raw_pi", "raw_watterson_theta", "tajima_d_stdev"]
    assert np.allclose(
        direct[float_columns].to_numpy(),
        posthoc[float_columns].to_numpy(),
        equal_nan=True,
    )


################################################################################
# Tests for pixy.main(): baseline input
################################################################################


@pytest.mark.regression
def test_pixy_main_valid_inputs(
    pixy_out_dir: Path,
    expected_outputs: Path,
    ag1000_pop_path: Path,
    ag1000_vcf_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """
    Given specific input data, assert that outputs do not change.

    Uses `filecmp` library to compare 2 files without opening them and reading line-by-line.
    `filecmp.cmp` returns True if 2 files are equal.
    """
    caplog.set_level(logging.INFO)
    run_pixy_helper(
        pixy_out_dir=pixy_out_dir,
        stats=["pi", "fst", "dxy"],
        window_size=10000,
        vcf_path=ag1000_vcf_path,
        populations_path=ag1000_pop_path,
    )
    assert (
        "[pixy] Data set contains 2 populations, 2 chromosome(s), and 36 sample(s)"
        in caplog.messages
    )

    expected_out_files: List[Path] = [
        Path("pixy_dxy.txt"),
        Path("pixy_fst.txt"),
        Path("pixy_pi.txt"),
    ]
    for file in expected_out_files:
        generated_data_path: Path = pixy_out_dir / file
        exp_data_path: Path = expected_outputs / "baseline" / file
        assert generated_data_path.exists()

        assert_files_are_consistent(generated_data_path, exp_data_path)


@pytest.mark.regression
def test_pixy_multicore_matches_single_core(
    pixy_out_dir: Path,
    tmp_path: Path,
    ag1000_pop_path: Path,
    ag1000_vcf_path: Path,
) -> None:
    """Smoke test for ``--n_cores > 1``: 2-core output must equal 1-core output byte-for-byte.

    Exercises the parts of ``pixy.__main__.main`` that are only reachable when
    ``--n_cores >= 2`` — the multiprocessing manager, the worker pool, the
    cross-process queue, and the temp-file listener. None of those are touched
    by the rest of the suite (which all runs at the default ``--n_cores 1``).

    The test pins ``--chunk_size`` to a value small enough that the input VCF
    produces several chunks; otherwise a 2-core pool only ever gets one chunk
    of work and the listener is never exercised. Output is compared directly
    to a 1-core run on the same inputs (not to a baseline file) so this stays
    decoupled from the output-format regression baselines.
    """
    sc_out = tmp_path / "single_core"
    mc_out = tmp_path / "multi_core"
    sc_out.mkdir()
    mc_out.mkdir()

    common = dict(
        stats=["pi", "fst", "dxy"],
        window_size=10000,
        vcf_path=ag1000_vcf_path,
        populations_path=ag1000_pop_path,
        # chunk_size must be >= window_size; this value yields ~5 chunks per chromosome
        # for the ag1000 test fixture, enough for a 2-core pool to see real parallelism.
        chunk_size=20000,
    )
    run_pixy_helper(pixy_out_dir=sc_out, cores=1, **common)
    run_pixy_helper(pixy_out_dir=mc_out, cores=2, **common)

    for name in ("pixy_pi.txt", "pixy_dxy.txt", "pixy_fst.txt"):
        sc = sc_out / name
        mc = mc_out / name
        assert sc.exists(), f"single-core run produced no {name}"
        assert mc.exists(), f"multi-core run produced no {name}"
        # `assert_files_are_consistent` sorts then line-compares, which is what we want:
        # multicore output rows can arrive in non-deterministic order across runs.
        assert_files_are_consistent(mc, sc)


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
        # shutil.copy(generated_data_path, exp_data_path)
        assert_files_are_consistent(generated_data_path, exp_data_path)


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
        # shutil.copy(generated_data_path, exp_data_path)
        assert_files_are_consistent(generated_data_path, exp_data_path)


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
        # assert generated_data_path.exists()
        assert_files_are_consistent(generated_data_path, exp_data_path)


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
        # shutil.copy(generated_data_path, exp_data_path)
        assert_files_are_consistent(generated_data_path, exp_data_path)


@pytest.mark.regression
@pytest.mark.parametrize(
    "fst_type, output_prefix",
    [
        ("wc", "wc_fst_components"),
        ("hudson", "hudson_fst_components"),
    ],
)
def test_pixy_fst_components_flag(
    pixy_out_dir: Path,
    expected_outputs: Path,
    ag1000_pop_path: Path,
    ag1000_vcf_path: Path,
    fst_type: str,
    output_prefix: str,
) -> None:
    """The --fst_components flag should append estimator components to FST output."""
    run_pixy_helper(
        pixy_out_dir=pixy_out_dir,
        stats=["fst"],
        window_size=10000,
        vcf_path=ag1000_vcf_path,
        populations_path=ag1000_pop_path,
        output_prefix=output_prefix,
        debug=True,
        cores=1,
        fst_type=fst_type,
        fst_components=True,
    )

    generated_data_path: Path = pixy_out_dir / f"{output_prefix}_fst.txt"
    exp_data_path: Path = expected_outputs / "fst_components" / f"{output_prefix}_fst.txt"
    assert generated_data_path.exists()
    assert_files_are_consistent(generated_data_path, exp_data_path)


@pytest.mark.regression
def test_pixy_fst_numeric_population_labels(
    pixy_out_dir: Path,
    expected_outputs: Path,
    simulated_pop_dxy_path: Path,
    missing50_vcf_path: Path,
) -> None:
    """Numeric population labels should still map to the correct VCF sample indices."""
    run_pixy_helper(
        pixy_out_dir=pixy_out_dir,
        stats=["fst"],
        window_size=10000,
        vcf_path=missing50_vcf_path,
        populations_path=simulated_pop_dxy_path,
        output_prefix="missing_genotypes",
        bypass_invariant_check=True,
        debug=True,
        cores=1,
    )

    generated_data_path: Path = pixy_out_dir / "missing_genotypes_fst.txt"
    exp_data_path: Path = expected_outputs / "missing_genotypes" / "missing_genotypes_fst.txt"
    assert generated_data_path.exists()
    assert_files_are_consistent(generated_data_path, exp_data_path)


################################################################################
# Tests for pixy.main(): multiallelic SNPs
################################################################################


@pytest.mark.regression
def test_pixy_multiallelic(
    pixy_out_dir: Path,
    expected_outputs: Path,
    multiallelic_pop_path: Path,
    multiallelic_vcf_path: Path,
) -> None:
    """Test that pixy produces deterministic stats with multiallelic SNPs included."""
    run_pixy_helper(
        pixy_out_dir=pixy_out_dir,
        stats=["pi", "dxy", "fst"],
        window_size=10000,
        vcf_path=multiallelic_vcf_path,
        populations_path=multiallelic_pop_path,
        include_multiallelic_snps=True,
        bypass_invariant_check=True,
        output_prefix="multiallelic",
    )

    expected_out_files: List[Path] = [
        Path("multiallelic_pi.txt"),
        Path("multiallelic_dxy.txt"),
        Path("multiallelic_fst.txt"),
    ]
    for file in expected_out_files:
        generated_data_path: Path = pixy_out_dir / file
        exp_data_path: Path = expected_outputs / "multiallelic" / file
        assert generated_data_path.exists()
        # shutil.copy(generated_data_path, exp_data_path)
        assert_files_are_consistent(generated_data_path, exp_data_path)


# Tests for pixy.main(): variable-ploidy VCFs
################################################################################


def _read_chromosomes_from_pixy_output(path: Path) -> List[str]:
    """Return the unique chromosome values from a pixy output file (tab-separated)."""
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        chrom_idx = header.index("chromosome")
        chroms = {line.rstrip("\n").split("\t")[chrom_idx] for line in f if line.strip()}
    return sorted(chroms)


@pytest.mark.regression
def test_mixed_ploidy_pi_dxy(
    pixy_out_dir: Path,
    mixed_ploidy_vcf_path: Path,
    mixed_ploidy_pop_path: Path,
) -> None:
    """
    pixy runs end-to-end on a VCF with variable ploidy across contigs.

    The simulated VCF has 500 diploid sites on ``chr1`` and 500 haploid sites on ``chrX``.
    With per-contig ploidy inference, pi and dxy should be computed for both contigs without
    needing the user to split the VCF.
    """
    run_pixy_helper(
        pixy_out_dir=pixy_out_dir,
        stats=["pi", "dxy"],
        window_size=100,
        vcf_path=mixed_ploidy_vcf_path,
        populations_path=mixed_ploidy_pop_path,
        output_prefix="mixed",
    )

    pi_path = pixy_out_dir / "mixed_pi.txt"
    dxy_path = pixy_out_dir / "mixed_dxy.txt"
    assert pi_path.exists()
    assert dxy_path.exists()

    # both contigs (one diploid, one haploid) must appear in the output
    assert _read_chromosomes_from_pixy_output(pi_path) == ["chr1", "chrX"]
    assert _read_chromosomes_from_pixy_output(dxy_path) == ["chr1", "chrX"]


@pytest.mark.regression
def test_mixed_ploidy_hudson_fst_runs_on_both_contigs(
    pixy_out_dir: Path,
    mixed_ploidy_vcf_path: Path,
    mixed_ploidy_pop_path: Path,
) -> None:
    """Hudson's FST is defined for any ploidy and should run on both contigs."""
    run_pixy_helper(
        pixy_out_dir=pixy_out_dir,
        stats=["fst"],
        window_size=100,
        vcf_path=mixed_ploidy_vcf_path,
        populations_path=mixed_ploidy_pop_path,
        output_prefix="mixed",
        fst_type="hudson",
    )

    fst_path = pixy_out_dir / "mixed_fst.txt"
    assert fst_path.exists()
    chroms = _read_chromosomes_from_pixy_output(fst_path)
    assert "chr1" in chroms
    assert "chrX" in chroms


@pytest.mark.regression
def test_mixed_ploidy_wc_fst_skips_haploid_contig(
    pixy_out_dir: Path,
    mixed_ploidy_vcf_path: Path,
    mixed_ploidy_pop_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """
    Weir-Cockerham FST is diploid-only.

    When WC is requested on a mixed-ploidy VCF, the haploid contig should be skipped (with a
    warning) and pi should still be computed for both contigs.
    """
    with caplog.at_level(logging.WARNING):
        run_pixy_helper(
            pixy_out_dir=pixy_out_dir,
            stats=["pi", "fst"],
            window_size=100,
            vcf_path=mixed_ploidy_vcf_path,
            populations_path=mixed_ploidy_pop_path,
            output_prefix="mixed",
            fst_type="wc",
        )

    # pi runs on both contigs regardless of FST estimator
    pi_path = pixy_out_dir / "mixed_pi.txt"
    assert pi_path.exists()
    assert _read_chromosomes_from_pixy_output(pi_path) == ["chr1", "chrX"]

    # FST output should only contain the diploid contig (chr1); chrX was skipped.
    fst_path = pixy_out_dir / "mixed_fst.txt"
    if fst_path.exists():
        fst_chroms = _read_chromosomes_from_pixy_output(fst_path)
        assert "chrX" not in fst_chroms, (
            "WC FST output unexpectedly contains the haploid contig chrX"
        )
        assert fst_chroms == ["chr1"]

    # a warning about WC FST being skipped on non-diploid contigs should have been emitted
    warning_msgs = " ".join(record.getMessage() for record in caplog.records)
    assert "Weir-Cockerham FST is not supported for non-diploid contigs" in warning_msgs
