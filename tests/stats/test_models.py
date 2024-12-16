from typing import Union

import pytest

from pixy.args_validation import PixyStat
from pixy.models import NA
from pixy.models import PixyTempResult


@pytest.mark.parametrize(
    "pixy_stat, pop1, pop2, chr, pos_1, pos_2, calculated_stat, shared, diffs, comps, missing,"
    "expected_str",
    [
        (
            PixyStat.DXY,
            "pop1",
            "pop2",
            "chr1",
            100,
            200,
            0.5,
            50,
            10,
            20,
            5,
            "dxy\tpop1\tpop2\tchr1\t100\t200\t0.5\t50\t10\t20\t5",
        ),  # no `None` -> no NA expected
        (
            PixyStat.FST.value,
            "pop1",
            "pop2",
            "chr1",
            100,
            200,
            0.5,
            50,
            "NA",
            20,
            5,
            "fst\tpop1\tpop2\tchr1\t100\t200\t0.5\t50\tNA\t20\t5",
        ),  # `total_differences` is `None`
        (
            PixyStat.FST,
            "pop1",
            "pop2",
            "chr1",
            100,
            200,
            0.5,
            50,
            10,
            "NA",
            5,
            "fst\tpop1\tpop2\tchr1\t100\t200\t0.5\t50\t10\tNA\t5",
        ),  # `total_comparisons` is `None`
        (
            PixyStat.FST,
            "pop1",
            "pop2",
            "chr1",
            100,
            200,
            0.5,
            50,
            10,
            20,
            "NA",
            "fst\tpop1\tpop2\tchr1\t100\t200\t0.5\t50\t10\t20\tNA",
        ),  # `total_missing` is `None`
        (
            PixyStat.PI,
            "pop1",
            "NA",
            "chr1",
            100,
            200,
            0.5,
            50,
            10,
            20,
            5,
            "pi\tpop1\tNA\tchr1\t100\t200\t0.5\t50\t10\t20\t5",
        ),  # `population_2` is `None`
    ],
)
def test_pixy_result_str(
    pixy_stat: PixyStat,
    pop1: str,
    pop2: str,
    chr: str,
    pos_1: int,
    pos_2: int,
    calculated_stat: float,
    shared: int,
    diffs: Union[int, NA],
    comps: Union[int, NA],
    missing: Union[int, NA],
    expected_str: str,
) -> None:
    """Assert that PixyResult.__str__() produces well-formed output."""
    generated_result: PixyTempResult = PixyTempResult(
        pixy_stat=pixy_stat,
        population_1=pop1,
        population_2=pop2,
        chromosome=chr,
        window_pos_1=pos_1,
        window_pos_2=pos_2,
        calculated_stat=calculated_stat,
        shared_sites_with_alleles=shared,
        total_differences=diffs,
        total_comparisons=comps,
        total_missing=missing,
    )
    assert str(generated_result) == expected_str
