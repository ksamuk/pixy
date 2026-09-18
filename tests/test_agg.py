import pytest

from pixy.agg import _final_stat


@pytest.mark.parametrize(
    "fst_type, diffs, comps, missing, expected",
    [
        # a negative numerator is a valid (negative) FST
        ("hudson", -0.5, 0.5, 0, -1.0),
        ("wc", -0.1, 0.3, 0.2, -0.25),
        # a non-positive denominator is not: `calc_fst` reports NA, so aggregation must too
        ("hudson", 0.1, 0.0, 0, "NA"),
        ("hudson", 0.1, -0.2, 0, "NA"),
        ("wc", 0.1, -0.3, 0.1, "NA"),
        ("wc", 0.0, 0.0, 0.0, "NA"),
    ],
)
def test_final_stat_fst_guards_match_calc_fst(
    fst_type: str, diffs: float, comps: float, missing: float, expected: object
) -> None:
    """Aggregated FST is defined under exactly the same condition as windowed FST."""
    result = _final_stat("fst", fst_type, 10, diffs, comps, missing)

    assert result == (expected if expected == "NA" else pytest.approx(expected))
