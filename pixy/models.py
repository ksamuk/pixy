from dataclasses import dataclass
from dataclasses import fields
from typing import Literal
from typing import Union

from typing_extensions import TypeAlias

from pixy.args_validation import PixyStat

NA: TypeAlias = Literal["NA"]


@dataclass(frozen=True)
class PixyTempResult:
    """
    Stores temporary `pixy` results.

    Currently, `pixy` computes results for all requested statistics ("pi", "dxy", or "fst")
    and writes them to the same temporary file before splitting this temp file into the final
    per-statistic output files.

    This dataclass represents a row in this temporary file.

    Attributes:
        pixy_stat: the genetic variance statistic that `pixy` calculated
            (one of "pi", "dxy", or "fst")
        population_1: one of two populations being compared (found in the `populations_path` file)
        population_2: the other of two populations being compared
        chromosome: the chromosome of interest
        window_pos_1: the start position of the window in which `pixy` is calculating variance
        window_pos_2: the end position of the window in which `pixy` is calculating variance
        calculated_stat: the result of the calculation of the `pixy_stat`
        shared_sites_with_alleles: the sites at which both populations have alleles
        total_differences: the total number of differences between populations summed across all
            sites
        total_comparisons: the total number of pairwise comparisons between sites
        total_missing: equal to (total count of possible pairwise comparisons at all sites
            - total_comparisons)

    `population_1` and `population_2` will both be populated for "dxy" and "fst" results.
    `population_2` will be "NA" when reporting "pi" results.

    """

    pixy_stat: PixyStat
    population_1: str
    population_2: Union[str, NA]
    chromosome: str
    window_pos_1: int
    window_pos_2: int
    calculated_stat: Union[float, NA]
    shared_sites_with_alleles: int
    total_differences: Union[int, float, NA]
    total_comparisons: Union[int, float, NA]
    total_missing: Union[int, float, NA]

    def __str__(self) -> str:
        """Returns a tab-delimited string representation for writing out to the temp file."""
        return "\t".join(str(getattr(self, field.name)) for field in fields(self))


@dataclass
class PiResult:
    """A result from calculating pi."""

    avg_pi: Union[float, NA]
    total_diffs: Union[int, NA]
    total_comps: Union[int, NA]
    total_missing: Union[int, NA]

    @classmethod
    def empty(cls) -> "PiResult":
        """An empty PiResult, with all fields set to NA."""
        return cls(avg_pi="NA", total_diffs="NA", total_comps="NA", total_missing="NA")


@dataclass
class DxyResult:
    """A result from calculating dxy."""

    avg_dxy: Union[float, NA]
    total_diffs: Union[int, NA]
    total_comps: Union[int, NA]
    total_missing: Union[int, NA]

    @classmethod
    def empty(cls) -> "DxyResult":
        """An empty DxyResult, with all fields set to NA."""
        return cls(avg_dxy="NA", total_diffs="NA", total_comps="NA", total_missing="NA")


@dataclass
class FstResult:
    """
    A result from calculating fst.

    If `fst_type` is `wc`, the following will be returned:
        fst: "NA" if no valid data
        a: variance between populations
        b: variance between individuals within populations
        c: variance between gametes within individuals
        n_sites: the number of sites over which variance was measured

    If `fst_type` is `hudson`, the following will be returned:
        fst: "NA" if no valid data.
        num: divergence between the two populations minus average of diversity within each
            population
        den: divergence between the two populations
        c: a placeholder of 0 to maintain the shape of the return Tuple
        n_sites: the number of sites over which variance was measured

    For simplicity and historical consistency, when `fst_type` is `hudson`, the `a` and `b` fields
    contain the `num` and `den` values.
    """

    fst: Union[float, NA]
    a: Union[float, NA]
    b: Union[float, NA]
    c: Union[float, int, NA]
    n_sites: int

    @classmethod
    def empty(cls) -> "FstResult":
        """An empty FstResult, with all fields set to NA."""
        return cls(fst="NA", a="NA", b="NA", c="NA", n_sites=0)
