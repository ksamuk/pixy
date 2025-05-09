from enum import Enum
from enum import unique


@unique
class PixyStat(Enum):
    """
    The genetic variance statistics that `pixy` can calculate.

    Attributes:
        PI: a measure of genetic diversity _within_ populations
        DXY: a measure of genetic diversity _between_ populations
        FST: the "fixation index"; represents a subpopulation-specific genetic variance
            relative to the total genetic variance
        TAJIMA_D: a measure of non-neutral evolution in a population
        WATTERSON_THETA: an estimator of genetic diversity
    """

    PI = "pi"
    DXY = "dxy"
    FST = "fst"
    TAJIMA_D = "tajima_d"
    WATTERSON_THETA = "watterson_theta"

    def __str__(self) -> str:
        return self.value


@unique
class FSTEstimator(Enum):
    """
    The mutually exclusive FST estimators that `pixy` offers.

    Attributes:
        WC: model from Weir and Cockerham 1984
        HUDSON: model from Hudson 1992 and Bhatia et al. 2013
    """

    WC = "wc"
    HUDSON = "hudson"

    def __str__(self) -> str:
        return self.value
