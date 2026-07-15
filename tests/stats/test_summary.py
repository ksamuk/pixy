import argparse

import numpy as np
import pytest
from allel import GenotypeArray
from allel import SortedIndex

from pixy.enums import FSTEstimator
from pixy.stats.summary import precompute_filtered_variant_array

# site 1 is invariant, site 2 is biallelic, and site 3 is triallelic
GT_ARRAY = GenotypeArray([
    [[0, 0], [0, 0], [0, 0], [0, 0]],
    [[0, 0], [0, 1], [0, 0], [1, 1]],
    [[0, 0], [0, 1], [2, 2], [1, 2]],
])

POS_ARRAY = SortedIndex(np.array([100, 200, 300]))


@pytest.mark.parametrize(
    "fst_type, fst_biallelic, expected_positions",
    [
        # Weir-Cockerham assumes biallelic data, so multiallelic sites are dropped
        (FSTEstimator.WC, False, [200]),
        # pixy's Hudson estimator handles any number of alleles, so only the invariant site is
        # dropped
        (FSTEstimator.HUDSON, False, [200, 300]),
        # --fst_biallelic opts Hudson back in to biallelic-only filtering
        (FSTEstimator.HUDSON, True, [200]),
        # ...and is a no-op for Weir-Cockerham, which is always biallelic
        (FSTEstimator.WC, True, [200]),
    ],
)
def test_precompute_filtered_variant_array_retains_multiallelic_sites_for_hudson_only(
    fst_type: FSTEstimator,
    fst_biallelic: bool,
    expected_positions: list,
) -> None:
    """Multiallelic sites should be retained only for Hudson FST without --fst_biallelic."""
    args = argparse.Namespace(
        populations="populations.txt",
        fst_type=fst_type.value,
        fst_biallelic=fst_biallelic,
    )

    _, gt_array_fst, pos_array_fst = precompute_filtered_variant_array(
        args=args,
        gt_array=GT_ARRAY,
        pos_array=POS_ARRAY,
        callset_is_none=False,
        window_size=1000,
        popindices={"pop1": np.array([0, 1]), "pop2": np.array([2, 3])},
        chromosome="chr1",
    )

    assert gt_array_fst is not None
    assert pos_array_fst is not None
    assert list(pos_array_fst) == expected_positions
    assert gt_array_fst.n_variants == len(expected_positions)
