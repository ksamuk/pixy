import allel
import numpy as np
from allel import GenotypeArray
from allel import HaplotypeArray
from allel import SortedIndex

from pixy.core import mask_non_target_sites


def test_mask_non_target_sites_diploid() -> None:
    """Non-target sites in a diploid GenotypeArray are replaced with missing rows."""
    gt_array = GenotypeArray([
        [[0, 0], [0, 1], [1, 1]],
        [[0, 1], [0, 0], [1, 1]],
        [[1, 1], [0, 1], [0, 0]],
    ])
    pos_array = SortedIndex(np.array([10, 20, 30]))
    sites_list_chunk = [20]

    masked = mask_non_target_sites(gt_array, pos_array, sites_list_chunk)

    assert np.all(masked[0] == -1)
    assert np.all(masked[2] == -1)
    assert np.array_equal(masked[1], np.array([[0, 1], [0, 0], [1, 1]]))


def test_mask_non_target_sites_haploid() -> None:
    """
    Non-target sites in a haploid HaplotypeArray are replaced with missing rows.

    Regression for https://github.com/ksamuk/pixy/issues/191: the previous implementation accessed
    ``gt_array.ploidy``/``n_samples``, which do not exist on a ``HaplotypeArray`` and raised an
    ``AttributeError`` whenever ``--sites_file`` was used on haploid input.
    """
    hap_array = HaplotypeArray([
        [0, 1, 0, 1],
        [1, 1, 0, 0],
        [0, 0, 1, 1],
    ])
    pos_array = SortedIndex(np.array([10, 20, 30]))
    sites_list_chunk = [20]

    masked = mask_non_target_sites(hap_array, pos_array, sites_list_chunk)

    assert isinstance(masked, allel.HaplotypeArray)
    assert np.all(np.asarray(masked[0]) == -1)
    assert np.all(np.asarray(masked[2]) == -1)
    assert np.array_equal(np.asarray(masked[1]), np.array([1, 1, 0, 0]))
