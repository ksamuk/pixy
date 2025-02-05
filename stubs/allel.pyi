"""
Minimal type hints for scikit-allel.

This stub module includes type hints for the `scikit-allel` classes, methods, and functions used by
`pixy`.

The structure of this module does not precisely mirror the structure of `scikit-allel`. For
simplicity, all type hints are provided in a single file. This does mirror the API of
`scikit-allel`, as all the public classes and functions typed within this stub module are available
as top-level imports from the `scikit-allel` module.
"""

from collections.abc import Sequence
from typing import Any
from typing import Dict
from typing import List
from typing import Literal
from typing import Optional
from typing import TextIO
from typing import Tuple
from typing import Union

import numpy as np
from allel.io.vcf_read import DEFAULT_ALT_NUMBER  # type: ignore[import-untyped]
from allel.io.vcf_read import DEFAULT_BUFFER_SIZE
from allel.io.vcf_read import DEFAULT_CHUNK_LENGTH
from numpy.typing import NDArray

####################################################################################################
# allel.abc
# https://github.com/cggh/scikit-allel/blob/master/allel/abc.py
####################################################################################################

# NB: This class does not technically inherit from `NDArray`, but it has an `NDArray` attribute
# (`values`) and a `__getattr__` definition that permits transparent access to member methods and
# attributes of this array, so it functionally behaves as a subclass. Including the `NDArray`
# inheritance in this type stub simplified some of the child stubs and corresponding type hints.
class ArrayWrapper(NDArray):
    def __init__(self, data: Union[NDArray, "ArrayWrapper"]) -> None: ...

####################################################################################################
# allel.model.ndarray
# https://github.com/cggh/scikit-allel/blob/master/allel/model/ndarray.py
####################################################################################################
class NumpyArrayWrapper(ArrayWrapper):
    def __init__(self, data: NDArray, copy: bool = False, **kwargs: Any) -> None: ...
    def __len__(self) -> int: ...

class Genotypes(NumpyArrayWrapper):
    def __init__(self, data: NDArray, copy: bool = False, **kwargs: Any) -> None: ...
    @property
    def ploidy(self) -> int: ...

class AlleleCountsArray(NumpyArrayWrapper):
    def __init__(self, data: NDArray, copy: bool = False, **kwargs: Any) -> None: ...
    def __getitem__(self, item: Any) -> Any: ...
    def to_frequencies(self, fill: float = np.nan) -> NDArray: ...

class GenotypeArray(Genotypes):
    def __init__(
        self,
        data: NDArray | List[List[List[int]]],
        copy: bool = False,
        **kwargs: Any,
    ) -> None: ...
    def __getitem__(self, item: Any) -> Any: ...
    @property
    def n_variants(self) -> int: ...
    @property
    def n_samples(self) -> int: ...
    def count_alleles(
        self,
        max_allele: Optional[int] = None,
        subpop: Optional[Sequence[int]] = None,
    ) -> AlleleCountsArray: ...
    def concatenate(
        self,
        others: Union["GenotypeArray", List["GenotypeArray"], Tuple["GenotypeArray"]],
        axis: int = 0,
    ) -> "GenotypeArray": ...

class GenotypeVector(Genotypes):
    def __init__(self, data: NDArray, copy: bool = False, **kwargs: Any) -> None: ...
    def __getitem__(self, item: Any) -> Any: ...
    def take(  # type: ignore[override]
        self,
        indices: NDArray,
        axis: int = 0,
        out: Optional[NDArray] = None,
        mode: Union[Literal["raise"], Literal["wrap"], Literal["clip"]] = "raise",
    ) -> GenotypeArray: ...

class SortedIndex(NumpyArrayWrapper):
    def __init__(self, data: NDArray, copy: bool = False, **kwargs: Any) -> None: ...
    def __getitem__(self, item: Any) -> Any: ...
    def locate_keys(self, keys: Sequence[int], strict: bool = True) -> NDArray: ...
    def locate_range(self, start: Optional[int] = None, stop: Optional[int] = None) -> slice: ...

####################################################################################################
# allel.model.dask
# https://github.com/cggh/scikit-allel/blob/master/allel/model/dask.py
####################################################################################################
class GenotypeDaskArray(ArrayWrapper):
    def __init__(
        self,
        data: NDArray,
        chunks: Optional[Tuple[int, ...]] = None,
        name: Optional[Union[str, bool]] = None,
        lock: bool = False,
    ) -> None: ...

####################################################################################################
# allel.stats.fst
# https://github.com/cggh/scikit-allel/blob/master/allel/stats/fst.py
####################################################################################################
def weir_cockerham_fst(
    g: GenotypeArray,
    subpops: Sequence[Sequence[int]],
    max_allele: Optional[int] = None,
    blen: Optional[int] = None,
) -> Tuple[NDArray, NDArray, NDArray]: ...
def hudson_fst(
    ac1: AlleleCountsArray,
    ac2: AlleleCountsArray,
    fill: float = np.nan,
) -> Tuple[NDArray, NDArray]: ...
def windowed_hudson_fst(
    pos: NDArray,
    ac1: NDArray,
    ac2: NDArray,
    size: Optional[int] = None,
    start: Optional[int] = None,
    stop: Optional[int] = None,
    step: Optional[int] = None,
    windows: Optional[NDArray] = None,
    fill: float = np.nan,
) -> Tuple[NDArray, NDArray, NDArray]: ...
def windowed_weir_cockerham_fst(
    pos: NDArray,
    g: NDArray,
    subpops: Sequence[Sequence[int]],
    size: Optional[int] = None,
    start: Optional[int] = None,
    stop: Optional[int] = None,
    step: Optional[int] = None,
    windows: Optional[NDArray] = None,
    fill: float = np.nan,
    max_allele: Optional[int] = None,
) -> Tuple[NDArray, NDArray, NDArray]: ...

####################################################################################################
# allel.io.vcf_read
# https://github.com/cggh/scikit-allel/blob/master/allel/io/vcf_read.py
####################################################################################################
class VCFHeaders:
    headers: List[str]
    filters: dict[str, Any]
    infos: dict[str, Any]
    formats: dict[str, Any]
    # I think this is technically optional in the reference implementation, but it's always
    # populated in the pixy VCFs
    samples: List[str]

def read_vcf(
    input: Union[str, TextIO],
    fields: Optional[List[str]] = None,
    exclude_fields: Optional[List[str]] = None,
    rename_fields: Optional[Dict[str, str]] = None,
    types: Optional[Dict[str, Any]] = None,
    numbers: Optional[Dict[str, Any]] = None,
    alt_number: int = DEFAULT_ALT_NUMBER,
    fills: Optional[Dict[str, Any]] = None,
    region: Optional[str] = None,
    tabix: Optional[str] = "tabix",
    samples: Optional[List[str]] = None,
    transformers: Optional[List[Any]] = None,
    buffer_size: int = DEFAULT_BUFFER_SIZE,
    chunk_length: int = DEFAULT_CHUNK_LENGTH,
    log: Optional[TextIO] = None,
) -> Optional[Dict[str, NDArray]]: ...
def read_vcf_headers(input: Any) -> VCFHeaders: ...
