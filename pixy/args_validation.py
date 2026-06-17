# `from __future__ import annotations` lazies all annotations into strings so referencing
# `allel.*` in type hints doesn't force that module to be imported at module load. The actual
# `import allel` is no longer needed here (we now use a tiny gzip-based VCF header reader),
# and pandas has been replaced wholesale with stdlib-based tables — see the dataclasses
# below. This keeps `pixy --help`, `--version`, and arg-parse error paths from paying the
# ~500 ms / ~130 MB cost of pandas import.
from __future__ import annotations

import argparse
import csv
import gzip
import logging
import os
import shutil
import subprocess
import uuid
from dataclasses import dataclass
from dataclasses import field
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING
from typing import Dict
from typing import List
from typing import Tuple
from typing import Union

import numpy as np
from numpy.typing import NDArray

from pixy.enums import FSTEstimator
from pixy.enums import PixyStat

if TYPE_CHECKING:
    from pixy.wisp import WispMask

# ---------------------------------------------------------------------------
# Pandas-free table types that replace the previous pandas.DataFrame fields on
# PixyArgs (populations_df / bed_df / sites_df). Each holds the same per-row
# information as the corresponding DataFrame did, with the access patterns used
# by downstream code exposed as small methods.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PopulationsTable:
    """
    Sample-id / population mapping derived from the --populations file.

    Replaces the previous `populations_df` DataFrame, whose columns were
    `ID`, `Population`, and (added after VCF lookup) `callset_index`.
    """

    ids: Tuple[str, ...]
    populations: Tuple[str, ...]
    # Per-sample index into the VCF's sample list. Populated after `_read_vcf_samples` runs.
    callset_index: Tuple[int, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if len(self.ids) != len(self.populations):
            raise ValueError("ids and populations must be the same length")
        if self.callset_index and len(self.callset_index) != len(self.ids):
            raise ValueError("callset_index, when set, must be the same length as ids")

    @cached_property
    def unique_populations(self) -> Tuple[str, ...]:
        """Distinct population names in first-appearance order."""
        return tuple(dict.fromkeys(self.populations))

    @property
    def num_populations(self) -> int:
        """Number of distinct populations in the table."""
        return len(self.unique_populations)

    def indices_for(self, population: str) -> NDArray[np.intp]:
        """Callset indices of all samples in `population`. Order is preserved."""
        return np.array(
            [
                ci
                for pop, ci in zip(self.populations, self.callset_index, strict=True)
                if pop == population
            ],
            dtype=np.intp,
        )

    def with_callset_index(self, callset_index: Tuple[int, ...]) -> PopulationsTable:
        """Return a new table with the supplied per-sample callset indices."""
        return PopulationsTable(
            ids=self.ids, populations=self.populations, callset_index=tuple(callset_index)
        )


@dataclass(frozen=True)
class BedTable:
    """
    List of windowing intervals derived from a user-supplied BED file.

    Stored coordinates are 1-based inclusive (pixy's internal convention), already
    converted from the BED 0-based half-open input by ``validate_bed_path``.
    """

    chroms: Tuple[str, ...]
    chrom_starts: Tuple[int, ...]
    chrom_ends: Tuple[int, ...]

    def __post_init__(self) -> None:
        n = len(self.chroms)
        if len(self.chrom_starts) != n or len(self.chrom_ends) != n:
            raise ValueError("chroms, chrom_starts, chrom_ends must be the same length")

    def intervals_for(self, chromosome: str) -> List[List[int]]:
        """
        Return `[[start, end], ...]` for the given chromosome.

        Parallel to the previous zip-of-dataframe-columns pattern used in `__main__.py`.
        """
        return [
            [s, e]
            for c, s, e in zip(self.chroms, self.chrom_starts, self.chrom_ends, strict=True)
            if c == chromosome
        ]

    def unique_chroms(self) -> List[str]:
        """Distinct chromosomes in input order."""
        return list(dict.fromkeys(self.chroms))


@dataclass(frozen=True)
class SitesTable:
    """
    Per-chromosome sorted position lists from the --sites_file.

    Replaces the previous `sites_df` DataFrame whose columns were `CHROM` and `POS`.
    The previous code repeatedly filtered the DataFrame by chromosome and sorted the
    positions; we pre-compute and cache that here.
    """

    by_chrom: Dict[str, Tuple[int, ...]]

    @cached_property
    def chromosomes(self) -> Tuple[str, ...]:
        """Tuple of chromosome names that have at least one site."""
        return tuple(self.by_chrom.keys())

    def positions_for(self, chromosome: str) -> List[int]:
        """Sorted list of positions for `chromosome`, or `[]` if absent."""
        return list(self.by_chrom.get(chromosome, ()))


@dataclass(frozen=True)
class PixyArgs:
    """
    Holds settings for running `pixy`.

    `pixy` has the following mutually exclusive, valid groups of inputs:
        * a `vcf_path`, a `populations_path`, and a `window_size`
        * a `vcf_path`, a `populations_path`, and a `bed_path`
        * a `vcf_path`, a `populations_path`, a `window_size`, `interval_start`, `interval_end`, and
            a single value for `chromosomes`

    Note that if a BED path is specified, the user is prohibited from specifying a window size or
    interval. If no BED path is specified, the user is required to specify a window size, and the
    interval is optional. When specifying an interval, both the start and end are required.

    Attributes:
        stats: list of statistics for `pixy` to calculate from the input VCF
        vcf_path: path to the input VCF (bgzipped and indexed)
        populations: a PopulationsTable derived from the user-specified populations file
        num_cores: number of CPUs to utilize for parallel processing (default = 1)
        include_multiallelic_snps: If True, include multiallelic sites in the analysis
        bypass_invariant_check: whether to allow computation of stats without invariant sites
            (this option is never recommended and defaults to False)
        bed: a BedTable derived from a user-specified BED file, or `None` if no BED file was
            given.
        output_dir: an optional path to which outputs will be written; default is current directory
        output_prefix: an optional prefix with which to prepend to `pixy` output (default is `pixy`)
        chromosomes: an optional comma-separated list of chroms over which stats will be calculated
            (defaults to all chromosomes in a given VCF)
        window_size: an optional length of base pairs over which to calculate stats
        interval_start: an optional 1-based position demarcating the start of an interval over which
            to calculate stats (only valid when calculating over a single chromosome)
        interval_end: an optional 1-based position demarcating the end of an interval over which
            to calculate stats (only valid when calculating over a single chromosome)
        sites: a SitesTable derived from a user-specified sites file, or `None` if no sites
            file was given.
        chunk_size: the approximate number of sites to read from VCF at any given time
        fst_type: the FST estimator to use, one of either 'WC' (Weir and Cockerham 1984) or
            'HUDSON' (Hudson 1992, Bhatia et al. 2013). Defaults to 'WC'.
        fst_components: whether to include FST estimator components in the final FST output table
        tajima_components: whether to include Tajima's D aggregation components in the final
            Tajima's D output table
        temp_file: a Path to which to write intermediate `pixy` results, assigned based on the value
            of `output_dir`
        ploidy_map: a mapping from contig name to inferred ploidy. Built from the first record of
            each analyzed contig in the input VCF, allowing per-contig variable ploidy (e.g.
            diploid autosomes alongside haploid sex chromosomes or organellar contigs).
        wisp_mask: an optional ``WispMask`` parsed from --wisp_bed. When provided, --vcf
            is treated as a variants-only callset and the per-window callable-site denominator
            for pi, dxy, Tajima's D, and Watterson's theta is sourced from the wisp mask
            rather than from invariant sites in the VCF. FST is unaffected.

    Raises:
        ValueError: if an interval is specified without both start and end positions
        ValueError: if a BED file is specified, and a window size or interval are also specified
        ValueError: if a BED file is not specified, and a window size is not specified
    """

    stats: List[PixyStat]
    vcf_path: Path
    populations: PopulationsTable
    output_dir: Path
    temp_file: Path
    chromosomes: List[str]
    bypass_invariant_check: bool
    include_multiallelic_snps: bool
    num_cores: int = 1
    fst_type: FSTEstimator = FSTEstimator.WC
    fst_components: bool = False
    tajima_components: bool = False
    output_prefix: str = "pixy"
    chunk_size: int = 100000
    bed: Union[BedTable, None] = None
    window_size: Union[int, None] = None
    interval_start: Union[int, None] = None
    interval_end: Union[int, None] = None
    sites: Union[SitesTable, None] = None
    ploidy_map: Union[Dict[str, int], None] = None
    # When set, indicates that --vcf is a variants-only VCF and the per-site callable
    # denominator for pi/dxy/Tajima's D / Watterson's theta should be analytically
    # supplied by the wisp mask. FST is unaffected (it uses variant sites only).
    wisp_mask: Union[WispMask, None] = None
    # When True, pi is computed from per-sample genotype likelihoods (PL or GL FORMAT field)
    # instead of hard-called genotypes. `likelihood_field` records which field was found in
    # the VCF header (PL preferred when both are present). v1 supports diploid biallelic only
    # and `--stats pi` only.
    use_likelihoods: bool = False
    likelihood_field: Union[str, None] = None

    def __post_init__(self) -> None:
        """Checks a subset of mutually exclusive `pixy` args to ensure compliance."""
        if self.interval_start is None != self.interval_end is None:
            raise ValueError(
                "interval_start and interval_end must be specified together or not at all"
            )

        if self.bed is None == self.window_size is None:
            raise ValueError("One but not both of a BED file or a window size must be specified")

        if self.bed is not None and self.interval_start is not None:
            raise ValueError("An interval cannot be specified with a BED file")

        if self.interval_start is not None and self.interval_end is not None:
            if self.interval_start > self.interval_end:
                raise ValueError(
                    f"The specified interval start {self.interval_start} exceeds the specified "
                    f"interval end {self.interval_end}"
                )

    @property
    def has_interval(self) -> bool:
        """True if the `pixy` args included an interval."""
        # The post-init check enforces that the start exists if the end exists
        return self.interval_start is not None

    @cached_property
    def pop_names(self) -> NDArray[np.str_]:
        """
        Unique population names in first-appearance order, as a numpy `str_` array.

        Kept as numpy (rather than a Python tuple) for compatibility with downstream code that
        iterates the result with numpy semantics. The underlying ordering comes from
        `PopulationsTable.unique_populations`.
        """
        return np.array(self.populations.unique_populations, dtype=np.str_)

    @cached_property
    def pop_ids(self) -> NDArray[np.str_]:
        """Sample IDs in the populations table, as a numpy `str_` array."""
        return np.array(self.populations.ids, dtype=np.str_)


def _read_tsv_rows(
    path: Path,
    expected_cols: int,
    missing_data_msg: str,
    too_few_columns_msg: str,
) -> List[List[str]]:
    """
    Read a tab-delimited file into a list of stripped, validated rows.

    Each row must contain exactly `expected_cols` non-empty fields after stripping. Empty
    lines are skipped. The two message arguments mirror the wording the previous
    pandas-based validators emitted, so the regression tests' regex matches still work.
    """
    rows: List[List[str]] = []
    with open(path, "r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for raw in reader:
            # `csv.reader` will give us a one-element list for blank lines; treat as skip.
            if not raw or (len(raw) == 1 and not raw[0].strip()):
                continue
            cells = [c.strip() for c in raw]
            if len(cells) < expected_cols:
                # Distinguish "row too narrow" (entire columns missing) from "row has empty cells"
                if len(cells) == 1:
                    raise ValueError(too_few_columns_msg)
                raise ValueError(missing_data_msg)
            if any(not c for c in cells[:expected_cols]):
                raise ValueError(missing_data_msg)
            rows.append(cells[:expected_cols])
    return rows


def validate_populations_path(populations_path: Path) -> PopulationsTable:
    """
    Read and validate the --populations file.

    A valid populations file has 2 tab-separated columns; column 1 is the sample identifier
    and column 2 is the population label.

    Raises:
        FileNotFoundError: if the path does not exist
        ValueError: if any row is missing fields

    Returns:
        A `PopulationsTable`. `callset_index` is empty here — it gets populated by
        `check_and_validate_args` after the VCF header has been parsed.
    """
    if not os.path.exists(populations_path):
        raise FileNotFoundError(f"The specified populations file {populations_path} does not exist")

    rows = _read_tsv_rows(
        populations_path,
        expected_cols=2,
        missing_data_msg=(
            "The specified populations file contains missing data, "
            "confirm all samples have population IDs and are assigned to a population."
        ),
        too_few_columns_msg=(
            "Too many columns specified: expected 2 and found 1 "
            f"in populations file {populations_path}"
        ),
    )
    return PopulationsTable(
        ids=tuple(r[0] for r in rows),
        populations=tuple(r[1] for r in rows),
    )


def validate_bed_path(bed_path: Path) -> BedTable:
    """
    Read and validate a BED3 file (chrom, chromStart, chromEnd).

    BED coordinates are 0-based half-open by the UCSC spec: ``chromStart`` is
    inclusive and ``chromEnd`` is exclusive. Pixy uses 1-based inclusive
    coordinates internally (matching VCF/tabix region strings and the
    ``window_pos_1`` / ``window_pos_2`` columns it emits), so we convert at
    read time: ``window_pos_1 = chromStart + 1`` and ``window_pos_2 = chromEnd``.

    Raises:
        FileNotFoundError: If the path does not exist.
        ValueError: If any row has fewer than 3 fields, empty values, non-integer
            coordinates, a negative ``chromStart``, or ``chromEnd <= chromStart``.
    """
    if not os.path.exists(bed_path):
        raise FileNotFoundError(f"The specified BED file {bed_path} does not exist")

    rows = _read_tsv_rows(
        bed_path,
        expected_cols=3,
        missing_data_msg=(
            "The specified BED file contains missing data, confirm all rows have all "
            "three fields (chrom, pos1, pos2)."
        ),
        too_few_columns_msg=(
            f"Too many columns specified: expected 3 and found 1 in BED file {bed_path}"
        ),
    )
    chroms: List[str] = []
    starts: List[int] = []
    ends: List[int] = []
    for r in rows:
        try:
            chrom_start = int(r[1])
            chrom_end = int(r[2])
        except ValueError as e:
            raise ValueError(
                f"{bed_path}: BED columns 2 and 3 must be integers; got {r[1]!r}, {r[2]!r}"
            ) from e
        if chrom_start < 0:
            raise ValueError(
                f"{bed_path}: BED chromStart must be >= 0 (got {chrom_start} on "
                f"{r[0]}:{chrom_start}-{chrom_end})"
            )
        if chrom_end <= chrom_start:
            raise ValueError(
                f"{bed_path}: BED chromEnd must be > chromStart "
                f"(got {r[0]}:{chrom_start}-{chrom_end})"
            )
        chroms.append(r[0])
        # Convert BED 0-based half-open to pixy's 1-based inclusive convention.
        starts.append(chrom_start + 1)
        ends.append(chrom_end)
    return BedTable(chroms=tuple(chroms), chrom_starts=tuple(starts), chrom_ends=tuple(ends))


def validate_sites_path(sites_path: Path) -> SitesTable:
    r"""
    Read and validate the --sites_file (`CHROM\tPOS`).

    The previous implementation kept the raw rows in a DataFrame and filtered/sorted them
    later in each consumer. We pre-group rows by chromosome and pre-sort positions here, so
    consumers can fetch them in O(1).

    Raises:
        FileNotFoundError: if the path does not exist
        ValueError: if any row has fewer than 2 fields or empty values, or POS is non-integer.
    """
    if not os.path.exists(sites_path):
        raise FileNotFoundError(f"The specified sites file {sites_path} does not exist")

    rows = _read_tsv_rows(
        sites_path,
        expected_cols=2,
        missing_data_msg=(
            "The specified sites file contains missing data, confirm all rows "
            "have two fields (chrom, pos)."
        ),
        too_few_columns_msg="Too many columns specified: expected 2 and found 1",
    )
    by_chrom: Dict[str, List[int]] = {}
    for r in rows:
        chrom = r[0]
        try:
            pos = int(r[1])
        except ValueError as e:
            raise ValueError(
                f"{sites_path}: sites column 2 (POS) must be an integer; got {r[1]!r}"
            ) from e
        by_chrom.setdefault(chrom, []).append(pos)
    return SitesTable(
        by_chrom={chrom: tuple(sorted(positions)) for chrom, positions in by_chrom.items()}
    )


def validate_vcf_path(vcf_path: str) -> None:
    """
    Validates user-specified path to a VCF file.

    Additional validation of VCF file contents is currently performed outside of this function.

    Args:
        vcf_path: the path to the VCF file of interest

    Raises:
        Exception, if the specific VCF path does not exist
        Exception, if the VCF has no .gz extension (e.g., is not compressed with `bgzip`)
        Exception, if the VCF is not indexed (e.g., lacking either a `.csi` or `.tbi` extension)

    """
    if not os.path.exists(vcf_path):
        raise FileNotFoundError(f"The specified VCF {vcf_path} does not exist.")

    if not vcf_path.endswith(".gz"):
        raise ValueError(
            "The vcf is not compressed with bgzip (or has no .gz extension). "
            'To fix this, run "bgzip [filename].vcf" first (and then index with '
            '"tabix [filename].vcf.gz" if necessary)'
        )

    # Pick whichever index is present (prefer .tbi to match tabix's own preference) and bump
    # its mtime/atime to now so it never appears older than the VCF on disk (some pipelines
    # rsync the VCF after the index and tabix then complains). `os.utime` is the in-process
    # equivalent of `touch -c` and avoids a fork+exec on every pixy run.
    if os.path.exists(vcf_path + ".tbi"):
        index_path = vcf_path + ".tbi"
    elif os.path.exists(vcf_path + ".csi"):
        index_path = vcf_path + ".csi"
    else:
        raise ValueError(
            "The vcf is not indexed. Please either use `tabix` or `bcftools` to"
            "produce a `.tbi` or `.csi` index."
        )
    os.utime(index_path, None)


def validate_output_path(output_folder: str, output_prefix: str) -> Tuple[str, str]:
    """
    Validates user-specified output paths for `pixy` output.

    Args:
        output_folder: the directory to which to write any `pixy` results
        output_prefix: the combination of a given `output_folder` and `output_prefix`

    Raises:
        Exception, if the output folder is not writeable
        Exception, if the output prefix contains slashes

    Returns:
        output_folder and output_prefix

    """
    # attempt to create the output folder
    Path(output_folder).mkdir(parents=True, exist_ok=True)

    # check if output folder is writable
    # if not os.access(re.sub(r"[^\/]+$", "", args.outfile_prefix), os.W_OK):
    if not os.access(output_folder, os.W_OK):
        raise OSError(f"The output folder {output_folder} is not writable")

    # check if output_prefix is correctly specified
    if "/" in str(output_prefix) or "\\" in str(output_prefix):
        raise ValueError(
            f"The output prefix {output_prefix} contains slashes. "
            f"Remove them and specify output folder structure "
            "with --output_folder if necessary."
        )
    if output_folder != "":
        output_folder = output_folder + "/"
    else:
        output_folder = os.path.expanduser(os.getcwd() + "/")
    output_prefix = output_folder + output_prefix

    return output_folder, output_prefix


def get_chrom_list(args: argparse.Namespace) -> List[str]:
    """
    Get the list of chromosomes for analysis.

    If `--chromosomes all` is specified, this will be all chromosomes in the provided VCF.
    Otherwise, it will be the list of chromosomes provided to `--chromosomes`.

    Raises:
        Exception: If any chromosomes specified are not found in the VCF.
    """
    # get the list of all chromosomes in the dataset. `shell=False` avoids the extra `/bin/sh`
    # fork+exec on each call and removes shell-injection risk from user-supplied VCF paths.
    chrom_all = subprocess.check_output(["tabix", "-l", args.vcf]).decode("utf-8").split()
    if args.chromosomes != "all":
        # If a subset of chromosomes were specified, limit our analysis to those, and ensure that
        # they are all present in the VCF
        chrom_list = list(str(args.chromosomes).split(","))

        missing = list(set(chrom_list) - set(chrom_all))
        if len(missing) > 0:
            raise ValueError(
                f"The following chromosomes were specified but do not occur in the VCF: {missing}"
            )

    else:
        # Otherwise return everything in the VCF
        chrom_list = chrom_all

    return chrom_list


def validate_window_and_interval_args(
    args: argparse.Namespace, chrom_list: Union[List[str], None] = None
) -> str:
    """
    Validate the window and interval arguments when a BED file is not provided.

    Args:
        args: The parsed command-line arguments.
        chrom_list: Pre-computed list of chromosomes to analyze. When `None`, the function
            falls back to calling `get_chrom_list(args)` itself. Callers that already have a
            chrom list should pass it in to skip the extra `tabix -l` subprocess.

    Returns:
        A "check message", which is "OK" if all conditions are met, or a "WARNING" if the specified
        interval is smaller than the specified window size.
    """
    assert args.bed_file is None, (
        "this function should only be invoked when a BED file is not specified"
    )
    logger: logging.Logger = logging.getLogger(__name__)
    check_message: str = "OK"

    if args.window_size is None:
        raise ValueError("In the absence of a BED file, a --window_size must be specified.")

    if args.interval_start is None and args.interval_end is not None:
        raise ValueError(
            "When specifying an interval, both --interval_start and --interval_end are required."
        )

    if args.interval_start is not None and args.interval_end is None:
        raise ValueError(
            "When specifying an interval, both --interval_start and --interval_end are required."
        )

    if chrom_list is None:
        chrom_list = get_chrom_list(args)
    if (args.interval_start is not None or args.interval_end is not None) and len(chrom_list) > 1:
        raise ValueError(
            "--interval_start and --interval_end are not valid when calculating over "
            "multiple chromosomes. Remove both arguments or specify a single chromosome."
        )

    if (args.interval_start is not None and args.interval_end is not None) and (
        (int(args.interval_end) - int(args.interval_start)) <= int(args.window_size)
    ):
        check_message = "WARNING"
        logger.warning(
            f"The specified interval {args.interval_start}-{args.interval_end} "
            f"is smaller than the window size ({args.window_size}). "
            "A single window will be returned."
        )

    return check_message


# helper for parsing ploidy from a single GT string (e.g. "0/0", "0|1", "1")
def _read_vcf_samples_and_alts(  # noqa: C901
    vcf_path: str, scan_alts: bool, max_alt_records: int = 100_000
) -> Tuple[List[str], Union[set, None]]:
    """
    Read the VCF header for the sample list and, optionally, sample ALT values from records.

    A single ``gzip.open`` covers both passes — initializing zlib twice on the same file
    used to cost ~5–50 ms per run depending on header size and disk speed.

    The returned ``alt_values`` set follows the original invariant-check semantics: scan up
    to ``max_alt_records`` data records, but exit early once both the invariant marker
    (``.``) and at least one variant ALT have been seen, since at that point the downstream
    decision is fully determined. When ``scan_alts`` is False the second pass is skipped
    entirely and ``alt_values`` is returned as ``None``.

    Replaces the previous `allel.read_vcf_headers(...).samples`-based reader, which pulled
    in dask/zarr/pyarrow/scipy (~50 MB RSS, ~300 ms wall) just to grab the sample list.

    Raises:
        ValueError: if the file ends before the `#CHROM` column line is found.
    """
    samples: Union[List[str], None] = None
    alt_values: Union[set, None] = set() if scan_alts else None
    seen_alt_records = 0
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if samples is None:
                if line.startswith("#CHROM"):
                    fields = line.rstrip("\n").split("\t")
                    # standard VCF: 9 fixed columns then per-sample columns
                    samples = fields[9:]
                    if not scan_alts:
                        break
                    continue
                if not line.startswith("#"):
                    # Past the header without seeing #CHROM — malformed VCF
                    break
                continue
            # samples already captured; we're now in the data section reading ALT values
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t", 5)
            if len(fields) >= 5:
                # `alt_values` is non-None here because `scan_alts` was True (we'd have
                # broken out above otherwise); the cast keeps mypy quiet without runtime cost.
                assert alt_values is not None
                alt_values.add(fields[4])
            seen_alt_records += 1
            # Once we've seen both an invariant ALT and a non-invariant ALT, the downstream
            # check is fully decided and there's no point reading further.
            assert alt_values is not None
            if "." in alt_values and (alt_values - {"."}):
                break
            if seen_alt_records >= max_alt_records:
                break

    if samples is None:
        raise ValueError(
            f"VCF {vcf_path!r} has no `#CHROM` header line; cannot determine sample names."
        )
    return samples, alt_values


def _read_vcf_format_ids(vcf_path: str) -> set:
    """
    Return the set of FORMAT IDs declared in the VCF header.

    Used by the `--use_likelihoods` validation to pick between PL and GL up front, and to
    fail fast when neither is declared. One small gzip pass over the header lines.
    """
    format_ids: set = set()
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("##FORMAT="):
                # Header line shape: ##FORMAT=<ID=PL,Number=G,Type=Integer,...>
                # Locate the ID= field without parsing the full bracketed expression.
                tag = "ID="
                idx = line.find(tag)
                if idx < 0:
                    continue
                rest = line[idx + len(tag) :]
                # ID value ends at the next comma or '>'.
                end = len(rest)
                for ch in (",", ">"):
                    j = rest.find(ch)
                    if 0 <= j < end:
                        end = j
                format_ids.add(rest[:end].strip())
            elif line.startswith("#CHROM"):
                break
            elif not line.startswith("#"):
                # Past the header without seeing #CHROM — malformed VCF, but not our concern here.
                break
    return format_ids


def _ploidy_from_gt(gt_str: str) -> int:
    """Infer ploidy from a single GT genotype string."""
    if "/" in gt_str:
        alleles = gt_str.split("/")
    elif "|" in gt_str:
        alleles = gt_str.split("|")
    else:
        # No separator: e.g., "0", ".", "1"
        alleles = [gt_str]

    # Count non-missing alleles
    non_missing = [a for a in alleles if a != "."]

    if len(alleles) == 1 or (len(non_missing) == 1 and len(alleles) > 1):
        return 1
    return len(alleles)


# function for inferring ploidy per contig from a VCF, using tabix to read the
# first record for each contig. Returns a dict mapping contig name -> ploidy.
# Supports VCFs with variable ploidy across contigs (e.g. autosomes vs. sex
# chromosomes or organellar contigs).
def infer_ploidy_per_contig(vcf_path: str, chrom_list: List[str]) -> "dict[str, int]":  # noqa: C901
    """
    Infer ploidy per contig by reading the first record of each contig via tabix.

    Args:
        vcf_path: path to a bgzipped, tabix-indexed VCF
        chrom_list: contigs over which to infer ploidy (typically the list of contigs
            that will be analyzed)

    Returns:
        A dict mapping each contig name to its inferred ploidy.

    Raises:
        RuntimeError: if a contig has no records in the VCF.
    """
    if not chrom_list:
        return {}

    # We spawn one tabix process and pass every contig as a region argument; tabix outputs
    # the records sequentially. The previous version spawned one tabix per contig — for
    # assemblies with hundreds of unplaced contigs that subprocess overhead dominated startup.
    #
    # For each contig we walk all sample columns (not just the first) and scan up to
    # max_sites records looking for a sample with ploidy > 1. Single-dot missing GTs
    # (which _ploidy_from_gt maps to 1) are naturally skipped by the p > 1 guard.
    max_sites = 100_000
    ploidy_map: dict[str, int] = {}
    remaining = set(chrom_list)
    site_count: dict[str, int] = {}
    with subprocess.Popen(
        ["tabix", vcf_path, *chrom_list],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    ) as proc:
        assert proc.stdout is not None
        for line in proc.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                # Not a usable record (no sample columns). Skip — error is raised below if no
                # usable record is ever seen for a contig.
                continue
            chrom = fields[0]
            if chrom not in remaining:
                continue
            found_diploid = False
            for sample_field in fields[9:]:
                p = _ploidy_from_gt(sample_field.split(":")[0])
                if p > 1:
                    ploidy_map[chrom] = p
                    remaining.discard(chrom)
                    found_diploid = True
                    break
            if not found_diploid:
                site_count[chrom] = site_count.get(chrom, 0) + 1
                if site_count[chrom] >= max_sites:
                    ploidy_map[chrom] = 1  # genuinely haploid
                    remaining.discard(chrom)
            if not remaining:
                break
        proc.terminate()

    # Chroms whose records were all haploid-looking but the stream ended before max_sites:
    # we've seen at least one record, so they're genuinely haploid rather than absent.
    for chrom in [c for c in remaining if c in site_count]:
        ploidy_map[chrom] = 1
        remaining.discard(chrom)

    if remaining:
        # Only chroms with NO records at all reach here.
        missing = [c for c in chrom_list if c in remaining]
        raise RuntimeError(
            f"No genotype records found in VCF for contig(s) {missing!r}; cannot infer ploidy."
        )

    return ploidy_map


def check_and_validate_args(  # noqa: C901
    args: argparse.Namespace,
) -> PixyArgs:
    """
    Checks whether user-specific arguments are valid.

    Args:
        args: parsed CLI args specified by the user

    Raises:
        Exception: if the output_folder is not writeable
        Exception: if any of the required input files are missing
        Exception: if the output_prefix contains either forward or backward slashes
        Exception: if any of the provided files do not exist
        Exception: if the VCF file is not compressed with `bgzip` or indexed with `tabix`
        Exception: if the VCF file does not contain variant sites
        Exception: if the provided `--chromosomes` do not occur in the specified VCF file
        Exception: if neither a `--bed-file` nor a `--window_size` is provided
        Exception: if only one of `--interval_start` and `--interval_end` is given
        Exception: if multiple `--chromosomes` and an interval are provided
        Exception: if any rows in either the `--bed-file` or `--sites-path` are missing data
        Exception: if any of the samples provided in the `populations_file` do not exist in the VCF
        Exception: if the `populations_file` does not contain at least 2 populations

    Returns:
        an instance of PixyArgs where each attribute is validated from user-specified input

    """
    # CHECK FOR TABIX
    tabix_path = shutil.which("tabix")
    logger = logging.getLogger(__name__)
    if tabix_path is None:
        raise ValueError(
            "`tabix` is not installed (or cannot be located in the path). "
            'Install tabix with "conda install -c bioconda htslib".'
        )

    if args.vcf is None:
        raise ValueError(f"The --vcf argument is missing or incorrectly specified: {args.vcf}")

    if args.populations is None:
        raise ValueError(
            f"The --populations argument is missing or incorrectly specified: {args.populations}."
        )

    # reformat file paths for compatibility

    populations_path: Path = Path(os.path.expanduser(args.populations))
    populations: PopulationsTable = validate_populations_path(populations_path)

    vcf_path: str = os.path.expanduser(args.vcf)  # we don't want a Path object just yet because
    # most of the downstream operations require a string
    validate_vcf_path(vcf_path)

    # Whether the invariant-ALT scan can be skipped is known up front from the user-supplied
    # flags: --bypass_invariant_check explicitly opts out, and supplying a --wisp_bed implies
    # it (the wisp mask supplies the invariant denominator analytically). Computing this here
    # lets the single-pass VCF reader below skip the data-section scan when it isn't needed.
    bypass_invariant_check: bool = args.bypass_invariant_check or (
        getattr(args, "wisp_bed", None) is not None
    )
    # One gzip pass for both the sample list (always needed) and the ALT-value scan (only
    # when invariant-checking is on). Previously these were two separate `gzip.open` calls
    # that each re-decompressed the VCF header.
    vcf_samples, alt_values = _read_vcf_samples_and_alts(
        vcf_path, scan_alts=not bypass_invariant_check
    )

    logger.info("Validating VCF and input parameters...")

    # CHECK OUTPUT FOLDER
    logger.info("Checking write access...")

    output_folder, output_prefix = validate_output_path(
        output_folder=args.output_folder, output_prefix=args.output_prefix
    )

    # CHECK CPU CONFIGURATION
    logger.info("Checking CPU configuration...")

    available_cores: int = os.cpu_count() or 1
    if args.n_cores > available_cores:
        logger.warning(
            f"{args.n_cores} CPU cores requested but only {available_cores} available. "
            f"Using {available_cores} cores."
        )
        args.n_cores = available_cores

    # CHECK FOR EXISTENCE OF INPUT FILES

    bed: Union[BedTable, None]
    if args.bed_file is not None:
        bed_path: Path = Path(os.path.expanduser(args.bed_file))
        bed = validate_bed_path(bed_path)
    else:
        bed = None

    # WISP MASK
    # If a wisp BED is supplied, --vcf is variants-only; the invariant-site denominator
    # comes from the wisp mask and the invariant-presence check on the VCF is moot
    # (`bypass_invariant_check` was already set True above when --wisp_bed was supplied).
    wisp_mask: Union["WispMask", None]
    wisp_bed_arg: Union[str, None] = getattr(args, "wisp_bed", None)
    if wisp_bed_arg is not None:
        # Defer wisp imports: they pull in `json` and the wisp dataclasses, which add ~14 ms
        # to startup but are only needed when --wisp_bed is supplied. Most pixy runs don't
        # use a wisp mask, so this keeps their args-validation step lighter.
        from pixy.wisp import WispMask
        from pixy.wisp import validate_wisp_against_populations
        from pixy.wisp import validate_wisp_data_rows

        wisp_path: Path = Path(os.path.expanduser(wisp_bed_arg))
        if not os.path.exists(wisp_path):
            raise FileNotFoundError(f"The specified wisp BED file {wisp_path} does not exist")
        # Tabix index required for per-window queries.
        tbi = Path(str(wisp_path) + ".tbi")
        csi = Path(str(wisp_path) + ".csi")
        if not (tbi.exists() or csi.exists()):
            raise ValueError(
                f"The wisp BED {wisp_path} is not indexed. "
                "Index it with `tabix -p bed [filename].bed.gz` first."
            )
        wisp_mask = WispMask.from_path(wisp_path)
        # Cross-check populations and per-pop sample counts against --populations.
        pop_to_sample_count: Dict[str, int] = {}
        for pop in populations.populations:
            pop_to_sample_count[pop] = pop_to_sample_count.get(pop, 0) + 1
        problem = validate_wisp_against_populations(
            wisp_mask=wisp_mask, pop_to_sample_count=pop_to_sample_count
        )
        if problem is not None:
            raise ValueError(problem)
        # Confirm the wisp BED actually carries data rows and that the first few
        # rows parse cleanly. Catches truncated / accidentally header-only wisp
        # outputs and malformed counts before pixy hits them mid-window.
        problem = validate_wisp_data_rows(wisp_mask=wisp_mask)
        if problem is not None:
            raise ValueError(problem)
        logger.info(
            f"Using wisp mask {wisp_path} "
            f"(threshold={wisp_mask.metadata.threshold}, "
            f"populations={list(wisp_mask.metadata.populations)})"
        )
    else:
        wisp_mask = None

    # VALIDATE THE VCF

    # check if the vcf contains any invariant sites
    # a very basic check: just looks for at least one invariant site in the alt field
    logger.info("Checking for invariant sites...")
    # `alt_values` was populated by the up-front single-pass VCF read above; it is None
    # exactly when `bypass_invariant_check` was already true, i.e. the scan was skipped.
    if alt_values is not None:
        if "." not in alt_values:
            raise ValueError(
                "The provided VCF appears to contain no invariant sites "
                '(ALT = "."). '
                "This check can be bypassed via --bypass_invariant_check 'yes'."
            )
        if alt_values == {"."}:
            logger.warning(
                "The provided VCF appears to contain no variable sites in the "
                "first 100 000 sites. It may have been filtered incorrectly, or genetic diversity "
                "may be extremely low. "
                "This warning can be suppressed via --bypass_invariant_check 'yes'.'"
            )
    else:
        # When the wisp path is in use, the invariant denominator is supplied by the
        # mask rather than by VCF invariants — the "incorrect estimates" warning would be
        # misleading there. Keep the existing warning for the manual --bypass case.
        if wisp_mask is None and not (len(args.stats) == 1 and (args.stats[0] == "fst")):
            logger.warning(
                "EXTREME WARNING: --bypass_invariant_check is set to True. Note that a "
                "lack of invariant sites will result in incorrect estimates."
            )

    # check if requested chromosomes exist in vcf
    logger.info("Checking chromosome data...")

    chrom_list: List[str] = get_chrom_list(args)

    # INTERVALS
    # check if intervals are correctly specified
    # validate the BED file (if present)

    logger.info("Checking intervals/sites...")

    if args.bed_file is None:
        check_message = validate_window_and_interval_args(args, chrom_list=chrom_list)
        logger.info(check_message)
    else:
        if (
            args.interval_start is not None
            or args.interval_end is not None
            or args.window_size is not None
        ):
            raise ValueError(
                "--interval_start, --interval_end, and --window_size are not valid "
                "when a BED file of windows is provided."
            )

        assert bed is not None  # narrow type for mypy
        bed_chrom: List[str] = list(bed.chroms)
        missing = list(set(bed_chrom) - set(chrom_list))
        chrom_list = list(set(chrom_list) & set(bed_chrom))

        if len(missing) > 0:
            logger.warning(
                "The following chromosomes are in the BED file but do not occur in the VCF "
                f"and will be ignored: {missing}"
            )
    sites: Union[SitesTable, None]
    if args.sites_file is None:
        sites = None
        chrom_sites: List[str] = []
        missing_sites: List[str] = []
    else:
        sites_path: Path = Path(os.path.expanduser(args.sites_file))
        sites = validate_sites_path(sites_path=sites_path)

        # all the chromosomes in the sites file
        chrom_sites = list(sites.chromosomes)

        # the difference between the chromosomes in the sites file and the VCF
        missing_sites = list(set(chrom_sites) - set(chrom_list))

        if len(missing_sites) > 0:
            logger.warning(
                "The following chromosomes occur in the sites file but do not occur in the "
                f"VCF and will be ignored: {missing_sites}"
            )

    # SAMPLES
    # check if requested samples exist in vcf

    logger.info("Checking sample data...")

    # - parse + validate the population file
    # - format is IND POP (tab separated)
    # - throws an error if individuals are missing from VCF

    # make sure every indiv in the pop file is in the VCF callset
    sample_ids: List[str] = list(populations.ids)
    missing = list(set(sample_ids) - set(vcf_samples))

    # find the samples in the callset index by matching up the order of samples between the
    # population file and the callset
    # also check if there are invalid samples in the popfile
    try:
        samples_callset_index = tuple(vcf_samples.index(s) for s in populations.ids)
    except ValueError as e:
        raise ValueError(
            f"The following samples are listed in the population file but not in the VCF: {missing}"
        ) from e
    # Re-build the populations table with the discovered callset indices attached.
    populations = populations.with_callset_index(samples_callset_index)

    if populations.num_populations == 1 and ("fst" in args.stats or "dxy" in args.stats):
        raise ValueError(
            "Calculation of fst and/or dxy requires at least two populations to be "
            "defined in the population file."
        )

    include_multiallelic_snps: bool = args.include_multiallelic_snps

    # GENOTYPE-LIKELIHOOD MODE
    # When --use_likelihoods is set, pi is estimated from per-sample PL/GL FORMAT fields
    # instead of hard calls. v1 supports diploid biallelic sites only and `--stats pi` only.
    use_likelihoods: bool = getattr(args, "use_likelihoods", False)
    likelihood_field: Union[str, None] = None
    if use_likelihoods:
        non_pi = [s for s in args.stats if s != "pi"]
        if non_pi:
            raise ValueError(
                "--use_likelihoods currently supports only --stats pi "
                "(dxy, fst, watterson_theta, tajima_d from likelihoods are planned for "
                f"a later release). Remove from --stats: {non_pi}"
            )
        if include_multiallelic_snps:
            raise ValueError(
                "--use_likelihoods supports diploid biallelic sites only in v1; "
                "remove --include_multiallelic_snps."
            )
        # Probe the VCF header for PL / GL FORMAT declarations and pick one.
        format_ids = _read_vcf_format_ids(vcf_path)
        if "PL" in format_ids:
            likelihood_field = "PL"
        elif "GL" in format_ids:
            likelihood_field = "GL"
        else:
            raise ValueError(
                f"--use_likelihoods requires the VCF {vcf_path!r} to declare a PL or GL "
                "FORMAT field; neither was found in the header."
            )
        logger.info(f"Using genotype likelihoods from FORMAT field: {likelihood_field}")

    # check ploidy per contig (supports VCFs with variable ploidy across contigs)
    ploidy_map: Dict[str, int] = infer_ploidy_per_contig(vcf_path, chrom_list)
    distinct_ploidies = sorted(set(ploidy_map.values()))
    if len(distinct_ploidies) == 1:
        logger.info(f"Inferred ploidy: {distinct_ploidies[0]} (uniform across contigs)")
    else:
        logger.info(f"Inferred variable ploidy across contigs: {ploidy_map}")

    # WC-FST is only supported for diploid data; warn up front if any non-diploid contig
    # will be encountered while WC-FST is requested (covers both mixed- and uniform-ploidy VCFs).
    if "fst" in args.stats and args.fst_type.upper() == "WC":
        non_diploid = [c for c, p in ploidy_map.items() if p != 2]
        if non_diploid:
            logger.warning(
                "Weir-Cockerham FST is not supported for non-diploid contigs. "
                f"FST will be skipped for: {non_diploid}. "
                "Use --fst_type hudson to compute FST on these contigs."
            )

    # --use_likelihoods is diploid-only in v1; reject up front if any analyzed contig is not.
    if use_likelihoods:
        non_diploid = [c for c, p in ploidy_map.items() if p != 2]
        if non_diploid:
            raise ValueError(
                "--use_likelihoods supports diploid contigs only in v1; "
                f"non-diploid contigs were detected: {non_diploid}."
            )

    logger.info("All initial checks passed!")
    stats: List[PixyStat] = [PixyStat[stat.upper()] for stat in args.stats]
    tmp_path: Path = _generate_tmp_path(output_dir=output_folder)
    _check_tmp_path(tmp_path)

    return PixyArgs(
        stats=stats,
        vcf_path=Path(vcf_path),
        populations=populations,
        num_cores=args.n_cores,
        bypass_invariant_check=bypass_invariant_check,
        include_multiallelic_snps=include_multiallelic_snps,
        bed=bed,
        output_dir=Path(output_folder),
        output_prefix=output_prefix,
        chromosomes=chrom_list,
        window_size=args.window_size,
        interval_start=args.interval_start,
        interval_end=args.interval_end,
        sites=sites,
        chunk_size=args.chunk_size,
        fst_type=FSTEstimator[args.fst_type.upper()],
        fst_components=getattr(args, "fst_components", False),
        tajima_components=getattr(args, "tajima_components", False),
        temp_file=tmp_path,
        ploidy_map=ploidy_map,
        wisp_mask=wisp_mask,
        use_likelihoods=use_likelihoods,
        likelihood_field=likelihood_field,
    )


def _generate_tmp_path(output_dir: str) -> Path:
    """Generates a temporary file path to which `pixy` will write intermediate results."""
    return Path(output_dir) / f"pixy_tmpfile_{uuid.uuid4().hex}.tmp"


def _check_tmp_path(temp_file: Path) -> None:
    # check if temp file is writable
    with open(temp_file, "w"):
        pass  # file is created and then closed
    assert os.access(temp_file, os.W_OK), "temp file is not writable"
