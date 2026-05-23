# `from __future__ import annotations` lazies all annotations into strings so referencing
# `pandas.DataFrame` / `allel.*` in type hints doesn't force those modules to be imported at
# module load. The actual `import pandas` / `import allel` calls are deferred to inside the
# helper functions that need them — this keeps `pixy --help`, `--version`, and arg-parse
# error paths from paying the ~500 ms / ~130 MB cost of those imports.
from __future__ import annotations

import argparse
import gzip
import logging
import os
import shutil
import subprocess
import uuid
from dataclasses import dataclass
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
    import pandas


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
        populations_df: a pandas DataFrame derived from a user-specified path to a headerless
            tab-separated populations file
        num_cores: number of CPUs to utilize for parallel processing (default = 1)
        include_multiallelic_snps: If True, include multiallelic sites in the analysis
        bypass_invariant_check: whether to allow computation of stats without invariant sites
            (this option is never recommended and defaults to False)
        bed_df: a pandas DataFrame derived from a user-specified path to a headerless BED file.
            Empty DataFrame if a path is not given.
        output_dir: an optional path to which outputs will be written; default is current directory
        output_prefix: an optional prefix with which to prepend to `pixy` output (default is `pixy`)
        chromosomes: an optional comma-separated list of chroms over which stats will be calculated
            (defaults to all chromosomes in a given VCF)
        window_size: an optional length of base pairs over which to calculate stats
        interval_start: an optional 1-based position demarcating the start of an interval over which
            to calculate stats (only valid when calculating over a single chromosome)
        interval_end: an optional 1-based position demarcating the end of an interval over which
            to calculate stats (only valid when calculating over a single chromosome)
        sites_df: a pandas DataFrame derived from a user-specified path to a headerless
            tab-separated file that defines specific sites for which summary stats should be
            calculated. Empty DataFrame if a path is not given.
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

    Raises:
        ValueError: if an interval is specified without both start and end positions
        ValueError: if a BED file is specified, and a window size or interval are also specified
        ValueError: if a BED file is not specified, and a window size is not specified
    """

    stats: List[PixyStat]
    vcf_path: Path
    populations_df: pandas.DataFrame
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
    bed_df: Union[pandas.DataFrame, None] = None
    window_size: Union[int, None] = None
    interval_start: Union[int, None] = None
    interval_end: Union[int, None] = None
    sites_df: Union[pandas.DataFrame, None] = None
    ploidy_map: Union[Dict[str, int], None] = None

    def __post_init__(self) -> None:
        """Checks a subset of mutually exclusive `pixy` args to ensure compliance."""
        if self.interval_start is None != self.interval_end is None:
            raise ValueError(
                "interval_start and interval_end must be specified together or not at all"
            )

        if self.bed_df is None == self.window_size is None:
            raise ValueError("One but not both of a BED file or a window size must be specified")

        if self.bed_df is not None and self.interval_start is not None:
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
    def pop_names(self) -> NDArray[np.bytes_]:
        """
        Returns the list of unique population names from the provided `populations_df`.

        NB, we cast to a datatype of np.str_ for `mypy`, which cannot infer the specific datatypes
        within the array.
        """
        return np.array(self.populations_df["Population"].unique(), dtype=np.str_)

    @cached_property
    def pop_ids(self) -> NDArray[np.bytes_]:
        """
        Returns the list of unique population identifiers from the provided `populations_df`.

        NB, we cast to a datatype of np.str_ for `mypy`, which cannot infer the specific datatypes
        within the array.
        """
        return np.array(self.populations_df["ID"].unique(), dtype=np.str_)


def validate_populations_path(populations_path: Path) -> pandas.DataFrame:
    """
    Validates user-specified path to a populations file and its contents.

    A valid populations file has at least 2 columns in it, the first of which must be a sample
    identifier. The second column must be the corresponding population to which that sample id
    belongs. The sample identifier is assumed to be alphanumeric.

    Raises:
        Exception, if the specific populations_path does not exist
        Exception, if any of the rows in the populations_path is missing data

    Returns:
        A pandas DataFrame containing the populations file contents, with one row for each row
        in the input populations file.
        This DataFrame includes two columns ("ID": str, "Population": str).
    """
    # read in the list of samples/populations
    if not os.path.exists(populations_path):
        raise FileNotFoundError(f"The specified populations file {populations_path} does not exist")

    # Local import — see module-level note. pandas is ~80 MB of RSS and ~500 ms to load; we
    # only pay that when the user actually has a populations file to parse.
    import pandas

    poppanel: pandas.DataFrame = pandas.read_csv(
        populations_path, sep="\t", usecols=[0, 1], names=["ID", "Population"]
    )

    # remove trailing whitespace from population names (if any)
    poppanel = poppanel.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

    poppanel["ID"] = poppanel["ID"].astype(str)

    # check for missing values

    if poppanel.isnull().values.any():
        raise ValueError(
            "The specified populations file contains missing data, "
            "confirm all samples have population IDs are are assigned to a population."
        )

    return poppanel


def validate_bed_path(bed_path: Path) -> pandas.DataFrame:
    """
    Validate and load the user-specified path to a BED file.

    The validation checks that the BED file exists and contains three columns. Chromosome names are
    coerced to string before returning.

    Args:
        bed_path: A path to a BED file.
            This file must be in BED3 format (three columns; chrom, start, and end).

    Returns:
        A pandas DataFrame containing the BED file contents, with one row for each row in the BED
        file.

        This DataFrame includes three columns ("chrom": str, "pos1": int, "pos2": int).

    Raises:
        Exception: If the provided filepath does not exist.
        Exception: If the BED file contains any rows with missing fields.
            Every row must include three fields - chrom, start, and stop.
        Exception: If the BED file is not in BED3 format; i.e. if the BED does not contain at
            least three columns.
    """
    if not os.path.exists(bed_path):
        raise FileNotFoundError(f"The specified BED file {bed_path} does not exist")

    import pandas  # deferred — see module-level note

    # read in the bed file, strip trailing whitespace, and extract the chromosome column
    bed_df = pandas.read_csv(bed_path, sep="\t", usecols=[0, 1, 2], names=["chrom", "pos1", "pos2"])
    bed_df = bed_df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

    # force chromosomes to strings
    bed_df["chrom"] = bed_df["chrom"].astype(str)

    if bed_df.isnull().values.any():
        raise ValueError(
            "The specified BED file contains missing data, confirm all rows have all "
            "three fields (chrom, pos1, pos2)."
        )

    if len(bed_df.columns) != 3:
        raise ValueError(f"The specified BED file has {len(bed_df.columns)} columns; expected 3.")

    else:
        bed_df.columns = ["chrom", "chromStart", "chromEnd"]

    return bed_df


def validate_sites_path(sites_path: Path) -> pandas.DataFrame:
    """
    Validates user-specified path to a sites file, if provided to `pixy`.

    A valid sites file has no header and at least 2 columns in it. The first column must be a
    chromosome and the second column must be a position. Position is assumed to be an int.
    Chromosome names are coerced to a string before returning the dataframe.

    Raises:
        Exception, if the specific sites path does not exist
        Exception, if any of the rows in the sites_path is missing data

    Returns:
        A pandas DataFrame containing the sites file contents, with one row for each row
        in the input sites file. The dataframe has two columns ("CHROM": str, "POS": str).
    """
    if not os.path.exists(sites_path):
        raise FileNotFoundError(f"The specified sites file {sites_path} does not exist")

    import pandas  # deferred — see module-level note

    sites_df = pandas.read_csv(sites_path, sep="\t", usecols=[0, 1], names=["chrom", "pos"])
    sites_df = sites_df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    sites_df["chrom"] = sites_df["chrom"].astype(str)

    if sites_df.isnull().values.any():
        raise ValueError(
            "The specified sites file contains missing data, confirm all rows each "
            "have two fields (chrom, pos)."
        )

    if len(sites_df.columns) != 2:
        raise ValueError(f"The specified BED file has {len(sites_df.columns)} columns; expected 2.")

    else:
        sites_df.columns = ["CHROM", "POS"]

    return sites_df


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

    if not (os.path.exists(vcf_path + ".tbi") or os.path.exists(vcf_path + ".csi")):
        raise ValueError(
            "The vcf is not indexed. Please either use `tabix` or `bcftools` to"
            "produce a `.tbi` or `.csi` index."
        )

    # update the date of the index
    # this apparently resolves issues of indexes being older than vcfs in some pipelines
    if os.path.exists(vcf_path + ".tbi"):
        subprocess.run(["touch", "-c", str(vcf_path + ".tbi")])
    else:
        subprocess.run(["touch", "-c", str(vcf_path + ".csi")])


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
    # get the list of all chromosomes in the dataset
    chrom_all = subprocess.check_output("tabix -l " + args.vcf, shell=True).decode("utf-8").split()
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


def validate_window_and_interval_args(args: argparse.Namespace) -> str:
    """
    Validate the window and interval arguments when a BED file is not provided.

    Args:
        args: The parsed command-line arguments.

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

    chrom_list: List[str] = get_chrom_list(args)
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
def _read_vcf_samples(vcf_path: str) -> List[str]:
    """Read the sample names from a (bgzipped) VCF's `#CHROM` header line.

    Replaces `allel.read_vcf_headers(...).samples` — pixy only needs the sample list, and
    going through scikit-allel forces the import of dask/zarr/pyarrow/scipy (~50 MB RSS,
    ~300 ms wall) every time the validator runs. This helper streams just the VCF header
    block and stops as soon as it finds the column line.
    """
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                # standard VCF: 9 fixed columns then per-sample columns
                return fields[9:]
            if not line.startswith("#"):
                # Past the header without seeing #CHROM — malformed VCF
                break
    raise ValueError(
        f"VCF {vcf_path!r} has no `#CHROM` header line; cannot determine sample names."
    )


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
def infer_ploidy_per_contig(vcf_path: str, chrom_list: List[str]) -> "dict[str, int]":
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
    # the records sequentially. As soon as we've seen the first record for each contig we
    # terminate. The previous version spawned one tabix per contig — for assemblies with
    # hundreds of unplaced contigs that subprocess overhead dominated startup time.
    ploidy_map: dict[str, int] = {}
    remaining = set(chrom_list)
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
            first_gt = fields[9].split(":")[0]
            ploidy_map[chrom] = _ploidy_from_gt(first_gt)
            remaining.discard(chrom)
            if not remaining:
                break
        proc.terminate()

    if remaining:
        # Order the error report by the user's original chrom_list for reproducible messages.
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
    populations_df: pandas.DataFrame = validate_populations_path(populations_path)

    vcf_path: str = os.path.expanduser(args.vcf)  # we don't want a Path object just yet because
    # most of the downstream operations require a string
    validate_vcf_path(vcf_path)

    if args.output_folder != "":
        output_folder = args.output_folder + "/"
    else:
        output_folder = os.path.expanduser(os.getcwd() + "/")

    output_prefix = output_folder + args.output_prefix

    # Read the VCF header to get the sample list. Previously this used
    # `allel.read_vcf_headers`, which transitively imports dask/zarr/pyarrow and adds ~300 ms
    # / ~50 MB to startup. The sample list is the only piece of header info we use here, and
    # the `#CHROM` column line that carries it is trivial to parse directly.
    vcf_samples: List[str] = _read_vcf_samples(vcf_path)

    logger.info("Validating VCF and input parameters...")

    # CHECK OUTPUT FOLDER
    logger.info("Checking write access...")

    output_folder, output_prefix = validate_output_path(
        output_folder=args.output_folder, output_prefix=args.output_prefix
    )

    # CHECK CPU CONFIGURATION
    logger.info("Checking CPU configuration...")
    check_message = "OK"

    available_cores: int = os.cpu_count() or 1
    if args.n_cores > available_cores:
        logger.warning(
            f"{args.n_cores} CPU cores requested but only {available_cores} available. "
            f"Using {available_cores} cores."
        )
        args.n_cores = available_cores

    # CHECK FOR EXISTENCE OF INPUT FILES

    if args.bed_file is not None:
        bed_path: Path = Path(os.path.expanduser(args.bed_file))
        bed_df: pandas.DataFrame = validate_bed_path(bed_path)

    else:
        bed_df = None

    # VALIDATE THE VCF

    # check if the vcf contains any invariant sites
    # a very basic check: just looks for at least one invariant site in the alt field
    logger.info("Checking for invariant sites...")
    check_message = "OK"
    bypass_invariant_check: bool = args.bypass_invariant_check
    if not bypass_invariant_check:
        # Read up to the first 100k data records and collect the distinct ALT values.
        # We bail out as soon as we have seen the invariant marker ('.') alongside at least one
        # variant ALT, since at that point both branches of the check below are already decided.
        # Implemented in pure Python (no shell pipeline) to avoid shell-injection risk from
        # user-supplied VCF paths, and to remove the dependency on system `gunzip`/`awk`/`sort`.
        alt_values: set[str] = set()
        max_records = 100_000
        seen_records = 0
        with gzip.open(vcf_path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t", 5)
                if len(fields) >= 5:
                    alt_values.add(fields[4])
                seen_records += 1
                # Once we've seen both an invariant ALT and a non-invariant ALT, the check below
                # is fully decided and there's no point reading further.
                if "." in alt_values and (alt_values - {"."}):
                    break
                if seen_records >= max_records:
                    break

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
        if not (len(args.stats) == 1 and (args.stats[0] == "fst")):
            logger.warning(
                "EXTREME WARNING: --bypass_invariant_check is set to True. Note that a "
                "lack of invariant sites will result in incorrect estimates."
            )

    # check if requested chromosomes exist in vcf
    # parses the whole CHROM column (!)

    logger.info("Checking chromosome data...")

    chrom_list: List[str] = get_chrom_list(args)

    # INTERVALS
    # check if intervals are correctly specified
    # validate the BED file (if present)

    logger.info("Checking intervals/sites...")
    check_message = "OK"

    if args.bed_file is None:
        check_message = validate_window_and_interval_args(args)
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

        bed_chrom: List[str] = list(bed_df["chrom"])
        missing = list(set(bed_chrom) - set(chrom_list))
        chrom_list = list(set(chrom_list) & set(bed_chrom))

        if len(missing) > 0:
            logger.warning(
                "The following chromosomes are in the BED file but do not occur in the VCF "
                f"and will be ignored: {missing}"
            )
    if args.sites_file is None:
        sites_df = None
        chrom_sites: List[str] = []
        missing_sites: List[str] = []

    else:
        sites_path: Path = Path(os.path.expanduser(args.sites_file))
        sites_df = validate_sites_path(sites_path=sites_path)

        # all the chromosomes in the sites file
        chrom_sites = list(sites_df["CHROM"])

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

    # get a list of samples from the callset
    samples_list = vcf_samples
    # make sure every indiv in the pop file is in the VCF callset
    sample_ids: List[str] = list(populations_df["ID"])
    missing = list(set(sample_ids) - set(samples_list))

    # find the samples in the callset index by matching up the order of samples between the
    # population file and the callset
    # also check if there are invalid samples in the popfile
    try:
        samples_callset_index = [samples_list.index(s) for s in populations_df["ID"]]
    except ValueError as e:
        raise ValueError(
            f"The following samples are listed in the population file but not in the VCF: {missing}"
        ) from e
    else:
        populations_df["callset_index"] = samples_callset_index

        # use the popindices dictionary to keep track of the indices for each population
        popindices = {}
        popnames = populations_df.Population.unique()
        for name in popnames:
            popindices[name] = populations_df[
                populations_df.Population == name
            ].callset_index.values

    if (populations_df["Population"].nunique() == 1) and (
        "fst" in args.stats or "dxy" in args.stats
    ):
        raise ValueError(
            "Calculation of fst and/or dxy requires at least two populations to be "
            "defined in the population file."
        )

    include_multiallelic_snps: bool = args.include_multiallelic_snps

    # check ploidy per contig (supports VCFs with variable ploidy across contigs)
    ploidy_map: Dict[str, int] = infer_ploidy_per_contig(vcf_path, chrom_list)
    distinct_ploidies = sorted(set(ploidy_map.values()))
    if len(distinct_ploidies) == 1:
        logger.info(f"Inferred ploidy: {distinct_ploidies[0]} (uniform across contigs)")
    else:
        logger.info(f"Inferred variable ploidy across contigs: {ploidy_map}")
        # WC-FST is only supported for diploid data; warn up front if a non-diploid contig
        # will be encountered while WC-FST is requested.
        if "fst" in args.stats and args.fst_type.upper() == "WC":
            non_diploid = [c for c, p in ploidy_map.items() if p != 2]
            if non_diploid:
                logger.warning(
                    "Weir-Cockerham FST is not supported for non-diploid contigs. "
                    f"FST will be skipped for: {non_diploid}. "
                    "Use --fst_type hudson to compute FST on these contigs."
                )

    logger.info("All initial checks passed!")
    stats: List[PixyStat] = [PixyStat[stat.upper()] for stat in args.stats]
    tmp_path: Path = _generate_tmp_path(output_dir=output_folder)
    _check_tmp_path(tmp_path)

    return PixyArgs(
        stats=stats,
        vcf_path=Path(vcf_path),
        populations_df=populations_df,
        num_cores=args.n_cores,
        bypass_invariant_check=bypass_invariant_check,
        include_multiallelic_snps=include_multiallelic_snps,
        bed_df=bed_df,
        output_dir=Path(output_folder),
        output_prefix=output_prefix,
        chromosomes=chrom_list,
        window_size=args.window_size,
        interval_start=args.interval_start,
        interval_end=args.interval_end,
        sites_df=sites_df,
        chunk_size=args.chunk_size,
        fst_type=FSTEstimator[args.fst_type.upper()],
        fst_components=getattr(args, "fst_components", False),
        tajima_components=getattr(args, "tajima_components", False),
        temp_file=tmp_path,
        ploidy_map=ploidy_map,
    )


def _generate_tmp_path(output_dir: str) -> Path:
    """Generates a temporary file path to which `pixy` will write intermediate results."""
    return Path(output_dir) / f"pixy_tmpfile_{uuid.uuid4().hex}.tmp"


def _check_tmp_path(temp_file: Path) -> None:
    # check if temp file is writable
    with open(temp_file, "w"):
        pass  # file is created and then closed
    assert os.access(temp_file, os.W_OK), "temp file is not writable"
