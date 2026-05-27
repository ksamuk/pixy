#!/usr/bin/env python
# coding: utf-8
#
# Heavy runtime imports (numpy, pandas, multiprocessing, pixy.calc, pixy.core) are deferred to
# inside main(), after argparse has finished. The reason: `pixy --version`, `pixy --help`, and
# argparse error paths previously paid ~0.9 s and ~160 MB just to discover that argparse was
# going to bail out. With deferred imports those paths return in ~0.2 s / ~30 MB.
#
# `from __future__ import annotations` makes all annotations below evaluate lazily, so we can
# annotate parameters with `PixyArgs` / `pandas.DataFrame` / etc. without importing them.
from __future__ import annotations

import argparse
import logging
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import TYPE_CHECKING
from typing import List
from typing import Optional
from typing import cast

from pixy.enums import FSTEstimator
from pixy.enums import PixyStat

if TYPE_CHECKING:
    from multiprocessing.context import BaseContext

    from pixy.args_validation import PixyArgs

# main pixy function


def main() -> None:  # noqa: C901
    """Parse arguments and execute pixy."""
    # argument parsing via argparse

    # the ascii help image
    help_image = "█▀▀█ ░▀░ █░█ █░░█\n█░░█ ▀█▀ ▄▀▄ █▄▄█\n█▀▀▀ ▀▀▀ ▀░▀ ▄▄▄█\n"

    help_text = (
        "pixy: unbiased estimates of pi, dxy, fst, Watterson's Theta, and Tajima's D from "
        "VCFs with invariant sites"
    )
    version = "2.0.0"
    citation = (
        "Korunes, KL and K Samuk. "
        "pixy: Unbiased estimation of nucleotide diversity and divergence in the presence of "
        "missing data. "
        "Mol Ecol Resour. 2021 Jan 16. doi: 10.1111/1755-0998.13326."
    )
    # initialize logger
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [pixy] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logger = logging.getLogger(__name__)
    # initialize arguments
    parser = argparse.ArgumentParser(
        description=help_image + help_text + "\n" + version,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")
    additional = parser.add_argument_group("in addition, one of")
    optional = parser.add_argument_group("optional arguments")

    required.add_argument(
        "--stats",
        nargs="+",
        choices=["pi", "dxy", "fst", "watterson_theta", "tajima_d"],
        help=(
            "List of statistics to calculate from the VCF, separated by spaces.\n"
            'e.g. "--stats pi fst" will result in pi and fst calculations.'
        ),
        required=True,
    )
    required.add_argument(
        "--vcf",
        type=str,
        nargs="?",
        help="Path to the input VCF (bgzipped and tabix indexed).",
        required=True,
    )
    required.add_argument(
        "--populations",
        type=str,
        nargs="?",
        help=(
            "Path to a headerless tab separated populations file with columns "
            "[SampleID Population]."
        ),
        required=True,
    )

    # NB: --window_size and --bed_file are mutually exclusive — that constraint is enforced in
    # pixy/args_validation.py rather than at argparse level so that the user gets a pixy-specific
    # ValueError (and so that the existing test suite, which matches on that message, keeps
    # working). Moving the check into argparse's `add_mutually_exclusive_group` is on the
    # roadmap; it requires updating the corresponding regression tests to match SystemExit.
    additional.add_argument(
        "--window_size",
        type=int,
        nargs="?",
        help=(
            "Window size in base pairs over which to calculate stats.\n"
            "Automatically determines window coordinates/bounds (see additional options below)."
        ),
        required=False,
    )
    additional.add_argument(
        "--bed_file",
        type=str,
        nargs="?",
        help=(
            "Path to a headerless .BED file with columns [chrom chromStart chromEnd].\n"
            "Manually defines window bounds, which can be heterogeneous in size."
        ),
        required=False,
    )

    optional.add_argument(
        "--n_cores",
        type=int,
        nargs="?",
        default=1,
        help="Number of CPUs to utilize for parallel processing (default=1).",
        required=False,
    )
    optional.add_argument(
        "--output_folder",
        type=str,
        nargs="?",
        default="",
        help=(
            "Folder where output will be written, e.g. path/to/output_folder.\n"
            "Defaults to current working directory."
        ),
        required=False,
    )
    optional.add_argument(
        "--output_prefix",
        type=str,
        nargs="?",
        default="pixy",
        help=(
            "Optional prefix for output file(s), with no slashes.\n"
            'e.g. "--output_prefix output" will result in [output folder]/output_pi.txt. \n'
            "Defaults to 'pixy'."
        ),
        required=False,
    )
    optional.add_argument(
        "--chromosomes",
        type=str,
        nargs="?",
        default="all",
        help=(
            "A single-quoted, comma separated list of chromosomes where stats will be calculated.\n"
            "e.g. --chromosomes 'X,1,2' will restrict stats to chromosomes X, 1, and 2.\n"
            "Defaults to all chromosomes in the VCF."
        ),
        required=False,
    )
    optional.add_argument(
        "--interval_start",
        type=str,
        nargs="?",
        help=(
            "The start of an interval over which to calculate stats.\n"
            "Only valid when calculating over a single chromosome.\n"
            "Defaults to 1."
        ),
        required=False,
    )
    optional.add_argument(
        "--interval_end",
        type=str,
        nargs="?",
        help=(
            "The end of the interval over which to calculate stats.\n"
            "Only valid when calculating over a single chromosome.\n"
            "Defaults to the last position for a chromosome."
        ),
        required=False,
    )
    optional.add_argument(
        "--sites_file",
        type=str,
        nargs="?",
        help=(
            "Path to a headerless tab separated file with columns [CHROM POS].\n"
            "This defines sites where summary stats should be calculated.\n"
            "Can be combined with the --bed_file and --window_size options."
        ),
        required=False,
    )
    optional.add_argument(
        "--chunk_size",
        type=int,
        nargs="?",
        default=100000,
        help=(
            "Approximate number of sites to read from VCF at any given time (default=100000).\n"
            "Larger numbers reduce I/O operations at the cost of memory."
        ),
        required=False,
    )

    optional.add_argument(
        "--fst_type",
        choices=["wc", "hudson"],
        default="wc",
        help=(
            "FST estimator to use, one of either: \n"
            "'wc' (Weir and Cockerham 1984) or\n"
            "'hudson' (Hudson 1992, Bhatia et al. 2013) \n"
            "Defaults to 'wc'."
        ),
        required=False,
    )
    optional.add_argument(
        "--fst_components",
        action="store_true",
        default=False,
        help=(
            "Include FST estimator components in the FST output table. "
            "For --fst_type wc, output wc_fst_a, wc_fst_b, and wc_fst_c. "
            "For --fst_type hudson, output hudson_fst_num and hudson_fst_den. "
            "Defaults to False."
        ),
        required=False,
    )
    optional.add_argument(
        "--tajima_components",
        action="store_true",
        default=False,
        help=(
            "Include Tajima's D aggregation components in the Tajima's D output table. "
            "This adds tajima_d_s_counts, a comma-separated list of "
            "observed_alleles:segregating_sites pairs that can be summed across windows. "
            "Defaults to False."
        ),
        required=False,
    )
    (
        optional.add_argument(
            "--include_multiallelic_snps",
            action="store_true",
            default=False,
            help=(
                "Multiallelic SNPs within the VCF will be included "
                "during calculation.(default=False)."
            ),
            required=False,
        ),
    )

    (
        optional.add_argument(
            "--bypass_invariant_check",
            action="store_true",
            default=False,
            help=(
                "Allow computation of stats without invariant sites (default=False).\n"
                "Will result in wildly incorrect estimates most of the time.\n"
                "Use with extreme caution!"
            ),
            required=False,
        ),
    )

    optional.add_argument(
        "--version",
        action="version",
        version=help_image + "version " + version,
        help="Print the version of pixy in use.",
    )
    optional.add_argument(
        "--citation",
        action="version",
        version=citation,
        help="Print the citation for pixy.",
    )
    optional.add_argument("--silent", action="store_true", help="Suppress all console output.")
    optional.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)
    optional.add_argument("--keep_temp_file", action="store_true", help=argparse.SUPPRESS)

    # catch arguments from the command line
    # automatically uncommented when a release is built
    args = parser.parse_args()

    # if not running in debug mode, suppress traceback
    if not args.debug:
        sys.tracebacklimit = 0

    # if running in silent mode, suppress output
    if args.silent:
        sys.stdout = open(os.devnull, "w")

    # ------------------------------------------------------------------
    # Heavy library imports (deferred from module top — see file header)
    # ------------------------------------------------------------------
    # By this point argparse has already exited for --version / --citation / --help / bad-arg
    # paths, so anything below only runs for an actual analysis invocation. Each `import x`
    # below binds `x` as a local in main(); later uses in the function find them via the
    # usual local-then-global lookup.
    import multiprocessing as mp

    # `BaseContext` is referenced at runtime by `cast(BaseContext, ...)` further down (not just
    # in annotations), so it has to be imported eagerly here — the TYPE_CHECKING-only import
    # at the top of the module covers static analysis but not the runtime cast.
    from multiprocessing.context import BaseContext

    import pixy.args_validation
    import pixy.calc  # noqa: F401  — workers reach it through pixy.core
    import pixy.core

    # validate arguments with the check_and_validate_args fuction
    # returns parsed populaion, chromosome, and sample info
    logger.info(f"pixy {version}")
    logger.info("See documentation at https://pixy.readthedocs.io/en/latest/")

    pixy_args: PixyArgs = pixy.args_validation.check_and_validate_args(args)
    # propagate per-contig ploidy map onto the raw args namespace so it is available to
    # worker functions (which currently receive `args`, not `pixy_args`).
    args.ploidy_map = pixy_args.ploidy_map
    popindices = {
        name: pixy_args.populations.indices_for(str(name)) for name in pixy_args.pop_names
    }
    chrom_list = pixy_args.chromosomes

    # For FST-only runs, scale up the chunk size to reduce VCF I/O overhead.
    # FST needs only variant sites; invariant sites are still read by allel.read_vcf (we cannot
    # skip them at parse time) but are filtered out immediately afterward. Larger chunks mean
    # fewer total allel.read_vcf calls, which is the dominant runtime cost for FST. The
    # multiplier trades proportionally more memory per chunk for far fewer I/O round-trips.
    # Users who are memory-constrained can override with an explicit --chunk_size.
    _stats_needing_invariants: set = {
        PixyStat.PI,
        PixyStat.DXY,
        PixyStat.TAJIMA_D,
        PixyStat.WATTERSON_THETA,
    }
    _needs_invariants: bool = bool(_stats_needing_invariants.intersection(pixy_args.stats))
    _fst_chunk_multiplier: int = 10
    effective_chunk_size: int = (
        pixy_args.chunk_size if _needs_invariants else pixy_args.chunk_size * _fst_chunk_multiplier
    )

    logger.info(
        f"Preparing for calculation of summary statistics: {', '.join(map(str, args.stats))}"
    )
    if not _needs_invariants:
        logger.info(
            f"FST-only run: using adaptive chunk size of {effective_chunk_size:,} bp "
            f"({_fst_chunk_multiplier}× the --chunk_size of {pixy_args.chunk_size:,} bp). "
            "Pass --chunk_size to override if memory is limited."
        )

    fst_cite: str
    if PixyStat.FST in pixy_args.stats:
        if pixy_args.fst_type is FSTEstimator.WC:
            fst_cite = "Weir and Cockerham (1984)"
        elif pixy_args.fst_type is FSTEstimator.HUDSON:
            fst_cite = "Hudson (1992)"
        logger.info(f"Using {fst_cite}'s estimator of FST.")

    logger.info(
        f"Data set contains {len(pixy_args.pop_names)} populations, "
        f"{len(chrom_list)} chromosome(s), "
        f"and {len(pixy_args.pop_ids)} sample(s)"
    )

    if pixy_args.window_size is not None:
        logger.info(f"Window size: {pixy_args.window_size} bp")

    if args.bed_file is not None:
        logger.info(f"Windows sourced from: {args.bed_file}")

    if args.sites_file is not None:
        logger.info(f"Calculations restricted to sites in {args.sites_file}")

    # time the calculations
    start_time = time.time()
    logger.info("Started calculations!")
    logger.info(f"Using {pixy_args.num_cores} out of {mp.cpu_count()} available CPU cores")
    # if in mc mode, set up multiprocessing
    if pixy_args.num_cores > 1:
        # Use `forkserver` on Linux and `spawn` elsewhere (macOS, Windows).
        # Plain `fork` would be marginally faster on Linux but is deprecated in
        # Python 3.14+ when the parent process has any threads — and numpy/scipy
        # BLAS initialization can start threads on import, which would trigger
        # the deprecation. `forkserver` keeps fork-style copy-on-write semantics
        # while running the actual fork from a dedicated single-threaded helper.
        ctx: BaseContext

        # stdlib `multiprocessing.get_context` returns `BaseContext | None`; cast pins the
        # narrower BaseContext for the type checker.
        if sys.platform == "linux":
            ctx = cast(BaseContext, mp.get_context("forkserver"))
        else:
            ctx = cast(BaseContext, mp.get_context("spawn"))

        # set up the multiprocessing manager, queue, and process pool. Types intentionally
        # left to inference — stdlib `multiprocessing` exposes `Pool`/`Queue` as factory
        # functions (not classes), and `manager.Queue()` returns a proxy whose stub type
        # doesn't match `multiprocessing.queues.Queue`. Annotating any of them invites mypy
        # noise without buying real safety; the call sites are short and locally obvious.
        manager = ctx.Manager()
        q = manager.Queue()
        pool = ctx.Pool(int(args.n_cores))

        # launch the listener for collecting output asynchronously. The listener lives in
        # `pixy.core` (not here) because worker processes started under forkserver/spawn
        # cannot resolve attributes of the `__main__` module — pickling a function from
        # `__main__` only works if the workers inherit the parent's `__main__` namespace,
        # which `fork` does but `forkserver`/`spawn` do not.
        watcher = pool.apply_async(  # noqa: F841
            pixy.core.temp_file_listener,
            args=(
                q,
                str(pixy_args.temp_file),
            ),
        )

    # begin processing each chromosome

    for chromosome in pixy_args.chromosomes:
        logger.info(f"Processing chromosome/contig {chromosome}")

        # if not using a bed file, build windows manually
        if pixy_args.bed is None:
            if pixy_args.window_size is not None:
                window_size: int = int(pixy_args.window_size)
            if pixy_args.has_interval:  # if an interval is specified, assign it
                # mypy cannot infer that by this point,
                # `interval_start` and `interval_end` are not `None`; we assert they are not `None`
                assert pixy_args.interval_start is not None
                assert pixy_args.interval_end is not None
                interval_start: int = int(pixy_args.interval_start)
                interval_end: int = int(pixy_args.interval_end)
            # otherwise, get the interval from the VCF's POS column
            else:
                if pixy_args.sites is None:
                    # Run `tabix VCF CHROM` (no shell, no piping through cut/tail) and read the
                    # POS column of the final line. Avoids shell-injection risk on the VCF and
                    # chromosome strings and removes a dependency on cut/tail being on PATH.
                    tabix_out = subprocess.check_output(["tabix", args.vcf, chromosome]).decode(
                        "utf-8"
                    )
                    last_pos: Optional[str] = None
                    for line in tabix_out.splitlines():
                        if line and not line.startswith("#"):
                            cols = line.split("\t", 2)
                            if len(cols) >= 2:
                                last_pos = cols[1]
                    if last_pos is None:
                        raise ValueError(
                            f"No records found in VCF for chromosome {chromosome!r}; cannot infer "
                            "interval end."
                        )
                    interval_start = 1
                    interval_end = int(last_pos)
                else:
                    sites_pre_list = pixy_args.sites.positions_for(chromosome)
                    interval_start = min(sites_pre_list)
                    interval_end = max(sites_pre_list)

            # final check if intervals are valid
            if interval_start > interval_end:
                raise ValueError(
                    f"The specified interval start {interval_start} exceeds the specified "
                    f"interval end {interval_end}"
                )

            targ_region = chromosome + ":" + str(interval_start) + "-" + str(interval_end)

            logger.info(f"Calculating statistics for region {targ_region}")

            # Determine list of windows over which to compute stats
            # in the case were window size = 1, AND there is a sites file, use the sites file as the
            # 'windows'
            if (
                pixy_args.sites is not None and window_size == 1
            ):  # TODO: dig into why `window_size` might be unbound
                # reference https://github.com/fulcrumgenomics/pixy-dev/issues/70
                window_list = [list(a) for a in zip(sites_pre_list, sites_pre_list, strict=True)]
            else:
                # if the interval is smaller than one window, make a list of length 1
                if (interval_end - interval_start) <= window_size:
                    window_pos_1_list = [interval_start]
                    window_pos_2_list = [interval_start + window_size - 1]
                else:
                    # if the interval_end is not a perfect multiple of the window size
                    # bump the interval_end up to the nearest multiple of window size
                    if not (interval_end % window_size == 0):
                        interval_end = interval_end + (window_size - (interval_end % window_size))

                    # create the start and stops for each window
                    window_pos_1_list = [*range(interval_start, int(interval_end), window_size)]
                    window_pos_2_list = [
                        *range(
                            interval_start + (window_size - 1),
                            int(interval_end) + window_size,
                            window_size,
                        )
                    ]

                # `strict=False` is intentional: when `window_size == 1`, `window_pos_2_list`
                # is one element longer than `window_pos_1_list` by construction (the trailing
                # entry is meaningful for larger window sizes but spurious here), and the loop
                # relies on `zip()` truncating to the shorter sequence.
                window_list = [
                    list(a) for a in zip(window_pos_1_list, window_pos_2_list, strict=False)
                ]

            # Set aggregate to true if
            # 1) the window size is larger than the chunk size OR
            # 2) the window size wasn't specified, but the chrom is longer than the cutoff
            if (window_size > effective_chunk_size) or (
                (pixy_args.window_size is None)
                and ((interval_end - interval_start) > effective_chunk_size)
            ):
                aggregate = True
            else:
                aggregate = False

        # if using a bed file, subset the bed file for the current chromosome
        else:
            aggregate = False
            window_list = pixy_args.bed.intervals_for(chromosome)

        if len(window_list) == 0:
            raise Exception(
                "ERROR: Window creation failed. Ensure that the POS column in the VCF is "
                "valid or change --window_size."
            )

        # if aggregating, break down large windows into smaller windows
        if aggregate:
            window_list = pixy.core.assign_subwindows_to_windows(window_list, effective_chunk_size)

        # using chunk_size, assign  windows to chunks
        window_list = pixy.core.assign_windows_to_chunks(window_list, effective_chunk_size)

        # if using a sites file, assign sites to chunks, as with windows above
        if pixy_args.sites is not None:
            sites_pre_list = pixy_args.sites.positions_for(chromosome)
            sites_list = pixy.core.assign_sites_to_chunks(sites_pre_list, effective_chunk_size)
        else:
            sites_list = None
        # obtain the list of chunks from the window list
        chunk_list = [i[2] for i in window_list]
        chunk_list = list(set(chunk_list))

        # if running in mc mode, send the summary stats function to the jobs pool
        if pixy_args.num_cores > 1:
            # the list of jobs to be launched
            jobs = []

            for chunk in chunk_list:
                # create a subset of the window list specific to this chunk
                window_list_chunk = [x for x in window_list if x[2] == chunk]

                # and for the site list (if it exists)
                sites_list_chunk: Optional[List[int]]
                if sites_list is None:
                    sites_list_chunk = None
                else:
                    chunk_sites_lists: List[List[int]] = [x for x in sites_list if x[1] == chunk]
                    sites_list_chunk = [x[0] for x in chunk_sites_lists]

                # determine the bounds of the chunk
                chunk_pos_1 = min(window_list_chunk, key=lambda x: x[1])[0]
                chunk_pos_2 = max(window_list_chunk, key=lambda x: x[1])[1]

                # launch a summary stats job for this chunk
                job = pool.apply_async(
                    pixy.core.compute_summary_stats,
                    args=(
                        args,
                        pixy_args.pop_names,
                        popindices,
                        pixy_args.temp_file,
                        chromosome,
                        chunk_pos_1,
                        chunk_pos_2,
                        window_list_chunk,
                        q,
                        sites_list_chunk,
                        aggregate,
                        args.window_size,
                    ),
                )
                jobs.append(job)

            # launch all the jobs onto the pool
            for job in jobs:
                job.get()

        # if running in single core mode, loop over the function manually
        elif pixy_args.num_cores == 1:
            for chunk in chunk_list:
                # create a subset of the window list specific to this chunk
                window_list_chunk = [x for x in window_list if x[2] == chunk]

                # and for the site list (if it exists)
                if pixy_args.sites is not None and sites_list is not None:
                    chunk_sites_lists = [x for x in sites_list if x[1] == chunk]
                    sites_list_chunk = [x[0] for x in chunk_sites_lists]
                else:
                    sites_list_chunk = None

                # determine the bounds of the chunk
                chunk_pos_1 = min(window_list_chunk, key=lambda x: x[1])[0]
                chunk_pos_2 = max(window_list_chunk, key=lambda x: x[1])[1]

                # don't use the queue (q) when running in single core mode; rebind to the
                # sentinel string `compute_summary_stats` checks for. The type changes here
                # (Queue proxy -> str) so the local is intentionally untyped.
                q = "NULL"  # type: ignore[assignment]

                # compute summary stats for all windows in the chunk window list
                pixy.core.compute_summary_stats(
                    args,
                    pixy_args.pop_names,
                    popindices,
                    pixy_args.temp_file,
                    chromosome,
                    chunk_pos_1,
                    chunk_pos_2,
                    window_list_chunk,
                    q,
                    sites_list_chunk,
                    aggregate,
                    args.window_size,
                )

    # clean up any remaining jobs and stop the listener
    if pixy_args.num_cores > 1:
        q.put("kill")
        pool.close()
        pool.join()

    # ------------------------------------------------------------------
    # Split and aggregate temp file to per-stat output files (streaming).
    # ------------------------------------------------------------------
    # Previously this loaded the entire tmp file with `pandas.read_csv` then did
    # `groupby + to_csv` per stat. Now we stream the tmp file once into a
    # `dict[(stat, chrom)] -> List[TempRow]` and let `pixy.agg.write_stat_file`
    # format/aggregate without pandas. Output is byte-identical to the prior path.
    import pixy.agg

    rows_by_stat_chrom: dict = {}
    successful_stat_names: set = set()
    found_any = False
    for trow in pixy.agg.iter_temp_rows(pixy_args.temp_file):
        found_any = True
        rows_by_stat_chrom.setdefault((trow.stat, trow.chrom), []).append(trow)
        successful_stat_names.add(trow.stat)
    if not found_any:
        raise Exception(
            "ERROR: pixy failed to write any output. Confirm that your bed/sites files and "
            "intervals refer to existing chromosomes and positions in the VCF."
        )

    successful_stats = [PixyStat(name) for name in successful_stat_names]

    if set(pixy_args.stats) != set(successful_stats):
        missing_stats = list(set(pixy_args.stats) - set(successful_stats))
        logger.warning(
            "WARNING: pixy failed to find any valid genotype data to calculate the "
            f"following summary statistics {', '.join([str(stat) for stat in missing_stats])}."
            " No output file will be created for these statistics."
        )

    # enforce chromosome IDs as strings
    chrom_list = list(map(str, pixy_args.chromosomes))

    # `window_size` is only bound when the per-chromosome loop above ran with no BED file
    # (BED mode never set it). For BED runs `aggregate` is always False, so the value is
    # never used by the writer — but it still has to be passed, so resolve it now.
    output_window_size: int = int(pixy_args.window_size) if pixy_args.window_size is not None else 0

    def _rows_for(stat_name: str) -> dict:
        """Return dict[chromosome] -> List[TempRow] for one statistic."""
        return {c: rows_by_stat_chrom.get((stat_name, c), []) for c in chrom_list}

    if PixyStat.PI in successful_stats:
        pixy.agg.write_stat_file(
            out_path=Path(pixy_args.output_prefix + "_pi.txt"),
            stat="pi",
            rows_by_chrom=_rows_for("pi"),
            chrom_list=chrom_list,
            aggregate=aggregate,
            window_size=output_window_size,
            fst_type=pixy_args.fst_type.value,
        )

    if PixyStat.DXY in successful_stats:
        pixy.agg.write_stat_file(
            out_path=Path(pixy_args.output_prefix + "_dxy.txt"),
            stat="dxy",
            rows_by_chrom=_rows_for("dxy"),
            chrom_list=chrom_list,
            aggregate=aggregate,
            window_size=output_window_size,
            fst_type=pixy_args.fst_type.value,
        )

    if PixyStat.FST in successful_stats:
        chroms_with_no_data = pixy.agg.write_stat_file(
            out_path=Path(pixy_args.output_prefix + "_fst.txt"),
            stat="fst",
            rows_by_chrom=_rows_for("fst"),
            chrom_list=chrom_list,
            aggregate=aggregate,
            window_size=output_window_size,
            fst_type=pixy_args.fst_type.value,
            fst_components=pixy_args.fst_components,
        )
        if chroms_with_no_data:
            logger.info(
                "NOTE: The following chromosomes/scaffolds did not have sufficient data "
                f"to estimate FST: {', '.join(chroms_with_no_data)}"
            )

    if PixyStat.WATTERSON_THETA in successful_stats:
        pixy.agg.write_stat_file(
            out_path=Path(f"{pixy_args.output_prefix}_{PixyStat.WATTERSON_THETA.value}.txt"),
            stat="watterson_theta",
            rows_by_chrom=_rows_for("watterson_theta"),
            chrom_list=chrom_list,
            aggregate=aggregate,
            window_size=output_window_size,
            fst_type=pixy_args.fst_type.value,
        )

    if PixyStat.TAJIMA_D in successful_stats:
        pixy.agg.write_stat_file(
            out_path=Path(f"{pixy_args.output_prefix}_{PixyStat.TAJIMA_D.value}.txt"),
            stat="tajima_d",
            rows_by_chrom=_rows_for("tajima_d"),
            chrom_list=chrom_list,
            aggregate=aggregate,
            window_size=output_window_size,
            fst_type=pixy_args.fst_type.value,
            tajima_components=pixy_args.tajima_components,
        )
    # remove the temp file(s)
    if args.keep_temp_file is not True:
        os.remove(pixy_args.temp_file)

    # confirm output was generated successfully
    outfolder_files = [
        f
        for f in os.listdir(pixy_args.output_dir)
        if os.path.isfile(os.path.join(pixy_args.output_dir, f))
    ]

    r = re.compile(".*_dxy.*|.*_pi.*|.*_fst.*|.*_watterson_theta.*|.*_tajima_d.*|")
    output_files = list(filter(r.match, outfolder_files))

    r = re.compile("pixy_tmpfile.*")
    leftover_tmp_files = list(filter(r.match, outfolder_files))

    if len(output_files) == 0:
        logger.warning(
            "WARNING: pixy failed to write any output files. Your VCF may not contain "
            "valid genotype data, or it was removed via filtering using the specified sites/bed "
            "file (if any)."
        )

    # print completion message
    logger.info("All calculations complete!")
    total_time = time.time() - start_time
    logger.info("Time elapsed: " + time.strftime("%H:%M:%S", time.gmtime(total_time)))
    logger.info(f"Output files written to {pixy_args.output_dir}")

    if len(leftover_tmp_files) > 0:
        logger.info(
            f"NOTE: There are pixy temp files in {pixy_args.output_dir}."
            "If these are not actively being used (e.g. by another running pixy process), "
            "they can be safely deleted."
        )

    logger.info(f"If you use pixy in your research, please cite the following paper: {citation}")

    # restore output
    if args.silent:
        sys.stdout = sys.__stdout__


if __name__ == "__main__":
    main()
