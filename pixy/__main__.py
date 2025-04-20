#!/usr/bin/env python
# coding: utf-8
import argparse
import logging
import os
import re
import subprocess
import sys
import time
from multiprocessing.context import BaseContext
from multiprocessing.managers import SyncManager
from multiprocessing.pool import ApplyResult
from typing import List
from typing import Optional
from typing import cast

import multiprocess as mp
import numpy as np
import pandas
from multiprocess import Pool
from multiprocess import Queue

import pixy.calc
import pixy.core
from pixy.args_validation import PixyArgs
from pixy.enums import FSTEstimator
from pixy.enums import PixyStat

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
    version = "2.0.0.beta4"
    citation = (
        "Korunes, KL and K Samuk. "
        "pixy: Unbiased estimation of nucleotide diversity and divergence in the presence of "
        "missing data. "
        "Mol Ecol Resour. 2021 Jan 16. doi: 10.1111/1755-0998.13326."
    )
    # initialize logger
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(name)s:%(funcName)s:%(lineno)s [%(levelname)s]: %(message)s",
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

    # validate arguments with the check_and_validate_args fuction
    # returns parsed populaion, chromosome, and sample info
    logger.info(f"[pixy] pixy {version}")
    logger.info("[pixy] See documentation at https://pixy.readthedocs.io/en/latest/")

    pixy_args: PixyArgs = pixy.args_validation.check_and_validate_args(args)
    popindices = {}
    for name in pixy_args.pop_names:
        popindices[name] = pixy_args.populations_df[
            pixy_args.populations_df.Population == name
        ].callset_index.values
    chrom_list = pixy_args.chromosomes

    logger.info(
        f"[pixy] Preparing for calculation of summary statistics: {', '.join(map(str, args.stats))}"
    )

    fst_cite: str
    if PixyStat.FST in pixy_args.stats:
        if pixy_args.fst_type is FSTEstimator.WC:
            fst_cite = "Weir and Cockerham (1984)"
        elif pixy_args.fst_type is FSTEstimator.HUDSON:
            fst_cite = "Hudson (1992)"
        logger.info(f"[pixy] Using {fst_cite}'s estimator of FST.")

    logger.info(
        f"[pixy] Data set contains {len(pixy_args.pop_names)} populations, "
        f"{len(chrom_list)} chromosome(s), "
        f"and {len(pixy_args.pop_ids)} sample(s)"
    )

    if pixy_args.window_size is not None:
        logger.info(f"[pixy] Window size: {pixy_args.window_size} bp")

    if args.bed_file is not None:
        logger.info(f"[pixy] Windows sourced from: {args.bed_file}")

    if args.sites_file is not None:
        logger.info(f"[pixy] Calculations restricted to sites in {args.sites_file}")

    # time the calculations
    start_time = time.time()
    logger.info(
        f"Started calculations at \
        {time.strftime('%H:%M:%S on %Y-%m-%d', time.localtime(start_time))}"
    )
    logger.info(f"[pixy] Using {pixy_args.num_cores} out of {mp.cpu_count()} available CPU cores")
    # if in mc mode, set up multiprocessing
    if pixy_args.num_cores > 1:
        # use forking context on linux, and spawn otherwise (macOS)
        ctx: BaseContext

        # `cast` is necessary because `multiprocess` is an untyped module, so mypy can only infer
        # that `get_context()` returns `Any`
        if sys.platform == "linux":
            ctx = cast(BaseContext, mp.get_context("fork"))
        else:
            ctx = cast(BaseContext, mp.get_context("spawn"))

        # set up the multiprocessing manager, queue, and process pool
        manager: SyncManager = ctx.Manager()
        q: Queue = manager.Queue()
        pool: Pool = ctx.Pool(int(args.n_cores))

        # a listener function for writing a temp file
        # used to write output in multicore mode
        def listener(q: Queue, temp_file: str) -> None:
            """
            Writes to a given `temp_file` in multicore mode.

            Args:
                q: the `Queue` from which to retrieve data
                temp_file: the file handle to which the data will be written
            """
            with open(temp_file, "a") as f:
                while 1:
                    m = q.get()
                    if m == "kill":  # we are done
                        break
                    f.write(str(m) + "\n")
                    f.flush()  # immediately write data, do not buffer

        # launch the watcher function for collecting output asynchronously
        watcher: ApplyResult = pool.apply_async(  # noqa: F841
            listener,
            args=(
                q,
                str(pixy_args.temp_file),
            ),
        )

    # begin processing each chromosome

    for chromosome in pixy_args.chromosomes:
        logger.info(f"[pixy] Processing chromosome/contig {chromosome}")

        # if not using a bed file, build windows manually
        if pixy_args.bed_df is None:
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
                if pixy_args.sites_df is None:
                    chrom_max: List[str] = (
                        subprocess.check_output(
                            "tabix " + args.vcf + " " + chromosome + " | cut -f 2 | tail -n 1",
                            shell=True,
                        )
                        .decode("utf-8")
                        .split()
                    )
                    interval_start = 1
                    interval_end = int(chrom_max[0])
                else:
                    sites_pre_list = pixy_args.sites_df[pixy_args.sites_df["CHROM"] == chromosome]
                    sites_pre_list = sorted(sites_pre_list["POS"].tolist())
                    interval_start = min(sites_pre_list)
                    interval_end = max(sites_pre_list)

            # final check if intervals are valid
            if interval_start > interval_end:
                raise ValueError(
                    f"The specified interval start {interval_start} exceeds the specified "
                    f"interval end {interval_end}"
                )

            targ_region = chromosome + ":" + str(interval_start) + "-" + str(interval_end)

            logger.info(f"[pixy] Calculating statistics for region {targ_region}")

            # Determine list of windows over which to compute stats
            # in the case were window size = 1, AND there is a sites file, use the sites file as the
            # 'windows'
            if (
                pixy_args.sites_df is not None and window_size == 1
            ):  # TODO: dig into why `window_size` might be unbound
                # reference https://github.com/fulcrumgenomics/pixy-dev/issues/70
                window_list = [list(a) for a in zip(sites_pre_list, sites_pre_list)]
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

                window_list = [list(a) for a in zip(window_pos_1_list, window_pos_2_list)]

            # Set aggregate to true if
            # 1) the window size is larger than the chunk size OR
            # 2) the window size wasn't specified, but the chrom is longer than the cutoff
            if (window_size > pixy_args.chunk_size) or (
                (pixy_args.window_size is None)
                and ((interval_end - interval_start) > pixy_args.chunk_size)
            ):
                aggregate = True
            else:
                aggregate = False

        # if using a bed file, subset the bed file for the current chromosome
        else:
            aggregate = False
            bed_df_chrom = pixy_args.bed_df.loc[pixy_args.bed_df["chrom"] == chromosome]
            window_list = [
                list(a) for a in zip(bed_df_chrom["chromStart"], bed_df_chrom["chromEnd"])
            ]

        if len(window_list) == 0:
            raise Exception(
                "[pixy] ERROR: Window creation failed. Ensure that the POS column in the VCF is "
                "valid or change --window_size."
            )

        # if aggregating, break down large windows into smaller windows
        if aggregate:
            window_list = pixy.core.assign_subwindows_to_windows(window_list, pixy_args.chunk_size)

        # using chunk_size, assign  windows to chunks
        window_list = pixy.core.assign_windows_to_chunks(window_list, pixy_args.chunk_size)

        # if using a sites file, assign sites to chunks, as with windows above
        if pixy_args.sites_df is not None:
            sites_pre_list = pixy_args.sites_df[pixy_args.sites_df["CHROM"] == chromosome]
            sites_pre_list = sites_pre_list["POS"].tolist()
            sites_list = pixy.core.assign_sites_to_chunks(sites_pre_list, pixy_args.chunk_size)
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
                if pixy_args.sites_df is not None and sites_list is not None:
                    chunk_sites_lists = [x for x in sites_list if x[1] == chunk]
                    sites_list_chunk = [x[0] for x in chunk_sites_lists]
                else:
                    sites_list_chunk = None

                # determine the bounds of the chunk
                chunk_pos_1 = min(window_list_chunk, key=lambda x: x[1])[0]
                chunk_pos_2 = max(window_list_chunk, key=lambda x: x[1])[1]

                # don't use the queue (q) when running in single core mode
                q = "NULL"

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

    # split and aggregate temp file to individual files

    # check if there is any output to process
    # halt execution if not
    try:
        outpanel: pandas.DataFrame = pandas.read_csv(pixy_args.temp_file, sep="\t", header=None)
    except pandas.errors.EmptyDataError as e:
        raise Exception(
            "[pixy] ERROR: pixy failed to write any output. Confirm that your bed/sites files and "
            "intervals refer to existing chromosomes and positions in the VCF."
        ) from e

    # check if particular stats failed to generate output
    # if not all requested stats were generated, produce a warning
    # and then remove the failed stats from the args list
    # TODO: add a unit-test for coverage here after closing
    # https://github.com/fulcrumgenomics/pixy-dev/issues/12
    successful_stats = [PixyStat(stat) for stat in np.unique(outpanel[0])]

    if set(pixy_args.stats) != set(successful_stats):
        missing_stats = list(set(pixy_args.stats) - set(successful_stats))
        logger.warning(
            "[pixy] WARNING: pixy failed to find any valid genotype data to calculate the "
            f"following summary statistics {', '.join([str(stat) for stat in missing_stats])}."
            " No output file will be created for these statistics."
        )

    outpanel[3] = outpanel[3].astype(str)  # force chromosome IDs to string
    outgrouped = outpanel.groupby([0, 3])  # groupby statistic, chromosome

    # enforce chromosome IDs as strings
    chrom_list = list(map(str, pixy_args.chromosomes))
    stat: str
    if PixyStat.PI in successful_stats:
        stat = "pi"
        pi_file: str = pixy_args.output_prefix + "_pi.txt"

        if os.path.exists(pi_file):
            os.remove(pi_file)

        outfile = open(pi_file, "a")
        outfile.write(
            "pop"
            + "\t"
            + "chromosome"
            + "\t"
            + "window_pos_1"
            + "\t"
            + "window_pos_2"
            + "\t"
            + "avg_pi"
            + "\t"
            + "no_sites"
            + "\t"
            + "count_diffs"
            + "\t"
            + "count_comparisons"
            + "\t"
            + "count_missing"
            + "\n"
        )

        if aggregate:  # put winsizes back together for each population to make final_window_size
            for chromosome in chrom_list:
                outpi = outgrouped.get_group(("pi", chromosome)).reset_index(
                    drop=True
                )  # get this statistic, this chrom only
                outpi.drop(
                    [0, 2], axis=1, inplace=True
                )  # get rid of "pi" and placeholder (NA) columns
                outsorted = pixy.core.aggregate_output(
                    outpi, stat, chromosome, window_size, pixy_args.fst_type.value
                )
                outsorted.to_csv(
                    outfile, sep="\t", mode="a", header=False, index=False, na_rep="NA"
                )  # write

        else:
            for chromosome in chrom_list:
                outpi = outgrouped.get_group(("pi", chromosome)).reset_index(
                    drop=True
                )  # get this statistic, this chrom only
                outpi.drop(
                    [0, 2], axis=1, inplace=True
                )  # get rid of "pi" and placeholder (NA) columns
                outsorted = outpi.sort_values([4])  # sort by position
                # make sure sites, comparisons, missing get written as integers
                cols = [7, 8, 9, 10]
                outsorted[cols] = outsorted[cols].astype("Int64")
                outsorted.to_csv(
                    outfile, sep="\t", mode="a", header=False, index=False, na_rep="NA"
                )  # write

        outfile.close()

    if PixyStat.DXY in successful_stats:
        stat = "dxy"
        dxy_file: str = pixy_args.output_prefix + "_dxy.txt"

        if os.path.exists(dxy_file):
            os.remove(dxy_file)

        outfile = open(dxy_file, "a")
        outfile.write(
            "pop1"
            + "\t"
            + "pop2"
            + "\t"
            + "chromosome"
            + "\t"
            + "window_pos_1"
            + "\t"
            + "window_pos_2"
            + "\t"
            + "avg_dxy"
            + "\t"
            + "no_sites"
            + "\t"
            + "count_diffs"
            + "\t"
            + "count_comparisons"
            + "\t"
            + "count_missing"
            + "\n"
        )

        if aggregate:  # put winsizes back together for each population to make final_window_size
            for chromosome in chrom_list:
                outdxy = outgrouped.get_group(("dxy", chromosome)).reset_index(
                    drop=True
                )  # get this statistic, this chrom only
                outdxy.drop([0], axis=1, inplace=True)  # get rid of "dxy"
                outsorted = pixy.core.aggregate_output(
                    outdxy, stat, chromosome, window_size, pixy_args.fst_type.value
                )
                outsorted.to_csv(
                    outfile, sep="\t", mode="a", header=False, index=False, na_rep="NA"
                )  # write

        else:
            for chromosome in chrom_list:
                outdxy = outgrouped.get_group(("dxy", chromosome)).reset_index(
                    drop=True
                )  # get this statistic, this chrom only
                outdxy.drop([0], axis=1, inplace=True)  # get rid of "dxy"
                outsorted = outdxy.sort_values([4])  # sort by position
                # make sure sites, comparisons, missing get written as integers
                cols = [7, 8, 9, 10]
                outsorted[cols] = outsorted[cols].astype("Int64")
                outsorted.to_csv(
                    outfile, sep="\t", mode="a", header=False, index=False, na_rep="NA"
                )  # write

        outfile.close()

    if PixyStat.FST in successful_stats:
        stat = "fst"
        fst_file: str = pixy_args.output_prefix + "_fst.txt"

        if os.path.exists(fst_file):
            os.remove(fst_file)

        outfile = open(fst_file, "a")
        outfile.write(
            "pop1"
            + "\t"
            + "pop2"
            + "\t"
            + "chromosome"
            + "\t"
            + "window_pos_1"
            + "\t"
            + "window_pos_2"
            + "\t"
            + "avg_"
            + args.fst_type
            + "_fst"
            + "\t"
            + "no_snps"
            + "\n"
        )

        # keep track of chrosomes with no fst data
        chroms_with_no_data = []

        if aggregate:  # put winsizes back together for each population to make final_window_size
            for chromosome in chrom_list:
                # logic to accommodate cases where pi/dxy have stats for a chromosome, but fst does
                # not
                chromosome_has_data = True

                # if there are no valid fst estimates, set chromosome_has_data = False
                try:
                    outfst = outgrouped.get_group(("fst", chromosome)).reset_index(
                        drop=True
                    )  # get this statistic, this chrom only
                except KeyError:
                    chroms_with_no_data.append(chromosome)
                    chromosome_has_data = False

                    pass

                if chromosome_has_data:
                    outfst.drop([0], axis=1, inplace=True)  # get rid of "fst"
                    outsorted = pixy.core.aggregate_output(
                        outfst, stat, chromosome, window_size, pixy_args.fst_type.value
                    )
                    outsorted = outsorted.iloc[:, :-3]  # drop components (for now)
                    outsorted.to_csv(
                        outfile,
                        sep="\t",
                        mode="a",
                        header=False,
                        index=False,
                        na_rep="NA",
                    )  # write

        else:
            for chromosome in chrom_list:
                # logic to accommodate cases where pi/dxy have stats for a chromosome, but fst does
                # not
                chromosome_has_data = True

                # if there are no valid fst estimates, set chromosome_has_data = False
                try:
                    outfst = outgrouped.get_group(("fst", chromosome)).reset_index(
                        drop=True
                    )  # get this statistic, this chrom only
                except KeyError:
                    chroms_with_no_data.append(chromosome)
                    chromosome_has_data = False
                    pass

                if chromosome_has_data:
                    outfst.drop([0], axis=1, inplace=True)  # get rid of "fst"
                    outsorted = outfst.sort_values([4])  # sort by position
                    # make sure sites (but not components like pi/dxy)
                    cols = [7]
                    outsorted[cols] = outsorted[cols].astype("Int64")
                    outsorted = outsorted.iloc[:, :-3]  # drop components (for now)
                    outsorted.to_csv(
                        outfile,
                        sep="\t",
                        mode="a",
                        header=False,
                        index=False,
                        na_rep="NA",
                    )

        outfile.close()

        if len(chroms_with_no_data) >= 1:
            logger.info(
                "[pixy] NOTE: The following chromosomes/scaffolds did not have sufficient data "
                f"to estimate FST: {', '.join(chroms_with_no_data)}"
            )
    if PixyStat.WATTERSON_THETA in successful_stats:
        watterson_theta_file = f"{pixy_args.output_prefix}_{PixyStat.WATTERSON_THETA.value}.txt"
        if os.path.exists(watterson_theta_file):
            os.remove(watterson_theta_file)
        outfile = open(watterson_theta_file, "a")
        outfile.write(
            "pop"
            + "\t"
            + "chromosome"
            + "\t"
            + "window_pos_1"
            + "\t"
            + "window_pos_2"
            + "\t"
            + "avg_watterson_theta"
            + "\t"
            + "no_sites"
            + "\t"
            + "raw_watterson_theta"
            + "\t"
            + "no_var_sites"
            + "\t"
            + "weighted_no_sites"
            + "\n"
        )

        if aggregate:  # put winsizes back together for each population to make final_window_size
            for chromosome in chrom_list:
                outwatterson_theta = outgrouped.get_group((
                    "watterson_theta",
                    chromosome,
                )).reset_index(drop=True)  # get this statistic, this chrom only
                outwatterson_theta.drop(
                    [0, 2], axis=1, inplace=True
                )  # get rid of "watterson_theta" and placeholder (NA) columns
                outsorted = pixy.core.aggregate_output(
                    outwatterson_theta,
                    PixyStat.WATTERSON_THETA.value,
                    chromosome,
                    window_size,
                    args.fst_type,
                )
                outsorted.to_csv(
                    outfile, sep="\t", mode="a", header=False, index=False, na_rep="NA"
                )  # write

        else:
            for chromosome in chrom_list:
                outwatterson_theta = outgrouped.get_group((
                    "watterson_theta",
                    chromosome,
                )).reset_index(drop=True)  # get this statistic, this chrom only
                outwatterson_theta.drop(
                    [0, 2], axis=1, inplace=True
                )  # get rid of "watterson_theta" and placeholder (NA) columns
                outsorted = outwatterson_theta.sort_values([4])  # sort by position
                # make sure sites, comparisons, missing get written as floats
                cols = [7, 8, 9, 10]
                outsorted[8] = outsorted[8].astype("Float64")  # raw_watterson_theta as float
                outsorted.to_csv(
                    outfile, sep="\t", mode="a", header=False, index=False, na_rep="NA"
                )  # write

        outfile.close()

    if PixyStat.TAJIMA_D in successful_stats:
        tajima_d_file = f"{pixy_args.output_prefix}_{PixyStat.TAJIMA_D.value}.txt"

        if os.path.exists(tajima_d_file):
            os.remove(tajima_d_file)

        outfile = open(tajima_d_file, "a")
        outfile.write(
            "pop"
            + "\t"
            + "chromosome"
            + "\t"
            + "window_pos_1"
            + "\t"
            + "window_pos_2"
            + "\t"
            + "tajima_d"
            + "\t"
            + "no_sites"
            + "\t"
            + "raw_pi"
            + "\t"
            + "raw_watterson_theta"
            + "\t"
            + "tajima_d_stdev"
            + "\n"
        )

        if aggregate:  # put winsizes back together for each population to make final_window_size
            for chromosome in chrom_list:
                outtajima_d = outgrouped.get_group(("tajima_d", chromosome)).reset_index(
                    drop=True
                )  # get this statistic, this chrom only
                outtajima_d.drop(
                    [0, 2], axis=1, inplace=True
                )  # get rid of "tajima_d" and placeholder (NA) columns
                outsorted = pixy.core.aggregate_output(
                    outtajima_d, PixyStat.TAJIMA_D.value, chromosome, window_size, args.fst_type
                )
                outsorted.to_csv(
                    outfile, sep="\t", mode="a", header=False, index=False, na_rep="NA"
                )  # write

        else:
            for chromosome in chrom_list:
                outtajima_d = outgrouped.get_group(("tajima_d", chromosome)).reset_index(
                    drop=True
                )  # get this statistic, this chrom only
                outtajima_d.drop(
                    [0, 2], axis=1, inplace=True
                )  # get rid of "theta" and placeholder (NA) columns
                outsorted = outtajima_d.sort_values([4])  # sort by position
                # make sure sites, comparisons, missing get written as floats
                cols = [7, 8, 9, 10]
                outsorted[cols] = outsorted[cols].astype("float64")
                outsorted.to_csv(
                    outfile, sep="\t", mode="a", header=False, index=False, na_rep="NA"
                )  # write

        outfile.close()
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
            "[pixy] WARNING: pixy failed to write any output files. Your VCF may not contain "
            "valid genotype data, or it was removed via filtering using the specified sites/bed "
            "file (if any)."
        )

    # print completion message
    end_time = time.time()
    logger.info(
        "[pixy] All calculations complete at "
        + time.strftime("%H:%M:%S on %Y-%m-%d", time.localtime(end_time))
    )
    total_time = time.time() - start_time
    logger.info("[pixy] Time elapsed: " + time.strftime("%H:%M:%S", time.gmtime(total_time)))
    logger.info(f"[pixy] Output files written to {pixy_args.output_dir}")

    if len(leftover_tmp_files) > 0:
        logger.info(
            f"[pixy] NOTE: There are pixy temp files in {pixy_args.output_dir}."
            "[pixy] If these are not actively being used (e.g. by another running pixy process), "
            "they can be safely deleted."
        )

    logger.info(
        f"[pixy] If you use pixy in your research, please cite the following paper: {citation}"
    )

    # restore output
    if args.silent:
        sys.stdout = sys.__stdout__


if __name__ == "__main__":
    main()
