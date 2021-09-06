#!/usr/bin/env python
# coding: utf-8



import allel
import numcodecs
import numpy as np
import sys
import os
import subprocess
import re
import operator
import pandas
import warnings
import time
import pathlib
import argparse
from argparse import RawTextHelpFormatter
import shutil

from multiprocess import Pool
from multiprocess import Queue
import multiprocess as mp

from scipy import special
from itertools import combinations
from collections import Counter

import pixy.core
import pixy.calc




# main pixy function

def main():
    # argument parsing via argparse

    # the ascii help image
    help_image = "█▀▀█ ░▀░ █░█ █░░█\n"     "█░░█ ▀█▀ ▄▀▄ █▄▄█\n"     "█▀▀▀ ▀▀▀ ▀░▀ ▄▄▄█\n"

    help_text = 'pixy: unbiased estimates of pi, dxy, and fst from VCFs with invariant sites'
    version = '1.2.4.beta1'
    citation = 'Korunes, KL and K Samuk. pixy: Unbiased estimation of nucleotide diversity and divergence in the presence of missing data. Mol Ecol Resour. 2021 Jan 16. doi: 10.1111/1755-0998.13326.'

    # initialize arguments
    parser = argparse.ArgumentParser(description=help_image+help_text+'\n'+version, formatter_class=argparse.RawTextHelpFormatter)

    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    additional = parser.add_argument_group('in addition, one of')
    optional = parser.add_argument_group('optional arguments')
    
    required.add_argument('--stats', nargs='+', choices=['pi', 'dxy', 'fst'], help='List of statistics to calculate from the VCF, separated by spaces.\ne.g. \"--stats pi fst\" will result in pi and fst calculations.', required=True)
    required.add_argument('--vcf', type=str, nargs='?', help='Path to the input VCF (bgzipped and tabix indexed).', required=True)
    required.add_argument('--populations', type=str, nargs='?', help='Path to a headerless tab separated populations file with columns [SampleID Population].', required=True)
    
    additional.add_argument('--window_size', type=int, nargs='?', help='Window size in base pairs over which to calculate stats.\nAutomatically determines window coordinates/bounds (see additional options below).', required=False)
    additional.add_argument('--bed_file', type=str, nargs='?', help='Path to a headerless .BED file with columns [chrom chromStart chromEnd].\nManually defines window bounds, which can be heterogenous in size.', required=False)

    optional.add_argument('--n_cores', type=int, nargs='?', default=1, help='Number of CPUs to utilize for parallel processing (default=1).', required=False)
    optional.add_argument('--output_folder', type=str, nargs='?', default='', help='Folder where output will be written, e.g. path/to/output_folder.\nDefaults to current working directory.', required=False)
    optional.add_argument('--output_prefix', type=str, nargs='?', default='pixy', help='Optional prefix for output file(s), with no slashes.\ne.g. \"--output_prefix output\" will result in [output folder]/output_pi.txt. \nDefaults to \'pixy\'.', required=False)
    optional.add_argument('--chromosomes', type=str, nargs='?', default='all', help='A single-quoted, comma separated list of chromosomes where stats will be calculated. \ne.g. --chromosomes \'X,1,2\' will restrict stats to chromosomes X, 1, and 2.\nDefaults to all chromosomes in the VCF.', required=False)
    optional.add_argument('--interval_start', type=str, nargs='?', help='The start of an interval over which to calculate stats.\nOnly valid when calculating over a single chromosome.\nDefaults to 1.', required=False)
    optional.add_argument('--interval_end', type=str, nargs='?', help='The end of the interval over which to calculate stats.\nOnly valid when calculating over a single chromosome.\nDefaults to the last position for a chromosome.', required=False)
    optional.add_argument('--sites_file', type=str, nargs='?', help='Path to a headerless tab separated file with columns [CHROM POS].\nThis defines sites where summary stats should be calculated.\nCan be combined with the --bed_file and --window_size options.', required=False)
    optional.add_argument('--chunk_size', type=int, nargs='?', default=100000, help='Approximate number of sites to read from VCF at any given time (default=100000).\nLarger numbers reduce I/O operations at the cost of memory.', required=False)
    
    optional.add_argument('--fst_type', choices=['wc', 'hudson'], default='wc', help='FST estimator to use, one of either: \n\'wc\' (Weir and Cockerham 1984) or\n\'hudson\' (Hudson 1992, Bhatia et al. 2013) \nDefaults to \'wc\'.', required=False)
    optional.add_argument('--bypass_invariant_check', choices=['yes', 'no'], default='no', help='Allow computation of stats without invariant sites (default=no).\nWill result in wildly incorrect estimates most of the time.\nUse with extreme caution!', required=False)
    optional.add_argument('--version', action='version', version=help_image+'version ' + version, help='Print the version of pixy in use.')
    optional.add_argument('--citation', action='version', version=citation, help='Print the citation for pixy.')
    optional.add_argument('--silent', action='store_true', help='Suppress all console output.')
    optional.add_argument('--debug', action='store_true', help=argparse.SUPPRESS)
    optional.add_argument('--keep_temp_file', action='store_true', help=argparse.SUPPRESS)


    
    args = parser.parse_args("--debug --chromosomes nDi.2.2.scaf00001,nDi.2.2.scaf00550 --window_size 10000 --stats pi fst --vcf testing/debugging_data/fst_bug/D_immitis.cohort.ALL.10K.vcf.gz --populations testing/debugging_data/fst_bug/aus_usa.txt --output_folder testing/notebook_output".split())
    
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
    print("[pixy] pixy " + version)
    print("[pixy] See documentation at https://pixy.readthedocs.io/en/latest/")
    popnames, popindices, chrom_list, IDs, temp_file, output_folder, output_prefix, bed_df, sites_df = pixy.core.check_and_validate_args(args)

    print("\n[pixy] Preparing for calculation of summary statistics: " + ', '.join(map(str, args.stats)))
    
    if 'fst' in args.stats:
        if args.fst_type == 'wc':
            fst_cite = 'Weir and Cockerham (1984)' 
        elif args.fst_type == 'hudson':
            fst_cite = 'Hudson (1992)' 
        print("[pixy] Using "+ fst_cite + "\'s estimator of FST.")
               
    print ("[pixy] Data set contains " + str(len(popnames)) + " population(s), " + str(len(chrom_list)) +" chromosome(s), and " + str(len(IDs)) + " sample(s)")

    if args.window_size is not None:
         print("[pixy] Window size: " + str(args.window_size) + " bp")
    
    if args.bed_file is not None:
         print("[pixy] Windows sourced from: " + args.bed_file)
    
    if args.sites_file is not None:
         print("[pixy] Calculations restricted to sites in: " + args.sites_file)
    
    print('')
    
    # time the calculations
    start_time = time.time()
    print("[pixy] Started calculations at " + time.strftime("%H:%M:%S on %Y-%m-%d", time.localtime(start_time)))
    print("[pixy] Using " + str(args.n_cores) + " out of " + str(mp.cpu_count()) + " available CPU cores\n")    
    # if in mc mode, set up multiprocessing 
    if (args.n_cores > 1):
        
        # use forking context on linux, and spawn otherwise (macOS)
        if sys.platform == "linux":
            ctx = mp.get_context("fork")
        else:
            ctx = mp.get_context("spawn")
        
        # set up the multiprocessing manager, queue, and process pool
        manager = ctx.Manager()
        q = manager.Queue()   
        pool = ctx.Pool(int(args.n_cores))
        
        # a listener function for writing a temp file
        # used to write output in multicore mode
        def listener(q, temp_file):

            with open(temp_file, 'a') as f:
                while 1:
                    m = q.get()
                    if m == 'kill':
                        break
                    f.write(str(m) + '\n')
                    f.flush()
        
        # launch the watcher function for collecting output
        watcher = pool.apply_async(listener, args = (q, temp_file,))
    
    # begin processing each chromosome

    for chromosome in chrom_list:

        print ("[pixy] Processing chromosome/contig " + chromosome + "...")
        
        # if not using a bed file, build windows manually
        if args.bed_file is None:
            
            # if an interval is specified, assign it
            if args.interval_start is not None:
                interval_start = int(args.interval_start)
                interval_end = int(args.interval_end)
                
            # otherwise, get the interval from the VCF's POS column
            else:
                chrom_max = subprocess.check_output("tabix " + args.vcf + " " + chromosome + " | cut -f 2 | tail -n 1", shell=True).decode("utf-8").split()
                interval_start = 1
                interval_end = int(chrom_max[0])

            # final check if intervals are valid
            try:     
                if (interval_start > interval_end):
                    raise ValueError()        
            except ValueError as e:
                raise Exception("[pixy] ERROR: The specified interval start ("+str(interval_start)+") exceeds the interval end ("+str(interval_end)+")") from e    
                
            targ_region = chromosome + ":" + str(interval_start) + "-" + str(interval_end)

            print("[pixy] Calculating statistics for region " + targ_region + "...")

            # Determine list of windows over which to compute stats

             # window size
            window_size = args.window_size
            
            # if the interval is smaller than one window, make a list of length 1
            if (interval_end - interval_start) <= window_size:
                window_pos_1_list = [interval_start]
                window_pos_2_list = [interval_start + window_size - 1]
            else:
                # if the interval_end is not a perfect multiple of the window size
                # bump the interval_end up to the nearest multiple of window size
                if not (interval_end % window_size == 0):
                    interval_end = interval_end + (window_size - (interval_end % window_size))
                
                # creat the start and stops for each window
                window_pos_1_list = [*range(interval_start, int(interval_end), window_size)]
                window_pos_2_list = [*range(interval_start + (window_size -1), int(interval_end) + window_size, window_size)]
            
            window_list = [list(a) for a in zip(window_pos_1_list,window_pos_2_list)]

            
            # Set aggregate to true if 1) the window size is larger than the chunk size OR 2) the window size wasn't specified, but the chrom is longer than the cutoff
            if (window_size > args.chunk_size) or ((args.window_size is None) and ((interval_end - interval_start) > args.chunk_size) ): 
                aggregate = True
            else:
                aggregate = False             
       
        # if using a bed file, subset the bed file for the current chromosome
        else:
            aggregate = False 
            bed_df_chrom = bed_df[bed_df['chrom'] == chromosome]
            window_list = [list(a) for a in zip(bed_df_chrom['chromStart'], bed_df_chrom['chromEnd'])]
        
        if(len(window_list) == 0):
            raise Exception("[pixy] ERROR: Window creation failed. Ensure that the POS column in the VCF is valid or change --window_size.")   
     
        # if aggregating, break down large windows into smaller windows
        if aggregate: 
            window_list = pixy.core.assign_subwindows_to_windows(window_list, args.chunk_size)
       
        # using chunk_size, assign  windows to chunks
        window_list = pixy.core.assign_windows_to_chunks(window_list, args.chunk_size)
            
        # if using a sites file, assign sites to chunks, as with windows above
        if args.sites_file is not None:
            sites_df = pandas.read_csv(args.sites_file, sep='\t', usecols=[0,1], names=['CHROM', 'POS'])
            sites_pre_list = sites_df[sites_df['CHROM'] == chromosome]
            sites_pre_list = sites_pre_list['POS'].tolist()
            sites_list = pixy.core.assign_sites_to_chunks(sites_pre_list, args.chunk_size)
        else:
            sites_df = None
        
        # obtain the list of chunks from the window list
        chunk_list = [i[2] for i in window_list]
        chunk_list = list(set(chunk_list))
        
        # if running in mc mode, send the summary stats funtion to the jobs pool
        if (args.n_cores > 1):
            
            # the list of jobs to be launched
            jobs = []
            
            for chunk in chunk_list:
                
                # create a subset of the window list specific to this chunk
                window_list_chunk = [x for x in window_list if x[2] == chunk]
                
                # and for the site list (if it exists)
                if sites_df is not None:
                    sites_list_chunk = [x for x in sites_list if x[1] == chunk]
                    sites_list_chunk = [x[0] for x in sites_list_chunk]
                else:
                    sites_list_chunk = None
                
                # determine the bounds of the chunk
                chunk_pos_1 = min(window_list_chunk, key=lambda x: x[1])[0]
                chunk_pos_2 = max(window_list_chunk, key=lambda x: x[1])[1]
                
                # launch a summary stats job for this chunk
                job = pool.apply_async(pixy.core.compute_summary_stats, args = (args, popnames, popindices, temp_file, chromosome, chunk_pos_1, chunk_pos_2, window_list_chunk, q, sites_list_chunk, aggregate, args.window_size))
                jobs.append(job)

            # launch all the jobs onto the pool
            for job in jobs: 
                job.get()
        
        # if running in single core mode, loop over the function manually
        elif (args.n_cores == 1):

            for chunk in chunk_list:
                
                 # create a subset of the window list specific to this chunk
                window_list_chunk = [x for x in window_list if x[2] == chunk]
                
                # and for the site list (if it exists)
                if sites_df is not None:
                    sites_list_chunk = [x for x in sites_list if x[1] == chunk]
                    sites_list_chunk = [x[0] for x in sites_list_chunk]
                else:
                    sites_list_chunk = None
                
                # determine the bounds of the chunk
                chunk_pos_1 = min(window_list_chunk, key=lambda x: x[1])[0]
                chunk_pos_2 = max(window_list_chunk, key=lambda x: x[1])[1]
                
                # don't use the queue (q) when running in single core mode
                q = "NULL"
                
                # compute summary stats for all windows in the chunk window list
                pixy.core.compute_summary_stats(args, popnames, popindices, temp_file, chromosome, chunk_pos_1, chunk_pos_2, window_list_chunk, q, sites_list_chunk, aggregate, args.window_size)
    
    # clean up any remaining jobs and stop the listener
    if (args.n_cores > 1):
        q.put('kill')
        pool.close()
        pool.join()
    
    # split and aggregate temp file to individual files
    
    # check if there is any output to process
    # halt execution if not
    try:
        outpanel = pandas.read_csv(temp_file, sep='\t', header=None)
    except pandas.errors.EmptyDataError:
        raise Exception('[pixy] ERROR: pixy failed to write any output. Confirm that your bed/sites files and intervals refer to existing chromosomes and positions in the VCF.')
        
    # check if particular stats failed to generate output     
    successful_stats = np.unique(outpanel[0])
    
    # if not all requested stats were generated, produce a warning 
    # and then remove the failed stats from the args list
    if set(args.stats) != set(successful_stats):
        missing_stats = list(set(args.stats) - set(successful_stats))
        print('\n[pixy] WARNING: pixy failed to find any valid gentoype data to calculate the following summary statistics: ' + ', '.join(missing_stats) + ". No output file will be created for these statistics.")
        args.stats = set(successful_stats)

    
    outpanel[3] = outpanel[3].astype(str) # force chromosome IDs to string
    outgrouped = outpanel.groupby([0,3]) #groupby statistic, chromosome

    
    # enforce chromosome IDs as strings
    chrom_list = list(map(str, chrom_list))
    
    if 'pi' in args.stats:
        stat = 'pi'
        pi_file = str(output_prefix) + "_pi.txt"

        if os.path.exists(pi_file):
            os.remove(pi_file)

        outfile = open(pi_file, 'a')
        outfile.write("pop" + "\t" + "chromosome" + "\t" + "window_pos_1" + "\t" + "window_pos_2" + "\t" + "avg_pi" + "\t" + "no_sites" + "\t" + "count_diffs" + "\t" + "count_comparisons" + "\t" + "count_missing" + "\n")
        
        if aggregate: #put winsizes back together for each population to make final_window_size
           
            for chromosome in chrom_list:
                outpi = outgrouped.get_group(("pi",chromosome)).reset_index(drop = True) #get this statistic, this chrom only
                outpi.drop([0,2], axis=1, inplace=True) #get rid of "pi" and placeholder (NA) columns
                outsorted = pixy.core.aggregate_output(outpi, stat, chromosome, window_size, args.fst_type)
                outsorted.to_csv(outfile, sep="\t", mode='a', header=False, index=False, na_rep='NA') #write
                
        else:
            for chromosome in chrom_list:
                outpi = outgrouped.get_group(("pi",chromosome)).reset_index(drop = True) #get this statistic, this chrom only
                outpi.drop([0,2], axis=1, inplace=True) #get rid of "pi" and placeholder (NA) columns
                outsorted = outpi.sort_values([4]) #sort by position
                # make sure sites, comparisons, missing get written as integers 
                cols = [7,8,9,10]
                outsorted[cols] = outsorted[cols].astype('Int64')
                outsorted.to_csv(outfile, sep="\t", mode='a', header=False, index=False, na_rep='NA') #write
        
        outfile.close()

    if 'dxy' in args.stats:
        stat = 'dxy'
        dxy_file = str(output_prefix) + "_dxy.txt"   

        if os.path.exists(dxy_file):
            os.remove(dxy_file)

        outfile = open(dxy_file, 'a')
        outfile.write("pop1" + "\t" + "pop2" + "\t" + "chromosome" + "\t" + "window_pos_1" + "\t" + "window_pos_2" + "\t" + "avg_dxy" + "\t" + "no_sites" + "\t" + "count_diffs" + "\t" + "count_comparisons" + "\t" + "count_missing" + "\n")

        if aggregate: # put winsizes back together for each population to make final_window_size
           
            for chromosome in chrom_list:
                outdxy = outgrouped.get_group(("dxy",chromosome)).reset_index(drop = True) #get this statistic, this chrom only
                outdxy.drop([0], axis=1, inplace=True) #get rid of "dxy"
                outsorted = pixy.core.aggregate_output(outdxy, stat, chromosome, window_size, args.fst_type)
                outsorted.to_csv(outfile, sep="\t", mode='a', header=False, index=False, na_rep='NA') #write
                
        else:
            for chromosome in chrom_list:
                outdxy = outgrouped.get_group(("dxy",chromosome)).reset_index(drop = True) #get this statistic, this chrom only
                outdxy.drop([0], axis=1, inplace=True) #get rid of "dxy"
                outsorted = outdxy.sort_values([4]) #sort by position
                # make sure sites, comparisons, missing get written as integers 
                cols = [7,8,9,10]
                outsorted[cols] = outsorted[cols].astype('Int64')
                outsorted.to_csv(outfile, sep="\t", mode='a', header=False, index=False, na_rep='NA') #write
        
        outfile.close()

    if 'fst' in args.stats:
        stat = 'fst'
        fst_file = str(output_prefix) + "_fst.txt"

        if os.path.exists(fst_file):
            os.remove(fst_file)

        outfile = open(fst_file, 'a')
        outfile.write("pop1" + "\t" + "pop2" + "\t" + "chromosome" + "\t" + "window_pos_1" + "\t" + "window_pos_2" + "\t" + "avg_" + args.fst_type + "_fst" + "\t" + "no_snps" + "\n")
        
        # keep track of chrosomes with no fst data
        chroms_with_no_data = []
 
        if aggregate: #put winsizes back together for each population to make final_window_size
           
            for chromosome in chrom_list:
                
                # logic to accodomate cases where pi/dxy have stats for a chromosome, but fst does not
                chromosome_has_data = True
                
                # if there are no valid fst estimates, set chromosome_has_data = False
                try:
                    outfst = outgrouped.get_group(("fst",chromosome)).reset_index(drop = True) #get this statistic, this chrom only
                except KeyError:
                    chroms_with_no_data.append(chromosome)
                    chromosome_has_data = False
                    
                    pass
                
                if chromosome_has_data:
                    outfst.drop([0], axis=1, inplace=True) #get rid of "fst"
                    outsorted = pixy.core.aggregate_output(outfst, stat, chromosome, window_size, args.fst_type)
                    outsorted = outsorted.iloc[:,:-3] #drop components (for now)
                    outsorted.to_csv(outfile, sep="\t", mode='a', header=False, index=False, na_rep='NA') #write
                
        else:
            for chromosome in chrom_list:
                
                # logic to accodomate cases where pi/dxy have stats for a chromosome, but fst does not
                chromosome_has_data = True
                
                # if there are no valid fst estimates, set chromosome_has_data = False
                try:
                    outfst = outgrouped.get_group(("fst",chromosome)).reset_index(drop = True) #get this statistic, this chrom only
                except KeyError:
                    chroms_with_no_data.append(chromosome)
                    chromosome_has_data = False
                    pass
                
                if chromosome_has_data:
                    outfst.drop([0], axis=1, inplace=True) #get rid of "fst" 
                    outsorted = outfst.sort_values([4]) #sort by position
                    # make sure sites (but not components like pi/dxy)
                    cols = [7]
                    outsorted[cols] = outsorted[cols].astype('Int64')
                    outsorted = outsorted.iloc[:,:-3] #drop components (for now)
                    outsorted.to_csv(outfile, sep="\t", mode='a', header=False, index=False, na_rep='NA')
        
        outfile.close()

        if len(chroms_with_no_data) >= 1:
            print('\n[pixy] NOTE: ' + 'The following chromosomes/scaffolds did not have sufficient data to estimate FST: ' + ', '.join(chroms_with_no_data))
    
    # remove the temp file(s)
    if (args.keep_temp_file is not True):
        os.remove(temp_file)
    
    # confirm output was generated successfully
    outfolder_files = [f for f in os.listdir(output_folder) if os.path.isfile(os.path.join(output_folder, f))]
    
    r = re.compile(".*_dxy.*|.*_pi.*|.*_fst.*")
    output_files = list(filter(r.match, outfolder_files)) 

    r = re.compile("pixy_tmpfile.*")
    leftover_tmp_files = list(filter(r.match, outfolder_files))

    if len(output_files) == 0:
        print('\n[pixy] WARNING: pixy failed to write any output files. Your VCF may not contain valid genotype data, or it was removed via filtering using the specified sites/bed file (if any).')
 
    # print completion message
    end_time = time.time()
    print("\n[pixy] All calculations complete at " + time.strftime("%H:%M:%S on %Y-%m-%d", time.localtime(end_time)))
    total_time = (time.time() - start_time)
    print("[pixy] Time elapsed: " + time.strftime("%H:%M:%S", time.gmtime(total_time)))
    print("[pixy] Output files written to: " + output_folder)
    
    if len(leftover_tmp_files) > 0:
        print("\n[pixy] NOTE: There are pixy temp files in " + output_folder)
        print("[pixy] If these are not actively being used (e.g. by another running pixy process), they can be safely deleted.")
    
    
    print("\n[pixy] If you use pixy in your research, please cite the following paper:\n[pixy] " + citation)
    
    # restore output
    if args.silent:
        sys.stdout = sys.__stdout__

if __name__ == "__main__":
            
   main()   
                






