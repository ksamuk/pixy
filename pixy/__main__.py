# -*- coding: utf-8 -*-

import allel
import zarr
import numcodecs
import numpy as np
import sys
import os
import shutil
import subprocess
import re
import operator
import pandas
from scipy import special
from itertools import combinations
from collections import Counter
import argparse


def main(args=None):

    if args is None:
        args = sys.argv[1:]    
    
    # argument parsing via argparse
    
    # the ascii help image
    help_image = "██████╗ ██╗██╗  ██╗██╗   ██╗\n" "██╔══██╗██║╚██╗██╔╝╚██╗ ██╔╝\n" "██████╔╝██║ ╚███╔╝  ╚████╔╝\n" "██╔═══╝ ██║ ██╔██╗   ╚██╔╝\n" "██║     ██║██╔╝ ██╗   ██║\n" "╚═╝     ╚═╝╚═╝  ╚═╝   ╚═╝\n" 
    
    help_text = 'pixy: senisble estimates of pi and dxy from a VCF'
    
    # initialize all the aruments
    parser = argparse.ArgumentParser(description=help_image+help_text, formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('--version', action='version', version='%(prog)s version 0.93')
    parser.add_argument('--stats', nargs='+', choices=['pi', 'dxy', 'fst'], help='Which statistics to calculate from the VCF (pi, dxy, and/or fst, separated by spaces)', required=True)
    parser.add_argument('--vcf', type=str, nargs='?', help='Path to the input VCF', required=True)
    parser.add_argument('--zarr_path', type=str, nargs='?', help='Folder in which to build the Zarr array(s)', required=True)
    parser.add_argument('--reuse_zarr', choices=['yes', 'no'], default='no', help='Use existing Zarr array(s) (saves time when re-running)')
    parser.add_argument('--populations', type=str, nargs='?', help='Path to the populations file', required=True)
    parser.add_argument('--window_size', type=int, nargs='?', help='Window size in base pairs over which to calculate pi/dxy')
    parser.add_argument('--chromosomes', type=str, nargs='?', default='all', help='A single-quoted, comma separated list of chromosome(s) (e.g. \'X,1,2\')', required=False)
    parser.add_argument('--interval_start', type=str, nargs='?', help='The start of the interval over which to calculate pi/dxy. Only valid when calculating over a single chromosome.')
    parser.add_argument('--interval_end', type=str, nargs='?', help='The end of the interval over which to calculate pi/dxy. Only valid when calculating over a single chromosome.')
    parser.add_argument('--variant_filter_expression', type=str, nargs='?', help='A single-quoted, comma separated list of genotype filters (e.g. \'DP>=10,GQ>=20\') to apply to SNPs', required=True)
    parser.add_argument('--invariant_filter_expression', type=str, nargs='?', help='A single-quoted, comma separated list of genotype filters (e.g. \'DP>=10,RGQ>=20\') to apply to invariant sites', required=True)
    parser.add_argument('--outfile_prefix', type=str, nargs='?', help='Path and prefix for the output file, e.g. path/to/outfile')
    parser.add_argument('--bypass_filtration', action='store_const', const='yes', default='no', help='Bypass all variant filtration (for data lacking FORMAT annotations, use with extreme caution)')
    parser.add_argument('--fst_maf_filter', default=0.0, type=float, nargs='?', help='Minor allele frequency filter for FST calculations, with value 0.0-1.0. Sites with MAF less than this value will be excluded.')
    
    ### test values for debugging
    
    # SIMULATED DATA
    #args = parser.parse_args('--interval_start 1000 --interval_end 8000 --bypass_filtration --stats pi fst dxy --vcf data/msprime_sim_invar/sim_dat_Ne=1.0e+06_mu=1e-08_samples=100_sites=10000_1_invar.vcf.gz --zarr_path data/msprime_sim_invar --chromosome 1 --window_size 1000 --populations data/msprime_sim_invar/populations.txt --regenerate_zarr yes --variant_filter_expression DP>=10,GQ>=20,RGQ>=20 --invariant_filter_expression DP>=10,RGQ>=20 --outfile_prefix output/pixy_out'.split())
    
    # ag1000g DATA
    #args = parser.parse_args('--stats pi --vcf data/vcf/multi_chr.vcf --zarr_path data/vcf/multi --window_size 10000 --populations data/vcf/ag1000/Ag1000_sampleIDs_popfile_3.txt --variant_filter_expression DP>=10,GQ>20 --invariant_filter_expression DP>=10,RGQ>20 --fst_maf_filter 0.05 --outfile_prefix output/pixy_out'.split())
    
    ###
    
    # catch arguments from the command line
    args = parser.parse_args()
    
    # CHECK FOR TABIX
    tabix_path = shutil.which("tabix")
    
    if tabix_path is None:
        raise Exception('ERROR: tabix is not installed! install with "conda install -c bioconda tabix"') 
    
    # VALIDATE ARGUMENTS
    
    # get vcf header info
    vcf_headers = allel.read_vcf_headers(args.vcf)
    
    # FILTER EXPRESSIONS
    # check if filter expressions are valid
    
    # get the list of format fields and requested filter fields
    format_fields = vcf_headers.formats.keys()
    filter_fields = list()
    for x in args.variant_filter_expression.split(","):
        filter_fields.append(re.sub("[^A-Za-z]+", "", x))
    
    missing = list(set(filter_fields)-set(format_fields))
        
    if len(missing) >0:
        raise Exception('ERROR: the following genotype filters were requested but not occur in the VCF:', missing) 
        
    # CHROMOSOMES
    # check if requested chromosomes exist in vcf
    
    # defaults to all the chromosomes contained in the VCF (first data column)
    if args.chromosomes == 'all':
        chrom_list = subprocess.check_output("grep -v '#' " + args.vcf + " | awk '{print $1}' | uniq", shell=True).decode("utf-8").split()
        chrom_all = chrom_list
    else:
        chrom_list = list(args.chromosomes.split(","))
        chrom_all = subprocess.check_output("grep -v '#' " + args.vcf + " | awk '{print $1}' | uniq", shell=True).decode("utf-8").split()
        missing = list(set(chrom_list)-set(chrom_all))
        if len(missing) >0:
            raise Exception('ERROR: the following chromosomes were requested but not occur in the VCF:', missing) 
    
    # INTERVALS
    # check if intervals are correctly specified
    
    if args.interval_start is not None and args.interval_end is None:
        raise Exception('ERROR: Both --interval_start and --interval_end must be specified') 
        
    if args.interval_start is None and args.interval_end is not None:
        raise Exception('ERROR: Both --interval_start and --interval_end must be specified') 
        
    if args.interval_start is not None and args.interval_end is not None and len(chrom_list) > 1:
        raise Exception('ERROR: --interval_start and --interval_end are not valid when calculating over multiple chromosomes. Remove both arguments or specify a single chromosome.') 
    
    # SAMPLES
    # check if requested samples exist in vcf
    
    # - parse + validate the population file
    # - format is IND POP (tab separated)
    # - throws an error if individuals are missing from VCF
    
    # read in the list of samples/populations
    poppanel = pandas.read_csv(args.populations, sep='\t', usecols=[0,1], names=['ID', 'Population'])
    poppanel.head()
    
    # get a list of samples from the callset
    samples_list = vcf_headers.samples
    
    # make sure every indiv in the pop file is in the VCF callset
    IDs = list(poppanel['ID'])
    missing = list(set(IDs)-set(samples_list))
    
    # find the samples in the callset index by matching up the order of samples between the population file and the callset
    # also check if there are invalid samples in the popfile
    try:
        samples_callset_index = [samples_list.index(s) for s in poppanel['ID']]
    except ValueError as e:
        raise Exception('ERROR: the following samples are listed in the population file but not in the VCF:', missing) from e
    else:   
        poppanel['callset_index'] = samples_callset_index
    
            # use the popindices dictionary to keep track of the indices for each population
        popindices={}
        popnames = poppanel.Population.unique()
        for name in popnames:
            popindices[name] = poppanel[poppanel.Population == name].callset_index.values
    
    print ("Preparing for calculation of summary statistics for " + str(len(chrom_list)) +" chromosome(s), and " + str(len(IDs)) + " sample(s)...")
    
    
    
    
    # initialize and remove any previous output files
    dxy_file = str(args.outfile_prefix) + "_dxy.txt"      
    if os.path.exists(dxy_file):
        os.remove(dxy_file)
        
    fst_file = str(args.outfile_prefix) + "_fst.txt"
    if os.path.exists(fst_file):
        os.remove(fst_file)
    
    pi_file = str(args.outfile_prefix) + "_pi.txt"
    if os.path.exists(pi_file):
        os.remove(pi_file)
        
    # initialize the output files for writing
    if 'pi' in args.stats:
        outfile = open(pi_file, 'a')
        outfile.write("pop" + "\t" + "chromosome" + "\t" + "window_pos_1" + "\t" + "window_pos_2" + "\t" + "avg_pi" + "\t" + "no_sites" + "\t" + "count_diffs" + "\t" + "count_comparisons" + "\t" + "count_missing" + "\n")
        outfile.close()
    
    if 'dxy' in args.stats:
        outfile = open(dxy_file, 'a')
        outfile.write("pop1" + "\t" + "pop2" + "\t" + "chromosome" + "\t" + "window_pos_1" + "\t" + "window_pos_2" + "\t" + "avg_dxy" + "\t" + "no_sites" + "\t" + "count_diffs" + "\t" + "count_comparisons" + "\t" + "count_missing" + "\n")
        outfile.close()
    
    if 'fst' in args.stats:
        outfile = open(fst_file, 'a')
        outfile.write("pop1" + "\t" + "pop2" + "\t" + "chromosome" + "\t" + "window_pos_1" + "\t" + "window_pos_2" + "\t" + "avg_wc_fst" + "\t" + "no_snps" + "\n")
        outfile.close()        
    
    
    
    
    # main loop for computing summary stats
    
    for chromosome in chrom_list:
        # Zarr array conversion
        print("chromosome " + chromosome + "...")
    
        # initialize the folder structure for the zarr array
        if os.path.exists(args.zarr_path) is not True:
            os.mkdir(args.zarr_path)
        
        # the chromosome specific zarr path
        zarr_path = args.zarr_path + "/" + chromosome
        
        # allow for resuse of previously calculated zarr arrays
        if os.path.exists(zarr_path) is not True:
            os.mkdir(zarr_path)
            allel.vcf_to_zarr(args.vcf, zarr_path, region=chromosome, fields='*', log=sys.stdout, overwrite=True, tabix=tabix_path)
        elif 'reuse_zarr' in args:
            if args.reuse_zarr == 'no':
                shutil.rmtree(zarr_path)
                os.mkdir(zarr_path)
                allel.vcf_to_zarr(args.vcf, zarr_path, region=chromosome, fields='*', log=sys.stdout, overwrite=True, tabix=tabix_path)
            if args.reuse_zarr == 'yes':
                print("Reusing existing zarr array for chromosome " + chromosome + "...")
    
        # inspect the structure of the zarr data
        callset = zarr.open_group(zarr_path, mode='r')
    
        # parse the filtration expression and build the boolean filter array
    
        # define an operator dictionary for parsing the operator strings
        ops = { "<": operator.lt, "<=": operator.le, ">": operator.gt, ">=": operator.ge, "==": operator.eq}
    
        # determine the complete list of available calldata fields usable for filtration
        calldata_fields = sorted(callset['/calldata/'].array_keys())
    
        # check if bypassing filtration, otherwise filter
        if args.bypass_filtration=='no':
            # intialize the filtration array (as a list)
            filters = []
    
            #print("Creating filters...")
    
            # split the filtration expressions by commas
            # parse out each component
            # apply the filter
            # TBD: handle cases where the specified filter isn't in the callset
            for x in args.variant_filter_expression.split(","):
                stat = re.sub("[^A-Za-z]+", "", x)
                value = int(re.sub("[^0-9]+", "", x))
                compare = re.sub("[A-Za-z0-9]+", "", x)
    
                #print(str(stat) + " " + str(compare) + " " + str(value))
    
                # check if the requested annotation exists in the VCF
                try: 
                    stat_index = calldata_fields.index(stat)
                except ValueError as e:
                    raise Exception("Error: The requested filter \'" + stat + "\' is not annotated in the input VCF") from e
                else: 
    
                    # check if this is the first run through the loop
                    # if so, either catch GQ and RGQ as separate filters or create the initial filter
                    # on subsequent runs ('filters' is not a list), only catch GQ and RGQ or update the filter (logical AND)
                    if type(filters) is list:
                        if stat == 'GQ':
                            GQ_filter = ops[compare](callset['/calldata/' + stat][:], value)
                        elif stat == 'RGQ':
                            RGQ_filter = ops[compare](callset['/calldata/' + stat][:], value)
                        else:
                            # this creates the initial filter array
                            filters = ops[compare](callset['/calldata/' + stat][:], value)
                    elif type(filters) is not list:
                        if stat == 'GQ':
                            GQ_filter = ops[compare](callset['/calldata/' + stat][:], value)
                        elif stat == 'RGQ':
                            RGQ_filter = ops[compare](callset['/calldata/' + stat][:], value)
                        else:
                            # this updates the filter array with additional filtration criteria
                            filters = np.logical_and(filters, ops[compare](callset['/calldata/' + stat][:], value))
    
            # check if GQ and RQG exist
            # if they both exist, perform a logical OR and join them into the filter
            # otherwise, perform a logical AND to join either one into the filter
    
            GQ_exists = 'GQ_filter' in args
            RGQ_exists = 'RGQ_filter' in args
    
            if GQ_exists & RGQ_exists:
                filters = np.logical_and(filters, np.logical_or(GQ_filter, RGQ_filter))
            elif GQ_exists:
                filters = np.logical_and(filters, GQ_filter)
            elif RGQ_exists:
                filters = np.logical_and(filters, RGQ_filter)
    
            # finally, invert the whole array 
            # this is for convenience/brevity in the next section
            filters = np.invert(filters)
    
        # applying the filter to the data
        # all the filters are in a boolean array ('filters') 
    
        # recode the gt matrix as a Dask array (saves memory)
        gt_dask = allel.GenotypeDaskArray(callset['/calldata/GT'])
    
        # create a packed genotype array 
        # this is a array with dims snps x samples
        # genotypes are represented by single byte codes 
        # critically, as the same dims as the filters array below
        gt_array = allel.GenotypeArray(gt_dask).to_packed()
    
        # only apply filter if not bypassing filtration
        if args.bypass_filtration=='no':
            # set all genotypes that fail filters to 'missing'
            # 239 = -1 (i.e. missing) for packed arrays
            gt_array[filters] = 239
    
        # remove sites with >1 alt allele from the filtered gt array
        gt_array = np.delete(gt_array, np.where(callset['/variants/numalt'][:] > 1), axis=0)
    
        # convert the packed array back to a GenotypeArray
        gt_array = allel.GenotypeArray.from_packed(gt_array)
    
        # build the position array
        pos_array = allel.SortedIndex(callset['/variants/POS'])
    
        # remove everything but biallelic snps and monomorphic sites from the position array
        pos_array = pos_array[callset['/variants/numalt'][:] < 2]
    
        #Basic functions for comparing the genotypes at each site in a region: counts differences out of sites with data
    
        #For the given region: return average pi, # of differences, # of comparisons, and # missing.
        # this function loops over every site in a region passed to it
    
        #Basic functions for comparing the genotypes at each site in a region: counts differences out of sites with data
    
        #For the given region: return average pi, # of differences, # of comparisons, and # missing.
        # this function loops over every site in a region passed to it
        def tallyRegion(gt_region):
            total_diffs = 0
            total_comps = 0
            total_missing = 0
            for site in gt_region:
                vec = site.flatten()
                #now we have an individual site as a numpy.ndarray, pass it to the comparison function
                site_diffs, site_comps, missing = compareGTs(vec)
                total_diffs += site_diffs
                total_comps += site_comps
                total_missing += missing
            if total_comps > 0:
                avg_pi = total_diffs/total_comps
            else:
                avg_pi = 0
            return(avg_pi, total_diffs, total_comps, total_missing)
    
        #For the given region: return average dxy, # of differences, # of comparisons, and # missing.
        # this function loops over every site in a region passed to it
        def dxyTallyRegion(pop1_gt_region, pop2_gt_region):
            total_diffs = 0
            total_comps = 0
            total_missing = 0
            for x in range(0,len(pop1_gt_region)):
                site1 = pop1_gt_region[x]
                site2 = pop2_gt_region[x]
                vec1 = site1.flatten()
                vec2 = site2.flatten()
                #now we have an individual site as 2 numpy.ndarrays, pass them to the comparison function
                site_diffs, site_comps, missing = dxyCompareGTs(vec1, vec2)
                total_diffs += site_diffs
                total_comps += site_comps
                total_missing += missing
            if total_comps > 0:
                avg_pi = total_diffs/total_comps
            else:
                avg_pi = 0
            return(avg_pi, total_diffs, total_comps, total_missing)
    
        #Return the number of differences, the number of comparisons, and missing data count.
        def compareGTs(vec): #for pi
            c = Counter(vec)
            diffs = c[1]*c[0]
            gts = c[1]+c[0]
            missing = (len(vec))-gts  #anything that's not 1 or 0 is ignored and counted as missing
            comps = int(special.comb(gts,2))
            return(diffs,comps,missing)
    
        def dxyCompareGTs(vec1, vec2): #for dxy
            c1 = Counter(vec1)
            c2 = Counter(vec2)
            gt1zeros = c1[0]
            gt1ones = c1[1]
            gts1 = c1[1]+c1[0]
            gt2zeros = c2[0]
            gt2ones = c2[1]
            gts2 = c2[1]+c2[0]
            missing = (len(vec1)+len(vec2))-(gts1+gts2)  #anything that's not 1 or 0 is ignored and counted as missing  
            diffs = (gt1zeros*gt2ones)+(gt1ones*gt2zeros)
            comps = gts1*gts2
            return(diffs,comps,missing)
    
        # Interval specification check
        # check if computing over specific intervals (otherwise, compute over whole chromosome)
    
        # window size
        window_size = args.window_size
    
        # set intervals based on args
        if (args.interval_end is None):
            interval_end = max(pos_array)
        else:
            interval_end = int(args.interval_end)
    
        if (args.interval_start is None):
                interval_start = 1
        else:
            interval_start = int(args.interval_start)
            
        try:    
            if (interval_start > interval_end):
                raise ValueError()        
        except ValueError as e:
            raise Exception("The specified interval start ("+str(interval_start)+") exceeds the interval end ("+str(interval_end)+")") from e
    
        # catch misspecified intervals
        try:    
            if (interval_end > max(pos_array)):
                raise ValueError()        
        except ValueError as e:
            raise Exception("The specified interval end ("+str(interval_end)+") exceeds the last position of the chromosome ("+str(max(pos_array))+")") from e
    
        try:           
            if (interval_start < min(pos_array)):
                raise ValueError()      
        except ValueError as e:
            raise Exception("The specified interval start ("+str(interval_start)+") begins before the first position of the chromosome ("+str(min(pos_array))+")") from e
    
        try:           
            if ((interval_end - interval_start + 1) < window_size):
                raise ValueError()      
        except ValueError as e:
            raise Exception("The specified interval ("+str(interval_start)+"-"+str(interval_end)+") is too small for the requested window size ("+str(window_size)+")") from e
    
        # PI:
        # AVERAGE NUCLEOTIDE VARIATION WITHIN POPULATIONS
    
        # Compute pi over a chosen interval and window size
    
        # TBD:
        # - write out summary of program parameters file* think about how this should work
    
    
        if (args.populations is not None) and ('pi' in args.stats):
    
            # open the pi output file for writing
            outfile = open(pi_file, 'a')
    
            for pop in popnames:
    
                # window size:
                window_size = args.window_size
    
                # initialize window_pos_2 
                window_pos_2 = (interval_start + window_size)-1
    
    
                # loop over populations and windows, compute stats and write to file
                for window_pos_1 in range(interval_start, interval_end, window_size):
    
                    # pull out the genotypes for the window
                    loc_region = pos_array.locate_range(window_pos_1, window_pos_2)
                    gt_region1 = gt_array[loc_region]
    
    
                    # subset the window for the individuals in each population 
                    gt_pop = gt_region1.take(popindices[pop], axis=1)
    
                    avg_pi, total_diffs, total_comps, total_missing = tallyRegion(gt_pop)
                    outfile.write(str(pop) + "\t" + str(chromosome) + "\t" + str(window_pos_1) + "\t" + str(window_pos_2) + "\t" + str(avg_pi) + "\t" + str(len(gt_region1)) + "\t" + str(total_diffs) + "\t" + str(total_comps) + "\t" + str(total_missing) + "\n")
                    window_pos_2 += window_size
    
                # close output file and print complete message
            outfile.close()
    
            print("Pi calculations for chromosome " + chromosome + " complete and written to " + args.outfile_prefix + "_pi.txt")
    
    
        # DXY:
        # AVERAGE NUCLEOTIDE VARIATION BETWEEN POPULATIONS
    
    
        if (args.populations is not None) and ('dxy' in args.stats):
    
            # create a list of all pairwise comparisons between populations in the popfile
            dxy_pop_list = list(combinations(popnames, 2))
    
            # open the dxy output file for writing
            outfile = open(dxy_file, 'a')
            
            # interate over all population pairs and compute dxy
            for pop_pair in dxy_pop_list:
                pop1 = pop_pair[0]
                pop2 = pop_pair[1]
    
                # window size:
                window_size = args.window_size
    
                # initialize window_pos_2 
                window_pos_2 = interval_start + window_size
    
                # perform the dxy calculation for all windows in the range
                for window_pos_1 in range (interval_start, interval_end, window_size):
                    loc_region = pos_array.locate_range(window_pos_1, window_pos_2)
                    gt_region1 = gt_array[loc_region]
                    # use the popGTs dictionary to keep track of this region's GTs for each population
                    popGTs={}
                    for name in pop_pair:
                        gt_pop = gt_region1.take(popindices[name], axis=1)
                        popGTs[name] = gt_pop
    
                    pop1_gt_region1 = popGTs[pop1]
                    pop2_gt_region1 = popGTs[pop2]
                    avg_dxy, total_diffs, total_comps, total_missing = dxyTallyRegion(pop1_gt_region1, pop2_gt_region1)
                    outfile.write(str(pop1) + "\t" + str(pop2) + "\t" + str(chromosome) + "\t" + str(window_pos_1) + "\t" + str(window_pos_2) + "\t" + str(avg_dxy) + "\t" + str(len(gt_region1)) + "\t" + str(total_diffs) + "\t" + str(total_comps) + "\t" + str(total_missing) + "\n")
    
                    window_pos_2 += window_size
    
            outfile.close()
            print("Dxy calculations chromosome " + chromosome + " complete and written to " + args.outfile_prefix + "_dxy.txt")
    
        # FST:
        # WEIR AND COCKERHAMS FST
        # This is just a plain wrapper for the scikit-allel fst function
    
        if (args.populations is not None) and ('fst' in args.stats):
            
            # open the fst output file for writing
            outfile = open(fst_file, 'a')
            
            # determine all the possible population pairings
            pop_names=list(popindices.keys())
            fst_pop_list = list(combinations(pop_names, 2))
    
            # for each pair, compute fst
            for pop_pair in fst_pop_list:
    
                # the indices for the individuals in each population
                fst_pop_indicies=[popindices[pop_pair[0]].tolist(), popindices[pop_pair[1]].tolist()]
    
                # rebuild the GT array for variant sites only
                # this is done for speed (also FST is undefined for monomorphic sites)
                gt_dask_fst = allel.GenotypeDaskArray(callset['/calldata/GT'])
                gt_array_fst = allel.GenotypeArray(gt_dask_fst).to_packed()
    
                # flag & remove non-biallelic sites
                non_bial_sites = np.where(np.logical_not(callset['/variants/numalt'][:] == 1))
                gt_array_fst = np.delete(gt_array_fst, non_bial_sites, axis=0)
                gt_array_fst = allel.GenotypeArray.from_packed(gt_array_fst)
            
                #apply the maf filter (default is 0, i.e. no filter)
                allele_counts= gt_array_fst.count_alleles()
                allele_freqs = allele_counts.to_frequencies()
                gt_array_fst = np.delete(gt_array_fst, np.where(allele_freqs < args.fst_maf_filter), axis=0)
    
                # also rebuild the position array for biallelic sites/maf filtered sites only
                pos_array_fst = allel.SortedIndex(callset['/variants/POS'])
                pos_array_fst = pos_array_fst[callset['/variants/numalt'][:] == 1]
                pos_array_fst = np.delete(pos_array_fst, np.where(allele_freqs < args.fst_maf_filter), axis=0)
         
                # compute FST
                # windowed_weir_cockerham_fst seems to generate (spurious?) warnings about div/0, so suppressing warnings
                # (this assumes that the scikit-allel function is working as intended)
                np.seterr(divide='ignore', invalid='ignore')
     
                a,b,c=allel.windowed_weir_cockerham_fst(pos_array_fst, gt_array_fst, subpops=fst_pop_indicies, size=args.window_size, start=interval_start, stop=interval_end)
    
                for fst,wind,snps in zip(a, b, c):
                    outfile.write(str(pop_pair[0]) + "\t" + str(pop_pair[1]) + "\t" + str(chromosome) + "\t" + str(wind[0]) + "\t" + str(wind[1]) + "\t" + str(fst) + "\t" + str(snps) +"\n")
            outfile.close()
            print("Fst calculations chromosome " + chromosome + " complete and written to " + args.outfile_prefix + "_fst.txt")
    
    
    
    
    
    
if __name__ == "__main__":
    main()