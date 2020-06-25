# -*- coding: utf-8 -*-

import allel
import zarr
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

from scipy import special
from itertools import combinations
from collections import Counter

def main(args=None):

    if args is None:
        args = sys.argv[1:]    version_text = 'version 0.95.0'
    
    # initialize arguments
    parser = argparse.ArgumentParser(description=help_image+help_text+'\n'+version_text, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--version', action='version', version=version_text)
    parser.add_argument('--stats', nargs='+', choices=['pi', 'dxy', 'fst'], help='Which statistics to calculate from the VCF (pi, dxy, and/or fst, separated by spaces)', required=True)
    parser.add_argument('--vcf', type=str, nargs='?', help='Path to the input VCF', required=True)
    parser.add_argument('--zarr_path', type=str, nargs='?', help='Folder in which to build the Zarr array(s)', required=True)
    parser.add_argument('--reuse_zarr', choices=['yes', 'no'], default='no', help='Use existing Zarr array(s) (saves time if re-running)')
    parser.add_argument('--populations', type=str, nargs='?', help='Path to the populations file', required=True)
    parser.add_argument('--window_size', type=int, nargs='?', help='Window size in base pairs over which to calculate pi/dxy')
    parser.add_argument('--chromosomes', type=str, nargs='?', default='all', help='A single-quoted, comma separated list of chromosome(s) (e.g. \'X,1,2\')', required=False)
    parser.add_argument('--interval_start', type=str, nargs='?', help='The start of the interval over which to calculate pi/dxy. Only valid when calculating over a single chromosome.')
    parser.add_argument('--interval_end', type=str, nargs='?', help='The end of the interval over which to calculate pi/dxy. Only valid when calculating over a single chromosome.')
    parser.add_argument('--variant_filter_expression', type=str, nargs='?', help='A single-quoted, comma separated list of genotype filters (e.g. \'DP>=10,GQ>=20\') to apply to SNPs', required=False)
    parser.add_argument('--invariant_filter_expression', type=str, nargs='?', help='A single-quoted, comma separated list of genotype filters (e.g. \'DP>=10,RGQ>=20\') to apply to invariant sites', required=False)
    parser.add_argument('--outfile_prefix', type=str, nargs='?', default='./pixy_output', help='Path and prefix for the output file, e.g. path/to/outfile')
    parser.add_argument('--bypass_filtration', choices=['yes', 'no'], default='no', help='Bypass all variant filtration (for data lacking FORMAT fields, use with caution)')
    parser.add_argument('--bypass_invariant_check', choices=['yes', 'no'], default='no', help='Allow computation of stats without invariant sites, will result in wildly incorrect estimates most of the time. Use with extreme caution.')
    parser.add_argument('--fst_maf_filter', default=0.05, type=float, nargs='?', help='Minor allele frequency filter for FST calculations, with value 0.0-1.0 (default 0.05).')
    
    # ag1000g test data
    # args = parser.parse_args('--stats fst --vcf data/vcf/multi_chr.vcf.gz --zarr_path data/vcf/multi --window_size 10000 --populations data/vcf/ag1000/Ag1000_sampleIDs_popfile_3.txt --variant_filter_expression DP>=10,GQ>20 --invariant_filter_expression DP>=10,RGQ>20 --outfile_prefix output/pixy_out'.split())
    
    # filter test data
    # args = parser.parse_args('--stats pi --vcf data/vcf/filter_test.vcf.gz --zarr_path data/vcf/filter_test --window_size 3 --populations data/vcf/ag1000/Ag1000_sampleIDs_popfile_3.txt --variant_filter_expression DP>=10,GQ>20 --invariant_filter_expression DP>=10,RGQ>20 --fst_maf_filter 0.05 --outfile_prefix output/pixy_out'.split())
    
    # catch arguments from the command line
    args = parser.parse_args()
    
    
    
    
    # CHECK FOR TABIX
    # (disabled until we implement site level and BED support)
    #tabix_path = shutil.which("tabix")
    
    #if tabix_path is None:
    #    warnings.warn('[pixy] WARNING: tabix is not installed (or cannot be located) -- this may reduce performance. install tabix with "conda install -c bioconda tabix"') 
    #if not os.path.exists(args.vcf + ".tbi") and tabix_path is not None:
    #    raise Exception('[pixy] ERROR: your vcf is not indexed with tabix, index the bgzipped vcf with "tabix your.vcf.gz"') 
    
    
    
    
    # VALIDATE ARGUMENTS
    
    print("[pixy] pixy " + version_text)
    print("[pixy] Validating VCF and input parameters (this may take some time)...")
    
    # expand all file paths
    args.vcf = os.path.expanduser(args.vcf)
    args.zarr_path = os.path.expanduser(args.zarr_path)
    args.populations = os.path.expanduser(args.populations)
    args.outfile_prefix = os.path.expanduser(args.outfile_prefix)
    
    # CHECK FOR EXISTANCE OF VCF AND POPFILES
    
    if os.path.exists(args.vcf) is not True:
        raise Exception('[pixy] ERROR: The specified VCF ' + str(args.vcf) + ' does not exist') 
    
    if os.path.exists(args.populations) is not True:
        raise Exception('[pixy] ERROR: The specified populations file ' + str(args.populations) + ' does not exist') 
    
    # VALIDATE FILTER EXPRESSIONS
    
    # get vcf header info
    vcf_headers = allel.read_vcf_headers(args.vcf)
    
    # skip invariant check if only asking for FST
    if len(args.stats) == 1 and (args.stats[0] == 'fst'):
        args.bypass_invariant_check = "yes"
    
    # if we are bypassing the invariant check, spoof in a invariant filter    
    if args.bypass_invariant_check=="yes":
        args.invariant_filter_expression="DP>=0"
    
    if args.bypass_filtration=='no' and (args.variant_filter_expression is None or args.invariant_filter_expression is None):
        raise Exception('[pixy] ERROR: One or more filter expression is missing. Provide two filter expressions, or set --bypass_filtration to \'yes\'') 
    
    if args.bypass_filtration=='no':
        # get the list of format fields and requested filter fields
        format_fields = vcf_headers.formats.keys()
        filter_fields = list()
    
        for x in args.variant_filter_expression.split(","):
            filter_fields.append(re.sub("[^A-Za-z]+", "", x))
    
        for x in args.invariant_filter_expression.split(","):
            filter_fields.append(re.sub("[^A-Za-z]+", "", x))
    
        missing = list(set(filter_fields)-set(format_fields))
    
        if len(missing) >0:
            raise Exception('[pixy] ERROR: the following genotype filters were requested but not occur in the VCF: ', missing) 
    else:
        print("[pixy] WARNING: --bypass_filtration is set to \'yes\', genotype filtration will be not be performed.")
    
    
    # VALIDATE THE VCF
    
    # check if the vcf is zipped
    if re.search(".gz", args.vcf):
        cat_prog = "gunzip -c "
    else:
        cat_prog = "cat "
    
    # check if the vcf contains any invariant sites
    # a very basic check: just looks for at least one invariant site in the alt field
    
    if args.bypass_invariant_check=='no':
        alt_list = subprocess.check_output(cat_prog + args.vcf + " | grep -v '#' | head -n 10000 | awk '{print $5}' | sort | uniq", shell=True).decode("utf-8").split()
        if "." not in alt_list:
            raise Exception('[pixy] ERROR: the provided VCF appears to contain no invariant sites (ALT = \".\"). This check can be bypassed via --bypass_invariant_check \'yes\'.') 
    else:
        if not (len(args.stats) == 1 and (args.stats[0] == 'fst')):
            print("[pixy] EXTREME WARNING: --bypass_invariant_check is set to \'yes\', which assumes that your VCF contains invariant sites. Lack of invariant sites will result in incorrect estimates.")
    
    
    # check if requested chromosomes exist in vcf
    # defaults to all the chromosomes contained in the VCF (first data column)
    
    if args.chromosomes == 'all':
        chrom_list = subprocess.check_output(cat_prog + args.vcf + " | grep -v '#' | awk '{print $1}' | uniq", shell=True).decode("utf-8").split()
        chrom_all = chrom_list
    
    if args.chromosomes == 'all':
        chrom_list = subprocess.check_output(cat_prog + args.vcf + " | grep -v '#' | awk '{print $1}' | uniq", shell=True).decode("utf-8").split()
        chrom_all = chrom_list
    else:
        chrom_list = list(args.chromosomes.split(","))
        chrom_all = subprocess.check_output(cat_prog + args.vcf + " | grep -v '#' | awk '{print $1}' | uniq", shell=True).decode("utf-8").split()
        missing = list(set(chrom_list)-set(chrom_all))
        if len(missing) >0:
            raise Exception('[pixy] ERROR: the following chromosomes were requested but not occur in the VCF: ', missing) 
    
    # INTERVALS
    # check if intervals are correctly specified
    
    if args.interval_start is not None and args.interval_end is None:
        raise Exception('[pixy] ERROR: Both --interval_start and --interval_end must be specified') 
        
    if args.interval_start is None and args.interval_end is not None:
        raise Exception('[pixy] ERROR: Both --interval_start and --interval_end must be specified') 
        
    if args.interval_start is not None and args.interval_end is not None and len(chrom_list) > 1:
        raise Exception('[pixy] ERROR: --interval_start and --interval_end are not valid when calculating over multiple chromosomes. Remove both arguments or specify a single chromosome.') 
    
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
        raise Exception('[pixy] ERROR: the following samples are listed in the population file but not in the VCF: ', missing) from e
    else:   
        poppanel['callset_index'] = samples_callset_index
    
            # use the popindices dictionary to keep track of the indices for each population
        popindices={}
        popnames = poppanel.Population.unique()
        for name in popnames:
            popindices[name] = poppanel[poppanel.Population == name].callset_index.values
    
    print ("[pixy] Preparing for calculation of summary statistics: " + ','.join(map(str, args.stats)))
    print ("[pixy] Data set contains " + str(len(popnames)) + " population(s), " + str(len(chrom_list)) +" chromosome(s), and " + str(len(IDs)) + " sample(s)")
    
    
    
    
    # initialize and remove any previous output files
    if os.path.exists(re.sub(r"[^\/]+$", "", args.outfile_prefix)) is not True:
        os.mkdir(re.sub(r"[^\/]+$", "", args.outfile_prefix))
    
            
    # initialize the output files for writing
    if 'pi' in args.stats:
        
        pi_file = str(args.outfile_prefix) + "_pi.txt"
        
        if os.path.exists(pi_file):
            os.remove(pi_file)
            
        outfile = open(pi_file, 'a')
        outfile.write("pop" + "\t" + "chromosome" + "\t" + "window_pos_1" + "\t" + "window_pos_2" + "\t" + "avg_pi" + "\t" + "no_sites" + "\t" + "count_diffs" + "\t" + "count_comparisons" + "\t" + "count_missing" + "\n")
        outfile.close()
    
    if 'dxy' in args.stats:
        
        dxy_file = str(args.outfile_prefix) + "_dxy.txt"     
        
        if os.path.exists(dxy_file):
            os.remove(dxy_file)
            
        outfile = open(dxy_file, 'a')
        outfile.write("pop1" + "\t" + "pop2" + "\t" + "chromosome" + "\t" + "window_pos_1" + "\t" + "window_pos_2" + "\t" + "avg_dxy" + "\t" + "no_sites" + "\t" + "count_diffs" + "\t" + "count_comparisons" + "\t" + "count_missing" + "\n")
        outfile.close()
    
    if 'fst' in args.stats:
        
        fst_file = str(args.outfile_prefix) + "_fst.txt"
        
        if os.path.exists(fst_file):
            os.remove(fst_file)
            
        outfile = open(fst_file, 'a')
        outfile.write("pop1" + "\t" + "pop2" + "\t" + "chromosome" + "\t" + "window_pos_1" + "\t" + "window_pos_2" + "\t" + "avg_wc_fst" + "\t" + "no_snps" + "\n")
        outfile.close()
    
    # initialize the folder structure for the zarr array
    if os.path.exists(args.zarr_path) is not True:
        pathlib.Path(args.zarr_path).mkdir(parents=True, exist_ok=True)
    
    
    
    
    # main loop for computing summary stats
    
    # time the calculations
    start_time = time.time()
    print("[pixy] Started calculations at " + time.strftime("%H:%M:%S", time.localtime(start_time)))
    
    for chromosome in chrom_list:
         
        # Zarr array conversion
        
        # the chromosome specific zarr path
        zarr_path = args.zarr_path + "/" + chromosome
        
        # determine the fields that will be included
        # TBD: just reading all fields currently
        # vcf_fields = ['variants/CHROM', 'variants/POS'] + ['calldata/' + s for s in np.unique(filter_fields)]
        
        # build region string (if using an interval)
        if args.interval_start is not None:
            targ_region = chromosome + ":" + str(args.interval_start) + "-" + str(args.interval_end)
        else:
            targ_region = chromosome
        
        # allow for resuse of previously calculated zarr arrays
        if args.reuse_zarr == 'yes' and os.path.exists(zarr_path):
            print("[pixy] If a zarr array exists, it will be reused for chromosome " + chromosome + "...")
        elif args.reuse_zarr == 'no' or os.path.exists(zarr_path) is not True:
            print("[pixy] Building zarr array for chromosome " + chromosome + "...")
            warnings.filterwarnings("ignore")
            allel.vcf_to_zarr(args.vcf, zarr_path, region=targ_region, fields='*', overwrite=True)
            warnings.resetwarnings()
            
        print("[pixy] Calculating statistics for chromosome " + targ_region + "...")
        
        # open the zarr
        callset = zarr.open_group(zarr_path, mode='r')
    
        # parse the filtration expression and build the boolean filter array
    
        # define an operator dictionary for parsing the operator strings
        ops = { "<": operator.lt, "<=": operator.le, ">": operator.gt, ">=": operator.ge, "==": operator.eq}
    
        # determine the complete list of available calldata fields usable for filtration
        calldata_fields = sorted(callset['/calldata/'].array_keys())
    
        # check if bypassing filtration, otherwise filter
        if args.bypass_filtration=='no':
            
            # VARIANT SITE FILTERS
            var_filters = []
            
            # iterate over each requested variant filter
            for x in args.variant_filter_expression.split(","):
                stat = re.sub("[^A-Za-z]+", "", x)
                value = int(re.sub("[^0-9]+", "", x))
                compare = re.sub("[A-Za-z0-9]+", "", x)
    
                # check if the requested filter/format exists in the VCF
                try: 
                    stat_index = calldata_fields.index(stat)
                except ValueError as e:
                    raise Exception("[pixy] ERROR: The requested filter \'" + stat + "\' is not annotated in the input VCF FORMAT field") from e
                else: 
                    if type(var_filters) is list:
                        var_filters = ops[compare](callset['/calldata/' + stat][:], value)
                    elif type(var_filters) is not list:
                        var_filters = np.logical_and(var_filters, ops[compare](callset['/calldata/' + stat][:], value))
            
            # create a mask for variants only
            # is snp is a site level (1d) array
            # np.tile below creates a column of "is_snp" once for each sample 
            # (i.e. makes it the same dimensions as the genotype table)
            is_snp = np.array([callset['/variants/is_snp'][:].flatten()]).transpose()
            snp_mask = np.tile(is_snp, (1,var_filters.shape[1]))
            
            # force only variant sites (snps, remember we ignore indels) to be included in the filter
            var_filters = np.logical_and(var_filters, snp_mask)
            
            # INVARIANT SITE FILTERS
            invar_filters = []
            
            for x in args.invariant_filter_expression.split(","):
                stat = re.sub("[^A-Za-z]+", "", x)
                value = int(re.sub("[^0-9]+", "", x))
                compare = re.sub("[A-Za-z0-9]+", "", x)
    
                # check if the requested filter/format exists in the VCF
                try: 
                    stat_index = calldata_fields.index(stat)
                except ValueError as e:
                    raise Exception("[pixy] ERROR: The requested filter \'" + stat + "\' is not annotated in the input VCF") from e
                else: 
                    if type(invar_filters) is list:
                        invar_filters = ops[compare](callset['/calldata/' + stat][:], value)
                    elif type(var_filters) is not list:
                        invar_filters = np.logical_and(invar_filters, ops[compare](callset['/calldata/' + stat][:], value))
            
            # create a mask for invariant sites by inverting the snp filter
            # join that to the invariant sites filter
            
            invar_filters = np.logical_and(invar_filters, np.invert(snp_mask))
    
            # join the variant and invariant filter masks (logical OR)
            filters = np.logical_or(invar_filters, var_filters)
                   
        # applying the filter to the data
        # all the filters are in a boolean array ('filters' above) 
    
        # first, recode the gt matrix as a Dask array (saves memory) -> packed
        # create a packed genotype array 
        # this is a array with dims snps x samples
        # genotypes are represented by single byte codes 
        # critically, as the same dims as the filters array below
        
        gt_array = allel.GenotypeArray(allel.GenotypeDaskArray(callset['/calldata/GT'])).to_packed()
    
        # apply filters
        # only if not bypassing filtration
        if args.bypass_filtration=='no':
            # set all genotypes that fail filters (the inversion of the array) 
            # to 'missing', 239 = -1 (i.e. missing) for packed arrays
            gt_array[np.invert(filters)] = 239
        
        # convert the packed array back to a GenotypeArray
        gt_array = allel.GenotypeArray.from_packed(gt_array)
    
        # build the position array
        pos_array = allel.SortedIndex(callset['/variants/POS'])
        
        # a mask for snps and invariant sites
        snp_invar_mask = np.logical_or(np.logical_and(callset['/variants/is_snp'][:] == 1, callset['/variants/numalt'][:] == 1), callset['/variants/numalt'][:] == 0)
        
        # remove rows that are NOT snps or invariant sites from the genotype array
        gt_array = np.delete(gt_array, np.where(np.invert(snp_invar_mask)), axis=0)
        gt_array = allel.GenotypeArray(gt_array)
    
        # select rows that ARE snps or invariant sites in the position array
        pos_array = pos_array[snp_invar_mask]
        
        # filter output for degbug
        # samples = callset['/samples'][:]
        # print(pos_array)
        # print(gt_array.shape)
        # print(gt_array)
        # np.savetxt("testing/filters/" + chromosome + '_snp_mask.csv', snp_mask, delimiter=',', fmt='%1d')
        # np.savetxt("testing/filters/" + chromosome + '_invar_mask.csv', np.invert(snp_mask), delimiter=',', fmt='%1d')
        # np.savetxt("testing/filters/" + chromosome + '_var_filter.csv', var_filters, delimiter=',', fmt='%1d')
        # np.savetxt("testing/filters/" + chromosome + '_invar_filter.csv', var_filters, delimiter=',', fmt='%1d')
        # np.savetxt("testing/filters/" + chromosome + '_complete_filter.csv', filters, delimiter=',', fmt='%1d')
        # np.savetxt("testing/filters/" + chromosome + '_filtered_genotypes.csv', gt_array.to_gt().decode('utf-8'), delimiter=',', fmt='%1s')
        # np.savetxt("testing/filters/" + chromosome + '_complete_positions.csv', pos_array, delimiter=',', fmt='%1d')
        # np.savetxt("testing/filters/" + chromosome + '_complete_samples.csv', samples, delimiter=',', fmt='%1s')
        
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
            interval_start = min(pos_array)
        else:
            interval_start = int(args.interval_start)
        
        try:     
            if (interval_start > interval_end):
                raise ValueError()        
        except ValueError as e:
            raise Exception("[pixy] ERROR: The specified interval start ("+str(interval_start)+") exceeds the interval end ("+str(interval_end)+")") from e
    
        # catch misspecified intervals
        # TBD: harmonize this with the new interval method for the zarr array
        if (interval_end > max(pos_array)):
            print("[pixy] WARNING: The specified interval end ("+str(interval_end)+") exceeds the last position of the chromosome and has been substituted with " + str(max(pos_array)))
            interval_end = max(pos_array)
            
        if (interval_start < min(pos_array)):
            print("[pixy] WARNING: The specified interval start ("+str(interval_start)+") begins before the first position of the chromosome and has been substituted with " + str(min(pos_array)))
            interval_start = min(pos_array)
    
        if ((interval_end - interval_start + 1) < window_size):
            print("[pixy] WARNING: The requested interval or total number of sites in the VCF ("+str(interval_start)+"-"+str(interval_end)+") is smaller than the requested window size ("+str(window_size)+")")
    
        # PI:
        # AVERAGE NUCLEOTIDE VARIATION WITHIN POPULATIONS
    
        # Compute pi over a chosen interval and window size
    
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
                    
                    # if the window has no sites, assign all NAs,
                    # otherwise calculate pi
                    if len(pos_array[(pos_array > window_pos_1 ) & (pos_array <window_pos_2)]) == 0:
                        avg_pi, total_diffs, total_comps, total_missing, no_sites = "NA", "NA", "NA", "NA", 0
                    else:
                    
                        # pull out the genotypes for the window
                        loc_region = pos_array.locate_range(window_pos_1, window_pos_2)
                        gt_region1 = gt_array[loc_region]
                        no_sites = len(gt_region1)
    
                        # subset the window for the individuals in each population 
                        gt_pop = gt_region1.take(popindices[pop], axis=1)
                        avg_pi, total_diffs, total_comps, total_missing = tallyRegion(gt_pop)
                        
                    outfile.write(str(pop) + "\t" + str(chromosome) + "\t" + str(window_pos_1) + "\t" + str(window_pos_2) + "\t" + str(avg_pi) + "\t" + str(no_sites) + "\t" + str(total_diffs) + "\t" + str(total_comps) + "\t" + str(total_missing) + "\n")
                    window_pos_2 += window_size
                    
                    if window_pos_2 > interval_end:
                        window_pos_2 = interval_end
    
                # close output file and print complete message
            outfile.close()
    
            print("[pixy] Pi calculations for chromosome " + chromosome + " complete and written to " + args.outfile_prefix + "_pi.txt")
    
    
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
                window_pos_2 = (interval_start + window_size)-1
    
                # perform the dxy calculation for all windows in the range
                for window_pos_1 in range (interval_start, interval_end, window_size):
                    
                    if len(pos_array[(pos_array > window_pos_1 ) & (pos_array <window_pos_2)]) == 0:
                        avg_dxy, total_diffs, total_comps, total_missing, no_sites = "NA", "NA", "NA", "NA", 0
                    else:
                        loc_region = pos_array.locate_range(window_pos_1, window_pos_2)
                        gt_region1 = gt_array[loc_region]
                        no_sites = len(gt_region1)
                        
                        # use the popGTs dictionary to keep track of this region's GTs for each population
                        popGTs={}
                        for name in pop_pair:
                            gt_pop = gt_region1.take(popindices[name], axis=1)
                            popGTs[name] = gt_pop
    
                        pop1_gt_region1 = popGTs[pop1]
                        pop2_gt_region1 = popGTs[pop2]
                        avg_dxy, total_diffs, total_comps, total_missing = dxyTallyRegion(pop1_gt_region1, pop2_gt_region1)
                        
                    outfile.write(str(pop1) + "\t" + str(pop2) + "\t" + str(chromosome) + "\t" + str(window_pos_1) + "\t" + str(window_pos_2) + "\t" + str(avg_dxy) + "\t" + str(no_sites) + "\t" + str(total_diffs) + "\t" + str(total_comps) + "\t" + str(total_missing) + "\n")
    
                    window_pos_2 += window_size
                    
                    if window_pos_2 > interval_end:
                        window_pos_2 = interval_end
    
            outfile.close()
            print("[pixy] Dxy calculations chromosome " + chromosome + " complete and written to " + args.outfile_prefix + "_dxy.txt")
    
        # FST:
        # WEIR AND COCKERHAMS FST
        # This is just a plain wrapper for the scikit-allel fst function
    
        if (args.populations is not None) and ('fst' in args.stats):
            
            # open the fst output file for writing
            outfile = open(fst_file, 'a')
            
            # determine all the possible population pairings
            pop_names=list(popindices.keys())
            fst_pop_list = list(combinations(pop_names, 2))
                      
            #calculate maf 
            allele_counts = gt_array.count_alleles()
            allele_freqs = allele_counts.to_frequencies()
            maf_array = allele_freqs[:,1] > args.fst_maf_filter
            
            # apply the maf filter to the genotype array]
            gt_array_fst = gt_array[maf_array]
            gt_array_fst = allel.GenotypeArray(gt_array_fst)
    
            # apply the maf filter to the position array
            pos_array_fst = pos_array[maf_array]
            
            # debug output
            # np.savetxt("testing/filters/" + chromosome + '_maf_array.csv', maf_array, delimiter=',', fmt='%1d')
            # np.savetxt("testing/filters/" + chromosome + '_fst_genotypes.csv', gt_array_fst.to_gt().decode('utf-8'), delimiter=',', fmt='%1s')
            # np.savetxt("testing/filters/" + chromosome + '_fst_positions.csv', pos_array_fst, delimiter=',', fmt='%1d')
    
            # for each pair, compute fst
            for pop_pair in fst_pop_list:
    
                # the indices for the individuals in each population
                fst_pop_indicies=[popindices[pop_pair[0]].tolist(), popindices[pop_pair[1]].tolist()]
                
                # compute FST
                # windowed_weir_cockerham_fst seems to generate (spurious?) warnings about div/0, so suppressing warnings
                # (this assumes that the scikit-allel function is working as intended)
                np.seterr(divide='ignore', invalid='ignore')
     
                a,b,c=allel.windowed_weir_cockerham_fst(pos_array_fst, gt_array_fst, subpops=fst_pop_indicies, size=args.window_size, start=interval_start, stop=interval_end)
    
                for fst,wind,snps in zip(a, b, c):
                    outfile.write(str(pop_pair[0]) + "\t" + str(pop_pair[1]) + "\t" + str(chromosome) + "\t" + str(wind[0]) + "\t" + str(wind[1]) + "\t" + str(fst) + "\t" + str(snps) +"\n")
            outfile.close()
            print("[pixy] Fst calculations chromosome " + chromosome + " complete and written to " + args.outfile_prefix + "_fst.txt")
    
    
    
    
    print("\n[pixy] All calculations complete at " + time.strftime("%H:%M:%S", time.localtime(start_time)))
    end_time = (time.time() - start_time)
    print("[pixy] Time elapsed: " + time.strftime("%H:%M:%S", time.gmtime(end_time)))
    
    
    
    
    
    
if __name__ == "__main__":
    main()