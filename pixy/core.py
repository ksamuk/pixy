import allel
import numcodecs
import numpy as np
import sys
import os
import subprocess
import re
import pandas
import operator
import warnings
import time
import pathlib
import argparse
import shutil
import uuid

import pixy.calc

from multiprocess import Pool
from multiprocess import Queue
import multiprocess as mp

from scipy import special
from itertools import combinations
from collections import Counter

# function for re-aggregating subwindows from the temporary output
def aggregate_output(df_for_stat, stat, chromosome, window_size, fst_type):
    outsorted = df_for_stat.sort_values([4]) #sort by position
    interval_start = df_for_stat[4].min()
    interval_end = df_for_stat[4].max()
    # create a subset of the window list specific to this chunk
    bins =  np.arange(interval_start-1, interval_end + window_size, window_size)
    #bin values into discrete intervals, both left and right inclusive
    assignments, edges = pandas.cut(outsorted[4], bins=bins, labels=False, retbins=True, include_lowest = True)
    outsorted['label'] = assignments
    outsorted['window_pos_1'] = (edges[assignments]+1)
    outsorted['window_pos_2'] = (edges[assignments+1])
 
    #group by population, window
    if stat == 'pi' or stat == 'watterson_theta' or stat == 'tajima_d': #pi, Watterson's theta, and Tajima's D only has one population field
    outsorted = outsorted.groupby([1,'window_pos_1','window_pos_2'], as_index=False, dropna=False).agg({7:'sum',8:'sum',9:'sum',10:'sum'}).reset_index()
    elif stat == 'dxy' or stat == 'fst': #dxy and fst have 2 population fields
        outsorted = outsorted.groupby([1,2,'window_pos_1','window_pos_2'], as_index=False, dropna=False).agg({7:'sum',8:'sum',9:'sum',10:'sum'}).reset_index()   
        
    if stat == 'pi' or stat == 'dxy' or stat == 'watterson_theta' or stat == 'tajima_d':
            outsorted[stat] = outsorted[8]/outsorted[9]
    elif stat == 'fst':
        if fst_type == 'wc':
            outsorted[stat] = outsorted[8]/(outsorted[8] + outsorted[9] + outsorted[10])
        elif fst_type == 'hudson':
            # 'a' is the numerator of hudson and 'b' is the denominator
            # (there is no 'c')
            outsorted[stat] = outsorted[8]/(outsorted[9])
    
    outsorted[stat].fillna('NA', inplace=True)
    outsorted["chromosome"]=chromosome
    
    #reorder columns
    if stat == 'pi' or stat == 'watterson_theta' or stat == 'tajima_d': #pi, Watterson's theta and Tajima's D only have one population field
        outsorted = outsorted[[1, 'chromosome', 'window_pos_1','window_pos_2',stat, 7, 8, 9, 10]] 
    else: #dxy and fst have 2 population fields
        outsorted = outsorted[[1,2,'chromosome', 'window_pos_1','window_pos_2',stat, 7, 8, 9, 10]]

    # make sure sites, comparisons, missing get written as floats
    if stat == 'pi' or stat == 'dxy' or stat == 'watterson_theta' or stat == 'tajima_d':
        cols = [7,8,9,10]
    elif stat == 'fst':
        cols = [7]    
    outsorted[cols] = outsorted[cols].astype('float64')
    return outsorted

# function for breaking down large windows into chunks
def assign_subwindows_to_windows(window_pre_list, chunk_size):
    
    # build list of subwindows  
    window_lst = []
    for i in range(len(window_pre_list)):
        original_start = window_pre_list[i][0]
        original_stop = window_pre_list[i][1]
        subwindow_pos_1_list = [*range(original_start, original_stop, chunk_size)]
        subwindow_pos_2_list = [(item + chunk_size - 1) for item in subwindow_pos_1_list]
        
        #end of last window is always the end of the original window
        subwindow_pos_2_list[-1] = original_stop 
        
        #label which window # each subwindow originated from
        greater_window = [i] * len(subwindow_pos_2_list) 
        
        sub_windows = [list(a) for a in zip(subwindow_pos_1_list, subwindow_pos_2_list)]
        window_lst.extend(sub_windows)
    
    return window_lst

# function for assinging windows to larger chunks
# goal here to reduce I/O bottlenecks by only doing a single read of 
# the VCF for a group of small windows
def assign_windows_to_chunks(window_pre_list, chunk_size):
    
    # the start and end positions of each window
    window_pos_1_list = [item[0] for item in window_pre_list]
    window_pos_2_list = [item[1] for item in window_pre_list]
    
    # assign starts and ends of windows to chunks using np.floor
    chunk_list_1 = [np.floor(x / chunk_size) for x in window_pos_1_list]
    chunk_list_2 = [np.floor(x / chunk_size) for x in window_pos_2_list]
    
    # bind the lists back together
    window_lst = [list(a) for a in zip(window_pos_1_list, window_pos_2_list, chunk_list_1, chunk_list_2)]

    # nudge windows that overlap two chunks into the first chunk
    for i in range(len(window_lst)):
        if not window_lst[i][2] == window_lst[i][3]:
            window_lst[i][3] = window_lst[i][2]
    
    # remove chunk placeholders
    for i in range(len(window_lst)):
        del window_lst[i][3]

    return window_lst

# function for assinging sites to larger chunks
def assign_sites_to_chunks(sites_pre_list, chunk_size):
    
    # assign sites to chunks using np.floor
    chunk_list = [np.floor(x / (chunk_size + 1)) for x in sites_pre_list]
    
    # bind the lists back together
    sites_list = [list(a) for a in zip(sites_pre_list, chunk_list)]

    return sites_list

# function for masking non-target sites in a genotype array
# used in conjuctions with the --sites_list command line option
def mask_non_target_sites(gt_array, pos_array, sites_list_chunk):

    # get the indexes of sites that are NOT the target sites
    # (these will be masked with missing rows to remove them from the calculations)
    masked_sites = set(pos_array) ^ set(sites_list_chunk) 
    masked_sites = sorted(list(set(pos_array).intersection(masked_sites)))
    
    gt_mask_indexes = list(np.flatnonzero(pos_array.locate_keys(masked_sites)))
    
    # a missing row of data to use as a mask
    missing_row = [[-1] * gt_array.ploidy] * gt_array.n_samples
    
    # apply the mask to all non-target sites
    for pos in gt_mask_indexes:
        gt_array[pos,:] = missing_row 
        
    return gt_array 


# function for reading in a genotype matrix from a VCF file
# also filters out all but biallelic SNPs and invariant sites
# returns a genotype matrix, and array of genomic coordinates and
# a logical indicating whether the array(s) are empty 
def read_and_filter_genotypes(args, chromosome, window_pos_1, window_pos_2, sites_list_chunk):
    
    # a string representation of the target region of the current window
    window_region = chromosome + ":" + str(window_pos_1) + "-" + str(window_pos_2)
    
    # read in data from the source VCF for the current window
    callset = allel.read_vcf(args.vcf, region = window_region, fields = ['CHROM', 'POS', 'calldata/GT', 'variants/is_snp', 'variants/numalt'])
    
    # keep track of whether the callset was empty (no sites for this range in the VCF)
    # used by compute_summary_stats to add info about completely missing sites
    if callset is None:
        callset_is_none = True
        gt_array = None
        pos_array = None
        
    else:
        # if the callset is NOT empty (None), continue with pipeline
        callset_is_none = False

        # convert to a genotype array object
        gt_array = allel.GenotypeArray(allel.GenotypeDaskArray(callset['calldata/GT']))

        # build an array of positions for the region
        pos_array = allel.SortedIndex(callset['variants/POS'])

        # create a mask for biallelic snps and invariant sites
        snp_invar_mask = np.logical_or(np.logical_and(callset['variants/is_snp'][:] == 1, callset['variants/numalt'][:] == 1), callset['variants/numalt'][:] == 0)

        # remove rows that are NOT snps or invariant sites from the genotype array
        gt_array = np.delete(gt_array, np.where(np.invert(snp_invar_mask)), axis=0)
        gt_array = allel.GenotypeArray(gt_array)

        # select rows that ARE snps or invariant sites in the position array
        pos_array = pos_array[snp_invar_mask]
        
        # if a list of target sites was specified, mask out all non-target sites
        if sites_list_chunk is not None:
            gt_array = mask_non_target_sites(gt_array, pos_array, sites_list_chunk)

        # extra 'none' check to catch cases where every site was removed by the mask
        if len(gt_array) == 0:
            callset_is_none = True
            gt_array = None
            pos_array = None

    return callset_is_none, gt_array, pos_array

# main pixy function for computing summary stats over a list of windows (for one chunk)

def compute_summary_stats(args, popnames, popindices, temp_file, chromosome, chunk_pos_1, chunk_pos_2, window_list_chunk, q, sites_list_chunk, aggregate, window_size):
    
    # read in the genotype data for the chunk
    callset_is_none, gt_array, pos_array = read_and_filter_genotypes(args, chromosome, chunk_pos_1, chunk_pos_2, sites_list_chunk)


    # if computing FST, pre-compute a filtered array of variants (only)
    if 'fst' in args.stats:
        
        if (not callset_is_none) and (args.populations is not None) and (len(gt_array) != 0):

            # compute allel freqs 
            allele_counts = gt_array.count_alleles()
            allele_freqs = allele_counts.to_frequencies()

            # remove invariant/polyallelic sites
            variants_array = [len(x) == 2 and x[0] < 1 for x in allele_freqs]
            
            # filter gt and position arrays for biallelic variable sites
            gt_array_fst = gt_array[variants_array]
            pos_array_fst = pos_array[variants_array]

        else:
            pos_array_fst = None
            gt_array_fst = None
            
            
        
        # if obtaining per-site esimates,
        # compute the FST values for the whole chunk
        # instead of looping over subwindows (below)
        
        if window_size == 1:
            
             # determine all the possible population pairings
            pop_names = list(popindices.keys())
            fst_pop_list = list(combinations(pop_names, 2))

            # for each pair, compute fst using the filtered gt_array
            for pop_pair in fst_pop_list:

                # the indices for the individuals in each population
                fst_pop_indicies = [popindices[pop_pair[0]].tolist(), popindices[pop_pair[1]].tolist()]

                # compute FST
                # windowed_weir_cockerham_fst seems to generate (spurious?) warnings about div/0, so suppressing warnings
                # (this assumes that the scikit-allel function is working as intended)
                np.seterr(divide='ignore', invalid='ignore')

                # if the genotype matrix is not empty, compute FST
                # other wise return NA


                if(not callset_is_none and gt_array_fst is not None and len(gt_array_fst) > 0) :
                    fst = pixy.calc.calc_fst_persite(gt_array_fst, fst_pop_indicies, args.fst_type)
                    window_positions = list(zip(pos_array_fst, pos_array_fst))
                    n_snps = [1]*len(pos_array_fst)
                    
                    for fst, wind, snps in zip(fst, window_positions, n_snps):
                        # append trailing NAs so that pi/dxy/fst have the same # of columns
                        pixy_result = "fst" + "\t" + str(pop_pair[0]) + "\t" + str(pop_pair[1]) + "\t" + str(chromosome) + "\t" + str(wind[0]) + "\t" + str(wind[1]) + "\t" + str(fst) + "\t" + str(snps)+ "\tNA\tNA\tNA"
                        
                        if 'pixy_output' in locals():
                            pixy_output = pixy_output + "\n" + pixy_result
                        else:
                            pixy_output = pixy_result

    # loop over the windows within the chunk and compute summary stats
    for window_index in range(0, len(window_list_chunk)):

        window_pos_1 = window_list_chunk[window_index][0]
        window_pos_2 = window_list_chunk[window_index][1]
 
        

        if pos_array is None:
            window_is_empty = True
        elif (len(pos_array[(pos_array >= window_pos_1 ) & (pos_array <= window_pos_2)]) == 0):
            window_is_empty = True
        else:
            window_is_empty = False
            # pull out the genotypes for the window
            loc_region = pos_array.locate_range(window_pos_1, window_pos_2)
            gt_region = gt_array[loc_region]
            
            # double check that the region is not empty after subsetting
            try:
                loc_region
            except:
                loc_region = None
                
            try:
                gt_region
            except:
                gt_region = None

            if (len(gt_region) == 0 or (loc_region is None) or (gt_region is None)):
                window_is_empty = True
            

        # PI:
        # AVERAGE NUCLEOTIDE DIFFERENCES WITHIN POPULATIONS

        if (args.populations is not None) and ('pi' in args.stats):

            for pop in popnames:
                # if the window has no sites in the VCF, assign all NAs,
                # otherwise calculate pi
                if window_is_empty:
                    avg_pi, total_diffs, total_comps, total_missing, no_sites = "NA", "NA", "NA", "NA", 0
                else:
                    
                    # subset the window for the individuals in each population 
                    gt_pop = gt_region.take(popindices[pop], axis=1)
                    
                    # if the population specific window for this region is empty, report it as such
                    if (len(gt_pop) == 0):
                        avg_pi, total_diffs, total_comps, total_missing, no_sites = "NA", "NA", "NA", "NA", 0
                    
                    # otherise compute pi as normal
                    else:
                        # number of sites genotyped in the population
                        # not directly used in the calculation
                        no_sites = np.count_nonzero(np.sum(gt_pop.count_alleles(max_allele = 1), 1))
                        avg_pi, total_diffs, total_comps, total_missing = pixy.calc.calc_pi(gt_pop)

                # create a string of the pi results to write to file
                #klk added NA so that pi/dxy/fst have the same # of columns
                pixy_result = "pi" + "\t" + str(pop) + "\tNA\t" + str(chromosome) + "\t" + str(window_pos_1) + "\t" + str(window_pos_2) + "\t" + str(avg_pi) + "\t" + str(no_sites) + "\t" + str(total_diffs) + "\t" + str(total_comps) + "\t" + str(total_missing)

                # append the result to the multiline output string
                if 'pixy_output' in locals():
                    pixy_output = pixy_output + "\n" + pixy_result
                else:
                    pixy_output = pixy_result

        # DXY:
        # AVERAGE NUCLEOTIDE DIFFERENCES BETWEEN POPULATIONS

        if (args.populations is not None) and ('dxy' in args.stats):

            # create a list of all pairwise comparisons between populations in the popfile
            dxy_pop_list = list(combinations(popnames, 2))

            # interate over all population pairs and compute dxy
            for pop_pair in dxy_pop_list:
                pop1 = pop_pair[0]
                pop2 = pop_pair[1]

                if window_is_empty:
                    avg_dxy, total_diffs, total_comps, total_missing, no_sites = "NA", "NA", "NA", "NA", 0
                
                else:
                    # use the popGTs dictionary to keep track of this region's GTs for each population
                    popGTs={}
                    for name in pop_pair:
                        gt_pop = gt_region.take(popindices[name], axis=1)
                        popGTs[name] = gt_pop

                    pop1_gt_region = popGTs[pop1]
                    pop2_gt_region = popGTs[pop2]
                    
                    # if either of the two population specific windows for this region are empty, report it missing
                    if (len(pop1_gt_region) == 0 or len(pop2_gt_region) == 0):
                        avg_dxy, total_diffs, total_comps, total_missing, no_sites = "NA", "NA", "NA", "NA", 0
                        
                    # otherwise compute dxy as normal
                    else:
                        # for number of sites (not used in calculation), report the 
                        # number of sites that have at least one genotype in BOTH populations
                        pop1_sites = np.sum(pop1_gt_region.count_alleles(max_allele = 1), 1) > 0
                        pop2_sites = np.sum(pop2_gt_region.count_alleles(max_allele = 1), 1) > 0
                        no_sites = np.sum(np.logical_and(pop1_sites, pop2_sites))
                        avg_dxy, total_diffs, total_comps, total_missing = pixy.calc.calc_dxy(pop1_gt_region, pop2_gt_region)

                    # create a string of for the dxy results
                pixy_result = "dxy" + "\t" + str(pop1) + "\t" + str(pop2) + "\t" + str(chromosome) + "\t" + str(window_pos_1) + "\t" + str(window_pos_2) + "\t" + str(avg_dxy) + "\t" + str(no_sites) + "\t" + str(total_diffs) + "\t" + str(total_comps) + "\t" + str(total_missing)

                # append the result to the multiline output string
                if 'pixy_output' in locals():
                    pixy_output = pixy_output + "\n" + pixy_result
                else:
                    pixy_output = pixy_result


        # FST:
        # WEIR AND COCKERHAMS FST
        # This is just a loose wrapper around the scikit-allel fst function
        # TBD: explicit fst when data is completely missing

        if (args.populations is not None) and ('fst' in args.stats) and window_size != 1:
            
            # check for valid sites in the FST (variant only) data position array
            if pos_array_fst is not None:
                
                if np.logical_and(pos_array_fst >= window_pos_1, pos_array_fst <= window_pos_2).any():
                    # if there are valid sites, determine all the possible population pairings
                    pop_names = list(popindices.keys())
                    fst_pop_list = list(combinations(pop_names, 2))

                    # for each pair, compute fst using the filtered gt_array
                    for pop_pair in fst_pop_list:

                        # the indices for the individuals in each population
                        fst_pop_indicies = [popindices[pop_pair[0]].tolist(), popindices[pop_pair[1]].tolist()]

                        # compute FST
                        # windowed_weir_cockerham_fst seems to generate (spurious?) warnings about div/0, so suppressing warnings
                        # (this assumes that the scikit-allel function is working as intended)
                        np.seterr(divide='ignore', invalid='ignore')

                        # if the genotype matrix is not empty, compute FST
                        # other wise return NA


                        if(not callset_is_none and gt_array_fst is not None and len(gt_array_fst) > 0 and not window_is_empty) :

                            # compute an ad-hoc window size
                            fst_window_size = (window_pos_2 - window_pos_1)

                            # otherwise, compute FST using the scikit-allel window fst functions or the pixy "direct" method if aggregating
                            if not aggregate:

                                if args.fst_type  == "wc":
                                    fst, window_positions, n_snps = allel.windowed_weir_cockerham_fst(pos_array_fst, gt_array_fst, subpops = fst_pop_indicies, size = fst_window_size, start = window_pos_1, stop = window_pos_2)

                                if args.fst_type  == "hudson":
                                    ac1 = gt_array_fst.count_alleles(subpop = fst_pop_indicies[0])
                                    ac2 = gt_array_fst.count_alleles(subpop = fst_pop_indicies[1])
                                    fst, window_positions, n_snps = allel.windowed_hudson_fst(pos_array_fst, ac1, ac2, size = fst_window_size, start = window_pos_1, stop = window_pos_2)
                            else:
                                fst, a, b, c, n_snps = pixy.calc.calc_fst(gt_array_fst, fst_pop_indicies, args.fst_type)
                                window_positions = [[window_pos_1, window_pos_2]]

                        else:
                            # if there are no variable sites in the window, output NA/0 
                            # edit: i think it actually makes more sense to just omit these sites
                            if not aggregate:
                                fst, window_positions, n_snps = ["NA"],[[window_pos_1,window_pos_2]],[0]
                            else:
                                fst, window_positions, n_snps, a, b, c = "NA",[[window_pos_1,window_pos_2]], 0, "NA", "NA", "NA"

                        # create an output string for the FST results

                        #print(fst)
                        #print(window_positions)
                        #print(n_snps)
                        if aggregate:
                            pixy_result = "fst" + "\t" + str(pop_pair[0]) + "\t" + str(pop_pair[1]) + "\t" + str(chromosome) + "\t" + str(window_pos_1) + "\t" + str(window_pos_2) + "\t" + str(fst) + "\t" + str(n_snps)+ "\t" + str(a) + "\t" + str(b) +"\t" + str(c)

                        else:
                            for fst,wind,snps in zip(fst, window_positions, n_snps):
                                # append trailing NAs so that pi/dxy/fst have the same # of columns
                                pixy_result = "fst" + "\t" + str(pop_pair[0]) + "\t" + str(pop_pair[1]) + "\t" + str(chromosome) + "\t" + str(wind[0]) + "\t" + str(wind[1]) + "\t" + str(fst) + "\t" + str(snps)+ "\tNA\tNA\tNA"

                        # append the result to the multiline output string

                        if 'pixy_output' in locals():
                            pixy_output = pixy_output + "\n" + pixy_result

                        else:
                            pixy_output = pixy_result
                            
                #else:
                #    # if there are no variable sites in the window, output NA/0 
                #    # edit: i think it actually makes more sense to just omit these sites
                #    if not aggregate:
                #        fst, window_positions, n_snps = ["NA"],[[window_pos_1,window_pos_2]],[0]
                #        pixy_result = "fst" + "\t" + str(pop_pair[0]) + "\t" + str(pop_pair[1]) + "\t" + str(chromosome) + "\t" + str(wind[0]) + "\t" + str(wind[1]) + "\t" + str(fst) + "\t" + str(snps)+ "\tNA\tNA\tNA"
                #    else:
                #        fst, window_positions, n_snps, a, b, c = "NA",[[window_pos_1,window_pos_2]], 0, "NA", "NA", "NA"
                #        pixy_result = "fst" + "\t" + str(pop_pair[0]) + "\t" + str(pop_pair[1]) + "\t" + str(chromosome) + "\t" + str(window_pos_1) + "\t" + str(window_pos_2) + "\t" + str(fst) + "\t" + str(n_snps)+ "\t" + str(a) + "\t" + str(b) +"\t" + str(c)
                #    
                #    if 'pixy_output' in locals():
                #        pixy_output = pixy_output + "\n" + pixy_result
                #
                #    else:
                #        pixy_output = pixy_result

        # WATTERSON'S THETA:
        # GENETIC DIVERSITY CALCULATED FROM NUMBER OF SEGREGATING (VARIANT) SITES

        if (args.populations is not None) and ('watterson_theta' in args.stats):

            for pop in popnames:
                # if the window has no sites in the VCF, assign all NAs,
                # otherwise calculate Watterson's theta
                if window_is_empty:
                    avg_watterson_theta, watterson_theta, weighted_sites, no_sites, no_var_sites = "NA", "NA", "NA", 0, 0
                else:

                    # subset the window for the individuals in each population 
                    gt_pop = gt_region.take(popindices[pop], axis=1)

                    # if the population specific window for this region is empty, report it as such
                    if (len(gt_pop) == 0):
                        avg_watterson_theta, watterson_theta, weighted_sites, no_sites, no_var_sites = "NA", "NA", "NA", 0, 0

                    # otherise compute Watterson's theta as normal
                    else:
                        # alleles are counted, variant alleles are extracted
                        # number (no) of sites is count of sites with more than zero alleles
                        # also give number of variant sites
                        allele_counts = gt_pop.count_alleles(max_allele = 1)
                        variant_counts = allele_counts[allele_counts[:,1] != 0]
                        no_sites = np.count_nonzero(np.sum(allele_counts, 1))
                        no_var_sites = np.count_nonzero(np.sum(variant_counts, 1))
                        avg_watterson_theta, watterson_theta, weighted_sites = pixy.calc.calc_watterson_theta(gt_pop)                        
                # create a string of the Watterson's theta results to write to file
                #klk added NA so that pi/dxy/fst have the same # of columns, npb has kept this for watterson theta
                pixy_result = "watterson_theta" + "\t" + str(pop) + "\tNA\t" + str(chromosome) + "\t" + str(window_pos_1) + "\t" + str(window_pos_2) + "\t" + str(avg_watterson_theta) + "\t" + str(no_sites) + "\t" + str(no_var_sites) + "\t" + str(watterson_theta) + "\t" + str(weighted_sites)
                if 'pixy_output' in locals():
                    pixy_output = pixy_output + "\n" + pixy_result
                else:
                    pixy_output = pixy_result                  

        # TAJIMA_D:
        # NEUTRALITY TEST TAKING NORMALIZED DIFFERENCE OF WATTERSON'S THETA AND PI

        if (args.populations is not None) and ('tajima_d' in args.stats):

            for pop in popnames:
                # if the window has no sites in the VCF, assign all NAs,
                # otherwise calculate Tajima's D
                if window_is_empty:
                    tajima_d, no_sites, pi, watterson_theta, d_stdev = "NA", 0, "NA", "NA", "NA"
                else:

                    # subset the window for the individuals in each population 
                    gt_pop = gt_region.take(popindices[pop], axis=1)

                    # if the population specific window for this region is empty, report it as such
                    if (len(gt_pop) == 0):
                        tajima_d, no_sites, pi, watterson_theta, d_stdev = "NA", 0, "NA", "NA", "NA"

                    # otherwise compute Tajima's D as normal
                    else:
                        # number of sites genotyped in the population
                        # not directly used in the calculation
                        allele_counts = gt_pop.count_alleles(max_allele = 1)
                        no_sites = np.count_nonzero(np.sum(allele_counts, 1))
                        tajima_d, pi, watterson_theta, d_stdev = pixy.calc.calc_tajima_d(gt_pop)

                # create a string of the Tajima's D results to write to file
                #klk added NA so that pi/dxy/fst have the same # of columns, npb has kept this for Tajima's D
                pixy_result = "tajima_d" + "\t" + str(pop) + "\tNA\t" + str(chromosome) + "\t" + str(window_pos_1) + "\t" + str(window_pos_2) + "\t" + str(tajima_d) + "\t" + str(no_sites) + "\t" + str(pi) + "\t" + str(watterson_theta) + "\t" + str(d_stdev)
                if 'pixy_output' in locals():
                    pixy_output = pixy_output + "\n" + pixy_result
                else:
                    pixy_output = pixy_result
    # OUTPUT
    # if in mc mode, put the results in the writing queue
    # otherwise just write to file
    
    # ensure the output variable exists in some form
    
    if 'pixy_output' in locals(): 
            
        if(args.n_cores > 1):
            q.put(pixy_output)

        elif(args.n_cores == 1):
            outfile = open(temp_file, 'a')
            outfile.write(pixy_output +"\n")
            outfile.close()
                
# function for checking & validating command line arguments
# when problems are detected, throws an error/warning + message 

def check_and_validate_args(args):
    
    # CHECK FOR TABIX
    tabix_path = shutil.which("tabix")

    if tabix_path is None:
        raise Exception('[pixy] ERROR: tabix is not installed (or cannot be located in the path). Install tabix with "conda install -c bioconda htslib".') 

    if args.vcf is None:
        raise Exception('[pixy] ERROR: The --vcf argument is missing or incorrectly specified.') 
        
    if args.populations is None:
        raise Exception('[pixy] ERROR: The --populations argument is missing or incorrectly specified.') 
    
    # reformat file paths for compatibility 
    args.vcf = os.path.expanduser(args.vcf)
    args.populations = os.path.expanduser(args.populations)
    
    if args.output_folder != '':
        output_folder = args.output_folder + "/"
    else:
        output_folder = os.path.expanduser(os.getcwd() + "/")
        
    output_prefix = output_folder + args.output_prefix
    
    # get vcf header info
    vcf_headers = allel.read_vcf_headers(args.vcf)
    
    print("\n[pixy] Validating VCF and input parameters...")
    
    # CHECK OUTPUT FOLDER 
    print("[pixy] Checking write access...", end = '')
    check_message = "OK"
    
    # attempt to create the output folder
    if os.path.exists(output_folder) is not True:
        os.makedirs(output_folder)
    
    # check if output folder is writable
    #if not os.access(re.sub(r"[^\/]+$", "", args.outfile_prefix), os.W_OK):
    if not os.access(output_folder, os.W_OK):
        raise Exception('[pixy] ERROR: The output folder ' + output_folder + ' is not writable')
        
    # check if output_prefix is correctly specified
    if "/" in str(args.output_prefix) or "\\" in str(args.output_prefix):
        raise Exception('[pixy] ERROR: The output prefix \'' + str(args.output_prefix) + '\' contains slashes. Remove them and specify output folder structure with --output_folder if necessary.')
        

    # generate a name for a unique temp file for collecting output
    temp_file = output_folder + "pixy_tmpfile_" + str(uuid.uuid4().hex) + ".tmp"
    
    # check if temp file is writable
    with open(temp_file, 'w') as f:
        pass

    if check_message == "OK":
        print(check_message)

    # CHECK CPU CONFIGURATION
    print("[pixy] Checking CPU configuration...", end = '')
    check_message = "OK"
    
    if (args.n_cores > mp.cpu_count()):
        check_message = "WARNING"
        print(check_message)
        print('[pixy] WARNING: ' + str(args.n_cores) + ' CPU cores requested but only ' + str(mp.cpu_count()) + ' are available. Using '+ str(mp.cpu_count()) +'.')
        args.n_cores = mp.cpu_count()
    
    if check_message == "OK":
        print(check_message)

    # CHECK FOR EXISTANCE OF INPUT FILES

    if os.path.exists(args.vcf) is not True:
        raise Exception('[pixy] ERROR: The specified VCF ' + str(args.vcf) + ' does not exist') 
        
    if not re.search(".gz", args.vcf):
        raise Exception('[pixy] ERROR: The vcf is not compressed with bgzip (or has no .gz extension). To fix this, run "bgzip [filename].vcf" first (and then index with "tabix [filename].vcf.gz" if necessary)') 
    
    if not os.path.exists(args.vcf + ".tbi"):
        raise Exception('[pixy] ERROR: The vcf is not indexed with tabix. To fix this, run "tabix [filename].vcf.gz" first') 

    if os.path.exists(args.populations) is not True:
        raise Exception('[pixy] ERROR: The specified populations file ' + str(args.populations) + ' does not exist') 
    
    if args.bed_file is not None:
        args.bed_file = os.path.expanduser(args.bed_file)
        
        if os.path.exists(args.bed_file) is not True:
            raise Exception('[pixy] ERROR: The specified BED file ' + str(args.bed_file) + ' does not exist') 
            
    else:
        bed_df = []
        
    if args.sites_file is not None:
        args.sites_file = os.path.expanduser(args.sites_file)
        
        if os.path.exists(args.sites_file) is not True:
            raise Exception('[pixy] ERROR: The specified sites file ' + str(args.sites_file) + ' does not exist') 
    else:
        sites_df = []
        
    # VALIDATE THE VCF

    # check if the vcf contains any invariant sites
    # a very basic check: just looks for at least one invariant site in the alt field
    print("[pixy] Checking for invariant sites...", end = '')
    check_message = "OK"

    if args.bypass_invariant_check=='no':
        alt_list = subprocess.check_output("gunzip -c " + args.vcf + " | grep -v '#' | head -n 10000 | awk '{print $5}' | sort | uniq", shell=True).decode("utf-8").split()
        if "." not in alt_list:
            raise Exception('[pixy] ERROR: the provided VCF appears to contain no invariant sites (ALT = \".\"). This check can be bypassed via --bypass_invariant_check \'yes\'.') 
        if "." in alt_list and len(alt_list) == 1 :
            raise Exception('[pixy] ERROR: the provided VCF appears to contain no variable sites. It may have been filtered incorrectly, or otherwise corrupted.') 
    else:
        if not (len(args.stats) == 1 and (args.stats[0] == 'fst')):
            check_message = "WARNING"
            print(check_message)
            print("[pixy] EXTREME WARNING: --bypass_invariant_check is set to \'yes\'. Note that a lack of invariant sites will result in incorrect estimates.")

    if check_message == "OK":
        print(check_message)

    # check if requested chromosomes exist in vcf
    # parses the whole CHROM column (!)

    print("[pixy] Checking chromosome data...", end = '')
    
    # get the list of all chromosomes in the dataset
    chrom_all = subprocess.check_output("tabix -l " + args.vcf, shell=True).decode("utf-8").split()
    
    if args.chromosomes != 'all': 

        chrom_list = list(args.chromosomes.split(","))
        # pretabix method, can remove
        # chrom_all = subprocess.check_output("gunzip -c " + args.vcf + " | grep -v '#' | awk '{print $1}' | uniq", shell=True).decode("utf-8").split()
        chrom_all = subprocess.check_output("tabix -l " + args.vcf, shell=True).decode("utf-8").split()
        missing = list(set(chrom_list)-set(chrom_all))
        if len(missing) >0:
            raise Exception('[pixy] ERROR: the following chromosomes were specified but not occur in the VCF: ', missing) 

    else: #added this else statement (klk)
        chrom_list = subprocess.check_output("tabix -l " + args.vcf, shell=True).decode("utf-8").split()
        chrom_all = chrom_list

    print("OK")

    # INTERVALS
    # check if intervals are correctly specified
    # validate the BED file (if present)

    print("[pixy] Checking intervals/sites...", end = '')
    check_message = "OK"

    if args.bed_file is None:
        
        if args.window_size is None:
            raise Exception('[pixy] ERROR: In the absence of a BED file, a --window_size must be specified.') 

        if args.interval_start is None and args.interval_end is not None:
            raise Exception('[pixy] ERROR: When specifying an interval, both --interval_start and --interval_end are required.') 
            
        if args.interval_start is not None and args.interval_end is None:
            raise Exception('[pixy] ERROR: When specifying an interval, both --interval_start and --interval_end are required.') 

        if (args.interval_start is not None or args.interval_end is not None) and len(chrom_list) > 1:
            raise Exception('[pixy] ERROR: --interval_start and --interval_end are not valid when calculating over multiple chromosomes. Remove both arguments or specify a single chromosome.')  

        if (args.interval_start is not None and args.interval_end is not None) and ((int(args.interval_end) - int(args.interval_start)) <= int(args.window_size)):
            check_message = "WARNING"
            print('[pixy] WARNING: The specified interval ' + str(args.interval_start) + '-' + str(args.interval_end) + ' is smaller than the window size (' + str(args.window_size) + '). A single window will be returned.')
            

    else:
        if args.interval_start is not None or args.interval_end is not None or args.window_size is not None:
            check_message = "ERROR"
            print(check_message)
            raise Exception('[pixy] ERROR: --interval_start, --interval_end, and --window_size are not valid when a BED file of windows is provided.') 

        # read in the bed file and extract the chromosome column
        bed_df = pandas.read_csv(args.bed_file, sep='\t', usecols=[0,1,2], names=['chrom', 'pos1', 'pos2'])
        bed_df['chrom'] = bed_df['chrom'].astype(str)
        
        # force chromosomes to strings
        
        
        if bed_df.isnull().values.any():
            check_message = "ERROR"
            print(check_message)
            raise Exception('[pixy] ERROR: your bed file contains missing data, confirm all rows have three fields (chrom, pos1, pos2).') 

        if len(bed_df.columns) != 3:
            check_message = "ERROR"
            print(check_message)
            raise Exception('[pixy] ERROR: The bed file has the wrong number of columns (should be 3, is ' + str(len(bed_df.columns)) + ')') 
        
        else:
            bed_df.columns = ['chrom', 'chromStart', 'chromEnd']
            bed_chrom = list(bed_df['chrom'])
            missing = list(set(bed_chrom)-set(chrom_all))
            chrom_all = list(set(chrom_all) & set(bed_chrom))
            chrom_list = list(set(chrom_all) & set(bed_chrom))

        if len(missing) >0:
            check_message = "WARNING"
            print(check_message)
            print('[pixy] WARNING: the following chromosomes in the BED file do not occur in the VCF and will be ignored: ' + str(missing))
        
            
    if args.sites_file is not None:
        sites_df = pandas.read_csv(args.sites_file, sep='\t', usecols=[0,1], names=['chrom', 'pos'])
        sites_df['chrom'] = sites_df['chrom'].astype(str)
        
        if sites_df.isnull().values.any():
            check_message = "ERROR"
            print(check_message)
            raise Exception('[pixy] ERROR: your sites file contains missing data, confirm all rows have two fields (chrom, pos).') 
        
        if len(sites_df.columns) != 2:
            raise Exception('[pixy] ERROR: The sites file has the wrong number of columns (should be 2, is ' + str(len(sites_df.columns)) + ')') 
            
        else:
            sites_df.columns = ['CHROM', 'POS']
            chrom_sites = list(sites_df['CHROM'])
            missing = list(set(chrom_sites)-set(chrom_all))
            chrom_list = list(set(chrom_all) & set(chrom_sites))
            
        if len(missing) >0:
            check_message = "WARNING"
            print(check_message)
            print('[pixy] WARNING: the following chromosomes in the sites file do not occur in the VCF and will be ignored: ' + str(missing))

    if check_message == "OK":
        print(check_message)
    
    # SAMPLES
    # check if requested samples exist in vcf

    print("[pixy] Checking sample data...", end = '')

    # - parse + validate the population file
    # - format is IND POP (tab separated)
    # - throws an error if individuals are missing from VCF

    # read in the list of samples/populations
    poppanel = pandas.read_csv(args.populations, sep='\t', usecols=[0,1], names=['ID', 'Population'])
    poppanel['ID'] = poppanel['ID'].astype(str)
    
    # check for missing values
    
    if poppanel.isnull().values.any():
        check_message = "ERROR"
        print(check_message)
        raise Exception('[pixy] ERROR: your populations file contains missing data, confirm all samples have population IDs (and vice versa).') 

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
        check_message = "ERROR"
        print(check_message)
        raise Exception('[pixy] ERROR: the following samples are listed in the population file but not in the VCF: ', missing) from e
    else:   
        poppanel['callset_index'] = samples_callset_index

        # use the popindices dictionary to keep track of the indices for each population
        popindices={}
        popnames = poppanel.Population.unique()
        for name in popnames:
            popindices[name] = poppanel[poppanel.Population == name].callset_index.values
            
    if len(popnames) == 1 and ("fst" in args.stats or "dxy" in args.stats):
        check_message = "ERROR"
        print(check_message)
        raise Exception('[pixy] ERROR: calcuation of fst and/or dxy requires at least two populations to be defined in the population file.') 
    
    print("OK")
    print("[pixy] All initial checks past!")
    
    return popnames, popindices, chrom_list, IDs, temp_file, output_folder, output_prefix, bed_df, sites_df
