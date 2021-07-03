import pixy.core
import allel
import numpy as np

from scipy import special
from itertools import combinations
from collections import Counter

# vectorized functions for calculating pi and dxy 
# these are reimplementations of the original functions

# helper function for calculation of pi
# for the given site (row of the count table) count # of differences, # of comparisons, and # missing.
# uses number of haploid samples (n_haps) to determine missing data
def count_diff_comp_missing(row, n_haps):
    
    diffs = row[1] * row[0] 
    gts = row[1] + row[0]
    comps = int(special.comb(gts, 2))
    missing =  n_haps - gts
    return diffs, comps, missing

# function for vectorized calculation of pi from a pre-filtered scikit-allel genotype matrix
def calc_pi(gt_array):
     
    # counts of each of the two alleles at each site
    allele_counts = gt_array.count_alleles(max_allele = 1)
    
    # the number of (haploid) samples in the population
    n_haps = gt_array.n_samples * gt_array.ploidy
    
    # compute the number of differences, comparisons, and missing data at each site
    diff_comp_missing_matrix = np.apply_along_axis(count_diff_comp_missing, 1, allele_counts, n_haps) 
    
    # sum up the above quantities for totals for the region
    diff_comp_missing_sums = np.sum(diff_comp_missing_matrix, 0)
    
    # extract the component values
    total_diffs = diff_comp_missing_sums[0]
    total_comps = diff_comp_missing_sums[1]
    total_missing = diff_comp_missing_sums[2]
    
    # if there are valid data (comparisons between genotypes) at the site, compute average dxy
    # otherwise return NA
    if total_comps > 0:
        avg_pi = total_diffs/total_comps
    else:
        avg_pi = "NA"
        
    return(avg_pi, total_diffs, total_comps, total_missing)

# function for vectorized calculation of dxy from a pre-filtered scikit-allel genotype matrix
def calc_dxy(pop1_gt_array, pop2_gt_array):
    
    # the counts of each of the two alleles in each population at each site
    pop1_allele_counts = pop1_gt_array.count_alleles(max_allele = 1)
    pop2_allele_counts = pop2_gt_array.count_alleles(max_allele = 1)
    
    # the number of (haploid) samples in each population
    pop1_n_haps = pop1_gt_array.n_samples * pop1_gt_array.ploidy
    pop2_n_haps = pop2_gt_array.n_samples * pop2_gt_array.ploidy
    
    # the total number of differences between populations summed across all sites
    total_diffs = (pop1_allele_counts[:,0] * pop2_allele_counts[:,1]) + (pop1_allele_counts[:,1] * pop2_allele_counts[:,0])
    total_diffs = np.sum(total_diffs, 0)
    
    # the total number of pairwise comparisons between sites
    total_comps = (pop1_allele_counts[:,0] + pop1_allele_counts[:,1]) * (pop2_allele_counts[:,0] + pop2_allele_counts[:,1])
    total_comps = np.sum(total_comps, 0)

    # the total count of possible pairwise comparisons at all sites
    total_possible = (pop1_n_haps * pop2_n_haps) * len(pop1_allele_counts)

    # the amount of missing is possible comps - actual ('total') comps
    total_missing = total_possible - total_comps
    
    # if there are valid data (comparisons between genotypes) at the site, compute average dxy
    # otherwise return NA
    if total_comps > 0:
        avg_dxy = total_diffs/total_comps 
    else:
        avg_dxy = "NA"
        
    return(avg_dxy, total_diffs, total_comps, total_missing)


# function for obtaining fst AND variance components via scikit allel function
# (need variance components for proper aggregation)
# for single sites, this is the final FST calculation
# in aggregation mode, we just want a,b,c and n_sites for aggregating and fst
def calc_fst(gt_array_fst, fst_pop_indicies, fst_type):
    
    # compute basic (multisite) FST via scikit allel
    
    # WC 84
    if fst_type == "wc":
        a, b, c = allel.weir_cockerham_fst(gt_array_fst, subpops = fst_pop_indicies)
        
        # compute variance component sums
        a = np.nansum(a).tolist()
        b = np.nansum(b).tolist()
        c = np.nansum(c).tolist()
        n_sites = len(gt_array_fst)
    
        # compute fst
        if (a + b + c) > 0:
            fst = a / (a + b + c)
        else:
            fst = "NA"
    
        return(fst, a, b, c, n_sites)
    
    # Hudson 92
    if fst_type == "hudson":
        
        # following scikit allel docs
        # allel counts for each population
        ac1 = gt_array_fst.count_alleles(subpop = fst_pop_indicies[0])
        ac2 = gt_array_fst.count_alleles(subpop = fst_pop_indicies[1])
        
        #hudson fst has two components (numerator & denominator)
        num, den = allel.hudson_fst(ac1, ac2)
        c = 0 # for compatibility with aggregation code for WC 84
        
        # compute variance component sums
        num = np.nansum(num).tolist()
        den = np.nansum(den).tolist()
        n_sites = len(gt_array_fst)
        
        # compute fst
        if (num + den) > 0:
            fst = num / den
        else:
            fst = "NA"
        
        # same abc format as WC84, where 'a' is the numerator and 
        # 'b' is the demoninator, and 'c' is a zero placeholder
        return(fst, num, den, c, n_sites)

# simplified version of above to handle the case 
# of per-site estimates of FST over whole chunks

def calc_fst_persite(gt_array_fst, fst_pop_indicies, fst_type):
    
    # compute basic (multisite) FST via scikit allel
    
    # WC 84
    if fst_type == "wc":
        a, b, c = allel.weir_cockerham_fst(gt_array_fst, subpops = fst_pop_indicies)

        fst = (np.sum(a, axis=1) / (np.sum(a, axis=1) + np.sum(b, axis=1) + np.sum(c, axis=1)))
    
        return(fst)
    
    # Hudson 92
    elif fst_type == "hudson":
        
        # following scikit allel docs
        # allel counts for each population
        ac1 = gt_array_fst.count_alleles(subpop = fst_pop_indicies[0])
        ac2 = gt_array_fst.count_alleles(subpop = fst_pop_indicies[1])
        
        #hudson fst has two components (numerator & denominator)
        num, den = allel.hudson_fst(ac1, ac2)
        
        fst = num/den

        return(fst)
    
