#!/usr/bin/env python
# coding: utf-8

# In[1]:


import allel
import zarr
import numcodecs
import numpy as np
import sys
import os
import re
import operator
import matplotlib.pyplot as plt
from itertools import combinations
from collections import Counter
from tqdm import tqdm
import argparse


# In[2]:


# parse the command line arguements
# - validate the arguements 
# - throw errors if requiremabents are missing
# - validate the filter strings

#pi_dxy --pi --dxy\ # must include at least 1 of pi and/dxy 
#--vcf allsites.vcf.gz \ # required
#--zarr path/to/zarr \  # default to the vcf folder
#--populations popfile.txt \ # only required for dxy
#--window_size 10000 \ # not required, defaults to whole genome
#--snp_filter_expression "DP>=10, GQ>=20, FS>2" \ # required
#--monomorphic_filter_expression "DP>=10, RGQ>=20" \ # required
#--out pi_dxy_out.txt # default to vcf path + suffix

# initialize and add arguments to the argparser

parser = argparse.ArgumentParser(description='pixy: senisbly calculate pi and/or dxy from a VCF containing variant and invariant sites')

parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('--stats', choices=['pi', 'dxy', 'pi_dxy'], help='Which stats to to calulate from the VCF (pi, dxy, or both)', required=True)
parser.add_argument('--vcf', type=str, nargs='?', help='Path to the input VCF', required=True)
parser.add_argument('--zarr', type=str, nargs='?', help='Folder in which to build the Zarr database', required=True)
parser.add_argument('--populations', type=str, nargs='?', help='Path to the populations file')
parser.add_argument('--window_size', type=int, nargs='?', help='Window size in base pairs over which to calculate pi/dxy')
parser.add_argument('--filter_expression', type=str, nargs='?', help='A comma separated list of filters (e.g. DP>=10,GQ>=20) to apply to SNPs', required=True)
parser.add_argument('--invariant_filter_expression', type=str, nargs='?', help='A comma separated list of filters (e.g. DP>=10,RGQ>=20) to apply to invariant sites', required=True)
parser.add_argument('--outfile', type=str, nargs='?', help='Path to the output file')

#parser.print_help()

# pull out varizbles from the parsed args
args = parser.parse_args('--stats pi_dxy --vcf test.vcf --zarr test/path/ --populations path/to/popfile.txt --filter_expression DP>=10,GQ>=20,RGQ>=20 --invariant_filter_expression DP>=10,RGQ>=20 --outfile pixy_out.txt'.split())
# sim_args = parser.parse_args()

if (args.populations is None) and ((args.stats == 'dxy' or args.stats == 'pi_dxy')):
    parser.error("--stats dxy and --stats pi_dxy requires --populations path/to/popfile.txt")


# In[3]:


print(args)

print(args.filter_expression)


# In[4]:


# validating inputs

# STEP 1 checking the vcf:
# - check if sci-kit allele can do this
# - check for contig info
# - alternatively use position vector
# - check for invariant sites and throw and error if they dont exist

# STEP 2 validate the population file
# - format is IND POP
# - throw an error if individuals are missing from VCF


# In[5]:


#Test data:
chromosome = "chrX"
vcf_path = 'data/vcf/ag1000/chrX_36Ag_allsites.vcf.gz'
zarr_path = 'data/vcf/ag1000/chrX_36Ag_allsites.zarr'
#vcf_path = '/Users/Katharine Korunes/Documents/Dxy_test_data/chrX_36Ag_allsites.vcf.gz'
#zarr_path = '/Users/Katharine Korunes/Documents/Dxy_test_data/chrX_36Ag_allsites.zarr'
# inspect the zarr data
callset = zarr.open_group(zarr_path, mode='r')
#callset.tree(expand=True)


# In[14]:


# parse the filtration expression and build the boolean filter array

# define an operator dictionary for parsing the operator strings
ops = { "<": operator.lt, "<=": operator.le, ">": operator.gt, ">=": operator.ge, "==": operator.eq}

filters = []

# split the filtration expressions by commas
# parse out each component
# apply the filter
# TBD: handle cases where the specified filter isn't in the callset
for x in args.filter_expression.split(","):
    stat = re.sub("[^A-Za-z]+", "", x)
    value = int(re.sub("[^0-9]+", "", x))
    compare = re.sub("[A-Za-z0-9]+", "", x)
    
    if type(filters) is list:
        filters = ops[compare](callset[chromosome + '/calldata/' + stat][:], value)
    else:
        if stat == 'GQ':
            GQ_filter = ops[compare](callset[chromosome + '/calldata/' + stat][:], value)
        elif stat == 'RGQ':
            RGQ_filter = ops[compare](callset[chromosome + '/calldata/' + stat][:], value)
        else :
            filters = np.logical_and(filters, ops[compare](callset[chromosome + '/calldata/' + stat][:], value))

# check if GQ and RQG exist
# if they both exist, perform a logical OR and join them into the filter
# otherwise, perform a logical AND to join either one into the filter

try:
    GQ_filter
except NameError:
    GQ_exists = False
else:
    GQ_exists = True
    
try:
    RGQ_filter
except NameError:
    RGQ_exists = False
else:
    RGQ_exists = True
    
if GQ_exists & RGQ_exists:
    filters = np.logical_and(filters, np.logical_or(GQ_filter, RGQ_filter))
elif GQ_exists:
    filters = np.logical_and(filters, GQ_filter)
elif RGQ_exists:
    filters = np.logical_and(filters, RGQ_filter)
    
# finally, invert the whole array 
# this is for convenience/brevity in the next section

filters = np.invert(filters)


# In[15]:


# applying the filter to the data
# all the filters are in a boolean array ('filters') 

# TBD: check to see there is any data left
# TBD: print warning for low data* think more about this

# the genotype calls
# recode the gt matrix as a Dask array (saves memory)
gt_dask = allel.GenotypeDaskArray(callset[chromosome + '/calldata/GT'])

# create a packed genotype array 
# this is a array with dims snps x samples
# genotypes are represented by single byte codes 
# critically, as the same dims as the filters array below
gt_array = allel.GenotypeArray(gt_dask).to_packed()

# set all genotypes that fail filters to 'missing'
# 239 = -1 (i.e. missing) for packed arrays
gt_array[filters] = 239

# remove non snps from the filtered gt array
gt_array = np.delete(gt_array, np.where(callset[chromosome + '/variants/numalt'][:] > 1), axis=0)

# convert the packed array back to a GenotypeArray
gt_array = allel.GenotypeArray.from_packed(gt_array)

# build the position array
pos_array = allel.SortedIndex(callset[chromosome + '/variants/POS'])

# remove non-snps from the position array
pos_array = pos_array[callset[chromosome + '/variants/numalt'][:] < 2]


# In[18]:


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

#Out of sites with data, count differences. 
#Return the number of differences, the number of comparisons, and missing data count.
def compareGTs(vec):
    # use gts to select only sites with data. 
    # counting missing data as a check, but it's implicit in the length of (gts)
    gts = []
    missing = 0
    for x in vec:
        if x in (0,1):
            gts.append(x)
        else:
            missing += 1
    c = Counter(gts)
    #print(c)
    diffs = c[1]*c[0]
    comps = len(list(combinations(gts, 2)))        
    return(diffs,comps,missing)


# In[19]:


# Compute pi over a chosen window size

# calculate pi

# TBD:
# - are any pis/dxys all zero?
# - check if # missing == (pos2 - pos1)
# - check if everything was either compared or missing
# - write out summary of program parameters file* think about how this should work

# total chromosome length
# chr_length = max(pos_array)
chr_length = 10000 # testing value

# window size:
window_size = 1000

# initialize window_pos_2 
window_pos_2 = window_size

# open the output file for writing
os.remove(args.outfile)
outfile = open(args.outfile, 'w')
outfile.write("chromosome" + "\t" + "window_pos_1" + "\t" + "window_pos_2" + "\t" + "avg_pi" + "\t" + "no_sites" + "\t" + "total_diffs" + "\t" + "total_comps" + "\t" + "total_missing" + "\n")

# loop over windows, compute stats and write to file
for window_pos_1 in tqdm(range (1, chr_length, window_size)):
    loc_region = pos_array.locate_range(window_pos_1, window_pos_2)
    gt_region1 = gt_array[loc_region]
    avg_pi, total_diffs, total_comps, total_missing = tallyRegion(gt_region1)
    outfile.write(str(chromosome) + "\t" + str(window_pos_1) + "\t" + str(window_pos_2) + "\t" + str(avg_pi) + "\t" + str(len(gt_region1)) + "\t" + str(total_diffs) + "\t" + str(total_comps) + "\t" + str(total_missing) + "\n")
    window_pos_2 += window_size

# close output file and print complete message
outfile.close()
print("Calculations complete and written to " + args.outfile)


# In[20]:


# convert this notebook to a .py script
get_ipython().system('jupyter nbconvert --to=python pixy.ipynb')


# In[ ]:




