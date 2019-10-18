#!/usr/bin/env python
# coding: utf-8

# In[3]:


import allel
import zarr
import numcodecs
import numpy as np
import sys
import os
import re
import operator
import pandas
from itertools import combinations
from collections import Counter
from tqdm import tqdm
import argparse


# In[4]:


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

help_image = "██████╗ ██╗██╗  ██╗██╗   ██╗\n" "██╔══██╗██║╚██╗██╔╝╚██╗ ██╔╝\n" "██████╔╝██║ ╚███╔╝  ╚████╔╝\n" "██╔═══╝ ██║ ██╔██╗   ╚██╔╝\n" "██║     ██║██╔╝ ██╗   ██║\n" "╚═╝     ╚═╝╚═╝  ╚═╝   ╚═╝\n" 

help_text = 'pixy: sensible estimates of pi and dxy from a VCF'


parser = argparse.ArgumentParser(description=help_image+help_text, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--version', action='version', version='%(prog)s version 1.0')
parser.add_argument('--stats', choices=['pi', 'dxy', 'pi_dxy'], help='Which stats to to calulate from the VCF (pi, dxy, or both)', required=True)
parser.add_argument('--vcf_path', type=str, nargs='?', help='Path to the input VCF', required=True)
parser.add_argument('--zarr_path', type=str, nargs='?', help='Folder in which to build the Zarr array', required=True)
parser.add_argument('--regenerate_zarr', choices=['yes', 'no'], help='Force regeneration of the Zarr array')
parser.add_argument('--populations', type=str, nargs='?', help='Path to the populations file', required=True)
parser.add_argument('--chromosome', type=str, nargs='?', help='Target chromosome (as annotated in the CHROM field)', required=True)
parser.add_argument('--window_size', type=int, nargs='?', help='Window size in base pairs over which to calculate pi/dxy')
parser.add_argument('--variant_filter_expression', type=str, nargs='?', help='A comma separated list of filters (e.g. DP>=10,GQ>=20) to apply to SNPs', required=True)
parser.add_argument('--invariant_filter_expression', type=str, nargs='?', help='A comma separated list of filters (e.g. DP>=10,RGQ>=20) to apply to invariant sites', required=True)
parser.add_argument('--outfile_prefix', type=str, nargs='?', help='Path and prefix for the output file, e.g. path/to/outfile')

#parser.print_help()

# test values
#args = parser.parse_args('--stats pi_dxy --vcf test.vcf --zarr test/path/ --populations data/vcf/ag1000/Ag1000_sampleIDs_popfile_3.txt --regenerate_zarr no --variant_filter_expression DP>=10,GQ>=20,RGQ>=20 --invariant_filter_expression DP>=10,RGQ>=20 --outfile_prefix output/pixy_out'.split())

args = parser.parse_args()


# In[3]:


#print(args)
#print(args.filter_expression)


# In[4]:


# validating inputs

# STEP 1 checking the vcf:
# - check if sci-kit allele can do this
# - check for contig info
# - alternatively use position vector
# - check for invariant sites and throw and error if they dont exist


# In[5]:


# Zarr array conversion

# perform the vcf to zarr conversion if the zarr array is missing, or regeneration has been requested

if os.path.exists(zarr_path) is not True:
    print("Zarr array does not exist, building...")
    allel.vcf_to_zarr(vcf_path, zarr_path, group=chromosome, fields='*', log=sys.stdout, overwrite=True)
elif 'regenerate_zarr' in args:
    if args.regenerate_zarr == 'yes':
        print("Regenerating Zarr array...")
        allel.vcf_to_zarr(vcf_path, zarr_path, group=chromosome', fields='*', log=sys.stdout, overwrite=True)

# inspect the structure of the zarr data
callset = zarr.open_group(zarr_path, mode='r')

#callset.tree(expand=True)


# In[6]:


# STEP 2 prase + validate the population file
# - format is IND POP (tab separated)
# - throws an error if individuals are missing from VCF

# read in the list of samples/populations
poppanel = pandas.read_csv(args.populations, sep='\t', usecols=[0,1], names=['ID', 'Population'])
poppanel.head()

# get a list of samples from the callset
samples = callset[chromosome + '/samples'][:]
samples_list = list(samples)
#print('VCF samples:', samples_list)

#make sure every indiv in the pop file is in the VCF callset
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
        #print(name)
        popindices[name] = poppanel[poppanel.Population == name].callset_index.values

    #print(popindices)


# In[7]:


# parse the filtration expression and build the boolean filter array

# define an operator dictionary for parsing the operator strings
ops = { "<": operator.lt, "<=": operator.le, ">": operator.gt, ">=": operator.ge, "==": operator.eq}

# determine the complete list of available calldata fields usable for filtration
calldata_fields = sorted(callset[chromosome + '/calldata/'].array_keys())

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
                GQ_filter = ops[compare](callset[chromosome + '/calldata/' + stat][:], value)
            elif stat == 'RGQ':
                RGQ_filter = ops[compare](callset[chromosome + '/calldata/' + stat][:], value)
            else:
                # this creates the initial filter array
                filters = ops[compare](callset[chromosome + '/calldata/' + stat][:], value)
        elif type(filters) is not list:
            if stat == 'GQ':
                GQ_filter = ops[compare](callset[chromosome + '/calldata/' + stat][:], value)
            elif stat == 'RGQ':
                RGQ_filter = ops[compare](callset[chromosome + '/calldata/' + stat][:], value)
            else:
                # this updates the filter array with additional filtration criteria
                filters = np.logical_and(filters, ops[compare](callset[chromosome + '/calldata/' + stat][:], value))

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


# In[8]:


# applying the filter to the data
# all the filters are in a boolean array ('filters') 

# TBD: check to see there is any data left
# TBD: print warning for low data* think more about this

#print("Applying filters...")

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

# remove sites with >1 alt allele from the filtered gt array
gt_array = np.delete(gt_array, np.where(callset[chromosome + '/variants/numalt'][:] > 1), axis=0)

# convert the packed array back to a GenotypeArray
gt_array = allel.GenotypeArray.from_packed(gt_array)

# build the position array
pos_array = allel.SortedIndex(callset[chromosome + '/variants/POS'])

# remove non-snps from the position array
pos_array = pos_array[callset[chromosome + '/variants/numalt'][:] < 2]


# In[9]:


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

def dxyCompareGTs(vec1, vec2):
    # use gts to select only sites with data. counting missing data as a check, but it's implicit in the length of (gts)
    #gts for each of the 2 populations being compared:
    gts1 = []
    gts2 = []
    missing = 0
    for x in vec1:
        if x in (0,1):
            gts1.append(x)
        else:
            missing += 1
    for y in vec2:
        if y in (0,1):
            gts2.append(y)
        else:
            missing += 1
     
    diffs = 0
    comps = 0
    length1 = len(gts1)
    length2 = len(gts2)
    for x in range(0,length1):
        i = gts1[x]
        for n in range(0,length2):
            j = gts2[n]
            comps += 1
            if i != j:
                diffs += 1    
        
    return(diffs,comps,missing)


# In[10]:


# PI:
# AVERAGE NUCLEOTIDE VARIATION WITHIN POPULATIONS

# Compute pi over a chosen window size

# calculate pi

# TBD:
# - are any pis/dxys all zero?
# - check if # missing == (pos2 - pos1)
# - check if everything was either compared or missing
# - write out summary of program parameters file* think about how this should work
# - if there are >1 populations, do this for each separately!

# total chromosome length
# chr_length = max(pos_array)

if (args.populations is not None) and ((args.stats == 'pi' or args.stats == 'pi_dxy')):

    chr_length = 1000000 # testing value

    # initialize the pi output file names

    for pop in popnames:

        # window size:
        window_size = 1000

        # initialize window_pos_2 
        window_pos_2 = window_size
        
        # create pi name via the prefix
        pi_file = str(args.outfile_prefix) + "_" + str(pop) +"_pi.txt"

        # remove any existing pop files
        if os.path.exists(pi_file):
            os.remove(pi_file)

        # open the dxy output file for writing
        outfile = open(pi_file, 'w')
        outfile.write("pop" + "\t" + "chromosome" + "\t" + "window_pos_1" + "\t" + "window_pos_2" + "\t" + "avg_pi" + "\t" + "no_sites" + "\t" + "count_diffs" + "\t" + "count_comparisons" + "\t" + "count_missing" + "\n")

        # loop over populations and windows, compute stats and write to file
        for window_pos_1 in tqdm(range (1, chr_length, window_size)):

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

    print("Pi calculations complete and written to " + args.outfile_prefix + "_[popname]_pi.txt")


# In[11]:


# DXY:
# AVERAGE NUCLEOTIDE VARIATION BETWEEN POPULATIONS

#come back to this later (parsing arguments to dictate pi/dxy behaviour)
#if (args.populations is not None) and ((args.stats == 'dxy' or args.stats == 'pi_dxy'))

# total chromosome length
# chr_length = max(pos_array)

if (args.populations is not None) and ((args.stats == 'dxy' or args.stats == 'pi_dxy')):

    chr_length = 10000 # testing value, total length of the chromosome

    # create a list of all pairwise comparisons between populations in the popfile
    dxy_pop_list = list(combinations(popnames, 2))

    # interate over all population pairs and compute dxy
    for pop_pair in dxy_pop_list:
        pop1 = pop_pair[0]
        pop2 = pop_pair[1]
        
        # window size:
        window_size = 1000

        # initialize window_pos_2 
        window_pos_2 = window_size
        
        # rename the dxy output file based on the prefix
        dxy_file = str(args.outfile_prefix) + "_" + str(pop1) + "_" + str(pop2) +"_dxy.txt"

        # remove any previous results
        if os.path.exists(dxy_file):
            os.remove(dxy_file)

        # open the dxy output file for writing
        outfile = open(dxy_file, 'w')
        outfile.write("pop1" + "\t" + "pop2" + "\t" + "chromosome" + "\t" + "window_pos_1" + "\t" + "window_pos_2" + "\t" + "avg_dxy" + "\t" + "no_sites" + "\t" + "count_diffs" + "\t" + "count_comparisons" + "\t" + "count_missing" + "\n")

        # perform the dxy calculation for all windows in the range
        for window_pos_1 in tqdm(range (1, chr_length, window_size)):
            loc_region = pos_array.locate_range(window_pos_1, window_pos_2)
            gt_region1 = gt_array[loc_region]
            # use the popGTs dictionary to keep track of this region's GTs for each population
            popGTs={}
            for name in pop_pair:
                #print(popindices[name])
                gt_pop = gt_region1.take(popindices[name], axis=1)
                popGTs[name] = gt_pop
            #print(popGTs)
            pop1_gt_region1 = popGTs[pop1]
            pop2_gt_region1 = popGTs[pop2]
            avg_dxy, total_diffs, total_comps, total_missing = dxyTallyRegion(pop1_gt_region1, pop2_gt_region1)
            outfile.write(str(pop1) + "\t" + str(pop2) + "\t" + str(chromosome) + "\t" + str(window_pos_1) + "\t" + str(window_pos_2) + "\t" + str(avg_dxy) + "\t" + str(len(gt_region1)) + "\t" + str(total_diffs) + "\t" + str(total_comps) + "\t" + str(total_missing) + "\n")

            #print("Region:", x , y, ",", "Region length:", len(pop1_gt_region1))
            #print("Average dxy:", avg_dxy)
            #print("Diffs, comps, missing:", total_diffs, total_comps, total_missing, "\n")
            window_pos_2 += window_size

        outfile.close()

    print("Dxy calculations complete and written to " + args.outfile_prefix + "_[pop1]_[pop2]_dxy.txt")


# In[1]:


# convert this notebook to a .py script
# get_ipython().system(u'jupyter nbconvert --to=python pixy.ipynb')


# In[ ]:




