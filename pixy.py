#!/usr/bin/env python
# coding: utf-8

# In[1]:


import allel
import zarr
import numcodecs
import numpy as np
import sys
import matplotlib.pyplot as plt
from itertools import combinations
from collections import Counter
import argparse, sys


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

# https://stackoverflow.com/questions/40001892/reading-named-command-arguments
# example:


parser = argparse.ArgumentParser(description='pixy: senisbly calculate pi and/or dxy from a VCF containing variant and invariant sites')

parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('--stats', choices=['pi', 'dxy', 'pi_dxy'], help='Which stats to to calulate from the VCF (pi, dxy, or both)', required=True)
parser.add_argument('--vcf', type=str, nargs='?', help='Path to the input VCF', required=True)
parser.add_argument('--zarr', type=str, nargs='?', help='Folder in which to build the Zarr database', required=True)
parser.add_argument('--populations', type=str, nargs='?', help='Path to the populations file')
parser.add_argument('--window_size', type=int, nargs='?', help='Window size in base pairs over which to calculate pi/dxy')
parser.add_argument('--snp_filter_expression', type=str, nargs='?', help='A comma separated list of filters (e.g. DP>=10,GQ>=20) to apply to SNPs', required=True)
parser.add_argument('--invariant_filter_expression', type=str, nargs='?', help='A comma separated list of filters (e.g. DP>=10,RGQ>=20) to apply to invariant sites', required=True)
parser.add_argument('--out', type=str, nargs='?', help='Path to the output file')

parser.print_help()

# pull out varizbles from the parsed args
args = parser.parse_args('--stats pi_dxy --vcf test.vcf --zarr test/path/ --populations path/to/popfile.txt --snp_filter_expression DP>=10,GQ>=20,FS>2 --invariant_filter_expression DP>=10,RGQ>=20 --out pi_out.txt'.split())
# sim_args = parser.parse_args()

if (args.populations is None) and ((args.stats == 'dxy' or args.stats == 'pi_dxy')):
    parser.error("--stats dxy and --stats pi_dxy requires --populations path/to/popfile.txt")


# In[4]:


print(args)

print(args.invariant_filter_expression)


# In[ ]:


# validating inputs

# STEP 1 checking the vcf:
# - check if sci-kit allele can do this
# - check for contig info
# - alternatively use position vector
# - check for invariant sites and throw and error if they dont exist

# STEP 2 validate the population file
# - format is IND POP
# - throw an error if individuals are missing from VCF


# In[ ]:


# filtering the data 

# - extracting the filtration values from the filter strings
# - applying the filter to the whole dataset
# - also filter out biallelic snps
# - check to see there is any data left
# - print warning for low data* think more about this


# In[ ]:


# calculate pi

# for loop wrapper for pi calcuation
# output looks like:

# chromosome pos1 pos2 no_differences no_comparisons no_missing pi_1 pi_2 dxy
# chromosome pos1 pos2 no_differences no_comparisons no_missing pi
# chromosome pos1 pos2 no_differences no_comparisons no_missing dxy

# check if pi output makes sense
# - are any pis/dxys all zero?
# - check if # missing == (pos2 - pos1)
# - check if everything was either compared or missing
# - write out summary of program parameters file* think about how this should work


# In[ ]:


# convert this notebook to a .py script
get_ipython().system('jupyter nbconvert --to=python pixy.ipynb')

