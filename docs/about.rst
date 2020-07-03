************
About pixy
************

**Authors: Kieran Samuk and Katharine Korunes, Duke University**

pixy is a command-line tool for painlessly and correctly estimating average nucleotide diversity within (π) and between (dxy) populations
from a VCF. In particular, pixy facilitates the use of VCFs containing invariant (AKA monomorphic) sites, which are essential for the 
correct computation of π and dxy:

Population geneticists are often interested in quantifying nucleotide diversity within and nucleotide differences between populations. 
The two most common summary statistics for these quantities were described by Nei and Li (1979), who discuss summarizing variation in the
case of two populations (denoted 'x' and 'y'):

* π - Average nucleotide diversity within populations, also sometimes denoted πx and πy to indicate the population of interest.
* dxy - Average nucleotide difference between populations, also sometimes denoted πxy (pixy, get it?), to indicate that the statistic is a 
  comparison between populations x and y.

.. note::
    pixy is currently in active development and is not ready for general use. 
    Our preprint has not yet undergone full peer-review but is available and citable here:
    https://www.biorxiv.org/content/10.1101/2020.06.27.175091v1

Common pitfalls in calculating π and dxy
#####################

Both π and dxy measure the average number of differences between sequences per nucleotide not per SNP. As such, one must include 
monomorphic/invariant sites when tallying differences between sequences. Prior to the genomic era, such sites were almost always explicitly 
included in datasets because sequence data was in FASTA format (e.g. Sanger reads). However, most modern genomics tools encode variants as 
VCFs which by design often omit invariant sites. With variants-only VCFs, there is often no way to distinguish missing sites from invariant 
sites. Further, when one does include invariant sites in a VCF, it generally results in very large files that are difficult to manipulate 
with standard tools.

Better π and dxy calculation with pixy
#####################

pixy provides the following solutions to problems inherent in computing π and dxy from a VCF:

1. Fast and efficient handing of invariant sites VCFs via conversion to on-disk chunked databases (Zarr format).
2. Standardized individual-level filtration of variant and invariant genotypes.
3. Computation of π and dxy for arbitrary numbers of populations
4. Computes all statistics in arbitrarily sized windows, and output contains all raw information for all computations (e.g. numerators and denominators).
5. The majority of this is made possible by extensive use of the existing data structures and functions found in the brilliant python library scikit-allel
