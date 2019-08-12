``pixy``
====================

A command-line tool for sensibly estimating average nucleotide diversity within (π) and between (d<sub>xy</sub>) populations from a VCF. In particular, pixy accomodates (and encourages) the use of VCFs containing invariant (AKA monomorphic) sites, which are **essential** for the correct computation of π and d<sub>xy</sub>.

***pixy is currently in active development and not intended for general use. The software will be fully described in a forthcoming publication.***

## Authors
Kieran Samuk and Katharine Korunes <p>
Duke University

## Installation
Currently: Clone from Github, with manual dependancy installation.<p>
Planned: Installation via Anaconda recipe.

## Background

Population geneticists are often interested in quantifying nucleotide diversity within and nucleotide differences between populations. The two most common summary statistics for these quantities were described by Nei and Li (1979), who discuss summarizing variation case of two populations (denoted 'x' and 'y'):

- π  - Average nucleotide diversity within populations, also sometimes denoted π<sub>x</sub> and π<sub>y</sub> to indicate the population of interest.
- d<sub>xy</sub> - Average nucleotide difference between populations, also sometimes denoted π<sub>xy</sub> (pixy, get it?), to indicate that the statistic is a comparison between populations x and y.

Both of these statistics use the same basic formula:

![](https://wikimedia.org/api/rest_v1/media/math/render/svg/be2956df9d2756a4f051f2516938d4831fcd3771)

x<sub>i</sub> and x<sub>j</sub> are the respective frequencies of the i<sub>th</sub>and j<sub>th</sub> sequences, π<sub>ij</sub> is the number of nucleotide differences per nucleotide site between the  i<sub>th</sub> and j<sub>th</sub> sequences, and n is the number of sequences in the sample. [(Source: Wikipedia)](https://en.wikipedia.org/wiki/Nucleotide_diversity)

In the case of π, all comparisons are made between sequences from *same* population, wherease for d<sub>xy</sub> all comparisons are made between sequences from two *different* populations.

## The problem

Both of π and d<sub>xy</sub> measure the average number of differences between sequences **per nucleotide not per SNP**. As such, one must include monomorphic/invariant sites when tallying differences between sequences. Prior to the genomic era,  such sites were almost always explicitly included in datasets because sequence data was in FASTA format (e.g. Sanger reads). However, modern genomics tools often encode variants as VCFs, which, by design, often omit invariant sites. As such, there is often no way to distinguish missing sites from invariant sites. When one does include invariant sites in a VCF, it generally results in very large files that are difficult to manipulate with standard tools. 

## The solution

pixy provides the following solutions to problems inherent in computing π and d<sub>xy</sub> from a VCF: 
 
1. Fast and efficient handing of invariant sites VCFs via conversion to on-disk chunked databases (Zarr format).
2. Standardized individual-level filtration of variant and invariant genotypes.
3. Computation of π and d<sub>xy</sub> for arbitrary numbers of populations 
4. Computes all statistics in arbitrarily sized windows, and output contains all raw information for all computations (e.g. numerators and denominators).

The majority of this is made possible by extensive use of the existing data structures and functions found in the brilliant python library [scikit-allel](https://github.com/cggh/scikit-allel). 
