************
Understanding pixy output
************

Output file contents
================

pixy outputs a slightly different file type for each summary statistic it calculates. The contents of the columns of this output file are detailed below.

Within population nucleotide diversity (pi)
-----------

``pop`` - The ID of the population from the population file

``chromosome`` - The chromosome/contig 

``window_pos_1`` - The first position of the genomic window

``window_pos_2`` - The last position of the genomic window

``avg_pi`` - Average per site nucleotide diversity for the window

``no_sites`` - The total number of sites in the window that have at least one valid genotype. This statistic is included for the user, and not directly used in any calculations.

``count_diffs`` - The raw number of pairwise differences between all genotypes in the window.

``count_comparisons`` - The raw number of non-missing pairwise comparisons between all genotypes in the window (i.e. cases where two genotypes were compared and both were valid). 

``count_missing`` - The raw number of missing pairwise comparisons between all genotypes in the window (i.e. cases where two genotypes were compared and at least one was missing). 

Between population nucleotide divergence (dxy)
-----------

``pop1`` - The ID of the first population from the population file.

``pop2`` - The ID of the second population from the population file.

``chromosome`` - The chromosome/contig.

``window_pos_1`` - The first position of the genomic window.

``window_pos_2`` - The last position of the genomic window.

``avg_dxy`` - Average per site nucleotide divergence for the window.

``no_sites`` - The total number of sites in the window that have at least one valid genotype in both populations. This statistic is included for the user, and not directly used in any calculations.

``count_diffs`` - The raw number of pairwise, cross-population differences between all genotypes.

``count_comparisons`` - The raw number of non-missing pairwise cross-population comparisons between all genotypes in the window (i.e. cases where two genotypes were compared and both were valid) .

``count_missing`` - The raw number of missing pairwise cross-population comparisons between all genotypes in the window (i.e. cases where two genotypes were compared and at least one was missing). This statistic is included for the user, and not directly used in any calculations. 

Weir and Cockerham's FST (fst)
-----------

``pop1`` - The ID of the first population from the population file.

``pop2`` - The ID of the second population from the population file.

``chromosome`` - The chromosome/contig.

``window_pos_1`` - The first position of the genomic window.

``window_pos_2`` - The last position of the genomic window.

``avg_wc_fst`` - Average Weir and Cockerham's FST for the window (per SNP, not per site).

``no_snps`` - Total number of variable sites (SNPs) in the window.

Working with pixy output data
================

Plotting results
------------------------

Post-hoc aggregating
------------------------

Note that if the user wishes to combine information across windows (e.g. by averaging) after the fact, they should use the raw counts, and not the average of the summary statistics themselves. 

For example, to get average pi or dxy for two windows, the correct forumla is: 

.. parsed-literal::

    (window 1 count_diffs + window 2 count_diffs) / (window 1 comparisons + window 2 comparisons)



 
