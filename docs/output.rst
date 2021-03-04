************
Understanding pixy output
************

Output file contents
================

pixy outputs a slightly different file type for each summary statistic it calculates. The contents of the columns of these output files are detailed below.

Within population nucleotide diversity (pi)
-----------

``pop`` - The ID of the population from the population file

``chromosome`` - The chromosome/contig 

``window_pos_1`` - The first position of the genomic window

``window_pos_2`` - The last position of the genomic window

``avg_pi`` - Average per site nucleotide diversity for the window

``no_sites`` - The total number of sites in the window that have at least one valid genotype. This statistic is included for the user, and not directly used in any calculations.

``count_diffs`` - The raw number of pairwise differences between all genotypes in the window. This is the numerator of avg_pi.

``count_comparisons`` - The raw number of non-missing pairwise comparisons between all genotypes in the window (i.e. cases where two genotypes were compared and both were valid). This is the denominator of avg_pi.

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

``count_diffs`` - The raw number of pairwise, cross-population differences between all genotypes. This is the numerator of avg_dxy.

``count_comparisons`` - The raw number of non-missing pairwise cross-population comparisons between all genotypes in the window (i.e. cases where two genotypes were compared and both were valid). This is the denominator of avg_dxy.

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

.. code:: r

    # Example R Script for simple output plots 
    # Here, we use pi and dxy output files directly from pixy.

    library(ggplot2)

    # Provide path to input. Can be pi or Dxy. 
    # NOTE: this is the only line you should have to edit to run this code:
    inp<-read.table("pixy_dxy.txt",sep="\t",header=T)

    # Find the chromosome names and order them: first numerical order, then any non-numerical chromosomes
    #   e.g., chr1, chr2, chr22, chrX
    chroms <- unique(inp$chromosome)
    chrOrder <- sort(chroms)
    inp$chrOrder <- factor(inp$chromosome,levels=chrOrder)

    # Plot pi for each population found in the input file
    # Saves a copy of each plot in the working directory
    if("avg_pi" %in% colnames(inp)){
        pops <- unique(inp$pop)
        for (p in pops){
            thisPop <- subset(inp, pop == p)
            # Plot stats along all chromosomes:
            popPlot <- ggplot(thisPop, aes(window_pos_1, avg_pi, color=chrOrder)) +
                geom_point()+
                facet_grid(. ~ chrOrder)+
                labs(title=paste("Pi for population", p))+
                labs(x="Position of window start", y="Pi")+
                scale_color_manual(values=rep(c("black","gray"),ceiling((length(chrOrder)/2))))+
                theme_classic()+
                theme(legend.position = "none")
            ggsave(paste("piplot_", p,".png", sep=""), plot = popPlot, device = "png", dpi = 300)
            }
    } else {
        print("Pi not found in this file")
    }

    # Plot Dxy for each combination of populations found in the input file
    # Saves a copy of each plot in the working directory
    if("avg_dxy" %in% colnames(inp)){
        # Get each unique combination of populations
        pops <- unique(inp[c("pop1", "pop2")])
        for (p in 1:nrow(pops)){
            combo <- pops[p,]
            thisPop <- subset(inp, pop1 == combo$pop1[[1]] & pop2 == combo$pop2[[1]])
            # Plot stats along all chromosomes:
            popPlot <- ggplot(thisPop, aes(window_pos_1, avg_dxy, color=chrOrder)) + 
                geom_point()+
                facet_grid(. ~ chrOrder)+
                labs(title=paste("Dxy for", combo$pop1[[1]], "&", combo$pop2[[1]]))+
                labs(x="Position of window start", y="Dxy")+
                theme(legend.position = "none")+
               scale_color_manual(values=rep(c("black","gray"),ceiling((length(chrOrder)/2))))+
               theme_classic()+
               theme(legend.position = "none") 
            ggsave(paste("dxyplot_", combo$pop1[[1]], "_", combo$pop2[[1]],".png", sep=""), plot = popPlot, device = "png", dpi = 300)
        }
    } else {
        print("Dxy not found in this file")
    }


Post-hoc aggregating
------------------------

Note that if the user wishes to combine information across windows (e.g. by averaging) after the fact, they should sum the raw counts and recompute the differences/comparisons ratios, and not take an average of the summary statistics themselves. 

For example, to get average pi or dxy for two windows, the correct forumla is: 

.. parsed-literal::

    (window 1 count_diffs + window 2 count_diffs) / (window 1 comparisons + window 2 comparisons)



 
