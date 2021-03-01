************
Companion files
************

Populations file
--------------------

*Description*

The populations file is used to assign individuals to populations. For pairwise statistics (dxy and fst), all pairwise comparisons between populations are generated. It is specified at the command line with the ``--populations`` argument.

*Format*

A headerless, tab separated file with the first column containing the sample name (as represented in the VCF), and the second column containing an arbitary population ID.

*Example*

.. parsed-literal::
    ERS223827	BFS
    ERS223759	BFS
    ERS223750	BFS
    ERS223967	AFS
    ERS223970	AFS
    ERS223924	AFS
    ERS224300	AFS
    ERS224168	KES
    ERS224314	KES

BED file
--------------------

*Description*

The BED file is used to manually defined genome regions over which to calculate summary statistics. It is specified at the command line with the ``--bed_file``.

*Format*

A headerless, tab separated file with the first column containing a chromosome ID, the second column containing the first position of a window, and the third column containing the last potion of the window. Each row = one window.

*Example*


.. parsed-literal::
    X	1	10000
    X	10001	20000
    X	20001	30000
    X	30001	40000
    X	40001	50000
    X	50001	60000
    1	1	10000
    1	10001	20000
    1	20001	30000
    1	30001	40000
    1	40001	50000
    1	50001	60000

Sites file
--------------------

*Description*

The sites is used to restrict calculations to specific sites (e.g. 4-fold degenerate sites, introns, etc.). It is specified at the command line with the ``--sites_file``.

*Format*

A headerless, tab separated file with the first column containing a chromosome ID (CHROM), the second column containing the position (POS) of a single site. Each row = one site.

*Example*

.. parsed-literal::
    X	1
    X	2
    X	3
    X	4
    X	5
    X	6
    X	7
    X	8
    X	9
    X	10
    X	11
    X	12
    X	13
    X	14
    X	15
    1	1
    1	2
    1	3
    1	4
    1	5
    1	6
    1	7
    1	8
    1	9
    1	10
    1	11
    1	12
    1	13
    1	14
    1	15