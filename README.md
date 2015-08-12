Implementation (C++) for gene expression with regulation.

There are several components in this project:

1. main.cpp: This is the entrance of the whole project. It will process some command line options, initialize some data and call other routines (mainly the routine for optimization).

2. genotype.cpp: This will read the genotype dosage data from prepared repository. It will also process the pruned SNPs according to the information we have at hand (dosage data of pruned SNPs, associated un-pruned SNPs/representative SNPs and the association signal R^2, and for the regularization we have the prior knowledge score).

3. expression.cpp: This will read and prepare the expression data.

4. optimization.cpp: This is the main routine of the optimization program. It will call other routine to get genotype data and work on expression data, for stochastic gradient descent. It may call other subroutines for detailed stochastic gradient descent algorithm.

5. basic.cpp: Some basic functions that may be used by the whole program.

6. xxx