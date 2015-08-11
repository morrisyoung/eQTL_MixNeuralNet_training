/*
the outline of the entire program:

what we should have at hand by now:
1. genotype: grouped by chromosomes, and split for different individuals; read on-demand;
2. expression: small enough to be fit in memory; should be loaded into memory immediately after initializing the program



1. pipeline for processing the genotype data and the expression date (querying and iterating, in a mini-batch manner);
2. after getting the data, do the stochastic gradient descent algorithm;
3. pay attention to the data structure used for storing all the parameters (coefficients);
4. find a way to terminate the optimization process;

*/


#include <iostream>
#include <sys/types.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include "genotype.h"


using namespace std;


int main()
{
	cout << "This is the entrance of the program...\n";


	int chr = 1;
	char individual[20] = "GTEX-TKQ1";
	dosage_load(chr, individual);
	// rpkm_load();








	//optimization();

	cout << "Optimization done! Please find the results in 'result' folder.\n";

	return 0;
}