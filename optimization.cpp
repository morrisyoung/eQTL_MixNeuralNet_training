// the main optimization routine

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unordered_map>
#include <string>
#include <vector>
#include "basic.h"
#include <forward_list>
#include <utility>
#include "genotype.h"
#include "optimization.h"
#include "global.h"
#include "main.h"



using namespace std;



void optimize()
{
	puts("============== entering the optimization routine...");


	//==================== get the actual genotype from (currently) files: for one individual on all chromosomes =====================
	puts("get the dosage data for one individual (test)...");
	// test sample (individual):
	string individual = "GTEX-TKQ1";  // for testing

	array<vector<float>, 22> snp_dosage_list;

	int i;
	for(i=0; i<22; i++)
	{
		int chr = i+1;
		vector<float> vec;
		snp_dosage_list[i] = vec;
		snp_dosage_load(&snp_dosage_list[i], chr, individual);
	}
	puts("snp dosage loading (test) done!");
	//===== simple test =====
	//cout << snp_dosage_list[2].size() << endl;
	//cout << snp_dosage_list[2][7] << endl;
	//===============================================================
	//============================================================================================================










	// 3 major parts for this routine:

	// 1. build the initial parameter space to be learned (in a trivial or non-trivial way)

	// 2. write the procedure to update these parameters in this program

	// 3. stop in some ways, and save the learned parameters to disk (interface); use another routine to check the predictive precision





	// part#1:
	// cis- parameter list for each gene (across different tissues); all SNP to cell_env parameter table; cell env to all 20k genes (across different tissues)























	puts("============== leaving the optimization routine...");

}