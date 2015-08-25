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









/// parameter space is as followed (used for main optimization routine):
// long int num_snp = 0;
// int num_cellenv = 400;
// long int num_gene = 0;
// int num_etissue = 0;


// // genotype relevant:
// array<vector<string>, 22> snp_name_list;
// array<vector<long>, 22> snp_pos_list;
// // TODO:
// //array<vector<float>, 22> snp_prior_list;


// // expression relevant:
// // what we need:
// // 1. list of eQTL tissues, hashing all samples with their rpkm value;
// // 2. hashed all eQTL samples, for convenience of reading relevant rpkm data from the course file
// // 3. array of all genes (assuming all genes in the source file are those to be used)
// unordered_map<string, unordered_map<string, vector<float>>> eQTL_tissue_rep;  // hashing all eTissues to their actual rep, in which all sample from that tissue is hashed to their rpkm array
// unordered_map<string, string> eQTL_samples;  // hashing all eQTL samples to their tissues
// vector<string> gene_list;  // all genes from the source file

// // information table:
// unordered_map<string, gene_pos> gene_tss;  // TSS for all genes (including those pruned genes)
// unordered_map<string, int> gene_xymt_rep;  // map all the X, Y, MT genes



// // parameter containers:
// vector<vector<float *> *> para_cis_gene;
// vector<float *> para_snp_cellenv;
// vector<vector<float *> *> para_cellenv_gene;

// // information table:
// unordered_map<string, tuple_long> gene_cis_index;  // mapping the gene to cis snp indices (start position and end position in the snp vector)























	puts("============== leaving the optimization routine...");

}