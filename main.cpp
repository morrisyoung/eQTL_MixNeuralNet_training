/*
the outline of the entire program:

what we should have at hand by now:
1. genotype: grouped by chromosomes, and split for different individuals; read in on-demand fashion;
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
#include <unordered_map>
#include <string.h>
#include <string>
#include <array>
#include <forward_list>
#include <utility>
#include "basic.h"
#include "expression.h"
#include "optimization.h"
#include "global.h"
#include "parameter_init.h"
#include "main.h"
#include <vector>
#include "parameter_save.h"
#include <sys/time.h>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */



using namespace std;



// global variables definition and initialization
//===========================================================
long int num_snp = 0;
int num_cellenv = 400;
long int num_gene = 0;
int num_etissue = 0;
// TODO: for the following two, we should later on prepared read data from file and automatically get their values
int num_batch = 100;
int num_batch_hidden = 100;
int num_individual = 0;


// genotype relevant:
array<vector<string>, 22> snp_name_list;
array<vector<long>, 22> snp_pos_list;
unordered_map<string, vector<vector<float>>> snp_dosage_rep;


// expression relevant:
// what we need:
// 1. list of eQTL tissues, hashing all samples with their rpkm value;
// 2. hashed all eQTL samples, for convenience of reading relevant rpkm data from the course file
// 3. array of all genes (assuming all genes in the source file are those to be used)
unordered_map<string, unordered_map<string, vector<float>>> eQTL_tissue_rep;  // hashing all eTissues to their actual rep, in which all sample from that tissue is hashed to their rpkm array
unordered_map<string, string> eQTL_samples;  // hashing all eQTL samples to their tissues
vector<string> gene_list;  // all genes from the source file
vector<string> etissue_list;  // eTissues in order
unordered_map<string, int> etissue_index_map;  // re-map those etissues into their order (reversed hashing above)
unordered_map<string, vector<string>> esample_tissue_rep;  // esample lists of all etissues

// information table:
unordered_map<string, gene_pos> gene_tss;  // TSS for all genes (including those pruned genes)
unordered_map<string, int> gene_xymt_rep;  // map all the X, Y, MT genes


// parameter containers:
vector<vector<float *>> para_cis_gene;
vector<float *> para_snp_cellenv;
vector<vector<float *>> para_cellenv_gene;
vector<float *> para_batch_batch_hidden;
vector<float *> para_batch_hidden_gene;

// information table:
unordered_map<string, tuple_long> gene_cis_index;  // mapping the gene to cis snp indices (start position and end position in the snp vector)

//===========================================================





int main()
{
	cout << "[now enter the program]\n";

	//============== timing starts ================
    struct timeval time_start;
    struct timeval time_end;
    double diff;
    gettimeofday(&time_start, NULL);


	// maybe here accept some command lines
	//
	//
	//
	//
	//
	//
	//
	//
	//


	//======================================= prepare the snp information ========================================
	puts("preparing the snp info (index --> snp name and chromosome positions)...");
	num_snp = snp_info_read();  // snp_name_list; snp_pos_list
	cout << "there are " << num_snp << " snps totally." << endl;

	/* temporarily (as there are no enough space locally in VM)
	// load the genotype for all individuals on all chromosomes
	puts("loading all dosage data for these snps for all individuals.");
	dosage_load();  // unordered_map<string, vector<vector<float>>> snp_dosage_rep;
	cout << "there are " << num_individual << " individuals." << endl;
	*/
	//============================================================================================================


	//===================================== prepare the expression matrix =======================================
	num_gene = gene_rpkm_load();  // eQTL_samples; gene_list; eQTL_tissue_rep
	num_etissue = eQTL_tissue_rep.size();
	cout << "there are " << num_gene << " genes totally." << endl;
	cout << "there are totally " << eQTL_samples.size() << " training samples from different eQTL tissues." << endl;
	cout << "there are " << num_etissue << " eTissues in the current framework." << endl;
	puts("number of training samples in each eTissue are as followed:");
	for(auto it=eQTL_tissue_rep.begin(); it != eQTL_tissue_rep.end(); ++it)
	{
		string etissue = it->first;
		cout << etissue << ":" << (it->second).size() << endl;
	}
	gene_tss_load();  // gene_tss
	gene_xymt_load();  // gene_xymt_rep
	//============================================================================================================




	//===================================== gene cis index data preparation ======================================
	puts("gene meta data (cis- index) preparation...");
	gene_cis_index_init();  // gene_cis_index
	//============================================================================================================
	//===================================== initialize all the parameters ========================================
	puts("parameter space initialization...");
	para_init();
	//============================================================================================================






	//======================================= main optimization routine ==========================================
	//optimize();
	//============================================================================================================


	cout << "test..." << endl;



	//============== save parameters and release memory ================
	para_save();  // para_cis_gene; para_snp_cellenv; para_cellenv_gene
	cout << "Optimization done! Please find the results in 'result' folder.\n";
	para_release();

	//============== timing ends ================
    gettimeofday(&time_end, NULL);
    diff = (double)(time_end.tv_sec-time_start.tv_sec) + (double)(time_end.tv_usec-time_start.tv_usec)/1000000;
    printf("Time used totally is %f seconds.\n", diff);

	cout << "[now leave the program]\n";
	return 0;
}

