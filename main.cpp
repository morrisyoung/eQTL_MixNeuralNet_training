/*
the outline of the entire program:

what we should have at hand by now:
1. genotype: grouped by chromosomes, and split for different individuals; read in on-demand fashion, or directly from memory;
2. expression: small enough to be fit in memory; should be loaded into memory immediately after initializing the program
3. batch variables


1. pipeline for processing the genotype data and the expression date (querying and iterating, in a mini-batch manner);
2. after getting the data, do the stochastic gradient descent algorithm;
3. pay attention to the data structure used for storing all the parameters (coefficients);
4. find a way to terminate the optimization process;
*/


// the information page of the data and the project is here:
//	https://github.com/morrisyoung/eQTL_script
// the project is here:
//	https://github.com/morrisyoung/eQTL_cplusplus


// something TODO in the main frame, e.g., some re-usable modules:
// 1. file operation module;
// 2. hashtable-in judgement module


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
#include "parameter_test.h"
#include "batch.h"



using namespace std;



// global variables definition and initialization
//===========================================================
long int num_snp = 0;		// TBD
int num_cellenv = 400;		// Specified
long int num_gene = 0;		// TBD
int num_etissue = 0;		// TBD
int num_batch = 0;			// TBD
int num_batch_hidden = 100;	// Specified
int num_individual = 0;		// TBD


//// genotype relevant:
array<vector<string>, 22> snp_name_list;
array<vector<long>, 22> snp_pos_list;
unordered_map<string, vector<vector<float>>> snp_dosage_rep;


//// expression relevant:
// what we need:
// 1. list of eQTL tissues, hashing all samples with their rpkm value;
// 2. hashed all eQTL samples, for convenience of reading relevant rpkm data from the course file
// 3. array of all genes (assuming all genes in the source file are those to be used)
unordered_map<string, unordered_map<string, vector<float>>> eQTL_tissue_rep;	// hashing all eTissues to their actual rep, in which all sample from that tissue is hashed to their rpkm array
unordered_map<string, string> eQTL_samples;										// hashing all eQTL samples to their tissues
vector<string> gene_list;														// all genes from the source file
unordered_map<string, int> gene_index_map;										// re-map those genes into their order (reversed hashing of above)
vector<string> etissue_list;													// eTissues in order
unordered_map<string, int> etissue_index_map;									// re-map those etissues into their order (reversed hashing of above)
unordered_map<string, vector<string>> esample_tissue_rep;						// esample lists of all etissues

// information table:
unordered_map<string, gene_pos> gene_tss;										// TSS for all genes (including those pruned genes)
unordered_map<string, int> gene_xymt_rep;										// map all the X, Y, MT genes


//// batch variables:
// batch variables are per genotype per sample (one individual may produce several RNA samples)
unordered_map<string, vector<float>> batch_individual;
unordered_map<string, vector<float>> batch_sample;


//// parameter containers:
vector<vector<float *>> para_cis_gene;
vector<float *> para_snp_cellenv;
vector<vector<float *>> para_cellenv_gene;
vector<float *> para_batch_batch_hidden;
vector<float *> para_batch_hidden_gene;

// information table:
unordered_map<string, tuple_long> gene_cis_index;  // mapping the gene to cis snp indices (start position and end position in the snp vector)



//// system control
// multi-threading mark
int MULTI_THREAD = 1;

//===========================================================





int main()
{
	cout << "[now enter the program]" << endl;



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
	puts("[xxx] preparing the snp info (index --> snp name and chromosome positions)...");
	num_snp = snp_info_read();  // snp_name_list; snp_pos_list
	cout << "there are " << num_snp << " snps totally." << endl;

	///* temporarily (as there are no enough space locally in VM)
	// load the genotype for all individuals on all chromosomes
	puts("[xxx] loading all dosage data for these snps for all individuals.");
	dosage_load();  // unordered_map<string, vector<vector<float>>> snp_dosage_rep;
	cout << "there are " << num_individual << " individuals." << endl;
	//*/





	//===================================== prepare the expression matrix =======================================
	puts("[xxx] loading the gene rpkm matrix...");
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
	puts("[xxx] loading the tss for genes...");
	gene_tss_load();  // gene_tss
	puts("[xxx] loading the X, Y, MT gene list...");
	gene_xymt_load();  // gene_xymt_rep





	//========================================== prepare the batch ==============================================
	// we know the num_batch and num_batch_hidden only after we read the batch source file
	puts("[xxx] batch variable values (for all individuals and samples) loading...");
	batch_load();  // load the batch variables
	cout << "there are totally " << num_batch << " batch variables (individuals and samples)..." << endl;






	// the following three must be in order
	//===================================== gene cis index data preparation ======================================
	puts("[xxx] gene meta data (cis- index) preparation...");
	gene_cis_index_init();  // gene_cis_index
	//==================================== initialize all other parameters =======================================
	puts("[xxx] parameter space initialization...");
	para_init();  // para_snp_cellenv; para_cellenv_gene; para_cis_gene
	//=========================================== set the beta prior =============================================
	puts("[xxx] beta prior values (from GTEx) loading...");
	beta_prior_fill();  // must happen after the above procedure







	//======================================= main optimization routine ==========================================
	optimize();






	//======================================= testing with the test subset =======================================
	puts("[xxx] testing the learned models...");
	para_test();  // currently this is empty




	//================================= save the parameters and release memory ===================================
	puts("[xxx] saving the models...");
	para_save();  // para_cis_gene; para_snp_cellenv; para_cellenv_gene
	cout << "Optimization done! Please find the results in 'result' folder." << endl;
	puts("[xxx] releasing the parameter space...");
	para_release();





	//============== timing ends ================
    gettimeofday(&time_end, NULL);
    diff = (double)(time_end.tv_sec-time_start.tv_sec) + (double)(time_end.tv_usec-time_start.tv_usec)/1000000;
    printf("Time used totally is %f seconds.\n", diff);

	cout << "[now leave the program]\n";
	return 0;
}

