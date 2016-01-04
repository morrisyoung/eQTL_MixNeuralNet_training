/*
the outline of the entire program:

what we should have at hand by now (to train the model):
1. genotype: grouped by chromosomes, and split for different individuals; read in on-demand fashion, or directly from memory;
2. expression: small enough to be fit in memory; should be loaded into memory immediately after initializing the program
3. batch variables

some notes:
1. pipeline for processing the genotype data and the expression data (querying and iterating, in a mini-batch manner);
2. after getting the data, do the stochastic gradient descent algorithm;
3. pay attention to the data structure used for storing all the parameters (coefficients);
4. find a way to terminate the optimization process;
*/


//// the information page of the data and the project is here:
//	https://github.com/morrisyoung/eQTL_script
// the project is here:
//	https://github.com/morrisyoung/eQTL_cplusplus


//// notes (Dec.30, 2015):
// 1. should try to make the training program as flexible as possible (operating on objects, other than specific data structure);
// 2. should add the single chromosome mode, as in simulation, we most probably simulate only one chromosome;
// 3. ...


//// source data re-organizing:
// NOTE: we will re-organize all the data source files, and make them standard (for convenience of working on both real data and simulated data)
// currently we need the following source files (data), not including the initialization of the parameter space:
/*
1. "./data_real/genotype/chrXXX/SNP_dosage_IndividualID.txt"		// the genotype of this individual
2. "./data_real/genotype/chrXXX/SNP_info.txt"						// SNP names and positions
3. "./data_real/list_individuals.txt"								// the list of individual IDs
4. "./data_real/list_samples_train.txt"								// the list of samples (from different tissues) used in training
5. "./data_real/expression.txt"										// the expression matrix
6. "./data_real/gene_tss.txt"										// the chr# and TSS of all genes
7. "./data_real/gene_xymt.txt"										// the genes to be excluded in the current model (X, Y, MT genes)
8. "./data_real/batch_individuals.txt"								// the individual batch variables
9. "./data_real/batch_samples.txt"									// the tissue sample batch variables (individuals x tissues)
*/
// where do they come from (the last processed data source files):
/*
1. "./genotype_185_dosage_matrix_qc/chrXXX/SNP_dosage_IndividualID.txt"
2. "./genotype_185_dosage_matrix_qc/chrXXX/SNP_info.txt"
3. "./list_individual.txt"
4. "./phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples_train"
5. "./GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_normalized"
6. "./gencode.v18.genes.patched_contigs.gtf_gene_tss"
7. "./gencode.v18.genes.patched_contigs.gtf_gene_xymt"
8. "./batch_var_individual.txt"
9. "./batch_var_sample.txt"
*/
// what's changed in this transmission:
/*
1. none
2. none
3. none
4. none
5. removed the first two useless lines
6. none
7. none
8. none
9. none
*/
// what needs to be further done in the simulation code:
/*
1. genotype: organize as the real data, even if there is only one chromosome (this is the case in simulated data)
2. as above
3. list_individuals.txt: ...
4. list_samples_train.txt: ...
5. expression.txt: re-organize the simulated data as the real data (a full matrix)
6. gene_tss.txt: ...
7. gene_xymt.txt: maybe make this empty?
8. batch_individuals.txt: will keep the first line -- variable name
9. batch_samples.txt: will keep the first line -- variable name
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
#include "batch.h"
#include "io_file.h"			// some Python-like file IO operations
#include "op_line.h"			// some Python-like line(string) operations




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


//// parameter containers: (we need to initialize all of them in an appropriate way)
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



//// file name space
// (note: if we standadize the source data format and name, we only need the upper folder name)
char filename_data_source[] = "../data_real/";
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
	char filename1[100];
	filename1[0] = '\0';
	strcat(filename1, filename_data_source);
	strcat(filename1, "list_samples_train.txt");
	char filename2[100];
	filename2[0] = '\0';
	strcat(filename2, filename_data_source);
	strcat(filename2, "expression.txt");

	num_gene = gene_rpkm_load(filename1, filename2);  // eQTL_samples; gene_list; eQTL_tissue_rep
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
	//
	// if we initialize all the parameters with pre-prepared data, we don't need the results from GTEx project
	//
	//puts("[xxx] beta prior values (from GTEx) loading...");
	//beta_prior_fill();  // must happen after the above procedure
	//








	//======================================= main optimization routine ==========================================
	optimize();







	//================================= save the parameters and release memory ===================================
	puts("[xxx] saving the models...");
	para_save();  // para_cis_gene; para_snp_cellenv; para_cellenv_gene; para_batch_batch_hidden; para_batch_hidden_gene
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

