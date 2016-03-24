/*
the testing routine:

1. read the saved parameters from the source files;
2. read the genotype data (and the batch variables) of testing dataset, to perform the prediction;
3. save the expected expression array for all the testing samples (from different tissue types) -- etissue list, and seperate files for each tissue;
4. [outside this program] calculate and plot the Pearson correlation plot for the testing dataset, for each tissue type;


several other issues:
1. for the testing, it's better if we can jump the process of saving all the learned parameters into file, as the numerical precision might be compromised;
2. for this testing, we know the sample batch values; while in practice, we don't actually know these sample batch variables (as we don't even have these samples);
	this to some extend limits the practical gene expression prediction of this model
3. ...

*/




/*
Mar.22:
I'm re-formating the files for testing set
what we need:
1. "./data_real/list_samples_test.txt"
where do they come from:
1. "./phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples_test"
what to be changed:
1. none





(Mar.23, 2016)
after checking the training set, I guess I still have the following to do:
//// [important] what to do to change the dataset from real dataset to simulated dataset?
// 1. change the file header (source data), from (char filename_data_source[] = "../data_real/";) to (char filename_data_source[] = "../data_simu/";)
// 2. [no need] change the file header (initial parameters), from (char file_para_init[] = "../result_init/";) to (char file_para_init[] = "../result_init_simu/";)
// 3. change (int num_cellenv = 400) and (int num_batch_hidden = 100), global variables, to (int num_cellenv = 400) and (int num_batch_hidden = 50)
// 4. change NUM_CHR (Macro), the number of chromosomes
// 5. change the indicator (of whether this is real data or not) from (int indicator_real = 1) to (int indicator_real = 0)
// 6. ...







TODO:
1. Matrix data structure, might not need to change;
2. output the testing errors/likelihood in the testing set (optional)





*/








// standard libraries:
#include <iostream>
#include <sys/types.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <unordered_map>
#include <string.h>
#include <string>
#include <array>
#include <forward_list>
#include <utility>
#include <vector>
#include <sys/time.h>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
// sub-routines:
#include "global.h"
#include "genotype.h"
#include "expression.h"
#include "batch.h"
#include "basic.h"
#include "test_para_read.h"
#include "test_predict.h"
#include "test_save.h"
#include "test_main.h"



using namespace std;



// global variables definition and initialization
//===========================================================
long int num_snp = 0;		// TBD
int num_cellenv = 400;		// Specified
long int num_gene = 0;		// TBD
int num_etissue = 0;		// TBD
int num_batch = 0;			// TBD
int num_batch_hidden = 50;	// Specified
int num_individual = 0;		// TBD


//// genotype relevant:
array<vector<string>, NUM_CHR> snp_name_list;
array<vector<long>, NUM_CHR> snp_pos_list;
unordered_map<string, vector<vector<float>>> snp_dosage_rep;


//// expression relevant: (I will redundently load also the training set)
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


// (Mar.22, 2016)
// I will replicate the "unordered_map<string, unordered_map<string, vector<float>>> eQTL_tissue_rep" and "unordered_map<string, string> eQTL_samples" and "unordered_map<string, vector<string>> esample_tissue_rep;" for the testing dataset:
unordered_map<string, unordered_map<string, vector<float>>> eQTL_tissue_rep_test;	// hashing all eTissues to their actual rep, in which all sample from that tissue is hashed to their rpkm array
unordered_map<string, string> eQTL_samples_test;	// hashing all eQTL samples to their tissues
unordered_map<string, vector<string>> esample_tissue_rep_test;	// esample lists of all etissues
// for prediction:
unordered_map<string, unordered_map<string, vector<float>>> eQTL_tissue_rep_predict;	// the predicted version of above one



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
//int MULTI_THREAD = 1;
// We won't use multi-threading for this testing program


//// file name space
// (note: if we standadize the source data format and name, we only need the upper folder name)
// data source:
//char filename_data_source[] = "../data_real/";
char filename_data_source[] = "../data_simu/";

// parameter source:
// we will always load parameters from the directory "../result/"
char file_para_init[] = "../result/";




//// indicator of whether this is working on real dataset
// this will be used to dicide parameters in some basic functions, like "sample ID to individual ID" function
int indicator_real = 0;
//===========================================================





int main()
{
	cout << "[now enter the testing program]" << endl;



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
	// loading the training dataset
	puts("[xxx] loading the gene rpkm matrix (training)...");
	char filename1[100];
	filename1[0] = '\0';
	strcat(filename1, filename_data_source);
	strcat(filename1, "list_samples_train.txt");
	char filename2[100];
	filename2[0] = '\0';
	strcat(filename2, filename_data_source);
	strcat(filename2, "expression.txt");

	num_gene = gene_train_load(filename1, filename2);  // eQTL_samples; gene_list; eQTL_tissue_rep
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


	// loading the testing dataset
	puts("[xxx] loading the gene rpkm matrix (testing)...");
	filename1[0] = '\0';
	strcat(filename1, filename_data_source);
	strcat(filename1, "list_samples_test.txt");
	filename2[0] = '\0';
	strcat(filename2, filename_data_source);
	strcat(filename2, "expression.txt");

	gene_test_load(filename1, filename2);  // eQTL_samples; gene_list; eQTL_tissue_rep


	// loading others
	puts("[xxx] loading the tss for genes...");
	gene_tss_load();  // gene_tss
	puts("[xxx] loading the X, Y, MT gene list...");
	gene_xymt_load();  // gene_xymt_rep



	// refine the following, as the reference are all built on the training dataset
	//vector<string> etissue_list;
	//unordered_map<string, int> etissue_index_map;
	char filename[100] = "../result/etissue_list.txt";
	FILE * file_in = fopen(filename, "r");
	if(file_in == NULL)
	{
		fputs("File error\n", stderr); exit (1);
	}
	int input_length = 100;
	char input[input_length];
	int count = 0;
	while(fgets(input, input_length, file_in) != NULL)
	{
		trim(input);

		const char * sep = "\t";
		char * p;
		p = strtok(input, sep);
		string etissue = p;
		etissue_list[count] = etissue;
		etissue_index_map[etissue] = count;
		count++;
	}
	fclose (file_in);


	// need to initialize the following (with "eQTL_tissue_rep_predict"):
	//unordered_map<string, unordered_map<string, vector<float>>> eQTL_tissue_rep_predict;
	for(int i=0; i<etissue_list.size(); i++)
	{
		string etissue = etissue_list[i];
		unordered_map<string, vector<float>> map;
		eQTL_tissue_rep_predict.emplace(etissue, map);
		for(auto it = eQTL_tissue_rep_test[etissue].begin(); it != eQTL_tissue_rep_test[etissue].end(); ++it)
		{
			string esample = it->first;
			vector<float> vec;
			eQTL_tissue_rep_predict[etissue].emplace(esample, vec);
			for(int j=0; j<num_gene; j++)
			{
				eQTL_tissue_rep_predict[etissue][esample].push_back(0);
			}
		}
	}







	//========================================== prepare the batch ==============================================
	// we know the num_batch and num_batch_hidden only after we read the batch source file
	puts("[xxx] batch variable values (for all individuals and samples) loading...");
	batch_load();  // load the batch variables
	cout << "there are totally " << num_batch << " batch variables (individuals and samples)..." << endl;






	// the following three must be in order
	//===================================== gene cis index data preparation ======================================
	puts("[xxx] gene meta data (cis- index) preparation...");
	gene_cis_index_init();  // gene_cis_index
	//==================================== initialize all parameters from learned results =======================================
	puts("[xxx] parameter space initialization and loading...");
	para_init();		// loading the learned parameters









	//======================================= main testing routine ==========================================
	predict();  // save unordered_map<string, unordered_map<string, vector<float>>> eQTL_tissue_rep_predict









	//================================= save the parameters and release memory ===================================
	puts("[xxx] saving the predicted expression arrays (samples) for all tissues...");
	predict_save();
	cout << "Optimization done! Please find the results in 'result_predict' folder." << endl;
	puts("[xxx] releasing the parameter space...");
	para_release();





	//============== timing ends ================
	gettimeofday(&time_end, NULL);
	diff = (double)(time_end.tv_sec-time_start.tv_sec) + (double)(time_end.tv_usec-time_start.tv_usec)/1000000;
	printf("Time used totally is %f seconds.\n", diff);

	cout << "[now leave the testing program]\n";
	return 0;
}


