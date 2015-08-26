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



//====================================== local global variables ========================================
// these variables are specially designed for this routine -- optimization
// need to initialize some local containers:
array<vector<float>, 22> snp_dosage_list;
vector<float> gene_rpkm_exp;  // with length "num_gene"
vector<float> cellenv_hidden_var;  // with length "num_cellenv"

// parameter derivative containers:
vector<vector<float *>> para_dev_cis_gene;
vector<float *> para_dev_snp_cellenv;
vector<vector<float *>> para_dev_cellenv_gene;


// learning control parameters:
int iter_learn_out = 5;  // iteration across all tissues
int iter_learn_in = 20;  // iteration across all samples from one tissue
int batch_size = 15;
int rate_learner = 1;  // the learning rate
// TODO: regularization coefficients (maybe most of them are pretty empirical)
//
//
//

//======================================================================================================



void opt_para_init()
{
	puts("opt_para_init..");

	//=============== snp_dosage_list ===============
	for(int i=0; i<22; i++)
	{
		for(int j=0; j<snp_name_list[i].size(); j++)
		{
			snp_dosage_list[i].push_back(0);
		}
	}

	//=============== gene_rpkm_exp ===============
	for(int i=0; i<num_gene; i++)
	{
		gene_rpkm_exp.push_back(0);
	}

	//=============== cellenv_hidden_var ===============
	for(int i=0; i<num_cellenv; i++)
	{
		cellenv_hidden_var.push_back(0);
	}

	//=============== para_dev_snp_cellenv ===============
	for(int i=0; i<num_cellenv; i++)
	{
		float * p = (float *)malloc( sizeof(float) * num_snp );
		para_dev_snp_cellenv.push_back(p);
	}

	//=============== para_dev_cellenv_gene ===============
	for(int j=0; j<num_etissue; j++)
	{
		vector<float *> vec;
		para_dev_cellenv_gene.push_back(vec);
		for(int i=0; i<num_gene; i++)
		{
			float * p = (float *)malloc( sizeof(float) * num_cellenv );
			para_dev_cellenv_gene[j].push_back(p);
		}
	}

	//=============== para_dev_cis_gene ===============
	for(int j=0; j<num_etissue; j++)
	{
		vector<float *> vec;
		para_dev_cis_gene.push_back(vec);
		for(int i=0; i<num_gene; i++)
		{
			string gene = gene_list[i];
			unordered_map<string, int>::const_iterator got = gene_xymt_rep.find(gene);
			if ( got != gene_xymt_rep.end() )
			{
				float * p = NULL;  // null pointer
				para_dev_cis_gene[j].push_back(p);
				continue;
			}
			else
			{
				long first = gene_cis_index[gene].first;  // index
				long second = gene_cis_index[gene].second;  // index
				long amount = second - first + 1;
				float * p = (float *)malloc( sizeof(float) * amount );
				para_dev_cis_gene[j].push_back(p);
			}
		}
	}


}



void opt_para_release()
{
	//=============== para_dev_snp_cellenv ===============
	for(int i=0; i<num_cellenv; i++)
	{
		free(para_dev_snp_cellenv[i]);
	}

	//=============== para_dev_cellenv_gene ===============
	for(int j=0; j<num_etissue; j++)
	{
		for(int i=0; i<num_gene; i++)
		{
			free(para_dev_cellenv_gene[j][i]);
		}
	}

	//=============== para_dev_cis_gene ===============
	for(int j=0; j<num_etissue; j++)
	{
		for(int i=0; i<num_gene; i++)
		{
			free(para_dev_cis_gene[j][i]);
		}
	}

}



// transform the sample ID (like "GTEX-R55E-0826-SM-2TC5M") into individual ID (here is the first 9 digits)
string sample_to_individual(string sample)
{
	string individual;
	for(int i=0; i<9; i++)  // TODO: in next stage, the individual ID may be longer
	{
		individual.push_back(sample[i]);
	}

	return individual;
}




// forward and backward propagation for one mini-batch
void forward_backward_prop_batch(string etissue, int pos_start, int num_esample)
{

	//=================================== our parameter space ===================================
	// long int num_snp = 0;
	// int num_cellenv = 400;
	// long int num_gene = 0;
	// int num_etissue = 0;
	// array<vector<string>, 22> snp_name_list;
	// array<vector<long>, 22> snp_pos_list;
	// unordered_map<string, unordered_map<string, vector<float>>> eQTL_tissue_rep;  // hashing all eTissues to their actual rep, in which all sample from that tissue is hashed to their rpkm array
	// unordered_map<string, string> eQTL_samples;  // hashing all eQTL samples to their tissues
	// vector<string> gene_list;  // all genes from the source file
	// vector<string> etissue_list;  // eTissues in order
	// vector<vector<string>> esample_tissue_list;  // esample lists of all etissues
	// unordered_map<string, gene_pos> gene_tss;  // TSS for all genes (including those pruned genes)
	// unordered_map<string, int> gene_xymt_rep;  // map all the X, Y, MT genes
	// vector<vector<float *>> para_cis_gene;
	// vector<float *> para_snp_cellenv;
	// vector<vector<float *>> para_cellenv_gene;
	// unordered_map<string, tuple_long> gene_cis_index;  // mapping the gene to cis snp indices (start position and end position in the snp vector)
	//============================================================================================


	for(int count=0; count<batch_size; count++)
	{
		int pos = (pos_start + count) % (num_esample);
		string esample = esample_tissue_rep[etissue][pos];
		string individual = sample_to_individual(esample);

		//=================================================== init ============================================================
		// get the: 0. esample and individual; 1. genotype; 2. expression data.
		// to: 1. forward_backward propagation



		// // genotype dosage data
		// cout << "getting the dosage data for individual #" << individual << endl;
		// //string individual = "GTEX-TKQ1";  // for testing sample (individual)
		// snp_dosage_load(&snp_dosage_list, individual);  // snp dosage data for one individual across all chromosomes
		// // expression rpkm data: eQTL_tissue_rep[etissue][esample]
		// cout << "we have this amount of genes expressed in this individual:" << eQTL_tissue_rep[etissue][esample].size() << endl;




		//========================================================================
		// two step: forward propagation; backward propagation (gradient descent)
		//========================================================================
		// step#1: ... (cis-; cell env)

		// ************** part1: cis- ***************

		// for cis-: 1. if this is a XYMT gene, jump; 2. we use "gene_cis_index" to get the length of the cis- parameter array



		// ************** part2: cell env relevant parameters **************


		//========================================================================
		// step#2: ... (cis- parameters;  cell env relevant parameters)

		// ************** part1: cis- ***************
		// for gene in "gene_list":
		// if gene in gene_xymt_rep: continue;
		// else:
		// int chr = gene_tss[gene]
		// long start = gene_cis_index[gene].first;
		// long end = gene_cis_index[gene].second;

		// for index= start; index<= end; index++:
		// dosage = snp_dosage_list[chr-1][index];


		// ************** part2: cell env relevant parameters **************
		// for gene in gene_list:


	}

	// aggragate the results from all samples in this mini- batch, in some way


}




void gradient_descent()
{
	// for all parameters in our scope, we do p = p - rate_learner * dp (we have all the components in the right hand, as followed)

	// parameter containers:
	//vector<vector<float *>> para_cis_gene;
	//vector<float *> para_snp_cellenv;
	//vector<vector<float *>> para_cellenv_gene;

	// parameter derivative containers:
	//vector<vector<float *>> para_dev_cis_gene;
	//vector<float *> para_dev_snp_cellenv;
	//vector<vector<float *>> para_dev_cellenv_gene;

	//====================== para_cis_gene ==========================
	for(int i=0; i<num_etissue; i++)
	{
		string etissue = etissue_list[i];
		for(int j=0; j<num_gene; j++)
		{
			string gene = gene_list[j];
			unordered_map<string, int>::const_iterator got = gene_xymt_rep.find(gene);
			if ( got != gene_xymt_rep.end() )
			{
				continue;
			}
			else
			{
				int num = gene_cis_index[gene].second - gene_cis_index[gene].first + 1;
				for(int k=0; k<num; k++)
				{
					para_cis_gene[i][j][k] = para_cis_gene[i][j][k] - rate_learner * para_dev_cis_gene[i][j][k];
				}
			}
		}
	}

	//====================== para_snp_cellenv ==========================
	for(int i=0; i<num_cellenv; i++)
	{
		for(long j=0; j<num_snp; j++)
		{
			para_snp_cellenv[i][j] = para_snp_cellenv[i][j] - rate_learner * para_dev_snp_cellenv[i][j];
		}
	}

	//====================== para_cellenv_gene ==========================
	for(int i=0; i<num_etissue; i++)
	{
		string etissue = etissue_list[i];
		for(int j=0; j<num_gene; j++)
		{
			string gene = gene_list[j];
			for(int k=0; k<num_cellenv; k++)
			{
				para_cellenv_gene[i][j][k] = para_cellenv_gene[i][j][k] - rate_learner * para_dev_cellenv_gene[i][j][k];
			}
		}
	}


}




// two (2) major parts for this routine:
// 1. write the procedure to update these parameters in this program (mini-batches gradient; gradient descent)
// 2. stop in some ways, and save the learned parameters to disk (interface; should be in a trivial way); use another routine to check the predictive precision
void optimize()
{
	puts("============== entering the optimization routine...");

	opt_para_init();

	for(int count1=0; count1<iter_learn_out; count1++)  // one count1 is for iteration across all tissues
	{
		for(int count2=0; count2<num_etissue; count2++)  // one count2 is for one tissue
		{
			string etissue = etissue_list[count2];
			int num_esample = eQTL_tissue_rep[etissue].size();

			for(int count3=0; count3<iter_learn_in; count3++)  // one count3 is for a 15-sized mini-batch in current tissue
			{
				int pos_start = (batch_size * count3) % (num_esample);

				printf("now we are working on %d iter_out, %s tissue (%d training samples in), #%d mini-batch (%d batch size, rounding all samples).\n", count1+1, etissue.c_str(), num_esample, count3+1, batch_size);

				forward_backward_prop_batch(etissue, pos_start, num_esample);

				gradient_descent();

			}
		}
	}

	opt_para_release();
	puts("============== leaving the optimization routine...");

}