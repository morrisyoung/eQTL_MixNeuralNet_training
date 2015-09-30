// calculate the expression level on testing dataset, from genotype data, with the learned parameters

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
#include <math.h>       /* exp */
// sub-routines:
#include "global.h"
#include "genotype.h"
#include "expression.h"
#include "batch.h"
#include "basic.h"
#include "test_predict.h"
#include "test_main.h"



using namespace std;




//====================================== local global variables ========================================
array<float *, 22> snp_dosage_list;
float * gene_rpkm_exp;  // with length "num_gene"
float * cellenv_hidden_var;  // with length "num_cellenv"
float * batch_var;  // with length "num_batch"
float * batch_hidden_var;  // with length "num_batch_hidden"
//======================================================================================================




void opt_para_init()
{
	//=============== snp_dosage_list ===============
	for(int i=0; i<22; i++)
	{
		long num_temp = snp_name_list[i].size();
		float * p = (float *)calloc( num_temp, sizeof(float) );
		snp_dosage_list[i] = p;
	}

	//=============== gene_rpkm_exp ===============
	gene_rpkm_exp = (float *)calloc( num_gene, sizeof(float) );

	//=============== cellenv_hidden_var ===============
	cellenv_hidden_var = (float *)calloc( num_cellenv, sizeof(float) );

	//=============== batch_var ===============
	batch_var = (float *)calloc( num_batch, sizeof(float) );

	//=============== batch_hidden_var ===============
	batch_hidden_var = (float *)calloc( num_batch_hidden, sizeof(float) );

}



void opt_para_release()
{
	//=============== snp_dosage_list ===============
	for(int i=0; i<22; i++)
	{
		free(snp_dosage_list[i]);
	}

	//=============== gene_rpkm_exp ===============
	free(gene_rpkm_exp);

	//=============== cellenv_hidden_var ===============
	free(cellenv_hidden_var);

	//=============== batch_var ===============
	free(batch_var);

	//=============== batch_hidden_var ===============
	free(batch_hidden_var);

}





void predict()
{
	opt_para_init();


	// for sample in eQTL_samples, get the individual, then get the predicted rpkm, then fill it into eQTL_tissue_rep_predict
	int count = 0;
	for(auto it = eQTL_samples.begin(); it != eQTL_samples.end(); ++it)
	{
		string esample = it->first;
		string individual = sample_to_individual(esample);
		string etissue = it->second;
		int etissue_index = etissue_index_map[etissue];

		cout << "current sample #" << count+1 << ": " << esample << endl;
		count++;
		//=================================================== init ============================================================
		// get the: 0. esample and individual; 1. genotype; (2. expression data;) 3. batch variables
		// to: 1. forward propagation;
		// genotype dosage data
		//cout << "getting the dosage data for individual #" << individual << endl;
		snp_dosage_load(&snp_dosage_list, individual);  // snp dosage data for one individual across all chromosomes
		// expression rpkm data: eQTL_tissue_rep[etissue][esample]
		//cout << "we have this amount of genes expressed in this individual:" << eQTL_tissue_rep[etissue][esample].size() << endl;
		// and the batch variable for this individual and this sample
		int num_batch_individual = batch_individual[individual].size();
		int index = 0;
		for(int i=0; i<num_batch_individual; i++)
		{
			float value = batch_individual[individual][i];
			batch_var[index] = value;
			index++;
		}
		int num_batch_sample = batch_sample[esample].size();
		for(int i=0; i<num_batch_sample; i++)
		{
			float value = batch_sample[esample][i];
			batch_var[index] = value;
			index++;
		}

		// ****************************** [part1] cis- *********************************
		// for cis-, two issues:
		// 1. if this is a XYMT gene, jump;
		// 2. we use (gene_cis_index[gene].second - gene_cis_index[gene].first + 1) as the length of the cis- parameter array
		for(int i=0; i<num_gene; i++)
		{
			string gene = gene_list[i];
			unordered_map<string, int>::const_iterator got = gene_xymt_rep.find(gene);
			if ( got != gene_xymt_rep.end() )
			{
				gene_rpkm_exp[i] = 0;
			}
			else
			{
				gene_rpkm_exp[i] = 0;
				int chr = gene_tss[gene].chr;
				int num = gene_cis_index[gene].second - gene_cis_index[gene].first + 1;
				for(int k=0; k<num; k++)
				{
					int pos = gene_cis_index[gene].first + k;
					float dosage = snp_dosage_list[chr-1][pos];  // dosage at k position
					gene_rpkm_exp[i] += dosage * para_cis_gene[etissue_index][i][k];
				}
			}
		}

		// ********************* [part2] cell env relevant parameters *********************
		// from snp to cell env variables
		for(int i=0; i<num_cellenv; i++)
		{
			cellenv_hidden_var[i] = 0;
			long count = 0;
			for(int j=0; j<22; j++)  // across all the chromosomes
			{
				int chr = j+1;
				for(long k=0; k<snp_name_list[j].size(); k++)
				{
					float dosage = snp_dosage_list[j][k];
					cellenv_hidden_var[i] += dosage * para_snp_cellenv[i][count];
					count ++;
				}
			}
		}
		//$$$$$$$$$$$ perform the activation function here (logistic or something else) $$$$$$$$$$$$
		for(int i=0; i<num_cellenv; i++)
		{
			cellenv_hidden_var[i] = 1 / ( 1 + exp( - cellenv_hidden_var[i] ));
		}
		// from cell env variables to genes
		for(int i=0; i<num_gene; i++)
		{
			string gene = gene_list[i];
			for(int j=0; j<num_cellenv; j++)
			{
				gene_rpkm_exp[i] += para_cellenv_gene[etissue_index][i][j] * cellenv_hidden_var[j];
			}
		}

		// ********************* [part3] linear or non-linear batches *********************
		// from original batch to hidden batch
		for(int i=0; i<num_batch_hidden; i++)
		{
			batch_hidden_var[i] = 0;
			for(int j=0; j<num_batch; j++)
			{
				batch_hidden_var[i] += batch_var[j] * para_batch_batch_hidden[i][j];
			}
		}
		//$$$$$$$$$$$ perform the activation function here (logistic or something else) $$$$$$$$$$$$
		for(int i=0; i<num_batch_hidden; i++)
		{
			batch_hidden_var[i] = 1 / ( 1 + exp( - batch_hidden_var[i] ));
		}
		// from hidden batch to genes
		for(int i=0; i<num_gene; i++)
		{
			string gene = gene_list[i];
			for(int j=0; j<num_batch_hidden; j++)
			{
				gene_rpkm_exp[i] += para_batch_hidden_gene[i][j] * batch_hidden_var[j];
			}
		}




		// fill in the eQTL_tissue_rep_predict[etissue][esample] for this esample
		for(int i=0; i<num_gene; i++)
		{
			eQTL_tissue_rep_predict[etissue][esample][i] = gene_rpkm_exp[i];
		}

		// leaving this esample
	}


	opt_para_release();	
}

