// subroutines of optimization procedure

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
#include "expression.h"
#include "global.h"
#include "main.h"  // typedef struct tuple_long
#include <math.h>       /* exp */
#include "opt_subroutine.h"
#include "optimization.h"



using namespace std;



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
	int etissue_index = etissue_index_map[etissue];

	//******************* initialize all the parameter derivatives (as 0) *******************
	//***************************************************************************************
	// vector<vector<float *>> para_dev_cis_gene;
	for(int i=0; i<num_gene; i++)
	{
		string gene = gene_list[i];
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
				para_dev_cis_gene[etissue_index][i][k] = 0;
			}
		}
	}

	// vector<float *> para_dev_snp_cellenv;
	for(int i=0; i<num_cellenv; i++)
	{
		for(long j=0; j<num_snp; j++)
		{
			para_dev_snp_cellenv[i][j] = 0;
		}
	}

	// vector<vector<float *>> para_dev_cellenv_gene;
	for(int i=0; i<num_gene; i++)
	{
		for(int j=0; j<num_cellenv; j++)
		{
			para_dev_cellenv_gene[etissue_index][i][j] = 0;
		}
	}

	// vector<float *> para_dev_batch_batch_hidden;
	for(int i=0; i<num_batch_hidden; i++)
	{
		for(int j=0; j<num_batch; j++)
		{
			para_dev_batch_batch_hidden[i][j] = 0;
		}
	}

	// vector<float *> para_dev_batch_hidden_gene;
	for(int i=0; i<num_gene; i++)
	{
		for(int j=0; j<num_batch_hidden; j++)
		{
			para_dev_batch_hidden_gene[i][j] = 0;
		}
	}


	//******************* enter the mini-batch *******************
	//************************************************************
	for(int count=0; count<batch_size; count++)
	{
		int pos = (pos_start + count) % (num_esample);
		string esample = esample_tissue_rep[etissue][pos];
		string individual = sample_to_individual(esample);

		//=================================================== init ============================================================
		// get the: 0. esample and individual; 1. genotype; 2. expression data.
		// to: 1. forward_backward propagation

		// genotype dosage data
		cout << "getting the dosage data for individual #" << individual << endl;
		snp_dosage_load(&snp_dosage_list, individual);  // snp dosage data for one individual across all chromosomes
		// expression rpkm data: eQTL_tissue_rep[etissue][esample]
		cout << "we have this amount of genes expressed in this individual:" << eQTL_tissue_rep[etissue][esample].size() << endl;


		//========================================================================
		// two step: forward propagation (get the function values); backward propagation (get the parameter derivatives)
		//========================================================================
		// step#1: ... (cis-; cell env)
		// ****************************** [part1] cis- *********************************
		// for cis-, two issues:
		// 1. if this is a XYMT gene, jump;
		// 2. we use (gene_cis_index[gene].second - gene_cis_index[gene].first + 1) as the length of the cis- parameter array
		for(int i=0; i<num_gene; i++)
		{
			string gene = gene_list[i];
			int chr = gene_tss[gene].chr;
			unordered_map<string, int>::const_iterator got = gene_xymt_rep.find(gene);
			if ( got != gene_xymt_rep.end() )
			{
				gene_rpkm_exp[i] = 0;
			}
			else
			{
				gene_rpkm_exp[i] = 0;
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
				for(long k=0; k<snp_dosage_list[j].size(); k++)
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


		//========================================================================
		// step#2: ... (cis- parameters;  cell env relevant parameters)
		// *********************** [part1] cis- ************************
		// pseudo: (expected rpkm - real rpkm) * genotype
		for(int i=0; i<num_gene; i++)
		{
			string gene = gene_list[i];
			int chr = gene_tss[gene].chr;
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
					int pos = gene_cis_index[gene].first + k;
					float dosage = snp_dosage_list[chr-1][pos];  // dosage at k position
					para_dev_cis_gene[etissue_index][i][k] += (gene_rpkm_exp[i] - eQTL_tissue_rep[etissue][esample][i]) * dosage;
				}
			}
		}

		// ***************** [part2] cell env relevant parameters *****************
		// from cell env to genes
		// pseudo: (expected rpkm - real rpkm) * cell_env
		for(int i=0; i<num_gene; i++)
		{
			string gene = gene_list[i];
			for(int j=0; j<num_cellenv; j++)
			{
				para_dev_cellenv_gene[etissue_index][i][j] += (gene_rpkm_exp[i] - eQTL_tissue_rep[etissue][esample][i]) * cellenv_hidden_var[j];
			}
		}

		// from snp to cell env
		// pseudo: [ \sum w3 * (expected rpkm - real rpkm) ] * g'(w2 * x1) * x1
		for(int i=0; i<num_cellenv; i++)
		{
			long count = 0;
			for(int j=0; j<22; j++)  // across all the chromosomes
			{
				int chr = j+1;
				for(long k=0; k<snp_dosage_list[j].size(); k++)
				{
					float dosage = snp_dosage_list[j][k];
					float temp = 0;
					for(int t=0; t<num_gene; t++)
					{
						temp += para_cellenv_gene[etissue_index][t][i] * (gene_rpkm_exp[t] - eQTL_tissue_rep[etissue][esample][t]);
					}
					temp *= cellenv_hidden_var[i] * ( 1 - cellenv_hidden_var[i] ) * dosage;
					para_dev_snp_cellenv[i][count] += temp;
					count ++;
				}
			}
		}


		// ********************* [part3] linear or non-linear batches *********************
		// from hidden batch to genes
		// pseudo: (expected rpkm - real rpkm) * hidden batch var
		for(int i=0; i<num_gene; i++)
		{
			string gene = gene_list[i];
			for(int j=0; j<num_batch_hidden; j++)
			{
				para_dev_batch_hidden_gene[i][j] += (gene_rpkm_exp[i] - eQTL_tissue_rep[etissue][esample][i]) * batch_hidden_var[j];
			}
		}

		// from original batch to hidden batch
		// pseudo: [ \sum w5 * (expected rpkm - real rpkm) ] * g'(w4 * x2) * x2
		for(int i=0; i<num_batch_hidden; i++)
		{
			for(int j=0; j<num_batch; j++)
			{
				float batch_value = batch_var[j];
				float temp = 0;
				for(int t=0; t<num_gene; t++)
				{
					temp += para_batch_hidden_gene[t][i] * (gene_rpkm_exp[t] - eQTL_tissue_rep[etissue][esample][t]);
				}
				temp *= batch_hidden_var[i] * ( 1 - batch_hidden_var[i] ) * batch_value;
				para_dev_batch_batch_hidden[i][j] += temp;
			}
		}

	}


	//******************* aggregation *******************
	//***************************************************
	// 1. average the derivatives calculated from previous steps;
	// 2. add the derivatives due to regularization;
	// vector<vector<float *>> para_dev_cis_gene;
	for(int i=0; i<num_gene; i++)
	{
		string gene = gene_list[i];
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
				para_dev_cis_gene[etissue_index][i][k] = para_dev_cis_gene[etissue_index][i][k] / batch_size;
			}
		}
	}

	// vector<float *> para_dev_snp_cellenv;
	for(int i=0; i<num_cellenv; i++)
	{
		for(long j=0; j<num_snp; j++)
		{
			para_dev_snp_cellenv[i][j] = para_dev_snp_cellenv[i][j] / batch_size;
		}
	}

	// vector<vector<float *>> para_dev_cellenv_gene;
	for(int i=0; i<num_gene; i++)
	{
		for(int j=0; j<num_cellenv; j++)
		{
			para_dev_cellenv_gene[etissue_index][i][j] = para_dev_cellenv_gene[etissue_index][i][j] / batch_size;
		}
	}

	// vector<float *> para_dev_batch_batch_hidden;
	for(int i=0; i<num_batch_hidden; i++)
	{
		for(int j=0; j<num_batch; j++)
		{
			para_dev_batch_batch_hidden[i][j] = para_dev_batch_batch_hidden[i][j] / batch_size;
		}
	}

	// vector<float *> para_dev_batch_hidden_gene;
	for(int i=0; i<num_gene; i++)
	{
		for(int j=0; j<num_batch_hidden; j++)
		{
			para_dev_batch_hidden_gene[i][j] = para_dev_batch_hidden_gene[i][j] / batch_size;
		}
	}


	// TODO add the regulation relevant items TODO
	// don't forget to add the distance prior to the cis- regulation
















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

	//====================== para_batch_batch_hidden ==========================
	for(int i=0; i<num_batch_hidden; i++)
	{
		for(int j=0; j<num_batch; j++)
		{
			para_batch_batch_hidden[i][j] = para_batch_batch_hidden[i][j] - rate_learner * para_dev_batch_batch_hidden[i][j];
		}
	}

	//====================== para_batch_hidden_gene ==========================
	for(int i=0; i<num_gene; i++)
	{
		for(int j=0; j<num_batch_hidden; j++)
		{
			para_batch_hidden_gene[i][j] = para_batch_hidden_gene[i][j] - rate_learner * para_dev_batch_hidden_gene[i][j];
		}
	}

}