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



	//===================================== Regularization in Regression =====================================
	//========================================================================================================
	// there are several classes of prior knowledge that we need to consider
	// 1. sparsity of cis- regulation, accompanied by ridge regression, achieved by elastic-net tuned by the prior number, and the distance prior
	// 2. sparsity (LASSO) for the coefficients from cell env to expression (with the assumption that one gene is only affected by several handful cell env)
	// 3.1. hierarchical regularization tuned by the learned tissue hierarchy
	// 3.2. or we can simply use group LASSO to encourage the tissue consistency
	// 4. penalize the batch variables hashly (from batch variables to batch_hidden, and from batch_hidden to genes)

	//===================================== part#0 =====================================
	// initialize some learning parameters

	// define the sigma here that may be used by L1 regularization
	float sigma = 0.0001;

	// regularization strength lambda:
	// path#1: add the lambda_{LASSO} and lambda_{ridge} for the cis- regularization
	// path#2: add the lambda for cellenv-gene regularization
	// path#3: and the lambda for batch-batch_hidden and batch_hidden-gene
	float lambda_lasso = 1.0;
	float lambda_ridge = 1.0;
	float lambda_cellenv_gene = 1.0;
	float lambda_batch_batch_hidden = 1.0;
	float lambda_batch_hidden_gene = 1.0;

	//===================================== part#1 =====================================
	// 1. sparsity of cis- regulation, accompanied by ridge regression, achieved by elastic-net tuned by the prior number, and the distance prior
	// TODO: not yet integrated the distance prior information
	for(int i=0; i<num_gene; i++)
	{
		string gene = gene_list[i];
		int chr = gene_tss[gene].chr;
		long int index_start = gene_cis_index[gene].first;
		long int index_end = gene_cis_index[gene].second;
		long int amount = index_end - index_start;
		for(int j=0; j<=amount; j++)
		{
			long int pos = j+index_start;  // the pos in snp list

			/// the value of current cis- beta:
			float beta = para_cis_gene[etissue_index][i][j];

			/// the prior that we need (if there is) for tuning the relative strength of L1 and L2 regularization:
			float prior = 0;
			unordered_map<string, vector<vector<float>>>::const_iterator got = prior_tissue_rep.find(etissue);
			if( got != prior_tissue_rep.end() )
			{
				prior = prior_tissue_rep[etissue][chr-1][pos];
			}
			else
			{
				prior = 1;
			}
			float alpha = 1 / ( 1 + exp(-(prior-1)) );

			/// the derivative of the beta:
			float derivative1 = beta / sqrt (beta * beta + sigma);  // this is an approximation of the LASSO regularization
			float derivative2 = 2 * beta;  // L2 regularization item is differentiable

			/// and the value of its derivative should be added with that derivative item from regularization:
			para_dev_cis_gene[etissue_index][i][j] += lambda_lasso * (1-alpha) * derivative1 + lambda_ridge * alpha * derivative2;
		}
	}

	//===================================== part#2 =====================================
	// 2. sparsity (LASSO) for the coefficients from cell env to expression (with the assumption that one gene is only affected by several handful cell env)
	for(int i=0; i<num_gene; i++)
	{
		string gene = gene_list[i];
		for(int j=0; j<num_cellenv; j++)
		{
			/// the value of current cellenv beta:
			float beta = para_cellenv_gene[etissue_index][i][j];
			/// the derivative of the beta:
			float derivative = beta / sqrt (beta * beta + sigma);  // this is an approximation of the LASSO regularization
			/// and the value of its derivative should be added with that derivative item from regularization:
			para_dev_cellenv_gene[etissue_index][i][j] +=  lambda_cellenv_gene * derivative;
		}
	}

	//===================================== part#3 =====================================
	// 3.2. or we can simply use group LASSO to encourage the tissue consistency
	// TODO: as this part is too un-stable (group LASSO, or hierarchical regularization), we now don't use them






	//===================================== part#4 =====================================
	// 4. penalize the batch variables hashly (from batch variables to batch_hidden, and from batch_hidden to genes)
	// from batch to batch_hidden:
	for(int i=0; i<num_batch_hidden; i++)
	{
		for(int j=0; j<num_batch; j++)
		{
			/// the value of current batch beta:
			float beta = para_batch_batch_hidden[i][j];
			/// the derivative of the beta:
			float derivative = beta / sqrt (beta * beta + sigma);  // this is an approximation of the LASSO regularization
			/// and the value of its derivative should be added with that derivative item from regularization:
			para_dev_batch_batch_hidden[i][j] += lambda_batch_batch_hidden * derivative;
		}
	}
	// from batch_hidden to gene:
	for(int i=0; i<num_gene; i++)
	{
		for(int j=0; j<num_batch_hidden; j++)
		{
			/// the value of current batch beta:
			float beta = para_batch_hidden_gene[i][j];
			/// the derivative of the beta:
			float derivative = beta / sqrt (beta * beta + sigma);  // this is an approximation of the LASSO regularization
			/// and the value of its derivative should be added with that derivative item from regularization:
			para_dev_batch_hidden_gene[i][j] += lambda_batch_hidden_gene * derivative;
		}
	}

	//===================================== end of Regularization ============================================
	//========================================================================================================

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