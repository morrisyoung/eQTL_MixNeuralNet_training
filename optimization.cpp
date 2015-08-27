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
#include "expression.h"
#include "optimization.h"
#include "global.h"
#include "main.h"  // typedef struct tuple_long
#include <math.h>       /* exp */
#include "opt_subroutine.h"



using namespace std;



//====================================== local global variables ========================================
// these variables are specially designed for this routine -- optimization
// need to initialize some local containers:
array<vector<float>, 22> snp_dosage_list;
vector<float> gene_rpkm_exp;  // with length "num_gene"
vector<float> cellenv_hidden_var;  // with length "num_cellenv"
vector<float> batch_var;  // with length "num_batch"
vector<float> batch_hidden_var;  // with length "num_batch_hidden"


// parameter derivative containers:
vector<vector<float *>> para_dev_cis_gene;
vector<float *> para_dev_snp_cellenv;
vector<vector<float *>> para_dev_cellenv_gene;
vector<float *> para_dev_batch_batch_hidden;
vector<float *> para_dev_batch_hidden_gene;


// some assistant components:
// the prior number for each un-pruned snp for regularization (from pruned snps and chromatin states)
// per etissue, per chromosome, for each snp
// we still need to integrate distance prior later on with the following prior information
vector<vector<vector<float>>> snp_prior_list;  // the prior number for each un-pruned snp for regularization (from pruned snps and chromatin states)
// pairwise phylogenetic distance between etissues
vector<vector<float>> tissue_hierarchical_pairwise;


// learning control parameters:
int iter_learn_out = 5;  // iteration across all tissues
int iter_learn_in = 20;  // iteration across all samples from one tissue
int batch_size = 15;
int rate_learner = 1;  // the learning rate
//
//
//

//======================================================================================================



// load all the cis- snp prior information (tissue specific) from prepared file outside
void opt_snp_prior_load()
{


}



// load the pairwise tissue hierarchy from prepared file outside
void opt_tissue_hierarch_load()
{


}



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

	//=============== batch_var ===============
	for(int i=0; i<num_batch; i++)
	{
		batch_var.push_back(0);
	}

	//=============== batch_hidden_var ===============
	for(int i=0; i<num_batch_hidden; i++)
	{
		batch_hidden_var.push_back(0);
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

	//=============== para_dev_batch_batch_hidden ===============
	for(int i=0; i<num_batch_hidden; i++)
	{
		float * p = (float *)malloc( sizeof(float) * num_batch );
		para_dev_batch_batch_hidden.push_back(p);
	}

	//=============== para_dev_batch_hidden_gene ===============
	for(int i=0; i<num_gene; i++)
	{
		float * p = (float *)malloc( sizeof(float) * num_batch_hidden );
		para_dev_batch_hidden_gene.push_back(p);
	}

	//=============== snp_prior_list ===============
	//vector<vector<vector<float>>> snp_prior_list;  // the prior number for each un-pruned snp for regularization (from pruned snps and chromatin states)
	// TODO: or fill this table with data from file
	for(int i=0; i<num_etissue; i++)
	{
		string etissue = etissue_list[i];
		vector<vector<float>> vec;
		snp_prior_list.push_back(vec);
		for(int j=0; j<22; j++)
		{
			int chr = j+1;
			vector<float> vec1;
			snp_prior_list[i].push_back(vec1);
			for(int k=0; k<snp_name_list.size(); k++)
			{
				snp_prior_list[i][j].push_back(1.0);  // 1.0 is used to cover the memory temporarily
			}
		}
	}

	//=============== tissue_hierarchical_pairwise ===============
	//vector<vector<float>> tissue_hierarchical_pairwise;
	// TODO: or fill this table with data from file
	for(int i=0; i<num_etissue; i++)
	{
		string etissue1 = etissue_list[i];
		vector<float> vec;
		tissue_hierarchical_pairwise.push_back(vec);
		for(int j=0; j<num_etissue; j++)
		{
			string etissue2 = etissue_list[j];
			tissue_hierarchical_pairwise[i].push_back(1.0);
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

	//=============== para_dev_batch_batch_hidden ===============
	for(int i=0; i<num_batch_hidden; i++)
	{
		free(para_dev_batch_batch_hidden[i]);
	}

	//=============== para_dev_batch_hidden_gene ===============
	for(int i=0; i<num_gene; i++)
	{
		free(para_dev_batch_hidden_gene[i]);
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