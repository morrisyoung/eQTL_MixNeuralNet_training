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
#include "opt_multi_thread.h"



using namespace std;



//====================================== local global variables ========================================
// these variables are specially designed for this routine -- optimization
// need to initialize some local containers:
array<float *, 22> snp_dosage_list;
float * gene_rpkm_exp;  // with length "num_gene"
float * cellenv_hidden_var;  // with length "num_cellenv"
float * batch_var;  // with length "num_batch"
float * batch_hidden_var;  // with length "num_batch_hidden"


// parameter derivative containers:
vector<vector<float *>> para_dev_cis_gene;
vector<float *> para_dev_snp_cellenv;
vector<vector<float *>> para_dev_cellenv_gene;
vector<float *> para_dev_batch_batch_hidden;
vector<float *> para_dev_batch_hidden_gene;


// some assistant components:
// the prior number for each un-pruned snp for regularization (from pruned snps and chromatin states); per etissue, per chromosome, for each snp
// TODO: we also still need to integrate distance prior later on with the following prior information
unordered_map<string, vector<vector<float>>> prior_tissue_rep;
// pairwise phylogenetic distance between etissues
vector<vector<float>> tissue_hierarchical_pairwise;


// learning control parameters:
int iter_learn_out = 1;  // iteration across all tissues
int iter_learn_in = 1;  // iteration across all samples from one tissue
int batch_size = 2;  // better be 20
int rate_learner = 1;  // the learning rate

//======================================================================================================



// load all the cis- snp prior information (tissue specific) from prepared file outside
// fill in the following: unordered_map<string, vector<vector<float>>> prior_tissue_rep
void opt_snp_prior_load()
{
	// get the eTissue-index map
	unordered_map<string, string> index_map;  // for temporary usage
	char filename[100] = "../prior.score.unpruned/prior_tissue_index.txt";
	FILE * file_in = fopen(filename, "r");
	if(file_in == NULL)
	{
		fputs("File error\n", stderr); exit (1);
	}
	int input_length = 1000;
	char input[input_length];
	while(fgets(input, input_length, file_in) != NULL)
	{
		trim(input);

		const char * sep = "\t";
		char * p;
		p = strtok(input, sep);
		string eTissue = p;

		int count = 0;
		while(p)
		{
			count++;
			if(count == 1)  // this is the eTissue
			{
				p = strtok(NULL, sep);
				continue;
			}
			if(count == 2)  // this is the index
			{
				string index = p;
				index_map[eTissue] = index;
				break;
			}
		}
	}
	fclose (file_in);

	// get the prior score for each eTissue, on all chromosomes
	for( auto it = index_map.begin(); it != index_map.end(); ++it )
	{
		string eTissue = it->first;
		string index = it->second;
		vector<vector<float>> vec;
		prior_tissue_rep[eTissue] = vec;

		int i;
		for(i=0; i<22; i++)
		{
			int chr = i+1;
			vector<float> vec;
			prior_tissue_rep[eTissue].push_back(vec);

			//======== get all SNPs with their snp_info (count, position) ========
			char filename[100] = "../prior.score.unpruned/etissue";
			char temp[10];
			StrToCharSeq(temp, index);
			strcat(filename, temp);
			strcat(filename, "/chr");
			sprintf(temp, "%d", chr);
			strcat(filename, temp);
			strcat(filename, ".score");
			//puts("the current file worked on is: ");
			//puts(filename);

			FILE * file_in = fopen(filename, "r");
			if(file_in == NULL)
			{
				fputs("File error\n", stderr); exit (1);
			}

			int input_length = 100;
			char input[input_length];
			while(fgets(input, input_length, file_in) != NULL)
			{
				trim(input);

				float prior = stof(input);
				prior_tissue_rep[eTissue][i].push_back(prior);
			}
			fclose(file_in);
			//======================================
		}
	}

}



// load the pairwise tissue hierarchy from prepared file outside
// TODO: maybe we should check whether this makes the results better
void opt_tissue_hierarchy_load()
{
	// target: vector<vector<float>> tissue_hierarchical_pairwise;
	// init
	for(int i=0; i<num_etissue; i++)
	{
		vector<float> vec;
		for(int j=0; j<num_etissue; j++)
		{
			vec.push_back(0);
		}
		tissue_hierarchical_pairwise.push_back(vec);
	}

	// load from data source
	char filename[100] = "../tissue_hierarchy_normalized.txt";
	FILE * file_in = fopen(filename, "r");
	if(file_in == NULL)
	{
		fputs("File error\n", stderr); exit (1);
	}
	int input_length = 100000;
	char input[input_length];
	while(fgets(input, input_length, file_in) != NULL)
	{
		trim(input);

		const char * sep = "\t";
		char * p;
		p = strtok(input, sep);
		string eTissue1 = p;
		int index1 = etissue_index_map[eTissue1];
		int index2 = 0;

		int count = 0;
		while(p)
		{
			count++;
			if(count == 1)  // this is the eTissue1
			{
				p = strtok(NULL, sep);
				continue;
			}
			if(count == 2)  // this is the eTissue2
			{
				string eTissue2 = p;
				int index2 = etissue_index_map[eTissue2];

				p = strtok(NULL, sep);
				continue;
			}
			if(count == 3)
			{
				float dist = stof(p);
				tissue_hierarchical_pairwise[index1][index2] = dist;
				tissue_hierarchical_pairwise[index2][index1] = dist;
				break;
			}
		}

	}
	fclose (file_in);

}



void opt_para_init()
{
	puts("opt_para_init..");

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


	//=============== para_dev_snp_cellenv ===============
	for(int i=0; i<num_cellenv; i++)
	{
		float * p = (float *)calloc( num_snp, sizeof(float) );
		para_dev_snp_cellenv.push_back(p);
	}

	//=============== para_dev_cellenv_gene ===============
	for(int j=0; j<num_etissue; j++)
	{
		vector<float *> vec;
		para_dev_cellenv_gene.push_back(vec);
		for(int i=0; i<num_gene; i++)
		{
			float * p = (float *)calloc( num_cellenv, sizeof(float) );
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
				float * p = (float *)calloc( amount, sizeof(float) );
				para_dev_cis_gene[j].push_back(p);
			}
		}
	}

	//=============== para_dev_batch_batch_hidden ===============
	for(int i=0; i<num_batch_hidden; i++)
	{
		float * p = (float *)calloc( num_batch, sizeof(float) );
		para_dev_batch_batch_hidden.push_back(p);
	}

	//=============== para_dev_batch_hidden_gene ===============
	for(int i=0; i<num_gene; i++)
	{
		float * p = (float *)calloc( num_batch_hidden, sizeof(float) );
		para_dev_batch_hidden_gene.push_back(p);
	}

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



//function: mini-batches gradient; gradient descent
void optimize()
{
	puts("============== entering the optimization routine...");
	puts("[xx] loading the tissue hierarchy...");
	opt_tissue_hierarchy_load();
	puts("[xx] initializing the parameter space in this optimization routine...");
	opt_para_init();
	puts("[xx] loading the prior information for cis- snps...");
	opt_snp_prior_load();


	for(int count1=0; count1<iter_learn_out; count1++)  // one count1 is for iteration across all tissues
	{
		for(int count2=0; count2<num_etissue; count2++)  // one count2 is for one tissue
		{
			string etissue = etissue_list[count2];
			int num_esample = eQTL_tissue_rep[etissue].size();

			for(int count3=0; count3<iter_learn_in; count3++)  // one count3 is for a batch_size mini-batch in current tissue
			{
				int pos_start = (batch_size * count3) % (num_esample);
				printf("[@@@] now we are working on %d iter_out (%d total), eTissue #%d (%d total) -- %s (%d training samples in), #%d mini-batch (%d batch size, rounding all samples).\n", count1+1, iter_learn_out, count2+1, num_etissue, etissue.c_str(), num_esample, count3+1, batch_size);
				if(MULTI_THREAD == 0)  // normal sequential program
				{
					forward_backward_prop_batch(etissue, pos_start, num_esample);
				}
				else  // multi-threading program
				{
					opt_mt_control(etissue, pos_start, num_esample);
				}
				// leaving this mini-batch
			}
			// leaving this etissue
		}
		//
	}


	opt_para_release();
	puts("============== leaving the optimization routine...");
}

