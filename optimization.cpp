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
#include "opt_para_save.h"
#include "opt_debugger.h"
#include "lib_matrix.h"




using namespace std;




//====================================== local global variables ========================================
// these variables are specially designed for this routine -- optimization
// need to initialize some local containers:
array<float *, NUM_CHR> snp_dosage_list;
float * gene_rpkm_exp;  // with length "num_gene"
float * cellenv_hidden_var;  // with length "num_cellenv"
float * batch_var;  // with length "num_batch"
float * batch_hidden_var;  // with length "num_batch_hidden"


// parameter derivative containers:
vector<Matrix_imcomp> cube_para_dev_cis_gene;
Matrix matrix_para_dev_snp_cellenv;
vector<Matrix> cube_para_dev_cellenv_gene;
Matrix matrix_para_dev_batch_batch_hidden;
Matrix matrix_para_dev_batch_hidden_gene;


// some assistant components:
// the prior number for each un-pruned snp for regularization (from pruned snps and chromatin states); per etissue, per chromosome, for each snp
// TODO: we also still need to integrate distance prior later on with the following prior information
vector<vector<vector<float>>> prior_tissue_vector;
// pairwise phylogenetic distance between etissues
vector<vector<float>> tissue_hierarchical_pairwise;


// learning control parameters:
int iter_learn_out = 1;  // iteration across all tissues
int iter_learn_in = 100;  // iteration across all samples from one tissue
int batch_size = 20;  // better be 20

// test different learning rate
//float rate_learner = 1.0;  // the learning rate; this doesn't work
//float rate_learner = 0.1;  // the learning rate; this doesn't work
//float rate_learner = 0.01;  // the learning rate; this doesn't work
//float rate_learner = 0.001;  // the learning rate; works!!!; bench#3
//float rate_learner = 0.0001;  // the learning rate; works!!!; bench#4
float rate_learner = 0.00001;  // the learning rate; works!!!; bench#5
//float rate_learner = 0.000001;  // the learning rate


//======================================================================================================







// TODO: we need to think more about this, say, how to use the data from other epigenomics projects
// load all the cis- snp prior information (tissue specific) from prepared file outside
// fill in the following: vector<vector<vector<float>>> prior_tissue_vector
void opt_snp_prior_load()
{

	for(int i=0; i<num_etissue; i++)
	{
		vector<vector<float>> matrix;
		prior_tissue_vector.push_back(matrix);
	}


	/*
	// TODO: there should always be prior information in the repo
	// TODO: in the simulating data, we don't have this prior, so temporarily stop this

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
		for(i=0; i<NUM_CHR; i++)
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
	*/

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
	fclose(file_in);

}




void opt_para_init()
{
	puts("opt_para_init..");

	//=============== snp_dosage_list ===============
	for(int i=0; i<NUM_CHR; i++)
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



	//=============== cube_para_dev_cis_gene ===============
	for(int j=0; j<num_etissue; j++)
	{
		Matrix_imcomp matrix_imcomp;
		matrix_imcomp.init(num_gene);
		for(long int i=0; i<num_gene; i++)
		{
			string gene = gene_list[i];
			unordered_map<string, int>::const_iterator got = gene_xymt_rep.find(gene);
			if ( got != gene_xymt_rep.end() )
			{
				continue;
			}
			else
			{
				long int first = gene_cis_index[gene].first;  // index
				long int second = gene_cis_index[gene].second;  // index
				long int amount = second - first + 1;
				matrix_imcomp.init_element(i, amount + 1);

				// assing the chr and the tss:
				matrix_imcomp.init_assign_chr(i, gene_tss[gene].chr);
				matrix_imcomp.init_assign_sst(i, gene_tss[gene].tss);
			}
		}
		cube_para_dev_cis_gene.push_back(matrix_imcomp);
	}

	//=============== matrix_para_dev_snp_cellenv ===============
	matrix_para_dev_snp_cellenv.init(num_cellenv, num_snp + 1);		// we do have the intercept term here

	//=============== cube_para_dev_cellenv_gene ===============
	for(int j=0; j<num_etissue; j++)
	{
		Matrix matrix;
		matrix.init(num_gene, num_cellenv + 1);						// we do have the intercept term here
		cube_para_dev_cellenv_gene.push_back(matrix);
	}

	//=============== matrix_para_dev_batch_batch_hidden ===============
	matrix_para_dev_batch_batch_hidden.init(num_batch_hidden, num_batch + 1);

	//=============== matrix_para_dev_batch_hidden_gene ===============
	matrix_para_dev_batch_hidden_gene.init(num_gene, num_batch_hidden + 1);


}




void opt_para_release()
{
	//=============== snp_dosage_list ===============
	for(int i=0; i<NUM_CHR; i++)
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



	//=============== cube_para_dev_cis_gene ===============
	for(int j=0; j<num_etissue; j++)
	{
		cube_para_dev_cis_gene[j].release();
	}


	//=============== matrix_para_dev_snp_cellenv ===============
	matrix_para_dev_snp_cellenv.release();


	//=============== cube_para_dev_cellenv_gene ===============
	for(int j=0; j<num_etissue; j++)
	{
		cube_para_dev_cellenv_gene[j].release();
	}


	//=============== matrix_para_dev_batch_batch_hidden ===============
	matrix_para_dev_batch_batch_hidden.release();


	//=============== matrix_para_dev_batch_hidden_gene ===============
	matrix_para_dev_batch_hidden_gene.release();


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




				// DEBUG
				// check nan after this mini-batch
				int flag = para_check_nan(etissue);
				if(flag == 1)
				{
					//
					cout << "we get nan..." << endl;
					cout << count3 << endl;
					break;
				}




			}
			// leaving this etissue

			//DEBUG
			break;  // won't consider other tissues

		}
		//
		// whenever we finish one iteration across all tissues, we should save the learned parameters
		//
		//para_inter_save(count1);
		//
		//

	}



	opt_para_release();
	puts("============== leaving the optimization routine...");
}

