// the control platform for multi-threading program

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <unordered_map>
#include <string>
#include "global.h"
#include "basic.h"
#include "main.h"				// typedef struct tuple_long
#include "optimization.h"		// optimization local global variables
#include "opt_multi_thread.h"	// declare functions in the current code
#include <math.h>				// /* exp */
#include <cstdlib>				// for mt
#include <pthread.h>			// for mt
#include "opt_subroutine.h"		// regularization and gradient descent
#include "genotype.h"
#include <string.h>



using namespace std;



// local global variables
int num_thread = 8;				// there are at maximum 8 cores in C2B2 cluster, but our main thread doesn't do extensive computation here
pthread_mutex_t mut;			// mutex used by all the threads
int * finish_table;				// finish table for all the samples in this batch ()



// initialize all the parameter space
void package_alloc(package_thread * package_pointer)
{
	//=============== (package_pointer->snp_dosage_list) ===============
	for(int i=0; i<22; i++)
	{
		long num_temp = snp_name_list[i].size();
		float * p = (float *)calloc( num_temp, sizeof(float) );
		(package_pointer->snp_dosage_list)[i] = p;
	}

	//=============== (package_pointer->gene_rpkm_exp) ===============
	(package_pointer->gene_rpkm_exp) = (float *)calloc( num_gene, sizeof(float) );

	//=============== (package_pointer->cellenv_hidden_var) ===============
	(package_pointer->cellenv_hidden_var) = (float *)calloc( num_cellenv, sizeof(float) );

	//=============== (package_pointer->batch_var) ===============
	(package_pointer->batch_var) = (float *)calloc( num_batch, sizeof(float) );

	//=============== (package_pointer->batch_hidden_var) ===============
	(package_pointer->batch_hidden_var) = (float *)calloc( num_batch_hidden, sizeof(float) );

	//=============== (package_pointer->para_dev_cis_gene) =============== TODO: we don't need all the tissues actually
	for(int j=0; j<num_etissue; j++)
	{
		vector<float *> vec;
		(package_pointer->para_dev_cis_gene).push_back(vec);
		for(int i=0; i<num_gene; i++)
		{
			string gene = gene_list[i];
			unordered_map<string, int>::const_iterator got = gene_xymt_rep.find(gene);
			if ( got != gene_xymt_rep.end() )
			{
				float * p = NULL;  // null pointer
				(package_pointer->para_dev_cis_gene)[j].push_back(p);
				continue;
			}
			else
			{
				long first = gene_cis_index[gene].first;  // index
				long second = gene_cis_index[gene].second;  // index
				long amount = second - first + 1;
				float * p = (float *)calloc( amount, sizeof(float) );
				(package_pointer->para_dev_cis_gene)[j].push_back(p);
			}
		}
	}

	//=============== (package_pointer->para_dev_snp_cellenv) ===============
	for(int i=0; i<num_cellenv; i++)
	{
		float * p = (float *)calloc( num_snp, sizeof(float) );
		(package_pointer->para_dev_snp_cellenv).push_back(p);
	}

	//=============== (package_pointer->para_dev_cellenv_gene) =============== TODO: we don't need all the tissues actually
	for(int j=0; j<num_etissue; j++)
	{
		vector<float *> vec;
		(package_pointer->para_dev_cellenv_gene).push_back(vec);
		for(int i=0; i<num_gene; i++)
		{
			float * p = (float *)calloc( num_cellenv, sizeof(float) );
			(package_pointer->para_dev_cellenv_gene)[j].push_back(p);
		}
	}

	//=============== (package_pointer->para_dev_batch_batch_hidden) ===============
	for(int i=0; i<num_batch_hidden; i++)
	{
		float * p = (float *)calloc( num_batch, sizeof(float) );
		(package_pointer->para_dev_batch_batch_hidden).push_back(p);
	}

	//=============== (package_pointer->para_dev_batch_hidden_gene) ===============
	for(int i=0; i<num_gene; i++)
	{
		float * p = (float *)calloc( num_batch_hidden, sizeof(float) );
		(package_pointer->para_dev_batch_hidden_gene).push_back(p);
	}

}



// to release the above space
void package_free(package_thread * package_pointer)
{
	//=============== (package_pointer->snp_dosage_list) ===============
	for(int i=0; i<22; i++)
	{
		free((package_pointer->snp_dosage_list)[i]);
	}

	//=============== (package_pointer->gene_rpkm_exp) ===============
	free(package_pointer->gene_rpkm_exp);

	//=============== (package_pointer->cellenv_hidden_var) ===============
	free(package_pointer->cellenv_hidden_var);

	//=============== (package_pointer->batch_var) ===============
	free(package_pointer->batch_var);

	//=============== (package_pointer->batch_hidden_var) ===============
	free(package_pointer->batch_hidden_var);

	//=============== (package_pointer->para_dev_cis_gene) ===============
	for(int j=0; j<num_etissue; j++)
	{
		for(int i=0; i<num_gene; i++)
		{
			free((package_pointer->para_dev_cis_gene)[j][i]);
		}
	}

	//=============== (package_pointer->para_dev_snp_cellenv) ===============
	for(int i=0; i<num_cellenv; i++)
	{
		free((package_pointer->para_dev_snp_cellenv)[i]);
	}

	//=============== (package_pointer->para_dev_cellenv_gene) ===============
	for(int j=0; j<num_etissue; j++)
	{
		for(int i=0; i<num_gene; i++)
		{
			free((package_pointer->para_dev_cellenv_gene)[j][i]);
		}
	}

	//=============== (package_pointer->para_dev_batch_batch_hidden) ===============
	for(int i=0; i<num_batch_hidden; i++)
	{
		free((package_pointer->para_dev_batch_batch_hidden)[i]);
	}

	//=============== (package_pointer->para_dev_batch_hidden_gene) ===============
	for(int i=0; i<num_gene; i++)
	{
		free((package_pointer->para_dev_batch_hidden_gene)[i]);
	}

}



// this is the working program for each thread
void * WorkPerThread(void * pointer)
{
	package_thread * package_pointer = (package_thread *)pointer;
	// allocate the memory for this thread (this is independent with other threads, but this should be visiable to the main thread -- aggregation)
	package_alloc(package_pointer);

	while(1)
	{
		int count = -1;
		//================ check whether there are still left samples to be processed ================
		// if there are, take into this thread; otherwise, terminate this thread
		pthread_mutex_lock(&mut);
		for(int i=0; i<batch_size; i++)
		{
			if(finish_table[i] == 0)
			{
				count = i;
				finish_table[i] = 1;
				break;
			}
		}
		pthread_mutex_unlock(&mut);

		if(count == -1)  // no left samples for processing
		{
			break;
		}

		//================ work on the current sample ================
		// count will determine the sample in this etissue
		// 1. get all the containers, and forward and backward propagation
		// 2. get the calculated parameters (derivatives) from this round, and fill that into something
		string etissue = package_pointer->etissue;
		int pos_start = package_pointer->pos_start;
		int num_esample = package_pointer->num_esample;
		int id = package_pointer->id;

		int pos = (pos_start + count) % (num_esample);
		string esample = esample_tissue_rep[etissue][pos];
		string individual = sample_to_individual(esample);
		cout << "[$$] thread #" << id+1 << " is working on eSample " << esample << " (#" << count+1 << " out of " << batch_size << ")" << endl;

		//=================================================== init ============================================================
		// get the: 0. esample and individual; 1. genotype; 2. expression data; 3. batch variables
		// to: 1. forward_backward propagation;
		// genotype dosage data
		//cout << "getting the dosage data for individual #" << individual << endl;
		snp_dosage_load(&(package_pointer->snp_dosage_list), individual);  // snp dosage data for one individual across all chromosomes
		// expression rpkm data: eQTL_tissue_rep[etissue][esample]
		//cout << "we have this amount of genes expressed in this individual:" << eQTL_tissue_rep[etissue][esample].size() << endl;
		// and the batch variable for this individual and this sample
		int num_batch_individual = batch_individual[individual].size();
		int index = 0;
		for(int i=0; i<num_batch_individual; i++)
		{
			float value = batch_individual[individual][i];
			(package_pointer->batch_var)[index] = value;
			index++;
		}
		int num_batch_sample = batch_sample[esample].size();
		for(int i=0; i<num_batch_sample; i++)
		{
			float value = batch_sample[esample][i];
			(package_pointer->batch_var)[index] = value;
			index++;
		}

		// no need to clean the derivative containers, as they are initially empty
		forward_backward(etissue,
						&(package_pointer->snp_dosage_list),
						&eQTL_tissue_rep[etissue][esample],
						(package_pointer->gene_rpkm_exp),
						(package_pointer->cellenv_hidden_var),
						(package_pointer->batch_var),
						(package_pointer->batch_hidden_var),
						&(package_pointer->para_dev_cis_gene),
						&(package_pointer->para_dev_cellenv_gene),
						&(package_pointer->para_dev_snp_cellenv),
						&(package_pointer->para_dev_batch_hidden_gene),
						&(package_pointer->para_dev_batch_batch_hidden));

	}

	pthread_exit(NULL);
}




void opt_mt_control(string etissue, int pos_start, int num_esample)
{
	cout << "[@@@] entering the current mini-batch (in multi-treading mode)..." << endl;

	//=============================== multi-threading parameter initialization ===============================
	// allocating all the other threads from here
	pthread_t threads[num_thread];
	void * status;
	pthread_attr_t attr;
	// Initialize and set thread joinable
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	//
	finish_table = (int *)calloc(batch_size, sizeof(int));
	//
	pthread_mutex_init(&mut, NULL);
	memset(&threads, 0, sizeof(threads));


	//=============================== thread local memory allocation ===============================
	package_thread para_array[num_thread];
	for(int i=0; i<num_thread; i++)
	{
		package_thread package;
		package.id = i;
		package.etissue = etissue;
		package.pos_start = pos_start;
		package.num_esample = num_esample;
		//package_alloc(&package);
		para_array[i] = package;
	}


	//=============================== thread initialization ===============================
	for(int i=0; i<num_thread; i++)
	{
		cout << "main() : creating thread#" << i+1 << endl;
		int rc = pthread_create(&threads[i], NULL, WorkPerThread, (void *)&para_array[i]);
		if(rc)
		{
			cout << "Error:unable to create thread," << rc << endl;
			exit(-1);
		}
	}


	//===================== waiting for all the threads to terminate =====================
	// free attribute and wait for the other threads
	pthread_attr_destroy(&attr);
	for(int i=0; i<num_thread; i++)
	{
		int rc = pthread_join(threads[i], &status);
		if(rc)
		{
			cout << "Error:unable to join " << rc << endl;
			exit(-1);
		}
		cout << "Main: completed thread#" << i+1;
		cout << " exiting with status: " << status << endl;
	}


	//===================== merge results, and release space =====================
	//// fill in the real para_dev_xxx space (aggregation)
	aggregation(para_array, etissue);
	//// add the regularization terms into the derivatives
	regularization(etissue);
	/// gradient descent
	gradient_descent(etissue);

	//// free the memory allocated for these threads
	for(int i=0; i<num_thread; i++)
	{
		package_free(&para_array[i]);
	}


	//===================== finish and quit =====================
	free(finish_table);
	cout << "[@@@] finishing the current mini-batch (in multi-treading mode)..." << endl;
}



void aggregation(package_thread * para_array_pointer, string etissue)
{
	cout << "[@@] entering the aggregation routine..." << endl;

	int etissue_index = etissue_index_map[etissue];

	//********************************* aggregation of this mini-batch *****************************************
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
				float temp = 0;
				for(int t=0; t<num_thread; t++)
				{
					temp += ((para_array_pointer+t)->para_dev_cis_gene)[etissue_index][i][k];
				}
				para_dev_cis_gene[etissue_index][i][k] = temp / batch_size;
			}
		}
	}

	// vector<float *> para_dev_snp_cellenv;
	for(int i=0; i<num_cellenv; i++)
	{
		for(long j=0; j<num_snp; j++)
		{
			float temp = 0;
			for(int t=0; t<num_thread; t++)
			{
				temp += ((para_array_pointer+t)->para_dev_snp_cellenv)[i][j];
			}
			para_dev_snp_cellenv[i][j] = temp / batch_size;
		}
	}

	// vector<vector<float *>> para_dev_cellenv_gene;
	for(int i=0; i<num_gene; i++)
	{
		for(int j=0; j<num_cellenv; j++)
		{
			float temp = 0;
			for(int t=0; t<num_thread; t++)
			{
				temp += ((para_array_pointer+t)->para_dev_cellenv_gene)[etissue_index][i][j];
			}
			para_dev_cellenv_gene[etissue_index][i][j] = temp / batch_size;
		}
	}

	// vector<float *> para_dev_batch_batch_hidden;
	for(int i=0; i<num_batch_hidden; i++)
	{
		for(int j=0; j<num_batch; j++)
		{
			float temp = 0;
			for(int t=0; t<num_thread; t++)
			{
				temp += ((para_array_pointer+t)->para_dev_batch_batch_hidden)[i][j];
			}
			para_dev_batch_batch_hidden[i][j] = temp / batch_size;
		}
	}

	// vector<float *> para_dev_batch_hidden_gene;
	for(int i=0; i<num_gene; i++)
	{
		for(int j=0; j<num_batch_hidden; j++)
		{
			float temp = 0;
			for(int t=0; t<num_thread; t++)
			{
				temp += ((para_array_pointer+t)->para_dev_batch_hidden_gene)[i][j];
			}
			para_dev_batch_hidden_gene[i][j] = temp / batch_size;
		}
	}

	cout << "[@@] leaving the aggregation routine..." << endl;
}

