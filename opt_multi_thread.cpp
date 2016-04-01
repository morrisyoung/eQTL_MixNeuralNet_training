// the control platform for multi-threading program
// this program will call sub routines in "opt_subroutine.cpp"

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
#include "expression.h"
#include "lib_matrix.h"



using namespace std;



// local global variables
//int num_thread = 12;				// there are at maximum 8 cores in C2B2 cluster, but our main thread doesn't do extensive computation here
// instead of declaring here, I will declare "num_thread" in the real global scope
pthread_mutex_t mut;			// mutex used by all the threads
int * finish_table;				// finish table for all the samples in this batch ()




// initialize all the parameter space
void package_alloc(package_dev * package_pointer)
{
	//=============== (package_pointer->snp_dosage_list) ===============
	for(int i=0; i<NUM_CHR; i++)
	{
		long num_temp = snp_name_list[i].size();						// TODO: maybe it's better that we can hide some global variables
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


	//=============== (package_pointer->matrix_imcomp_para_dev_cis_gene) ===============
	(package_pointer->matrix_imcomp_para_dev_cis_gene).init(num_gene);
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
			long int first = gene_cis_index[gene].first;
			long int second = gene_cis_index[gene].second;
			long int amount = second - first + 1;
			(package_pointer->matrix_imcomp_para_dev_cis_gene).init_element(i, amount + 1);			// we do have the intercept term

			// assing the chr and the tss:
			(package_pointer->matrix_imcomp_para_dev_cis_gene).init_assign_chr(i, gene_tss[gene].chr);
			//(package_pointer->matrix_imcomp_para_dev_cis_gene).init_assign_sst(i, gene_tss[gene].tss);	// here is a bug!!! sst != tss
			(package_pointer->matrix_imcomp_para_dev_cis_gene).init_assign_sst(i, gene_cis_index[gene].first);	// Here is a BUG: sst != tss
		}
	}

	//=============== (package_pointer->matrix_para_dev_snp_cellenv) ===============
	(package_pointer->matrix_para_dev_snp_cellenv).init(num_cellenv, num_snp + 1);					// we do have the intercept term here

	//=============== (package_pointer->matrix_para_dev_cellenv_gene) ===============
	(package_pointer->matrix_para_dev_cellenv_gene).init(num_gene, num_cellenv + 1);				// we do have the intercept term here

	//=============== (package_pointer->matrix_para_dev_batch_batch_hidden) ===============
	(package_pointer->matrix_para_dev_batch_batch_hidden).init(num_batch_hidden, num_batch + 1);	// we do have the intercept term here

	//=============== (package_pointer->matrix_para_dev_batch_hidden_gene) ===============
	(package_pointer->matrix_para_dev_batch_hidden_gene).init(num_gene, num_batch_hidden + 1);		// we do have the intercept term here

}



// to release the above space
void package_free(package_dev * package_pointer)
{
	//=============== (package_pointer->snp_dosage_list) ===============
	for(int i=0; i<NUM_CHR; i++)
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



	//=============== (package_pointer->matrix_imcomp_para_dev_cis_gene) ===============
	(package_pointer->matrix_imcomp_para_dev_cis_gene).release();

	//=============== (package_pointer->matrix_para_dev_snp_cellenv) ===============
	(package_pointer->matrix_para_dev_snp_cellenv).release();

	//=============== (package_pointer->matrix_para_dev_cellenv_gene) ===============
	(package_pointer->matrix_para_dev_cellenv_gene).release();

	//=============== (package_pointer->matrix_para_dev_batch_batch_hidden) ===============
	(package_pointer->matrix_para_dev_batch_batch_hidden).release();

	//=============== (package_pointer->matrix_para_dev_batch_hidden_gene) ===============
	(package_pointer->matrix_para_dev_batch_hidden_gene).release();

}



// this is the working program for each thread
void * WorkPerThread(void * pointer)
{
	package_dev * package_pointer = (package_dev *)pointer;
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
		// 2. get the calculated parameters (derivatives) from this round, and fill that into para_dev repo
		string etissue = package_pointer->etissue;
		int etissue_index = package_pointer->etissue_index;
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
		forward_backward(etissue_index,
						&(package_pointer->snp_dosage_list),
						&eQTL_tissue_rep[etissue][esample],

						(package_pointer->gene_rpkm_exp),
						(package_pointer->cellenv_hidden_var),
						(package_pointer->batch_var),
						(package_pointer->batch_hidden_var),

						(package_pointer->matrix_imcomp_para_dev_cis_gene),
						(package_pointer->matrix_para_dev_snp_cellenv),
						(package_pointer->matrix_para_dev_cellenv_gene),
						(package_pointer->matrix_para_dev_batch_batch_hidden),
						(package_pointer->matrix_para_dev_batch_hidden_gene)
						);

	}

	pthread_exit(NULL);
}




void opt_mt_control(string etissue, int pos_start, int num_esample)
{
	cout << "[@@@] entering the current mini-batch (in multi-treading mode)..." << endl;

	int etissue_index = etissue_index_map[etissue];				// why we transform this information back and forth???

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

	//=============================== thread local memory preparation (not yet allocation heap) ===============================
	package_dev para_array[num_thread];
	for(int i=0; i<num_thread; i++)
	{
		package_dev package;
		package.id = i;
		package.etissue = etissue;
		package.etissue_index = etissue_index;
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

	//===================== waiting for all the threads to terminate, and aggregate =====================
	// free attribute and wait for the other threads
	aggregation_init(etissue_index);
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
		cout << " exiting with status: " << status;
		cout << ", aggregate it's results..." << endl;

		// aggregate the results from this current thread
		aggregation_add(&para_array[i], etissue_index);
	}
	aggregation_ave(etissue_index);


	//
	//
	// can debug the derivatives from here
	//
	//


	//===================== regularize, gradient descent, and release space =====================
	//// add the regularization terms into the derivatives
	regularization(etissue_index);
	/// gradient descent
	gradient_descent(etissue_index);

	//// free the memory allocated for these threads
	for(int i=0; i<num_thread; i++)
	{
		package_free(&para_array[i]);
	}


	//===================== finish and quit =====================
	free(finish_table);
	cout << "[@@@] finishing the current mini-batch (in multi-treading mode)..." << endl;
}




// set the deravetives to 0
void aggregation_init(int etissue_index)
{
	//******************* initialize all the parameter derivatives (as 0) *******************
	// vector<Matrix_imcomp> cube_para_dev_cis_gene;
	cube_para_dev_cis_gene[etissue_index].clean();

	// Matrix matrix_para_dev_snp_cellenv;
	matrix_para_dev_snp_cellenv.clean();

	// vector<Matrix> matrix_para_dev_cellenv_gene;
	cube_para_dev_cellenv_gene[etissue_index].clean();

	// Matrix matrix_para_dev_batch_batch_hidden;
	matrix_para_dev_batch_batch_hidden.clean();

	// Matrix matrix_para_dev_batch_hidden_gene;
	matrix_para_dev_batch_hidden_gene.clean();

}




// add the results from one thread into the final repo
void aggregation_add(package_dev * para_array_pointer, int etissue_index)
{
	long int dimension1;
	long int dimension2;


	//********************************* aggregation of this mini-batch *****************************************
	// vector<Matrix_imcomp> cube_para_dev_cis_gene;
	dimension1 = cube_para_dev_cis_gene[etissue_index].get_dimension1();
	for(long int i=0; i<dimension1; i++)
	{
		dimension2 = cube_para_dev_cis_gene[etissue_index].get_dimension2(i);
		for(long int j=0; j<dimension2; j++)
		{
			float dev = (para_array_pointer->matrix_imcomp_para_dev_cis_gene).get(i, j);
			cube_para_dev_cis_gene[etissue_index].add_on(i, j, dev);
		}
	}

	// Matrix matrix_para_dev_snp_cellenv;
	dimension1 = matrix_para_dev_snp_cellenv.get_dimension1();
	for(long int i=0; i<dimension1; i++)
	{
		dimension2 = matrix_para_dev_snp_cellenv.get_dimension2();
		for(long int j=0; j<dimension2; j++)
		{
			float dev = (para_array_pointer->matrix_para_dev_snp_cellenv).get(i, j);
			matrix_para_dev_snp_cellenv.add_on(i, j, dev);

		}
	}

	// vector<Matrix> cube_para_dev_cellenv_gene;
	dimension1 = cube_para_dev_cellenv_gene[etissue_index].get_dimension1();
	for(long int i=0; i<dimension1; i++)
	{
		dimension2 = cube_para_dev_cellenv_gene[etissue_index].get_dimension2();
		for(long int j=0; j<dimension2; j++)
		{
			float dev = (para_array_pointer->matrix_para_dev_cellenv_gene).get(i, j);
			cube_para_dev_cellenv_gene[etissue_index].add_on(i, j, dev);
		}
	}

	// Matrix matrix_para_dev_batch_batch_hidden;
	dimension1 = matrix_para_dev_batch_batch_hidden.get_dimension1();
	for(long int i=0; i<dimension1; i++)
	{
		dimension2 = matrix_para_dev_batch_batch_hidden.get_dimension2();
		for(long int j=0; j<dimension2; j++)
		{
			float dev = (para_array_pointer->matrix_para_dev_batch_batch_hidden).get(i, j);
			matrix_para_dev_batch_batch_hidden.add_on(i, j, dev);
		}
	}

	// Matrix matrix_para_dev_batch_hidden_gene;
	dimension1 = matrix_para_dev_batch_hidden_gene.get_dimension1();
	for(long int i=0; i<dimension1; i++)
	{
		dimension2 = matrix_para_dev_batch_hidden_gene.get_dimension2();
		for(long int j=0; j<dimension2; j++)
		{
			float dev = (para_array_pointer->matrix_para_dev_batch_hidden_gene).get(i, j);
			matrix_para_dev_batch_hidden_gene.add_on(i, j, dev);
		}
	}

}




// average the summed results from all threads
void aggregation_ave(int etissue_index)
{
	//********************************* averaging of this mini-batch *****************************************
	cout << "aggregation of this mini-batch..." << endl;

	// vector<Matrix_imcomp> cube_para_dev_cis_gene;
	cube_para_dev_cis_gene[etissue_index].scale( 1.0 / batch_size );

	// Matrix matrix_para_dev_snp_cellenv;
	matrix_para_dev_snp_cellenv.scale( 1.0 / batch_size );

	// vector<Matrix> cube_para_dev_cellenv_gene;
	cube_para_dev_cellenv_gene[etissue_index].scale( 1.0 / batch_size );

	// Matrix matrix_para_dev_batch_batch_hidden;
	matrix_para_dev_batch_batch_hidden.scale( 1.0 / batch_size );

	// Matrix matrix_para_dev_batch_hidden_gene;
	matrix_para_dev_batch_hidden_gene.scale( 1.0 / batch_size );

}




// this is the working program for each thread, for calculating the log-likelihood
void * WorkPerThread_loglike_testerror(void * pointer)
{
	float result = 0;

	package_loglike_testerror * package_pointer = (package_loglike_testerror *)pointer;
	int indicator = package_pointer->indicator;
	string etissue = package_pointer->etissue;
	//=================================================================================================================
	//******************************************* loglike or testerror ************************************************
	//=================================================================================================================
	int amount_sample;
	if(indicator == 0)
	{
		amount_sample = esample_tissue_rep[etissue].size();
	}
	else
	{
		amount_sample = esample_tissue_rep_test[etissue].size();
	}


	//======== allocating memory for local containers ========
	array<float *, NUM_CHR> snp_dosage_list;
	float * gene_rpkm_exp;  // with length "num_gene"
	float * cellenv_hidden_var;  // with length "num_cellenv"
	float * batch_var;  // with length "num_batch"
	float * batch_hidden_var;  // with length "num_batch_hidden"
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


	while(1)
	{
		int count = -1;
		//================ check whether there are still left samples to be processed ================
		// if there are, take into this thread; otherwise, terminate this thread
		pthread_mutex_lock(&mut);
		for(int i=0; i<amount_sample; i++)
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

		// NOTE: debugging
		//cout << package_pointer->id << ": " << count << endl;

		//================ work on the current sample ================
		//=================================================================================================================
		//******************************************* loglike or testerror ************************************************
		//=================================================================================================================
		if(indicator == 0)
		{
			result += forward_loglike_testerror(0, etissue, count, &snp_dosage_list, gene_rpkm_exp, cellenv_hidden_var, batch_var, batch_hidden_var);	// indicator=0 --> loglike; indicator=1 --> testerror
		}
		else
		{
			result += forward_loglike_testerror(1, etissue, count, &snp_dosage_list, gene_rpkm_exp, cellenv_hidden_var, batch_var, batch_hidden_var);	// indicator=0 --> loglike; indicator=1 --> testerror
		}
	}

	package_pointer->result = result;


	// free memory
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




	pthread_exit(NULL);
}




// calculate the log-likelihood in a multithreading fashion
float cal_loglike_multithread(string etissue)
{
	cout << "[@@@] now calculating the log-likelihood (multithreading)..." << endl;

	float loglike = 0;
	int etissue_index = etissue_index_map[etissue];
	int amount_sample = esample_tissue_rep[etissue].size();


	//=============================== multi-threading parameter initialization ===============================
	// allocating all the other threads from here
	pthread_t threads[num_thread];
	void * status;
	pthread_attr_t attr;
	// Initialize and set thread joinable
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	//
	finish_table = (int *)calloc(amount_sample, sizeof(int));		// calloc --> initialize as 0
	//
	pthread_mutex_init(&mut, NULL);
	memset(&threads, 0, sizeof(threads));

	//=============================== thread local memory preparation (not yet allocation heap) ===============================
	package_loglike_testerror para_array[num_thread];
	for(int i=0; i<num_thread; i++)
	{
		package_loglike_testerror package;
		package.id = i;
		package.indicator = 0;		// indicating loglike
		package.etissue = etissue;
		package.result = 0;		// to be filled with the thread work
		para_array[i] = package;
	}

	//=============================== thread initialization ===============================
	for(int i=0; i<num_thread; i++)
	{
		cout << "main() : creating thread#" << i+1 << endl;
		int rc = pthread_create(&threads[i], NULL, WorkPerThread_loglike_testerror, (void *)&para_array[i]);
		if(rc)
		{
			cout << "Error:unable to create thread," << rc << endl;
			exit(-1);
		}
	}


	//===================== waiting for all the threads to terminate, and aggregate =====================
	loglike = 0;
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
		cout << " exiting with status: " << status;
		cout << ", aggregate it's results..." << endl;

		// aggregate the results from this current thread
		loglike += para_array[i].result;
	}


	//===================== finish and quit =====================
	free(finish_table);
	cout << "[@@@] finishing the current loglike cal (in multi-treading mode)..." << endl;


	return loglike;
}




float cal_testerror_multithread(string etissue)
{
	cout << "[@@@] now calculating the testing error (multithreading)..." << endl;

	float testerror = 0;
	int etissue_index = etissue_index_map[etissue];
	int amount_sample = esample_tissue_rep_test[etissue].size();


	//=============================== multi-threading parameter initialization ===============================
	// allocating all the other threads from here
	pthread_t threads[num_thread];
	void * status;
	pthread_attr_t attr;
	// Initialize and set thread joinable
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	//
	finish_table = (int *)calloc(amount_sample, sizeof(int));		// calloc --> initialize as 0
	//
	pthread_mutex_init(&mut, NULL);
	memset(&threads, 0, sizeof(threads));

	//=============================== thread local memory preparation (not yet allocation heap) ===============================
	package_loglike_testerror para_array[num_thread];
	for(int i=0; i<num_thread; i++)
	{
		package_loglike_testerror package;
		package.id = i;
		package.indicator = 1;		// indicating testerror
		package.etissue = etissue;
		package.result = 0;		// to be filled with the thread work
		para_array[i] = package;
	}

	//=============================== thread initialization ===============================
	for(int i=0; i<num_thread; i++)
	{
		cout << "main() : creating thread#" << i+1 << endl;
		int rc = pthread_create(&threads[i], NULL, WorkPerThread_loglike_testerror, (void *)&para_array[i]);
		if(rc)
		{
			cout << "Error:unable to create thread," << rc << endl;
			exit(-1);
		}
	}


	//===================== waiting for all the threads to terminate, and aggregate =====================
	testerror = 0;
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
		cout << " exiting with status: " << status;
		cout << ", aggregate it's results..." << endl;

		// aggregate the results from this current thread
		testerror += para_array[i].result;
	}


	//===================== finish and quit =====================
	free(finish_table);
	cout << "[@@@] finishing the current loglike cal (in multi-treading mode)..." << endl;


	return testerror;
}


