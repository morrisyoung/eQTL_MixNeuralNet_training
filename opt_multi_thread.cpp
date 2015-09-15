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



using namespace std;



// TODO: fill them into the .h file
int num_thread = 8;				// there are at maximum 8 cores in C2B2 cluster, but our main thread doesn't do extensive computation here
pthread_mutex_t mut;			// mutex used by all the threads
int * finish_table;				// finish table for all the samples in this batch ()


// each thread should have such a local parameter space
typedef struct package_dev
{
	array<float *, 22> snp_dosage_list;
	float * gene_rpkm_exp;  	// with length "num_gene"
	float * cellenv_hidden_var; // with length "num_cellenv"
	float * batch_var;			// with length "num_batch"
	float * batch_hidden_var;	// with length "num_batch_hidden"

	vector<vector<float *>> para_dev_cis_gene;
	vector<float *> para_dev_snp_cellenv;
	vector<vector<float *>> para_dev_cellenv_gene;
	vector<float *> para_dev_batch_batch_hidden;
	vector<float *> para_dev_batch_hidden_gene;

}package_dev;


// to fill in the above space
void package_alloc(package_dev * package_pointer)
{
	//

}

// to release the above space
void package_free(package_dev * package_pointer)
{
	//

}




void * WorkPerThread(void * threadid)
{
   long tid;
   tid = (long)threadid;
   cout << "Thread ID " << tid << " is working..." << endl;

   int mark = 1;

   while(mark)
   {
   	//
   	
   }


   //================ check whether there are still left samples to be processed ================
   // if there are, take into this thread; otherwise, terminate this thread
	pthread_mutex_lock(&mut);
	// add the present segment to the reporting list
	if(report_start == NULL && entrance == NULL)
	{
		report_start = ibd_seg;
		entrance = ibd_seg;
	}
	else
	{
		entrance->next = ibd_seg;
		entrance = ibd_seg;
	}
	pthread_mutex_unlock(&mut);



	pthread_exit(NULL);

}


void opt_mt_control(string etissue, int pos_start, int num_esample)
{
	//=============================== multi-threading parameter initialization ===============================
	// allocating all the other threads from here
	pthread_t threads[num_thread];
	void *status;
	pthread_attr_t attr;
	// Initialize and set thread joinable
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	finish_table = (int *)calloc(batch_size, sizeof(int));

	pthread_mutex_init(&mut, NULL);
	memset(&thread, 0, sizeof(thread));


	//=============================== thread local memory allocation ===============================
	package_dev para_array[num_thread];
	for(int i=0; i<num_thread; i++)
	{
		package_dev package;
		package_alloc(&package);
		// append
	}



	//=============================== thread initialization ===============================
	for(int i=0; i<num_thread; i++)
	{
		cout << "main() : creating thread, " << i << endl;
		int rc = pthread_create(&threads[i], NULL, WorkPerThread, (void *)&para_array[i]);
		if(rc)
		{
			cout << "Error:unable to create thread," << rc << endl;
			exit(-1);
		}
	}



	//===================== waiting for all the threads to terminate =====================
	for(int i=0; i<num_thread; i++)
	{
		int rc = pthread_join(threads[i], &status);
		if (rc)
		{
			cout << "Error:unable to join," << rc << endl;
			exit(-1);
		}
		cout << "Main: completed thread id :" << i ;
		cout << "  exiting with status :" << status << endl;
	}


	//===================== merge results, and release space =====================
	//// fill in the true para_dev_xxx space (aggregation)
	aggregation(&para_array);
	//// add the regularization terms into the derivatives
	regualrization(etissue);
	/// gradient descent
	gradient_descent(etissue);

	//// free the memory allocated for these threads
	for(int i=0; i<num_thread; i++)
	{
		package_free(&para_array[i]);
	}


	//===================== finish and quit =====================
	free(finish_table);
	cout << "[@@@] finishing the current mini-batch..." << endl;
}


void aggregation(package_dev * para_array_pointer)
{
	//
	//

}

