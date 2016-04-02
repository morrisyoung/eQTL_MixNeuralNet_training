// this contains the mt versions of some functions of operations on the matrix type
//	1. regularization
//	2. gradient descent

#include "lib_matrix.h"
#include <math.h>
#include <string>
#include <vector>
#include <unordered_map>
#include "global.h"
#include <array>
#include "opt_nn_acfunc.h"
#include <math.h>       /* exp */
#include "libfunc_matrix_mt.h"



using namespace std;



// local global variables
pthread_mutex_t mut_regu;				// mutex used by all the threads
pthread_mutex_t mut_gd;				// mutex used by all the threads
int * finish_table_regu;			// finish table for all the samples in this batch ()
int * finish_table_gd;				// finish table for all the samples in this batch ()




// Apr.1: multi-threading versions of the following functions:
//	void para_penalty_lasso_approx(Matrix matrix_para, Matrix matrix_para_dev, float lambda, float sigma)
//	void para_penalty_cis(Matrix_imcomp matrix_imcomp_para, Matrix_imcomp matrix_imcomp_para_dev, vector<vector<float>> repo_prior, float lambda_lasso, float lambda_ridge, float sigma)
//	void para_gradient_descent(Matrix matrix_para, Matrix matrix_para_dev, float rate)
//	void para_gradient_descent_cis(Matrix_imcomp matrix_imcomp_para, Matrix_imcomp matrix_imcomp_para_dev, float rate)
// I can do this because this happens after the whole mini-batch, so there are no multi-threading any more
// I will write the multi-threading version of these routines in "libfunc_matrix_mt.cpp"





//=============================================================================================================
//============================================ regularization =================================================
//=============================================================================================================
//==== local global variables ====
Matrix matrix_para_regu_global;
Matrix matrix_para_dev_regu_global;
float lambda_global;
float sigma_global;


// this is the working program for each thread, for gradient descent
void * WorkPerThread_penalty_lasso_approx(void * pointer)
{
	int * package_pointer = (int *)pointer;
	int thread_index = (* package_pointer);
	long int dimension1 = matrix_para_regu_global.get_dimension1();
	long int dimension2 = matrix_para_regu_global.get_dimension2();


	while(1)
	{
		int count = -1;
		//================ check whether there are still left samples to be processed ================
		// if there are, take into this thread; otherwise, terminate this thread
		pthread_mutex_lock(&mut_regu);
		for(long int i=0; i<dimension1; i++)
		{
			if(finish_table_regu[i] == 0)
			{
				count = i;
				finish_table_regu[i] = 1;
				break;
			}
		}
		pthread_mutex_unlock(&mut_regu);

		if(count == -1)  // no left samples for processing
		{
			break;
		}

		//================ work on the current sample ================
		for(long int j=0; j<dimension2; j++)
		{
			/// the value of current para beta:
			float beta = matrix_para_regu_global.get(count, j);
			/// the derivative of the beta:
			float derivative = beta / sqrt (beta * beta + sigma_global);  // this is an approximation of the LASSO regularization
			/// and the value of its derivative should be added with that derivative term for regularization:
			matrix_para_dev_regu_global.add_on(count, j, lambda_global * derivative);
		}

	}// all samples finished


	pthread_exit(NULL);
}


// func: add the LASSO penalty term in the matrix_para_dev; destroy type
//		the matrix_para is used to calculate the derivative term, that will be added to the matrix_para_dev
void para_penalty_lasso_approx_mt(Matrix matrix_para, Matrix matrix_para_dev, float lambda, float sigma)
{
	cout << "[@@@] now do the regularization for Matrix (multithreading)..." << endl;
	//==== local global variables assignment
	matrix_para_regu_global = matrix_para;
	matrix_para_dev_regu_global = matrix_para_dev;
	lambda_global = lambda;
	sigma_global = sigma;

	long int dimension1 = matrix_para.get_dimension1();

	//=============================== multi-threading parameter initialization ===============================
	// allocating all the other threads from here
	pthread_t threads[num_thread];
	void * status;
	pthread_attr_t attr;
	// Initialize and set thread joinable
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	//
	finish_table_regu = (int *)calloc(dimension1, sizeof(int));		// calloc --> initialize as 0
	//
	pthread_mutex_init(&mut_regu, NULL);
	memset(&threads, 0, sizeof(threads));

	//=============================== thread local memory preparation (not yet allocation heap) ===============================
	int para_array[num_thread];			// we only need the thread# in this program
	for(int i=0; i<num_thread; i++)
	{
		para_array[i] = i;
	}

	//=============================== thread initialization ===============================
	for(int i=0; i<num_thread; i++)
	{
		cout << "main() : creating thread#" << i+1 << endl;
		int rc = pthread_create(&threads[i], NULL, WorkPerThread_penalty_lasso_approx, (void *)&para_array[i]);
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

	//===================== finish and quit =====================
	free(finish_table_regu);
	cout << "[@@@] finishing the gradient descent for Matrix (multi-treading mode)..." << endl;

	return;
}























/*
// func:
//		regularization for the cis- parameters, with both lasso and ridge (composed as elastic net)
//		for sparsity and pruning signal
void para_penalty_cis_mt(Matrix_imcomp matrix_imcomp_para, Matrix_imcomp matrix_imcomp_para_dev, vector<vector<float>> repo_prior, float lambda_lasso, float lambda_ridge, float sigma)
{

	long int dimension1 = matrix_imcomp_para.get_dimension1();
	for(long int i=0; i<dimension1; i++)
	{
		int chr = matrix_imcomp_para.get_chr(i);
		long int sst = matrix_imcomp_para.get_sst(i);

		long int dimension2 = matrix_imcomp_para.get_dimension2(i);
		for(long int j=0; j<dimension2; j++)
		{
			long int pos = j + sst;  // the pos in snp list

			/// the value of current cis- beta:
			float beta = matrix_imcomp_para.get(i, j);

			/// the prior that we need (if there is) for tuning the relative strength of L1 and L2 regularization:
			// TODO: we don't load this currently; NOTE that we don't have the prior information for the intercept term (the last term)
			//prior = prior_tissue_rep[etissue][chr-1][pos];
			float prior = 1.0;

			float alpha = 1 / ( 1 + exp(-(prior-1)) );

			/// the derivative of the beta:
			float derivative1 = beta / sqrt (beta * beta + sigma);  // this is an approximation of the LASSO regularization
			float derivative2 = 2 * beta;  // L2 regularization item is differentiable

			/// and the value of its derivative should be added with that derivative item from regularization:
			float value = lambda_lasso * (1 - alpha) * derivative1 + lambda_ridge * alpha * derivative2;
			matrix_imcomp_para_dev.add_on(i, j, value);
		}

	}

	return;
}
*/










//===============================================================================================================
//============================================ gradient descent =================================================
//===============================================================================================================
//==== local global variables ====
Matrix matrix_para_global;
Matrix matrix_para_dev_global;
float rate_global;


// this is the working program for each thread, for gradient descent
void * WorkPerThread_gradient_descent(void * pointer)
{
	int * package_pointer = (int *)pointer;
	int thread_index = (* package_pointer);
	long int dimension1 = matrix_para_global.get_dimension1();
	long int dimension2 = matrix_para_global.get_dimension2();

	while(1)
	{
		int count = -1;
		//================ check whether there are still left samples to be processed ================
		// if there are, take into this thread; otherwise, terminate this thread
		pthread_mutex_lock(&mut_gd);
		for(long int i=0; i<dimension1; i++)
		{
			if(finish_table_gd[i] == 0)
			{
				count = i;
				finish_table_gd[i] = 1;
				break;
			}
		}
		pthread_mutex_unlock(&mut_gd);

		if(count == -1)  // no left samples for processing
		{
			break;
		}

		//================ work on the current sample ================
		for(long int j=0; j<dimension2; j++)
		{
			float dev = matrix_para_dev_global.get(count, j);
			matrix_para_global.add_on(count, j, - rate_global * dev);
		}

	}// all samples finished


	pthread_exit(NULL);
}


// func: do the gradient descent for this parameter matrix
//		with: beta = beta - rate * dev
void para_gradient_descent_mt(Matrix matrix_para, Matrix matrix_para_dev, float rate)
{
	cout << "[@@@] now do the gradient descent for Matrix (multithreading)..." << endl;
	//==== local global variables assignment
	matrix_para_global = matrix_para;
	matrix_para_dev_global = matrix_para_dev;
	rate_global = rate;

	long int dimension1 = matrix_para.get_dimension1();

	//=============================== multi-threading parameter initialization ===============================
	// allocating all the other threads from here
	pthread_t threads[num_thread];
	void * status;
	pthread_attr_t attr;
	// Initialize and set thread joinable
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	//
	finish_table_gd = (int *)calloc(dimension1, sizeof(int));		// calloc --> initialize as 0
	//
	pthread_mutex_init(&mut_gd, NULL);
	memset(&threads, 0, sizeof(threads));

	//=============================== thread local memory preparation (not yet allocation heap) ===============================
	int para_array[num_thread];			// we only need the thread# in this program
	for(int i=0; i<num_thread; i++)
	{
		para_array[i] = i;
	}

	//=============================== thread initialization ===============================
	for(int i=0; i<num_thread; i++)
	{
		cout << "main() : creating thread#" << i+1 << endl;
		int rc = pthread_create(&threads[i], NULL, WorkPerThread_gradient_descent, (void *)&para_array[i]);
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

	//===================== finish and quit =====================
	free(finish_table_gd);
	cout << "[@@@] finishing the gradient descent for Matrix (multi-treading mode)..." << endl;

	return;
}





//==== local global variables ====
Matrix_imcomp matrix_imcomp_para_global;
Matrix_imcomp matrix_imcomp_para_dev_global;
//float rate_global;		// have this before


// this is the working program for each thread, for gradient descent
void * WorkPerThread_gradient_descent_cis(void * pointer)
{
	int * package_pointer = (int *)pointer;
	int thread_index = (* package_pointer);
	long int dimension1 = matrix_imcomp_para_global.get_dimension1();

	while(1)
	{
		int count = -1;
		//================ check whether there are still left samples to be processed ================
		// if there are, take into this thread; otherwise, terminate this thread
		pthread_mutex_lock(&mut_gd);
		for(long int i=0; i<dimension1; i++)
		{
			if(finish_table_gd[i] == 0)
			{
				count = i;
				finish_table_gd[i] = 1;
				break;
			}
		}
		pthread_mutex_unlock(&mut_gd);

		if(count == -1)  // no left samples for processing
		{
			break;
		}

		//================ work on the current sample ================
		long int dimension2 = matrix_imcomp_para_global.get_dimension2(count);
		for(long int j=0; j<dimension2; j++)
		{
			float dev = matrix_imcomp_para_dev_global.get(count, j);
			matrix_imcomp_para_global.add_on(count, j, - rate_global * dev);
		}

	}// all samples finished


	pthread_exit(NULL);
}


// func: do the gradient descent for this parameter matrix_incomp
//		with: beta = beta - rate * dev
void para_gradient_descent_cis_mt(Matrix_imcomp matrix_imcomp_para, Matrix_imcomp matrix_imcomp_para_dev, float rate)
{
	cout << "[@@@] now do the gradient descent for Matrix_imcomp (multithreading)..." << endl;
	//==== local global variables assignment
	matrix_imcomp_para_global = matrix_imcomp_para;
	matrix_imcomp_para_dev_global = matrix_imcomp_para_dev;
	rate_global = rate;

	long int dimension1 = matrix_imcomp_para.get_dimension1();

	//=============================== multi-threading parameter initialization ===============================
	// allocating all the other threads from here
	pthread_t threads[num_thread];
	void * status;
	pthread_attr_t attr;
	// Initialize and set thread joinable
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	//
	finish_table_gd = (int *)calloc(dimension1, sizeof(int));		// calloc --> initialize as 0
	//
	pthread_mutex_init(&mut_gd, NULL);
	memset(&threads, 0, sizeof(threads));

	//=============================== thread local memory preparation (not yet allocation heap) ===============================
	int para_array[num_thread];			// we only need the thread# in this program
	for(int i=0; i<num_thread; i++)
	{
		para_array[i] = i;
	}

	//=============================== thread initialization ===============================
	for(int i=0; i<num_thread; i++)
	{
		cout << "main() : creating thread#" << i+1 << endl;
		int rc = pthread_create(&threads[i], NULL, WorkPerThread_gradient_descent_cis, (void *)&para_array[i]);
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

	//===================== finish and quit =====================
	free(finish_table_gd);
	cout << "[@@@] finishing the gradient descent for Matrix_imcomp (multi-treading mode)..." << endl;

	return;
}



