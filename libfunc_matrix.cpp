// this contains the operations on the matrix type



#include "libfunc_matrix.h"
#include "lib_matrix.h"
#include <math.h>
#include <string>
#include <vector>
#include <unordered_map>





using namespace std;




//============================================ task specific functions =================================================
// func: add the LASSO penalty term in the matrix_para_dev; destroy type
//		the matrix_para is used to calculate the derivative term, that will be added to the matrix_para_dev
void para_penalty_lasso_approx(Matrix matrix_para, Matrix matrix_para_dev, float lambda, float sigma)
{
	long int dimension1 = matrix_para.get_dimension1();
	long int dimension2 = matrix_para.get_dimension2();



	for(int i=0; i<dimension1; i++)
	{
		for(int j=0; j<dimension2; j++)
		{
			/// the value of current para beta:
			float beta = matrix_para.get(i, j);
			/// the derivative of the beta:
			float derivative = beta / sqrt (beta * beta + sigma);  // this is an approximation of the LASSO regularization
			/// and the value of its derivative should be added with that derivative term for regularization:
			matrix_para_dev.add_on(i, j, lambda * derivative);
		}
	}

	return;
}



// func:
//		regularization for the cis- parameters, with both lasso and ridge (composed as elastic net)
//		for sparsity and pruning signal
void para_penalty_cis(Matrix_imcomp matrix_imcomp_para, Matrix_imcomp matrix_imcomp_para_dev, vector<vector<float>>> repo_prior, float lambda_lasso, float lambda_ridge, float sigma)
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
			//prior = prior_tissue_rep[etissue][chr][pos];
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



// func: do the gradient descent for this parameter matrix
//		with: beta = beta - rate * dev
void para_gradient_descent(Matrix matrix_para, Matrix matrix_para_dev, float rate)
{
	long int dimension1 = matrix_para.get_dimension1();
	long int dimension2 = matrix_para.get_dimension2();

	for(long int i=0; i<dimension1; i++)
	{
		for(long int j=0; j<dimension2; j++)
		{
			float dev = matrix_para_dev.get(i, j);
			matrix_para.add_on(i, j, - rate * dev);
		}
	}

	return;
}


// func: do the gradient descent for this parameter matrix_incomp
//		with: beta = beta - rate * dev
void para_gradient_descent_cis(Matrix_imcomp matrix_imcomp_para, Matrix_imcomp matrix_imcomp_para_dev, float rate)
{
	long int dimension1 = matrix_imcomp_para.get_dimension1();
	for(long int i=0; i<dimension1; i++)
	{
		long int dimension2 = matrix_imcomp_para.get_dimension2(i);
		for(long int j=0; j<dimension2; j++)
		{
			float dev = matrix_imcomp_para_dev.get(i, j);
			matrix_imcomp_para.add_on(i, j, - rate * dev);
		}
	}

	return;
}



//============================================ abstract functions =================================================
// func: multiplay the var array with the parameter matrix to get the result
// input: variable array; parameter matrix
// output: result variable array
void multi_array_matrix(float * input, Matrix matrix_para, float * result)
{
	long int dimension1 = matrix_para.get_dimension1();
	long int dimension2 = matrix_para.get_dimension2();

	for(long int i=0; i< dimension1; i++)
	{
		result[i] = 0;
		for(long int j=0; j<dimension2; j++)
		{
			float para = matrix_para.get(i, j);
			result[i] += input[j] * para;
		}
	}

	return;
}


