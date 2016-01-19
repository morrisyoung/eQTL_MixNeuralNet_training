// this contains the operations on the matrix type



#include "libfunc_matrix.h"
#include "lib_matrix.h"
#include <math.h>
#include <string>
#include <vector>
#include <unordered_map>
#include "global.h"
#include <array>





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
		for(long int j=0; j<dimension2-1; j++)
		{
			float para = matrix_para.get(i, j);
			result[i] += input[j] * para;
		}
		result[i] += matrix_para.get(i, dimension2 - 1);		// here we do have the regression intercept term
	}

	return;
}



// func: multiply the SNP array with the cis- SNP parameter matrix, to get the expression list
// TODO: this might need to be changed, as I don't want to bring the global variables in this routine
void multi_array_matrix_imcomp(array<float *, NUM_CHR> * input_pointer, Matrix_imcomp matrix_imcomp_para, float * result)
{
	long int dimension1 = matrix_imcomp_para.get_dimension1();
	for(long int i=0; i<dimension1; i++)
	{
		int chr = matrix_imcomp_para.get_chr(i);

		long int dimension2 = matrix_imcomp_para.get_dimension2(i);
		for(long int j=0; j<dimension2; j++)
		{
			if(j == dimension2 - 1)
			{
				result[i] += 1 * par;			// the last one in the parameter list is for the intercept term
			}
			else
			{
				int pos = matrix_imcomp_para.get_sst(i) + j;
				float var = (* input_pointer)[chr][pos];
				float par = matrix_imcomp_para.get(i, j);
				result[i] += var * par;
			}
		}
	}

	return;
}



// func: multiply the SNP array list with the parameter matrix
void multi_array_list_matrix(array<float *, NUM_CHR> * input_pointer, Matrix matrix_para, float * result)
{
	long int dimension1 = matrix_para.get_dimension1();
	long int dimension2 = matrix_para.get_dimension2();
	for(int i=0; i<dimension1; i++)
	{
		result[i] = 0;
		long int count = 0;
		for(int j=0; j<NUM_CHR; j++)  // across all the chromosomes
		{
			for(long k=0; k<snp_name_list[j].size(); k++)			// TODO: this is to be corrected, as we don't want to see global variables here
			{
				float var = (*dosage_list_pointer)[j][k];
				float par = matrix_para[i][count];
				cellenv_con_pointer[i] += var * par;
				count ++;
			}
		}
		float par = matrix_para[i][dimension2 - 1];					// we do have the intercept term here
		cellenv_con_pointer[i] += par;
	}

	return;
}


