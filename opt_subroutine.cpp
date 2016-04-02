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
#include "opt_nn_acfunc.h"
#include "opt_debugger.h"
#include "libfunc_matrix.h"
#include "libfunc_matrix_mt.h"
#include <cmath>        // std::abs




using namespace std;




// forward and backward propagation for one mini-batch
void forward_backward_prop_batch(string etissue, int pos_start, int num_esample)
{
	cout << "[@@] entering the forward-backward propagation..." << endl;

	int etissue_index = etissue_index_map[etissue];

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



	//****************************** enter the mini-batch ***********************************
	cout << "we are entering a new mini-batch..." << endl;
	for(int count=0; count<batch_size; count++)
	{
		int pos = (pos_start + count) % (num_esample);
		string esample = esample_tissue_rep[etissue][pos];
		string individual = sample_to_individual(esample);
		cout << "current sample #" << pos+1 << ": " << esample << endl;

		//=================================================== init ============================================================
		// get the: 0. esample and individual; 1. genotype; 2. expression data; 3. batch variables
		// to: 1. forward_backward propagation;
		// genotype dosage data
		//cout << "getting the dosage data for individual #" << individual << endl;
		snp_dosage_load(&snp_dosage_list, individual);  // snp dosage data for one individual across all chromosomes
		// expression rpkm data: eQTL_tissue_rep[etissue][esample]
		//cout << "we have this amount of genes expressed in this individual:" << eQTL_tissue_rep[etissue][esample].size() << endl;
		// and the batch variable for this individual and this sample
		int num_batch_individual = batch_individual[individual].size();
		int index = 0;
		for(int i=0; i<num_batch_individual; i++)
		{
			float value = batch_individual[individual][i];
			batch_var[index] = value;
			index++;
		}
		int num_batch_sample = batch_sample[esample].size();
		for(int i=0; i<num_batch_sample; i++)
		{
			float value = batch_sample[esample][i];
			batch_var[index] = value;
			index++;
		}

		forward_backward(etissue_index,
						&snp_dosage_list,
						&eQTL_tissue_rep[etissue][esample],

						gene_rpkm_exp,
						cellenv_hidden_var,
						batch_var,
						batch_hidden_var,

						cube_para_dev_cis_gene[etissue_index],
						matrix_para_dev_snp_cellenv,
						cube_para_dev_cellenv_gene[etissue_index],
						matrix_para_dev_batch_batch_hidden,
						matrix_para_dev_batch_hidden_gene
						);



		// DEBUG
		// for this sample (or for this individual as we fix the tissue type currenlty):
		// I want to check the genotype of this individual, as the problem seems to come up there!!!
		// save: individual ID; genotype
		// debug here and also in the subroutine
		// for(long int i=0; i<snp_name_list[0].size(); i++)
		// {
		// 	float dosage = snp_dosage_list[0][i];
		// 	if(isnan(dosage))
		// 	{
		// 		cout << individual << " " << i << " " << dosage << endl;
		// 	}
		// }



		// iterating in this mini-batch
	}



	// DEBUG: debug the parameter dev, and the current cellenv and hiddenbatch
	para_temp_save_dev(etissue_index);

	// DEBUG: save some variables (is this so necessary especially after a whole batch?)
	char filename[100] = "../result_tempdata/var_cellenv.txt";
	para_temp_save_var(cellenv_hidden_var, num_cellenv, filename);
	sprintf(filename, "%s", "../result_tempdata/var_batch_hidden.txt");
	para_temp_save_var(batch_hidden_var, num_batch_hidden, filename);



	//********************************* aggregation of this mini-batch *****************************************
	// 1. average the derivatives calculated from previous steps
	// 2. will add the derivatives due to regularization in the next part
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




	//===================================== Regularization in Regression =====================================
	regularization(etissue_index);



	//=========================================== Gradient Descent ===========================================
	gradient_descent(etissue_index);




	cout << "[@@] leaving the forward-backward propagation..." << endl;

}




// NOTE: ADDITIVE
// property of this function: additive to the total derivative after this round (additive)
// what we need for the following routine:
// dosage list; expression value list; expression list; cellenv list; batch list; batch hidden list; ALL parameter (derivative) containers
void forward_backward(int etissue_index,
	array<float *, NUM_CHR> * dosage_list_pointer,
	vector<float> * expr_list_pointer,

	float * expr_con_pointer,
	float * cellenv_con_pointer,
	float * batch_list_pointer,
	float * batch_hidden_con_pointer,
	// the new Matrix/Matrix_imcomp classes:
	Matrix_imcomp matrix_imcomp_para_dev_cis_gene,				// drop the Matrix_imcomp object, other than the full cube. TODO: not sure whether there are problems
	Matrix matrix_para_dev_snp_cellenv,
	Matrix matrix_para_dev_cellenv_gene,						// drop the Matrix object, other than the full cube
	Matrix matrix_para_dev_batch_batch_hidden,
	Matrix matrix_para_dev_batch_hidden_gene
	)
{
	// how to map these pointers:
	/*
	dosage_list_pointer --> &snp_dosage_list
	expr_list_pointer --> &eQTL_tissue_rep[etissue][esample]
	batch_list_pointer --> batch_var
	expr_con_pointer  --> gene_rpkm_exp
	cellenv_con_pointer --> cellenv_hidden_var
	batch_hidden_con_pointer --> batch_hidden_var

	// all the other Matrix/Matrix_imcomp are directly called
	*/




	// DEBUG: check the Nan in genotype
	/*
	cout << "$$$" << endl;
	for(long int i=0; i<snp_name_list[0].size(); i++)
	{
		float dosage = snp_dosage_list[0][i];
		if(isnan(dosage))
		{
			cout << i << " " << dosage << endl;
		}
	}
	*/







	// DEBUG: this is used to debug the program
	char filename[100];



	/*
	// DEBUG: check the genotyoe list
	for(int j=0; j<NUM_CHR; j++)  // across all the chromosomes
	{
		for(long k=0; k<snp_name_list[j].size(); k++)			// TODO: this is to be corrected, as we don't want to see global variables here
		{
			float var = (*dosage_list_pointer)[j][k];
			cout << var << "\t";
		}
	}
	cout << endl;
	// DEBUG done: there is no problem
	*/




	//========================================================================
	// two step: forward propagation (get the function values); backward propagation (get the parameter derivatives)
	//========================================================================
	//========================================================================
	// step#1: forward-propogation (cis-; cell env; batch)
	//========================================================================
	//========================================================================

	// ****************************** [part1] cis- *********************************
	// for cis-, two issues:
	// 1. if this is a XYMT gene, we don't have signal from it's cis- SNPs (not consider currently);
	// 2. we use (gene_cis_index[gene].second - gene_cis_index[gene].first + 1) as the length of the cis- parameter array
	float * expr_con_pointer_cis = (float *)calloc( num_gene, sizeof(float) );
	multi_array_matrix_imcomp(dosage_list_pointer, cube_para_cis_gene[etissue_index], expr_con_pointer_cis);


	// // DEBUG mode: let's save the signal from all the three pathways
	// sprintf(filename, "%s", "../result_tempdata/var_expr_cis.txt");
	// para_temp_save_var(expr_con_pointer_cis, num_gene, filename);


	// ********************* [part2] cell env relevant parameters *********************
	// from snp to cell env variables
	float * expr_con_pointer_cellenv = (float *)calloc( num_gene, sizeof(float) );
	multi_array_list_matrix(dosage_list_pointer, matrix_para_snp_cellenv, cellenv_con_pointer);


	// // // DEBUG
	// sprintf(filename, "%s", "../result_tempdata/var_cellenv_before.txt");
	// para_temp_save_var(cellenv_con_pointer, num_cellenv, filename);


	//$$$$$$$$$$$ perform the activation function here (logistic or something else) $$$$$$$$$$$$
	neuralnet_ac_func(cellenv_con_pointer, num_cellenv);



	// // // DEBUG
	// sprintf(filename, "%s", "../result_tempdata/var_cellenv_after.txt");
	// para_temp_save_var(cellenv_con_pointer, num_cellenv, filename);



	// from cell env variables to genes
	multi_array_matrix(cellenv_con_pointer, cube_para_cellenv_gene[etissue_index], expr_con_pointer_cellenv);



	// // // DEBUG
	// sprintf(filename, "%s", "../result_tempdata/var_expr_cellenv.txt");
	// para_temp_save_var(expr_con_pointer_cellenv, num_gene, filename);




	// ********************* [part3] linear or non-linear batches *********************
	float * expr_con_pointer_batch = (float *)calloc( num_gene, sizeof(float) );
	// from original batch to hidden batch
	multi_array_matrix(batch_list_pointer, matrix_para_batch_batch_hidden, batch_hidden_con_pointer);


	// // DEBUG
	// sprintf(filename, "%s", "../result_tempdata/var_batch_hidden_before.txt");
	// para_temp_save_var(batch_hidden_con_pointer, num_batch_hidden, filename);



	//$$$$$$$$$$$ perform the activation function here (logistic or something else) $$$$$$$$$$$$
	neuralnet_ac_func(batch_hidden_con_pointer, num_batch_hidden);



	// // DEBUG
	// sprintf(filename, "%s", "../result_tempdata/var_batch_hidden_after.txt");
	// para_temp_save_var(batch_hidden_con_pointer, num_batch_hidden, filename);



	// from hidden batch to genes
	multi_array_matrix(batch_hidden_con_pointer, matrix_para_batch_hidden_gene, expr_con_pointer_batch);



	// // DEBUG
	// sprintf(filename, "%s", "../result_tempdata/var_expr_batch.txt");
	// para_temp_save_var(expr_con_pointer_batch, num_gene, filename);




	// ********************* [end] merge the signal from three pathways here, to expr_con_pointer *********************
	for(long int i=0; i<num_gene; i++)
	{
		expr_con_pointer[i] = expr_con_pointer_cis[i] + expr_con_pointer_cellenv[i] + expr_con_pointer_batch[i];
	}
	free(expr_con_pointer_cis);
	free(expr_con_pointer_cellenv);
	free(expr_con_pointer_batch);


	// error is the thing actually needed
	float * error_list = (float *)calloc(num_gene, sizeof(float));
	for(long int i=0; i<num_gene; i++)
	{
		error_list[i] = expr_con_pointer[i] - (*expr_list_pointer)[i];
	}



	// DEBUG
	sprintf(filename, "%s", "../result_tempdata/var_expr_exp.txt");
	para_temp_save_var(expr_con_pointer, num_gene, filename);
	// DEBUG: check the error
	sprintf(filename, "%s", "../result_tempdata/var_expr_error.txt");
	para_temp_save_var(error_list, num_gene, filename);




	//========================================================================
	//========================================================================
	// step#2: back-propogation (cis-;  cell env; batch)
	//========================================================================
	//========================================================================
	// *********************** [part1] cis- ************************
	backward_error_prop_direct_imcomp(matrix_imcomp_para_dev_cis_gene, error_list, dosage_list_pointer);


	// ***************** [part2] cell env relevant parameters *****************
	//// from cell env to genes
	backward_error_prop_last_layer(matrix_para_dev_cellenv_gene, error_list, cellenv_con_pointer);

	//// from snp to cell env
	backward_error_prop_inter_layer_1(error_list, cube_para_cellenv_gene[etissue_index], matrix_para_dev_snp_cellenv, cellenv_con_pointer, dosage_list_pointer);


	// ********************* [part3] linear or non-linear batches *********************
	//// from hidden batch to genes
	backward_error_prop_last_layer(matrix_para_dev_batch_hidden_gene, error_list, batch_hidden_con_pointer);

	// from original batch to hidden batch
	backward_error_prop_inter_layer_2(error_list, matrix_para_batch_hidden_gene, matrix_para_dev_batch_batch_hidden, batch_hidden_con_pointer, batch_list_pointer);



	free(error_list);

	return;
}







void regularization(int etissue_index)
{


	// Apr.2: I'm testing the speedup of mt code of this routine
	int marker_mt = 1;




	cout << "[@@] entering the regularization routine..." << endl;

	// there are several classes of prior knowledge that we need to consider
	// 1. sparsity of cis- regulation, accompanied by ridge regression, achieved by elastic-net tuned by the prior number, and the distance prior
	// 2. sparsity (LASSO) for the coefficients from cell env to expression (with the assumption that one gene is only affected by several handful cell env)
	// 3.1.[TODO] hierarchical regularization tuned by the learned tissue hierarchy
	// 3.2.[TODO] or we can simply use group LASSO to encourage the tissue consistency
	// 4. penalize the batch variables hashly (from batch variables to batch_hidden, and from batch_hidden to genes)
	cout << "adding the regularization items to the derivatives..." << endl;

	//===================================== part#0 =====================================
	// initialize some learning parameters

	// define the sigma here that may be used by L1 regularization
	float sigma = 0.0001;

	// regularization strength lambda:
	// path#1: add the lambda_{LASSO} and lambda_{ridge} for the cis- regularization
	// path#2: add the lambda for cellenv-gene regularization
	// path#3: and the lambda for batch-batch_hidden and batch_hidden-gene

	// TODO: DEBUG, see below
	//float lambda_lasso = 1.0;
	//float lambda_ridge = 1.0;
	//float lambda_snp_cellenv = 1.0;
	//float lambda_cellenv_gene = 1.0;
	//float lambda_batch_batch_hidden = 1.0;
	//float lambda_batch_hidden_gene = 1.0;


	//========================================================================================================
	// DEBUG
	// Feb.12: (let's test different lambda for one path -- batch path)
	// different lambda strength (listed below) x different batch iterations (200, 1000)
	// Feb.14: it seems the group#2 and group#3 are good

	// group#1
	//lambda = 0.1;
	// float lambda_lasso = 0.1;
	// float lambda_ridge = 0.1;
	// float lambda_snp_cellenv = 0.1;
	// float lambda_cellenv_gene = 0.1;
	// float lambda_batch_batch_hidden = 0.1;
	// float lambda_batch_hidden_gene = 0.1;

	// group#2
	//lambda = 0.001;
	float lambda_lasso = 0.001;
	float lambda_ridge = 0.001;
	float lambda_snp_cellenv = 0.001;
	float lambda_cellenv_gene = 0.001;
	float lambda_batch_batch_hidden = 0.001;
	float lambda_batch_hidden_gene = 0.001;

	// group#3
	// float lambda_lasso = 0.00001;
	// float lambda_ridge = 0.00001;
	// float lambda_snp_cellenv = 0.00001;
	// float lambda_cellenv_gene = 0.00001;
	// float lambda_batch_batch_hidden = 0.00001;
	// float lambda_batch_hidden_gene = 0.00001;

	// group#4
	// float lambda_lasso = 0.0000001;
	// float lambda_ridge = 0.0000001;
	// float lambda_snp_cellenv = 0.0000001;
	// float lambda_cellenv_gene = 0.0000001;
	// float lambda_batch_batch_hidden = 0.0000001;
	// float lambda_batch_hidden_gene = 0.0000001;




	// group#11
	// float lambda_lasso = 10;
	// float lambda_ridge = 10;
	// float lambda_snp_cellenv = 10;
	// float lambda_cellenv_gene = 10;
	// float lambda_batch_batch_hidden = 10;
	// float lambda_batch_hidden_gene = 10;

	// group#12
	// float lambda_lasso = 1000;
	// float lambda_ridge = 1000;
	// float lambda_snp_cellenv = 1000;
	// float lambda_cellenv_gene = 1000;
	// float lambda_batch_batch_hidden = 1000;
	// float lambda_batch_hidden_gene = 1000;

	// group#13
	// float lambda_lasso = 100000;
	// float lambda_ridge = 100000;
	// float lambda_snp_cellenv = 100000;
	// float lambda_cellenv_gene = 100000;
	// float lambda_batch_batch_hidden = 100000;
	// float lambda_batch_hidden_gene = 100000;


	//========================================================================================================




	if(marker_mt)
	{
		//===================================== part#1 =====================================
		// 1. sparsity of cis- regulation, accompanied by ridge regression, achieved by elastic-net tuned by the prior number, (and the distance prior)
		// TODO: not yet integrated the distance prior information
		para_penalty_cis_mt(cube_para_cis_gene[etissue_index], cube_para_dev_cis_gene[etissue_index], prior_tissue_vector[etissue_index], lambda_lasso, lambda_ridge, sigma);


		//===================================== part#2 =====================================
		// 2.1. snp to cellenv
		para_penalty_lasso_approx_mt(matrix_para_snp_cellenv, matrix_para_dev_snp_cellenv, lambda_snp_cellenv, sigma);
		// 2.2. sparsity (LASSO) for the coefficients from cell env to expression (with the assumption that one gene is only affected by several handful cell env)
		para_penalty_lasso_approx_mt(cube_para_cellenv_gene[etissue_index], cube_para_dev_cellenv_gene[etissue_index], lambda_cellenv_gene, sigma);


		//===================================== part#3 =====================================
		// 3.1. tissue hierarchy regularization;
		// 3.2. or we can simply use group LASSO to encourage the tissue consistency;
		// TODO: test later on
		//
		//
		//


		//===================================== part#4 =====================================
		// 4. penalize the batch variables hashly
		// 4.1. from batch to batch_hidden:
		para_penalty_lasso_approx_mt(matrix_para_batch_batch_hidden, matrix_para_dev_batch_batch_hidden, lambda_batch_batch_hidden, sigma);
		// 4.2. from batch_hidden to gene:
		para_penalty_lasso_approx_mt(matrix_para_batch_hidden_gene, matrix_para_dev_batch_hidden_gene, lambda_batch_hidden_gene, sigma);
	}
	else
	{
		//===================================== part#1 =====================================
		// 1. sparsity of cis- regulation, accompanied by ridge regression, achieved by elastic-net tuned by the prior number, (and the distance prior)
		// TODO: not yet integrated the distance prior information
		para_penalty_cis(cube_para_cis_gene[etissue_index], cube_para_dev_cis_gene[etissue_index], prior_tissue_vector[etissue_index], lambda_lasso, lambda_ridge, sigma);


		//===================================== part#2 =====================================
		// 2.1. snp to cellenv
		para_penalty_lasso_approx(matrix_para_snp_cellenv, matrix_para_dev_snp_cellenv, lambda_snp_cellenv, sigma);
		// 2.2. sparsity (LASSO) for the coefficients from cell env to expression (with the assumption that one gene is only affected by several handful cell env)
		para_penalty_lasso_approx(cube_para_cellenv_gene[etissue_index], cube_para_dev_cellenv_gene[etissue_index], lambda_cellenv_gene, sigma);


		//===================================== part#3 =====================================
		// 3.1. tissue hierarchy regularization;
		// 3.2. or we can simply use group LASSO to encourage the tissue consistency;
		// TODO: test later on
		//
		//
		//


		//===================================== part#4 =====================================
		// 4. penalize the batch variables hashly
		// 4.1. from batch to batch_hidden:
		para_penalty_lasso_approx(matrix_para_batch_batch_hidden, matrix_para_dev_batch_batch_hidden, lambda_batch_batch_hidden, sigma);
		// 4.2. from batch_hidden to gene:
		para_penalty_lasso_approx(matrix_para_batch_hidden_gene, matrix_para_dev_batch_hidden_gene, lambda_batch_hidden_gene, sigma);
	}





	cout << "[@@] leaving the regularization routine..." << endl;
}





// for all parameters in our scope, we do p = p - rate_learner * dp (we have all the components in the right hand, as followed)
void gradient_descent(int etissue_index)
{

	// Apr.2: I'm testing the speedup of mt code of this routine
	int marker_mt = 1;


	cout << "[@@] entering the gradient descent..." << endl;


	if(marker_mt)
	{
		//============================================ pathway#1 ================================================
		//====================== cube_para_cis_gene ==========================
		para_gradient_descent_cis_mt(cube_para_cis_gene[etissue_index], cube_para_dev_cis_gene[etissue_index], rate_learner);


		//============================================ pathway#2 ================================================
		//====================== matrix_para_snp_cellenv ==========================
		para_gradient_descent_mt(matrix_para_snp_cellenv, matrix_para_dev_snp_cellenv, rate_learner);


		//====================== cube_para_cellenv_gene ==========================
		para_gradient_descent_mt(cube_para_cellenv_gene[etissue_index], cube_para_dev_cellenv_gene[etissue_index], rate_learner);


		//============================================ pathway#3 ================================================
		//====================== matrix_para_batch_batch_hidden ==========================
		para_gradient_descent_mt(matrix_para_batch_batch_hidden, matrix_para_dev_batch_batch_hidden, rate_learner);


		//====================== matrix_para_batch_hidden_gene ==========================
		para_gradient_descent_mt(matrix_para_batch_hidden_gene, matrix_para_dev_batch_hidden_gene, rate_learner);

	}
	else
	{
		//============================================ pathway#1 ================================================
		//====================== cube_para_cis_gene ==========================
		para_gradient_descent_cis(cube_para_cis_gene[etissue_index], cube_para_dev_cis_gene[etissue_index], rate_learner);


		//============================================ pathway#2 ================================================
		//====================== matrix_para_snp_cellenv ==========================
		para_gradient_descent(matrix_para_snp_cellenv, matrix_para_dev_snp_cellenv, rate_learner);


		//====================== cube_para_cellenv_gene ==========================
		para_gradient_descent(cube_para_cellenv_gene[etissue_index], cube_para_dev_cellenv_gene[etissue_index], rate_learner);


		//============================================ pathway#3 ================================================
		//====================== matrix_para_batch_batch_hidden ==========================
		para_gradient_descent(matrix_para_batch_batch_hidden, matrix_para_dev_batch_batch_hidden, rate_learner);


		//====================== matrix_para_batch_hidden_gene ==========================
		para_gradient_descent(matrix_para_batch_hidden_gene, matrix_para_dev_batch_hidden_gene, rate_learner);
	}


	cout << "[@@] leaving the gradient descent..." << endl;
}






// calculate the loglikelihood for all the samples in one tissue
float cal_loglike(string etissue)
{
	cout << "now calculating the log-likelihood..." << endl;

	float loglike = 0;
	int etissue_index = etissue_index_map[etissue];

	for(int pos = 0; pos < esample_tissue_rep[etissue].size(); pos++)
	{
		loglike += forward_loglike_testerror(0, etissue, pos, &snp_dosage_list, gene_rpkm_exp, cellenv_hidden_var, batch_var, batch_hidden_var);	// indicator=0 --> loglike; indicator=1 --> testerror
	}

	return loglike;
}




// forward process, accumulated the errors
float forward_loglike_testerror(int indicator, string etissue, int pos, array<float *, NUM_CHR> * dosage_list_pointer, float * expr_con_pointer, float * cellenv_con_pointer, float * batch_list_pointer, float * batch_hidden_con_pointer)
{
	float result = 0;
	int etissue_index = etissue_index_map[etissue];


	//========================================================================
	// prepare all the containers, and the input variables
	//========================================================================
	//========================================================================
	//=================================================================================================================
	//******************************************* loglike or testerror ************************************************
	//=================================================================================================================
	string esample;
	if(indicator == 0)	// training set, loglike
	{
		esample = esample_tissue_rep[etissue][pos];
	}
	else 				// testing set, testerror
	{
		esample = esample_tissue_rep_test[etissue][pos];
	}
	string individual = sample_to_individual(esample);
	cout << "loglike: current sample #" << pos+1 << ": " << esample << endl;


	// make it compatible with the old code
	//array<float *, NUM_CHR> * dosage_list_pointer = &snp_dosage_list;
	//=================================================================================================================
	//******************************************* loglike or testerror ************************************************
	//=================================================================================================================
	vector<float> * expr_list_pointer;
	if(indicator == 0)	// training set, loglike
	{
		expr_list_pointer = &eQTL_tissue_rep[etissue][esample];
	}
	else 				// testing set, testerror
	{
		expr_list_pointer = &eQTL_tissue_rep_test[etissue][esample];
	}
	//float * expr_con_pointer = gene_rpkm_exp;
	//float * cellenv_con_pointer = cellenv_hidden_var;
	//float * batch_list_pointer = batch_var;
	//float * batch_hidden_con_pointer = batch_hidden_var;



	//=================================================== init ============================================================
	// get the: 0. esample and individual; 1. genotype; 2. expression data; 3. batch variables
	// to: 1. forward_backward propagation;
	// genotype dosage data
	//cout << "getting the dosage data for individual #" << individual << endl;
	snp_dosage_load(dosage_list_pointer, individual);  // snp dosage data for one individual across all chromosomes
	// expression rpkm data: eQTL_tissue_rep[etissue][esample]
	//cout << "we have this amount of genes expressed in this individual:" << eQTL_tissue_rep[etissue][esample].size() << endl;
	// and the batch variable for this individual and this sample
	int num_batch_individual = batch_individual[individual].size();
	int index = 0;
	for(int i=0; i<num_batch_individual; i++)
	{
		float value = batch_individual[individual][i];
		batch_list_pointer[index] = value;
		index++;
	}
	int num_batch_sample = batch_sample[esample].size();
	for(int i=0; i<num_batch_sample; i++)
	{
		float value = batch_sample[esample][i];
		batch_list_pointer[index] = value;
		index++;
	}


	//========================================================================
	// forward-propogation (cis-; cell env; batch)
	//========================================================================
	//========================================================================

	// ****************************** [part1] cis- *********************************
	// for cis-, two issues:
	// 1. if this is a XYMT gene, we don't have signal from it's cis- SNPs (not consider currently);
	// 2. we use (gene_cis_index[gene].second - gene_cis_index[gene].first + 1) as the length of the cis- parameter array
	float * expr_con_pointer_cis = (float *)calloc( num_gene, sizeof(float) );
	multi_array_matrix_imcomp(dosage_list_pointer, cube_para_cis_gene[etissue_index], expr_con_pointer_cis);


	// ********************* [part2] cell env relevant parameters *********************
	// from snp to cell env variables
	float * expr_con_pointer_cellenv = (float *)calloc( num_gene, sizeof(float) );
	multi_array_list_matrix(dosage_list_pointer, matrix_para_snp_cellenv, cellenv_con_pointer);

	//$$$$$$$$$$$ perform the activation function here (logistic or something else) $$$$$$$$$$$$
	neuralnet_ac_func(cellenv_con_pointer, num_cellenv);

	// from cell env variables to genes
	multi_array_matrix(cellenv_con_pointer, cube_para_cellenv_gene[etissue_index], expr_con_pointer_cellenv);


	// ********************* [part3] linear or non-linear batches *********************
	float * expr_con_pointer_batch = (float *)calloc( num_gene, sizeof(float) );
	// from original batch to hidden batch
	multi_array_matrix(batch_list_pointer, matrix_para_batch_batch_hidden, batch_hidden_con_pointer);

	//$$$$$$$$$$$ perform the activation function here (logistic or something else) $$$$$$$$$$$$
	neuralnet_ac_func(batch_hidden_con_pointer, num_batch_hidden);

	// from hidden batch to genes
	multi_array_matrix(batch_hidden_con_pointer, matrix_para_batch_hidden_gene, expr_con_pointer_batch);



	// ********************* [end] merge the signal from three pathways here, to expr_con_pointer *********************
	for(long int i=0; i<num_gene; i++)
	{
		expr_con_pointer[i] = expr_con_pointer_cis[i] + expr_con_pointer_cellenv[i] + expr_con_pointer_batch[i];
	}
	free(expr_con_pointer_cis);
	free(expr_con_pointer_cellenv);
	free(expr_con_pointer_batch);



	// error is the thing actually needed
	//=================================================================================================================
	//******************************************* loglike or testerror ************************************************
	//=================================================================================================================
	if(indicator == 0)		// training set, loglike
	{
		float loglike = 0;
		for(long int i=0; i<num_gene; i++)
		{
			float error = expr_con_pointer[i] - (*expr_list_pointer)[i];
			// save the loglike
			loglike -= pow(error, 2.0);
		}
		result = loglike;
	}
	else 					// testing error, testerror
	{
		float testerror = 0;
		for(long int i=0; i<num_gene; i++)
		{
			float error = abs(expr_con_pointer[i] - (*expr_list_pointer)[i]);
			testerror += error;
		}
		result = testerror;
	}


	return result;
}



// testing error for one tissue
// type: MAE or AE (mean absolute error or absolute error)
float cal_testerror(string etissue)
{
	float testerror = 0;

	cout << "now calculating the testing error (for the current tissue)..." << endl;

	int etissue_index = etissue_index_map[etissue];

	for(int pos = 0; pos < esample_tissue_rep_test[etissue].size(); pos++)
	{
		testerror += forward_loglike_testerror(1, etissue, pos, &snp_dosage_list, gene_rpkm_exp, cellenv_hidden_var, batch_var, batch_hidden_var);	// indicator=0 --> loglike; indicator=1 --> testerror
	}


	return testerror;
}


