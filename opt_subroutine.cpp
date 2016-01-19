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


		forward_backward(etissue,
						&snp_dosage_list,
						&eQTL_tissue_rep[etissue][esample],
						gene_rpkm_exp,
						cellenv_hidden_var,
						batch_var,
						batch_hidden_var,
						&(para_dev_cis_gene[etissue_index]),
						&(para_dev_cellenv_gene[etissue_index]),
						&para_dev_snp_cellenv,
						&para_dev_batch_hidden_gene,
						&para_dev_batch_batch_hidden);

		// leaving the mini-batch
	}




	// DEBUG: debug the parameter dev, and the current cellenv and hiddenbatch
	para_temp_save_dev(etissue_index);
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
	regularization(etissue);



	//=========================================== Gradient Descent ===========================================
	gradient_descent(etissue);



	cout << "[@@] leaving the forward-backward propagation..." << endl;

}




// property of this function: additive to the total derivative after this round (additive)
// what we need for the following routine:
// dosage list; expression value list; expression list; cellenv list; batch list; batch hidden list; ALL parameter (derivative) containers
void forward_backward(string etissue,
	array<float *, NUM_CHR> * dosage_list_pointer,
	vector<float> * expr_list_pointer,
	float * expr_con_pointer,
	float * cellenv_con_pointer,
	float * batch_list_pointer,
	float * batch_hidden_con_pointer,
	vector<float *> * para_dev_cis_gene_pointer,
	vector<float *> * para_dev_cellenv_gene_pointer,
	vector<float *> * para_dev_snp_cellenv_pointer,
	vector<float *> * para_dev_batch_hidden_gene_pointer,
	vector<float *> * para_dev_batch_batch_hidden_pointer)
{
	// how to map these pointers:
	/*
	dosage_list_pointer --> &snp_dosage_list
	expr_list_pointer --> &eQTL_tissue_rep[etissue][esample]
	batch_list_pointer --> batch_var
	expr_con_pointer  --> gene_rpkm_exp
	cellenv_con_pointer --> cellenv_hidden_var
	batch_hidden_con_pointer --> batch_hidden_var

	para_dev_cis_gene_pointer --> &para_dev_cis_gene[etissue_index]
	para_dev_cellenv_gene_pointer --> &para_dev_cellenv_gene[etissue_index]
	para_dev_snp_cellenv_pointer --> &para_dev_snp_cellenv
	para_dev_batch_hidden_gene_pointer --> &para_dev_batch_hidden_gene
	para_dev_batch_batch_hidden_pointer --> &para_dev_batch_batch_hidden
	*/

	int etissue_index = etissue_index_map[etissue];

	//========================================================================
	// two step: forward propagation (get the function values); backward propagation (get the parameter derivatives)
	//========================================================================
	//========================================================================
	// step#1: forward-propogation (cis-; cell env; batch)
	//========================================================================
	//========================================================================



	// TODO: to get the new data structure
	// ****************************** [part1] cis- *********************************
	// for cis-, two issues:
	// 1. if this is a XYMT gene, jump;
	// 2. we use (gene_cis_index[gene].second - gene_cis_index[gene].first + 1) as the length of the cis- parameter array
	float * expr_con_pointer_cis = (float *)calloc( num_gene, sizeof(float) );
	for(int i=0; i<num_gene; i++)
	{
		//expr_con_pointer_cis[i] = 0;
		string gene = gene_list[i];
		unordered_map<string, int>::const_iterator got = gene_xymt_rep.find(gene);
		if ( got != gene_xymt_rep.end() )
		{
			continue;
		}
		else
		{
			int chr = gene_tss[gene].chr;
			int num = gene_cis_index[gene].second - gene_cis_index[gene].first + 1;
			for(int k=0; k<num; k++)
			{
				int pos = gene_cis_index[gene].first + k;
				float dosage = (*dosage_list_pointer)[chr-1][pos];  // dosage at k position
				expr_con_pointer_cis[i] += dosage * para_cis_gene[etissue_index][i][k];
			}
		}
	}



	// ********************* [part2] cell env relevant parameters *********************
	// from snp to cell env variables
	// TODO: to get the new data structure
	float * expr_con_pointer_cellenv = (float *)calloc( num_gene, sizeof(float) );
	for(int i=0; i<num_cellenv; i++)
	{
		cellenv_con_pointer[i] = 0;
		long count = 0;
		for(int j=0; j<NUM_CHR; j++)  // across all the chromosomes
		{
			int chr = j+1;
			for(long k=0; k<snp_name_list[j].size(); k++)
			{
				float dosage = (*dosage_list_pointer)[j][k];
				cellenv_con_pointer[i] += dosage * para_snp_cellenv[i][count];
				count ++;
			}
		}
	}

	// // DEBUG
	// char filename[100] = "../result_tempdata/var_cellenv_before.txt";
	// para_temp_save_var(cellenv_con_pointer, num_cellenv, filename);

	//$$$$$$$$$$$ perform the activation function here (logistic or something else) $$$$$$$$$$$$
	neuralnet_ac_func(cellenv_con_pointer, num_cellenv);

	// from cell env variables to genes
	multi_array_matrix(cellenv_con_pointer, cube_para_cellenv_gene[etissue_index], expr_con_pointer_cellenv);




	// ********************* [part3] linear or non-linear batches *********************
	float * expr_con_pointer_batch = (float *)calloc( num_gene, sizeof(float) );
	// from original batch to hidden batch
	multi_array_matrix(batch_list_pointer, matrix_para_batch_batch_hidden, batch_hidden_con_pointer);

	// // DEBUG
	// sprintf(filename, "%s", "../result_tempdata/var_batch_hidden_before.txt");
	// para_temp_save_var(batch_hidden_con_pointer, num_batch_hidden, filename);

	//$$$$$$$$$$$ perform the activation function here (logistic or something else) $$$$$$$$$$$$
	neuralnet_ac_func(batch_hidden_con_pointer, num_batch_hidden);

	// from hidden batch to genes
	multi_array_matrix(batch_hidden_con_pointer, matrix_para_batch_hidden_gene, expr_con_pointer_batch);


	//// merge the signal from three pathways here, to expr_con_pointer
	for(long int i=0; i<num_gene; i++)
	{
		expr_con_pointer[i] = expr_con_pointer_cis[i] + expr_con_pointer_cellenv[i] + expr_con_pointer_batch[i];
	}


	// // DEBUG
	// sprintf(filename, "%s", "../result_tempdata/var_expr_exp.txt");
	// para_temp_save_var(expr_con_pointer, num_gene, filename);



	//========================================================================
	//========================================================================
	// step#2: back-propogation (cis-;  cell env; batch)
	//========================================================================
	//========================================================================
	// *********************** [part1] cis- ************************
	// pseudo: (expected rpkm - real rpkm) * genotype
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
			int chr = gene_tss[gene].chr;
			int num = gene_cis_index[gene].second - gene_cis_index[gene].first + 1;
			float diff = expr_con_pointer[i] - (*expr_list_pointer)[i];
			for(int k=0; k<num; k++)
			{
				int pos = gene_cis_index[gene].first + k;
				float dosage = (*dosage_list_pointer)[chr-1][pos];  // dosage at k position
				(*para_dev_cis_gene_pointer)[i][k] += diff * dosage;
			}
		}
	}

	// ***************** [part2] cell env relevant parameters *****************
	// from cell env to genes
	// pseudo: (expected rpkm - real rpkm) * cell_env
	for(int i=0; i<num_gene; i++)
	{
		string gene = gene_list[i];
		float diff = expr_con_pointer[i] - (*expr_list_pointer)[i];
		for(int j=0; j<num_cellenv; j++)
		{
			(*para_dev_cellenv_gene_pointer)[i][j] += diff * cellenv_con_pointer[j];
		}
	}
	// from snp to cell env
	// pseudo: [ \sum w3 * (expected rpkm - real rpkm) ] * g'(w2 * x1) * x1







	// // DEBUG
	// //===========================================================================================
	// //===========================================================================================
 //    FILE * file_out1 = fopen("../temp_data/error_cellenv_1.txt", "w+");
 //    if(file_out1 == NULL)
 //    {
 //        fputs("File error\n", stderr); exit(1);
 //    }
	// FILE * file_out2 = fopen("../temp_data/error_cellenv_2.txt", "w+");
	// if(file_out2 == NULL)
	// {
	// 	fputs("File error\n", stderr); exit(1);
	// }
	// //===========================================================================================
	// //===========================================================================================

	// // DEBUG
	// char buf[1024];

	for(int i=0; i<num_cellenv; i++)
	{
		//
		float temp = 0;
		for(int t=0; t<num_gene; t++)
		{
			temp += para_cellenv_gene[etissue_index][t][i] * (expr_con_pointer[t] - (*expr_list_pointer)[t]);
		}

		// // DEBUG
		// // save the back-propogated errors to the file, and check
		// sprintf(buf, "%f\n", temp);
		// fwrite(buf, sizeof(char), strlen(buf), file_out1);

		//
		temp *= neuralnet_ac_func_dev(cellenv_con_pointer[i]);
		//

		// // DEBUG
		// // save the back-propogated errors to the file, and check
		// sprintf(buf, "%f\n", temp);
		// fwrite(buf, sizeof(char), strlen(buf), file_out2);

		long count = 0;
		for(int j=0; j<NUM_CHR; j++)  // across all the chromosomes
		{
			int chr = j+1;
			for(long k=0; k<snp_name_list[j].size(); k++)
			{
				float dosage = (*dosage_list_pointer)[j][k];
				(*para_dev_snp_cellenv_pointer)[i][count] += temp * dosage;
				count ++;
			}
		}
	}

	// // DEBUG
	// //===========================================================================================
	// //===========================================================================================
	// fclose(file_out1);
	// fclose(file_out2);
	// //===========================================================================================
	// //===========================================================================================




	// ********************* [part3] linear or non-linear batches *********************
	// from hidden batch to genes
	// pseudo: (expected rpkm - real rpkm) * hidden batch var
	for(int i=0; i<num_gene; i++)
	{
		string gene = gene_list[i];
		float diff = expr_con_pointer[i] - (*expr_list_pointer)[i];
		for(int j=0; j<num_batch_hidden; j++)
		{
			(*para_dev_batch_hidden_gene_pointer)[i][j] += diff * batch_hidden_con_pointer[j];
		}
	}
	// from original batch to hidden batch
	// pseudo: [ \sum w5 * (expected rpkm - real rpkm) ] * g'(w4 * x2) * x2






	// // DEBUG
	// //===========================================================================================
	// //===========================================================================================
 //    file_out1 = fopen("../temp_data/error_batch_1.txt", "w+");
 //    if(file_out1 == NULL)
 //    {
 //        fputs("File error\n", stderr); exit(1);
 //    }
	// file_out2 = fopen("../temp_data/error_batch_2.txt", "w+");
	// if(file_out2 == NULL)
	// {
	// 	fputs("File error\n", stderr); exit(1);
	// }
	// FILE * file_out3 = fopen("../temp_data/error_batch_3.txt", "w+");
	// if(file_out3 == NULL)
	// {
	// 	fputs("File error\n", stderr); exit(1);
	// }
	// //===========================================================================================
	// //===========================================================================================





	for(int i=0; i<num_batch_hidden; i++)
	{
		//
		float temp = 0;
		for(int t=0; t<num_gene; t++)
		{
			temp += para_batch_hidden_gene[t][i] * (expr_con_pointer[t] - (*expr_list_pointer)[t]);
		}




		// // DEBUG
		// // save the back-propogated errors to the file, and check
		// sprintf(buf, "%f\n", temp);
		// fwrite(buf, sizeof(char), strlen(buf), file_out1);





		//
		temp *= neuralnet_ac_func_dev(batch_hidden_con_pointer[i]);
		//


		// // DEBUG
		// // save the back-propogated errors to the file, and check
		// sprintf(buf, "%f\n", temp);
		// fwrite(buf, sizeof(char), strlen(buf), file_out2);


		// // DEBUG: save all the batch var and the multiplied error
		// for(int j=0; j<num_batch; j++)
		// {
		// 	float batch_value = batch_list_pointer[j];
		// 	(*para_dev_batch_batch_hidden_pointer)[i][j] += temp * batch_value;


		// 	sprintf(buf, "%f\t%f\t", batch_value, temp * batch_value);
		// 	fwrite(buf, sizeof(char), strlen(buf), file_out3);
		// }
		// fwrite("\n", sizeof(char), 1, file_out3);

		for(int j=0; j<num_batch; j++)
		{
			float batch_value = batch_list_pointer[j];
			(*para_dev_batch_batch_hidden_pointer)[i][j] += temp * batch_value;
		}
	}

	// // DEBUG
	// //===========================================================================================
	// //===========================================================================================
	// fclose(file_out1);
	// fclose(file_out2);
	// fclose(file_out3);
	// //===========================================================================================
	// //===========================================================================================




}







void regularization(string etissue)
{
	cout << "[@@] entering the regularization routine..." << endl;

	// there are several classes of prior knowledge that we need to consider
	// 1. sparsity of cis- regulation, accompanied by ridge regression, achieved by elastic-net tuned by the prior number, and the distance prior
	// 2. sparsity (LASSO) for the coefficients from cell env to expression (with the assumption that one gene is only affected by several handful cell env)
	// 3.1.[TODO] hierarchical regularization tuned by the learned tissue hierarchy
	// 3.2.[TODO] or we can simply use group LASSO to encourage the tissue consistency
	// 4. penalize the batch variables hashly (from batch variables to batch_hidden, and from batch_hidden to genes)
	cout << "adding the regularization items to the derivatives..." << endl;

	int etissue_index = etissue_index_map[etissue];

	//===================================== part#0 =====================================
	// initialize some learning parameters

	// define the sigma here that may be used by L1 regularization
	float sigma = 0.0001;

	// regularization strength lambda:
	// path#1: add the lambda_{LASSO} and lambda_{ridge} for the cis- regularization
	// path#2: add the lambda for cellenv-gene regularization
	// path#3: and the lambda for batch-batch_hidden and batch_hidden-gene
	float lambda_lasso = 1.0;
	float lambda_ridge = 1.0;
	float lambda_snp_cellenv = 1.0;
	float lambda_cellenv_gene = 1.0;
	float lambda_batch_batch_hidden = 1.0;
	float lambda_batch_hidden_gene = 1.0;


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


	cout << "[@@] leaving the regularization routine..." << endl;
}







void gradient_descent(string etissue)
{
	cout << "[@@] entering the gradient descent..." << endl;

	// for all parameters in our scope, we do p = p - rate_learner * dp (we have all the components in the right hand, as followed)

	// parameter containers:
	//vector<vector<float *>> para_cis_gene;
	//vector<float *> para_snp_cellenv;
	//vector<vector<float *>> para_cellenv_gene;
	//vector<float *> para_batch_batch_hidden;
	//vector<float *> para_batch_hidden_gene;

	// parameter derivative containers:
	//vector<vector<float *>> para_dev_cis_gene;
	//vector<float *> para_dev_snp_cellenv;
	//vector<vector<float *>> para_dev_cellenv_gene;
	//vector<float *> para_dev_batch_batch_hidden;
	//vector<float *> para_dev_batch_hidden_gene;

	int etissue_index = etissue_index_map[etissue];


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


	cout << "[@@] leaving the gradient descent..." << endl;
}


