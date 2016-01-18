// some debugger routines for the "optimization" routine (including all the involved subroutines)

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
//#include <string.h>
#include <unordered_map>
#include <string>
#include <vector>
//#include "basic.h"
//#include <forward_list>
//#include <utility>
//#include "genotype.h"
//#include "expression.h"
//#include "optimization.h"
#include "global.h"
#include "main.h"  // typedef struct tuple_long
#include <math.h>       /* exp */
#include "opt_debugger.h"
#include "optimization.h"




using namespace std;




// func:
//		0. this is for one specific tissue
//		1. check whether there are nan for any of the parameter (the actual parameter space, not the parameter dev space); if there is, return 1, otherwise return 0
//		2. also return where is the nan
int para_check_nan(string etissue)
{
	int etissue_index = etissue_index_map[etissue];

	int flag = 0;


	//================================ vector<vector<float *>> para_cis_gene ================================
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
			//for(int k=0; k<num; k++)
			// add the intercept:
			for(int k=0; k<num+1; k++)
			{
				float parameter = para_cis_gene[etissue_index][i][k];
				// check nan
				if(isnan(parameter))
				{
					flag = 1;
					break;
				}
			}
		}
	}
	if(flag == 1)
	{
		cout << "I find Nan in para_cis_gene ..." << endl;
		return flag;
	}



	//================================== vector<float *> para_snp_cellenv ===================================
	for(int i=0; i<num_cellenv; i++)
	{
		//for(long j=0; j<num_snp; j++)
		// add the intercept:
		for(long j=0; j<num_snp+1; j++)
		{
			float parameter = para_snp_cellenv[i][j];
			// check nan
			if(isnan(parameter))
			{
				flag = 1;
				break;
			}
		}
	}
	if(flag == 1)
	{
		cout << "I find Nan in para_snp_cellenv ..." << endl;
		return flag;
	}



	//============================== vector<vector<float *>> para_cellenv_gene ==============================
	for(int i=0; i<num_gene; i++)
	{
		string gene = gene_list[i];
		//for(int j=0; j<num_cellenv; j++)
		// add the intercept:
		for(int j=0; j<num_cellenv+1; j++)
		{
			float parameter = para_cellenv_gene[etissue_index][i][j];
			// check nan
			if(isnan(parameter))
			{
				flag = 1;
				break;
			}
		}
	}
	if(flag == 1)
	{
		cout << "I find Nan in para_cellenv_gene ..." << endl;
		return flag;
	}



	//=============================== vector<float *> para_batch_batch_hidden ===============================
	for(int i=0; i<num_batch_hidden; i++)
	{
		//for(int j=0; j<num_batch; j++)
		// add the intercept:
		for(int j=0; j<num_batch+1; j++)
		{
			float parameter = para_batch_batch_hidden[i][j];
			// check nan
			if(isnan(parameter))
			{
				flag = 1;
				break;
			}
		}
	}
	if(flag == 1)
	{
		cout << "I find Nan in para_batch_batch_hidden ..." << endl;
		return flag;
	}




	//=============================== vector<float *> para_batch_hidden_gene ================================
	for(int i=0; i<num_gene; i++)
	{
		//for(int j=0; j<num_batch_hidden; j++)
		// add the intercept:
		for(int j=0; j<num_batch_hidden+1; j++)
		{
			float parameter = para_batch_hidden_gene[i][j];
			// check nan
			if(isnan(parameter))
			{
				flag = 1;
				break;
			}
		}
	}
	if(flag == 1)
	{
		cout << "I find Nan in para_batch_hidden_gene ..." << endl;
		return flag;
	}


	return flag;
}





// func:
//		save all the parameter dev (five parts, for one specific tissue) into "../result_tempdata/"
void para_temp_save_dev(int etissue_index)
{
	//================================ vector<vector<float *>> para_dev_cis_gene ================================
	char filename[100] = "../result_tempdata/para_dev_cis_gene.txt";
	//puts("the current file worked on is: ");
	//puts(filename);

    FILE * file_out = fopen(filename, "w+");
    if(file_out == NULL)
    {
        fputs("File error\n", stderr); exit(1);
    }

	for(int i=0; i<num_gene; i++)
	{
		string gene = gene_list[i];
		unordered_map<string, int>::const_iterator got = gene_xymt_rep.find(gene);
		if ( got != gene_xymt_rep.end() )
		{
			fwrite("\n", sizeof(char), 1, file_out);
		}
		else
		{
			int num = gene_cis_index[gene].second - gene_cis_index[gene].first + 1;
			for(int k=0; k<num; k++)
			{
				float parameter = para_dev_cis_gene[etissue_index][i][k];
				char buf[1024];
				sprintf(buf, "%f\t", parameter);
				fwrite(buf, sizeof(char), strlen(buf), file_out);
			}
			fwrite("\n", sizeof(char), 1, file_out);
		}
	}
	fclose(file_out);


	//================================== vector<float *> para_dev_snp_cellenv ===================================
	sprintf(filename, "%s", "../result_tempdata/para_dev_snp_cellenv.txt");
	//puts("the current file worked on is: ");
	//puts(filename);

    file_out = fopen(filename, "w+");
    if(file_out == NULL)
    {
        fputs("File error\n", stderr); exit(1);
    }

	for(int i=0; i<num_cellenv; i++)
	{
		for(long j=0; j<num_snp; j++)
		{
			float parameter = para_dev_snp_cellenv[i][j];
			char buf[1024];
			sprintf(buf, "%f\t", parameter);
			fwrite(buf, sizeof(char), strlen(buf), file_out);
		}
		fwrite("\n", sizeof(char), 1, file_out);
		// leaving this cellenv
	}
	fclose(file_out);


	//============================== vector<vector<float *>> para_dev_cellenv_gene ==============================
	sprintf(filename, "%s", "../result_tempdata/para_dev_cellenv_gene.txt");
	//puts("the current file worked on is: ");
	//puts(filename);

	file_out = fopen(filename, "w+");
    if(file_out == NULL)
    {
        fputs("File error\n", stderr); exit(1);
    }

	for(int i=0; i<num_gene; i++)
	{
		string gene = gene_list[i];
		for(int j=0; j<num_cellenv; j++)
		{
			float parameter = para_dev_cellenv_gene[etissue_index][i][j];
			char buf[1024];
			sprintf(buf, "%f\t", parameter);
			fwrite(buf, sizeof(char), strlen(buf), file_out);
			// or:
			// fprintf(file_out, "%s", str);
		}
		fwrite("\n", sizeof(char), 1, file_out);
		// leaving this gene
	}
	fclose(file_out);


	//=============================== vector<float *> para_dev_batch_batch_hidden ===============================
	sprintf(filename, "%s", "../result_tempdata/para_dev_batch_batch_hidden.txt");
	//puts("the current file worked on is: ");
	//puts(filename);

    file_out = fopen(filename, "w+");
    if(file_out == NULL)
    {
        fputs("File error\n", stderr); exit(1);
    }

	for(int i=0; i<num_batch_hidden; i++)
	{
		for(int j=0; j<num_batch; j++)
		{
			float parameter = para_dev_batch_batch_hidden[i][j];
			char buf[1024];
			sprintf(buf, "%f\t", parameter);
			fwrite(buf, sizeof(char), strlen(buf), file_out);
		}
		fwrite("\n", sizeof(char), 1, file_out);
		// leaving this batch_hidden
	}
	fclose(file_out);


	//=============================== vector<float *> para_dev_batch_hidden_gene ================================
	sprintf(filename, "%s", "../result_tempdata/para_dev_batch_hidden_gene.txt");
	//puts("the current file worked on is: ");
	//puts(filename);

    file_out = fopen(filename, "w+");
    if(file_out == NULL)
    {
        fputs("File error\n", stderr); exit(1);
    }

	for(int i=0; i<num_gene; i++)
	{
		for(int j=0; j<num_batch_hidden; j++)
		{
			float parameter = para_dev_batch_hidden_gene[i][j];
			char buf[1024];
			sprintf(buf, "%f\t", parameter);
			fwrite(buf, sizeof(char), strlen(buf), file_out);
		}
		fwrite("\n", sizeof(char), 1, file_out);
		// leaving this gene
	}
	fclose(file_out);


	return;
}



void para_temp_save_cellenv()
{
	//======================== cellenv variables ========================
	char filename[100];
	sprintf(filename, "%s", "../result_tempdata/var_cellenv.txt");
	//puts("the current file worked on is: ");
	//puts(filename);

    FILE * file_out = fopen(filename, "w+");
    if(file_out == NULL)
    {
        fputs("File error\n", stderr); exit(1);
    }
	for(int i=0; i<num_cellenv; i++)
	{
		float parameter = cellenv_hidden_var[i];
		char buf[1024];
		sprintf(buf, "%f\n", parameter);
		fwrite(buf, sizeof(char), strlen(buf), file_out);
	}
	fclose(file_out);

	return;
}



void para_temp_save_hiddenbatch()
{
	//======================== batch hidden variables ========================
	char filename[100];
	sprintf(filename, "%s", "../result_tempdata/var_batch_hidden.txt");
	//puts("the current file worked on is: ");
	//puts(filename);

    FILE * file_out = fopen(filename, "w+");
    if(file_out == NULL)
    {
        fputs("File error\n", stderr); exit(1);
    }
	for(int i=0; i<num_batch_hidden; i++)
	{
		float parameter = batch_hidden_var[i];
		char buf[1024];
		sprintf(buf, "%f\n", parameter);
		fwrite(buf, sizeof(char), strlen(buf), file_out);
	}
	fclose(file_out);

	return;
}


