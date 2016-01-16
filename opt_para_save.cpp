// function: save all the parameters into file (after this iteration), and test them or interpret them later on
// this is mainly used for testing the training code

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unordered_map>
#include <string>
#include <vector>
#include <utility>
#include "global.h"
#include <cstring>
#include "expression.h"
#include "opt_para_save.h"




using namespace std;




void para_inter_save(int iteration)
{
	cout << "saving the parameters into file (after this iteration)..." << endl;
	char temp[10];


	// write the etissue_list into a file
	//================================ vector<string> etissue_list ================================
	char filename[100] = "../result_inter/etissue_list.txt";
	//puts("the current file worked on is: ");
	//puts(filename);

    FILE * file_out = fopen(filename, "w+");
    if(file_out == NULL)
    {
        fputs("File error\n", stderr); exit(1);
    }

    int count = 0;
    for(int i=0; i<etissue_list.size(); i++)
    {
    	count ++;
    	string etissue = etissue_list[i];
		char buf[100];
    	sprintf(buf, "%s\t%d\n", etissue.c_str(), count);
    	fwrite(buf, sizeof(char), strlen(buf), file_out);
    }
	fclose(file_out);



	//================================ vector<vector<float *>> para_cis_gene ================================
	// this is tissue specific
	for(int i=0; i<num_etissue; i++)
	{
		string etissue = etissue_list[i];
		int etissue_index = i;

		char filename[100] = "../result_inter/para_cis_gene/";
		strcat(filename, "etissue");

		// eTissue #
		sprintf(temp, "%d", i+1);
		strcat(filename, temp);
		//

		strcat(filename, "_");

		// iteration #
		sprintf(temp, "%d", iteration);
		strcat(filename, temp);
		//

		strcat(filename, ".txt");
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
				//for(int k=0; k<num; k++)
				// add the intercept:
				for(int k=0; k<num+1; k++)
				{
					float parameter = para_cis_gene[etissue_index][i][k];
					char buf[1024];
					sprintf(buf, "%f\t", parameter);
					fwrite(buf, sizeof(char), strlen(buf), file_out);
				}
				fwrite("\n", sizeof(char), 1, file_out);
			}
		}
		fclose(file_out);
		// leaving this etissue
	}



	//================================== vector<float *> para_snp_cellenv ===================================
	sprintf(filename, "%s", "../result_inter/para_snp_cellenv_");

	// iteration #
	sprintf(temp, "%d", iteration);
	strcat(filename, temp);
	//

	strcat(filename, ".txt");
	//puts("the current file worked on is: ");
	//puts(filename);


    file_out = fopen(filename, "w+");
    if(file_out == NULL)
    {
        fputs("File error\n", stderr); exit(1);
    }

	for(int i=0; i<num_cellenv; i++)
	{
		//for(long j=0; j<num_snp; j++)
		// add the intercept:
		for(long j=0; j<num_snp+1; j++)
		{
			float parameter = para_snp_cellenv[i][j];
			char buf[1024];
			sprintf(buf, "%f\t", parameter);
			fwrite(buf, sizeof(char), strlen(buf), file_out);
		}
		fwrite("\n", sizeof(char), 1, file_out);
		// leaving this cellenv
	}
	fclose(file_out);


	//============================== vector<vector<float *>> para_cellenv_gene ==============================
	// this is tissue specific
	for(int i=0; i<num_etissue; i++)
	{
		string etissue = etissue_list[i];
		int etissue_index = i;

		char filename[100] = "../result_inter/para_cellenv_gene/";
		strcat(filename, "etissue");

		// eTissue #
		sprintf(temp, "%d", i+1);
		strcat(filename, temp);
		//

		strcat(filename, "_");

		// iteration #
		sprintf(temp, "%d", iteration);
		strcat(filename, temp);
		//

		strcat(filename, ".txt");
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
			//for(int j=0; j<num_cellenv; j++)
			// add the intercept:
			for(int j=0; j<num_cellenv+1; j++)
			{
				float parameter = para_cellenv_gene[etissue_index][i][j];
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
		// leaving this etissue
	}


	//=============================== vector<float *> para_batch_batch_hidden ===============================
	sprintf(filename, "%s", "../result_inter/para_batch_batch_hidden_");

	// iteration #
	sprintf(temp, "%d", iteration);
	strcat(filename, temp);
	//
	
	strcat(filename, ".txt");
	//puts("the current file worked on is: ");
	//puts(filename);


    file_out = fopen(filename, "w+");
    if(file_out == NULL)
    {
        fputs("File error\n", stderr); exit(1);
    }

	for(int i=0; i<num_batch_hidden; i++)
	{
		//for(int j=0; j<num_batch; j++)
		// add the intercept
		for(int j=0; j<num_batch+1; j++)
		{
			float parameter = para_batch_batch_hidden[i][j];
			char buf[1024];
			sprintf(buf, "%f\t", parameter);
			fwrite(buf, sizeof(char), strlen(buf), file_out);
		}
		fwrite("\n", sizeof(char), 1, file_out);
		// leaving this batch_hidden
	}
	fclose(file_out);


	//=============================== vector<float *> para_batch_hidden_gene ================================
	sprintf(filename, "%s", "../result_inter/para_batch_hidden_gene_");

	// iteration #
	sprintf(temp, "%d", iteration);
	strcat(filename, temp);
	//

	strcat(filename, ".txt");
	//puts("the current file worked on is: ");
	//puts(filename);


    file_out = fopen(filename, "w+");
    if(file_out == NULL)
    {
        fputs("File error\n", stderr); exit(1);
    }

	for(int i=0; i<num_gene; i++)
	{
		//for(int j=0; j<num_batch_hidden; j++)
		// add the intercept:
		for(int j=0; j<num_batch_hidden+1; j++)
		{
			float parameter = para_batch_hidden_gene[i][j];
			char buf[1024];
			sprintf(buf, "%f\t", parameter);
			fwrite(buf, sizeof(char), strlen(buf), file_out);
		}
		fwrite("\n", sizeof(char), 1, file_out);
		// leaving this gene
	}
	fclose(file_out);



	cout << "all parameters have been saved into files (after this iteration)..." << endl;

}


