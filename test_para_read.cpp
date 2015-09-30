// initializing all the parameters from the learned results

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unordered_map>
#include <string>
#include <vector>
#include <forward_list>
#include <utility>
#include "global.h"
#include "test_para_read.h"
#include "expression.h"
#include "genotype.h"
#include "basic.h"
#include "test_main.h"




using namespace std;


// initializing the parameter space
void para_init()
{
	//
	// initializing the memory
	//
	//==================================== cellular factor pathway =====================================
	//=============== from snp to cell env variables ===============
	for(int i=0; i<num_cellenv; i++)
	{
		float * p = (float *)calloc( num_snp, sizeof(float) );
		para_snp_cellenv.push_back(p);
	}

	//=============== from cell env variables to genes ===============
	for(int j=0; j<num_etissue; j++)
	{
		vector<float *> vec;
		para_cellenv_gene.push_back(vec);
		for(int i=0; i<num_gene; i++)
		{
			float * p = (float *)calloc( num_cellenv, sizeof(float) );
			para_cellenv_gene[j].push_back(p);
		}
	}

	//==================================== cis- association pathway =====================================
	//=============== initialize: vector<float *> para_cis_gene ===============
	for(int j=0; j<num_etissue; j++)
	{
		vector<float *> vec;
		para_cis_gene.push_back(vec);
		for(long i=0; i<gene_list.size(); i++)
		{
			string gene = gene_list[i];
			unordered_map<string, int>::const_iterator got = gene_xymt_rep.find(gene);
			if ( got != gene_xymt_rep.end() )
			{
				float * p = NULL;  // null pointer
				para_cis_gene[j].push_back(p);
				continue;
			}
			else
			{
				long first = gene_cis_index[gene].first;  // index
				long second = gene_cis_index[gene].second;  // index
				long amount = second - first + 1;
				float * p = (float *)calloc( amount, sizeof(float) );
				para_cis_gene[j].push_back(p);
			}
		}
	}

	//==================================== batch effect pathway =====================================
	//=============== from original batch to hidden batch ===============
	for(int i=0; i<num_batch_hidden; i++)
	{
		float * p = (float *)calloc( num_batch, sizeof(float) );
		para_batch_batch_hidden.push_back(p);
	}

	//=============== from hidden batch to genes ===============
	for(int i=0; i<num_gene; i++)
	{
		float * p = (float *)calloc( num_batch_hidden, sizeof(float) );
		para_batch_hidden_gene.push_back(p);
	}





	//
	// fill in the parameter values
	//
	//==================================== cellular factor pathway =====================================
	//=============== from snp to cell env variables ===============
	// vector<float *> para_snp_cellenv
	char filename[100] = "../result/para_snp_cellenv.txt";
	FILE * file_in = fopen(filename, "r");
	if(file_in == NULL)
	{
		fputs("File error\n", stderr); exit (1);
	}
	long input_length = 5000000000;
	char * input = (char *)malloc( sizeof(char) * input_length );
	int count = 0;
	while(fgets(input, input_length, file_in) != NULL)
	{
		trim(input);

		const char * sep = "\t";
		char * p;
		p = strtok(input, sep);

		long count1 = 0;
		while(p)
		{
			float para = stof(p);
			para_snp_cellenv[count][count1] = para;
			p = strtok(NULL, sep);
			count1++;
		}

		count++;
	}
	fclose(file_in);
	//=============== from cell env variables to genes ===============
	// vector<vector<float *>> para_cellenv_gene
	for(int i=0; i<num_etissue; i++)
	{
		string etissue = etissue_list[i];
		int etissue_index = i;

		char filename[100] = "../result/para_cellenv_gene/";
		char temp[10];
		sprintf(temp, "%d", i+1);
		strcat(filename, "etissue");
		strcat(filename, temp);
		strcat(filename, ".txt");
		//puts("the current file worked on is: ");
		//puts(filename);

		FILE * file_in = fopen(filename, "r");
		if(file_in == NULL)
		{
			fputs("File error\n", stderr); exit (1);
		}
		int input_length = 100000;
		char input[input_length];
		int count = 0;
		while(fgets(input, input_length, file_in) != NULL)
		{
			trim(input);

			const char * sep = "\t";
			char * p;
			p = strtok(input, sep);

			int count1 = 0;
			while(p)
			{
				float para = stof(p);
				para_cellenv_gene[i][count][count1] = para;
				p = strtok(NULL, sep);
				count1++;
			}

			count++;
		}
		fclose(file_in);
	}

	//==================================== cis- association pathway =====================================
	// vector<vector<float *>> para_cis_gene
	for(int i=0; i<num_etissue; i++)
	{
		string etissue = etissue_list[i];
		int etissue_index = i;

		char filename[100] = "../result/para_cis_gene/";
		char temp[10];
		sprintf(temp, "%d", i+1);
		strcat(filename, "etissue");
		strcat(filename, temp);
		strcat(filename, ".txt");
		//puts("the current file worked on is: ");
		//puts(filename);

	    FILE * file_out = fopen(filename, "r");
	    if(file_out == NULL)
	    {
	        fputs("File error\n", stderr); exit(1);
	    }
		int input_length = 100000;
		char input[input_length];
		int count = 0;
		while(fgets(input, input_length, file_in) != NULL)
		{
			trim(input);

			const char * sep = "\t";
			char * p;
			p = strtok(input, sep);

			int count1 = 0;
			while(p)
			{
				float para = stof(p);
				para_cis_gene[etissue_index][count][count1] = para;
				count1++;
				p = strtok(NULL, sep);
			}

			count++;
		}
		fclose(file_in);
		// leaving this etissue
	}

	//==================================== batch effect pathway =====================================
	//=============== from original batch to hidden batch ===============
	// vector<float *> para_batch_batch_hidden
	char filename1[100] = "../result/para_batch_batch_hidden.txt";
	file_in = fopen(filename1, "r");
	if(file_in == NULL)
	{
		fputs("File error\n", stderr); exit (1);
	}
	//int input_length = 100000;
	//char input[input_length];
	count = 0;
	while(fgets(input, input_length, file_in) != NULL)
	{
		trim(input);

		const char * sep = "\t";
		char * p;
		p = strtok(input, sep);

		int count1 = 0;
		while(p)
		{
			float para = stof(p);
			para_batch_batch_hidden[count][count1] = para;
			p = strtok(NULL, sep);
			count1++;
		}

		count++;
	}
	fclose(file_in);
	//=============== from hidden batch to genes ===============
	// vector<float *> para_batch_hidden_gene
	char filename2[100] = "../result/para_batch_hidden_gene.txt";
	file_in = fopen(filename2, "r");
	if(file_in == NULL)
	{
		fputs("File error\n", stderr); exit (1);
	}
	//int input_length = 100000;
	//char input[input_length];
	count = 0;
	while(fgets(input, input_length, file_in) != NULL)
	{
		trim(input);

		const char * sep = "\t";
		char * p;
		p = strtok(input, sep);

		int count1 = 0;
		while(p)
		{
			float para = stof(p);
			para_batch_hidden_gene[count][count1] = para;
			p = strtok(NULL, sep);
			count1++;
		}

		count++;
	}
	fclose(file_in);



	// release the huge memory used as buff
	free(input);

}




void para_release()
{
	//=============== from snp to cell env variables ===============
	for(int i=0; i<num_cellenv; i++)
	{
		free(para_snp_cellenv[i]);
	}

	//=============== from cell env variables to genes ===============
	for(int j=0; j<num_etissue; j++)
	{
		for(int i=0; i<num_gene; i++)
		{
			free(para_cellenv_gene[j][i]);
		}
	}

	//=============== initialize: vector<float *> para_cis_gene ===============
	for(int j=0; j<num_etissue; j++)
	{
		for(long i=0; i<gene_list.size(); i++)
		{
			free(para_cis_gene[j][i]);
		}
	}

	//=============== from original batch to hidden batch ===============
	for(int i=0; i<num_batch_hidden; i++)
	{
		free(para_batch_batch_hidden[i]);
	}

	//=============== from hidden batch to genes ===============
	for(int i=0; i<num_gene; i++)
	{
		free(para_batch_hidden_gene[i]);
	}

}



// loading and preparing some gene (cis- relevant) mate data
void gene_cis_index_init()
{
	// initialize the following one:
	// unordered_map<string, tuple_long> gene_cis_index;  // mapping the gene to cis snp indices (start position and end position)
	//
	// with the following:
	// genotype relevant:
	// array<vector<long>, 22> snp_pos_list;
	// expression relevant:
	// vector<string> gene_list;
	// unordered_map<string, gene_pos> gene_tss;  // TSS for all genes (including those pruned genes)

	int i;
	int j;
	for(i=0; i<gene_list.size(); i++)
	{
		string gene = gene_list[i];
		tuple_long tuple;
		gene_cis_index[gene] = tuple;

		int chr = gene_tss[gene].chr;
		long tss = gene_tss[gene].tss;

		int flag1 = 0;
		int flag2 = 0;
		for(j=0; j<snp_pos_list[chr-1].size(); j++)
		{
			if(flag1 == 0)
			{
				if((snp_pos_list[chr-1][j] - tss > -1000000 && snp_pos_list[chr-1][j] - tss < 1000000) || (tss - snp_pos_list[chr-1][j] > -1000000 && tss - snp_pos_list[chr-1][j] < 1000000))
				{
					gene_cis_index[gene].first = j;
					flag1 = 1;
				}
			}

			if(flag1 == 1 && flag2 == 0)
			{
				if(snp_pos_list[chr-1][j] - tss > -1000000 && snp_pos_list[chr-1][j] - tss < 1000000 || (tss - snp_pos_list[chr-1][j] > -1000000 && tss - snp_pos_list[chr-1][j] < 1000000))
				{
					gene_cis_index[gene].second = j;
				}
				else
				{
					flag2 = 1;
				}
			}

			if(flag1 == 1 && flag2 == 1)break;

		}

	}

}
