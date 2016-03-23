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
// NOTE: we have ordered the tissues, and we need to load the parameter files for the ordered tissues
void para_init()
{
	//
	// initializing the memory (pre-fill with 0)
	//
	//==================================== cellular factor pathway =====================================
	//=============== from snp to cell env variables ===============
	for(int i=0; i<num_cellenv; i++)
	{
		//float * p = (float *)calloc( num_snp, sizeof(float) );
		// add the intercept:
		float * p = (float *)calloc( num_snp + 1, sizeof(float) );
		para_snp_cellenv.push_back(p);
	}

	//=============== from cell env variables to genes ===============
	for(int j=0; j<num_etissue; j++)
	{
		vector<float *> vec;
		para_cellenv_gene.push_back(vec);
		for(int i=0; i<num_gene; i++)
		{
			//float * p = (float *)calloc( num_cellenv, sizeof(float) );
			// add the intercept:
			float * p = (float *)calloc( num_cellenv + 1, sizeof(float) );
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
				//float * p = (float *)calloc( amount, sizeof(float) );
				// add the intercept:
				float * p = (float *)calloc( amount + 1, sizeof(float) );
				para_cis_gene[j].push_back(p);
			}
		}
	}

	//==================================== batch effect pathway =====================================
	//=============== from original batch to hidden batch ===============
	for(int i=0; i<num_batch_hidden; i++)
	{
		//float * p = (float *)calloc( num_batch, sizeof(float) );
		// add the intercept:
		float * p = (float *)calloc( num_batch + 1, sizeof(float) );
		para_batch_batch_hidden.push_back(p);
	}

	//=============== from hidden batch to genes ===============
	for(int i=0; i<num_gene; i++)
	{
		//float * p = (float *)calloc( num_batch_hidden, sizeof(float) );
		// add the intercept:
		float * p = (float *)calloc( num_batch_hidden + 1, sizeof(float) );
		para_batch_hidden_gene.push_back(p);
	}






	//
	// initializing their value of those parameters
	//
	// TODO: Don't forget to divide the signal from each pathway into their share
	//		temporarily don't, cause there are some non-linearality in the hidden layer (maybe the initialization is also not so good)
	//
	//==================================== cellular factor pathway =====================================
	//=============== from snp to cell env variables ===============
	// vector<float *> para_snp_cellenv
	char filename[100];
	filename[0] = '\0';
	strcat(filename, file_para_init);
	strcat(filename, "para_snp_cellenv.txt");
	//puts("the current file worked on is: ");
	//puts(filename);

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
	// get the temp rep first
	unordered_map<string, string> map_temp;

	//char filename[100];
	filename[0] = '\0';
	strcat(filename, file_para_init);
	strcat(filename, "etissue_list.txt");
	//puts("the current file worked on is: ");
	//puts(filename);

	//FILE * file_in = fopen(filename, "r");
	file_in = fopen(filename, "r");
	if(file_in == NULL)
	{
		fputs("File error\n", stderr); exit (1);
	}
	while(fgets(input, input_length, file_in) != NULL)
	{
		trim(input);

		const char * sep = "\t";
		char * p;
		p = strtok(input, sep);
		string etissue = p;
		
		p = strtok(NULL, sep);
		string etissue_index = p;

		map_temp[etissue] = etissue_index;
	}
	fclose(file_in);

	// fill the parameter space for each tissue then
	for(int j=0; j<num_etissue; j++)
	{
		string etissue = etissue_list[j];
		string etissue_index = map_temp[etissue];

		char filename[100];
		filename[0] = '\0';
		strcat(filename, file_para_init);
		strcat(filename, "para_cellenv_gene/etissue");
		strcat(filename, etissue_index.c_str());
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
				para_cellenv_gene[j][count][count1] = para;
				p = strtok(NULL, sep);
				count1++;
			}

			count++;
		}
		fclose(file_in);
	}


	//==================================== cis- association pathway =====================================
	// vector<vector<float *>> para_cis_gene
	// map_temp can be used here
	for(int j=0; j<num_etissue; j++)
	{
		// build the rep first, then use it to fill all the parameter space
		unordered_map<string, vector<float>> rep_para_cis_gene;

		string etissue = etissue_list[j];
		string etissue_index = map_temp[etissue];

		char filename[100];
		filename[0] = '\0';
		strcat(filename, file_para_init);
		strcat(filename, "para_cis_gene/etissue");
		strcat(filename, etissue_index.c_str());
		strcat(filename, ".txt");
		//puts("the current file worked on is: ");
		//puts(filename);

		FILE * file_in = fopen(filename, "r");
		if(file_in == NULL)
		{
			fputs("File error\n", stderr); exit(1);
		}
		//int input_length = 100000;
		//char input[input_length];
		while(fgets(input, input_length, file_in) != NULL)
		{
			trim(input);

			const char * sep = "\t";
			char * p;
			p = strtok(input, sep);
			string gene = p;
			vector<float> vec;
			rep_para_cis_gene.emplace(gene, vec);

			int count = 0;
			while(p)
			{
				count++;
				if(count == 1)  // this is the gene
				{
					p = strtok(NULL, sep);
					continue;
				}
				// append this para, and iterate across all samples
				float para = stof(p);
				// NOTE: there may be "nan" in the cis- parameters, as these parameters are from regression analysis
				if(isnan(para))
				{
					para = 0;
				}
				rep_para_cis_gene[gene].push_back(para);

				p = strtok(NULL, sep);
			}
		}
		fclose(file_in);


		// then fill in the same parameters for all tissue types
		for(long i=0; i<gene_list.size(); i++)
		{
			string gene = gene_list[i];
			// will not consider xymt genes
			unordered_map<string, int>::const_iterator got = gene_xymt_rep.find(gene);
			if ( got != gene_xymt_rep.end() )
			{
				continue;
			}
			// for security
			unordered_map<string, vector<float>>::const_iterator got1 = rep_para_cis_gene.find(gene);
			if ( got1 == rep_para_cis_gene.end() )
			{
				continue;
			}

			long first = gene_cis_index[gene].first;  // index
			long second = gene_cis_index[gene].second;  // index
			long amount = second - first + 1;
			//for(int k=0; k<amount; k++)  // NOTE: we will drop the last item, which is the intercept of the multi-linear model
			// add the intercept:
			for(int k = 0; k < amount + 1; k++)  // NOTE: we will drop the last item, which is the intercept of the multi-linear model
			{
				para_cis_gene[j][i][k] = rep_para_cis_gene[gene][k];
			}
			// leaving this gene
		}
	}


	//==================================== batch effect pathway =====================================
	//=============== from original batch to hidden batch ===============
	// vector<float *> para_batch_batch_hidden
	//char filename[100];
	filename[0] = '\0';
	strcat(filename, file_para_init);
	strcat(filename, "para_batch_batch_hidden.txt");
	//puts("the current file worked on is: ");
	//puts(filename);

	file_in = fopen(filename, "r");
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
	//char filename[100];
	filename[0] = '\0';
	strcat(filename, file_para_init);
	strcat(filename, "para_batch_hidden_gene.txt");
	//puts("the current file worked on is: ");
	//puts(filename);

	file_in = fopen(filename, "r");
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


	// release the huge memory used as buffer
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
	for(int i=0; i<num_etissue; i++)
	{
		for(int j=0; j<num_gene; j++)
		{
			free(para_cellenv_gene[i][j]);
		}
	}

	//=============== initialize: vector<float *> para_cis_gene ===============
	for(int i=0; i<num_etissue; i++)
	{
		for(long j=0; j<num_gene; j++)
		{
			free(para_cis_gene[i][j]);
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
	// array<vector<long>, NUM_CHR> snp_pos_list;
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

