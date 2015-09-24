// initializing all the parameters

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
#include "parameter_init.h"
#include "main.h"
#include "expression.h"
#include "basic.h"


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
	// initializing their value of those parameters
	//
	// TODO: Don't forget to divide the signal from each pathway into their share
	//		temporarily don't, cause there are some non-linearality in the hidden layer (maybe the initialization is also not so good)
	//
	//==================================== cellular factor pathway =====================================
	//=============== from snp to cell env variables ===============
	// vector<float *> para_snp_cellenv
	char filename[100] = "../result_init/para_init_train_snp_cellenv.txt";
	FILE * file_in = fopen(filename, "r");
	if(file_in == NULL)
	{
		fputs("File error\n", stderr); exit (1);
	}
	int input_length = 50000000;
	char input[input_length];
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
	for(int j=0; j<num_etissue; j++)
	{
		char filename[100] = "../result_init/para_init_train_cellenv_gene.txt";
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
	// build the rep first, then use it to fill all the parameter space
	unordered_map<string, vector<float>> rep_para_cis_gene;
	char filename0[100] = "../result_init/para_init_train_cis.txt";
	file_in = fopen(filename0, "r");
	if(file_in == NULL)
	{
		fputs("File error\n", stderr); exit (1);
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
			rep_para_cis_gene[gene].push_back(para);

			p = strtok(NULL, sep);
		}
	}
	fclose(file_in);
	// then fill in the same parameters for all tissue types
	for(int j=0; j<num_etissue; j++)
	{
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
			for(int k=0; k<amount; k++)
			{
				para_cis_gene[j][i][k] = rep_para_cis_gene[gene][k];
			}
			// leaving this gene
		}
	}
	//==================================== batch effect pathway =====================================
	//=============== from original batch to hidden batch ===============
	// vector<float *> para_batch_batch_hidden
	char filename1[100] = "../result_init/para_init_train_batch_batch_hidden.txt";
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
	char filename2[100] = "../result_init/para_init_train_batch_hidden_gene.txt";
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


void beta_prior_fill()
{
	//===================================== part#2: get the beta table from files =====================================
	// we can simply build a double hashing table to map the genes to snps to beta's, the later on waling through them and fill in those prior beta's
	// see below:
	unordered_map<string, unordered_map<string, unordered_map<string, float>>> beta_rep;  // we have tissue specificity for this
	//========= get the index map
	unordered_map<string, string> index_map;  // for temporary usage
	char filename[100] = "../GTEx_Analysis_V4_eQTLs/etissue_list.txt";
	FILE * file_in = fopen(filename, "r");
	if(file_in == NULL)
	{
		fputs("File error\n", stderr); exit (1);
	}
	int input_length = 1000;
	char input[input_length];
	while(fgets(input, input_length, file_in) != NULL)
	{
		trim(input);

		const char * sep = "\t";
		char * p;
		p = strtok(input, sep);
		string eTissue = p;

		int count = 0;
		while(p)
		{
			count++;
			if(count == 1)  // this is the eTissue
			{
				p = strtok(NULL, sep);
				continue;
			}
			if(count == 2)  // this is the index
			{
				string index = p;
				index_map[eTissue] = index;
				break;
			}
		}
	}
	fclose (file_in);

	//========= get the index map
	// get the prior beta for each eTissue (if there is), for all gene cis- snps
	for( auto it = index_map.begin(); it != index_map.end(); ++it )
	{
		string eTissue = it->first;
		string index = it->second;

		unordered_map<string, unordered_map<string, float>> map;
		beta_rep[eTissue] = map;

		char filename[100] = "../GTEx_Analysis_V4_eQTLs/etissue";
		char temp[10];
		StrToCharSeq(temp, index);
		strcat(filename, temp);
		strcat(filename, ".beta");
		//puts("the current file worked on is: ");
		//puts(filename);

		FILE * file_in = fopen(filename, "r");
		if(file_in == NULL)
		{
			fputs("File error\n", stderr); exit (1);
		}

		int input_length = 1000;
		char input[input_length];
		while(fgets(input, input_length, file_in) != NULL)
		{
			trim(input);

			const char * sep = "\t";
			char * p;
			p = strtok(input, sep);
			string gene = p;
			// check whether to add this gene to the rep
			unordered_map<string, unordered_map<string, float>>::const_iterator got = beta_rep[eTissue].find(gene);
			if( got == beta_rep[eTissue].end() )
			{
				unordered_map<string, float> map;
				beta_rep[eTissue][gene] = map;
			}
			string snp;

			int count = 0;
			while(p)
			{
				count++;
				if(count == 1)  // this is the gene
				{
					p = strtok(NULL, sep);
					continue;
				}
				if(count == 2)  // this is the snp
				{
					snp = p;
					p = strtok(NULL, sep);
					continue;
				}
				if(count == 3)  // this is the beta
				{
					char temp[100];
					strcpy(temp, p);
					float beta = stof(temp);
					beta_rep[eTissue][gene][snp] = beta;
					break;
				}
			}
		}
		fclose(file_in);
	}


	//===================================== part#2: filling the beta's =====================================
	//unordered_map<string, unordered_map<string, unordered_map<string, float>>> beta_rep;
	for( auto it = beta_rep.begin(); it != beta_rep.end(); ++it )
	{
		string eTissue = it->first;
		int eTissue_index = etissue_index_map[eTissue];  // re-map those etissues into their order (reversed hashing above)

		for(int i=0; i<num_gene; i++)
		{
			string gene = gene_list[i];

			// check whether this is a xymt gene
			unordered_map<string, int>::const_iterator got = gene_xymt_rep.find(gene);
			if( got != gene_xymt_rep.end() )
			{
				continue;
			}

			// otherwise, we should fill in the prior cis- coefficients
			int chr = gene_tss[gene].chr;
			long start = gene_cis_index[gene].first;
			long end = gene_cis_index[gene].second;
			long amount = end - start + 1;
			for(long j=0; j<amount; j++)
			{
				long pos = j + start;
				string snp = snp_name_list[chr-1][pos];
				// check whether snp has prior beta
				unordered_map<string, float>::const_iterator got = beta_rep[eTissue][gene].find(snp);
				if( got != beta_rep[eTissue][gene].end() )
				{
					float beta = beta_rep[eTissue][gene][snp];
					para_cis_gene[eTissue_index][i][j] = beta;
				}
				else
				{
					continue;
				}

			}

			// end filling this gene
		}

		// end filling this tissue
	}

}

