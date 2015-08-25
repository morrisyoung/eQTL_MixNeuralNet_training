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


using namespace std;


// initializing the parameter space
// TODO we should actually initialize the value of these parameters in some way (from prior knowledge, or some other ways)
void para_init()
{
	// from snp to cell env variables
	long int i;
	for(i=0; i<num_cellenv; i++)
	{
		float * p = (float *)malloc( sizeof(float) * num_snp );
		para_snp_cellenv.push_back(p);
	}
	// from cell env variables to genes
	for(i=0; i<num_gene; i++)
	{
		float * p = (float *)malloc( sizeof(float) * num_cellenv );
		para_cellenv_gene.push_back(p);
	}

	// initialize: vector<float *> para_cis_gene
	for(i=0; i<gene_list.size(); i++)
	{
		string gene = gene_list[i];
		long first = gene_cis_index[gene].first;  // index
		long second = gene_cis_index[gene].second;  // index
		long amount = second - first + 1;
		float * p = (float *)malloc( sizeof(float) * amount );
		para_cis_gene.push_back(p);
	}

}


void para_release()
{
	// from snp to cell env variables
	long int i;
	for(i=0; i<num_cellenv; i++)
	{
		free(para_snp_cellenv[i]);
	}
	// from cell env variables to genes
	for(i=0; i<num_gene; i++)
	{
		free(para_cellenv_gene[i]);
	}
}


// loading and preparing some gene (cis- relevant) mate data
void gene_meta_init()
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
