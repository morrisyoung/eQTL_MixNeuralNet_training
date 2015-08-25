// the main optimization routine

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
#include "optimization.h"
#include "global.h"
#include "main.h"



using namespace std;


//====================================== local global variables ========================================
// these variables are specially designed for this routine -- optimization
// need to initialize some local containers:
array<vector<float>, 22> snp_dosage_list;
vector<float> gene_rpkm_exp;
vector<float> cellenv_hidden_var;



//======================================================================================================



void opt_para_init()
{
	puts("opt_para_init..");
	for(int i=0; i<22; i++)
	{
		for(int j=0; j<snp_name_list[i].size(); j++)
		{
			snp_dosage_list[i].push_back(0);
		}
	}

	for(int i=0; i<gene_list.size(); i++)
	{
		gene_rpkm_exp.push_back(0);
	}

	for(int i=0; i<num_cellenv; i++)
	{
		cellenv_hidden_var.push_back(0);
	}

}



void opt_para_release()
{

}




void optimize()
{
	puts("============== entering the optimization routine...");

	opt_para_init();




	// puts("get the dosage data for one individual (test)...");
	// string individual = "GTEX-TKQ1";  // for testing sample (individual)
	// snp_dosage_load(&snp_dosage_list, individual);  // snp dosage data for one individual across all chromosomes
	// puts("snp dosage loading (test) done!");









	// two (2) major parts for this routine:

	// 1. write the procedure to update these parameters in this program (mini-batches)

	// 2. stop in some ways, and save the learned parameters to disk (interface; should be in a trivial way); use another routine to check the predictive precision





//=================================== our parameter space ===================================
// long int num_snp = 0;
// int num_cellenv = 400;
// long int num_gene = 0;
// int num_etissue = 0;
// array<vector<string>, 22> snp_name_list;
// array<vector<long>, 22> snp_pos_list;
// unordered_map<string, unordered_map<string, vector<float>>> eQTL_tissue_rep;  // hashing all eTissues to their actual rep, in which all sample from that tissue is hashed to their rpkm array
// unordered_map<string, string> eQTL_samples;  // hashing all eQTL samples to their tissues
// vector<string> gene_list;  // all genes from the source file
// unordered_map<string, gene_pos> gene_tss;  // TSS for all genes (including those pruned genes)
// unordered_map<string, int> gene_xymt_rep;  // map all the X, Y, MT genes
// vector<vector<float *> *> para_cis_gene;
// vector<float *> para_snp_cellenv;
// vector<vector<float *> *> para_cellenv_gene;
// unordered_map<string, tuple_long> gene_cis_index;  // mapping the gene to cis snp indices (start position and end position in the snp vector)
//============================================================================================




// we probably need two more extra space for * cell env *



//========================================================================
// two step: forward propagation; backward propagation (gradient descent)
//========================================================================
// step#1: ... (cis-; cell env)

// ************** part1: cis- ***************

// ************** part2: cell env relevant parameters **************





//========================================================================
// step#2: ... (cis- parameters;  cell env relevant parameters)

// ************** part1: cis- ***************
// for gene in "gene_list":
// if gene in gene_xymt_rep: continue;
// else:
// int chr = gene_tss[gene]
// long start = gene_cis_index[gene].first;
// long end = gene_cis_index[gene].second;

// for index= start; index<= end; index++:
// dosage = snp_dosage_list[chr-1][index];


// ************** part2: cell env relevant parameters **************
// for gene in gene_list:























	opt_para_release();
	puts("============== leaving the optimization routine...");

}