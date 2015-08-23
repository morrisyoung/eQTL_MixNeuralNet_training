/*
the outline of the entire program:

what we should have at hand by now:
1. genotype: grouped by chromosomes, and split for different individuals; read in on-demand fashion;
2. expression: small enough to be fit in memory; should be loaded into memory immediately after initializing the program



1. pipeline for processing the genotype data and the expression date (querying and iterating, in a mini-batch manner);
2. after getting the data, do the stochastic gradient descent algorithm;
3. pay attention to the data structure used for storing all the parameters (coefficients);
4. find a way to terminate the optimization process;

*/


#include <iostream>
#include <sys/types.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include "genotype.h"
#include <unordered_map>
#include <string.h>
#include <string>
#include <array>
#include <forward_list>
#include <utility>
#include "basic.h"
#include "expression.h"



using namespace std;


int main()
{
	cout << "This is the entrance of the program...\n";


	// sub-routine to be constructed: get the actual genotype from (currently) files
	//int chr = 1;
	//char individual[20] = "GTEX-TKQ1";
	//dosage_load(chr, individual);










	/* temporarily
	// yes we need this information to characterize the cis- snps or not, in practical computation


	//==================== prepare the snp information (hashtable: (snp, (count, position))) =====================
	puts("preparing the snp info...");
	array<unordered_map<string, snp_info>, 22> snp_info_list;

	int i;
	for(i=0; i<22; i++)
	{
		int chr = i+1;
		unordered_map<string, snp_info> hashtable;
		snp_info_list[i] = hashtable;
		snp_info_read(&snp_info_list[i], chr);
	}
	puts("snp info preparation done!");

	//===== test the hashtable here: 1. the length; 2. a sample =====
	//cout << (snp_info_list[0]).size() << endl;
	//cout << (snp_info_list[0]).at("rs2073814").count << "+" << (snp_info_list[0]).at("rs2073814").position << endl;
	//cout << (snp_info_list[5]).size() << endl;
	// for (auto& x: hashtable)
	// {
	// 	cout << x.first << ":" << x.second.count << "+" << x.second.position << endl;
	// }
	//===============================================================
	//============================================================================================================
	*/










	/* temporarily
	// feedback: we don't need to load all these information into the main program to process; we can simply pre-process them, and
	//			directly load the pre-processed prior information into this main program (to keep coding simple)
	// so: we need another source file -- each snp with its prior score value (integrated from other associated pruned snps, and chromatin states prior knowledge)



	//===================================== prepare the pruning information ======================================
	puts("preparing the pruning (association of pruned with un-pruned snps) info...");
	array<vector<string>, 22> snp_unpruned_list;
	array<unordered_map<string, forward_list<pair<string, float>>>, 22> snp_prune_assoc_list;

	int i;  // remvoe later on
	for(i=0; i<22; i++)
	{
		int chr = i+1;
		vector<string> snp_unpruned;
		unordered_map<string, forward_list<pair<string, float>>> snp_prune_assoc;

		snp_unpruned_list[i] = snp_unpruned;
		snp_prune_assoc_list[i] = snp_prune_assoc;

		prune_info_read(&snp_unpruned_list[i], &snp_prune_assoc_list[i], chr);
	}
	puts("pruning info preparation done!");

	//===== test the hashtable here: 1. the length; 2. a sample =====
	//forward_list<pair<string, float>> list = snp_prune_assoc_list[2]["rs83616"];
	//for ( auto it = list.begin(); it != list.end(); ++it )
	//	cout << (*it).first << "+" << (*it).second << endl;
	//cout << snp_prune_assoc_list[5].size() << endl;

	//cout << (snp_unpruned_list[3]).size() << endl;
	//cout << (snp_unpruned_list[4])[2] << endl;

	// for (auto& x: hashtable)
	// {
	// 	cout << x.first << ":" << x.second.count << "+" << x.second.position << endl;
	// }
	//===============================================================
	//============================================================================================================


	*/









	/* temporarily

	//===================================== prepare the expression matrix ======================================
	// what we need:
	// 1. list of eQTL tissues, hashing all samples with their rpkm value;
	// 2. hashed all eQTL samples, for convenience of reading relevant rpkm data from the course file; (so we'll need a index list to pick up relevant rpkm values)
	// 3. array of all genes (assuming all genes in the source file are those to be used)
	// specifically:
	unordered_map<string, unordered_map<string, vector<float>>> eQTL_tissue_rep;  // hashing all eTissues to their actual rep, in which all sample from that tissue is hashed to their rpkm array
	unordered_map<string, string> eQTL_samples;  // hashing all eQTL samples to their tissues
	vector<string> gene_list;  // all genes from the source file

	// then we need to initialize some of them before reading the rpkm file
	// eQTL_tissue_rep, eQTL_samples --> "phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples"
	char filename[100] = "../phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples";
	FILE * file_in = fopen(filename, "r");
	if(file_in == NULL)
	{
		fputs("File error\n", stderr); exit (1);
	}
	int input_length = 10000;
	char input[input_length];
	while(fgets(input, input_length, file_in) != NULL)
	{
		trim(input);

		const char * sep = "\t";
		char * p;
		p = strtok(input, sep);
		string eTissue = p;
		unordered_map<string, vector<float>> rep;
		eQTL_tissue_rep.emplace(eTissue, rep);

		int count = 0;
		while(p)
		{
			count++;
			if(count == 1)  // this is the eTissue
			{
				p = strtok(NULL, sep);
				continue;
			}

			// append this sample, and iterate across all samples
			string sample = p;
			vector<float> list;
			eQTL_tissue_rep[eTissue].emplace(sample, list);
			eQTL_samples.emplace(sample, eTissue);

			p = strtok(NULL, sep);
		}

	}
	fclose (file_in);


	rpkm_load(&eQTL_tissue_rep, &eQTL_samples, &gene_list);
	// // simple testing:
	// //gene_list_pointer: size of gene_list; a sample
	// cout << gene_list.size() << endl;
	// cout << gene_list[2] << endl;
	// // eQTL_tissue_rep_pointer:
	// cout << eQTL_tissue_rep["Thyroid"]["GTEX-SN8G-1526-SM-4DM79"].size() << endl;
	// cout << eQTL_tissue_rep["Thyroid"]["GTEX-SN8G-1526-SM-4DM79"][3] << endl;
	// //============================================================================================================

	*/








	//optimization();

	cout << "Optimization done! Please find the results in 'result' folder.\n";

	return 0;
}