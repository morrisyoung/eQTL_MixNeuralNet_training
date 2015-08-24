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
#include "optimization.h"
#include "global.h"
#include "parameter_init.h"



using namespace std;



// global variables definition and initialization
//===========================================================
long int num_snp = 0;
int cell_env = 400;
long int num_gene = 0;


// genotype relevant:
array<vector<string>, 22> snp_name_list;
array<vector<long>, 22> snp_pos_list;

// expression relevant:
// what we need:
// 1. list of eQTL tissues, hashing all samples with their rpkm value;
// 2. hashed all eQTL samples, for convenience of reading relevant rpkm data from the course file
// 3. array of all genes (assuming all genes in the source file are those to be used)
unordered_map<string, unordered_map<string, vector<float>>> eQTL_tissue_rep;  // hashing all eTissues to their actual rep, in which all sample from that tissue is hashed to their rpkm array
unordered_map<string, string> eQTL_samples;  // hashing all eQTL samples to their tissues
vector<string> gene_list;  // all genes from the source file
unordered_map<string, gene_pos> gene_tss;  // TSS for all genes (including those pruned genes)

// parameter space:





//===========================================================





int main()
{
	cout << "[now enter the program]\n";


	// maybe here accept some command lines
	//
	//
	//
	//
	//
	//
	//
	//
	//


	// yes we need this information to characterize the cis- snps or not, in practical computation
	//==================== prepare the snp information (hashtable: (snp, (count, position))) =====================
	puts("preparing the snp info (index --> snp name and chromosome positions)...");
	num_snp = snp_info_read();
	cout << "there are " << num_snp << " snps totally." << endl;
	//============================================================================================================







	/* permanently

	// thought (Aug.23): we don't need to load all these information into the main program to process; we can simply pre-process them, and
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






	//===================================== prepare the expression matrix ======================================
	num_gene = rpkm_load();
	tss_load();  // gene_tss
	//============================================================================================================



	// initialize all the parameters
	para_init();




	optimize();



	cout << "Optimization done! Please find the results in 'result' folder.\n";
	cout << "[now leave the program]\n";

	return 0;
}