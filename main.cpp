/*
the outline of the entire program:

what we should have at hand by now:
1. genotype: grouped by chromosomes, and split for different individuals; read on-demand;
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
#include <string>
#include <array>
#include <forward_list>
#include <utility>



using namespace std;


int main()
{
	cout << "This is the entrance of the program...\n";


	//int chr = 1;
	//char individual[20] = "GTEX-TKQ1";
	//dosage_load(chr, individual);






	/* temporarily

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














	// rpkm_load();








	//optimization();

	cout << "Optimization done! Please find the results in 'result' folder.\n";

	return 0;
}