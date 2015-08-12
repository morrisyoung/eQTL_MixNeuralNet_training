// function: processing genotype relevant data (dosage data)

#include "genotype.h"
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


using namespace std;



void snp_info_read(unordered_map<string, snp_info> * hashtable_pointer, int chr)
{
	//======== get all SNPs with their snp_info (count, position) ========
	char filename[100] = "../genotype_185_dosage_matrix_qc/chr";
	char chrom[10];
	sprintf(chrom, "%d", chr);
	strcat(filename, chrom);
	strcat(filename, "/SNP_info.txt");
	//puts("the current file worked on is: ");
	//puts(filename);

	FILE * file_in = fopen(filename, "r");
	if(file_in == NULL)
	{
		fputs("File error\n", stderr); exit (1);
	}

	int input_length = 100;
	char input[input_length];
	char input2[input_length];
	long count = 0;
	while(fgets(input, input_length, file_in) != NULL)
	{
		// the target from this code section: snp (string); count (long); position (long)
		count++;

		strcpy(input2, input);
		char * pos = strstr(input2, " ");
		pos++;
		string snp = strtok(input, " ");
		long position = strtol(pos, NULL, 10);

		snp_info SNP_info;
		SNP_info.count = count;
		SNP_info.position = position;

		(* hashtable_pointer).emplace(snp, SNP_info);
	}

	fclose(file_in);
	//======================================

}





void prune_info_read(vector<string> * pointer1, unordered_map<string, forward_list<pair<string, float>>> * pointer2, int chr)
{
	
	//======== fill (*pointer1) ========
	char filename[100] = "../genotype_185_dosage_matrix_qc/post_prune/chr";
	char chrom[10];
	sprintf(chrom, "%d", chr);
	strcat(filename, chrom);
	strcat(filename, ".prune.in");
	//puts(filename);

	FILE * file_in = fopen(filename, "r");
	if(file_in == NULL)
	{
		fputs("File error\n", stderr); exit (1);
	}

	int input_length = 100;
	char input[input_length];
	while(fgets(input, input_length, file_in) != NULL)
	{
		string snp = input;
		(*pointer1).push_back(snp);
	}

	fclose (file_in);
	//======================================


	//======== fill (*pointer2) ========
	strcpy(filename, "../genotype_185_dosage_matrix_qc/post_prune/chr");
	sprintf(chrom, "%d", chr);
	strcat(filename, chrom);
	strcat(filename, ".post_prune.txt");
	//puts(filename);

	file_in = fopen(filename, "r");
	if(file_in == NULL)
	{
		fputs("File error\n", stderr); exit (1);
	}

	input_length = 10000;
	char input1[input_length];
	while(fgets(input1, input_length, file_in) != NULL)
	{
		trim(input1);

		const char * sep = "\t";
		char * p;
		p = strtok(input1, sep);
		string snp_unpruned = p;

		// initialize the hashtable item
		forward_list<pair<string, float>> list;
		(* pointer2).emplace(snp_unpruned, list);

		int count = 0;
		while(p)
		{
			count++;
			if(count == 1)
			{
				p = strtok(NULL, sep);
				continue;
			}
			// append this pair into the hashed forward_list
			pair<string, float> snp_pair;
			pair_split(p, &snp_pair);
			(* pointer2)[snp_unpruned].push_front(snp_pair);

			p = strtok(NULL, sep);
		}
	}

	fclose (file_in);
	//======================================

}





void dosage_load(int chr, char * individual)
{

	/*
	//======== file reading module ========
	char filename[100] = "../genotype_185_dosage_matrix_qc/chr";
	char chrom[10];
	sprintf(chrom, "%d", chr);
	strcat(filename, chrom);
	strcat(filename, "/SNP_dosage_");
	strcat(filename, individual);
	strcat(filename, ".txt");
	puts(filename);

	FILE * file_in = fopen(filename, "r");
	if(file_in == NULL)
	{
		fputs("File error\n", stderr); exit (1);
	}

	int input_length = 100;
	char input[input_length];
	while(fgets(input, input_length, file_in) != NULL)
	{
		cout << input;
		// do something here
	}

	fclose (file_in);
	//======================================
	*/


}