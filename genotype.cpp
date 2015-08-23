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



void snp_info_read(vector<string> * vec_name_pointer, vector<long> * vec_pos_pointer, int chr)
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
	while(fgets(input, input_length, file_in) != NULL)
	{
		// the target from this code section: snp (string); position (long)
		char input2[input_length];
		strcpy(input2, input);
		char * pos = strstr(input2, " ");
		pos++;
		string snp = strtok(input, " ");
		long position = strtol(pos, NULL, 10);

		(* vec_name_pointer).push_back(snp);
		(* vec_pos_pointer).push_back(position);
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




void snp_dosage_load(vector<float> * vec_pointer, int chr, string individual)
{

	//======== get all SNPs with their snp_info (count, position) ========
	char filename[100] = "../genotype_185_dosage_matrix_qc/chr";
	char chrom[10];
	sprintf(chrom, "%d", chr);
	strcat(filename, chrom);
	strcat(filename, "/SNP_dosage_");
	char individual1[20];
	StrToCharSeq(individual1, individual);
	strcat(filename, individual1);
	strcat(filename, ".txt");
	//puts("the current file worked on is: ");
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
		trim(input);

		float dosage = stof(input);
		(* vec_pointer).push_back(dosage);

	}

	fclose(file_in);
	//======================================

}