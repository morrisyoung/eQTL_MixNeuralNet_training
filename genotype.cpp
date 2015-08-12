// function: processing genotype relevant data (dosage data)

#include "genotype.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unordered_map>
#include <string>
#include <vector>

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





void prune_info_read(vector<string> * pointer1, unordered_map<string, snp_assoc> * pointer2, int chr)
{
	// fill (*pointer1)
	
	//======== file reading module ========
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



	// fill (*pointer2)











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