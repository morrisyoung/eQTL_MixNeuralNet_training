// function: processing genotype relevant data (dosage data)

#include "genotype.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unordered_map>
#include <string>

using namespace std;

// other routines: processing the information from pruning




int dosage_load(int chr, char * individual)
{

	/*
	//======== file reading module ========
	char filename[100] = "../genotype_185_dosage_matrix/chr";
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



	//======== get all SNPs with their snp_info (count, position) ========
	char filename[100] = "../genotype_185_dosage_matrix/chr";
	char chrom[10];
	sprintf(chrom, "%d", chr);
	strcat(filename, chrom);
	strcat(filename, "/SNP_info.txt");
	puts("the current file worked on is: ");
	puts(filename);

	FILE * file_in = fopen(filename, "r");
	if(file_in == NULL)
	{
		fputs("File error\n", stderr); exit (1);
	}

	int input_length = 100;
	char input[input_length];
	char input2[input_length];
	long count = 0;
	unordered_map<string, snp_info> hashtable;
	while(fgets(input, input_length, file_in) != NULL)
	{
		count++;

		strcpy(input2, input);
		char * pos = strstr(input2, " ");
		pos++;
		string snp = strtok(input, " ");
		long position = strtol(pos, NULL, 10);

		snp_info SNP_info;
		SNP_info.count = count;
		SNP_info.position = position;
		pair<string, snp_info> pair (snp, SNP_info);

		hashtable.insert(pair);

		// snp (string); count (long); position (long)
		//cout << snp << "+" << count << "+" << position << "\n";

	}

	fclose(file_in);
	//======================================

	// test the hashtable here
	// for (auto& x: hashtable)
	// {
	// 	cout << x.first << ":" << x.second.count << "+" << x.second.position << endl;
	// }











	return 0;
}