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
#include "global.h"
#include "main.h"



using namespace std;


long int snp_info_read()
{
	long int num = 0;

	int i;
	for(i=0; i<22; i++)
	{
		int chr = i+1;
		vector<string> vec1;
		vector<long> vec2;
		snp_name_list[i] = vec1;
		snp_pos_list[i] = vec2;

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

			snp_name_list[i].push_back(snp);
			snp_pos_list[i].push_back(position);
		}

		fclose(file_in);
		//======================================

		num += snp_name_list[i].size();
	}

	return num;
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