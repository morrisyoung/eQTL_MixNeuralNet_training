// expression data relevant operations, like loading the expressin matrix into the memory

#include "expression.h"
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


void rpkm_load(unordered_map<string, unordered_map<string, vector<float>>> * eQTL_tissue_rep_pointer, unordered_map<string, string> * eQTL_samples_pointer, vector<string> * gene_list_pointer)
{

	int index = 0;
	unordered_map<int, int> index_rep;



	/* to work on
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
	*/








}