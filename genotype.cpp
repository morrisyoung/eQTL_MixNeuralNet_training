// function: processing genotype relevant data (dosage data)

#include "genotype.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

// other routines: processing the information from pruning




int dosage_load(int chr, char * individual)
{

	int input_length = 100;
	char input[input_length];
	char filename[100] = "../genotype_185_dosage_matrix/chr1/SNP_dosage_GTEX-TKQ1.txt";

	FILE * file_in = fopen(filename, "r");
	if(file_in == NULL){fputs("File error\n", stderr); exit (1);}

	while(fgets(input, input_length, file_in) != NULL)
	{
		cout << input;
		cout << "haha";
	}
	fclose (file_in);





	return 0;
}