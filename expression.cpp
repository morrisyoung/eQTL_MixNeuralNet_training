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
	unordered_map<int, string> index_rep;

	char filename[100] = "../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_normalized";
	FILE * file_in = fopen(filename, "r");
	if(file_in == NULL)
	{
		fputs("File error\n", stderr); exit (1);
	}
	int input_length = 100000;
	char input[input_length];
	int count = 0;
	while(fgets(input, input_length, file_in) != NULL)
	{
		count++;
		switch(count)
		{
			case 1:
			{
				break;
			}

			case 2:
			{
				break;
			}

			case 3:
			{
				// // fill the index_rep, with (* eQTL_samples_pointer)
				// all samples are from eTissues, as we have preprocessed the expression file
				index = 0;
				trim(input);

				const char * sep = "\t";
				char * p;
				p = strtok(input, sep);

				while(p)
				{
					index++;
					if(index == 1 || index == 2)
					{
						p = strtok(NULL, sep);
						continue;
					}

					string sample = p;
 					//unordered_map<string, string>::const_iterator got = (* eQTL_samples_pointer).find(sample);
  					//if ( got != (* eQTL_samples_pointer).end() )
  					//{
  					//	index_rep.emplace(index, sample);
  					//}
					index_rep.emplace(index, sample);

					p = strtok(NULL, sep);
				}
				break;
			}

			default:
			{
				// fill the (* gene_list_pointer), and (* eQTL_tissue_rep_pointer)
				index = 0;
				trim(input);

				const char * sep = "\t";
				char * p;
				p = strtok(input, sep);
				string gene = p;
				(* gene_list_pointer).push_back(gene);

				while(p)
				{
					index++;
					if(index == 1 || index == 2)
					{
						p = strtok(NULL, sep);
						continue;
					}

 					//unordered_map<int, string>::const_iterator got = index_rep.find(index);
  					//if ( got != index_rep.end() )
  					//{
					//	char rpkm[100];
					//	strcpy(rpkm, p);
					//	float expression = stof(rpkm);
					//	string sample = index_rep[index];
					//	string eTissue = (* eQTL_samples_pointer)[sample];
					//	(* eQTL_tissue_rep_pointer)[eTissue][sample].push_back(expression);
  					//}
					char rpkm[100];
					strcpy(rpkm, p);
					float expression = stof(rpkm);
					string sample = index_rep[index];
					string eTissue = (* eQTL_samples_pointer)[sample];
					(* eQTL_tissue_rep_pointer)[eTissue][sample].push_back(expression);

					p = strtok(NULL, sep);
				}
				break;
			}
		}


	}
	fclose (file_in);


}