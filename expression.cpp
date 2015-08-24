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
#include "global.h"


using namespace std;


long int rpkm_load()  // fill in: eQTL_samples; gene_list; eQTL_tissue_rep
{
	puts("load rpkm matrix...");
	//===================================== eQTL_samples ===========================================
	char filename[100] = "../phs000424.v4.pht002743.v4.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_60_samples_train";
	FILE * file_in = fopen(filename, "r");
	if(file_in == NULL)
	{
		fputs("File error\n", stderr); exit (1);
	}
	int input_length = 100000;
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


	//===================================== gene_list; eQTL_tissue_rep ===========================================
	unordered_map<int, string> index_rep;

	filename[0] = '\0';
	strcpy(filename, "../GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_normalized");
	file_in = fopen(filename, "r");
	if(file_in == NULL)
	{
		fputs("File error\n", stderr); exit (1);
	}
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
				// // fill the index_rep, with eQTL_samples
				int index = 0;
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
 					unordered_map<string, string>::const_iterator got = eQTL_samples.find(sample);
  					if ( got != eQTL_samples.end() )
  					{
  						index_rep.emplace(index, sample);
  					}

					p = strtok(NULL, sep);
				}
				break;
			}

			default:
			{
				// fill the gene_list, and eQTL_tissue_rep
				int index = 0;
				trim(input);

				const char * sep = "\t";
				char * p;
				p = strtok(input, sep);
				string gene = p;
				gene_list.push_back(gene);

				while(p)
				{
					index++;
					if(index == 1 || index == 2)
					{
						p = strtok(NULL, sep);
						continue;
					}

 					unordered_map<int, string>::const_iterator got = index_rep.find(index);
  					if ( got != index_rep.end() )
  					{
						char rpkm[100];
						strcpy(rpkm, p);
						float expression = stof(rpkm);
						string sample = index_rep[index];
						string eTissue = eQTL_samples[sample];
						eQTL_tissue_rep[eTissue][sample].push_back(expression);
  					}

					p = strtok(NULL, sep);
				}
				break;
			}
		}


	}
	fclose (file_in);

	return gene_list.size();
}



void tss_load()
{
	puts("loading the tss for genes...");

	int index = 0;

	char filename[100] = "../gencode.v18.genes.patched_contigs.gtf_gene_tss";
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
		// fill in this: unordered_map<string, gene_pos> gene_tss
		trim(input);

		const char * sep = "\t";
		char * p;
		p = strtok(input, sep);
		string gene = p;
		gene_pos tuple;
		gene_tss[gene] = tuple;

		index = 0;
		while(p)
		{
			index++;
			if(index == 1)
			{
				p = strtok(NULL, sep);
				continue;
			}
			if(index == 2)
			{
				// chr
				//char temp[100];
				//strcpy(temp, p);
				long temp = strtol(p, NULL, 10);
				gene_tss[gene].chr = temp;
				p = strtok(NULL, sep);
				continue;
			}
			if(index == 3)
			{
				// tss
				//char temp[100];
				//strcpy(temp, p);
				long temp = strtol(p, NULL, 10);
				gene_tss[gene].tss = temp;
				break;
			}
		}


	}
	fclose (file_in);



}


