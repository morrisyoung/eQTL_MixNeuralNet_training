// save the predicted rpkm
// we don't need to split all the samples into their tissues; we only need to save all the predicted expression arrays

// standard libraries:
#include <iostream>
#include <sys/types.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <unordered_map>
#include <string.h>
#include <string>
#include <array>
#include <forward_list>
#include <utility>
#include <vector>
#include <sys/time.h>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
// sub-routines:
#include "global.h"
#include "genotype.h"
#include "expression.h"
#include "batch.h"
#include "basic.h"
#include "test_main.h"
#include "test_save.h"




using namespace std;



void predict_save()
{
	cout << "saving the predicted gene expression into file (per eTissue)..." << endl;


	//================================ vector<string> etissue_list ================================
	char filename[100] = "../result_predict/etissue_list.txt";
	//puts("the current file worked on is: ");
	//puts(filename);

    FILE * file_out = fopen(filename, "w+");
    if(file_out == NULL)
    {
        fputs("File error\n", stderr); exit(1);
    }

    for(int i=0; i<etissue_list.size(); i++)
    {
    	string etissue = etissue_list[i];
		char buf[100];
    	sprintf(buf, "%s\t%d\n", etissue.c_str(), i+1);
    	fwrite(buf, sizeof(char), strlen(buf), file_out);
    }
	fclose(file_out);


	//========================================== save samples by tissues ==========================================
	for(int i=0; i<num_etissue; i++)
	{
		string etissue = etissue_list[i];
		// open the file
		char filename[100] = "../result_predict/";
		char temp[10];
		sprintf(temp, "%d", i+1);
		strcat(filename, "etissue");
		strcat(filename, temp);
		strcat(filename, ".txt");
		//puts("the current file worked on is: ");
		//puts(filename);

		FILE * file_out = fopen(filename, "w+");
    	if(file_out == NULL)
    	{
	        fputs("File error\n", stderr); exit(1);
	    }

	    // one line one sample
		for(auto it = eQTL_tissue_rep_predict[etissue].begin(); it != eQTL_tissue_rep_predict[etissue].end(); ++it)
		{
			string esample = it->first;
			char buf[100];
    		sprintf(buf, "%s\t", esample.c_str());
			fwrite(buf, sizeof(char), strlen(buf), file_out);
			for(int j=0; j<num_gene; j++)
			{
				float rpkm = eQTL_tissue_rep_predict[etissue][esample][j];
    			sprintf(buf, "%f\t", rpkm);
				fwrite(buf, sizeof(char), strlen(buf), file_out);
			}
			fwrite("\n", sizeof(char), 1, file_out);
			// leaving this sample
		}
		fclose(file_out);
		// leaving this tissue
	}


}

