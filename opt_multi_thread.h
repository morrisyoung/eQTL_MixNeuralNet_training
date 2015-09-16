// opt_multi_thread.h
// function:

#ifndef OPT_MULTI_THREAD_H
#define OPT_MULTI_THREAD_H


#include <vector>
#include <array>
#include <string>



using namespace std;



// each thread should have such a local parameter space
typedef struct package_thread
{
	// meta data
	int id;
	string etissue;
	int pos_start;
	int num_esample;

	// variable containers
	array<float *, 22> snp_dosage_list;
	float * gene_rpkm_exp;  	// with length "num_gene"
	float * cellenv_hidden_var; // with length "num_cellenv"
	float * batch_var;			// with length "num_batch"
	float * batch_hidden_var;	// with length "num_batch_hidden"

	// parameter containers
	vector<float *> para_dev_cis_gene;
	vector<float *> para_dev_snp_cellenv;
	vector<float *> para_dev_cellenv_gene;
	vector<float *> para_dev_batch_batch_hidden;
	vector<float *> para_dev_batch_hidden_gene;

}package_dev;


void package_alloc(package_dev *);


void package_free(package_dev *);


void * WorkPerThread(void *);


// the control platform for multi-threading program
void opt_mt_control(string, int, int);



// aggregation
void aggregation_init(string);
void aggregation_add(package_dev *, string);
void aggregation_ave(string);




#endif

// end of opt_multi_thread.h

