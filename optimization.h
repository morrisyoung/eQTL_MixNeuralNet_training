// optimization.h
// function: the main optimization routine

#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H


#include <string>


using namespace std;



extern array<vector<float>, 22> snp_dosage_list;
extern vector<float> gene_rpkm_exp;  // with length "num_gene"
extern vector<float> cellenv_hidden_var;  // with length "num_cellenv"
extern vector<float> batch_var;  // with length "num_batch"
extern vector<float> batch_hidden_var;  // with length "num_batch_hidden"


// parameter derivative containers:
extern vector<vector<float *>> para_dev_cis_gene;
extern vector<float *> para_dev_snp_cellenv;
extern vector<vector<float *>> para_dev_cellenv_gene;
extern vector<float *> para_dev_batch_batch_hidden;
extern vector<float *> para_dev_batch_hidden_gene;


// some assistant components:
// the prior number for each un-pruned snp for regularization (from pruned snps and chromatin states)
// per etissue, per chromosome, for each snp
// we still need to integrate distance prior later on with the following prior information
extern vector<vector<vector<float>>> snp_prior_list;  // the prior number for each un-pruned snp for regularization (from pruned snps and chromatin states)
// pairwise phylogenetic distance between etissues
extern vector<vector<float>> tissue_hierarchical_pairwise;


// learning control parameters:
extern int iter_learn_out;  // iteration across all tissues
extern int iter_learn_in;  // iteration across all samples from one tissue
extern int batch_size;
extern int rate_learner;  // the learning rate




// load all the cis- snp prior information (tissue specific) from file outside
void opt_snp_prior_load();


// load the pairwise tissue hierarchy from prepared file outside
void opt_tissue_hierarchy_load();


// initialize some local parameter containers
void opt_para_init();


// release the memory for some dynamically allocated space (if there is)
void opt_para_release();



// main entrance
void optimize();




#endif

// end of optimization.h