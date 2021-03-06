// optimization.h
// function: the main optimization routine

#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H


#include <string>
#include "global.h"
#include "lib_matrix.h"


using namespace std;






extern array<float *, NUM_CHR> snp_dosage_list;
extern float * gene_rpkm_exp;  // with length "num_gene"
extern float * cellenv_hidden_var;  // with length "num_cellenv"
extern float * batch_var;  // with length "num_batch"
extern float * batch_hidden_var;  // with length "num_batch_hidden"


extern vector<Matrix_imcomp> cube_para_dev_cis_gene;
extern Matrix matrix_para_dev_snp_cellenv;
extern vector<Matrix> cube_para_dev_cellenv_gene;
extern Matrix matrix_para_dev_batch_batch_hidden;
extern Matrix matrix_para_dev_batch_hidden_gene;




typedef struct hierarchy_neighbor
{
	string node;
    float branch;
}hierarchy_neighbor;


extern vector<Matrix_imcomp> cube_para_cis_gene_parent;
extern vector<Matrix> cube_para_cellenv_gene_parent;


extern unordered_map<string, hierarchy_neighbor> hash_leaf_parent;
extern unordered_map<string, vector< hierarchy_neighbor >> hash_internode_neighbor;
extern vector<string> internode_list;
extern unordered_map<string, int> internode_index_map;
extern int num_internode;
extern vector<float> etissue_dis_par_list;




// some assistant components:
// the prior number for each un-pruned snp for regularization (from pruned snps and chromatin states); per etissue, per chromosome, for each snp
// TODO: we also still need to integrate distance prior later on with the following prior information
extern vector<vector<vector<float>>> prior_tissue_vector;
// pairwise phylogenetic distance between etissues
extern vector<vector<float>> tissue_hierarchical_pairwise;


// learning control parameters:
extern int iter_learn_out;  // iteration across all tissues
extern int iter_learn_in;  // iteration across all samples from one tissue
extern int batch_size;
extern float rate_learner;  // the learning rate




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