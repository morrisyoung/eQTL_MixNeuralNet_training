// global.h
// function:

#ifndef GLOBAL_H
#define GLOBAL_H

#include "expression.h"
#include "main.h"


using namespace std;


// number of cell environment variables; default value is 400
extern int num_cellenv;

// number of un-pruned snps
extern long int num_snp;

// number of genes
extern long int num_gene;

// mapping the gene to cis snp indices (start position and end position)
extern unordered_map<string, tuple_long> gene_cis_index;

// cis parameter table
extern vector<float *> para_cis_gene;

// snp to cell env variable parameter container
extern vector<float *> para_snp_cellenv;

// cell env variable to gene parameter container
extern vector<float *> para_cellenv_gene;;


// genotype relevant:
extern array<vector<string>, 22> snp_name_list;
extern array<vector<long>, 22> snp_pos_list;

// expression relevant:
extern unordered_map<string, unordered_map<string, vector<float>>> eQTL_tissue_rep;  // hashing all eTissues to their actual rep, in which all sample from that tissue is hashed to their rpkm array
extern unordered_map<string, string> eQTL_samples;  // hashing all eQTL samples to their tissues
extern vector<string> gene_list;  // all genes from the source file
extern unordered_map<string, gene_pos> gene_tss;  // TSS for all genes (including those pruned genes)



#endif

// end of global.h