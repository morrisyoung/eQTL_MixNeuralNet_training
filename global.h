// global.h
// function:

#ifndef GLOBAL_H
#define GLOBAL_H


using namespace std;


// number of cell environment variables; default value is 400
extern int cell_env;

// genotype relevant:
extern array<vector<string>, 22> snp_name_list;
extern array<vector<long>, 22> snp_pos_list;

// expression relevant:
extern unordered_map<string, unordered_map<string, vector<float>>> eQTL_tissue_rep;  // hashing all eTissues to their actual rep, in which all sample from that tissue is hashed to their rpkm array
extern unordered_map<string, string> eQTL_samples;  // hashing all eQTL samples to their tissues
extern vector<string> gene_list;  // all genes from the source file



#endif

// end of global.h