// genotype.h
// function: read dosage chunk (for one individual on one chromosome), according to the information from pruning

#ifndef GENOTYPE_H
#define GENOTYPE_H

#include <unordered_map>
#include <string>
#include <vector>
#include <forward_list>
#include <utility>

using namespace std;


// this can also be substituted by pair<long, long>
typedef struct snp_info
{
	long count;
    long position;
}snp_info;


// read snp_info into the specified hashtable
long int snp_info_read();


// read the un-pruned snp list and the association file (pruned snps with the un-pruned snps) into specified vector or hashtable
void prune_info_read(vector<string> *, unordered_map<string, forward_list<pair<string, float>>> *, int);


// read dosage chunk (for one individual on all chromosomes)
void snp_dosage_load(vector<float> *, int, string);



#endif

// end of genotype.h