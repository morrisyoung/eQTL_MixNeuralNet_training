// genotype.h
// function: read dosage chunk (for one individual on one chromosome), according to the information from pruning

#ifndef GENOTYPE_H
#define GENOTYPE_H

#include <unordered_map>
#include <string>
#include <vector>

using namespace std;


// the MRCA struct for the h1/h11, storing the tMRCA value and the pointed position in list1
typedef struct snp_info
{
	long count;
    long position;
}snp_info;

// the snp (pruned) with its association signal (R^2) with the representative un-pruned snp
typedef struct snp_assoc
{
	string snp;
    float assoc;
}snp_assoc;


// read snp_info into the specified hashtable
void snp_info_read(unordered_map<string, snp_info> *, int);


// read the un-pruned snp list and the association file (pruned snps with the un-pruned snps) into specified vector or hashtable
void prune_info_read(vector<string> *, unordered_map<string, snp_assoc> *, int);


// read dosage chunk (for one individual on one chromosome)
void dosage_load(int, char *);



#endif

// end of genotype.h