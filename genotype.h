// genotype.h
// function: read dosage chunk (for one individual on one chromosome), according to the information from pruning

#ifndef GENOTYPE_H
#define GENOTYPE_H


using namespace std;


// the MRCA struct for the h1/h11, storing the tMRCA value and the pointed position in list1
typedef struct snp_info
{
	long count;
    long position;
}snp_info;


// read dosage chunk (for one individual on one chromosome)
int dosage_load(int, char *);



#endif

// end of genotype.h