// expression.h
// function: read the expression matrix (rpkm values) into the specified data structure

#ifndef EXPRESSION_H
#define EXPRESSION_H


using namespace std;


// load the expression matrix into memory
void rpkm_load();


// load tss for all genes, from the annotation file
void tss_load();




#endif

// end of expression.h