// expression.h
// function: read the expression matrix (rpkm values) into the specified data structure

#ifndef EXPRESSION_H
#define EXPRESSION_H


using namespace std;


typedef struct gene_pos
{
	int chr;
	long tss;
}gene_pos;



// load the expression matrix into memory
long int gene_rpkm_load(char *, char *);  // fill in: eQTL_samples; gene_list; eQTL_tissue_rep


// load tss for all genes, from the annotation file
void gene_tss_load();


// map all the X, Y, MT genes to unordered_map<string, int> gene_xymt_rep
void gene_xymt_load();



#endif

// end of expression.h