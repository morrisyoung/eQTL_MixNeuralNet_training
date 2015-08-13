// expression.h
// function: read the expression matrix (rpkm values) into the specified data structure

#ifndef EXPRESSION_H
#define EXPRESSION_H


#include <unordered_map>
#include <string>
#include <vector>


using namespace std;


// load the expression matrix into memory
void rpkm_load(unordered_map<string, unordered_map<string, vector<float>>> *, unordered_map<string, string> *, vector<string> *);




#endif

// end of expression.h