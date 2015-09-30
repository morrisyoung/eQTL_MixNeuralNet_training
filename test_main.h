// test_main.h
// function:

#ifndef TEST_H
#define TEST_H

#include <unordered_map>
#include <string>
#include <vector>


using namespace std;


// lock the following as a global
extern unordered_map<string, unordered_map<string, vector<float>>> eQTL_tissue_rep_predict;	// the predicted version of above one


// the main function
int main();



#endif

// end of test_main.h