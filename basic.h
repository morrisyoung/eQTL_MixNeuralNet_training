// basic.h
// function: some basic functions that maybe used by some parts of the whole program

#ifndef BASIC_H
#define BASIC_H

#include <utility>
#include <string>


using namespace std;


// trimming functions
void rtrim(char *);
void ltrim(char *);
void trim(char *);


// split the pair of (snp, R^2) and return pair variable
void pair_split(char *, pair<string, float> *);


// transform the string to char sequence
void StrToCharSeq(char *, string);



#endif

// end of basic.h