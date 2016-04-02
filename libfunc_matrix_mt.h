// libfunc_matrix_mt.h

#ifndef LIBFUNC_MATRIX_MT_H
#define LIBFUNC_MATRIX_MT_H


#include <iostream>
//#include <sys/types.h>
//#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <vector>
#include <unordered_map>
#include "lib_matrix.h"
#include <array>
#include "global.h"




using namespace std;




//==== regularization ====
void * WorkPerThread_penalty_lasso_approx(void *);
void para_penalty_lasso_approx_mt(Matrix, Matrix, float, float);

void * WorkPerThread_penalty_cis(void *);
void para_penalty_cis_mt(Matrix_imcomp, Matrix_imcomp, vector<vector<float>> &, float, float, float);




//==== gradient descent ====
void * WorkPerThread_gradient_descent(void *);
void para_gradient_descent_mt(Matrix, Matrix, float);

void * WorkPerThread_gradient_descent_cis(void *);
void para_gradient_descent_cis_mt(Matrix_imcomp, Matrix_imcomp, float);





#endif

// end of libfunc_matrix_mt.h


