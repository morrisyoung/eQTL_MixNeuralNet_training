// libfunc_matrix.h

#ifndef LIBFUNC_MATRIX_H
#define LIBFUNC_MATRIX_H


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




using namespace std;




void para_penalty_lasso_approx(Matrix, Matrix, float, float);

void para_penalty_cis(Matrix_imcomp, Matrix_imcomp, vector<vector<float>>>, float, float, float);


void para_gradient_descent(Matrix, Matrix, float);

void para_gradient_descent_cis(Matrix_imcomp, Matrix_imcomp, rate);


void multi_array_matrix(float *, Matrix, float *);





#endif

// end of libfunc_matrix.h

