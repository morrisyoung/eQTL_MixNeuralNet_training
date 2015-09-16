// opt_subroutine.h
// function: subroutines of optimization

#ifndef OPT_SUBROUTINE_H
#define OPT_SUBROUTINE_H


#include <string>


using namespace std;



// transform the sample ID (like "GTEX-R55E-0826-SM-2TC5M") into individual ID (here is the first 9 digits)
string sample_to_individual(string);


// the main optimization routine: the forward_backward propagation, and the gradient descent
void forward_backward_prop_batch(string, int, int);
void forward_backward(string, array<float *, 22> *, vector<float> *, float *, float *, float *, float *, vector<float *> *, vector<float *> *, vector<float *> *, vector<float *> *, vector<float *> *);
void regularization(string);
void gradient_descent(string);




#endif

// end of opt_subroutine.h
