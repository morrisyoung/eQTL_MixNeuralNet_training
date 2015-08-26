// optimization.h
// function: the main optimization routine

#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H


#include <string>


using namespace std;



// initialize some local parameter containers
void opt_para_init();


// release the memory for some dynamically allocated space (if there is)
void opt_para_release();



// transform the sample ID (like "GTEX-R55E-0826-SM-2TC5M") into individual ID (here is the first 9 digits)
string sample_to_individual(string);



// the main optimization routine: the forward_backward propagation, and the gradient descent
void forward_backward_prop_batch(string, int, int);
void gradient_descent();



// main entrance
void optimize();




#endif

// end of optimization.h