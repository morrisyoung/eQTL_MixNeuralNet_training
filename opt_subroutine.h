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
void gradient_descent();



#endif

// end of opt_subroutine.h