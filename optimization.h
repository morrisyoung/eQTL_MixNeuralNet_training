// optimization.h
// function: the main optimization routine

#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H


using namespace std;



// initialize some local parameter containers
void opt_para_init();


// release the memory for some dynamically allocated space (if there is)
void opt_para_release();



void forward_backward_prop_batch(string, int, int);



void gradient_descent();



// main routine
void optimize();



#endif

// end of optimization.h