// test_para_read.h
// function:

#ifndef TEST_PARA_READ_H
#define TEST_PARA_READ_H



using namespace std;



// initializing the parameter space; and fill them with the prepared parameters (training), or learned parameters (testing)
void para_init();


// release all the dynamically allocated memory at the end of the program
void para_release();


// loading and preparing some gene (cis- relevant) mate data
void gene_cis_index_init();




#endif

// end of test_para_read.h