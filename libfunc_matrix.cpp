// this contains the operations on the matrix type



#include "libfunc_matrix.h"
#include "lib_matrix.h"
#include <math.h>





using namespace std;





// func: add the LASSO penalty term in the matrix_para_dev; destroy type
//		the matrix_para is used to calculate the derivative term, that will be added to the matrix_para_dev
void para_penalty_lasso_approx(Matrix matrix_para, Matrix matrix_para_dev, float lambda, float sigma)
{
	long int dimension1 = matrix_para.get_dimension1();
	long int dimension2 = matrix_para.get_dimension2();



	for(int i=0; i<dimension1; i++)
	{
		for(int j=0; j<dimension2; j++)
		{
			/// the value of current para beta:
			float beta = matrix_para.get(i, j);
			/// the derivative of the beta:
			float derivative = beta / sqrt (beta * beta + sigma);  // this is an approximation of the LASSO regularization
			/// and the value of its derivative should be added with that derivative term for regularization:
			matrix_para_dev.add_on(i, j, lambda * derivative);
		}
	}

	return;
}




// func: do the gradient descent for this parameter matrix
//		with: beta = beta - rate * dev
void para_gradient_descent(Matrix matrix_para, Matrix matrix_para_dev, float rate)
{

	long int dimension1 = matrix_para.get_dimension1();
	long int dimension2 = matrix_para.get_dimension2();

	for(long int i=0; i<dimension1; i++)
	{
		for(long int j=0; j<dimension2; j++)
		{
			float dev = matrix_para_dev.get(i, j);
			matrix_para.add_on(i, j, - rate * dev);
		}
	}


	return;
}


