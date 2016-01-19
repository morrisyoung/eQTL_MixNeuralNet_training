// lib_matrix.h
// function: the matrix class and some basic operations

#ifndef LIB_MATRIX_H
#define LIB_MATRIX_H


#include <iostream>
//#include <sys/types.h>
//#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
//#include <string.h>
#include <string>
#include <string.h>
#include <vector>
#include <math.h>       /* exp */




using namespace std;




// matrix class
class Matrix
{
	long int dimension1;
	long int dimension2;
	vector<float *> matrix;			// element as "float *" other than "array", as we don't pre-know the length

	public:
		//================================ constructor =======================================
		void init(long int value1, long int value2)
		{
			dimension1 = value1;
			dimension2 = value2;

			// the matrix will be initialized as 0 matrix, with "calloc"
			for(long int i=0; i<dimension1; i++)
			{
				float * pointer = (float *)calloc( dimension2, sizeof(float) );
				matrix.push_back(pointer);
			}
			return;
		}
		void init(long int value1, long int value2, vector<float *> vec)
		{
			dimension1 = value1;
			dimension2 = value2;

			// the matrix will be initialized as 0 matrix, with "calloc"
			for(long int i=0; i<dimension1; i++)
			{
				float * pointer = (float *)calloc( dimension2, sizeof(float) );
				matrix.push_back(pointer);
			}

			for(long int i=0; i<dimension1; i++)
			{
				for(long int j=0; j<dimension2; j++)
				{
					matrix[i][j] = vec[i][j];
				}
			}

			return;
		}


		//================================ operations =======================================
		// dimension evaluation
		int get_dimension1()
		{
			return dimension1;
		}
		int get_dimension2()
		{
			return dimension2;
		}

		// clean: set all the values in this matrix to 0
		void clean()
		{
			for(long int i=0; i<dimension1; i++)
			{
				for(long int j=0; j<dimension2; j++)
				{
					matrix[i][j] = 0;
				}
			}

			return;
		}

		// scale: scale the full matrix according to a given factor
		void scale(float factor)
		{
			for(long int i=0; i<dimension1; i++)
			{
				for(long int j=0; j<dimension2; j++)
				{
					matrix[i][j] = matrix[i][j] * factor;
				}
			}

			return;
		}

		// get
		float get(long int pos1, long int pos2)
		{
			return matrix[pos1][pos2];
		}

		// assign
		void assign(long int pos1, long int pos2, float value)
		{
			matrix[pos1][pos2] = value;

			return;
		}

		// add on
		void add_on(long int pos1, long int pos2, float value)
		{
			matrix[pos1][pos2] = matrix[pos1][pos2] + value;

			return;
		}

		// check nan: check whether there is Nan in this Matrix
		int check_nan()
		{
			int flag = 0;

			for(long int i=0; i<dimension1; i++)
			{
				for(long int j=0; j<dimension2; j++)
				{
					float value = matrix[i][j];
					// check nan
					if(isnan(value))
					{
						flag = 1;
						break;
					}
				}
			}

			return flag;
		}


		//================================ destructor =======================================
		void release()
		{
			for(long int i=0; i<dimension1; i++)
			{
				free(matrix[i]);
			}
			return;
		}


};




// matrix (imcomplete, for cis- association parameters) class
class Matrix_imcomp
{
	long int dimension;
	vector<long int> list_length;
	vector<float *> matrix;			// element as "float *" other than "array", as we don't pre-know the length

	public:
		//================================ constructor =======================================
		void init(long int value)
		{
			dimension = value;
			for(long int i=0; i<dimension; i++)
			{
				list_length.push_back(0);
				float * pointer = NULL;
				matrix.push_back(pointer);
			}

			return;
		}

		// initialize each element (empty) with specified length
		void init_element(long int pos, long int length)
		{
			list_length[pos] = length;
			matrix[pos] = (float *)calloc( length, sizeof(float) );

			return;
		}

		// initialize each element (with values) with specified length
		void fill_element(long int pos, long int length, float * list)
		{
			list_length[pos] = length;
			matrix[pos] = (float *)calloc( length, sizeof(float) );
			for(long int i=0; i<length; i++)
			{
				matrix[pos][i] = list[i];
			}

			return;
		}



		//================================ operations =======================================







		//================================ destructor =======================================
		void release()
		{
			for(long int i=0; i<dimension; i++)
			{
				free(matrix[i]);
			}
			return;
		}


};




#endif

// end of lib_matrix.h

