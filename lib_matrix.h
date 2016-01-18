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
//m#include <string>
#include <string.h>
#include <vector>



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






#endif

// end of lib_matrix.h

