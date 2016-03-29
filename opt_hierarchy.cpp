// the MLE (maximum likelihood estimate) program given fixed hierarchy
// comments: or you can treat this as expectation of the internal nodes, given the branch lengths; the thing is that we have close-form solution of expectation when it's a joint Gaussian

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unordered_map>
#include <string>
#include <vector>
#include "basic.h"
#include <forward_list>
#include <utility>
#include "genotype.h"
#include "expression.h"
#include "optimization.h"
#include "global.h"
#include "main.h"  // typedef struct tuple_long
#include <math.h>       /* exp */
#include "lib_matrix.h"
#include "opt_hierarchy.h"



using namespace std;



/*
//// input:
//		vector<Matrix_imcomp> cube_para_cis_gene;
//		vector<Matrix> cube_para_cellenv_gene;



//// output: two cubes with all parental variables filled
//		vector<Matrix_imcomp> cube_para_cis_gene_parent;
//		vector<Matrix> cube_para_cellenv_gene_parent;



// INPUT: the hierarchy (with branch length), should be able to fill hashtable2 and hashtable4 below
//		three files (or two practically):
//			1. all leaves hashing all parents (directly intermediate nodes) with the branch lengths;
//			2. all internal nodes hashing children and parents with the branch lengths;
//			3. the root hashing it's child (top node of the tree)

# this is the reference for tree file; detailed data structures see below



// what to build for the hierarchy computing (prepared and renewed):
//		[1]. (not any more) hashing all the leaves and the root ([0]^n) to their variable array
//		[2]. hashing all the leaves to their parents (in order to retrieve the parental variable array)
//		[3]. (not any more) hashing all the internal nodes to variable array to be filled
//		[4]. hashing all the internal nodes to their children and parent (with length to them, or the variance), in order to build the computational matrices
//		[5]. having a bi-directional list for the internal nodes (in order to build and fill in computtional matrices)
//		(6). building the tissue distance list (to its parent) based on etissue_list, to be used by the actual regularization

unordered_map<string, hierarchy_neighbor> hash_leaf_parent;		// --> [2]
unordered_map<string, vector< hierarchy_neighbor >> hash_internode_neighbor;
																// --> [4]
vector<string> internode_list;									// --> [5]
unordered_map<string, int> internode_index_map;					// --> [5]
int num_internode;
vector<float> etissue_dis_par_list;								// --> (6)

with:

typedef struct hierarchy_neighbor
{
	string node;
    float branch;
}hierarchy_neighbor;



// the computing matrices:
// assuming there are n internal nodes:

k1_p1	k1_p2	k1_p3	...	k1_pn		p1 		f1(leaves, root)
...										p2 		f2(leaves, root)
...										.		.
...									x 	.	=	.
...										.		.
...										.		.
kn_p1	k1_p2	k1_p3	...	k1_pn		pn 		fn(leaves, root)

-->
A x P = B

-->
P = A^{-1} x B

// I will make the following variables local
vector<vector<float>> matrix_computation;						// --> A
vector<float> array_computation;								// --> B
vector<float> P;												// --> P



//pseudocode: (some of the parameters are re-usable)
do the following for both cis- regulator (G) and cellular regulator (C):
	for each gene i:
		retrieve [1], with input and INPUT; build [6] with INPUT
		fill in [2], with INPUT
		initialize [3], with INPUT
		fill in [4], with INPUT
		for each position j in len(G) or len(C):
			go over [4], fill in the A and B matrices above
			P = A^{-1} x B
		fill in cube_para_cis_gene_parent[][i][] or cube_para_cellenv_gene_parent[][i][] with [3] and etissue_list

return cube_para_cis_gene_parent, cube_para_cellenv_gene_parent and [6]
*/







// func: inverse A, and do A^{-1} x B^{T} to get P
// input involves the following:
//	vector<vector<float>> *  matrix_computation;						// --> A
//	vector<float> array_computation;								// --> B
//	vector<float> P;												// --> P
void matrix_inv_multiply(vector<vector<float>> * matrix_computation_pointer, vector<float> * array_computation_pointer, vector<float> * P_pointer, int index1, int index2)
{
	// do something
	for(int i=0; i<num_internode; i++)
	{
		string internode = internode_list[i];

		string neighbor1 = hash_internode_neighbor[internode][0].node;
		float branch1 = hash_internode_neighbor[internode][0].branch;
		string neighbor2 = hash_internode_neighbor[internode][1].node;
		float branch2 = hash_internode_neighbor[internode][1].branch;
		string neighbor3 = hash_internode_neighbor[internode][2].node;
		float branch3 = hash_internode_neighbor[internode][2].branch;
		

		//==== check the three, fill in A and B

		// how to get the value for one leaf:
		string etissue;
		int etissue_index = etissue_index_map[etissue];
		float value1 = cube_para_cis_gene[etissue_index].get(index1, index2);
		float value2 = cube_para_cellenv_gene[etissue_index].get(index1, index2);



		//==== inverse A, and multiply B (getting P')


		//==== use P' to fill in P



	}


	return;
}






void hierarchy()
{
	//============= variable preparation =============
	// initialize as 0:
	for(int i=0; i<num_etissue; i++)
	{
		cube_para_cis_gene_parent[i].clean();
		cube_para_cellenv_gene_parent[i].clean();
	}

	//==== what we can utilize:
	//unordered_map<string, hierarchy_neighbor> hash_leaf_parent;		// --> [2]
	//unordered_map<string, vector< hierarchy_neighbor >> hash_internode_neighbor;
																		// --> [4]
	//vector<string> internode_list;									// --> [5]
	//unordered_map<string, int> internode_index_map;					// --> [5]
	//int num_internode;
	//vector<float> etissue_dis_par_list;								// --> (6)


	//==== and the following local ones:
	//vector<vector<float>> matrix_computation;						// --> A
	//vector<float> array_computation;								// --> B
	//vector<float> P;												// --> P
	vector<vector<float>> matrix_computation;						// --> A
	for(int i=0; i<num_internode; i++)
	{
		vector<float> vec;
		matrix_computation.push_back(vec);
		for(int j=0; j<num_internode; j++)
		{
			matrix_computation[i].push_back(0);
		}
	}
	vector<float> array_computation;								// --> B
	for(int i=0; i<num_internode; i++)
	{
		array_computation.push_back(0);
	}
	vector<float> P;												// --> P
	for(int i=0; i<num_internode; i++)
	{
		P.push_back(0);
	}


	//=================================
	//==== for cis- regulator cube ====
	//=================================
	for(int i = 0; i < num_gene; i++)
	{
		int num_regulator = cube_para_cis_gene_parent[0].get_dimension2(i);		// can use any tissue to retrieve the dimension
		for(int j=0; j<num_regulator; j++)
		{
			//============= variable preparation =============
			//vector<vector<float>> matrix_computation;						// --> A
			//vector<float> array_computation;								// --> B
			//vector<float> P;												// --> P
			for(int count1=0; count1<matrix_computation.size(); count1++)
			{
				for(int count2=0; count2<matrix_computation[count1].size(); count2++)
				{
					matrix_computation[count1][count2] = 0;
				}
			}
			for(int count1=0; count1<array_computation.size(); count1++)
			{
				array_computation[count1] = 0;
			}
			for(int count1=0; count1<P.size(); count1++)
			{
				P[count1] = 0;
			}


  			//==== inverse "xxx", do A^{-1} x B^{T} (probably through libraries)
			matrix_inv_multiply(&matrix_computation, &array_computation, &P, i, j);


  			//==== fill in cube_para_cis_gene_parent[x][i][j], as the computation matrices will be re-used
  			// with: hash_leaf_parent, internode_index_map, P
  			for(int count=0; count<num_etissue; count++)
  			{
				string etissue = etissue_list[count];
				string parent = hash_leaf_parent[etissue].node;
				int parent_index = internode_index_map[parent];

				cube_para_cis_gene_parent[count].assign(i, j, P[parent_index]);
			}

		}// end j, the current regulator

	}//end i, the current gene




	//============================================
	//==== for cellular factor regulator cube ====
	//============================================
	for(int i = 0; i < num_gene; i++)
	{
		int num_regulator = cube_para_cellenv_gene_parent[0].get_dimension2();		// can use any tissue to retrieve the dimension
		for(int j=0; j<num_regulator; j++)
		{
			//============= variable preparation =============
			//vector<vector<float>> matrix_computation;						// --> A
			//vector<float> array_computation;								// --> B
			//vector<float> P;												// --> P
			for(int count1=0; count1<matrix_computation.size(); count1++)
			{
				for(int count2=0; count2<matrix_computation[count1].size(); count2++)
				{
					matrix_computation[count1][count2] = 0;
				}
			}
			for(int count1=0; count1<array_computation.size(); count1++)
			{
				array_computation[count1] = 0;
			}
			for(int count1=0; count1<P.size(); count1++)
			{
				P[count1] = 0;
			}


  			//==== inverse "xxx", do A^{-1} x B^{T} (probably through libraries)
			matrix_inv_multiply(&matrix_computation, &array_computation, &P, i, j);


  			//==== fill in cube_para_cellenv_gene_parent[x][i][j], as the computation matrices will be re-used
  			// with: hash_leaf_parent, internode_index_map, P
  			for(int count=0; count<num_etissue; count++)
  			{
				string etissue = etissue_list[count];
				string parent = hash_leaf_parent[etissue].node;
				int parent_index = internode_index_map[parent];

				cube_para_cellenv_gene_parent[count].assign(i, j, P[parent_index]);


			}

		}// end j, the current regulator

	}//end i, the current gene



	return;
}


