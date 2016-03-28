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
//		[1]. hashing all the leaves and the root ([0]^n) to their variable array
//		[2]. hashing all the leaves to their parents (in order to retrieve the parental variable array)
//		[3]. hashing all the internal nodes to variable array to be filled
//		[4]. hashing all the internal nodes to their children and parent (with length to them, or the variance), in order to build the computational matrices
//		[5]. having a bi-directional list for the internal nodes (in order to build and fill in computtional matrices)
//		(6). building the tissue distance list (to its parent) based on etissue_list, to be used by the actual regularization

unordered_map<string, vector<float>> hash_array_leaf_root;		// --> [1]
unordered_map<string, hierarchy_neighbor> hash_leaf_parent;		// --> [2]
unordered_map<string, vector<float>> hash_array_internode;		// --> [3]
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

vector<vector<float>> matrix_computation;						// --> A
vector<float> array_computation;								// --> B



//pseudocode:
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







void hierarchy()
{

	//============= variable =============
	// initialize as 0:
	cube_para_cis_gene_parent.clean();
	cube_para_cellenv_gene_parent.clean();


// --> has prepared
	// --> need to renew

	//unordered_map<string, vector<float>> hash_array_leaf_root;		// --> [1]
//	unordered_map<string, hierarchy_neighbor> hash_leaf_parent;		// --> [2]
	//unordered_map<string, vector<float>> hash_array_internode;		// --> [3]
//	unordered_map<string, vector< hierarchy_neighbor >> hash_internode_neighbor;
																	// --> [4]
//	vector<string> internode_list;									// --> [5]
//	unordered_map<string, int> internode_index_map;					// --> [5]
//	int num_internode;
//	vector<float> etissue_dis_par_list;								// --> (6)

	//vector<vector<float>> matrix_computation;						// --> A
	//vector<float> array_computation;								// --> B



	for(int i = 0; i < num_gene; i++)
	{
		string gene = gene_list[i];

		for( auto it = hash_array_leaf_root.begin(); it != hash_array_leaf_root.end(); ++it )
		{
			string etissue = it->first;
			for(int i = 0; i < it->second.size(); i++)
			{
				(it->second)[i] = 0;
			}
		}

		for( auto it = hash_array_internode.begin(); it != hash_array_internode.end(); ++it )
		{
			for(int i = 0; i < it->second.size(); i++)
			{
				(it->second)[i] = 0;
			}
		}






	}













}

