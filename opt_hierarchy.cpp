// the MLE (maximum likelihood estimate) program given fixed hierarchy

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







void hierarchy()
{
	/*
	// input:
	//		vector<Matrix_imcomp> cube_para_cis_gene;
	//		vector<Matrix> cube_para_cellenv_gene;



	// output: two cubes with all parental variables filled
	//		vector<Matrix_imcomp> cube_para_cis_gene_parent;
	//		vector<Matrix> cube_para_cellenv_gene_parent;
	//		vector<float> hierarchy_dist_parent;			// the distance of this tissue to its parent



	// INPUT: the hierarchy (with branch length), should be able to fill hashtable2 and hashtable4 below
	//		two files:
	//			1. all leaves hashing all parents (directly intermediate nodes) with the branch lengths;
	//			2. all internal nodes hashing children and parents with the branch lengths



	// what to build for the hierarchy computing:
	//		[1]. hashing all the leaves and the root ([0]^n) to their variable array
	//		[2]. hashing all the leaves to their parents (in order to retrieve the parental variable array)
	//		[3]. hashing all the internal nodes to variable array to be filled
	//		[4]. hashing all the internal nodes to their children and parent (with length to them, or the variance), in order to build the computational matrices
	//		[5]. having a bi-directional list for the internal nodes (in order to build and fill in computtional matrices)
	//		[6]. building the tissue distance list (to its parent) based on etissue_list, to return after this program


	// the computing matrices:
	// assuming there are n parental nodes:

	k1_p1	k1_p2	k1_p3	...	k1_pn		p1 		f1(leaves)
	...										p2 		f2(leaves)
	...										.		.
	...									x 	.	=	.
	...										.		.
	...										.		.
	kn_p1	k1_p2	k1_p3	...	k1_pn		pn 		fn(leaves)

	-->
	A x P = B

	-->
	P = A^{-1} x B



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








}


