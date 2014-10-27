/*  Copyright 2014 - UVSQ
    Authors list: Loïc Thébault, Eric Petit

    This file is part of the D&C library.

    D&C library is free software: you can redistribute it and/or modify it under the
    terms of the GNU Lesser General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later version.

    D&C library is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
    PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License along with
    the D&C library. If not, see <http://www.gnu.org/licenses/>. */

#include <cilk/cilk.h>
#include <pthread.h>

#include "tools.h"
#include "coloring.h"
#include "permutations.h"
#include "partitioning.h"
#include "tree.h"

// The D&C tree and the permutations are global in order to persist from one call
// to the library to another
tree_t *treeHead = nullptr;
int *elemPerm = nullptr, *nodePerm = nullptr;

#ifdef TREE_CREATION

// Mutex to avoid race condition in merge permutations
pthread_mutex_t mergeMutex = PTHREAD_MUTEX_INITIALIZER;

// Compute the edge interval for each leaf of the D&C tree
void compute_edge_intervals (tree_t &tree, int *nodeToNodeRow, int *elemToNode,
                             int nbNodes)
{
    // If current node is a leaf
    if (tree.left == nullptr && tree.right == nullptr) {
        // Get the first and last edges of the leaf
        int firstNode = tree.firstCSR, lastNode = tree.lastCSR;
        tree.firstCSR = nodeToNodeRow[firstNode];
        tree.lastCSR  = nodeToNodeRow[lastNode+1] - 1;
    }
    else {
        // Left & right recursion
        cilk_spawn
        compute_edge_intervals (*tree.left, nodeToNodeRow, elemToNode, nbNodes);
        compute_edge_intervals (*tree.right, nodeToNodeRow, elemToNode, nbNodes);

        // Synchronization
        cilk_sync;
    }
}

// Wrapper used to get the root of the D&C tree before computing the edge intervals
// for CSR reset
void DC_finalize_tree (int *nodeToNodeRow, int *elemToNode, int nbNodes)
{
    compute_edge_intervals (*treeHead, nodeToNodeRow, elemToNode, nbNodes);
}

// Initialize a node of the D&C tree
void init_dc_tree (tree_t &tree, int firstElem, int lastElem, int nbSepElem,
				   int firstNode, int lastNode, bool isLeaf)
{
	tree.firstElem = firstElem;
	tree.lastElem  = lastElem - nbSepElem;
	tree.lastSep   = lastElem;
    tree.firstCSR  = firstNode;
    tree.lastCSR   = lastNode;
    tree.vecOffset = 0;
	tree.left   = nullptr;
	tree.right  = nullptr;
	tree.sep    = nullptr;

	if (isLeaf == false) {
  		tree.left  = new tree_t;
  		tree.right = new tree_t;
		if (nbSepElem > 0) {
			tree.sep = new tree_t;
		}
	}
}

// Create element partition & count left & separator elements
void create_elem_part (int *elemPart, int *nodePart, int *elemToNode, int nbElem,
                       int dimElem, int separator, int offset, int *nbLeftElem,
                       int *nbSepElem)
{
	for (int i = 0; i < nbElem; i++) {
		int node, leftCtr = 0, rightCtr = 0;

		for (int j = 0; j < dimElem; j++) {
			node = elemToNode[(i+offset)*dimElem+j];
			if (nodePart[node] <= separator) leftCtr++;
			else                             rightCtr++;
		}
		if (leftCtr == dimElem) {
			elemPart[i] = 0;
			(*nbLeftElem)++;
		}
		else if (rightCtr == dimElem) {
			elemPart[i] = 1;
		}
		else {
			elemPart[i] = 2;
			(*nbSepElem)++;
		}
	}
}

// Compute the interval of nodes and elements at each node of the D&C tree
void recursive_tree_creation (tree_t &tree, int *elemToNode, int *sepToNode,
                              int *nodePart, int *nodePartSize, int globalNbElem,
                              int dimElem, int firstPart, int lastPart, int firstElem,
                              int lastElem, int firstNode, int lastNode, int sepOffset
#ifdef STATS
                              , ofstream &dcFile, int curNode, int LRS)
#else
                              )
#endif
{
    // If current node is a leaf, initialize it & exit
    int nbPart = lastPart - firstPart + 1;
    int localNbElem = lastElem - firstElem + 1;
    if (nbPart < 2 || localNbElem <= MAX_ELEM_PER_PART) {
        init_dc_tree (tree, firstElem, lastElem, 0, firstNode, lastNode, true);
        #ifdef STATS
            fill_dc_file_leaves (dcFile, curNode, firstElem, lastElem, LRS);
        #endif
        return;
	}

	// Else, prepare next left, right & separator recursion
	int nbLeftElem = 0, nbSepElem = 0, nbLeftNodes = 0;
	int separator = firstPart + (lastPart - firstPart) / 2;

	// Create local element partition & count left & separator elements
	int *localElemPart = new int [localNbElem];
	if (sepToNode == nullptr) {
        create_elem_part (localElemPart, nodePart, elemToNode, localNbElem, dimElem,
                          separator, firstElem, &nbLeftElem, &nbSepElem);
	}
	else {
        create_elem_part (localElemPart, nodePart, sepToNode, localNbElem, dimElem,
                          separator, sepOffset, &nbLeftElem, &nbSepElem);
	}

    // Count the left nodes if this not a separator
    if (nodePartSize != nullptr) {
        for (int i = firstPart; i <= separator; i++) {
            nbLeftNodes += nodePartSize[i];
        }
    }

	// Create local element permutation
	int *localElemPerm = new int [localNbElem];
	DC_create_permutation (localElemPerm, localElemPart, localNbElem, 3);
	delete[] localElemPart;

    // Execution is correct without mutex although cilkscreen detects a race condition
    pthread_mutex_lock (&mergeMutex);
	// Apply local element permutation to global element permutation
	merge_permutations (localElemPerm, globalNbElem, localNbElem, firstElem,
						lastElem);
    pthread_mutex_unlock (&mergeMutex);

	// Permute elemToNode and sepToNode with local element permutation
	DC_permute_int_2d_array (elemToNode, localElemPerm, localNbElem, dimElem,
                             firstElem);
	if (sepToNode != nullptr) {
		DC_permute_int_2d_array (sepToNode, localElemPerm, localNbElem, dimElem,
                                 sepOffset);
	}
	delete[] localElemPerm;

	// Initialize current node
	init_dc_tree (tree, firstElem, lastElem, nbSepElem, firstNode, lastNode, false);
    #ifdef STATS
        fill_dc_file_nodes (dcFile, curNode, firstElem, lastElem, nbSepElem);
    #endif

	// Left & right recursion
	cilk_spawn
	recursive_tree_creation (*tree.left, elemToNode, sepToNode, nodePart, nodePartSize,
                             globalNbElem, dimElem, firstPart, separator, firstElem,
                             firstElem+nbLeftElem-1, firstNode,
                             firstNode+nbLeftNodes-1, sepOffset
    #ifdef STATS
    					     , dcFile, 3*curNode+1, 1);
    #else
                             );
    #endif
	recursive_tree_creation (*tree.right, elemToNode, sepToNode, nodePart,
                             nodePartSize, globalNbElem, dimElem, separator+1,
                             lastPart, firstElem+nbLeftElem, lastElem-nbSepElem,
                             firstNode+nbLeftNodes, lastNode, sepOffset+nbLeftElem
    #ifdef STATS
                             , dcFile, 3*curNode+2, 2);
    #else
                             );
    #endif

	// D&C partitioning of separator elements
	cilk_sync;
	if (nbSepElem > 0) {
		sep_partitioning (*tree.sep, elemToNode, globalNbElem, dimElem,
        #ifdef STATS
						  lastElem-nbSepElem+1, lastElem, dcFile, 3*curNode+3);
        #else
						  lastElem-nbSepElem+1, lastElem);
        #endif
	}
}

// Create the D&C tree and the permutations
void DC_create_tree (int *elemToNode, int nbElem, int dimElem, int nbNodes)
{
    // Allocate the D&C tree & the permutation functions
    treeHead = new tree_t;
    elemPerm = new int [nbElem];
    nodePerm = new int [nbNodes];

    // Create the D&C tree & the permutation functions
    partitioning (elemToNode, nbElem, dimElem, nbNodes);

    // Hybrid version with coloring of the leaves of the D&C tree
    #ifdef HYBRID
    	coloring (elemToNode, nbElem, dimElem, nbNodes);
    #endif
}

#endif
