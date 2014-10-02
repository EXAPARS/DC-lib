#ifdef CREATE_PERM_AND_TREE

#include <iostream>
#include <cilk/cilk.h>

#include "tools.h"
#include "permutations.h"
#include "partitioning.h"
#include "tree.h"

#include <pthread.h>
pthread_mutex_t mergeMutex = PTHREAD_MUTEX_INITIALIZER;

// Compare two different D&C trees and display the differences
void compare_dc_trees (tree_t &tree1, tree_t &tree2, int curNode)
{
	// If trees have not the same topology
	if ((tree1.left == NULL && tree1.right == NULL) &&
		(tree2.left != NULL || tree2.right != NULL)) {
		cout << "Node " << curNode << " -> Tree1 is a leaf, Tree2 is not.\n";
	}
	else if ((tree2.left == NULL && tree2.right == NULL)
		  && (tree1.left != NULL || tree1.right != NULL)) {
		cout << "Node " << curNode << " -> Tree2 is a leaf, Tree1 is not.\n";
	}
	if (tree1.sep == NULL && tree2.sep != NULL) {
		cout << "Node " << curNode << " -> Tree1 has a sep, Tree2 has not.\n";
    }
	else if (tree1.sep != NULL && tree2.sep == NULL) {
		cout << "Node " << curNode << " -> Tree2 has a sep, Tree1 has not.\n";
    }

	// If trees have not the same intervals
	if (tree1.firstElem != tree2.firstElem ||
		tree1.lastElem  != tree2.lastElem  || tree1.lastSep != tree2.lastSep) {
		cout << "Node " << curNode << " -> "
			 << "Tree 1: [" << tree1.firstElem << ", " << tree1.lastElem
			 << ", " << tree1.lastSep << "], "
			 << "Tree 2: [" << tree2.firstElem << ", " << tree2.lastElem
			 << ", " << tree2.lastSep << "]\n";
	}

	// Left, right & separator recursion
	if ((tree1.left != NULL && tree1.right != NULL) &&
		(tree2.left != NULL && tree2.right != NULL)) {
		compare_dc_trees (*tree1.left,  *tree2.left,  3*curNode+1);
		compare_dc_trees (*tree1.right, *tree2.right, 3*curNode+2);
	}
	if (tree1.sep != NULL && tree2.sep != NULL) {
		compare_dc_trees (*tree1.sep,   *tree2.sep,   3*curNode+3);
    }
}

// Compute the interval of edges of each leaf of the D&C tree
void edge_intervals (tree_t &tree, int *nodeToNodeRow, int *elemToNode, int nbNodes)
{
     // If current node is a leaf
    if (tree.left == NULL && tree.right == NULL) {
        // Get the first and last edges of the leaf
        int firstNode = tree.firstCSR, lastNode = tree.lastCSR;
        tree.firstCSR = nodeToNodeRow[firstNode];
        tree.lastCSR  = nodeToNodeRow[lastNode+1] - 1;
    }
    else {
        // Left & right recursion
        cilk_spawn
        edge_intervals (*tree.left, nodeToNodeRow, elemToNode, nbNodes);
        edge_intervals (*tree.right, nodeToNodeRow, elemToNode, nbNodes);

        // Synchronization
        cilk_sync;
    }
}

// Fill the element intervals of the D&C tree
void fill_dc_tree (tree_t &tree, int firstElem, int lastElem, int nbSepElem,
				   int firstNode, int lastNode, bool isLeaf)
{
	tree.firstElem = firstElem;
	tree.lastElem  = lastElem - nbSepElem;
	tree.lastSep   = lastElem;
    tree.firstCSR  = firstNode;
    tree.lastCSR   = lastNode;
    tree.lastFullColor = 0;
	tree.left   = NULL;
	tree.right  = NULL;
	tree.sep    = NULL;

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
                       int separator, int offset, int *nbLeftElem, int *nbSepElem)
{
	for (int i = 0; i < nbElem; i++) {
		int node, leftCtr = 0, rightCtr = 0;

		for (int j = 0; j < DIM_ELEM; j++) {
			node = elemToNode[(i+offset)*DIM_ELEM+j];
			if (nodePart[node] <= separator) leftCtr++;
			else                             rightCtr++;
		}
		if (leftCtr == DIM_ELEM) {
			elemPart[i] = 0;
			(*nbLeftElem)++;
		}
		else if (rightCtr == DIM_ELEM) {
			elemPart[i] = 1;
		}
		else {
			elemPart[i] = 2;
			(*nbSepElem)++;
		}
	}
}

// Compute the interval of nodes and elements at each node of the D&C tree
void create_dc_tree (tree_t &tree, int *elemToNode, int *sepToNode, int *nodePart,
                     int *nodePartSize, int globalNbElem, int firstPart, int lastPart,
                     int firstElem, int lastElem, int firstNode, int lastNode,
#ifdef STATS
					 int sepOffset, ofstream &dcFile, int curNode, int LRS)
#else
					 int sepOffset)
#endif
{
    // If current node is a leaf, fill D&C tree & exit
    int nbPart = lastPart - firstPart + 1;
    int localNbElem = lastElem - firstElem + 1;
    if (nbPart < 2 || localNbElem <= MAX_ELEM_PER_PART) {
        fill_dc_tree (tree, firstElem, lastElem, 0, firstNode, lastNode, true);
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
	if (sepToNode == NULL) {
        create_elem_part (localElemPart, nodePart, elemToNode, localNbElem, separator,
                          firstElem, &nbLeftElem, &nbSepElem);
	}
	else {
        create_elem_part (localElemPart, nodePart, sepToNode, localNbElem, separator,
                          sepOffset, &nbLeftElem, &nbSepElem);
	}

    // Count the left nodes if this not a separator
    if (nodePartSize != NULL) {
        for (int i = firstPart; i <= separator; i++) {
            nbLeftNodes += nodePartSize[i];
        }
    }

	// Create local element permutation
	int *localElemPerm = new int [localNbElem];
	create_perm_array (localElemPerm, localElemPart, localNbElem, 3);
	delete[] localElemPart;

    // Execution is correct without mutex although cilkscreen detects a race condition
    pthread_mutex_lock (&mergeMutex);
	// Apply local element permutation to global element permutation
	merge_permutations (localElemPerm, globalNbElem, localNbElem, firstElem,
						lastElem);
    pthread_mutex_unlock (&mergeMutex);

	// Permute elemToNode and sepToNode with local element permutation
	permute_int_2d_array (elemToNode, localElemPerm, localNbElem, DIM_ELEM,
                          firstElem);
	if (sepToNode != NULL) {
		permute_int_2d_array (sepToNode, localElemPerm, localNbElem, DIM_ELEM,
                              sepOffset);
	}
	delete[] localElemPerm;

	// Fill D&C tree
	fill_dc_tree (tree, firstElem, lastElem, nbSepElem, firstNode, lastNode, false);
#ifdef STATS
	fill_dc_file_nodes (dcFile, curNode, firstElem, lastElem, nbSepElem);
#endif

	// Left & right recursion
	cilk_spawn
	create_dc_tree (*tree.left, elemToNode, sepToNode, nodePart, nodePartSize,
                    globalNbElem, firstPart, separator, firstElem,
                    firstElem+nbLeftElem-1, firstNode, firstNode+nbLeftNodes-1,
#ifdef STATS
					sepOffset, dcFile, 3*curNode+1, 1);
#else
					sepOffset);
#endif
	create_dc_tree (*tree.right, elemToNode, sepToNode, nodePart, nodePartSize,
                    globalNbElem, separator+1, lastPart, firstElem+nbLeftElem,
                    lastElem-nbSepElem, firstNode+nbLeftNodes, lastNode,
#ifdef STATS
                    sepOffset+nbLeftElem, dcFile, 3*curNode+2, 2);
#else
                    sepOffset+nbLeftElem);
#endif

	// D&C partitioning of separator elements
	cilk_sync;
	if (nbSepElem > 0) {
		sep_partitioning (*tree.sep, elemToNode, globalNbElem, lastElem-nbSepElem+1,
#ifdef STATS
						  lastElem, dcFile, 3*curNode+3);
#else
						  lastElem);
#endif
	}
}

#endif
