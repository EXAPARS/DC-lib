/*  Copyright 2014 - UVSQ
    Authors list: Loïc Thébault, Eric Petit

    This file is part of the DC-lib.

    DC-lib is free software: you can redistribute it and/or modify it under the
    terms of the GNU Lesser General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later version.

    DC-lib is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
    PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License along with
    the DC-lib. If not, see <http://www.gnu.org/licenses/>. */

#ifdef CILK
    #include <cilk/cilk.h>
#endif
#include <pthread.h>
#include <limits.h>

#include "tools.h"
#include "coloring.h"
#include "permutations.h"
#include "partitioning.h"
#include "tree_creation.h"

// The D&C tree and the permutations are global in order to persist from one call
// to the library to another
tree_t *treeHead = nullptr;
int *elemPerm = nullptr, *nodePerm = nullptr, *nodeOwner = nullptr;

#ifdef TREE_CREATION

// Mutex to avoid race condition in merge permutations
pthread_mutex_t mergeMutex = PTHREAD_MUTEX_INITIALIZER;

/*
// Split the D&C tree into two D&C trees. One for the interfaces, the other for the
// inner domain
bool intf_tree_creation ()
{
    // Look if current D&C node contains at least one node on one interface
    bool hasIntfNode = false;
    for (int i = firstElem * dimElem; i < (lastElem+1) * dimElem; i++) {
        for (int j = 0; j < nbIntf; j++) {
            for (int k = intfIndex[j]; k < intfIndex[j+1]; k++) {
                int intfNode = intfNodes[k];
                if (elemToNode[i] == intfNode) {
                    hasIntfNode = true;
                    break;
                }
            }
            if (hasIntfNode) break;
        }
        if (hasIntfNode) break;
    }

    // If current node is a leaf
    if (tree.left == nullptr && tree.right == nullptr) {
        // If leaf contains interface node(s)
        if (hasIntfNode) {
            // new feuille in DC intf tree containing only current leaf
        }
    }
    else {
        bool leftHasIntfNode, rightHasIntfNode, sepHasIntfNode = false;

        // Left & right recursion
        cilk_spawn
        rightHasIntfNode = intf_tree_creation (*tree.right);
        leftHasIntfNode  = intf_tree_creation (*tree.left);

        // Synchronization
        cilk_sync;

        // Separator recursion, if it is not empty
        if (tree.sep != nullptr) {
            sepHasIntfNode = intf_tree_creation (*tree.sep);
        }

        // If both left & right sons have interface node(s)
        if (leftHasIntfNode && rightHasIntfNode) {
            // new feuille in DC intf tree containing only current node
        }
        // If exclusively left or right son has interface node(s)
        if ((leftHasIntfNode && !rightHasIntfNode) ||
           (!leftHasIntfNode &&  rightHasIntfNode)) {
            // new feuille in DC intf tree containing current node & the left or right
        }
    }

    return hasIntfNode;
}
*/

// Compute the edge interval and the list of nodes owned by each leaf of the D&C tree
void compute_intervals (tree_t &tree, int *nodeToNodeRow, int *elemToNode, int curNode)
{
    // If current node is a leaf
    if (tree.left == nullptr && tree.right == nullptr) {

        // Compute the number of nodes owned by current leaf and fill the node list
        int nbOwnedNodes = 0, ctr = 0;
        for (int i = tree.firstNode; i <= tree.lastNode; i++) {
            if (nodeOwner[i] == curNode) nbOwnedNodes++;
        }
        tree.nbOwnedNodes = nbOwnedNodes;
        if (nbOwnedNodes > 0) tree.ownedNodes = new int [nbOwnedNodes];
        for (int i = tree.firstNode; i <= tree.lastNode; i++) {
            if (nodeOwner[i] == curNode) {
                tree.ownedNodes[ctr] = i;
                ctr++;
            }
            if (ctr == nbOwnedNodes) break;
        }

        // Get the first and last edges of the leaf if current leaf is not a separator
        if (tree.isSep) {
            tree.firstNode = -1;
            tree.lastNode  = -1;
            tree.firstEdge = -1;
            tree.lastEdge  = -1;
        }
        else {
            tree.firstEdge = nodeToNodeRow[tree.firstNode];
            tree.lastEdge  = nodeToNodeRow[tree.lastNode+1] - 1;
        }
    }
    else {
        #ifdef OMP
            // Left & right recursion
            #pragma omp task default (shared)
            compute_intervals (*tree.right, nodeToNodeRow, elemToNode, 3*curNode+2);
            #pragma omp task default (shared)
            compute_intervals (*tree.left, nodeToNodeRow, elemToNode, 3*curNode+1);
            // Synchronization
            #pragma omp taskwait
        #elif CILK
            // Left & right recursion
            cilk_spawn
            compute_intervals (*tree.right, nodeToNodeRow, elemToNode, 3*curNode+2);
            compute_intervals (*tree.left, nodeToNodeRow, elemToNode, 3*curNode+1);
            // Synchronization
            cilk_sync;
        #endif
        // Separator recursion, if it is not empty
        if (tree.sep != nullptr) {
            compute_intervals (*tree.sep, nodeToNodeRow, elemToNode, 3*curNode+3);
        }
    }
}

// Wrapper used to get the root of the D&C tree before computing the edge intervals and
// the list of nodes owned by each leaf of the D&C tree
void DC_finalize_tree (int *nodeToNodeRow, int *elemToNode)
{
    #ifdef OMP
        #pragma omp parallel
        #pragma omp single nowait
    #endif
    compute_intervals (*treeHead, nodeToNodeRow, elemToNode, 0);
    delete[] nodeOwner;
}

// Initialize a node of the D&C tree
void init_dc_tree (tree_t &tree, int *elemToNode, int *intfIndex, int *intfNodes,
                   int firstElem, int lastElem, int nbSepElem, int dimElem,
                   int firstNode, int lastNode, int nbIntf, int commDepth,
                   int curDepth, bool isSep, bool isLeaf, bool *hasIntfNode)
{
    tree.ownedNodes   = nullptr;
    tree.nbOwnedNodes = -1;
    tree.firstElem    = firstElem;
    tree.lastElem     = lastElem - nbSepElem;
    tree.lastSep      = lastElem;
    tree.firstNode    = firstNode;
    tree.lastNode     = lastNode;
    tree.firstEdge    = -1;
    tree.lastEdge     = -1;
    #ifdef DC_VEC
        tree.vecOffset = 0;
    #endif
    tree.isSep = isSep;
    tree.left  = nullptr;
    tree.right = nullptr;
    tree.sep   = nullptr;
    if (isLeaf == false) {
        tree.left  = new tree_t;
        tree.right = new tree_t;
        if (nbSepElem > 0) {
            tree.sep = new tree_t;
        }
    }

    #ifdef MULTI_THREADED_COMM
        tree.intfIndex = nullptr;
        tree.intfNodes = nullptr;

        // If current node is a leaf upper or at the communication level or an internal
        // node at communication level, initialize its interface
        if ((isLeaf && curDepth <= commDepth) || (!isLeaf && curDepth == commDepth)) {
            int nbIntfNodes = 0, ctr = 0;
            // Count the number of nodes on the interface
            for (int i = firstElem * dimElem; i < (lastElem+1) * dimElem; i++) {
                for (int j = 0; j < nbIntf; j++) {
                    for (int k = intfIndex[j]; k < intfIndex[j+1]; k++) {
                        int intfNode = intfNodes[k] - 1;
                        if (elemToNode[i] == intfNode) {
                            *hasIntfNode = true;
                            nbIntfNodes++;
                        }
                    }
                }
            }
            // If there is at least one node on the interface
            if (nbIntfNodes > 0) {
                // Allocate the interface index
                tree.intfIndex = new int [nbIntf+1];
                tree.intfNodes = new int [nbIntfNodes];

                // Initialize the interface index
                for (int j = 0; j < nbIntf; j++) {
                    tree.intfIndex[j] = ctr;
                    for (int k = intfIndex[j]; k < intfIndex[j+1]; k++) {
                        int intfNode = intfNodes[k] - 1;
                        for (int i = firstElem*dimElem; i < (lastElem+1)*dimElem; i++){
                            if (elemToNode[i] == intfNode) {
                                tree.intfNodes[ctr] = intfNode + 1;
                                ctr++;
                            }
                        }
                    }
                }
                tree.intfIndex[nbIntf] = ctr;
            }
        }
    #endif
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

// Create the D&C tree and the element permutation, and compute the intervals of nodes
// and elements at each node of the tree
void tree_creation (tree_t &tree, int *elemToNode, int *sepToNode, int *nodePart,
                    int *nodePartSize, int *intfIndex, int *intfNodes,
                    int globalNbElem, int dimElem, int firstPart, int lastPart,
                    int firstElem, int lastElem, int firstNode, int lastNode,
                    int sepOffset, int nbIntf, int curNode, int commDepth,
#ifdef STATS
                    int curDepth, bool isSep, ofstream &dcFile, int LRS)
#else
                    int curDepth, bool isSep)
#endif
{
    // If current node is a leaf
    int nbPart = lastPart - firstPart + 1;
    int localNbElem = lastElem - firstElem + 1;
    bool hasIntfNode = false;
    if (nbPart < 2 || localNbElem <= MAX_ELEM_PER_PART) {

        // Initialize the leaf
        init_dc_tree (tree, elemToNode, intfIndex, intfNodes, firstElem, lastElem, 0,
                      dimElem, firstNode, lastNode, nbIntf, commDepth, curDepth, isSep,
                      true, &hasIntfNode);

        // Set the last updater of each node
        if (isSep) {
            for (int i = firstElem * dimElem; i < (lastElem+1) * dimElem; i++) {
                int node = nodePerm[elemToNode[i]];
                nodeOwner[node] = curNode;
            }
        }
        else {
            for (int i = firstNode; i <= lastNode; i++) {
                nodeOwner[i] = curNode;
            }
        }

        #ifdef STATS
            fill_dc_file_leaves (dcFile, curNode, firstElem, lastElem, LRS,
                                 hasIntfNode);
            count_intf_stats (hasIntfNode);
        #endif

        // End of recursion
        return;
    }

    // Else, prepare next left, right & separator recursion
    int nbLeftElem = 0, nbSepElem = 0, nbLeftNodes = 0, nbRightNodes = 0;
    int separator = firstPart + (lastPart - firstPart) / 2;

    // Create local element partition & count left & separator elements
    int *localElemPart = new int [localNbElem];
    if (isSep) {
        create_elem_part (localElemPart, nodePart, sepToNode, localNbElem, dimElem,
                          separator, sepOffset, &nbLeftElem, &nbSepElem);
    }
    else {
        create_elem_part (localElemPart, nodePart, elemToNode, localNbElem, dimElem,
                          separator, firstElem, &nbLeftElem, &nbSepElem);
    }

    // Count the left & right nodes if current D&C node is not a separator
    if (isSep) {
        nbLeftNodes  = lastNode - firstNode + 1;
        nbRightNodes = lastNode - firstNode + 1;
    }
    else {
        for (int i = firstPart; i <= separator; i++) {
            nbLeftNodes += nodePartSize[i];
        }
        nbRightNodes = (lastNode - firstNode + 1) - nbLeftNodes;
    }

    // Create local element permutation
    int *localElemPerm = new int [localNbElem];
    DC_create_permutation (localElemPerm, localElemPart, localNbElem, 3);
    delete[] localElemPart;

    // Execution is correct without mutex although cilkscreen detects a race condition
    pthread_mutex_lock (&mergeMutex);
    // Apply local element permutation to global element permutation
    merge_permutations (localElemPerm, globalNbElem, localNbElem, firstElem, lastElem);
    pthread_mutex_unlock (&mergeMutex);

    // Permute elemToNode and sepToNode with local element permutation
    DC_permute_int_2d_array (elemToNode, localElemPerm, localNbElem, dimElem,
                             firstElem);
    if (isSep) {
        DC_permute_int_2d_array (sepToNode, localElemPerm, localNbElem, dimElem,
                                 sepOffset);
    }
    delete[] localElemPerm;

    // Initialize current node
    init_dc_tree (tree, elemToNode, intfIndex, intfNodes, firstElem, lastElem,
                  nbSepElem, dimElem, firstNode, lastNode, nbIntf, commDepth, curDepth,
                  isSep, false, &hasIntfNode);
    #ifdef STATS
        fill_dc_file_nodes (dcFile, curNode, firstElem, lastElem, nbSepElem,
                            hasIntfNode);
    #endif

    // Left & right recursion
    #ifdef OMP
        #pragma omp task default(shared)
    #elif CILK
        cilk_spawn
    #endif
    tree_creation (*tree.right, elemToNode, sepToNode, nodePart, nodePartSize,
                   intfIndex, intfNodes, globalNbElem, dimElem, separator+1, lastPart,
                   firstElem+nbLeftElem, lastElem-nbSepElem, lastNode-nbRightNodes+1,
                   lastNode, sepOffset+nbLeftElem, nbIntf, 3*curNode+2, commDepth,
    #ifdef STATS
                   curDepth+1, isSep, dcFile, 2);
    #else
                   curDepth+1, isSep);
    #endif
    #ifdef OMP
        #pragma omp task default(shared)
    #endif
    tree_creation (*tree.left, elemToNode, sepToNode, nodePart, nodePartSize,
                   intfIndex, intfNodes, globalNbElem, dimElem, firstPart, separator,
                   firstElem, firstElem+nbLeftElem-1, firstNode, firstNode+nbLeftNodes
                   -1, sepOffset, nbIntf, 3*curNode+1, commDepth, curDepth+1, isSep
    #ifdef STATS
                   , dcFile, 1);
    #else
                   );
    #endif

    // Synchronization
    #ifdef OMP
        #pragma omp taskwait
    #elif CILK
        cilk_sync;
    #endif

    // D&C partitioning of separator elements
    if (nbSepElem > 0) {
        sep_partitioning (*tree.sep, elemToNode, intfIndex, intfNodes, globalNbElem,
                          dimElem, lastElem-nbSepElem+1, lastElem, firstNode, lastNode,
        #ifdef STATS
                          nbIntf, 3*curNode+3, commDepth, curDepth+1, dcFile);
        #else
                          nbIntf, 3*curNode+3, commDepth, curDepth+1);
        #endif
    }
}

// Create the D&C tree and the permutations
void DC_create_tree (double *coord, int *elemToNode, int *intfIndex, int *intfNodes,
                     int nbElem, int dimElem, int nbNodes, int dimNode, int nbIntf,
                     int rank)
{
    // Allocate the D&C tree & the permutation functions
    treeHead = new tree_t;
    elemPerm = new int [nbElem];
    nodePerm = new int [nbNodes];

    // Initialize the owner of the nodes to MAX_INT (lower is owner)
    nodeOwner = new int [nbNodes];
    #ifdef OMP
        #pragma omp parallel for
        for (int i = 0; i < nbNodes; i++) {
    #elif CILK
        cilk_for (int i = 0; i < nbNodes; i++) {
    #endif
        nodeOwner[i] = INT_MAX;
    }

    // Create the D&C tree & the permutation functions
    partitioning (elemToNode, intfIndex, intfNodes, nbElem, dimElem, nbNodes, nbIntf,
                  rank);

    // Vectorial version with coloring of the leaves of the D&C tree
    #ifdef DC_VEC
        coloring (elemToNode, nbElem, dimElem, nbNodes);
    #endif
}

#endif
