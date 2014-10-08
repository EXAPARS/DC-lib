#include <mpi.h>
#include <iostream>
#include <cilk/cilk.h>

#include "permutations.h"
#include "IO.h"
#include "coloring.h"
#include "partitioning.h"
#include "DC_lib.h"
 
// The D&C tree is a global variable in order to persist from one call to the library
// to another
tree_t *treeHead = NULL;

// Follow the D&C tree to execute the assembly step in parallel
void assembly (tree_t &tree, void (*userSeqFct) (DCArgs_t *, void *),
               void (*userVecFct) (DCArgs_t *, void *), void *userArgs,
               double *nodeToNodeValue, int operatorDim)
{
    // If current node is a leaf, call the appropriate assembly function
    if (tree.left == NULL && tree.right == NULL) {
        DCArgs_t DCArgs;

        // If leaf is not a separator, reset locally the CSR matrix
        if (tree.firstCSR != -1) {
            int firstCSR = tree.firstCSR * operatorDim;
            int lastCSR  = (tree.lastCSR + 1) * operatorDim - firstCSR;
            nodeToNodeValue[firstCSR:lastCSR] = 0;
        }

#ifdef HYBRID
        // Call user vectorial function on full colors
        DCArgs = {tree.firstElem, tree.vecOffset};
        userVecFct (&DCArgs, userArgs);

        // Call user sequential function on other colors
        DCArgs = {tree.vecOffset+1, tree.lastElem};
        userSeqFct (&DCArgs, userArgs);
#else
        // Call user sequential function
        DCArgs = {tree.firstElem, tree.lastElem};
        userSeqFct (&DCArgs, userArgs);
#endif
    }
    else {
        // Left & right recursion
        cilk_spawn
        assembly (*tree.left, userSeqFct, userVecFct, userArgs, nodeToNodeValue,
                  operatorDim);
        assembly (*tree.right, userSeqFct, userVecFct, userArgs, nodeToNodeValue,
                  operatorDim);

        // Synchronization
        cilk_sync;

        // Separator recursion, if it is not empty
        if (tree.sep != NULL) {
            assembly (*tree.sep, userSeqFct, userVecFct, userArgs, nodeToNodeValue,
                      operatorDim);
        }
    }
}

// Wrapper used to get the root of the D&C tree before calling the real assembly
// function
void DC_assembly (void (*userSeqFct) (DCArgs_t *, void *),
                  void (*userVecFct) (DCArgs_t *, void *),
                  void *userArgs, double *nodeToNodeValue, int operatorDim)
{
    assembly (*treeHead, userSeqFct, userVecFct, userArgs, nodeToNodeValue,
              operatorDim);
}

// Compute edge intervals for CSR reset and store the D&C tree and the permutations
void DC_finalize (int *nodeToNodeRow, int *elemToNode, int nbElem, int nbNodes,
                  int nbBlocks, int mpiRank)
{
    #ifdef CREATE_PERM_AND_TREE
        double t1, t2;
        // Compute the interval of edges of each leaf of the D&C tree
        if (mpiRank == 0) {
            cout << "Computing edge intervals...          ";
    	    t1 = MPI_Wtime ();
        }
        edge_intervals (*treeHead, nodeToNodeRow, elemToNode, nbNodes);
        if (mpiRank == 0) {
            t2 = MPI_Wtime ();
            cout << "done  (" << t2 - t1 << " seconds)\n";
        }

        // Store the D&C tree & the permutation functions
    	if (mpiRank == 0) {
            cout << "Storing D&C tree and permutations... ";
            t1 = MPI_Wtime ();
        }
    	store_perm_and_tree (nbElem, nbNodes, nbBlocks, mpiRank);
        if (mpiRank == 0) {
            t2 = MPI_Wtime ();
        	cout << "done  (" << t2 - t1 << " seconds)\n";
        }
        delete[] nodePerm, delete[] elemPerm; 
    #endif
}

// Apply the permutations
void DC_permutation (double *coord, int *elemToNode, int *intfNodes, int *dispList,
                     int *boundNodesCode, int nbElem, int nbNodes, int nbIntfNodes,
                     int nbDispNodes, int mpiRank)
{
    double t1, t2;
    if (mpiRank == 0) {
        cout << "Applying permutations...             ";
        t1 = MPI_Wtime ();
    }
    apply_transformations (coord, elemToNode, intfNodes, dispList, boundNodesCode,
                           nbElem, nbNodes, nbIntfNodes, nbDispNodes);
    if (mpiRank == 0) {
        t2 = MPI_Wtime ();
    	cout << "done  (" << t2 - t1 << " seconds)\n";
    }
}

// Create or load from file the D&C tree and the permutations
void DC_start (int *elemToNode, int nbElem, int nbNodes, int nbBlocks, int mpiRank)
{
    double t1, t2;

    // Allocate the D&C tree & the permutation functions
    treeHead = new tree_t;
    elemPerm = new int [nbElem];
    nodePerm = new int [nbNodes];

    // Create the D&C tree & the permutation functions
    #ifdef CREATE_PERM_AND_TREE
        if (mpiRank == 0) {
            cout << "D&C partitioning...                  ";
            t1 = MPI_Wtime ();
        }
        partitioning (elemToNode, nbElem, nbNodes);
        if (mpiRank == 0) {
            t2 = MPI_Wtime ();
        	cout << "done  (" << t2 - t1 << " seconds)\n";
        }

        // Hybrid version with coloring of the leaves of the D&C tree
        #ifdef HYBRID
        	if (mpiRank == 0) {
            	cout << "Coloring D&C subdomains...           ";
                t1 = MPI_Wtime ();
            }
        	coloring (elemToNode, nbElem, nbNodes);
        	if (mpiRank == 0) {
                t2 = MPI_Wtime ();
            	cout << "done  (" << t2-t1 << " seconds)\n";
            }
        #endif

    // Read the D&C tree & the permutation functions
    #else
    	if (mpiRank == 0) {
            cout << "Reading D&C tree and permutations... ";
            t1 = MPI_Wtime ();
        }
    	read_perm_and_tree (nbElem, nbNodes, nbBlocks, mpiRank);
        if (mpiRank == 0) {
            t2 = MPI_Wtime ();
        	cout << "done  (" << t2 - t1 << " seconds)\n";
        }
    #endif
}
