#ifndef DC_LIB_H
#define DC_LIB_H

#include <stdint.h>

typedef struct {
    int firstElem, lastElem;
} DCArgs_t;

// Get CPU cycles
uint64_t DC_get_cycles ();

// Wrapper used to get the root of the D&C tree before calling the real assembly
// function
void DC_assembly (void (*userSeqFct) (DCArgs_t *, void *),
                  void (*userVecFct) (DCArgs_t *, void *),
                  void *userArgs, double *nodeToNodeValue, int operatorDim);

// Apply the permutations
void DC_permutation (double *coord, int *elemToNode, int *intfNodes, int *dispList,
                     int *boundNodesCode, int nbElem, int dimElem, int nbNodes,
                     int dimNode, int nbIntfNodes, int nbDispNodes);

// Read the D&C tree and the permutation functions
void DC_read_tree (int nbElem, int nbNodes, int nbBlocks, int mpiRank);

// Store the D&C tree and the permutation functions to a binary file
void DC_store_tree (int nbElem, int nbNodes, int nbBlocks, int mpiRank);

// Wrapper used to get the root of the D&C tree before computing the edge intervals
// for CSR reset
void DC_finalize_tree (int *nodeToNodeRow, int *elemToNode, int nbNodes);

// Create the D&C tree and the permutations
void DC_create_tree (int *elemToNode, int nbElem, int dimElem, int nbNodes);

#endif
