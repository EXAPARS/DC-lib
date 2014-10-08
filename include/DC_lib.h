#ifndef DC_LIB_H
#define DC_LIB_H

#include "globals.h"

typedef struct {
    int firstElem, lastElem;
} DCArgs_t;

// Follow the D&C tree to execute the assembly step in parallel
void assembly (tree_t &tree, void (*userSeqFct) (DCArgs_t *, void *),
               void (*userVecFct) (DCArgs_t *, void *), void *userArgs,
               double *nodeToNodeValue, int operatorDim);

// Wrapper used to get the root of the D&C tree before calling the real assembly
// function
void DC_assembly (void (*userSeqFct) (DCArgs_t *, void *),
                  void (*userVecFct) (DCArgs_t *, void *),
                  void *userArgs, double *nodeToNodeValue, int operatorDim);

// Compute edge intervals for CSR reset and store the D&C tree and the permutations
void DC_finalize (int *nodeToNodeRow, int *elemToNode, int nbElem, int nbNodes,
                  int nbBlocks, int mpiRank);

// Apply the permutations
void DC_permutation (double *coord, int *elemToNode, int *intfNodes, int *dispList,
                     int *boundNodesCode, int nbElem, int nbNodes, int nbIntfNodes,
                     int nbDispNodes, int mpiRank);

// Create or load from file the D&C tree and the permutations
void DC_start (int *elemToNode, int nbElem, int nbNodes, int nbBlocks, int mpiRank);

#endif
