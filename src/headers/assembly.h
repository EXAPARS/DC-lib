#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include "DC.h"

// Follow the D&C tree to execute the assembly step in parallel
void recursive_assembly (void (*userSeqFct) (void *, int, int),
                         void (*userVecFct) (void *, int, int), void *userArgs,
                         double *nodeToNodeValue, int operatorDim, tree_t &tree);

// Wrapper used to get the root of the D&C tree before calling the real assembly
// function
void DC_assembly (void (*userSeqFct) (void *, int, int),
                  void (*userVecFct) (void *, int, int),
                  void *userArgs, double *nodeToNodeValue, int operatorDim);

#endif
