#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include "DC_lib.h"
#include "globals.h"

// Follow the D&C tree to execute the assembly step in parallel
void assembly (void (*userSeqFct) (DCArgs_t *, void *),
               void (*userVecFct) (DCArgs_t *, void *), void *userArgs,
               double *nodeToNodeValue, int operatorDim, tree_t &tree);

// Wrapper used to get the root of the D&C tree before calling the real assembly
// function
void DC_assembly (void (*userSeqFct) (DCArgs_t *, void *),
                  void (*userVecFct) (DCArgs_t *, void *),
                  void *userArgs, double *nodeToNodeValue, int operatorDim);

#endif
