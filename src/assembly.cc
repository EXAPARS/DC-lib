#include <cilk/cilk.h>
#include "assembly.h"

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
