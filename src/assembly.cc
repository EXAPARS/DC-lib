#include <cilk/cilk.h>
#include <stdlib.h>

#include "assembly.h"

extern tree_t *treeHead;

// Follow the D&C tree to execute the assembly step in parallel
void recursive_assembly (void (*userSeqFct) (void *, int, int),
                         void (*userVecFct) (void *, int, int), void *userArgs,
                         double *nodeToNodeValue, int operatorDim, tree_t &tree)
{
    // If current node is a leaf, call the appropriate assembly function
    if (tree.left == nullptr && tree.right == nullptr) {

        // If leaf is not a separator, reset locally the CSR matrix
        if (tree.firstCSR != -1) {
            int firstCSR = tree.firstCSR * operatorDim;
            int lastCSR  = (tree.lastCSR + 1) * operatorDim - firstCSR;
            nodeToNodeValue[firstCSR:lastCSR] = 0;
        }

        #ifdef HYBRID
            // Call user vectorial function on full colors
            userVecFct (userArgs, tree.firstElem, tree.vecOffset);

            // Call user sequential function on other colors
            userSeqFct (userArgs, tree.vecOffset+1, tree.lastElem);
        #else
            // Call user sequential function
            userSeqFct (userArgs, tree.firstElem, tree.lastElem);
        #endif
    }
    else {
        // Left & right recursion
        cilk_spawn
        recursive_assembly (userSeqFct, userVecFct, userArgs, nodeToNodeValue,
                            operatorDim, *tree.left);
        recursive_assembly (userSeqFct, userVecFct, userArgs, nodeToNodeValue,
                            operatorDim, *tree.right);

        // Synchronization
        cilk_sync;

        // Separator recursion, if it is not empty
        if (tree.sep != nullptr) {
            recursive_assembly (userSeqFct, userVecFct, userArgs, nodeToNodeValue,
                                operatorDim, *tree.sep);
        }
    }
}

// Wrapper used to get the root of the D&C tree before calling the real assembly
// function
void DC_assembly (void (*userSeqFct) (void *, int, int),
                  void (*userVecFct) (void *, int, int),
                  void *userArgs, double *nodeToNodeValue, int operatorDim)
{
    recursive_assembly (userSeqFct, userVecFct, userArgs, nodeToNodeValue, operatorDim,
                        *treeHead);
}
