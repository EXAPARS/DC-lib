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
#include <stdlib.h>

#include "tree_traversal.h"

extern tree_t *treeHead;

// Follow the D&C tree to execute the given function in parallel
void tree_traversal (void (*userSeqFct) (void *, int, int),
                     void (*userVecFct) (void *, int, int), void *userArgs,
                     double *nodeToNodeValue, int operatorDim, tree_t &tree)
{
    // If current node is a leaf, call the appropriate function
    if (tree.left == nullptr && tree.right == nullptr) {

        // If leaf is not a separator, reset locally the CSR matrix
        if (tree.firstEdge != -1) {
            int firstEdge = tree.firstEdge * operatorDim;
            int lastEdge  = (tree.lastEdge + 1) * operatorDim - firstEdge;
            nodeToNodeValue[firstEdge:lastEdge] = 0;
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
        tree_traversal (userSeqFct, userVecFct, userArgs, nodeToNodeValue,
                        operatorDim, *tree.right);
        tree_traversal (userSeqFct, userVecFct, userArgs, nodeToNodeValue,
                        operatorDim, *tree.left);

        // Synchronization
        cilk_sync;

        // Separator recursion, if it is not empty
        if (tree.sep != nullptr) {
            tree_traversal (userSeqFct, userVecFct, userArgs, nodeToNodeValue,
                            operatorDim, *tree.sep);
        }
    }
}

// Wrapper used to get the root of the D&C tree before calling the real tree traversal
void DC_tree_traversal (void (*userSeqFct) (void *, int, int),
                        void (*userVecFct) (void *, int, int),
                        void *userArgs, double *nodeToNodeValue, int operatorDim)
{
    tree_traversal (userSeqFct, userVecFct, userArgs, nodeToNodeValue, operatorDim,
                    *treeHead);
}
