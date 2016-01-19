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
#include <stdlib.h>

#include "tree_traversal.h"

extern tree_t *treeHead;

// Follow the D&C tree to execute the given function in parallel
void tree_traversal (void (*userSeqFct)  (void *, DCargs_t *),
                     void (*userVecFct)  (void *, DCargs_t *),
                     void (*userCommFct) (void *, DCcommArgs_t *),
                     void *userArgs, void *userCommArgs, tree_t &tree)
{
    // If current node is a leaf, call the appropriate function
    if (tree.left == nullptr && tree.right == nullptr) {

        // Initialize the D&C arguments
        DCargs_t DCargs;
        DCargs.firstNode = tree.firstNode;
        DCargs.lastNode  = tree.lastNode;
        DCargs.firstEdge = tree.firstEdge;
        DCargs.lastEdge  = tree.lastEdge;
        DCargs.isSep     = tree.isSep;
        #ifdef MULTITHREADED_COMM
            DCargs.nbOwnedNodes = tree.nbOwnedNodes;
            DCargs.ownedNodes   = tree.ownedNodes;
        #endif

        #ifdef DC_VEC
            // Call user vectorial function on full colors
            DCargs.firstElem = tree.firstElem;
            DCargs.lastElem  = tree.vecOffset;
            userVecFct (userArgs, &DCargs);

            // Call user sequential function on other colors
            DCargs.firstElem = tree.vecOffset+1;
            DCargs.lastElem  = tree.lastElem;
            userSeqFct (userArgs, &DCargs);
        #else
            // Call user sequential function
            DCargs.firstElem = tree.firstElem;
            DCargs.lastElem  = tree.lastElem;
            userSeqFct (userArgs, &DCargs);
        #endif

        // Call user communication function
        #ifdef MULTITHREADED_COMM
            if (tree.nbIntfNodes > 0) {
                DCcommArgs_t DCcommArgs;
                DCcommArgs.intfIndex = tree.intfIndex;
                DCcommArgs.intfNodes = tree.intfNodes;
                DCcommArgs.intfDst   = tree.intfDst;
                userCommFct (userCommArgs, &DCcommArgs);
            }
        #endif
    }
    else {
        #ifdef OMP
            // Left & right recursion
            #pragma omp task default(shared)
            tree_traversal (userSeqFct, userVecFct, userCommFct, userArgs,
                            userCommArgs, *tree.right);
            #pragma omp task default(shared)
            tree_traversal (userSeqFct, userVecFct, userCommFct, userArgs,
                            userCommArgs, *tree.left);
            // Synchronization
            #pragma omp taskwait
        #elif CILK
            // Left & right recursion
            cilk_spawn
            tree_traversal (userSeqFct, userVecFct, userCommFct, userArgs,
                            userCommArgs, *tree.right);
            tree_traversal (userSeqFct, userVecFct, userCommFct, userArgs,
                            userCommArgs, *tree.left);
            // Synchronization
            cilk_sync;
        #endif

        // Separator recursion, if it is not empty
        if (tree.sep != nullptr) {
            tree_traversal (userSeqFct, userVecFct, userCommFct, userArgs,
                            userCommArgs, *tree.sep);
        }
    }
}

// Wrapper used to get the root of the D&C tree before calling the real tree traversal
void DC_tree_traversal (void (*userSeqFct)  (void *, DCargs_t *),
                        void (*userVecFct)  (void *, DCargs_t *),
                        void (*userCommFct) (void *, DCcommArgs_t *),
                        void *userArgs, void *userCommArgs)
{
    #ifdef OMP
        #pragma omp parallel
        #pragma omp single nowait
    #endif
    tree_traversal (userSeqFct, userVecFct, userCommFct, userArgs, userCommArgs,
                    *treeHead);
}
