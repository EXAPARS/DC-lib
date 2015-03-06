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

#ifndef TREE_TRAVERSAL_H
#define TREE_TRAVERSAL_H

#include "DC.h"

// Follow the D&C tree to execute the given function in parallel
void tree_traversal (void (*userSeqFct) (void *, int, int),
                     void (*userVecFct) (void *, int, int), void *userArgs,
                     double *nodeToNodeValue, int operatorDim, tree_t &tree);

// Wrapper used to get the root of the D&C tree before calling the real tree traversal
void DC_tree_traversal (void (*userSeqFct) (void *, int, int),
                        void (*userVecFct) (void *, int, int),
                        void *userArgs, double *nodeToNodeValue, int operatorDim);

#endif
