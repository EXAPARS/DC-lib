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

#ifndef TREE_H
#define TREE_H

#ifdef STATS
    #include <fstream>
#endif
#include "DC.h"

// Compute the edge interval for each leaf of the D&C tree
void compute_edge_intervals (tree_t &tree, int *nodeToNodeRow, int *elemToNode,
                             int nbNodes);

// Wrapper used to get the root of the D&C tree before computing the edge intervals
// for CSR reset
void DC_finalize_tree (int *nodeToNodeRow, int *elemToNode, int nbNodes);

// Initialize a node of the D&C tree
void init_dc_tree (tree_t &tree, int firstElem, int lastElem, int nbSepElem,
				   int firstNode, int lastNode, bool isLeaf);

// Create element partition & count left & separator elements
void create_elem_part (int *elemPart, int *nodePart, int *elemToNode, int nbElem,
                       int dimElem, int separator, int offset, int *nbLeftElem,
                       int *nbSepElem);

// Compute the interval of nodes and elements at each node of the D&C tree
void recursive_tree_creation (tree_t &tree, int *elemToNode, int *sepToNode,
                              int *nodePart, int *nodePartSize, int globalNbElem,
                              int dimElem, int firstPart, int lastPart, int firstElem,
                              int lastElem, int firstNode, int lastNode, int sepOffset
#ifdef STATS
                              , ofstream &dcFile, int curNode, int LRS);
#else
                              );
#endif

// Create the D&C tree and the permutations
void DC_create_tree (int *elemToNode, int nbElem, int dimElem, int nbNodes);

#endif
