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

#ifndef TREE_CREATION_H
#define TREE_CREATION_H

#ifdef STATS
    #include <fstream>
#endif
#include "DC.h"

// Compute the edge interval and the list of nodes owned by each leaf of the D&C tree
void compute_intervals (tree_t &tree, int *nodeToNodeRow, int *elemToNode,int curNode);

// Wrapper used to get the root of the D&C tree before computing the edge intervals and
// the list of nodes owned by each leaf of the D&C tree
void DC_finalize_tree (int *nodeToNodeRow, int *elemToNode);

// Initialize a node of the D&C tree
void init_dc_tree (tree_t &tree, int *elemToNode, int *intfIndex, int *intfNodes,
                   int firstElem, int lastElem, int nbSepElem, int dimElem,
                   int firstNode, int lastNode, int nbIntf, int commDepth,
                   int curDepth, bool isSep, bool isLeaf, bool *hasIntfNode);

// Create element partition & count left & separator elements
void create_elem_part (int *elemPart, int *nodePart, int *elemToNode, int nbElem,
                       int dimElem, int separator, int offset, int *nbLeftElem,
                       int *nbSepElem);

// Create the D&C tree and the element permutation, and compute the intervals of nodes
// and elements at each node of the tree
void tree_creation (tree_t &tree, int *elemToNode, int *sepToNode, int *nodePart,
                    int *nodePartSize, int *intfIndex, int *intfNodes,
                    int globalNbElem, int dimElem, int firstPart, int lastPart,
                    int firstElem, int lastElem, int firstNode, int lastNode,
                    int sepOffset, int nbIntf, int curNode, int commDepth,
#ifdef STATS
                    int curDepth, bool isSep, ofstream &dcFile, int LRS);
#else
                    int curDepth, bool isSep);
#endif

// Create the D&C tree and the permutations
void DC_create_tree (double *coord, int *elemToNode, int *intfIndex, int *intfNodes,
                     int nbElem, int dimElem, int nbNodes, int dimNode, int nbIntf,
                     int rank);

#endif
