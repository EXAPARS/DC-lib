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

#ifdef MULTITHREADED_COMM

// Determine if nodeID is a descendant of current node
bool isDescendant (int curNode, int nodeID);

// Initialize D&C tree interface for multithreaded communication
void create_multithreaded_intf (tree_t &tree, int *elemToNode, int *intfIndex,
                                int *intfNodes, int dimElem, int nbIntf, int nbBlocks,
                                int curNode, int curLevel, bool isLeaf);

// Compute the number of nodes owned by current leaf and fill the list
void create_owned_nodes_list (tree_t &tree, int *elemToNode, int dimElem, int curNode);

#endif

// Compute the edge interval, the list of nodes owned by each leaf of the D&C tree,
// and the interface for multithreaded communication
void tree_finalize (tree_t &tree, int *nodeToNodeRow, int *elemToNode, int *intfIndex,
                    int *intfNodes, int dimElem, int nbBlocks, int nbIntf, int curNode,
                    int curLevel);

// Wrapper used to get the root of the D&C tree before calling the real tree finalize
void DC_finalize_tree (int *nodeToNodeRow, int *elemToNode, int *intfIndex,
                       int *intfNodes, int dimElem, int nbBlocks, int nbIntf);

// Initialize the content of D&C tree nodes
void init_dc_tree (tree_t &tree, int firstElem, int lastElem, int nbSepElem,
                   int firstNode, int lastNode, bool isSep, bool isLeaf);

#ifdef MULTITHREADED_COMM
// Set the last updater of each node. The owners are the closer leaves to the root.
void fill_node_owner (int *elemToNode, int firstElem, int lastElem, int dimElem,
                      int firstNode, int lastNode, int curNode, bool isSep);

// Initialize the owner of the nodes to MAX_INT (lower is owner)
void init_node_owner (int nbNodes);
#endif

// Create element partition & count left & separator elements
void create_elem_part (int *elemPart, int *nodePart, int *elemToNode, int nbElem,
                       int dimElem, int separator, int offset, int *nbLeftElem,
                       int *nbSepElem);

// Create the D&C tree and the element permutation, and compute the intervals of nodes
// and elements at each node of the tree
void tree_creation (tree_t &tree, int *elemToNode, int *sepToNode, int *nodePart,
                    int *nodePartSize, int globalNbElem, int dimElem, int firstPart,
                    int lastPart, int firstElem, int lastElem, int firstNode,
                    int lastNode, int sepOffset, int curNode, bool isSep
#ifdef STATS
                    , ofstream &dcFile, int LRS);
#else
                    );
#endif

// Create the D&C tree and the permutations
void DC_create_tree (int *elemToNode, int nbElem, int dimElem, int nbNodes, int rank);

#endif
