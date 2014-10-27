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
