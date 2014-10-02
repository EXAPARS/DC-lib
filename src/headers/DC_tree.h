#ifndef TREE_H
#define TREE_H

#ifdef STATS
#include <fstream>
#endif
#include "globals.h"

// Compare two different D&C trees and display the differences
void compare_dc_trees (tree_t &tree1, tree_t &tree2, int curNode);

// Compute the interval of edges of each leaf of the D&C tree
void edge_intervals (tree_t &tree, int *nodeToNodeRow, int *elemToNode, int nbNodes);

// Fill the element intervals of the D&C tree
void fill_dc_tree (tree_t &tree, int firstElem, int lastElem, int nbSepElem,
				   int firstNode, int lastNode, bool isLeaf);

// Create element partition & count left & separator elements
void create_elem_part (int *elemPart, int *nodePart, int *elemToNode, int nbElem,
                       int separator, int offset, int *nbLeftElem, int *nbSepElem);

// Compute the interval of nodes and elements at each node of the D&C tree
void create_dc_tree (tree_t &tree, int *elemToNode, int *sepToNode, int *nodePart,
                     int *nodePartSize, int globalNbElem, int firstPart, int lastPart,
                     int firstElem, int lastElem, int firstNode, int lastNode,
#ifdef STATS
					 int sepOffset, ofstream &dcFile, int curNode, int LRS);
#else
					 int sepOffset);
#endif

#endif
