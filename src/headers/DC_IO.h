#ifndef DC_IO_H
#define DC_IO_H

#include <fstream>
#include "globals.h"

// Read recursively the intervals of each node of the D&C tree
void read_dc_tree (tree_t &tree, ifstream &permAndTree);

// Read the permutation functions and the D&C tree
void read_perm_and_tree (int nbElem, int nbNodes, int nbBlocks, int mpiRank);

// Store recursively the intervals of each node of the D&C tree
void store_dc_tree (tree_t &tree, ofstream &permAndTree);

// Store the permutation functions and the D&C tree to a binary file
void store_perm_and_tree (int nbElem, int nbNodes, int nbBlocks, int mpiRank);

#endif
