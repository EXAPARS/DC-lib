#ifndef IO_H
#define IO_H

#include <fstream>
#include "globals.h"

// Read recursively each node of the D&C tree
void recursive_reading (tree_t &tree, ifstream &treeFile);

// Read the D&C tree and the permutation functions
void DC_read_tree (int nbElem, int nbNodes, int nbBlocks, int mpiRank);

// Store recursively each node of the D&C tree
void recursive_storing (tree_t &tree, ofstream &treeFile);

// Store the D&C tree and the permutation functions to a binary file
void DC_store_tree (int nbElem, int nbNodes, int nbBlocks, int mpiRank);

#endif
