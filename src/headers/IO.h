#ifndef IO_H
#define IO_H

#include <fstream>
#include "DC.h"

using namespace std;

// Read recursively each node of the D&C tree
void recursive_reading (tree_t &tree, ifstream &treeFile);

// Read the D&C tree and the permutation functions
void DC_read_tree (string &treePath, int nbElem, int nbNodes);

// Store recursively each node of the D&C tree
void recursive_storing (tree_t &tree, ofstream &treeFile);

// Store the D&C tree and the permutation functions to a binary file
void DC_store_tree (string &treePath, int nbElem, int nbNodes);

#endif
