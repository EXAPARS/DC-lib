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
