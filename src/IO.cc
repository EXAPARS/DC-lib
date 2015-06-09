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

#include <iostream>
#include "IO.h"

extern tree_t *treeHead;
extern int *elemPerm, *nodePerm;

// Read recursively each node of the D&C tree
void recursive_reading (tree_t &tree, ifstream &treeFile)
{
    bool isLeaf;
    treeFile.read ((char*)&tree.firstElem, sizeof (int));
    treeFile.read ((char*)&tree.lastElem,  sizeof (int));
    treeFile.read ((char*)&tree.lastSep,   sizeof (int));
    treeFile.read ((char*)&tree.firstNode, sizeof (int));
    treeFile.read ((char*)&tree.lastNode,  sizeof (int));
    treeFile.read ((char*)&tree.firstEdge, sizeof (int));
    treeFile.read ((char*)&tree.lastEdge,  sizeof (int));
    #ifdef DC_VEC
        treeFile.read ((char*)&tree.vecOffset, sizeof (int));
    #endif
    treeFile.read ((char*)&tree.isSep, sizeof (bool));
    treeFile.read ((char*)&isLeaf, sizeof (bool));

    tree.nbOwnedNodes = -1;
    tree.ownedNodes   = nullptr;
    tree.left         = nullptr;
    tree.right        = nullptr;
    tree.sep          = nullptr;

    if (isLeaf) {
        treeFile.read ((char*)&tree.nbOwnedNodes, sizeof (int));
        if (tree.nbOwnedNodes > 0) {
            tree.ownedNodes = new int [tree.nbOwnedNodes];
            treeFile.read ((char*)tree.ownedNodes, tree.nbOwnedNodes * sizeof (int));
        }
    }
    else {
        tree.left  = new tree_t;
        tree.right = new tree_t;
        if (tree.lastSep > tree.lastElem) {
            tree.sep = new tree_t;
        }

        // Left, right & separator recursion
        recursive_reading (*tree.left,  treeFile);
        recursive_reading (*tree.right, treeFile);
        if (tree.sep != nullptr) {
            recursive_reading (*tree.sep, treeFile);
        }
    }
}

// Read the D&C tree and the permutation functions
void DC_read_tree (string &treePath, int nbElem, int nbNodes)
{
    // Allocate the D&C tree & the permutation functions
    treeHead = new tree_t;
    elemPerm = new int [nbElem];
    nodePerm = new int [nbNodes];

	ifstream treeFile (treePath, ios::in | ios::binary);
    if (!treeFile.is_open ()) {
        cerr << "Error: cannot read DC tree & permutations!\n";
        exit (EXIT_FAILURE);
    }
	treeFile.read ((char*)elemPerm, nbElem  * sizeof (int));
	treeFile.read ((char*)nodePerm, nbNodes * sizeof (int));
	recursive_reading (*treeHead, treeFile);
	treeFile.close ();
}

// Store recursively each node of the D&C tree
void recursive_storing (tree_t &tree, ofstream &treeFile)
{
    bool isLeaf;
    treeFile.write ((char*)&tree.firstElem, sizeof (int));
    treeFile.write ((char*)&tree.lastElem,  sizeof (int));
    treeFile.write ((char*)&tree.lastSep,   sizeof (int));
    treeFile.write ((char*)&tree.firstNode, sizeof (int));
    treeFile.write ((char*)&tree.lastNode,  sizeof (int));
    treeFile.write ((char*)&tree.firstEdge, sizeof (int));
    treeFile.write ((char*)&tree.lastEdge,  sizeof (int));
    #ifdef DC_VEC
        treeFile.write ((char*)&tree.vecOffset, sizeof (int));
    #endif
    treeFile.write ((char*)&tree.isSep, sizeof (bool));

    if (tree.left == nullptr && tree.right == nullptr) {
        isLeaf = true;
        treeFile.write ((char*)&isLeaf, sizeof (bool));
        treeFile.write ((char*)&tree.nbOwnedNodes, sizeof (int));
        if (tree.nbOwnedNodes > 0) {
            treeFile.write ((char*)tree.ownedNodes, tree.nbOwnedNodes * sizeof (int));
        }
    }
    else {
        isLeaf = false;
        treeFile.write ((char*)&isLeaf, sizeof (bool));

        // Left, right & separator recursion
        recursive_storing (*tree.left,  treeFile);
        recursive_storing (*tree.right, treeFile);
        if (tree.sep != nullptr) {
            recursive_storing (*tree.sep, treeFile);
        }
    }
}

// Store the D&C tree and the permutation functions to a binary file
void DC_store_tree (string &treePath, int nbElem, int nbNodes)
{
    ofstream treeFile (treePath, ios::out | ios::trunc | ios::binary);
    if (!treeFile.is_open ()) {
        cerr << "Error: cannot store DC tree & permutations!\n";
        exit (EXIT_FAILURE);
    }
    treeFile.write ((char*)elemPerm, nbElem  * sizeof (int));
    treeFile.write ((char*)nodePerm, nbNodes * sizeof (int));
    recursive_storing (*treeHead, treeFile);
    treeFile.close ();

    delete[] nodePerm, delete[] elemPerm; 
}
