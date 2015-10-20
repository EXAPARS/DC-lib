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
void recursive_reading (tree_t &tree, ifstream &treeFile, int nbIntf)
{
    bool isLeaf;
    tree.nbOwnedNodes = -1;
    tree.intfIndex    = nullptr;
    tree.intfNodes    = nullptr;
    tree.ownedNodes   = nullptr;
    tree.left         = nullptr;
    tree.right        = nullptr;
    tree.sep          = nullptr;

    treeFile.read ((char*)&tree.firstElem,    sizeof (int));
    treeFile.read ((char*)&tree.lastElem,     sizeof (int));
    treeFile.read ((char*)&tree.lastSep,      sizeof (int));
    treeFile.read ((char*)&tree.firstNode,    sizeof (int));
    treeFile.read ((char*)&tree.lastNode,     sizeof (int));
    treeFile.read ((char*)&tree.firstEdge,    sizeof (int));
    treeFile.read ((char*)&tree.lastEdge,     sizeof (int));
    treeFile.read ((char*)&tree.nbIntfNodes,  sizeof (int));
    if (tree.nbIntfNodes > 0) {
        tree.intfIndex = new int [nbIntf+1];
        tree.intfNodes = new int [tree.nbIntfNodes];
        treeFile.read ((char*)tree.intfIndex,     (nbIntf + 1) * sizeof (int));
        treeFile.read ((char*)tree.intfNodes, tree.nbIntfNodes * sizeof (int));
        treeFile.read ((char*)tree.intfDest,  tree.nbIntfNodes * sizeof (int));
    }
    treeFile.read ((char*)&tree.vecOffset, sizeof (int));
    treeFile.read ((char*)&tree.isSep,     sizeof (bool));
    treeFile.read ((char*)&isLeaf,         sizeof (bool));

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
        recursive_reading (*tree.left,  treeFile, nbIntf);
        recursive_reading (*tree.right, treeFile, nbIntf);
        if (tree.sep != nullptr) {
            recursive_reading (*tree.sep, treeFile, nbIntf);
        }
    }
}

// Read the D&C tree and the permutation functions
void DC_read_tree (string &treePath, int nbElem, int nbNodes, int nbIntf)
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
	recursive_reading (*treeHead, treeFile, nbIntf);
	treeFile.close ();
}

// Store recursively each node of the D&C tree
void recursive_storing (tree_t &tree, ofstream &treeFile, int nbIntf)
{
    bool isLeaf;
    treeFile.write ((char*)&tree.firstElem,   sizeof (int));
    treeFile.write ((char*)&tree.lastElem,    sizeof (int));
    treeFile.write ((char*)&tree.lastSep,     sizeof (int));
    treeFile.write ((char*)&tree.firstNode,   sizeof (int));
    treeFile.write ((char*)&tree.lastNode,    sizeof (int));
    treeFile.write ((char*)&tree.firstEdge,   sizeof (int));
    treeFile.write ((char*)&tree.lastEdge,    sizeof (int));
    treeFile.write ((char*)&tree.nbIntfNodes, sizeof (int));
    if (tree.nbIntfNodes > 0) {
        treeFile.write ((char*)tree.intfIndex,     (nbIntf + 1) * sizeof (int));
        treeFile.write ((char*)tree.intfNodes, tree.nbIntfNodes * sizeof (int));
        treeFile.write ((char*)tree.intfDest,  tree.nbIntfNodes * sizeof (int));
    }
    treeFile.write ((char*)&tree.vecOffset, sizeof (int));
    treeFile.write ((char*)&tree.isSep,     sizeof (bool));

    if (tree.left == nullptr && tree.right == nullptr) {
        isLeaf = true;
        treeFile.write ((char*)&isLeaf,            sizeof (bool));
        treeFile.write ((char*)&tree.nbOwnedNodes, sizeof (int));
        if (tree.nbOwnedNodes > 0) {
            treeFile.write ((char*)tree.ownedNodes, tree.nbOwnedNodes * sizeof (int));
        }
    }
    else {
        isLeaf = false;
        treeFile.write ((char*)&isLeaf, sizeof (bool));

        // Left, right & separator recursion
        recursive_storing (*tree.left,  treeFile, nbIntf);
        recursive_storing (*tree.right, treeFile, nbIntf);
        if (tree.sep != nullptr) {
            recursive_storing (*tree.sep, treeFile, nbIntf);
        }
    }
}

// Store the D&C tree and the permutation functions to a binary file
void DC_store_tree (string &treePath, int nbElem, int nbNodes, int nbIntf)
{
    ofstream treeFile (treePath, ios::out | ios::trunc | ios::binary);
    if (!treeFile.is_open ()) {
        cerr << "Error: cannot store DC tree & permutations!\n";
        exit (EXIT_FAILURE);
    }
    treeFile.write ((char*)elemPerm, nbElem  * sizeof (int));
    treeFile.write ((char*)nodePerm, nbNodes * sizeof (int));
    recursive_storing (*treeHead, treeFile, nbIntf);
    treeFile.close ();

    delete[] nodePerm, delete[] elemPerm; 
}
