#include "dc_IO.h"

// Read recursively the intervals of each node of the D&C tree
void read_dc_tree (tree_t &tree, ifstream &permAndTree)
{
	bool isLeaf;
	permAndTree.read ((char*)&tree.firstElem, sizeof (int));
	permAndTree.read ((char*)&tree.lastElem,  sizeof (int));
	permAndTree.read ((char*)&tree.lastSep,   sizeof (int));
    permAndTree.read ((char*)&tree.firstCSR,  sizeof (int));
    permAndTree.read ((char*)&tree.lastCSR,   sizeof (int));
#ifdef HYBRID
    permAndTree.read ((char*)&tree.lastFullColor,  sizeof (int));
#else
    tree.lastFullColor = 0;
#endif
	permAndTree.read ((char*)&isLeaf, sizeof (bool));

	tree.left   = NULL;
	tree.right  = NULL;
	tree.sep    = NULL;

	if (isLeaf == false) {
		tree.left  = new tree_t;
		tree.right = new tree_t;
        if (tree.lastSep > tree.lastElem) {
			tree.sep = new tree_t;
		}

		// Left, right & separator recursion
		read_dc_tree (*tree.left,  permAndTree);
		read_dc_tree (*tree.right, permAndTree);
		if (tree.sep != NULL) {
			read_dc_tree (*tree.sep, permAndTree);
		}
	}
}

// Read the permutation functions and the D&C tree
void read_perm_and_tree (int nbElem, int nbNodes, int nbBlocks, int mpiRank)
{
	string fileName = meshName + "/permAndTree/" +
#ifdef HYBRID
                      "Hybrid_" +
#else
                      "DC_" +
#endif
                      to_string ((long long)MAX_ELEM_PER_PART) + "_" +
                      to_string ((long long)nbBlocks) + "_" +
                      to_string ((long long)mpiRank);
	ifstream permAndTree (fileName, ios::in | ios::binary);
	permAndTree.read ((char*)elemPerm, nbElem  * sizeof (int));
	permAndTree.read ((char*)nodePerm, nbNodes * sizeof (int));
	read_dc_tree (*treeHead, permAndTree);
	permAndTree.close ();
}

// Store recursively the intervals of each node of the D&C tree
void store_dc_tree (tree_t &tree, ofstream &permAndTree)
{
    bool isLeaf;
    permAndTree.write ((char*)&tree.firstElem, sizeof (int));
    permAndTree.write ((char*)&tree.lastElem,  sizeof (int));
    permAndTree.write ((char*)&tree.lastSep,   sizeof (int));
    permAndTree.write ((char*)&tree.firstCSR,  sizeof (int));
    permAndTree.write ((char*)&tree.lastCSR,   sizeof (int));
#ifdef HYBRID
    permAndTree.write ((char*)&tree.lastFullColor,  sizeof (int));
#endif

	// If current node is a leaf, dump local leaf interval
	if (tree.left == NULL && tree.right == NULL) {
		isLeaf = true;
		permAndTree.write ((char*)&isLeaf, sizeof (bool));
	}
	else {
        isLeaf = false;
        permAndTree.write ((char*)&isLeaf, sizeof (bool));

		// Left, right & separator recursion
		store_dc_tree (*tree.left,  permAndTree);
		store_dc_tree (*tree.right, permAndTree);
		if (tree.sep != NULL) {
			store_dc_tree (*tree.sep, permAndTree);
		}
	}
}

// Store the permutation functions and the D&C tree to a binary file
void store_perm_and_tree (int nbElem, int nbNodes, int nbBlocks, int mpiRank)
{
    string fileName = meshName + "/permAndTree/" +
#ifdef HYBRID
                      "Hybrid_" +
#else
                      "DC_" +
#endif
                      to_string ((long long)MAX_ELEM_PER_PART) + "_" +
                      to_string ((long long)nbBlocks) + "_" +
                      to_string ((long long)mpiRank);
    ofstream permAndTree (fileName, ios::out | ios::trunc | ios::binary);
    permAndTree.write ((char*)elemPerm, nbElem  * sizeof (int));
    permAndTree.write ((char*)nodePerm, nbNodes * sizeof (int));
    store_dc_tree (*treeHead, permAndTree);
    permAndTree.close ();
}
