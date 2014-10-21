#include "IO.h"

// Read recursively each node of the D&C tree
void recursive_reading (tree_t &tree, ifstream &treeFile)
{
	bool isLeaf;
	treeFile.read ((char*)&tree.firstElem, sizeof (int));
	treeFile.read ((char*)&tree.lastElem,  sizeof (int));
	treeFile.read ((char*)&tree.lastSep,   sizeof (int));
    treeFile.read ((char*)&tree.firstCSR,  sizeof (int));
    treeFile.read ((char*)&tree.lastCSR,   sizeof (int));
    #ifdef HYBRID
        treeFile.read ((char*)&tree.vecOffset, sizeof (int));
    #else
        tree.vecOffset = 0;
    #endif
	treeFile.read ((char*)&isLeaf, sizeof (bool));

	tree.left  = NULL;
	tree.right = NULL;
	tree.sep   = NULL;

	if (isLeaf == false) {
		tree.left  = new tree_t;
		tree.right = new tree_t;
        if (tree.lastSep > tree.lastElem) {
			tree.sep = new tree_t;
		}

		// Left, right & separator recursion
		recursive_reading (*tree.left,  treeFile);
		recursive_reading (*tree.right, treeFile);
		if (tree.sep != NULL) {
			recursive_reading (*tree.sep, treeFile);
		}
	}
}

// Read the D&C tree and the permutation functions
void DC_read_tree (int nbElem, int nbNodes, int nbBlocks, int mpiRank)
{
	string fileName = "../data" + meshName + "/DC_tree/" +
    #ifdef HYBRID
                      "Hybrid_" +
    #else
                      "DC_" +
    #endif
                      to_string ((long long)MAX_ELEM_PER_PART) + "_" +
                      to_string ((long long)nbBlocks) + "_" +
                      to_string ((long long)mpiRank);

	ifstream treeFile (fileName, ios::in | ios::binary);
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
    treeFile.write ((char*)&tree.firstCSR,  sizeof (int));
    treeFile.write ((char*)&tree.lastCSR,   sizeof (int));
    #ifdef HYBRID
        treeFile.write ((char*)&tree.vecOffset, sizeof (int));
    #endif

	if (tree.left == NULL && tree.right == NULL) {
		isLeaf = true;
		treeFile.write ((char*)&isLeaf, sizeof (bool));
	}
	else {
        isLeaf = false;
        treeFile.write ((char*)&isLeaf, sizeof (bool));

		// Left, right & separator recursion
		recursive_storing (*tree.left,  treeFile);
		recursive_storing (*tree.right, treeFile);
		if (tree.sep != NULL) {
			recursive_storing (*tree.sep, treeFile);
		}
	}
}

// Store the D&C tree and the permutation functions to a binary file
void DC_store_tree (int nbElem, int nbNodes, int nbBlocks, int mpiRank)
{
    string fileName = "../data" + meshName + "/DC_tree/" +
    #ifdef HYBRID
                      "Hybrid_" +
    #else
                      "DC_" +
    #endif
                      to_string ((long long)MAX_ELEM_PER_PART) + "_" +
                      to_string ((long long)nbBlocks) + "_" +
                      to_string ((long long)mpiRank);

    ofstream treeFile (fileName, ios::out | ios::trunc | ios::binary);
    treeFile.write ((char*)elemPerm, nbElem  * sizeof (int));
    treeFile.write ((char*)nodePerm, nbNodes * sizeof (int));
    recursive_storing (*treeHead, treeFile);
    treeFile.close ();

    delete[] nodePerm, delete[] elemPerm; 
}
