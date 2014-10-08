#include <iostream>
#include <cilk/cilk.h>

#include "tools.h"
#include "permutations.h"
#include "coloring.h"

#ifdef STATS
    #include <pthread.h>
    pthread_mutex_t statMutex = PTHREAD_MUTEX_INITIALIZER;
#endif

// Assign a color to the elements of a given leaf & return the number of colors
int create_color_part (int *colorPart, int *colorCard, list_t *elemToElem, int nbElem)
{
	__uint128_t *elemToColor = new __uint128_t [nbElem] ();
	int nbColors = 0;

    // For each element of local interval
	for (int i = 0; i < nbElem; i++) {
		__uint128_t mask = 1, neighborColor = 0;
		int color = 0;

        // Get the color of all neigbor elements
        for (int j = 0; j < elemToElem[i].size; j++) {
            int neighbor = elemToElem[i].list[j];
            neighborColor |= elemToColor[neighbor];
		}
        // Get the first free color (position of the first 0 bit)
		while (neighborColor & mask || colorCard[color] >= VEC_SIZE) {
		    neighborColor = neighborColor >> 1;
		    color++;
		}
		if (color >= MAX_COLOR) {
		    cerr << "Error: Not enough colors !\n";
		    exit (EXIT_FAILURE);
        }
        // Assign the first free color to current element
		elemToColor[i] = (mask << color);
        colorPart[i] = color;
        colorCard[color]++;

        // Compute the total number of colors
        if (color > nbColors) nbColors = color;
	}
    nbColors++;

	delete[] elemToColor;
    return nbColors;
}

// Create a coloring permutation for each leaf of the D&C tree
void leaves_coloring (tree_t &tree, index_t &nodeToElem, int *elemToNode,
#ifdef STATS
					  int *elemPerColor, ofstream &colorPerLeaf,
#endif
					  int globalNbElem)
{
	// If current node is a leaf
	if (tree.left == NULL && tree.right == NULL) {
        // List the neighbor elements of each element of the leaf
		int nbColors, localNbElem = tree.lastSep - tree.firstElem + 1;
        list_t *elemToElem = new list_t [localNbElem];
        elem_to_elem (elemToElem, nodeToElem, elemToNode, tree.firstElem,
                      tree.lastSep);

		// Assign a color to each element of the leaf
		int *colorPart = new int [localNbElem];
        int colorCard[MAX_COLOR] = {0};
    	nbColors = create_color_part (colorPart, colorCard, elemToElem, localNbElem);
        delete[] elemToElem;

#ifdef STATS
		// Compute statistics on the coloring of the leaf
        pthread_mutex_lock (&statMutex);
		leaf_coloring_stat (colorPerLeaf, elemPerColor, colorPart, localNbElem,
                            nbColors);
        pthread_mutex_unlock (&statMutex);
#endif

		// Create a local permutation on the elements of the leaf & get the index
        // of the last element in a full vectorial color
		int *colorPerm = new int [localNbElem];
        int vecOffset  = create_coloring_perm (colorPerm, colorPart, colorCard,
                                               localNbElem, nbColors);
        tree.vecOffset = vecOffset + tree.firstElem;
		delete[] colorPart;

		// Apply local permutation on elemToNode & on the global element permutation
		permute_int_2d_array (elemToNode, colorPerm, localNbElem, DIM_ELEM,
                              tree.firstElem);
		merge_permutations (colorPerm, globalNbElem, localNbElem, tree.firstElem,
                            tree.lastSep);
		delete[] colorPerm;
	}
	else {
		cilk_spawn
#ifdef STATS
		leaves_coloring (*tree.left, nodeToElem, elemToNode, elemPerColor,
						 colorPerLeaf, globalNbElem);
		leaves_coloring (*tree.right, nodeToElem, elemToNode, elemPerColor,
						 colorPerLeaf, globalNbElem);
#else
		leaves_coloring (*tree.left, nodeToElem, elemToNode, globalNbElem);
		leaves_coloring (*tree.right, nodeToElem, elemToNode, globalNbElem);
#endif
		cilk_sync;
		if (tree.sep != NULL) {
#ifdef STATS
			leaves_coloring (*tree.sep, nodeToElem, elemToNode, elemPerColor,
							 colorPerLeaf, globalNbElem);
#else
            leaves_coloring (*tree.sep, nodeToElem, elemToNode, globalNbElem);
#endif
		}
	}
}

// Coloring of the D&C tree
void coloring (int *elemToNode, int nbElem, int nbNodes)
{
    // List the neighbor elements of each node
    index_t nodeToElem;
    nodeToElem.index = new int [nbNodes + 1];
    nodeToElem.value = new int [nbElem * DIM_ELEM];
    node_to_elem (nodeToElem, elemToNode, nbElem, nbNodes);

    #ifdef STATS
        string fileName = "colorPerLeaf_" + meshName + "_" +
                          to_string ((long long)MAX_ELEM_PER_PART) + ".csv";
        ofstream colorPerLeaf (fileName, ios::out | ios::trunc);
        colorPerLeaf << "leafNb nbColors\n";
        int *elemPerColor = new int [MAX_ELEM_PER_PART] ();
        leaves_coloring (*treeHead, nodeToElem, elemToNode, elemPerColor, colorPerLeaf,
                         nbElem);
        coloring_stat (elemPerColor, nbElem);
        delete[] elemPerColor;
        colorPerLeaf.close ();
    #else
        leaves_coloring (*treeHead, nodeToElem, elemToNode, nbElem);
    #endif

    delete[] nodeToElem.value, delete[] nodeToElem.index;
}
