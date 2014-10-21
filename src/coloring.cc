#ifdef HYBRID

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
    int nbColors = 0;

    // For each element of local interval
    for (int i = 0; i < nbElem; i++) {
        int neighborColor[NB_BLOCKS] = {0}, color = 0, mask = 1, block = 0;

        // Get the color of all neigbor elements
        for (int j = 0; j < elemToElem[i].size; j++) {
            int neighbor = elemToElem[i].list[j];
            neighborColor[colorPart[neighbor] / BLOCK_SIZE] |=
                 mask << (colorPart[neighbor] % BLOCK_SIZE);
        }
        // Get the first free color (position of the first 0 bit)
        while (neighborColor[block] & mask || colorCard[color] >= VEC_SIZE) {
            neighborColor[block] = neighborColor[block] >> 1;
            color++;
            if (color % BLOCK_SIZE == 0) block++;
        }
        if (color >= MAX_COLOR) {
            cerr << "Error: Not enough colors !\n";
            exit (EXIT_FAILURE);
        }
        // Assign the first free color to current element
        colorPart[i] = color;
        colorCard[color]++;

        // Compute the total number of colors
        if (color > nbColors) nbColors = color;
    }
    nbColors++;

    return nbColors;
}

// Construct element to element array from element to node and node to element
// Two elements are neighbors if they share a node
void elem_to_elem (list_t *elemToElem, index_t &nodeToElem, int *elemToNode,
                   int firstElem, int lastElem, int dimElem)
{
	// For each node of given element interval
    cilk_for (int i = firstElem; i <= lastElem; i++) {
        int nbNeighbors = 0;
        for (int j = 0; j < dimElem; j++) {
			int node = elemToNode[i*dimElem+j] - 1;
			// For each neighbor element of current node
            for (int k = nodeToElem.index[node];
                     k < nodeToElem.index[node+1]; k++) {
                int elemNeighbor = nodeToElem.value[k];
				// If neighbor element is in the element interval and is not
                // current element
                if (elemNeighbor >= firstElem && elemNeighbor <= lastElem &&
                    elemNeighbor != i) {
                    bool isNew = true;
					// If neighbor element is not already stored
                    for (int l = 0; l < nbNeighbors; l++) {
                        if ((elemNeighbor - firstElem) ==
                            elemToElem[i-firstElem].list[l]) {
							isNew = false;
                            break;
						}
					}
					// If neighbor is not in the list, add it
                    if (isNew) {
                       	elemToElem[i-firstElem].list[nbNeighbors] =
                            elemNeighbor - firstElem;
                        nbNeighbors++;
					}
				}
			}
		}
        elemToElem[i-firstElem].size = nbNeighbors;
	}
}

// Construct node to element structure from element to node
void node_to_elem (index_t &nodeToElem, int *elemToNode, int nbElem, int dimElem,
                   int nbNodes)
{
    couple_t *couple = new couple_t [nbElem * dimElem];
    cilk_for (int i = 0; i < nbElem; i++) {
        for (int j = 0; j < dimElem; j++) {
            couple[i*dimElem+j].elem = i;
            couple[i*dimElem+j].node = elemToNode[i*dimElem+j] - 1;
        }
    }
    quick_sort (couple, 0, nbElem * dimElem - 1);

    int ctr = 0;
    for (int i = 0; i < nbNodes; i++) {
        nodeToElem.index[i] = ctr;
        while (couple[ctr].node == i) {
            nodeToElem.value[ctr] = couple[ctr].elem;
            ctr++;
			if (ctr == nbElem * dimElem) break;
        }
    }
    nodeToElem.index[nbNodes] = ctr;
    delete[] couple;
}

// Create a coloring permutation for each leaf of the D&C tree
void leaves_coloring (tree_t &tree, index_t &nodeToElem, int *elemToNode,
#ifdef STATS
					  int *elemPerColor, ofstream &colorPerLeaf,
#endif
					  int globalNbElem, int dimElem)
{
	// If current node is a leaf
	if (tree.left == NULL && tree.right == NULL) {
        // List the neighbor elements of each element of the leaf
		int nbColors, localNbElem = tree.lastSep - tree.firstElem + 1;
        list_t *elemToElem = new list_t [localNbElem];
        elem_to_elem (elemToElem, nodeToElem, elemToNode, tree.firstElem,
                      tree.lastSep, dimElem);

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
		permute_int_2d_array (elemToNode, colorPerm, localNbElem, dimElem,
                              tree.firstElem);
		merge_permutations (colorPerm, globalNbElem, localNbElem, tree.firstElem,
                            tree.lastSep);
		delete[] colorPerm;
	}
	else {
		cilk_spawn
#ifdef STATS
		leaves_coloring (*tree.left, nodeToElem, elemToNode, elemPerColor,
						 colorPerLeaf, globalNbElem, dimElem);
		leaves_coloring (*tree.right, nodeToElem, elemToNode, elemPerColor,
						 colorPerLeaf, globalNbElem, dimElem);
#else
		leaves_coloring (*tree.left, nodeToElem, elemToNode, globalNbElem, dimElem);
		leaves_coloring (*tree.right, nodeToElem, elemToNode, globalNbElem, dimElem);
#endif
		cilk_sync;
		if (tree.sep != NULL) {
#ifdef STATS
			leaves_coloring (*tree.sep, nodeToElem, elemToNode, elemPerColor,
							 colorPerLeaf, globalNbElem, dimElem);
#else
            leaves_coloring (*tree.sep, nodeToElem, elemToNode, globalNbElem, dimElem);
#endif
		}
	}
}

// Coloring of the D&C tree
void coloring (int *elemToNode, int nbElem, int dimElem, int nbNodes)
{
    // List the neighbor elements of each node
    index_t nodeToElem;
    nodeToElem.index = new int [nbNodes + 1];
    nodeToElem.value = new int [nbElem * dimElem];
    node_to_elem (nodeToElem, elemToNode, nbElem, dimElem, nbNodes);

    #ifdef STATS
        string fileName = "colorPerLeaf_" + meshName + "_" +
                          to_string ((long long)MAX_ELEM_PER_PART) + ".csv";
        ofstream colorPerLeaf (fileName, ios::out | ios::trunc);
        colorPerLeaf << "leafNb nbColors\n";
        int *elemPerColor = new int [MAX_ELEM_PER_PART] ();
        leaves_coloring (*treeHead, nodeToElem, elemToNode, elemPerColor, colorPerLeaf,
                         nbElem, dimElem);
        coloring_stat (elemPerColor, nbElem);
        delete[] elemPerColor;
        colorPerLeaf.close ();
    #else
        leaves_coloring (*treeHead, nodeToElem, elemToNode, nbElem, dimElem);
    #endif

    delete[] nodeToElem.value, delete[] nodeToElem.index;
}

#endif
