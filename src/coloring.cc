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

#ifdef CILK
    #include <cilk/cilk.h>
#endif
#include <iostream>
#include <string.h>

#include "tools.h"
#include "permutations.h"
#include "coloring.h"

#ifdef STATS
    #include <pthread.h>
    pthread_mutex_t statMutex = PTHREAD_MUTEX_INITIALIZER;
#endif

extern tree_t *treeHead;

// Create element to element array from element to node and node to element
// Two elements are neighbors if they share a node
void DC_create_elemToElem (list_t *elemToElem, index_t &nodeToElem, int *elemToNode,
                           int firstElem, int lastElem, int dimElem)
{
	// For each node of given element interval
    #ifdef OMP
        #pragma omp parallel for
        for (int i = firstElem; i <= lastElem; i++) {
    #elif CILK
        cilk_for (int i = firstElem; i <= lastElem; i++) {
    #endif
        int nbNeighbors = 0;
        for (int j = 0; j < dimElem; j++) {
		    int node = elemToNode[i*dimElem+j] - 1;
		    // For each neighbor element of current node
            for (int k = nodeToElem.index[node]; k < nodeToElem.index[node+1]; k++) {
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

// Create node to element structure from element to node
void DC_create_nodeToElem (index_t &nodeToElem, int *elemToNode, int nbElem,
                          int dimElem, int nbNodes)
{
    couple_t *couple = new couple_t [nbElem * dimElem];
    
    #ifdef OMP
        #pragma omp parallel for
        for (int i = 0; i < nbElem; i++) {
    #elif CILK
        cilk_for (int i = 0; i < nbElem; i++) {
    #endif
        for (int j = 0; j < dimElem; j++) {
            couple[i*dimElem+j].elem = i;
            couple[i*dimElem+j].node = elemToNode[i*dimElem+j] - 1;
        }
    }
    #ifdef OMP
        #pragma omp parallel
        #pragma omp single nowait
    #endif
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

#ifdef DC_VEC

// Assign a color to the elements of a given leaf following the bounded colors strategy
// & return the number of colors
int create_bounded_color_part (int *colorPart, int *colorCard, list_t *elemToElem,
                               int nbElem)
{
    int nbColors = 0;

    // For each element of local interval
    for (int i = 0; i < nbElem; i++) {
        int neighborsColor[NB_BLOCKS] = {0}, color = 0, mask = 1, block = 0;

        // Get the color of all neigbor elements
        for (int j = 0; j < elemToElem[i].size; j++) {
            int neighborColor = colorPart[elemToElem[i].list[j]];
            // If neighbor is initialized
            if (neighborColor != -1) {
                neighborsColor[neighborColor / BLOCK_SIZE] |=
                      mask << (neighborColor % BLOCK_SIZE);
            }
        }
        // Get the first free color (position of the first 0 bit)
        while (neighborsColor[block] & mask || colorCard[color] >= VEC_SIZE) {
            neighborsColor[block] = neighborsColor[block] >> 1;
            color++;
            if (color % BLOCK_SIZE == 0) block++;
        }
        if (color >= MAX_COLOR) {
            cerr << "Error: Not enough colors.\n";
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

// Create a coloring permutation for each leaf of the D&C tree
void leaves_coloring (tree_t &tree, index_t &nodeToElem, int *elemToNode,
#ifdef STATS
					  int *elemPerColor, ofstream &colorPerLeaf,
#endif
					  int globalNbElem, int dimElem)
{
	// If current node is a leaf
	if (tree.left == nullptr && tree.right == nullptr) {
        // List the neighbor elements of each element of the leaf
		int nbColors, localNbElem = tree.lastSep - tree.firstElem + 1;
        list_t *elemToElem = new list_t [localNbElem];
        DC_create_elemToElem (elemToElem, nodeToElem, elemToNode, tree.firstElem,
                              tree.lastSep, dimElem);

		// Assign a color to each element of the leaf
		int *colorPart = new int [localNbElem];
        int colorCard[MAX_COLOR] = {0};
        memset (colorPart, -1, localNbElem * sizeof (int));
    	nbColors = create_bounded_color_part (colorPart, colorCard, elemToElem,
                                              localNbElem);
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
        int vecOffset  = create_coloring_permutation (colorPerm, colorPart, colorCard,
                                                      localNbElem, nbColors);
        tree.vecOffset = vecOffset + tree.firstElem;
		delete[] colorPart;

		// Apply local permutation on elemToNode & on the global element permutation
		DC_permute_int_2d_array (elemToNode, colorPerm, localNbElem, dimElem,
                                 tree.firstElem);
		merge_permutations (colorPerm, globalNbElem, localNbElem, tree.firstElem,
                            tree.lastSep);
		delete[] colorPerm;
	}
	else {
        #ifdef OMP
            #ifdef STATS
                #pragma omp task default (shared)
                leaves_coloring (*tree.right, nodeToElem, elemToNode, elemPerColor,
                                 colorPerLeaf, globalNbElem, dimElem);
                #pragma omp task default (shared)
                leaves_coloring (*tree.left, nodeToElem, elemToNode, elemPerColor,
                                 colorPerLeaf, globalNbElem, dimElem);
            #else
                #pragma omp task default (shared)
                leaves_coloring (*tree.right, nodeToElem, elemToNode, globalNbElem,
                                 dimElem);
                #pragma omp task default (shared)
                leaves_coloring (*tree.left, nodeToElem, elemToNode, globalNbElem,
                                 dimElem);
            #endif
            #pragma omp taskwait
        #else
		    cilk_spawn
            #ifdef STATS
                leaves_coloring (*tree.right, nodeToElem, elemToNode, elemPerColor,
                				 colorPerLeaf, globalNbElem, dimElem);
                leaves_coloring (*tree.left, nodeToElem, elemToNode, elemPerColor,
                                 colorPerLeaf, globalNbElem, dimElem);
            #else
                leaves_coloring (*tree.right, nodeToElem, elemToNode, globalNbElem,
                                 dimElem);
                leaves_coloring (*tree.left, nodeToElem, elemToNode, globalNbElem,
                                 dimElem);
            #endif
            cilk_sync;
        #endif

		if (tree.sep != nullptr) {
            #ifdef STATS
			    leaves_coloring (*tree.sep, nodeToElem, elemToNode, elemPerColor,
			    				 colorPerLeaf, globalNbElem, dimElem);
            #else
                leaves_coloring (*tree.sep, nodeToElem, elemToNode, globalNbElem,
                                 dimElem);
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
    DC_create_nodeToElem (nodeToElem, elemToNode, nbElem, dimElem, nbNodes);

    #ifdef STATS
        string fileName = "colorPerLeaf_" +
                          to_string ((long long)MAX_ELEM_PER_PART) + ".csv";
        ofstream colorPerLeaf (fileName, ios::out | ios::trunc);
        colorPerLeaf << "leafNb nbColors\n";
        int *elemPerColor = new int [MAX_ELEM_PER_PART] ();
        #ifdef OMP
            #pragma omp parallel
            #pragma omp single nowait
        #endif
        leaves_coloring (*treeHead, nodeToElem, elemToNode, elemPerColor,
                         colorPerLeaf, nbElem, dimElem);
        coloring_stat (elemPerColor, nbElem);
        delete[] elemPerColor;
        colorPerLeaf.close ();
    #else
        #ifdef OMP
            #pragma omp parallel
            #pragma omp single nowait
        #endif
        leaves_coloring (*treeHead, nodeToElem, elemToNode, nbElem, dimElem);
    #endif

    delete[] nodeToElem.value, delete[] nodeToElem.index;
}

#endif
