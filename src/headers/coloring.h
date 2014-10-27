#ifndef COLORING_H
#define COLORING_H

#include "DC_lib.h"

#define NB_BLOCKS 32
#define BLOCK_SIZE sizeof (int)
#define MAX_COLOR (NB_BLOCKS * BLOCK_SIZE)

// Create element to element array from element to node and node to element
// Two elements are neighbors if they share a node
void DC_create_elemToElem (list_t *elemToElem, index_t &nodeToElem, int *elemToNode,
                          int firstElem, int lastElem, int dimElem);

// Create node to element structure from element to node
void DC_create_nodeToElem (index_t &nodeToElem, int *elemToNode, int nbElem,
                          int dimElem, int nbNodes);

// Assign a color to the elements of a given leaf following the bounded colors strategy
// & return the number of colors
int create_bounded_color_part (int *colorPart, int *colorCard, list_t *elemToElem,
                               int nbElem);

// Create a coloring permutation for each leaf of the D&C tree
void leaves_coloring (tree_t &tree, index_t &nodeToElem, int *elemToNode,
#ifdef STATS
                      int *elemPerColor, ofstream &colorPerLeaf,
#endif
					  int globalNbElem, int dimElem);

// Coloring of the D&C tree
void coloring (int *elemToNode, int nbElem, int dimElem, int nbNodes);

#endif
