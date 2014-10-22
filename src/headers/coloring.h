#ifndef COLORING_H
#define COLORING_H

#include "globals.h"

#define NB_BLOCKS 32
#define BLOCK_SIZE sizeof (int)
#define MAX_COLOR (NB_BLOCKS * BLOCK_SIZE)

typedef struct {
    int elem, node;
} couple_t;
typedef struct {
    int *index, *value;
} index_t;
typedef struct {
    int list[MAX_ELEM_NEIGHBORS];
    int size;
} list_t;

// Assign a color to the elements of a given leaf & return the number of colors
int create_color_part (int *colorPart, int *colorCard, list_t *elemToElem, int nbElem);

// Construct element to element array from element to node and node to element
// Two elements are neighbors if they share a node
void elem_to_elem (list_t *elemToElem, index_t &nodeToElem, int *elemToNode,
                   int firstElem, int lastElem, int dimElem);

// Construct node to element structure from element to node
void node_to_elem (index_t &nodeToElem, int *elemToNode, int nbElem, int dimElem,
                   int nbNodes);

// Create a coloring permutation for each leaf of the D&C tree
void leaves_coloring (tree_t &tree, index_t &nodeToElem, int *elemToNode,
#ifdef STATS
                      int *elemPerColor, ofstream &colorPerLeaf,
#endif
					  int globalNbElem, int dimElem);

// Coloring of the D&C tree
void coloring (int *elemToNode, int nbElem, int dimElem, int nbNodes);

#endif
