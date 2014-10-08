#ifndef COLORING_H
#define COLORING_H

#include "globals.h"
#include "matrix.h"

#define MAX_COLOR 128

// Assign a color to the elements of a given leaf & return the number of colors
int create_color_part (int *colorPart, int *colorCard, list_t *elemToElem, int nbElem);

// Create a coloring permutation for each leaf of the D&C tree
void leaves_coloring (tree_t &tree, index_t &nodeToElem, int *elemToNode,
#ifdef STATS
                      int *elemPerColor, ofstream &colorPerLeaf,
#endif
					  int globalNbElem);

// Coloring of the D&C tree
void coloring (int *elemToNode, int nbElem, int nbNodes);

#endif
