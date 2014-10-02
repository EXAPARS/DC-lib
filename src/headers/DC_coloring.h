#ifndef DC_COLORING_H
#define DC_COLORING_H

#include "matrix.h"
#include "globals.h"

// Create a coloring permutation for each leaf of the D&C tree
void leaves_coloring (tree_t &tree, index_t &nodeToElem, int *elemToNode,
#ifdef STATS
                      int *elemPerColor, ofstream &colorPerLeaf,
#endif
					  int globalNbElem);

// Coloring of the D&C tree
void dc_coloring (int *elemToNode, int nbElem, int nbNodes);

#endif
