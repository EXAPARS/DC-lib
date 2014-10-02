#ifndef DC_ASSEMBLY_H
#define DC_ASSEMBLY_H

#include "globals.h"

// Follow the D&C tree to execute the assembly step in parallel
void dc_assembly (tree_t &tree, double *coord, double *nodeToNodeValue,
                  int *nodeToNodeRow, int *nodeToNodeColumn, int *elemToNode,
                  int *elemToEdge, int operatorID);

// Iterate over the colors & execute the assembly step on the elements of a same color
// in parallel
void coloring_assembly (double *coord, double *nodeToNodeValue, int *nodeToNodeRow,
                        int *nodeToNodeColumn, int *elemToNode, int *elemToEdge,
                        int operatorID);

// Call the appropriate function to perform the assembly step
void assembly (double *coord, double *nodeToNodeValue, int *nodeToNodeRow,
               int *nodeToNodeColumn, int *elemToNode, int *elemToEdge, int nbElem,
               int nbEdges, int operatorDim, int operatorID);

#endif
