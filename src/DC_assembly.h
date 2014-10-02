#include <cilk/cilk.h>

#include "assembly.h"
#include "dc_assembly.h"

// Follow the D&C tree to execute the assembly step in parallel
void dc_assembly (tree_t &tree, double *coord, double *nodeToNodeValue,
                  int *nodeToNodeRow, int *nodeToNodeColumn, int *elemToNode,
                  int *elemToEdge, int operatorID, int operatorDim)
{
    // If current node is a leaf, call the appropriate assembly function
    if (tree.left == NULL && tree.right == NULL) {

        // If leaf is not a separator, reset locally the CSR matrix
        if (tree.firstCSR != -1) {
            int firstCSR = tree.firstCSR * operatorDim;
            int lastCSR  = (tree.lastCSR + 1) * operatorDim - firstCSR;
            nodeToNodeValue[firstCSR:lastCSR] = 0;
        }

#ifdef HYBRID
        // Call assembly function for all full colors in vectorial and other colors
        // in sequential using laplacian operator
        if (operatorID == 0) {
            assembly_lap_vec (coord, nodeToNodeValue, nodeToNodeRow, nodeToNodeColumn,
                              elemToNode, elemToEdge, tree.firstElem,
                              tree.lastFullColor);
            assembly_lap_seq (coord, nodeToNodeValue, nodeToNodeRow, nodeToNodeColumn,
                              elemToNode, elemToEdge, tree.lastFullColor+1,
                              tree.lastElem);
        }
        // Using elasticity operator
        else {
            assembly_ela_vec (coord, nodeToNodeValue, nodeToNodeRow, nodeToNodeColumn,
                              elemToNode, elemToEdge, tree.firstElem,
                              tree.lastFullColor);
            assembly_ela_seq (coord, nodeToNodeValue, nodeToNodeRow, nodeToNodeColumn,
                              elemToNode, elemToEdge, tree.lastFullColor+1,
                              tree.lastElem);
        }
#else
        // Laplacian sequential assembly
        if (operatorID == 0) {
            assembly_lap_seq (coord, nodeToNodeValue, nodeToNodeRow, nodeToNodeColumn,
                              elemToNode, elemToEdge, tree.firstElem, tree.lastElem);
        }
        // Elasticity sequential assembly
        else {
            assembly_ela_seq (coord, nodeToNodeValue, nodeToNodeRow, nodeToNodeColumn,
                              elemToNode, elemToEdge, tree.firstElem, tree.lastElem);
        }
#endif
    }
    else {
        // Left & right recursion
        cilk_spawn
        dc_assembly (*tree.left, coord, nodeToNodeValue, nodeToNodeRow,
                     nodeToNodeColumn, elemToNode, elemToEdge, operatorID,
                     operatorDim);
        dc_assembly (*tree.right, coord, nodeToNodeValue, nodeToNodeRow,
                     nodeToNodeColumn, elemToNode, elemToEdge, operatorID,
                     operatorDim);

        // Synchronization
        cilk_sync;

        // Separator recursion, if it is not empty
        if (tree.sep != NULL) {
            dc_assembly (*tree.sep, coord, nodeToNodeValue, nodeToNodeRow,
                         nodeToNodeColumn, elemToNode, elemToEdge, operatorID,
                         operatorDim);
        }
    }
}

// Iterate over the colors & execute the assembly step on the elements of a same color
// in parallel
void coloring_assembly (double *coord, double *nodeToNodeValue, int *nodeToNodeRow,
                        int *nodeToNodeColumn, int *elemToNode, int *elemToEdge,
                        int operatorID)
{
    // For each color
    for (int color = 0; color < nbTotalColors; color++) {

        // Get the interval of elements of the current color
        int firstElem = colorToElem[color];
        int lastElem  = colorToElem[color+1] - 1;

        // Call assembly function using laplacian operator
        if (operatorID == 0) {
            assembly_lap_seq (coord, nodeToNodeValue, nodeToNodeRow, nodeToNodeColumn,
                              elemToNode, elemToEdge, firstElem, lastElem);
        }
        // Using elasticity operator
        else {
            assembly_ela_seq (coord, nodeToNodeValue, nodeToNodeRow, nodeToNodeColumn,
                              elemToNode, elemToEdge, firstElem, lastElem);
        }
    }
}

// Call the appropriate function to perform the assembly step
void assembly (double *coord, double *nodeToNodeValue, int *nodeToNodeRow,
               int *nodeToNodeColumn, int *elemToNode, int *elemToEdge, int nbElem,
               int nbEdges, int operatorDim, int operatorID)
{
#ifdef REF
    // Sequential reset of CSR matrix
    for (int i = 0; i < nbEdges * operatorDim; i++) {
        nodeToNodeValue[i] = 0;
    }
    // Laplacian sequential assembly
    if (operatorID == 0) {
        assembly_lap_seq (coord, nodeToNodeValue, nodeToNodeRow, nodeToNodeColumn,
                          elemToNode, elemToEdge, 0, nbElem-1);
    }
    // Elasticity sequential assembly
    else {
        assembly_ela_seq (coord, nodeToNodeValue, nodeToNodeRow, nodeToNodeColumn,
                          elemToNode, elemToEdge, 0, nbElem-1);
    }
#elif COLORING
    // Parallel reset of CSR matrix
    //#pragma omp parallel for
    cilk_for (int i = 0; i < nbEdges * operatorDim; i++) {
        nodeToNodeValue[i] = 0;
    }
    // Coloring parallel assembly
    coloring_assembly (coord, nodeToNodeValue, nodeToNodeRow, nodeToNodeColumn,
                       elemToNode, elemToEdge, operatorID);
#elif DC
    // D&C parallel assembly
    dc_assembly (*treeHead, coord, nodeToNodeValue, nodeToNodeRow, nodeToNodeColumn,
                 elemToNode, elemToEdge, operatorID, operatorDim);
#endif
}
