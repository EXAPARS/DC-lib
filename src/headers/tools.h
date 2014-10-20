#ifndef TOOLS_H
#define TOOLS_H

#include <fstream>
#include "globals.h"

// Get CPU cycles
uint64_t DC_get_cycles ();

// Sort by ascending node values couple arrays using parallel quick sort
void quick_sort (couple_t *tab, int begin, int end);

// Return the euclidean norm of given array
double compute_double_norm (double *tab, int size);

/*****************************************************************************/
/***********                    2D matrix plot                     ***********/
/*****************************************************************************/

// Determine for each tree node if it's a left, right or separator node
void rec_2d_map (tree_t &tree, int *checkNode, int *elemToNode, int LRS);

// Fill 2D map file for D&C version
void create_dc_2d_map (int *elemToNode, int *nodeToNodeRow,
                       int *nodeToNodeColumn, int nbNodes);

// Fill 2D map file for Ref version
extern "C"
void create_ref_2d_map_ (int *nodeToNodeRow, int *nodeToNodeColumn,
                         int *nbNodes);

/*****************************************************************************/
/***********                  Vectorization stats                  ***********/
/*****************************************************************************/

// Compute the number of colors used per leaf of the D&C tree
// and the number of elements per color
void leaf_coloring_stat (ofstream &colorPerLeaf, int *elemPerColor,
						 int *colorPart, int nbElem, int nbColors);

// Fill the elemToElem, colorPerLeaf & elemPerColor files
void coloring_stat (int *elemPerColor, int nbElem);

// Compute the number of elements per leaf of the D&C tree
void leaf_dc_stat (tree_t &tree, ofstream &elemPerLeaf);

// Fill the elemPerLeaf file
void dc_stat ();

/*****************************************************************************/
/***********                   D&C tree dot file                   ***********/
/*****************************************************************************/

// Fill the leaves of the D&C tree dot file
void fill_dc_file_leaves (ofstream &dcFile, int curNode, int firstElem,
						  int lastElem, int LRS);

// Fill the nodes of the D&C tree dot file
void fill_dc_file_nodes (ofstream &dcFile, int curNode, int firstElem,
                         int lastElem, int nbSepElem);

// Close the D&C tree dot file
void close_dc_file (ofstream &dcFile);

// Initialize the D&C tree dot file with default layout
void init_dc_file (ofstream &dcFile, int nbPart);

#endif
