#ifndef TOOLS_H
#define TOOLS_H

#include <fstream>

#include "globals.h"
#include "coloring.h"

// Get CPU cycles
uint64_t DC_get_cycles ();

// Sort by ascending node values couple arrays using parallel quick sort
void quick_sort (couple_t *tab, int begin, int end);

/*****************************************************************************/
/***********                  Vectorization stats                  ***********/
/*****************************************************************************/

// Compute the number of colors used per leaf and the number of elements per color
void leaf_coloring_stat (ofstream &colorPerLeaf, int *elemPerColor, int *colorPart,
                         int nbElem, int nbColors);

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
void fill_dc_file_leaves (ofstream &dcFile, int curNode, int firstElem, int lastElem,
                          int LRS);

// Fill the nodes of the D&C tree dot file
void fill_dc_file_nodes (ofstream &dcFile, int curNode, int firstElem, int lastElem,
                         int nbSepElem);

// Close the D&C tree dot file
void close_dc_file (ofstream &dcFile);

// Initialize the D&C tree dot file with default layout
void init_dc_file (ofstream &dcFile, int nbPart);

#endif
