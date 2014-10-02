#ifndef COLORING_H
#define COLORING_H

#include "matrix.h"
#include "globals.h"

#define MAX_COLOR 128

// Fill the index of elements per color
void fill_color_index (int *colors, int *colorPart, int nbElem, int nbColors,
                       int offset);

// Assign a color to the elements of a given leaf & return the number of colors
int create_color_part (int *colorPart, int *colorCard, list_t *elemToElem, int nbElem);

// Mesh coloring of the whole mesh (elemToNode)
void coloring (int *elemToNode, int nbElem, int nbNodes);

#endif
