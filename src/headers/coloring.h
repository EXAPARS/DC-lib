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

#ifndef COLORING_H
#define COLORING_H

#include "DC.h"

#define NB_BLOCKS 32
#define BLOCK_SIZE 32
#define MAX_COLOR (NB_BLOCKS * BLOCK_SIZE)

// Create element to element array from element to node and node to element
// Two elements are neighbors if they share a node
void DC_create_elemToElem (list_t *elemToElem, index_t &nodeToElem, int *elemToNode,
                          int firstElem, int lastElem, int dimElem);

// Create node to element structure from element to node
void DC_create_nodeToElem (index_t &nodeToElem, int *elemToNode, int nbElem,
                          int dimElem, int nbNodes);

#ifdef DC_VEC

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
#endif
