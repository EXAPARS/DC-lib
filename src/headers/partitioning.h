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

#ifndef PARTITIONING_H
#define PARTITIONING_H

#ifdef STATS
    #include <fstream>
#endif
#include "DC.h"

// Create a nodal graph from a tetrahedron mesh (created from METIS)
void mesh_to_nodal (int *graphIndex, int *graphValue, int *elemToNode, int nbElem,
                    int dimElem, int nbNodes);

// Create temporal sepToNode array containing separator elements indexed
// contiguously from 0 to nbSepElem and return the number of separator nodes
int create_sepToNode (int *sepToNode, int *elemToNode, int firstSepElem,
                      int lastSepElem, int dimElem);

// D&C partitioning of separators with more than MAX_ELEM_PER_PART elements
void sep_partitioning (tree_t &tree, int *elemToNode, int globalNbElem, int dimElem,
                       int firstSepElem, int lastSepElem, int firstNode, int lastNode,
                       int nbIntf, int *intfIndex, int *intfNodes,
#ifdef STATS
                       int curNode, ofstream &dcFile);
#else
                       int curNode);
#endif

// Divide & Conquer partitioning
void partitioning (int *elemToNode, int nbElem, int dimElem, int nbNodes,
                   int nbIntf, int *intfIndex, int *intfNodes, int rank);

int intf_partitioning (double *coord, int *elemToNode, int *intfIndex, int *intfNodes,
                       int nbElem, int dimElem, int nbNodes, int dimNode, int nbIntf);

#endif
