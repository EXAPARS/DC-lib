#ifndef PARTITIONING_H
#define PARTITIONING_H

#ifdef STATS
#include <fstream>
#endif
#include "globals.h"

// Create a nodal graph from a tetrahedron mesh (created from METIS)
void mesh_to_nodal (int *graphIndex, int *graphValue, int *elemToNode,
                    int nbElem, int nbNodes);

// Create temporal sepToNode array containing separator elements indexed
// contiguously from 0 to nbSepElem and return the number of separator nodes
int create_sepToNode (int *sepToNode, int *elemToNode, int firstSepElem,
                      int lastSepElem);

// D&C partitioning of separators with more than MAX_ELEM_PER_PART elements
void sep_partitioning (tree_t &tree, int *elemToNode, int globalNbElem,
                       int firstSepElem, int lastSepElem
#ifdef STATS
					   , ofstream &dcFile, int curNode);
#else
                       );
#endif

// Divide & Conquer partitioning of elemToNode array
void partitioning (int *elemToNode, int nbElem, int nbNodes);

#endif
