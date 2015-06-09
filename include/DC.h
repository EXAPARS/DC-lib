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

#ifndef DC_H
#define DC_H

#include <stdint.h>
#include <string>

#define MAX_ELEM_NEIGHBORS 250
#define MAX_ELEM_PER_PART 200

// D&C tree structure
typedef struct tree_s {
    int *ownedNodes;
    int nbOwnedNodes,
        firstElem, lastElem, lastSep,
        firstNode, lastNode,
        firstEdge, lastEdge;
    #ifdef DC_VEC
        int vecOffset;
    #endif
    bool isSep;
    struct tree_s *left, *right, *sep;
} tree_t;

// D&C arguments structure
typedef struct DCargs_s {
    int *ownedNodes;
    int nbOwnedNodes,
        firstElem, lastElem,
        firstNode, lastNode,
        firstEdge, lastEdge,
        isSep;
} DCargs_t;

typedef struct {
    int *index, *value;
} index_t;
typedef struct {
    int list[MAX_ELEM_NEIGHBORS];
    int size;
} list_t;
typedef struct {
    int elem, node;
} couple_t;

// Compute the average time/cycles of each (start, stop) intervals until reset call
class DC_timer
{
    public:
        // Constructor
        DC_timer ();

        // Get time of day
        double get_avg_time ();
        void reset_time ();
        void stop_time ();
        void start_time ();

        // RDTSC
        uint64_t get_avg_cycles ();
        void reset_cycles ();
        void stop_cycles ();
        void start_cycles ();

    private:
        // Get time of day
        double startTime, avgTime;
        int timeCtr;

        // RDTSC
        uint64_t startCycles, avgCycles;
        int cyclesCtr;
};

// Get time of day
double DC_get_time ();

// RDTSC
uint64_t DC_get_cycles ();

// Wrapper used to get the root of the D&C tree before calling the real tree traversal
void DC_tree_traversal (void (*userSeqFct) (void *, DCargs_t *),
                        void (*userVecFct) (void *, DCargs_t *), void *userArgs);

// Create element to element array from element to node and node to element
// Two elements are neighbors if they share a node
void DC_create_elemToElem (list_t *elemToElem, index_t &nodeToElem, int *elemToNode,
                           int firstElem, int lastElem, int dimElem);

// Create node to element structure from element to node
void DC_create_nodeToElem (index_t &nodeToElem, int *elemToNode, int nbElem,
                           int dimElem, int nbNodes);

// Permute "tab" 2D array of double using node permutation
void DC_permute_double_2d_array (double *tab, int nbItem, int dimItem);

// Permute "tab" 2D array of int using "perm"
void DC_permute_int_2d_array (int *tab, int *perm, int nbItem, int dimItem,
                              int offset);

// Permute "tab" 1D array of int using node permutation
void DC_permute_int_1d_array (int *tab, int size);

// Renumber "tab" array of int using node permutation
void DC_renumber_int_array (int *tab, int size, bool isFortran);

// Create permutation array from partition array
void DC_create_permutation (int *perm, int *part, int size, int nbPart);

// Read the D&C tree and the permutation functions
void DC_read_tree (std::string &treePath, int nbElem, int nbNodes);

// Store the D&C tree and the permutation functions to a binary file
void DC_store_tree (std::string &treePath, int nbElem, int nbNodes);

// Wrapper used to get the root of the D&C tree before computing the edge intervals
// for CSR reset
void DC_finalize_tree (int *nodeToNodeRow, int *elemToNode);

// Create the D&C tree and the permutations
void DC_create_tree (int *elemToNode, int nbElem, int dimElem, int nbNodes);

#endif
