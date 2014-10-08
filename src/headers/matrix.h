#ifndef MATRIX_H
#define MATRIX_H

#include "globals.h"

typedef struct {
    int elem, node;
} couple_t;
typedef struct {
    int *index, *value;
} index_t;
typedef struct {
    int list[MAX_ELEM_NEIGHBORS];
    int size;
} list_t;

// Sort by ascending node values couple_t arrays using parallel quick sort
void quick_sort (couple_t *tab, int begin, int end);

// Create elem to edge array giving the index of each edge of each element
void elem_to_edge (int *nodeToNodeRow, int *nodeToNodeColumn, int *elemToNode,
                   int *elemToEdge, int nbElem);

// Construct node to node arrays from node to element and element to node
void node_to_node (int *nodeToNodeRow, int *nodeToNodeColumn,
                   index_t &nodeToElem, int *elemToNode, int nbNodes);

// Construct element to element array from element to node and node to element
// Two elements are neighbors if they share a node
void elem_to_elem (list_t *elemToElem, index_t &nodeToElem, int *elemToNode,
                   int firstElem, int lastElem);

// Construct node to element structure from element to node
void node_to_elem (index_t &nodeToElem, int *elemToNode, int nbElem,
                   int nbNodes);

#endif
