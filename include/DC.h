#ifndef DC_H
#define DC_H

#include <stdint.h>
#include <string>

#define MAX_ELEM_NEIGHBORS 250
#define MAX_ELEM_PER_PART 200

// D&C tree structure
typedef struct tree_s {
    int firstElem, lastElem, lastSep;
    int firstNode, lastNode;
    int firstEdge, lastEdge;
    int vecOffset;
    struct tree_s *left, *right, *sep;
} tree_t;

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

// Get CPU cycles
uint64_t DC_get_cycles ();

// Wrapper used to get the root of the D&C tree before calling the real tree traversal
void DC_tree_traversal (void (*userSeqFct) (void *, int, int),
                        void (*userVecFct) (void *, int, int),
                        void *userArgs, double *nodeToNodeValue, int operatorDim);

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
void DC_finalize_tree (int *nodeToNodeRow, int *elemToNode, int nbNodes);

// Create the D&C tree and the permutations
void DC_create_tree (int *elemToNode, int nbElem, int dimElem, int nbNodes);

#endif
