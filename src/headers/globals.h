#ifndef GLOBALS_H
#define GLOBALS_H

#include <stdlib.h>
#include <string>

#define DIM_ELEM 4
#define DIM_NODE 3
#define VALUES_PER_ELEM 16
#define MAX_ELEM_NEIGHBORS 250
#define MAX_ELEM_PER_PART 200
//extern int MAX_ELEM_PER_PART;

using namespace std;

typedef struct tree_s {
    int firstElem, lastElem, lastSep;
    int firstCSR, lastCSR;
    int lastFullColor;
    struct tree_s *left, *right, *sep;
} tree_t;

extern string meshName, operatorName;
extern tree_t *treeHead;
extern int *nodePerm, *elemPerm, *colorToElem;
extern int nbTotalColors;

#endif
