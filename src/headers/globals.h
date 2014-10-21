#ifndef GLOBALS_H
#define GLOBALS_H

#include <stdlib.h>
#include <string>

#define MAX_ELEM_NEIGHBORS 250
#define MAX_ELEM_PER_PART 200

using namespace std;

typedef struct tree_s {
    int firstElem, lastElem, lastSep;
    int firstCSR, lastCSR;
    int vecOffset;
    struct tree_s *left, *right, *sep;
} tree_t;

extern string meshName;
extern tree_t *treeHead;
extern int *nodePerm, *elemPerm;

#endif
