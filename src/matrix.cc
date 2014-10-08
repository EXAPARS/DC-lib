#include <cilk/cilk.h>
#include "matrix.h"

// Sort by ascending node values couple_t arrays using parallel quick sort
void quick_sort (couple_t *tab, int begin, int end)
{
    int left = begin - 1, right = end + 1, pivot = tab[begin].node;

    if (begin < end) {
        while (1)
        {
            do right--; while(tab[right].node > pivot);
            do left++;  while(tab[left].node  < pivot);

            if (left < right) {
                couple_t tmp = tab[left];
                tab[left]    = tab[right];
                tab[right]   = tmp;
            }
            else break;
        }

        cilk_spawn
        quick_sort (tab, begin, right);
        quick_sort (tab, right+1, end);
        cilk_sync;
    }
}

// Create elem to edge array giving the index of each edge of each element
void elem_to_edge (int *nodeToNodeRow, int *nodeToNodeColumn, int *elemToNode,
                   int *elemToEdge, int nbElem)
{
    // For each element
    cilk_for (int i = 0; i < nbElem; i++) {
        int ctr = 0;
        // For each edge of current element
        for (int j = 0; j < DIM_ELEM; j++) {
            int node1 = elemToNode[i*DIM_ELEM+j] - 1;
            for (int k = 0; k < DIM_ELEM; k++) {
                int node2 = elemToNode[i*DIM_ELEM+k] - 1;
                // Get the index of current edge from nodeToNode
                for (int l = nodeToNodeRow[node1];
                         l < nodeToNodeRow[node1+1]; l++) {
                    if (nodeToNodeColumn[l] == (node2 + 1)) {
                        elemToEdge[i*VALUES_PER_ELEM+ctr] = l;
                        ctr++;
                        break;
                    }
                }
            }
        }
    }
}

// Construct node to node arrays from node to element and element to node
void node_to_node (int *nodeToNodeRow, int *nodeToNodeColumn,
                   index_t &nodeToElem, int *elemToNode, int nbNodes)
{
    int nodeToNodeCtr = 0;

	// For each node
    for (int i = 0; i < nbNodes; i++) {
        int nbNeighbors = 0,
            nbMaxNeighbors = (nodeToElem.index[i+1] - nodeToElem.index[i]) * DIM_ELEM;
        int *neighbors = new int [nbMaxNeighbors];
        nodeToNodeRow[i] = nodeToNodeCtr;
		// For each neighbor element of current node
        for (int j = nodeToElem.index[i]; j < nodeToElem.index[i+1]; j++) {
            int elemNeighbor = nodeToElem.value[j];
		    // For each node of neighbor element
            for (int k = 0; k < DIM_ELEM; k++) {
                int nodeNeighbor = elemToNode[elemNeighbor*DIM_ELEM+k];
                bool isNew = true;
                // Check if current node neighbor is already stored
                for (int l = 0; l < nbNeighbors; l++) {
                    if (nodeNeighbor == neighbors[l]) {
                        isNew = false;
                    }
                }
				// If current neighbor is not in the list, add it
                if (isNew) {
                    nodeToNodeColumn[nodeToNodeCtr] = nodeNeighbor;
                    neighbors[nbNeighbors] = nodeNeighbor;
                    nodeToNodeCtr++;
                    nbNeighbors++;
                }
            }
        }
        delete[] neighbors;
    }
    nodeToNodeRow[nbNodes] = nodeToNodeCtr;
}

// Construct element to element array from element to node and node to element
// Two elements are neighbors if they share a node
void elem_to_elem (list_t *elemToElem, index_t &nodeToElem, int *elemToNode,
                   int firstElem, int lastElem)
{
	// For each node of given element interval
    cilk_for (int i = firstElem; i <= lastElem; i++) {
        int nbNeighbors = 0;
        for (int j = 0; j < DIM_ELEM; j++) {
			int node = elemToNode[i*DIM_ELEM+j] - 1;
			// For each neighbor element of current node
            for (int k = nodeToElem.index[node];
                     k < nodeToElem.index[node+1]; k++) {
                int elemNeighbor = nodeToElem.value[k];
				// If neighbor element is in the element interval and is not
                // current element
                if (elemNeighbor >= firstElem && elemNeighbor <= lastElem &&
                    elemNeighbor != i) {
                    bool isNew = true;
					// If neighbor element is not already stored
                    for (int l = 0; l < nbNeighbors; l++) {
                        if ((elemNeighbor - firstElem) ==
                            elemToElem[i-firstElem].list[l]) {
							isNew = false;
                            break;
						}
					}
					// If neighbor is not in the list, add it
                    if (isNew) {
                       	elemToElem[i-firstElem].list[nbNeighbors] =
                            elemNeighbor - firstElem;
                        nbNeighbors++;
					}
				}
			}
		}
        elemToElem[i-firstElem].size = nbNeighbors;
	}
}

// Construct node to element structure from element to node
void node_to_elem (index_t &nodeToElem, int *elemToNode, int nbElem,
                   int nbNodes)
{
    couple_t *couple = new couple_t [nbElem * DIM_ELEM];
    cilk_for (int i = 0; i < nbElem; i++) {
        for (int j = 0; j < DIM_ELEM; j++) {
            couple[i*DIM_ELEM+j].elem = i;
            couple[i*DIM_ELEM+j].node = elemToNode[i*DIM_ELEM+j] - 1;
        }
    }
    quick_sort (couple, 0, nbElem * DIM_ELEM - 1);

    int ctr = 0;
    for (int i = 0; i < nbNodes; i++) {
        nodeToElem.index[i] = ctr;
        while (couple[ctr].node == i) {
            nodeToElem.value[ctr] = couple[ctr].elem;
            ctr++;
			if (ctr == nbElem * DIM_ELEM) break;
        }
    }
    nodeToElem.index[nbNodes] = ctr;
    delete[] couple;
}
