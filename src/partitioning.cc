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

#ifdef TREE_CREATION

#include <cstring>
#include <cmath>
#ifdef CILK
    #include <cilk/cilk.h>
#endif
#include <pthread.h>
#include <metis.h>

#include "tools.h"
#include "permutations.h"
#include "tree_creation.h"
#include "partitioning.h"

extern tree_t *treeHead;
extern int *elemPerm, *nodePerm, *nodeOwner;

pthread_mutex_t metisMutex = PTHREAD_MUTEX_INITIALIZER;

// Create a nodal graph from a tetrahedron mesh (created from METIS)
void mesh_to_nodal (int *graphIndex, int *graphValue, int *elemToNode, int nbElem,
                    int dimElem, int nbNodes)
{
    int nEdges, *nPtr, *nInd, *marker;

    nPtr = new int [nbNodes + 1] ();
    for (int i = 0; i < dimElem * nbElem; i++) {
        nPtr[elemToNode[i]]++;
    }
    for (int i = 1; i < nbNodes; i++) {
        nPtr[i] += nPtr[i-1];
    }
    for (int i = nbNodes; i > 0; i--) {
        nPtr[i] = nPtr[i-1];
    }
    nPtr[0] = 0;

    nInd = new int [nPtr[nbNodes]];
    for (int i = 0; i < nbElem; i++) {
        for (int j = 0; j < dimElem; j++) {
            nInd[nPtr[elemToNode[i*dimElem+j]]++] = i;
        }
    }
    for (int i = nbNodes; i > 0; i--) {
        nPtr[i] = nPtr[i-1];
    }
    nPtr[0] = 0;

    marker = new int [nbNodes];
    memset (marker, -1, nbNodes * sizeof (int));
    nEdges = graphIndex[0] = 0;
    for (int i = 0; i < nbNodes; i++) {
        marker[i] = i;
        for (int j = nPtr[i]; j < nPtr[i+1]; j++) {
            int jj = dimElem * nInd[j];
            for (int k = 0; k < dimElem; k++) {
                int kk = elemToNode[jj];
                if (marker[kk] != i) {
                    marker[kk] = i;
                    graphValue[nEdges++] = kk;
                }
                jj++;
            }
        }
        graphIndex[i+1] = nEdges;
    }
    delete[] marker, delete[] nInd, delete[] nPtr;
}

// Create local elemToNode array containing elements indexed contiguously from 0 to
// nbElem and return the number of nodes accessed
int create_local_elemToNode (int *localElemToNode, int *elemToNode, int firstElem,
                             int lastElem, int dimElem)
{
    int nbNodes = 0;
    int *tmp = new int [(lastElem - firstElem + 1) * dimElem];

    for (int i = firstElem * dimElem, j = 0; i < (lastElem+1)*dimElem; i++, j++) {
        int newNode = 0, oldNode = elemToNode[i];
        bool isNew = true;
        for (newNode = 0; newNode < nbNodes; newNode++) {
            if (oldNode == tmp[newNode]) {
                isNew = false;
                break;
            }
        }
        if (isNew) {
            tmp[nbNodes] = oldNode;
            nbNodes++;
        }
        localElemToNode[j] = newNode;
    }
    delete[] tmp;
    return nbNodes;
}

// D&C partitioning of separators with more than MAX_ELEM_PER_PART elements
void sep_partitioning (tree_t &tree, int *elemToNode, int *intfIndex, int *intfNodes,
                       int globalNbElem, int dimElem, int firstSepElem,
                       int lastSepElem, int firstNode, int lastNode, int nbIntf,
                       int nbBlocks, int curNode, int commDepth, int curDepth
#ifdef STATS
                       , ofstream &dcFile)
#else
                       )
#endif
{
    // If there is not enough element in the separator
    int nbSepElem = lastSepElem - firstSepElem + 1;
    int nbSepPart = ceil (nbSepElem / (double)MAX_ELEM_PER_PART);
    if (nbSepPart < 2 || nbSepElem <= MAX_ELEM_PER_PART) {

        // Initialize the leaf
        bool hasIntfNode = false;
        init_dc_tree (tree, elemToNode, intfIndex, intfNodes, firstSepElem,
                      lastSepElem, 0, dimElem, firstNode, lastNode, nbIntf, nbBlocks,
                      commDepth, curDepth, true, true, &hasIntfNode);

        // Set the last updater of each node
        for (int i = firstSepElem * dimElem; i < (lastSepElem+1) * dimElem; i++){
            int node = nodePerm[elemToNode[i]];
            nodeOwner[node] = curNode;
        }

        #ifdef STATS
            fill_dc_file_leaves (dcFile, curNode, firstSepElem, lastSepElem, 3,
                                 hasIntfNode);
            count_intf_stats (hasIntfNode);
        #endif

        // End of recursion
        return;
    }

    // Create temporal elemToNode containing the separator elements
    int *sepToNode = new int [nbSepElem * dimElem];
    int nbSepNodes = create_local_elemToNode (sepToNode, elemToNode, firstSepElem,
                                              lastSepElem, dimElem);

    // Configure METIS & compute the node partitioning of the separators
    int constraint = 1, objVal;
    int *graphIndex = new int [nbSepNodes + 1];
    int *graphValue = new int [nbSepNodes * 15];
    int *nodePart   = new int [nbSepNodes];
    mesh_to_nodal (graphIndex, graphValue, sepToNode, nbSepElem, dimElem, nbSepNodes);

    // Execution is correct without mutex although cilkscreen detects many race
    // conditions. Check if the problem is solved with future version of METIS (5.0)
    pthread_mutex_lock (&metisMutex);
    METIS_PartGraphRecursive (&nbSepNodes, &constraint, graphIndex, graphValue,
                              nullptr, nullptr, nullptr, &nbSepPart, nullptr, nullptr,
                              nullptr, &objVal, nodePart);
    pthread_mutex_unlock (&metisMutex);
    delete[] graphValue, delete[] graphIndex;

    // Create the separator D&C tree
    tree_creation (tree, elemToNode, sepToNode, nodePart, nullptr, intfIndex,
                   intfNodes, globalNbElem, dimElem, 0, nbSepPart-1, firstSepElem,
                   lastSepElem, firstNode, lastNode, 0, nbIntf, nbBlocks, curNode,
    #ifdef STATS
                   commDepth, curDepth, true, dcFile, -1);
    #else
                   commDepth, curDepth, true);
    #endif
    delete[] nodePart, delete[] sepToNode;
}

// Divide & Conquer partitioning
void partitioning (int *elemToNode, int *intfIndex, int *intfNodes, int nbElem,
                   int dimElem, int nbNodes, int nbIntf, int nbBlocks, int rank)
{
    // Fortran to C elemToNode conversion
    #ifdef OMP
        #pragma omp parallel for
        for (int i = 0; i < nbElem * dimElem; i++) {
    #elif CILK
        cilk_for (int i = 0; i < nbElem * dimElem; i++) {
    #endif
        elemToNode[i]--;
    }

    // Configure METIS & compute the node partitioning of the mesh
    int nbPart = ceil (nbElem / (double)MAX_ELEM_PER_PART);
    int commDepth = ceil ((double)log2 (nbPart) / 2.);
	int constraint = 1, objVal;
    int *graphIndex = new int [nbNodes + 1];
    int *graphValue = new int [nbNodes * 15];
	int *nodePart   = new int [nbNodes];
    mesh_to_nodal (graphIndex, graphValue, elemToNode, nbElem, dimElem, nbNodes);
    METIS_PartGraphRecursive (&nbNodes, &constraint, graphIndex, graphValue,
                              nullptr, nullptr, nullptr, &nbPart, nullptr, nullptr,
                              nullptr, &objVal, nodePart);
    delete[] graphValue, delete[] graphIndex;

    // Create node permutation from node partition
    DC_create_permutation (nodePerm, nodePart, nbNodes, nbPart);

    // Compute the number of nodes per partition
    int *nodePartSize = new int [nbPart] ();
    for (int i = 0; i < nbNodes; i++) {
        nodePartSize[nodePart[i]]++;
    }

    // Create D&C tree dot file
    #ifdef STATS
        string fileName = "dcTree_" + to_string ((long long)rank) + "_" +
                           to_string ((long long)MAX_ELEM_PER_PART) + ".dot";
        ofstream dcFile (fileName, ios::out | ios::trunc);
        init_dc_file (dcFile, nbPart);
    #endif

	// Create D&C tree
	#ifdef OMP
		#pragma omp parallel
		#pragma omp single nowait
    #endif
    tree_creation (*treeHead, elemToNode, nullptr, nodePart, nodePartSize, intfIndex,
                   intfNodes, nbElem, dimElem, 0, nbPart-1, 0, nbElem-1, 0, nbNodes-1,
    #ifdef STATS
                   0, nbIntf, nbBlocks, 0, commDepth, 0, false, dcFile, -1);
    #else
                   0, nbIntf, nbBlocks, 0, commDepth, 0, false);
    #endif
    delete[] nodePartSize, delete[] nodePart;

    #ifdef STATS
    	close_dc_file (dcFile);
        store_intf_stats (nbElem, rank);
    	dc_stat ();
    #endif

	// C to Fortran elemToNode conversion
    #ifdef OMP
        #pragma omp parallel for
	    for (int i = 0; i < nbElem * dimElem; i++) {
    #elif CILK
    	cilk_for (int i = 0; i < nbElem * dimElem; i++) {
    #endif	
	    elemToNode[i]++;
	}
}
/*
// Store the elements & the nodes on the interface at the beginning and create a
// permutation array for each of them
int intf_partitioning (double *coord, int *elemToNode, int *intfIndex, int *intfNodes,
                       int nbElem, int dimElem, int nbNodes, int dimNode, int nbIntf)
{
    // Store the elements on the interfaces at the beginning
    int *elemIntfPart = new int [nbElem];
    int nbIntfElem = 0;
    cilk_for (int i = 0; i < nbElem; i++) elemIntfPart[i] = 1;
    for (int i = 0; i < nbIntf; i++) {
        for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
            int intfNode = intfNodes[j];
            for (int k = 0; k < nbElem; k++) {
                for (int l = 0; l < dimElem; l++) {
                    if (elemToNode[k*dimElem+l] == intfNode) {
                        elemIntfPart[k] = 0;
                        nbIntfElem++;
                        break;
                    }
                }
            }
        }
    }
    DC_create_permutation (elemPerm, elemIntfPart, nbElem, 2);
    DC_permute_int_2d_array (elemToNode, elemPerm, nbElem, dimElem, 0);
    delete[] elemIntfPart;

    // Store the nodes on the interfaces at the beginning
    int *nodeIntfPart = new int [nbNodes];
    cilk_for (int i = 0; i < nbNodes; i++) nodeIntfPart[i] = 1;
    for (int i = 0; i < nbIntf; i++) {
        for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
            int intfNode = intfNodes[j] - 1;
            nodeIntfPart[intfNode] = 0;
        }
    }
    DC_create_permutation (nodePerm, nodeIntfPart, nbNodes, 2);
    DC_permute_double_2d_array (coord, nbNodes, dimNode);
    delete[] nodeIntfPart;

    // Return the number of elements on the interface
    return nbIntfElem;
}
*/
#endif
