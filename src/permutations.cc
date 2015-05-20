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

#include "permutations.h"

extern int *elemPerm, *nodePerm;

// Permute "tab" 2D array of double using node permutation
void DC_permute_double_2d_array (double *tab, int nbItem, int dimItem)
{
	char   *checkPerm = new char   [nbItem] ();
	double *tmpSrc    = new double [dimItem];
	double *tmpDest   = new double [dimItem];
	for (int i = 0; i < nbItem; i++) {
		if (checkPerm[i] == 1) continue;

		int init = i, src = i, dest;
		for (int j = 0; j < dimItem; j++) {
			tmpSrc[j] = tab[i*dimItem+j];
		}
		do {
			dest = nodePerm[src];
			for (int j = 0; j < dimItem; j++) {
				tmpDest[j] = tab[dest*dimItem+j];
				tab[dest*dimItem+j] = tmpSrc[j];
				tmpSrc[j] = tmpDest[j];
			}
			src = dest;
			checkPerm[src] = 1;
		}
		while (src != init);
	}
	delete[] tmpDest, delete[] tmpSrc, delete[] checkPerm;
}

// Permute "tab" 2D array of int using "perm"
void DC_permute_int_2d_array (int *tab, int *perm, int nbItem, int dimItem, int offset)
{
	char *checkPerm = new char [nbItem] ();
	int  *tmpSrc    = new int  [dimItem];
	int  *tmpDest   = new int  [dimItem];

    // If no permutation is given, default behavior is to use D&C elemPerm
    if (perm == nullptr) perm = elemPerm;

	for (int i = 0; i < nbItem; i++) {
		if (checkPerm[i] == 1) continue;

		int init = i, src = i, dest;
		for (int j = 0; j < dimItem; j++) {
			tmpSrc[j] = tab[(i+offset)*dimItem+j];
		}
		do {
			dest = perm[src];
			for (int j = 0; j < dimItem; j++) {
				tmpDest[j] = tab[(dest+offset)*dimItem+j];
				tab[(dest+offset)*dimItem+j] = tmpSrc[j];
				tmpSrc[j] = tmpDest[j];
			}
			src = dest;
			checkPerm[src] = 1;
		}
		while (src != init);
	}
	delete[] tmpDest, delete[] tmpSrc, delete[] checkPerm;
}

// Permute "tab" 1D array of int using node permutation
void DC_permute_int_1d_array (int *tab, int size)
{
	char *checkPerm = new char [size] ();
	for (int i = 0; i < size; i++) {
		if (checkPerm[i] == 1) continue;

		int init = i, src = i, dest;
		int tmpSrc = tab[i], tmpDest;
		checkPerm[i] = 1;

		do {
			dest      = nodePerm[src];
			tmpDest   = tab[dest];
			tab[dest] = tmpSrc;
			tmpSrc    = tmpDest;
			src       = dest;
			checkPerm[src] = 1;
		}
		while (src != init);
	}
	delete[] checkPerm;
}

// Renumber "tab" array of int using node permutation
void DC_renumber_int_array (int *tab, int size, bool isFortran)
{
	for (int i = 0; i < size; i++) {
		int tmp = tab[i];
		if (isFortran) tmp--;
		tab[i] = nodePerm[tmp];
		if (isFortran) tab[i]++;
	}
}

// Apply local element permutation to global element permutation
void merge_permutations (int *localElemPerm, int globalNbElem, int localNbElem,
						 int firstElem, int lastElem)
{
	int ctr = 0;
	for (int i = 0; i < globalNbElem; i++) {
		int dest = elemPerm[i];
		if (dest >= firstElem && dest <= lastElem) {
			elemPerm[i] = localElemPerm[dest-firstElem] + firstElem;
			ctr++;
		}
		if (ctr == localNbElem)	break;
	}
}

#ifdef DC_VEC
// Create coloring permutation array with full vectorial colors stored first &
// return the index of the last element in a full vectorial color
int create_coloring_permutation (int *perm, int *part, int *card, int size,
                                 int nbColors)
{
    int ptr = 0, lastFullColor;

    // Full colors
    for (int i = 0; i < nbColors; i++) {
        if (card[i] == VEC_SIZE) {
            for (int j = 0; j < size; j++) {
                if (part[j] == i) {
                    perm[j] = ptr;
                    ptr++;
                }
            }
        }
    }
    lastFullColor = ptr - 1;

    // Other colors
    for (int i = 0; i < nbColors; i++) {
        if (card[i] < VEC_SIZE) {
            for (int j = 0; j < size; j++) {
                if (part[j] == i) {
                    perm[j] = ptr;
                    ptr++;
                }
            }
        }
    }
    return lastFullColor;
}
#endif

// Create permutation array from partition array
void DC_create_permutation (int *perm, int *part, int size, int nbPart)
{
    int ptr = 0;
    for (int i = 0; i < nbPart; i++) {
        for (int j = 0; j < size; j++) {
            if (part[j] == i) {
                perm[j] = ptr;
                ptr++;
            }
        }
    }
}
