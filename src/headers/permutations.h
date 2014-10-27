#ifndef PERMUTATIONS_H
#define PERMUTATIONS_H

// Permute "tab" 2D array of double using node permutation
void DC_permute_double_2d_array (double *tab, int nbItem, int dimItem);

// Permute "tab" 2D array of int using "perm"
void DC_permute_int_2d_array (int *tab, int *perm, int nbItem, int dimItem,
                              int offset);

// Permute "tab" 1D array of int using node permutation
void DC_permute_int_1d_array (int *tab, int size);

// Renumber "tab" array of int using node permutation
void DC_renumber_int_array (int *tab, int size, bool isFortran);

// Apply local element permutation to global element permutation
void merge_permutations (int *localElemPerm, int globalNbElem, int localNbElem,
						 int firstElem, int lastElem);

#ifdef HYBRID
    // Create coloring permutation array with full vectorial colors stored first &
    // return the index of the last element in a full vectorial color
    int create_coloring_permutation (int *perm, int *part, int *card, int size,
                                     int nbColors);
#endif

// Create permutation array from partition array
void DC_create_permutation (int *perm, int *part, int size, int nbPart);

#endif
