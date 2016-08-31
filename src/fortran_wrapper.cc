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

#include <string>
#include "DC_f.h"
/*
// Get the list of user arguments
void dc_get_args_ (uint64_t *argsList, int *size, ...)
{
    va_list userArgs;
    va_start (userArgs, *size);

    for (int i = 0; i < (*size); i++) {
        uint64_t **tmpArg = (uint64_t**)va_arg (userArgs, void*);

fprintf (stderr, "C recoit l'adr %lu\n", tmpArg);

        *tmpArg = (uint64_t*)argsList[i];

fprintf (stderr, "C envoie l'arg %lu\n", *tmpArg);

    }

    va_end (userArgs);
}

// Set the list of user arguments
void dc_set_args_ (uint64_t *argsList, int *size, ...)
{
    va_list userArgs;
    va_start (userArgs, *size);

    for (int i = 0; i < (*size); i++) {
        argsList[i] = (uint64_t)va_arg (userArgs, void*);

fprintf (stderr, "C genere l'arg %lu\n", argsList[i]);

    }

    va_end (userArgs);
}

// Initialize the list of user arguments
uint64_t* dc_init_args_list_ (int *size)
{
    uint64_t *argsList = new uint64_t [*size];
    return argsList;
}
*/
// Get CPU cycles
uint64_t dc_get_cycles_ ()
{
    return DC_get_cycles ();
}

// Get time of day
double dc_get_time_ ()
{
    return DC_get_time ();
}

// Wrapper used to get the root of the D&C tree before calling the real tree traversal
void dc_tree_traversal_ (void (*userSeqFct)  (void *, DCargs_t *),
                         void (*userVecFct)  (void *, DCargs_t *),
                         void (*userCommFct) (void *, DCcommArgs_t *),
                         void *userArgs, void *userCommArgs)
{
    DC_tree_traversal (userSeqFct, userVecFct, userCommFct, userArgs, userCommArgs);
}

// Permute "tab" 2D array of double using node permutation
void dc_permute_double_2d_array_ (double *tab, int *nbItem, int *dimItem)
{
    DC_permute_double_2d_array (tab, *nbItem, *dimItem);
}

// Permute "tab" 2D array of int using "perm"
void dc_permute_int_2d_array_ (int *tab, int *nbItem, int *dimItem, int *offset)
{
    DC_permute_int_2d_array (tab, nullptr, *nbItem, *dimItem, *offset);
}

// Permute "tab" 1D array of int using node permutation
void dc_permute_int_1d_array_ (int *tab, int *size)
{
    DC_permute_int_1d_array (tab, *size);
}

// Renumber "tab" array of int using node permutation
void dc_renumber_int_array_ (int *tab, int *size)
{
    DC_renumber_int_array (tab, *size, true);
}

// Create permutation array from partition array
void dc_create_permutation_ (int *perm, int *part, int *size, int *nbPart)
{
    DC_create_permutation (perm, part, *size, *nbPart);
}

// Read the D&C tree and the permutation functions
void dc_read_tree_ (char *treePath, int *nbElem, int *nbNodes, int *nbIntf,
                    int *nbMaxComm)
{
    std::string treePath_c (treePath);
    DC_read_tree (treePath_c, *nbElem, *nbNodes, *nbIntf, nbMaxComm);
}

// Store the D&C tree and the permutation functions to a binary file
void dc_store_tree_ (char *treePath, int *nbElem, int *nbNodes, int *nbIntf,
                     int *nbMaxComm)
{
    std::string treePath_c (treePath);
    DC_store_tree (treePath_c, *nbElem, *nbNodes, *nbIntf, *nbMaxComm);
}

#ifdef TREE_CREATION

// Wrapper used to get the root of the D&C tree before calling the real tree finalize
void dc_finalize_tree_ (int *nodeToNodeRow, int *elemToNode, int *intfIndex,
                        int *intfNodes, int *intfDstIndex, int *nbDCcomm,
                        int *nbElem, int *dimElem, int *nbBlocks, int *nbIntf,
                        int *rank)
{
    DC_finalize_tree (nodeToNodeRow, elemToNode, intfIndex, intfNodes, intfDstIndex,
                      nbDCcomm, *nbElem, *dimElem, *nbBlocks, *nbIntf, *rank);
}

// Create the D&C tree and the permutations
void dc_create_tree_ (int *elemToNode, int *nbElem, int *dimElem, int *nbNodes)
{
    DC_create_tree (elemToNode, *nbElem, *dimElem, *nbNodes);
}

#endif
