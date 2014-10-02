#include <iostream>

#include "permutations.h"
#include "coloring.h"

// Fill the index of elements per color
void fill_color_index (int *colors, int *colorPart, int nbElem, int nbColors,
                       int offset)
{
    int *colorCount = new int [nbColors] ();

    // Count the number of elements of each color
    for (int i = 0; i < nbElem; i++) {
        colorCount[colorPart[i]]++;
    }
    // Create the coloring index
    colors[0] = offset;
    for (int i = 1; i <= nbColors; i++) {
        colors[i] = colors[i-1] + colorCount[i-1];
    }
    delete[] colorCount;
}

// Assign a color to the elements of a given leaf & return the number of colors
int create_color_part (int *colorPart, int *colorCard, list_t *elemToElem, int nbElem)
{
	__uint128_t *elemToColor = new __uint128_t [nbElem] ();
	int nbColors = 0;

    // For each element of local interval
	for (int i = 0; i < nbElem; i++) {
		__uint128_t mask = 1, neighborColor = 0;
		int color = 0;

        // Get the color of all neigbor elements
        for (int j = 0; j < elemToElem[i].size; j++) {
            int neighbor = elemToElem[i].list[j];
            neighborColor |= elemToColor[neighbor];
		}
        // Get the first free color (position of the first 0 bit)
		while (neighborColor & mask || colorCard[color] >= VEC_SIZE) {
		    neighborColor = neighborColor >> 1;
		    color++;
		}
		if (color >= MAX_COLOR) {
		    cerr << "Error: Not enough colors !\n";
		    exit (EXIT_FAILURE);
        }
        // Assign the first free color to current element
		elemToColor[i] = (mask << color);
        colorPart[i] = color;
        colorCard[color]++;

        // Compute the total number of colors
        if (color > nbColors) nbColors = color;
	}
    nbColors++;

	delete[] elemToColor;
    return nbColors;
}

void coloring_permutation ()
{
    // Apply element permutation function
    if (mpiRank == 0) {
        cout << "Applying element permutation...      ";
        t1 = MPI_Wtime ();
    }
    permute_int_2d_array (elemToNode, elemPerm, nbElem, DIM_ELEM, 0);
    delete[] elemPerm;
    if (mpiRank == 0) {
        t2 = MPI_Wtime ();
    	cout << "done  (" << t2 - t1 << " seconds)\n";
    }
}

// Mesh coloring of the whole mesh (elemToNode)
void coloring_creation (int *elemToNode, int nbElem, int nbNodes)
{
    if (mpiRank == 0) {
        cout << "Coloring of the mesh...              ";
        t1 = MPI_Wtime ();
    }
    elemPerm = new int [nbElem];

    // List the neighbor elements of each node
    index_t nodeToElem;
    nodeToElem.index = new int [nbNodes + 1];
    nodeToElem.value = new int [nbElem * DIM_ELEM];
    node_to_elem (nodeToElem, elemToNode, nbElem, nbNodes);

    // List the neighbor elements of each element
    list_t *elemToElem = new list_t [nbElem];
    elem_to_elem (elemToElem, nodeToElem, elemToNode, 0, nbElem-1);
    delete[] nodeToElem.value, delete[] nodeToElem.index;

    // Assign a color to each element
    delete[] elemToElem;

    int *colorPart = new int [nbElem];
    int colorCard[MAX_COLOR] = {0};
    nbTotalColors = create_color_part (colorPart, colorCard, elemToElem, nbElem);
    delete[] elemToElem;

    // Fill the index of elements per color
    colorToElem = new int [nbTotalColors+1];
    fill_color_index (colorToElem, colorPart, nbElem, nbTotalColors, 0);

    // Create a permutation array to sort the elements per color
    elemPerm = new int [nbElem];
    create_perm_array (elemPerm, colorPart, nbElem, MAX_COLOR);
    delete[] colorPart;

    if (mpiRank == 0) {
        t2 = MPI_Wtime ();
    	cout << "done  (" << t2 - t1 << " seconds)\n";
    }
}
