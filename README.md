DC-lib
======

Overview
--------

The DC-lib is a parallelization library based on the Divide & Conquer approach.
It is used to exploit efficiently shared memory parallelism in 3D unstuctured applications.
The DC-lib uses the Intel Cilk Plus tasks based runtime.
Two versions of the library are available:
- A pure D&C version.
- An hybrid version combining D&C with mesh coloring to enable vectorization at task level.

How to compile
--------------

DC-lib requires CMake 2.6 or later and METIS 5.1.0.
The build directory is "DC-lib/build".
The path to the METIS library can be se at the beginning of the iMake file.

A single command is required to compile a new binary:
  ./iMake [hybrid] [$VECTOR_LENGTH] [tree]

- The hybrid option is used to compile the library using the D&C + coloring version.
  Default uses the pure D&C version.

- The $VECTOR_LENGTH variable must be specified when using the hybrid version.
  It can be either SSE, AVX or MIC depending on the target architecture.

- The tree option is used to create a new D&C tree and new permutation functions.
  If not specified, the application will try to read the existing tree and permutations.
  The created tree and permutations are stored in "Mini-FEM/data/$USE_CASE/DC_tree".
  They are associated to the code version (D&C or D&C Hybrid), to the partition size,
  and to the number of MPI processes. The files name can be read this way:
    $VERSION_$PARTITION_SIZE_$NB_PROCESS_$PROCESS_RANK

If you want to build binaries for Xeon Phi, you need to use the MIC option whatever
the code version.

You can vary the size of the D&C partitions by modifying the "DC.h" header file. It is
located in "DC-lib/include/". The value corresponds to the "MAX_ELEM_PER_PART" define.

API
---

- Creation of the D&C tree and the permutation functions:

    void DC_create_tree (int *elemToNode, int nbElem, int dimElem, int nbNodes);

This function is required to create the D&C tree and the permutation functions.
Using this function requires the tree option at compile time. 

Inputs:
- int *elemToNode : Connectivity table (contains the list of nodes used per element)
- int nbElem      : Number of elements
- int dimElem     : Dimension of the elements
- int nbNodes     : Number of nodes

---

- Finalization of the D&C tree:

    void DC_finalize_tree (int *nodeToNodeRow, int *elemToNode, int nbNodes);

This function is used to compute the edges intervals for CSR reset.

Inputs:
- int *nodeToNodeRow : Row indexes of the CSR matrix
- int *elemToNode    : Connectivity table (contains the list of nodes used per element)
- int nbNodes        : Number of nodes

---

- Storing the D&C tree:

    void DC_store_tree (std::string &treePath, int nbElem, int nbNodes);

This function stores the current D&C tree and permutation functions to a binary file.

Inputs:
- std::string &treePath : Path to the destination of the binary file
- int nbElem            : Number of elements
- int nbNodes           : Number of nodes

---

- Reading the D&C tree:

    void DC_read_tree (std::string &treePath, int nbElem, int nbNodes);

This function reads previously created D&C trees and permutation functions depending of
the input parameters.

Inputs:
- std::string &treePath : Path to the targeted binary file
- int nbElem            : Number of elements
- int nbNodes           : Number of nodes

---

- Permutations:

    void DC_permute_double_2d_array (double *tab, int nbItem, int dimItem);

    void DC_permute_int_2d_array (int *tab, int *perm, int nbItem, int dimItem, int offset);

    void DC_permute_int_1d_array (int *tab, int size);

These functions permute a given 1D/2D array of int/double using node/perm permutation

Inputs:
- double/int *tab : Array that will be permuted
- int *perm       : Permutation array (default is node permutation)
- int nbItem/size : Number of elements/nodes in tab
- int dimItem     : Dimension of the elements/nodes of tab
- int offset      : Usefull to permute a subpart of tab (default is 0)

---

- Parallel assembly function:

    void DC_assembly (void (*userSeqFct) (void *, int, int),
                      void (*userVecFct) (void *, int, int),
                      void *userArgs, double *nodeToNodeValue, int operatorDim);

This function is used to parallelized user functions iterating over mesh elements.
It uses the precomputed D&C tree to split the element interval in many parallel tasks.
If the hybrid option is used at compile time, users can call both a function iterating
sequentially over the elements, and a function iterating over the elements with a stride
equal to the vector length specified at compile time.

Inputs:
- void (*userSeqFct) (void *, int, int) : Function pointer to the user function
                                          iterating sequentially over elements
- void (*userVecFct) (void *, int, int) : Function pointer to the user function
                                          iterating by vectors of size $VECTOR_LENGTH
                                          (default is NULL)
- void *userArgs                        : Structure containing all the arguments used
                                          in the seq/vec user functions
- double *nodeToNodeValue               : Values of the CSR matrix
