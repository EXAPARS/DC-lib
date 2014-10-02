// Compute edge intervals for CSR reset and store the D&C tree and the permutations
void DC_finalize ()
{
    #ifdef CREATE_PERM_AND_TREE
        // Compute the interval of edges of each leaf of the D&C tree
        if (mpiRank == 0) {
            cout << "Computing edge intervals...          ";
    	    t1 = MPI_Wtime ();
        }
        edge_intervals (*treeHead, nodeToNodeRow, elemToNode, nbNodes);
        if (mpiRank == 0) {
            t2 = MPI_Wtime ();
            cout << "done  (" << t2 - t1 << " seconds)\n";
        }

        // Store the D&C tree & the permutation functions
    	if (mpiRank == 0) {
            cout << "Storing D&C tree and permutations... ";
            t1 = MPI_Wtime ();
        }
    	store_perm_and_tree (nbElem, nbNodes, nbBlocks, mpiRank);
        if (mpiRank == 0) {
            t2 = MPI_Wtime ();
        	cout << "done  (" << t2 - t1 << " seconds)\n";
        }
        delete[] nodePerm, delete[] elemPerm; 
    #endif
}

// Apply the permutations
void DC_permutation ()
{
    if (mpiRank == 0) {
        cout << "Applying permutations...             ";
        t1 = MPI_Wtime ();
    }
    apply_transformations (coord, elemToNode, intfNodes, dispList, boundNodesCode,
                           nbElem, nbNodes, nbIntfNodes, nbDispNodes);
    if (mpiRank == 0) {
        t2 = MPI_Wtime ();
    	cout << "done  (" << t2 - t1 << " seconds)\n";
    }
}

// Create or load from file the D&C tree and the permutations
void DC_creation ()
{
    // Allocate the D&C tree & the permutation functions
    treeHead = new tree_t;
    elemPerm = new int [nbElem];
    nodePerm = new int [nbNodes];

    // Create the D&C tree & the permutation functions
    #ifdef CREATE_PERM_AND_TREE
        if (mpiRank == 0) {
            cout << "D&C partitioning...                  ";
            t1 = MPI_Wtime ();
        }
        partitioning (elemToNode, nbElem, nbNodes);
        if (mpiRank == 0) {
            t2 = MPI_Wtime ();
        	cout << "done  (" << t2 - t1 << " seconds)\n";
        }

        // Hybrid version with coloring of the leaves of the D&C tree
        #ifdef HYBRID
        	if (mpiRank == 0) {
            	cout << "Coloring D&C subdomains...           ";
                t1 = MPI_Wtime ();
            }
        	dc_coloring (elemToNode, nbElem, nbNodes);
        	if (mpiRank == 0) {
                t2 = MPI_Wtime ();
            	cout << "done  (" << t2-t1 << " seconds)\n";
            }
        #endif

    // Read the D&C tree & the permutation functions
    #else
    	if (mpiRank == 0) {
            cout << "Reading D&C tree and permutations... ";
            t1 = MPI_Wtime ();
        }
    	read_perm_and_tree (nbElem, nbNodes, nbBlocks, mpiRank);
        if (mpiRank == 0) {
            t2 = MPI_Wtime ();
        	cout << "done  (" << t2 - t1 << " seconds)\n";
        }
    #endif
}
