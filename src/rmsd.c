#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "common.h"

#include "utils.h"
#include "io.h"

int check_for_basic_matches(MolecularData* mol_data_1, MolecularData* mol_data_2) {
    // Check if atom counts are the same
    if (mol_data_1->num_atoms != mol_data_2->num_atoms) {
        fprintf(stderr, "Number of atoms in the two molecules are different.\n"); return 0;
    }

    // Check if atom types are the same
    if (!check_array_unordered_identity(mol_data_1->num_atoms, mol_data_1->atomic_numbers, mol_data_2->atomic_numbers)) {
        fprintf(stderr, "Atoms in the two molecules are different.\n"); return 0;
    }

    // Check if the layer data is the same
    if (!check_array_unordered_identity(mol_data_1->num_atoms, mol_data_1->layer_data, mol_data_2->layer_data)) {
        fprintf(stderr, "Atom descriptors used for matching in the two molecules are different.\n"); return 0;
    }

    return 1;
}

double rmsd(MolecularData* mol_data_1, MolecularData* mol_data_2, int* best_assignment) {
    // Check for basic matches
    if (!check_for_basic_matches(mol_data_1, mol_data_2)) { 
        fprintf(stderr, "Basic matches failed. Molecules are not isomorphic.\n"); return INFINITY; 
    }

    int num_atoms = mol_data_1->num_atoms;

    // Get candidate data
    int** candidate_data = get_candidate_data(num_atoms, mol_data_1->layer_data, mol_data_2->layer_data);
    if (!candidate_data) { fprintf(stderr, "Failed to get candidate data.\n"); return INFINITY; }

    double** candidate_squared_distances = get_candidate_squared_distances(num_atoms, mol_data_1->coordinates, mol_data_2->coordinates, candidate_data);
    if (!candidate_squared_distances) { 
        for (int i = 0; i < num_atoms; ++i) free(candidate_data[i]); free(candidate_data);
        fprintf(stderr, "Failed to get candidate squared distances.\n"); return INFINITY; 
    }

    const SearchData* data = init_search_data(num_atoms, candidate_data, candidate_squared_distances, mol_data_1->adjacency_list, mol_data_2->adjacency_list);
    
    if (!data) { 
        for (int i = 0; i < num_atoms; ++i) free(candidate_data[i]); free(candidate_data);
        for (int i = 0; i < num_atoms; ++i) free(candidate_squared_distances[i]); free(candidate_squared_distances);
        fprintf(stderr, "Failed to initialize recursive state.\n"); return INFINITY; 
    }

    // Also iterative function is available, but this one is more efficient
    double best_rmsd = search_assignment_recurse(data, best_assignment);
    
    for (int i = 0; i < num_atoms; ++i) free(candidate_data[i]); free(candidate_data);
    for (int i = 0; i < num_atoms; ++i) free(candidate_squared_distances[i]); free(candidate_squared_distances);

    return best_rmsd;
}
