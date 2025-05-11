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


void rmsd_naive_rec_helper(MolecularData* mol_data_1, MolecularData* mol_data_2, 
                             int** candidate_data, double** candidate_squared_distances,
                             int* current_assignment, double current_sum, int minimize,
                             int* best_assignment, double* best_sum, 
                             int current_group, int group_count, int** groups) {

    if (current_group == group_count) {
        if (minimize)
            current_sum = kabsch_squared_dist_sum(mol_data_1->num_atoms, mol_data_1->coordinates, mol_data_2->coordinates, current_assignment);

        if (current_sum < *best_sum) {
            *best_sum = current_sum;
            memcpy(best_assignment, current_assignment, mol_data_1->num_atoms * sizeof(int));
        }
        return;
    }

    int group_size = groups[current_group][0];
    int* mol1_indices = groups[current_group] + 1;
    int* mol2_indices = candidate_data[mol1_indices[0]] + 1;

    int* permutation = (int*)malloc(group_size * sizeof(int));
    for (int i = 0; i < group_size; i++) {
        permutation[i] = i;
    }

    int** adjacency_list_1 = mol_data_1->adjacency_list;
    int** adjacency_list_2 = mol_data_2->adjacency_list;

    // Generate all permutations of the group
    do {
        for (int i = 0; i < group_size; i++) {
            int id = mol1_indices[i];
            current_assignment[id] = -1;
        }

        int is_permutation_valid = 1;

        double group_sum = 0.0;
        for (int i = 0; i < group_size; i++) {
            int id = mol1_indices[i];
            int mapped_choice = mol2_indices[permutation[i]];

            // Check if the current assignment is valid
            int is_mapping_valid = 1;
            for (int j = 1; j <= adjacency_list_1[id][0]; j++) {
                int neighbor_id_1 = get_bonded_atom_id(adjacency_list_1[id][j]);
                int bond_type_mask = get_bond_type(adjacency_list_1[id][j]) << BOND_TYPE_SHIFT;
                int mapped_neighbor_id = current_assignment[neighbor_id_1];

                if (mapped_neighbor_id == -1) continue;

                int found = 0;
                for (int k = 1; k <= adjacency_list_2[mapped_choice][0]; k++)
                    if ((mapped_neighbor_id | bond_type_mask) == adjacency_list_2[mapped_choice][k]) {
                        found = 1; break;
                    }

                if (!found) { is_mapping_valid = 0; break; }
            }

            if (!is_mapping_valid) {
                is_permutation_valid = 0;
                break;
            }

            current_assignment[id] = mapped_choice;
            group_sum += get_sq_distance(mol_data_1->coordinates, id, mol_data_2->coordinates, mapped_choice);
        }

        if (!is_permutation_valid) {
            continue;
        }

        // Recursively call for the next group
        rmsd_naive_rec_helper(mol_data_1, mol_data_2, candidate_data, candidate_squared_distances,
                              current_assignment, current_sum + group_sum, minimize, best_assignment, best_sum,
                              current_group + 1, group_count, groups);
    } while (next_permutation(permutation, group_size));

    free(permutation);
    for (int i = 0; i < group_size; i++) {
        int id = mol1_indices[i];
        current_assignment[id] = -1;
    }
}

double rmsd_naive(MolecularData* mol_data_1, MolecularData* mol_data_2, int* best_assignment, int minimize) {
    // Check for basic matches
    if (!check_for_basic_matches(mol_data_1, mol_data_2)) { 
        fprintf(stderr, "Basic matches failed. Molecules are not isomorphic.\n"); return INFINITY; 
    }


    int num_atoms = mol_data_1->num_atoms;
    double best_rmsd = INFINITY;

    int** candidate_data = get_candidate_data(num_atoms, mol_data_1->layer_data, mol_data_2->layer_data);
    if (!candidate_data) { fprintf(stderr, "Failed to get candidate data.\n"); return INFINITY; }

    double** candidate_squared_distances = get_candidate_squared_distances(num_atoms, mol_data_1->coordinates, mol_data_2->coordinates, candidate_data);
    if (!candidate_squared_distances) { 
        for (int i = 0; i < num_atoms; ++i) free(candidate_data[i]); free(candidate_data);
        fprintf(stderr, "Failed to get candidate squared distances.\n"); return INFINITY; 
    }

    double current_sum = 0.0;
    int* current_assignment = (int*)malloc(num_atoms * sizeof(int));
    if (!current_assignment) { 
        for (int i = 0; i < num_atoms; ++i) free(candidate_data[i]); free(candidate_data);
        fprintf(stderr, "Failed to allocate memory for current mapping.\n"); return INFINITY; 
    }


    for (int i = 0; i < num_atoms; i++) {
        if (candidate_data[i][0] == 1) {
            current_assignment[i] = candidate_data[i][1];
            current_sum += candidate_squared_distances[i][0];
        } else {
            current_assignment[i] = -1;
        }
    }

    int** groups = NULL;
    int groups_count = 0;

    for (int i = 0; i < num_atoms; i++) {
        if (current_assignment[i] != -1) continue;

        int before_count = 0, after_count = 0;
        for (int j = 0; j < num_atoms; j++) {
            if (mol_data_1->layer_data[i] == mol_data_1->layer_data[j]) {
                if (j <= i) before_count++;
                else after_count++;
            }
        }

        if (after_count == 0) {
            groups = (int**)realloc(groups, (groups_count + 1) * sizeof(int*));
            groups[groups_count] = (int*)malloc((before_count + 1) * sizeof(int));
            groups[groups_count][0] = before_count;

            for (int j = 0, k = 1; j < num_atoms; j++) {
                if (mol_data_1->layer_data[i] == mol_data_1->layer_data[j]) {
                    groups[groups_count][k++] = j;
                }
            }

            groups_count++;
        }
    }

    if (groups_count == 0) {
        if (minimize)
            current_sum = kabsch_squared_dist_sum(mol_data_1->num_atoms, mol_data_1->coordinates, mol_data_2->coordinates, current_assignment);
        best_rmsd = sqrt(current_sum / num_atoms);
        memcpy(best_assignment, current_assignment, num_atoms * sizeof(int));
    } else {
        double best_sum = INFINITY;
        rmsd_naive_rec_helper(mol_data_1, mol_data_2, candidate_data, candidate_squared_distances,
                              current_assignment, current_sum, minimize, best_assignment, &best_sum,
                              0, groups_count, groups);
        best_rmsd = sqrt(best_sum / num_atoms);
    }
    free(current_assignment);
    for (int i = 0; i < groups_count; i++) {
        free(groups[i]);
    }
    free(groups);
    for (int i = 0; i < num_atoms; ++i) free(candidate_data[i]); free(candidate_data);
    for (int i = 0; i < num_atoms; ++i) free(candidate_squared_distances[i]); free(candidate_squared_distances);
    
    return best_rmsd;
}
