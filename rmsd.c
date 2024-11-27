#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

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
    if (!check_for_basic_matches(mol_data_1, mol_data_2)) { fprintf(stderr, "Basic matches failed.\n"); return INFINITY; }
    int num_atoms = mol_data_1->num_atoms;

    // Get candidate data
    int** candidate_data = get_candidate_data(num_atoms, mol_data_1->layer_data, mol_data_2->layer_data);
    if (!candidate_data) { fprintf(stderr, "Failed to get candidate data.\n"); return INFINITY; }

    // double total_assignments = 1;
    // for (int i = 0; i < num_atoms; ++i)
    //     if (candidate_data[i][0] > 1) {
    //         printf("Atom %d: ", i);
    //         for (int j = 1; j <= candidate_data[i][0]; ++j)
    //             printf("%d ", candidate_data[i][j]);
    //         puts("");
    //     }

    double** candidate_squared_distances = get_candidate_squared_distances(num_atoms, mol_data_1->coordinates, mol_data_2->coordinates, candidate_data);
    if (!candidate_squared_distances) { 
        for (int i = 0; i < num_atoms; ++i) free(candidate_data[i]); free(candidate_data);
        fprintf(stderr, "Failed to get candidate squared distances.\n"); return INFINITY; 
    }

    // Search for the best assignment
    
    if (!best_assignment) { 
        for (int i = 0; i < num_atoms; ++i) free(candidate_data[i]); free(candidate_data);
        for (int i = 0; i < num_atoms; ++i) free(candidate_squared_distances[i]); free(candidate_squared_distances);
        fprintf(stderr, "Memory allocation failed for assignment.\n"); return INFINITY; 
    }

    const SearchData* data = init_search_data(num_atoms, candidate_data, candidate_squared_distances, mol_data_1->adjacency_list, mol_data_2->adjacency_list);
    
    if (!data) { 
        for (int i = 0; i < num_atoms; ++i) free(candidate_data[i]); free(candidate_data);
        for (int i = 0; i < num_atoms; ++i) free(candidate_squared_distances[i]); free(candidate_squared_distances);
        free(best_assignment);
        fprintf(stderr, "Failed to initialize recursive state.\n"); return INFINITY; 
    }
    
    double best_rmsd = search_assignment_recurse(data, best_assignment);
    // double best_rmsd = search_assignment(data, best_assignment);
    
    for (int i = 0; i < num_atoms; ++i) free(candidate_data[i]); free(candidate_data);
    for (int i = 0; i < num_atoms; ++i) free(candidate_squared_distances[i]); free(candidate_squared_distances);

    return best_rmsd;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <template_file_path> <query_file_path>\n", argv[0]);
        return 1;
    }

    const char* filename_1 = argv[1];
    const char* filename_2 = argv[2];

    // const char* filename_1 = "data/1azx_ligand.sdf";
    // const char* filename_2 = "data/1azx_ligand_2.sdf";

    int start = clock();
    MolecularData template_mol_data;

    if (read_file(filename_1, &template_mol_data, 0, MAX_CONF_COUNT) != 1) {
        fprintf(stderr, "Failed to read template file. Single molecule conformation should be present in template file.\n");
        return 1;
    }

    MolecularData* query_mol_data = (MolecularData*)malloc(MAX_CONF_COUNT * sizeof(MolecularData));

    int conf_count = read_file(filename_2, query_mol_data, 0, MAX_CONF_COUNT);
    if (!conf_count) {
        fprintf(stderr, "Failed to read query file.\n");
        return 1;
    }
    int end = clock();
    printf("Time taken to read files: %lf\n", (double)(end - start) / CLOCKS_PER_SEC);

    // Output the molecular data
    query_mol_data = (MolecularData*)realloc(query_mol_data, conf_count * sizeof(MolecularData));

    for (int i = 0; i < conf_count; ++i) {
        start = clock();
        
        int* best_assignment = (int*)malloc(template_mol_data.num_atoms * sizeof(int));
        if (!best_assignment) {
            fprintf(stderr, "Memory allocation failed for best assignment.\n");
            continue;
        }
        double best_rmsd = rmsd(&template_mol_data, query_mol_data + i, best_assignment);

        printf("Best RMSD for conformer %d: %lf\n", i, best_rmsd);

        end = clock();

        printf("Time taken: %lf\n", (double)(end - start) / CLOCKS_PER_SEC);
    }
    
    free_all_mol_data_allocations(&template_mol_data, query_mol_data, conf_count);

    return 0;
}
