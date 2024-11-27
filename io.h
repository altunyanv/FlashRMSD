#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "utils.h"

typedef struct {
    int num_atoms;

    int* atomic_numbers;
    double* coordinates;

    int** adjacency_list;

    int* layer_data;
} MolecularData;

int handle_error(const char* message, FILE* infile) {
    if (infile) fclose(infile);
    fprintf(stderr, "%s\n", message);
    return 0;
}

const int MAX_CONF_COUNT = 100;

int read_sdf_block_from_file_ptr(FILE* infile, MolecularData* mol_data, int read_hydrogens) {
    char line[256];

    // Skip the first three header lines
    if (!fgets(line, sizeof(line), infile)) return 0;
    if (!fgets(line, sizeof(line), infile)) return handle_error("Unexpected end of file while reading header line 2.", infile); // Line 2
    if (!fgets(line, sizeof(line), infile)) return handle_error("Unexpected end of file while reading header line 3.", infile); // Line 3

    // Read counts line
    if (!fgets(line, sizeof(line), infile)) return handle_error("Unexpected end of file while reading counts line.", infile);

    int num_atoms_total = 0, num_bonds_total = 0, chiral_flag = 0;
    char version[7] = "";

    // Parse counts line
    if (sscanf(line, "%3d%3d%*3d%*3d%*3d%3d%*10c%6s", &num_atoms_total, &num_bonds_total, &chiral_flag, version) < 2) return handle_error("Failed to parse counts line.", infile);

    int* atomic_numbers_all = (int*)malloc(num_atoms_total * sizeof(int));
    double* coordinates_all = (double*)malloc(num_atoms_total * 3 * sizeof(double));
    int* atom_index_map = (int*)malloc(num_atoms_total * sizeof(int));

    if (!atomic_numbers_all || !coordinates_all || !atom_index_map) return handle_error("Memory allocation failed for atomic numbers/coordinates/atom index map.", infile);

    int num_atoms = 0; // Number of atoms after possibly excluding hydrogens

    // Read atom block
    for (int i = 0; i < num_atoms_total; ++i) {
        if (!fgets(line, sizeof(line), infile)) return handle_error("Unexpected end of file while reading atom block.", infile);

        // Extract x, y, z coordinates
        char x_str[11] = {0}, y_str[11] = {0}, z_str[11] = {0};
        strncpy(x_str, line + 0, 10);
        strncpy(y_str, line + 10, 10);
        strncpy(z_str, line + 20, 10);

        double x = atof(x_str);
        double y = atof(y_str);
        double z = atof(z_str);

        // Extract atom symbol
        char atom_symbol[4] = {0};
        strncpy(atom_symbol, line + 31, 3);

        // Trim leading and trailing spaces
        int start = 0, len = strlen(atom_symbol);
        while (len > 0 && isspace(atom_symbol[len - 1])) atom_symbol[--len] = '\0';
        while (isspace(atom_symbol[start])) ++start;

        // Map atom symbol to atomic number
        int atomic_number = get_atomic_number(atom_symbol + start);

        // Decide whether to include this atom
        if (!read_hydrogens && atomic_number == 1) { atom_index_map[i] = -1; continue; }

        // Store atom data
        atomic_numbers_all[num_atoms] = atomic_number;
        coordinates_all[3 * num_atoms] = x;
        coordinates_all[3 * num_atoms + 1] = y;
        coordinates_all[3 * num_atoms + 2] = z;

        atom_index_map[i] = num_atoms; // Map old index to new index
        num_atoms++;
    }

    mol_data->num_atoms = num_atoms;

    // Allocate memory for atomic_numbers and coordinates in mol_data
    mol_data->atomic_numbers = (int*)malloc(num_atoms * sizeof(int));
    mol_data->coordinates = (double*)malloc(num_atoms * 3 * sizeof(double));
    if (!mol_data->atomic_numbers || !mol_data->coordinates) return handle_error("Memory allocation failed for atomic numbers/coordinates in mol_data.", infile);

    // Copy the included atoms to mol_data
    memcpy(mol_data->atomic_numbers, atomic_numbers_all, num_atoms * sizeof(int));
    memcpy(mol_data->coordinates, coordinates_all, num_atoms * 3 * sizeof(double));

    free(atomic_numbers_all);   free(coordinates_all);

    // Temporary storage for adjacency list
    int** adjacency_list_temp = (int**)malloc(num_atoms * sizeof(int*));
    int* adjacency_sizes = (int*)malloc(num_atoms * sizeof(int));

    if (!adjacency_list_temp || !adjacency_sizes) return handle_error("Memory allocation failed for adjacency list.", infile);

    // Initialize adjacency list pointers to NULL
    for (int i = 0; i < num_atoms; ++i) adjacency_list_temp[i] = NULL;

    // Read bond block
    for (int i = 0; i < num_bonds_total; ++i) {
        if (!fgets(line, sizeof(line), infile)) return handle_error("Unexpected end of file while reading bond block.", infile);

        char atom1_str[4] = {0}, 
             atom2_str[4] = {0};
        strncpy(atom1_str, line + 0, 3);
        strncpy(atom2_str, line + 3, 3);

        int atom1 = atoi(atom1_str) - 1, 
            atom2 = atoi(atom2_str) - 1;

        int new_idx1 = atom_index_map[atom1],
            new_idx2 = atom_index_map[atom2];

        // If either atom is excluded, skip this bond
        if (new_idx1 == -1 || new_idx2 == -1)   continue;

        // Increment adjacency sizes
        adjacency_sizes[new_idx1]++;    adjacency_sizes[new_idx2]++;

        // Reallocate adjacency lists
        adjacency_list_temp[new_idx1] = (int*)realloc(adjacency_list_temp[new_idx1], adjacency_sizes[new_idx1] * sizeof(int));
        adjacency_list_temp[new_idx2] = (int*)realloc(adjacency_list_temp[new_idx2], adjacency_sizes[new_idx2] * sizeof(int));

        if (!adjacency_list_temp[new_idx1] || !adjacency_list_temp[new_idx2]) return handle_error("Memory allocation failed for adjacency list.", infile);

        // Add each atom to the other's adjacency list
        adjacency_list_temp[new_idx1][adjacency_sizes[new_idx1] - 1] = new_idx2;
        adjacency_list_temp[new_idx2][adjacency_sizes[new_idx2] - 1] = new_idx1;
    }

    free(atom_index_map);

    // Build adjacency_list in MolecularData
    mol_data->adjacency_list = (int**)malloc(num_atoms * sizeof(int*));
    if (!mol_data->adjacency_list) return handle_error("Memory allocation failed for adjacency list in mol_data.", infile);

    for (int i = 0; i < num_atoms; ++i) {
        int num_neighbors = adjacency_sizes[i];
        mol_data->adjacency_list[i] = (int*)malloc((num_neighbors + 1) * sizeof(int)); // First element is num_neighbors
        if (!mol_data->adjacency_list[i]) return handle_error("Memory allocation failed for adjacency list in mol_data.", infile);

        mol_data->adjacency_list[i][0] = num_neighbors;
        memcpy(mol_data->adjacency_list[i] + 1, adjacency_list_temp[i], num_neighbors * sizeof(int));

        free(adjacency_list_temp[i]);
    }

    // Free temporary adjacency sizes and list
    free(adjacency_list_temp);  free(adjacency_sizes);

    mol_data->layer_data = (int*)malloc(num_atoms * sizeof(int));
    if (!mol_data->layer_data) return handle_error("Memory allocation failed for layer data.", infile);

    if (!generate_layer_data(mol_data->num_atoms, mol_data->atomic_numbers, mol_data->adjacency_list, mol_data->layer_data))
        return handle_error("Failed to generate layer data.", infile);
    
    while (fgets(line, sizeof(line), infile)) {
        if (strstr(line, "$$$$")) break;
    }

    return 1;
}

int read_mol2_block_from_file_ptr(FILE* infile, MolecularData* mol_data, int read_hydrogens) {
    char line[256];

    int found_molecule = 0;
    while (fgets(line, sizeof(line), infile))
        if (strstr(line, "@<TRIPOS>MOLECULE")) {
            found_molecule = 1;
            break;
        }

    if (found_molecule == 0)    return 0;
    
    // Read molecule name
    if (!fgets(line, sizeof(line), infile)) return handle_error("Unexpected end of file while reading molecule name.", infile);

    // Read counts line
    if (!fgets(line, sizeof(line), infile)) return handle_error("Unexpected end of file while reading counts line.", infile);

    int num_atoms_total = 0, num_bonds_total = 0, num_substructures = 0, num_features = 0, num_sets = 0;
    if (sscanf(line, "%d%d%d%d%d", &num_atoms_total, &num_bonds_total, &num_substructures, &num_features, &num_sets) < 2) return handle_error("Failed to parse counts line.", infile);

    int* atomic_numbers_all = (int*)malloc(num_atoms_total * sizeof(int));
    double* coordinates_all = (double*)malloc(num_atoms_total * 3 * sizeof(double));
    int* atom_index_map = (int*)malloc(num_atoms_total * sizeof(int));

    if (!atomic_numbers_all || !coordinates_all || !atom_index_map) return handle_error("Memory allocation failed for atomic numbers/coordinates/atom index map.", infile);

    int num_atoms = 0; // Number of atoms after possibly excluding hydrogens

    // Read atom block
    int found_atom_block = 0;
    while (fgets(line, sizeof(line), infile))
        if (strstr(line, "@<TRIPOS>ATOM")) {
            found_atom_block = 1;
            break;
        }

    if (found_atom_block == 0) return handle_error("Failed to find atom block.", infile);

    for (int i = 0; i < num_atoms_total; ++i) {
        if (!fgets(line, sizeof(line), infile)) return handle_error("Unexpected end of file while reading atom block.", infile);

        // Extract atom data
        char atom_name[6] = {0}, x_str[11] = {0}, y_str[11] = {0}, z_str[11] = {0};
        if (sscanf(line, "%*d%*s%lf%lf%lf%5s", coordinates_all + 3 * num_atoms, coordinates_all + 3 * num_atoms + 1, coordinates_all + 3 * num_atoms + 2, atom_name) < 4) return handle_error("Failed to parse atom block.", infile);
        char* atom_symbol = strtok(atom_name, ".");

        // Trim leading and trailing spaces
        int start = 0, len = strlen(atom_name);
        while (len > 0 && isspace(atom_name[len - 1])) atom_name[--len] = '\0';
        while (isspace(atom_name[start])) ++start;
        
        // Map atom symbol to atomic number
        int atomic_number = get_atomic_number(atom_name + start);

        // Decide whether to include this atom
        if (!read_hydrogens && atomic_number == 1) { atom_index_map[i] = -1; continue; }

        // Store atom data
        atomic_numbers_all[num_atoms] = atomic_number;
        atom_index_map[i] = num_atoms; 
        num_atoms++;
    }

    mol_data->num_atoms = num_atoms;

    // Allocate memory for atomic_numbers and coordinates in mol_data
    mol_data->atomic_numbers = (int*)malloc(num_atoms * sizeof(int));
    mol_data->coordinates = (double*)malloc(num_atoms * 3 * sizeof(double));

    if (!mol_data->atomic_numbers || !mol_data->coordinates) return handle_error("Memory allocation failed for atomic numbers/coordinates in mol_data.", infile);

    // Copy the included atoms to mol_data
    memcpy(mol_data->atomic_numbers, atomic_numbers_all, num_atoms * sizeof(int));
    memcpy(mol_data->coordinates, coordinates_all, num_atoms * 3 * sizeof(double));

    free(atomic_numbers_all);   free(coordinates_all);

    // Temporary storage for adjacency list
    int** adjacency_list_temp = (int**)malloc(num_atoms * sizeof(int*));
    int* adjacency_sizes = (int*)malloc(num_atoms * sizeof(int));

    if (!adjacency_list_temp || !adjacency_sizes) return handle_error("Memory allocation failed for adjacency list.", infile);

    // Initialize adjacency list pointers to NULL
    for (int i = 0; i < num_atoms; ++i) adjacency_list_temp[i] = NULL;

    // Read bond block
    int found_bond_block = 0;
    while (fgets(line, sizeof(line), infile))
        if (strstr(line, "@<TRIPOS>BOND")) {
            found_bond_block = 1;
            break;
        }

    if (found_bond_block == 0) return handle_error("Failed to find bond block.", infile);

    for (int i = 0; i < num_bonds_total; ++i) {
        if (!fgets(line, sizeof(line), infile)) return handle_error("Unexpected end of file while reading bond block.", infile);

        int atom1 = 0, atom2 = 0;
        char bond_type[3] = {0};
        if (sscanf(line, "%*d%d%d%2s", &atom1, &atom2, bond_type) < 3) return handle_error("Failed to parse bond block.", infile);

        int new_idx1 = atom_index_map[atom1 - 1],
            new_idx2 = atom_index_map[atom2 - 1];

        // If either atom is excluded, skip this bond
        if (new_idx1 == -1 || new_idx2 == -1)   continue;

        // Increment adjacency sizes
        adjacency_sizes[new_idx1]++;    adjacency_sizes[new_idx2]++;

        // Reallocate adjacency lists
        adjacency_list_temp[new_idx1] = (int*)realloc(adjacency_list_temp[new_idx1], adjacency_sizes[new_idx1] * sizeof(int));
        adjacency_list_temp[new_idx2] = (int*)realloc(adjacency_list_temp[new_idx2], adjacency_sizes[new_idx2] * sizeof(int));

        if (!adjacency_list_temp[new_idx1] || !adjacency_list_temp[new_idx2]) return handle_error("Memory allocation failed for adjacency list.", infile);

        // Add each atom to the other's adjacency list
        adjacency_list_temp[new_idx1][adjacency_sizes[new_idx1] - 1] = new_idx2;
        adjacency_list_temp[new_idx2][adjacency_sizes[new_idx2] - 1] = new_idx1;
    }

    free(atom_index_map);

    // Build adjacency_list in MolecularData
    mol_data->adjacency_list = (int**)malloc(num_atoms * sizeof(int*));

    if (!mol_data->adjacency_list) return handle_error("Memory allocation failed for adjacency list in mol_data.", infile);

    for (int i = 0; i < num_atoms; ++i) {
        int num_neighbors = adjacency_sizes[i];
        mol_data->adjacency_list[i] = (int*)malloc((num_neighbors + 1) * sizeof(int)); // First element is num_neighbors
        if (!mol_data->adjacency_list[i]) return handle_error("Memory allocation failed for adjacency list in mol_data.", infile);

        mol_data->adjacency_list[i][0] = num_neighbors;
        memcpy(mol_data->adjacency_list[i] + 1, adjacency_list_temp[i], num_neighbors * sizeof(int));

        free(adjacency_list_temp[i]);
    }

    // Free temporary adjacency sizes and list
    free(adjacency_list_temp);  free(adjacency_sizes);

    mol_data->layer_data = (int*)malloc(num_atoms * sizeof(int));

    if (!mol_data->layer_data) return handle_error("Memory allocation failed for layer data.", infile);

    if (!generate_layer_data(mol_data->num_atoms, mol_data->atomic_numbers, mol_data->adjacency_list, mol_data->layer_data))
        return handle_error("Failed to generate layer data.", infile);

    return 1;
}

int read_sdf_file(const char* filename, MolecularData* mol_data, int read_hydrogens, int max_conf_count) {
    FILE* infile = fopen(filename, "r");
    if (!infile) return handle_error("Cannot open file", infile);

    int i = 0;
    for (; i < max_conf_count; ++i) {
        if (!read_sdf_block_from_file_ptr(infile, mol_data + i, read_hydrogens)) {
            if (i == 0) { fprintf(stderr, "Failed to read SDF file %s.\n", filename); return 0; }
            return i;
        }
    }
    
    if (i == max_conf_count) fclose(infile);

    fprintf(stdout, "[Warning] Reached maximum number of conformers in %s. Returning first %d conformers\n", filename, max_conf_count);
    return max_conf_count;
}

int read_mol2_file(const char* filename, MolecularData* mol_data, int read_hydrogens, int max_conf_count) {
    FILE* infile = fopen(filename, "r");
    if (!infile) return handle_error("Cannot open file", infile);

    int i = 0;
    for (; i < max_conf_count; ++i) {
        if (!read_mol2_block_from_file_ptr(infile, mol_data + i, read_hydrogens)) {
            if (i == 0) { fprintf(stderr, "Failed to read MOL2 file %s.\n", filename); return 0; }
            return i;
        }
    }
    
    if (i == max_conf_count) fclose(infile);

    fprintf(stdout, "[Warning] Reached maximum number of conformers in %s. Returning first %d conformers\n", filename, max_conf_count);
    return max_conf_count;
}

int read_file(const char* filename, MolecularData* mol_data, int read_hydrogens, int max_conf_count) {
    FILE* infile = fopen(filename, "r");
    if (!infile) return handle_error("Cannot open file", infile);

    // Check file extension
    const char* extension = strrchr(filename, '.');
    if (!extension) return handle_error("Failed to get file extension.", infile);

    if (strcmp(extension, ".sdf") == 0) {
        int num_confs = read_sdf_file(filename, mol_data, read_hydrogens, max_conf_count);
        fclose(infile);
        return num_confs;
    } else if (strcmp(extension, ".mol2") == 0) {
        int num_confs = read_mol2_file(filename, mol_data, read_hydrogens, max_conf_count);
        fclose(infile);
        return num_confs;
    } else {
        fprintf(stderr, "Unsupported file format: %s\n", extension);
        fclose(infile);
        return 0;
    }
}

void free_mol_data(MolecularData* mol_data) {
    free(mol_data->atomic_numbers);
    free(mol_data->coordinates);
    for (int i = 0; i < mol_data->num_atoms; ++i) {
        free(mol_data->adjacency_list[i]);
    }
    free(mol_data->adjacency_list);
    free(mol_data->layer_data);
}

void free_all_mol_data_allocations(MolecularData* template_mol_data, MolecularData* query_mol_data, int query_conf_count) {
    free_mol_data(template_mol_data);
    for (int i = 0; i < query_conf_count; ++i) {
        free_mol_data(query_mol_data + i);
    }
    free(query_mol_data);
}

#endif // IO_H