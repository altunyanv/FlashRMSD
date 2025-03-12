#include "io.h"

int get_atomic_number(const char* symbol) {
    int num_elements = sizeof(atom_table) / sizeof(AtomEntry);
    for (int i = 0; i < num_elements; ++i) {
        if (strcmp(atom_table[i].symbol, symbol) == 0) {
            return atom_table[i].atomic_number;
        }
    }
    return 0;
}

int get_bond_number(const char* bond_symbol) {
    if (strcmp(bond_symbol, "1") == 0) return 1;
    if (strcmp(bond_symbol, "2") == 0) return 2;
    if (strcmp(bond_symbol, "3") == 0) return 3;
    if (strcmp(bond_symbol, "ar") == 0) return 4;
    if (strcmp(bond_symbol, "am") == 0) return 5;
    if (strcmp(bond_symbol, "du") == 0) return 6;
    if (strcmp(bond_symbol, "un") == 0) return 7;
    return 0;
}

int handle_error(const char* message, FILE* infile) {
    if (infile) fclose(infile);
    fprintf(stderr, "%s\n", message);
    return 0;
}

int is_end_of_file(FILE *file) {
    int ch;
    
    while ((ch = fgetc(file)) != EOF) {
        if (!isspace(ch)) {
            ungetc(ch, file);
            return 0;
        }
    }
    
    return 1;
}

int read_sdf_block_from_file_ptr(FILE* infile, MolecularData* mol_data, int read_hydrogens, int read_bonds) {
    char line[256];

    // Skip the first three header lines
    if (!fgets(line, sizeof(line), infile)) return 0;
    if (!fgets(line, sizeof(line), infile)) return handle_error("Unexpected end of file while reading header line 2.", infile);
    if (!fgets(line, sizeof(line), infile)) return handle_error("Unexpected end of file while reading header line 3.", infile);

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

        char x_str[11] = {0}, y_str[11] = {0}, z_str[11] = {0}, atom_symbol[4] = {0};
        double* coordinates_shifted = coordinates_all + 3 * num_atoms;
        strncpy(x_str, line + 0, 10);   
        strncpy(y_str, line + 10, 10); 
        strncpy(z_str, line + 20, 10); 
        strncpy(atom_symbol, line + 31, 3);
        
        // Convert coordinates to double
        coordinates_shifted[0] = atof(x_str);
        coordinates_shifted[1] = atof(y_str);
        coordinates_shifted[2] = atof(z_str);

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

        char atom1_str[4] = {0}, atom2_str[4] = {0}, bond_type[4] = {0};  
        int atom1, atom2, bond_type_int;      
        
        strncpy(atom1_str, line, 3);        atom1 = atoi(atom1_str);
        strncpy(atom2_str, line + 3, 3);    atom2 = atoi(atom2_str);

        // Trim leading and trailing spaces, because of length
        strncpy(bond_type, line + 6, 3);
        int start = 0, len = strlen(bond_type);
        while (len > 0 && isspace(bond_type[len - 1])) bond_type[--len] = '\0';
        while (isspace(bond_type[start])) ++start;

        bond_type_int = read_bonds * get_bond_number(bond_type + start);
        
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
        adjacency_list_temp[new_idx1][adjacency_sizes[new_idx1] - 1] = (new_idx2 | (bond_type_int << BOND_TYPE_SHIFT));
        adjacency_list_temp[new_idx2][adjacency_sizes[new_idx2] - 1] = (new_idx1 | (bond_type_int << BOND_TYPE_SHIFT));
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
    
    while (fgets(line, sizeof(line), infile))
        if (strstr(line, "$$$$")) break;

    return 1;
}

int read_mol2_block_from_file_ptr(FILE* infile, MolecularData* mol_data, int read_hydrogens, int read_bonds) {
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
        char atom_name[6] = {0};
        double* coordinates_shifted = coordinates_all + 3 * num_atoms;
        if (sscanf(line, 
                   "%*d%*s%lf%lf%lf%5s",
                   coordinates_shifted, 
                   coordinates_shifted + 1,
                   coordinates_shifted + 2,
                   atom_name) < 4) return handle_error("Failed to parse atom block.", infile);
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

        int atom1, atom2, bond_type_int;
        char bond_type[3] = {0};
        if (sscanf(line, 
                   "%*d%d%d%2s", 
                   &atom1, 
                   &atom2, 
                   bond_type) < 3) return handle_error("Failed to parse bond block.", infile);

        bond_type_int = read_bonds * get_bond_number(bond_type);

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
        adjacency_list_temp[new_idx1][adjacency_sizes[new_idx1] - 1] = (new_idx2 | (bond_type_int << BOND_TYPE_SHIFT));
        adjacency_list_temp[new_idx2][adjacency_sizes[new_idx2] - 1] = (new_idx1 | (bond_type_int << BOND_TYPE_SHIFT));
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

    while (!feof(infile) && fgets(line, sizeof(line), infile))
        if (strstr(line, "@<TRIPOS>MOLECULE")) break;

    if (!feof(infile)) fseek(infile, -strlen(line), SEEK_CUR);

    return 1;
}

int read_sdf_file(const char* filename, MolecularData* mol_data, int read_hydrogens, int read_bonds, int max_conf_count) {
    FILE* infile = fopen(filename, "r");
    if (!infile) return handle_error("Cannot open file", infile);

    int i = 0;
    for (; i < max_conf_count; ++i) {
        if (!read_sdf_block_from_file_ptr(infile, mol_data + i, read_hydrogens, read_bonds)) {
            if (i == 0) { fprintf(stderr, "Failed to read SDF file %s.\n", filename); return 0; }
            return i;
        }
    }
    
    if (!is_end_of_file(infile))
        fprintf(stdout, "[Warning] Reached maximum number of conformers in %s. Returning first %d conformers\n", filename, max_conf_count);

    fclose(infile);
    return max_conf_count;
}

int read_mol2_file(const char* filename, MolecularData* mol_data, int read_hydrogens, int read_bonds, int max_conf_count) {
    FILE* infile = fopen(filename, "r");
    if (!infile) return handle_error("Cannot open file", infile);

    int i = 0;
    for (; i < max_conf_count; ++i) {
        if (!read_mol2_block_from_file_ptr(infile, mol_data + i, read_hydrogens, read_bonds)) {
            if (i == 0) { fprintf(stderr, "Failed to read MOL2 file %s.\n", filename); return 0; }
            return i;
        }
    }
    
    if (!is_end_of_file(infile))
        fprintf(stdout, "[Warning] Reached maximum number of conformers in %s. Returning first %d conformers\n", filename, max_conf_count);

    fclose(infile);
    return max_conf_count;
}

int read_file(const char* filename, MolecularData* mol_data, int read_hydrogens, int read_bonds, int max_conf_count) {
    FILE* infile = fopen(filename, "r");
    if (!infile) return handle_error("Cannot open file", infile);

    // Check file extension
    const char* extension = strrchr(filename, '.');
    if (!extension) return handle_error("Failed to get file extension.", infile);

    if (strcmp(extension, ".sdf") == 0) {
        int num_confs = read_sdf_file(filename, mol_data, read_hydrogens, read_bonds, max_conf_count);
        fclose(infile);
        return num_confs;
    } else if (strcmp(extension, ".mol2") == 0) {
        int num_confs = read_mol2_file(filename, mol_data, read_hydrogens, read_bonds, max_conf_count);
        fclose(infile);
        return num_confs;
    } else {
        fprintf(stderr, "Unsupported file format: %s\n", extension);
        fclose(infile);
        return 0;
    }
}

MolecularData* read_input_files(const char* filename_1, const char* filename_2, int read_hydrogens, int read_bonds, int* total_conf_count) {
    MolecularData* mol_data = (MolecularData*)malloc((1 + MAX_CONF_COUNT) * sizeof(MolecularData));

    if (!mol_data) { fprintf(stderr, "Memory allocation failed for mol_data.\n"); return NULL; }

    if (!read_file(filename_1, mol_data, read_hydrogens, read_bonds, 1)) {
        fprintf(stderr, "Failed to read template file. Single molecule conformation should be present in template file.\n");
        exit(1);
    }

    int conf_count = read_file(filename_2, mol_data + 1, read_hydrogens, read_bonds, MAX_CONF_COUNT);
    if (!conf_count) {
        fprintf(stderr, "Failed to read query file.\n");
        exit(1);
    }

    *total_conf_count = conf_count + 1;

    return mol_data;
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

void free_all_mol_data_allocations(MolecularData* mol_data, int total_conf_count) {
    for (int i = 0; i < total_conf_count; ++i) {
        free_mol_data(mol_data + i);
    }
    free(mol_data);
}
