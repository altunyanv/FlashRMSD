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

int read_sdf_block_from_file_ptr(FILE* infile, MolecularData* mol_data, int read_hydrogens, int read_bonds)
{
    char line[256];

    /* 1) Skip the three header lines */
    if (!fgets(line, sizeof(line), infile)) return 0;
    if (!fgets(line, sizeof(line), infile)) return handle_error(
            "Unexpected end of file while reading header line 2.", infile);
    if (!fgets(line, sizeof(line), infile)) return handle_error(
            "Unexpected end of file while reading header line 3.", infile);

    /* 2) Peek at the next line to decide V2000 vs. V3000 (MDL CTAB) */
    if (!fgets(line, sizeof(line), infile))
        return handle_error("Unexpected end of file while reading counts.", infile);

    int use_v3000 = (strstr(line, "V3000") != NULL);

    int num_atoms_total = 0, num_bonds_total = 0;
    /* for V2000 we already have line; for V3000 we must skip to "M  V30 COUNTS" */
    if (use_v3000) {
        /* read V3000 COUNTS line */
        while (!strstr(line, "M  V30 COUNTS")) {
            if (!fgets(line, sizeof(line), infile))
                return handle_error("Failed to find V3000 COUNTS.", infile);
        }
        /* format: "M  V30 COUNTS natoms nbonds ..." */
        if (sscanf(line+14, "%d %d", &num_atoms_total, &num_bonds_total) < 2)
            return handle_error("Failed to parse V3000 COUNTS.", infile);
    } else {
        /* V2000, parse fixed columns */
        char buf_atoms[4]  = {0},
             buf_bonds[4]  = {0};
        if ((int)strlen(line) < 6)
            return handle_error("Counts line too short.", infile);
        memcpy(buf_atoms, line+0, 3);
        memcpy(buf_bonds, line+3, 3);
        if (sscanf(buf_atoms, "%d", &num_atoms_total) != 1 ||
            sscanf(buf_bonds, "%d", &num_bonds_total) != 1)
            return handle_error("Failed to parse V2000 counts.", infile);
    }

    /* allocate workspace */
    int*    atomic_numbers_all = malloc(num_atoms_total * sizeof(int));
    double* coordinates_all    = malloc(num_atoms_total * 3 * sizeof(double));
    int*    atom_index_map     = malloc(num_atoms_total * sizeof(int));
    if (!atomic_numbers_all || !coordinates_all || !atom_index_map)
        return handle_error("Memory allocation failed.", infile);

    int num_atoms = 0;

    /* 3) Read atom block */
    if (use_v3000) {
        /* skip to "M  V30 BEGIN ATOM" */
        while (!strstr(line, "M  V30 BEGIN ATOM")) {
            if (!fgets(line, sizeof(line), infile))
                return handle_error("Failed to find V3000 BEGIN ATOM.", infile);
        }
        /* read until "M  V30 END ATOM" */
        while (fgets(line, sizeof(line), infile)) {
            if (strstr(line, "M  V30 END ATOM")) break;
            /* parse: M  V30 <id> <x> <y> <z> <symbol> ... */
            int    idx;
            double x, y, z;
            char   sym[16];
            if (sscanf(line+6, "%d %s %lf %lf %lf", &idx, sym, &x, &y, &z) < 5)
                return handle_error("Failed to parse V3000 atom.", infile);
            int Z = get_atomic_number(sym);
            atom_index_map[idx-1] = (Z==1 && !read_hydrogens) ? -1 : num_atoms;
            if (!(Z==1 && !read_hydrogens)) {
                atomic_numbers_all[num_atoms] = Z;
                coordinates_all[3*num_atoms+0] = x;
                coordinates_all[3*num_atoms+1] = y;
                coordinates_all[3*num_atoms+2] = z;
                ++num_atoms;
            }
        }
    } else {
        /* V2000 atom block: next num_atoms_total lines */
        for (int i = 0; i < num_atoms_total; ++i) {
            if (!fgets(line, sizeof(line), infile))
                return handle_error("Unexpected EOF in atom block.", infile);

            char xs[11]={0}, ys[11]={0}, zs[11]={0}, symb[4]={0};
            strncpy(xs, line+0, 10);
            strncpy(ys, line+10,10);
            strncpy(zs, line+20,10);
            strncpy(symb,line+31,3);

            double x = atof(xs), y = atof(ys), z = atof(zs);
            /* trim symb */
            int s = 0, l = strlen(symb);
            while (l>0 && isspace(symb[l-1])) symb[--l]=0;
            while (isspace(symb[s])) ++s;

            int Z = get_atomic_number(symb+s);
            atom_index_map[i] = (Z==1 && !read_hydrogens) ? -1 : num_atoms;
            if (!(Z==1 && !read_hydrogens)) {
                atomic_numbers_all[num_atoms] = Z;
                coordinates_all[3*num_atoms+0]=x;
                coordinates_all[3*num_atoms+1]=y;
                coordinates_all[3*num_atoms+2]=z;
                ++num_atoms;
            }
        }
    }

    mol_data->num_atoms = num_atoms;
    mol_data->atomic_numbers =
        malloc(num_atoms * sizeof(int));
    mol_data->coordinates =
        malloc(num_atoms * 3 * sizeof(double));
    if (!mol_data->atomic_numbers || !mol_data->coordinates)
        return handle_error("Memory allocation failed in mol_data.", infile);

    memcpy(mol_data->atomic_numbers,
           atomic_numbers_all,
           num_atoms * sizeof(int));
    memcpy(mol_data->coordinates,
           coordinates_all,
           num_atoms * 3 * sizeof(double));
    free(atomic_numbers_all);
    free(coordinates_all);

    /* 4) Build adjacency list temp */
    int** adjacency_list_temp = calloc(num_atoms, sizeof(int*));
    int*  adjacency_sizes     = calloc(num_atoms, sizeof(int));
    if (!adjacency_list_temp || !adjacency_sizes)
        return handle_error("Memory allocation failed for adjacency.", infile);

    /* 5) Read bond block */
    if (use_v3000) {
        /* skip to "M  V30 BEGIN BOND" */
        while (!strstr(line, "M  V30 BEGIN BOND")) {
            if (!fgets(line, sizeof(line), infile))
                return handle_error("Failed to find V3000 BEGIN BOND.", infile);
        }
        /* read until "M  V30 END BOND" */
        while (fgets(line, sizeof(line), infile)) {
            if (strstr(line, "M  V30 END BOND")) break;
            /* parse: M  V30 <id> <a1> <a2> <type> ... */
            int    bidx, a1, a2;
            char   btype[8];
            if (sscanf(line+6, "%d %s %d %d",
                       &bidx,btype,&a1,&a2) < 4)
                return handle_error("Failed to parse V3000 bond.", infile);
            int bt = read_bonds ? get_bond_number(btype) : 0;
            int i1 = atom_index_map[a1-1],
                i2 = atom_index_map[a2-1];
            if (i1<0||i2<0) continue;
            /* append to temp lists */
            adjacency_list_temp[i1] = realloc(
                adjacency_list_temp[i1],
                ++adjacency_sizes[i1] * sizeof(int));
            adjacency_list_temp[i2] = realloc(
                adjacency_list_temp[i2],
                ++adjacency_sizes[i2] * sizeof(int));
            adjacency_list_temp[i1][adjacency_sizes[i1]-1] =
                i2 | (bt<<BOND_TYPE_SHIFT);
            adjacency_list_temp[i2][adjacency_sizes[i2]-1] =
                i1 | (bt<<BOND_TYPE_SHIFT);
        }
    } else {
        /* V2000 bond block: next num_bonds_total lines */
        for (int i = 0; i < num_bonds_total; ++i) {
            if (!fgets(line, sizeof(line), infile))
                return handle_error("Unexpected EOF in bond block.", infile);
            char a1s[4]={0}, a2s[4]={0}, bt[4]={0};
            strncpy(a1s,line+0,3);
            strncpy(a2s,line+3,3);
            strncpy(bt, line+6,3);
            int a1=atoi(a1s), a2=atoi(a2s);
            /* trim bt */
            int si=0, li=strlen(bt);
            while (li>0&&isspace(bt[li-1])) bt[--li]=0;
            while (isspace(bt[si])) ++si;
            int bti = read_bonds ? get_bond_number(bt+si) : 0;
            int i1 = atom_index_map[a1-1],
                i2 = atom_index_map[a2-1];
            if (i1<0||i2<0) continue;
            adjacency_list_temp[i1] = realloc(
                adjacency_list_temp[i1],
                ++adjacency_sizes[i1] * sizeof(int));
            adjacency_list_temp[i2] = realloc(
                adjacency_list_temp[i2],
                ++adjacency_sizes[i2] * sizeof(int));
            adjacency_list_temp[i1][adjacency_sizes[i1]-1] =
                i2 | (bti<<BOND_TYPE_SHIFT);
            adjacency_list_temp[i2][adjacency_sizes[i2]-1] =
                i1 | (bti<<BOND_TYPE_SHIFT);
        }
    }

    free(atom_index_map);

    /* 6) Fold temp adjacency into mol_data */
    mol_data->adjacency_list =
        malloc(num_atoms * sizeof(int*));
    if (!mol_data->adjacency_list)
        return handle_error("Memory allocation failed.", infile);
    for (int i = 0; i < num_atoms; ++i) {
        int n = adjacency_sizes[i];
        mol_data->adjacency_list[i] =
            malloc((n+1) * sizeof(int));
        if (!mol_data->adjacency_list[i])
            return handle_error("Memory allocation failed.", infile);
        mol_data->adjacency_list[i][0] = n;
        memcpy(mol_data->adjacency_list[i]+1,
               adjacency_list_temp[i],
               n * sizeof(int));
        free(adjacency_list_temp[i]);
    }
    free(adjacency_list_temp);
    free(adjacency_sizes);

    /* 7) Layer data */
    mol_data->layer_data = malloc(num_atoms * sizeof(int));
    if (!mol_data->layer_data)
        return handle_error("Memory allocation failed.", infile);
    if (!generate_layer_data(mol_data->num_atoms,
                             mol_data->atomic_numbers,
                             mol_data->adjacency_list,
                             mol_data->layer_data))
        return handle_error("Failed to generate layer data.", infile);

    /* 8) Consume record terminator */
    while (fgets(line, sizeof(line), infile)) {
        if (strstr(line, "$$$$"))  /* MDL Molfile end */
        {
            break;
        }
    }

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
    for (int i = 0; i < num_atoms; ++i) { adjacency_list_temp[i] = NULL; adjacency_sizes[i] = 0; }

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
