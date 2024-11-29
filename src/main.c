#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "rmsd.h"

int main(int argc, char* argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <template_file_path> <query_file_path> [-h <if passed H atoms will be used as well>] [-b <if passed bond types will be used as well>]\n", argv[0]);
        return 1;
    }

    const char* filename_1 = argv[1];
    const char* filename_2 = argv[2];

    const char* hydrogen_flag = "-h";
    const char* bond_flag = "-b";

    int read_hydrogens = 0;
    int read_bonds = 0;

    for (int i = 3; i < argc; ++i) {
        if (strcmp(argv[i], hydrogen_flag) == 0)    read_hydrogens = 1;
        else if (strcmp(argv[i], bond_flag) == 0) read_bonds = 1;
    }

    int total_conf_count = 0;
    MolecularData* mol_data = read_input_files(filename_1, filename_2, read_hydrogens, read_bonds, &total_conf_count);

    MolecularData* template_mol_data = mol_data;

    int start, end;

    // Output the molecular data
    for (int i = 1; i < total_conf_count; ++i) {
        start = clock();
        
        int* best_assignment = (int*)malloc(template_mol_data->num_atoms * sizeof(int));
        if (!best_assignment) {
            fprintf(stderr, "Memory allocation failed for best assignment.\n");
            continue;
        }

        double best_rmsd = rmsd(template_mol_data, mol_data + i, best_assignment);

        printf("Best RMSD for conformer %d: %lf\n", i, best_rmsd);

        end = clock();

        printf("Time taken: %lf\n", (double)(end - start) / CLOCKS_PER_SEC);
    }

    return 0;
}
