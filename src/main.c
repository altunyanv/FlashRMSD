#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "rmsd.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <query_file_path> [<template_file_path>] \n\
            [-n <if passed naive version of algorithm will be run] \n\
            [-x <if passed the cross RMSD will be calculated for all conformations of query file] \n\
            [-h <if passed H atoms will be used as well>] \n\
            [-b <if passed bond types will be used as well>] \n\
            [-v <if passed verbosity is set to true (i.e. runtimes and additional comments)>] \n\
            [-a <if passed atom assignments will be printed]\n", argv[0]);
        return 1;
    }

    int argi = 1;    
    if (argv[argi][0] == '-') {
        fprintf(stderr, "Expected query file path, found option '%s'\n", argv[argi]);
        return 1;
    }

    const char* query_file_path = argv[argi++];
    const char* template_file_path = NULL;

    // Optional template file
    if (argi < argc && argv[argi][0] != '-') {
        template_file_path = argv[argi++];
    }

    int naive_rmsd = 0;
    int cross_rmsd = 0;
    int read_hydrogens = 0;
    int read_bonds = 0;
    int verbosity = 0;
    int print_assignments = 0;

    for (; argi < argc; ++argi) {
        if      (strcmp(argv[argi], "-n") == 0) naive_rmsd = 1;
        else if (strcmp(argv[argi], "-x") == 0) cross_rmsd = 1;
        else if (strcmp(argv[argi], "-h") == 0) read_hydrogens = 1;
        else if (strcmp(argv[argi], "-b") == 0) read_bonds = 1;
        else if (strcmp(argv[argi], "-v") == 0) verbosity = 1;
        else if (strcmp(argv[argi], "-a") == 0) print_assignments = 1;
        else {
            fprintf(stderr, "Unrecognized flag %s. Exiting...\n", argv[argi]);
            return 1;
        }
    }

    if (cross_rmsd) {
        MolecularData* mol_data = (MolecularData*)malloc(MAX_CONF_COUNT * sizeof(MolecularData));

        int conf_count = read_file(query_file_path, mol_data, read_hydrogens, read_bonds, MAX_CONF_COUNT);
        if (!conf_count) {
            fprintf(stderr, "Failed to read query file.\n");
            exit(1);
        }

        int* best_assignment = (int*)malloc(mol_data->num_atoms * sizeof(int));
        double* rmsd_table = (double*)malloc(conf_count * conf_count * sizeof(double));

        int start, end;
        start = clock();

        for (int i = 0; i < conf_count; ++i) {
            rmsd_table[i * conf_count + i] = 0.0;
            for (int j = i + 1; j < conf_count; ++j) {
                double rmsd_value = naive_rmsd ? rmsd_naive(mol_data + i, mol_data + j, best_assignment) : rmsd(mol_data + i, mol_data + j, best_assignment);
                rmsd_table[i * conf_count + j] = rmsd_value;
                rmsd_table[j * conf_count + i] = rmsd_value;
            }
        }

        end = clock();

        for (int i = 0; i < conf_count; ++i) {
            for (int j = 0; j < conf_count; ++j) {
                printf("%lf ", rmsd_table[i * conf_count + j]);
            }
            printf("\n");
        }

        if (verbosity) {
            printf("Cross RMSD calculation completed in %lf seconds.\n", (double)(end - start) / CLOCKS_PER_SEC);
        }
    } else {
        int total_conf_count = 0;
        MolecularData* mol_data = read_input_files(template_file_path, query_file_path, read_hydrogens, read_bonds, &total_conf_count);

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

            double best_rmsd = naive_rmsd ? rmsd_naive(template_mol_data, mol_data + i, best_assignment) : rmsd(template_mol_data, mol_data + i, best_assignment);

            end = clock();

            if (verbosity) {
                printf("Conformer %d: \n", i - 1);
                printf("RMSD: %lf\n", best_rmsd);
                if (print_assignments) {
                    printf("Assignment: \n");
                    for (int j = 0; j < template_mol_data->num_atoms; ++j) {
                        printf("%d -> %d\n", j, best_assignment[j]);
                    }
                }
                printf("Time taken: %lf\n", (double)(end - start) / CLOCKS_PER_SEC);
            } else {
                printf("%lf\n", best_rmsd);
                if (print_assignments) {
                    for (int j = 0; j < template_mol_data->num_atoms; ++j) {
                        printf("%d -> %d\n", j, best_assignment[j]);
                    }
                }
            }
        }
    }

    return 0;
}
