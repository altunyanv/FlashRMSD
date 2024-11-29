#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "common.h"

#include "utils.h"

typedef struct {
    const char* symbol;
    int atomic_number;
} AtomEntry;

static AtomEntry atom_table[] = {
    {"H", 1},   {"He", 2},  {"Li", 3},  {"Be", 4},  {"B", 5},
    {"C", 6},   {"N", 7},   {"O", 8},   {"F", 9},   {"Ne", 10},
    {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15},
    {"S", 16},  {"Cl", 17}, {"Ar", 18}, {"K", 19},  {"Ca", 20},
    {"Sc", 21}, {"Ti", 22}, {"V", 23},  {"Cr", 24}, {"Mn", 25},
    {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30},
    {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Se", 34}, {"Br", 35},
    {"Kr", 36}, {"Rb", 37}, {"Sr", 38}, {"Y", 39},  {"Zr", 40},
    {"Nb", 41}, {"Mo", 42}, {"Tc", 43}, {"Ru", 44}, {"Rh", 45},
    {"Pd", 46}, {"Ag", 47}, {"Cd", 48}, {"In", 49}, {"Sn", 50},
    {"Sb", 51}, {"Te", 52}, {"I", 53},  {"Xe", 54}, {"Cs", 55},
    {"Ba", 56}, {"La", 57}, {"Ce", 58}, {"Pr", 59}, {"Nd", 60},
    {"Pm", 61}, {"Sm", 62}, {"Eu", 63}, {"Gd", 64}, {"Tb", 65},
    {"Dy", 66}, {"Ho", 67}, {"Er", 68}, {"Tm", 69}, {"Yb", 70},
    {"Lu", 71}, {"Hf", 72}, {"Ta", 73}, {"W", 74},  {"Re", 75},
    {"Os", 76}, {"Ir", 77}, {"Pt", 78}, {"Au", 79}, {"Hg", 80},
    {"Tl", 81}, {"Pb", 82}, {"Bi", 83}, {"Po", 84}, {"At", 85},
    {"Rn", 86}, {"Fr", 87}, {"Ra", 88}, {"Ac", 89}, {"Th", 90},
    {"Pa", 91}, {"U", 92},  {"Np", 93}, {"Pu", 94}, {"Am", 95},
    {"Cm", 96}, {"Bk", 97}, {"Cf", 98}, {"Es", 99}, {"Fm", 100},
    {"Md", 101},{"No", 102},{"Lr", 103},{"Rf", 104},{"Db", 105},
    {"Sg", 106},{"Bh", 107},{"Hs", 108},{"Mt", 109},{"Ds", 110},
    {"Rg", 111},{"Cn", 112},{"Nh", 113},{"Fl", 114},{"Mc", 115},
    {"Lv", 116},{"Ts", 117},{"Og", 118}
};

int get_atomic_number(const char*);

typedef struct {
    int num_atoms;

    int* atomic_numbers;
    double* coordinates;

    int** adjacency_list;

    int* layer_data;
} MolecularData;

int handle_error(const char*, FILE*);

#define MAX_CONF_COUNT 100

int read_sdf_block_from_file_ptr(FILE*, MolecularData*, int, int);
int read_mol2_block_from_file_ptr(FILE*, MolecularData*, int, int);

int read_sdf_file(const char*, MolecularData*, int, int, int);
int read_mol2_file(const char*, MolecularData*, int, int, int);

int read_file(const char*, MolecularData*, int, int, int);

MolecularData* read_input_files(const char*, const char*, int, int, int*);
void free_mol_data(MolecularData*);
void free_all_mol_data_allocations(MolecularData*, int);

#endif // IO_H