#ifndef COMMON_H
#define COMMON_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define MAX_BLOCK_LENGTH 10     // Used both in adjacency_list items packing and layer value packing
#define LAYER_ID_SHIFT 20       // 2 * MAX_BLOCK_LENGTH
#define ATOMIC_NUMBER_SHIFT 10  // MAX_BLOCK_LENGTH
#define BLOCK_MASK 1023         // 2^MAX_BLOCK_LENGTH - 1

#define BOND_TYPE_SHIFT 10      // MAX_BLOCK_LENGTH

// ! Note: Maximum supported atoms count is 1023

int get_layer_id(int);
int get_atom_id(int);

int get_bond_type(int);
int get_bonded_atom_id(int);

int compare_ints(const void*, const void*);
int check_array_unordered_identity(int, int*, int*);

#endif // COMMON_H