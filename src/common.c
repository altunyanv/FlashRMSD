#include "common.h"

int get_layer_id(int value) { return (value >> LAYER_ID_SHIFT) & BLOCK_MASK;}
int get_atom_id(int value) { return value & BLOCK_MASK; }

int get_bond_type(int value) { return (value >> BOND_TYPE_SHIFT) & BLOCK_MASK; }
int get_bonded_atom_id(int value) { return value & BLOCK_MASK; }

int compare_ints(const void* a, const void* b) {
    int ai = *(const int*)a, bi = *(const int*)b;   return (ai > bi) - (ai < bi);
}

int check_array_unordered_identity(int num_atoms, int* arr1, int* arr2) {
    int* arr1_sorted = (int*)malloc(num_atoms * sizeof(int));
    int* arr2_sorted = (int*)malloc(num_atoms * sizeof(int));

    memcpy(arr1_sorted, arr1, num_atoms * sizeof(int));
    memcpy(arr2_sorted, arr2, num_atoms * sizeof(int));

    qsort(arr1_sorted, num_atoms, sizeof(int), compare_ints);
    qsort(arr2_sorted, num_atoms, sizeof(int), compare_ints);

    for (int i = 0; i < num_atoms; i++)
        if (arr1_sorted[i] != arr2_sorted[i]) {
            free(arr1_sorted);
            free(arr2_sorted);

            return 0;
        }
    
    free(arr1_sorted);
    free(arr2_sorted);
    
    return 1;
}

