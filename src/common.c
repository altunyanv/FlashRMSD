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

void swap(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

void reverse_subarray(int arr[], int start, int end) {
    while (start < end) {
        swap(&arr[start], &arr[end]);
        start++;
        end--;
    }
}

int next_permutation(int arr[], int n) {
    if (n <= 1) return 0;

    int i = n - 2;
    while (i >= 0 && arr[i] >= arr[i + 1]) {
        i--;
    }

    if (i < 0) {
        reverse_subarray(arr, 0, n - 1);
        return 0;
    }

    int j = n - 1;
    while (arr[j] <= arr[i]) {
        j--;
    }

    swap(&arr[i], &arr[j]);

    reverse_subarray(arr, i + 1, n - 1);

    return 1;
}

double get_sq_distance(double* coords1, int id1, double* coords2, int id2) {
    double dx = coords1[3 * id1] - coords2[3 * id2];
    double dy = coords1[3 * id1 + 1] - coords2[3 * id2 + 1];
    double dz = coords1[3 * id1 + 2] - coords2[3 * id2 + 2];

    return dx * dx + dy * dy + dz * dz;
}
