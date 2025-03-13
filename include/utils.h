#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#include "common.h"
#include "data_structs.h"


int layer_queue_hash(int*, int);
int generate_layer_data(int, int*, int**, int*);

int** get_candidate_data(int, int*, int*);
double** get_candidate_squared_distances(int, double*, double*, int**);

typedef struct {
    int num_atoms;
    const int** candidates;
    const double** squared_distances;
    
    const int** adjacency_list_1;
    const int** adjacency_list_2;
} SearchData;

const SearchData* init_search_data(int, int**, double**, int**, int**);

typedef struct {
    int* used_mask_1;
    int* used_mask_2;

    int* assignment;
    int assigned_count;
    double current_sum;

    int* best_assignment;
    double best_sum;
} RecursionState;

double search_assignment(const SearchData*, int*);

double search_assignment_recurse(const SearchData*, int*);
void search_assignment_recurse_helper(const SearchData*, RecursionState*, int, int*, int);

#endif