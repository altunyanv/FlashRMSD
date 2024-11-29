#ifndef DATA_STRUCTS_H
#define DATA_STRUCTS_H

#include <stdlib.h>
#include <stdio.h>

#include "common.h"

typedef struct {
    int size;
    int* parent;
} DisjointSetUnion;

int find_parent(DisjointSetUnion*, int);
int merge(DisjointSetUnion*, int, int);

DisjointSetUnion* init_disjoint_set_union(int, int*, const int**, const int**);
int** get_disjoint_sets_and_free_up(DisjointSetUnion*);

#endif // DATA_STRUCTS_H