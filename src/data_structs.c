#include "data_structs.h"

int find_parent(DisjointSetUnion* dsu, int x) {
    if (dsu->parent[x] == x) return x;
    return dsu->parent[x] = find_parent(dsu, dsu->parent[x]);
}

int merge(DisjointSetUnion* dsu, int x, int y) {
    x = find_parent(dsu, x);
    y = find_parent(dsu, y);

    if (x != y) {
        dsu->parent[x] = y;
        return 1;
    }

    return 0;
}

DisjointSetUnion* init_disjoint_set_union(int num_atoms, int* look_up_ids, const int** adjacency_list, const int** candidates_data) {
    DisjointSetUnion* dsu = (DisjointSetUnion*)malloc(sizeof(DisjointSetUnion));
    if (!dsu) { fprintf(stderr, "Memory allocation failed for DisjointSetUnion.\n"); return NULL; }

    dsu->size = look_up_ids[0];

    int* map_to_look_up_ids = (int*)malloc(num_atoms * sizeof(int));
    int* kernel_atom = (int*)malloc(num_atoms * sizeof(int));
    dsu->parent = (int*)malloc(look_up_ids[0] * sizeof(int));

    if (!map_to_look_up_ids || !kernel_atom || !dsu->parent) {
        fprintf(stderr, "Memory allocation failed for map_to_look_up_ids or kernel_atom or parent.\n");
        free(map_to_look_up_ids); free(kernel_atom); free(dsu->parent); free(dsu);
        return NULL;
    }

    for (int i = 0; i < num_atoms; i++) {
        kernel_atom[i] = -1;
        map_to_look_up_ids[i] = -1;
    }
    
    for (int i = 0; i < look_up_ids[0]; i++) {
        dsu->parent[i] = i;
        map_to_look_up_ids[look_up_ids[i + 1]] = i;
    }

    for (int i = 1; i <= look_up_ids[0]; i++) {
        int id = look_up_ids[i];

        for (int j = 1; j <= adjacency_list[id][0]; j++) {
            int neighbor_id = get_bonded_atom_id(adjacency_list[id][j]);
            if (map_to_look_up_ids[neighbor_id] != -1)
                merge(dsu, i - 1, map_to_look_up_ids[neighbor_id]);
        }

        for (int j = 1; j <= candidates_data[id][0]; j++) {
            int candidate_id = candidates_data[id][j];
            if (kernel_atom[candidate_id] == -1) kernel_atom[candidate_id] = i - 1;
            else merge(dsu, i - 1, kernel_atom[candidate_id]);
        }
    }

    return dsu;
}

int** get_disjoint_sets_and_free_up(DisjointSetUnion* dsu) {
    int** disjoint_sets = (int**)malloc((dsu->size + 1) * sizeof(int*));
    for (int i = 0; i <= dsu->size; i++) disjoint_sets[i] = NULL;
    
    for (int i = 0; i < dsu->size; i++) find_parent(dsu, i);

    int count = 0;
    for (int i = 0; i < dsu->size; i++) {
        if (dsu->parent[i] != -1) {
            int cnt = 0, root = dsu->parent[i];
            for (int j = i; j < dsu->size; j++)
                if (dsu->parent[j] == root) cnt++;
            
            disjoint_sets[count] = (int*)malloc((cnt + 1) * sizeof(int));
            if (!disjoint_sets[count]) {
                fprintf(stderr, "Memory allocation failed for disjoint_sets[%d].\n", count);
                for (int j = 0; j < count; j++) free(disjoint_sets[j]);
                free(disjoint_sets); free(dsu->parent); free(dsu);
                return NULL;
            }

            for (int j = i, k = 1; j < dsu->size; j++)
                if (dsu->parent[j] == root) {
                    disjoint_sets[count][k++] = j;
                    dsu->parent[j] = -1;
                }
            
            disjoint_sets[count][0] = cnt;
            count++;
        }
    }
    
    disjoint_sets = (int**)realloc(disjoint_sets, (count + 1) * sizeof(int*));
    free(dsu->parent); free(dsu);

    return disjoint_sets;
}
