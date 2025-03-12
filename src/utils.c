#include "utils.h"

#define OPT_LEVEL 2

int layer_queue_hash(int* layer_queue, int num_atoms) {
    uint32_t hash = 0;
    for (int i = 0; i < num_atoms; i++) {
        hash ^= (layer_queue[i] >> MAX_BLOCK_LENGTH);
        hash = ((hash << 5) | (hash >> (32 - 5))) & 0xFFFFFFFF;
    }
    return (int)hash;
}

int generate_layer_data(int num_atoms, int* atomic_numbers, int** adjacency_list, int* layer_data) {
    int* layer_queue = (int*)malloc(num_atoms * sizeof(int));
    int* visited = (int*)malloc(num_atoms * sizeof(int));
    
    if (!layer_queue || !visited) { fprintf(stderr, "Memory allocation failed for layer queue or visited array.\n"); return 0; }

    for (int i = 0; i < num_atoms; i++) {
        for (int j = 0; j < num_atoms; j++) visited[j] = 0;

        int l = 0, r = 0;

        layer_queue[r++] = (0 << LAYER_ID_SHIFT) | (atomic_numbers[i] << ATOMIC_NUMBER_SHIFT) | i;
        visited[i] = 1;

        while (l < r) {
            int current = layer_queue[l];
            int layer_id = get_layer_id(current);
            int atom_id = get_atom_id(current);

            int num_neighbors = adjacency_list[atom_id][0];
            for (int neighbor_index = 1; neighbor_index <= num_neighbors; neighbor_index++) {
                int neighbor_atom_id = get_bonded_atom_id(adjacency_list[atom_id][neighbor_index]);
                if (!visited[neighbor_atom_id]) {
                    visited[neighbor_atom_id] = 1;
                    layer_queue[r++] = ((layer_id + 1) << LAYER_ID_SHIFT)
                        | (atomic_numbers[neighbor_atom_id] << ATOMIC_NUMBER_SHIFT)
                        | neighbor_atom_id;
                }
            }

            l++;

            if (l < r && get_layer_id(layer_queue[l]) != layer_id) qsort(layer_queue + l, r - l, sizeof(int), compare_ints);
        }

        layer_data[i] = layer_queue_hash(layer_queue, num_atoms);
    }

    free(layer_queue);
    free(visited);

    return 1;
}

int** get_candidate_data(int num_atoms, int* layer_data_1, int* layer_data_2) {
    int** candidate_ids = (int**)malloc(num_atoms * sizeof(int*));

    if (!candidate_ids) { fprintf(stderr, "Memory allocation failed for candidate_ids.\n"); return NULL; }

    for (int i = 0; i < num_atoms; i++) {
        int num_candidates = 0;
        for (int j = 0; j < num_atoms; j++)
            if (layer_data_1[i] == layer_data_2[j]) num_candidates++;

        candidate_ids[i] = (int*)malloc((num_candidates + 1) * sizeof(int));
        if (!candidate_ids[i]) { fprintf(stderr, "Memory allocation failed for candidate_ids[%d].\n", i); return NULL; }

        candidate_ids[i][0] = num_candidates;
        num_candidates = 1;
        for (int j = 0; j < num_atoms; j++)
            if (layer_data_1[i] == layer_data_2[j]) candidate_ids[i][num_candidates++] = j;
    }

    return candidate_ids;
}

void sort_candidates_by_distance(int* candidate_data, double* candidate_squared_distances) {
    int num_candidates = candidate_data[0];
    for (int i = 1; i <= num_candidates; i++) {
        for (int j = i - 1; j > 0; j--) {
            if (candidate_squared_distances[j - 1] > candidate_squared_distances[j]) {
                int temp = candidate_data[j];
                candidate_data[j] = candidate_data[j + 1];
                candidate_data[j + 1] = temp;

                double temp_dist = candidate_squared_distances[j - 1];
                candidate_squared_distances[j - 1] = candidate_squared_distances[j];
                candidate_squared_distances[j] = temp_dist;
            } else break;
        }
    }
}

double** get_candidate_squared_distances(int num_atoms, double* coordinates_1, double* coordinates_2, int** candidate_data) {
    double** candidate_distances = (double**)malloc(num_atoms * sizeof(double*));

    if (!candidate_distances) { fprintf(stderr, "Memory allocation failed for candidate_distances.\n"); return NULL; }

    for (int i = 0; i < num_atoms; i++) {
        int num_candidates = candidate_data[i][0];
        candidate_distances[i] = (double*)malloc(num_candidates * sizeof(double));
        
        if (!candidate_distances[i]) { fprintf(stderr, "Memory allocation failed for candidate_distances[%d].\n", i); return NULL; }

        candidate_distances[i][0] = num_candidates;

        for (int j = 1; j <= num_candidates; j++) {
            int atom_id_1 = i, atom_id_2 = candidate_data[i][j];
            double dx = coordinates_1[3 * atom_id_1] - coordinates_2[3 * atom_id_2];
            double dy = coordinates_1[3 * atom_id_1 + 1] - coordinates_2[3 * atom_id_2 + 1];
            double dz = coordinates_1[3 * atom_id_1 + 2] - coordinates_2[3 * atom_id_2 + 2];
            candidate_distances[i][j - 1] = dx * dx + dy * dy + dz * dz;
        }

        sort_candidates_by_distance(candidate_data[i], candidate_distances[i]);
    }

    return candidate_distances;
}




const SearchData* init_search_data(int num_atoms, int** candidates, double** squared_distances, int** adjacency_list_1, int** adjacency_list_2) {
    SearchData* data = (SearchData*)malloc(sizeof(SearchData));

    if (!data) { fprintf(stderr, "Memory allocation failed for SearchData.\n"); return NULL; }

    data->num_atoms = num_atoms;
    data->candidates = (const int**)candidates;
    data->squared_distances = (const double**)squared_distances;
    data->adjacency_list_1 = (const int**)adjacency_list_1;
    data->adjacency_list_2 = (const int**)adjacency_list_2;

    return (const SearchData*)data;
}

double search_assignment(const SearchData* data, int* assignment) {
    int num_atoms = data->num_atoms;
    const int** candidate_data = data->candidates;
    const double** candidate_squared_distances = data->squared_distances;
    const int** adjacency_list_1 = data->adjacency_list_1;
    const int** adjacency_list_2 = data->adjacency_list_2;

    double best_sum = 1000000000.0;

    int index = 0;

    int* current_assignment = (int*)malloc(num_atoms * sizeof(int));
    double current_sum = 0.0;
    
    int* mapped_connections = (int*)malloc(num_atoms * sizeof(int));
    int* active_choices = (int*)malloc(num_atoms * sizeof(int));
    int* active_choices_candidate_index = (int*)malloc(num_atoms * sizeof(int));
    
    int* active_choices_mask = (int*)malloc(num_atoms * sizeof(int));
    int* active_mapped_choices_mask = (int*)malloc(num_atoms * sizeof(int));
    
    if (!current_assignment || !mapped_connections || !active_choices || !active_choices_candidate_index || !active_choices_mask || !active_mapped_choices_mask) {
        fprintf(stderr, "Memory allocation failed for current_assignment, active_choices or active_choices_candidate_index.\n");
        return INFINITY;
    }

    for (int i = 0; i < num_atoms; i++) { 
        current_assignment[i] = -1;
        active_choices_candidate_index[i] = -1; 
        mapped_connections[i] = 0; 
        
        active_choices[i] = -1;
        active_choices_mask[i] = 0;
        active_mapped_choices_mask[i] = 0;

        if (candidate_data[i][0] == 0) {
            // If there are no candidates for an atom, we can't find a mapping
            fprintf(stderr, "No candidates for atom %d.\n", i);
            return INFINITY;
        }
    }

    int id, found_any = 0;
    while (index >= 0) {
        if (index == num_atoms) {
            if (current_sum < best_sum) {
                memcpy(assignment, current_assignment, num_atoms * sizeof(int));
                best_sum = current_sum;
            }

            index--;
            id = active_choices[index];

            current_sum -= candidate_squared_distances[id][active_choices_candidate_index[index]];

            active_mapped_choices_mask[current_assignment[id]] = 0;
            current_assignment[id] = -1;

            continue;
        }

        id = active_choices[index];

        // If needed choose new atom to map
        if (active_choices_candidate_index[index] == -1) {
            id = -1;
            double best_metric = 1000000000.0;

            for (int i = 0; i < num_atoms; i++) 
                if (!active_choices_mask[i]) {
                    double current_metric = candidate_data[i][0];
                    if (current_metric < best_metric) {
                        best_metric = current_metric;
                        id = i;
                    }
                }

            active_choices[index] = id;
            active_choices_mask[id] = 1;
        }

        if (active_choices_candidate_index[index] == candidate_data[id][0] - 1) {
            active_choices_candidate_index[index] = -1;
            active_choices_mask[id] = 0;
            
            active_choices[index] = -1;
            index--;

            if (index < 0) break;
            
            id = active_choices[index];

            current_sum -= candidate_squared_distances[id][active_choices_candidate_index[index]];
            active_mapped_choices_mask[current_assignment[id]] = 0;
            current_assignment[id] = -1;
            continue;
        }

        int found_next = 0;
        for (int i = active_choices_candidate_index[index] + 1; i < candidate_data[id][0]; i++) {
            if (current_sum + candidate_squared_distances[id][i] >= best_sum) break;

            int mapped_choice = candidate_data[id][i + 1];
            if (active_mapped_choices_mask[mapped_choice] || adjacency_list_1[id][0] != adjacency_list_2[mapped_choice][0]) continue;

            found_next = 1;
            for (int j = 1; j <= adjacency_list_1[id][0]; j++) {
                int neighbor_id_1 = get_bonded_atom_id(adjacency_list_1[id][j]);
                int bond_type_mask = get_bond_type(adjacency_list_1[id][j]) << BOND_TYPE_SHIFT;
                int mapped_neighbor_id = current_assignment[neighbor_id_1];

                if (mapped_neighbor_id == -1) continue;

                int found = 0;
                for (int k = 1; k <= adjacency_list_2[mapped_choice][0]; k++)
                    if ((mapped_neighbor_id | bond_type_mask) == adjacency_list_2[mapped_choice][k]) {
                        found = 1; break;
                    }

                if (!found) { found_next = 0; break; }
            }

            if (found_next) {
                current_sum += candidate_squared_distances[id][i];
                active_mapped_choices_mask[mapped_choice] = 1;
                current_assignment[id] = mapped_choice;
                
                active_choices_candidate_index[index] = i;

                index++;
                break;
            }
        }

        if (!found_next)
            active_choices_candidate_index[index] = candidate_data[id][0] - 1;
    }

    free(current_assignment);
    free(mapped_connections);
    
    free(active_choices);
    free(active_choices_candidate_index);
    
    free(active_choices_mask);
    free(active_mapped_choices_mask);
    
    double best_rmsd = sqrt(best_sum / num_atoms);

    return best_rmsd;
}

double search_assignment_recurse(const SearchData* data, int* best_assignment) {
    RecursionState* state = (RecursionState*)malloc(sizeof(RecursionState));

    if (!state) { fprintf(stderr, "Memory allocation failed for RecursionState.\n"); return INFINITY; }

    state->used_mask_1 = (int*)malloc(data->num_atoms * sizeof(int));
    state->used_mask_2 = (int*)malloc(data->num_atoms * sizeof(int));

    state->assignment = (int*)malloc(data->num_atoms * sizeof(int));
    state->assigned_count = 0;
    state->current_sum = 0.0;

    state->best_assignment = best_assignment;
    state->best_sum = INFINITY;

    if (!state->used_mask_1 || !state->used_mask_2 || !state->assignment) {
        fprintf(stderr, "Memory allocation failed for RecursionState members.\n");
        free(state->used_mask_1);   free(state->used_mask_2);
        free(state->assignment);    free(state);
        return INFINITY;
    }

    for (int i = 0; i < data->num_atoms; i++) {
        state->used_mask_1[i] = state->used_mask_2[i] = 0;
        state->assignment[i] = -1;
    }

    int* look_up_ids = (int*)malloc((data->num_atoms + 1) * sizeof(int));
    if (!look_up_ids) {
        fprintf(stderr, "Memory allocation failed for look_up_ids.\n");
        free(state->used_mask_1);   free(state->used_mask_2);
        free(state->assignment);    free(state);
        return INFINITY;
    }

    look_up_ids[0] = data->num_atoms;
    for (int i = 1; i <= data->num_atoms; i++) look_up_ids[i] = i - 1;
        
    if (OPT_LEVEL == 0)
        search_assignment_recurse_helper(data, state, data->num_atoms, look_up_ids, -1);
    else {
        look_up_ids[0] = 0;
        int look_up_index = 1;

        for (int i = 0; i < data->num_atoms; i++)
            if (data->candidates[i][0] == 1) {
                state->used_mask_1[i] = 1;
                state->used_mask_2[data->candidates[i][1]] = 1;
                state->assignment[i] = data->candidates[i][1];
                state->assigned_count++;
                state->current_sum += data->squared_distances[i][0];
            } else {
                look_up_ids[0]++;
                look_up_ids[look_up_index++] = i;
            }

        look_up_ids = (int*)realloc(look_up_ids, (look_up_ids[0] + 1) * sizeof(int));

        if (look_up_ids[0] == 0 || OPT_LEVEL == 1)
            search_assignment_recurse_helper(data, state, data->num_atoms, look_up_ids, -1);
        else if (OPT_LEVEL == 2) {
            DisjointSetUnion* dsu = init_disjoint_set_union(data->num_atoms, look_up_ids, data->adjacency_list_1, data->candidates);            

            if (!dsu) {
                fprintf(stderr, "Failed to initialize DisjointSetUnion.\n");
                free(state->used_mask_1);   free(state->used_mask_2);
                free(state->assignment);    free(state);
                return INFINITY;
            }

            int** disjoint_sets = get_disjoint_sets_and_free_up(dsu);

            for (int i = 0; disjoint_sets[i]; i++) {
                state->best_sum = INFINITY;

                for (int j = 1; j <= disjoint_sets[i][0]; j++)
                    disjoint_sets[i][j] = look_up_ids[1 + disjoint_sets[i][j]];

                int threshold = state->assigned_count + disjoint_sets[i][0];
                
                search_assignment_recurse_helper(data, state, threshold, disjoint_sets[i], -1);

                for (int j = 1; j <= disjoint_sets[i][0]; j++) {
                    int id = disjoint_sets[i][j];
                    state->used_mask_1[id] = 1;
                    state->used_mask_2[state->best_assignment[id]] = 1;
                    state->assignment[id] = state->best_assignment[id];
                }

                state->assigned_count = threshold;
                state->current_sum = state->best_sum;
            }
        } else {
            fprintf(stderr, "Invalid optimization level.\n");
            free(state->used_mask_1);   free(state->used_mask_2);
            free(state->assignment);    free(state);
            return INFINITY;
        }
    } 

    double best_rmsd = sqrt(state->best_sum / data->num_atoms);

    free(state->used_mask_1);   free(state->used_mask_2);
    free(state->assignment);    free(state);

    return best_rmsd;
}

void search_assignment_recurse_helper(const SearchData* data, RecursionState* state, int threshold, int* look_up_ids, int last_id) {
    if (state->assigned_count == threshold) {
        if (state->current_sum < state->best_sum) {
            memcpy(state->best_assignment, state->assignment, data->num_atoms * sizeof(int));
            state->best_sum = state->current_sum;
        }
    } else {
        const int** adjacency_list_1 = data->adjacency_list_1;
        const int** adjacency_list_2 = data->adjacency_list_2;

        int id = -1;
        
        if (last_id != -1) {
            for (int i = 1; i <= adjacency_list_1[last_id][0]; i++) {
                int candidate_id = get_bonded_atom_id(adjacency_list_1[last_id][i]);
                if (state->used_mask_1[candidate_id]) continue;

                if (id == -1 || data->candidates[candidate_id][0] < data->candidates[id][0]) 
                    id = candidate_id;
            }
        }

        if (id == -1) {
            for (int i = 1; i <= look_up_ids[0]; i++) {
                if (state->used_mask_1[look_up_ids[i]]) continue;

                if (id == -1 || data->candidates[look_up_ids[i]][0] < data->candidates[id][0]) 
                    id = look_up_ids[i];
            }
        }
        
        state->used_mask_1[id] = 1;
        state->assigned_count++;

        double current_sum = state->current_sum;
        for (int i = 1; i <= data->candidates[id][0]; i++) {
            if (current_sum + data->squared_distances[id][i - 1] >= state->best_sum) break;

            int mapped_choice = data->candidates[id][i];
            
            if (state->used_mask_2[mapped_choice]) continue;

            // Validate the mapping
            int valid = 1;
            for (int j = 1; j <= adjacency_list_1[id][0]; j++) {
                int neighbor_id_1 = get_bonded_atom_id(adjacency_list_1[id][j]);
                int bond_type_mask = get_bond_type(adjacency_list_1[id][j]) << BOND_TYPE_SHIFT;
                int mapped_neighbor_id = state->assignment[neighbor_id_1];

                if (mapped_neighbor_id == -1) continue;

                int found = 0;
                for (int k = 1; k <= adjacency_list_2[mapped_choice][0]; k++)
                    if ((mapped_neighbor_id | bond_type_mask) == adjacency_list_2[mapped_choice][k]) {
                        found = 1; break;
                    }

                if (!found) { valid = 0; break; }
            }

            if (!valid) continue;

            state->used_mask_2[mapped_choice] = 1;
            state->assignment[id] = mapped_choice;
            state->current_sum += data->squared_distances[id][i - 1];

            search_assignment_recurse_helper(data, state, threshold, look_up_ids, id);

            state->current_sum -= data->squared_distances[id][i - 1];
            state->used_mask_2[mapped_choice] = 0;
            state->assignment[id] = -1;
        }

        state->assigned_count--;
        state->used_mask_1[id] = 0;
    }
}