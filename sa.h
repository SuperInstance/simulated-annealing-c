#ifndef SA_H
#define SA_H

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef double (*sa_cost_fn)(const void *state);
typedef void (*sa_neighbor_fn)(void *state, void *next, double temp);
typedef int (*sa_accept_fn)(double old_cost, double new_cost, double temp);

typedef enum {
    SA_SCHEDULE_LINEAR,
    SA_SCHEDULE_EXPONENTIAL,
    SA_SCHEDULE_LOGARITHMIC
} sa_schedule_type;

typedef struct {
    sa_cost_fn cost;
    sa_neighbor_fn neighbor;
    size_t state_size;
    double temp_init;
    double temp_min;
    double cooling_rate;
    int max_iterations;
    sa_schedule_type schedule;
    int reheat_interval;
    double reheat_factor;
} sa_config;

typedef struct {
    void *best_state;
    double best_cost;
    void *final_state;
    double final_cost;
    int iterations;
    int reheat_count;
    double *cost_history;
    int history_len;
    double convergence_temp;
} sa_result;

static double sa_temperature(const sa_config *cfg, int iteration) {
    double t;
    switch (cfg->schedule) {
        case SA_SCHEDULE_LINEAR:
            t = cfg->temp_init - (cfg->temp_init - cfg->temp_min) * 
                ((double)iteration / cfg->max_iterations);
            break;
        case SA_SCHEDULE_EXPONENTIAL:
            t = cfg->temp_init * pow(cfg->cooling_rate, iteration);
            break;
        case SA_SCHEDULE_LOGARITHMIC:
            t = cfg->temp_init / (1.0 + cfg->cooling_rate * iteration);
            break;
        default:
            t = cfg->temp_init;
    }
    return t < cfg->temp_min ? cfg->temp_min : t;
}

static int sa_default_accept(double old_cost, double new_cost, double temp) {
    if (new_cost <= old_cost) return 1;
    if (temp <= 0) return 0;
    double delta = new_cost - old_cost;
    double prob = exp(-delta / temp);
    return ((double)rand() / RAND_MAX) < prob;
}

static sa_result sa_run(const sa_config *cfg, const void *initial_state) {
    sa_result result = {0};
    
    void *current = malloc(cfg->state_size);
    void *next = malloc(cfg->state_size);
    result.best_state = malloc(cfg->state_size);
    result.final_state = malloc(cfg->state_size);
    
    memcpy(current, initial_state, cfg->state_size);
    double current_cost = cfg->cost(current);
    
    memcpy(result.best_state, current, cfg->state_size);
    result.best_cost = current_cost;
    
    result.cost_history = malloc(sizeof(double) * (cfg->max_iterations + 1));
    result.cost_history[0] = current_cost;
    result.history_len = 1;
    
    int stale_count = 0;
    double last_best = result.best_cost;
    
    for (int i = 0; i < cfg->max_iterations; i++) {
        double temp = sa_temperature(cfg, i);
        
        // Reheat if stuck
        if (cfg->reheat_interval > 0 && i > 0 && i % cfg->reheat_interval == 0) {
            if (fabs(result.best_cost - last_best) < 1e-10) {
                stale_count++;
                if (stale_count >= 2) {
                    // Reheat
                    temp = cfg->temp_init * cfg->reheat_factor;
                    stale_count = 0;
                    result.reheat_count++;
                }
            } else {
                stale_count = 0;
            }
            last_best = result.best_cost;
        }
        
        cfg->neighbor(current, next, temp);
        double next_cost = cfg->cost(next);
        
        if (sa_default_accept(current_cost, next_cost, temp)) {
            memcpy(current, next, cfg->state_size);
            current_cost = next_cost;
        }
        
        if (current_cost < result.best_cost) {
            memcpy(result.best_state, current, cfg->state_size);
            result.best_cost = current_cost;
        }
        
        result.cost_history[result.history_len++] = result.best_cost;
        result.iterations = i + 1;
    }
    
    memcpy(result.final_state, current, cfg->state_size);
    result.final_cost = current_cost;
    result.convergence_temp = sa_temperature(cfg, cfg->max_iterations);
    
    free(current);
    free(next);
    
    return result;
}

static void sa_result_free(sa_result *r) {
    if (r->best_state) { free(r->best_state); r->best_state = NULL; }
    if (r->final_state) { free(r->final_state); r->final_state = NULL; }
    if (r->cost_history) { free(r->cost_history); r->cost_history = NULL; }
}

#endif
