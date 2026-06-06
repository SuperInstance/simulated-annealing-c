#include "sa.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>

/* Test 1: Minimize sphere function f(x) = x^2 */
typedef struct { double x; } sphere_state;

double sphere_cost(const void *s) {
    return pow(((const sphere_state*)s)->x, 2);
}

void sphere_neighbor(void *s, void *n, double temp) {
    sphere_state *next = (sphere_state*)n;
    next->x = ((sphere_state*)s)->x + (2.0 * rand() / RAND_MAX - 1.0) * sqrt(temp) * 0.5;
}

void test_sphere() {
    sphere_state init = {10.0}; // Start far from optimum
    sa_config cfg = {
        .cost = sphere_cost,
        .neighbor = sphere_neighbor,
        .state_size = sizeof(sphere_state),
        .temp_init = 100.0,
        .temp_min = 0.001,
        .cooling_rate = 0.999,
        .max_iterations = 10000,
        .schedule = SA_SCHEDULE_EXPONENTIAL,
        .reheat_interval = 0,
        .reheat_factor = 1.0
    };
    sa_result r = sa_run(&cfg, &init);
    sphere_state *best = (sphere_state*)r.best_state;
    printf("Sphere: x=%.6f, cost=%.6f (expected ~0)\n", best->x, r.best_cost);
    assert(fabs(best->x) < 1.0);
    assert(r.best_cost < 1.0);
    sa_result_free(&r);
    printf("  PASS: sphere function\n");
}

/* Test 2: Rastrigin function */
double rastrigin_cost(const void *s) {
    double x = ((const sphere_state*)s)->x;
    return x*x - 10*cos(2*M_PI*x) + 10;
}

void rastrigin_neighbor(void *s, void *n, double temp) {
    sphere_state *next = (sphere_state*)n;
    next->x = ((sphere_state*)s)->x + (2.0 * rand() / RAND_MAX - 1.0) * sqrt(temp) * 0.3;
}

void test_rastrigin() {
    sphere_state init = {5.0};
    sa_config cfg = {
        .cost = rastrigin_cost,
        .neighbor = rastrigin_neighbor,
        .state_size = sizeof(sphere_state),
        .temp_init = 200.0,
        .temp_min = 0.01,
        .cooling_rate = 0.9995,
        .max_iterations = 20000,
        .schedule = SA_SCHEDULE_EXPONENTIAL,
        .reheat_interval = 0,
        .reheat_factor = 1.0
    };
    sa_result r = sa_run(&cfg, &init);
    printf("Rastrigin: x=%.4f, cost=%.4f (expected ~0)\n", 
           ((sphere_state*)r.best_state)->x, r.best_cost);
    assert(r.best_cost < 5.0); // Should get close to a minimum
    sa_result_free(&r);
    printf("  PASS: rastrigin function\n");
}

/* Test 3: TSP (small, 5 cities) */
typedef struct { int cities[5]; int n; } tsp_state;

double tsp_cost(const void *s) {
    const tsp_state *t = (const tsp_state*)s;
    // Simple 5-city distances (predefined coordinates)
    double coords[5][2] = {{0,0},{1,0},{1,1},{0,1},{0.5,0.5}};
    double total = 0;
    for (int i = 0; i < t->n; i++) {
        int a = t->cities[i];
        int b = t->cities[(i+1) % t->n];
        total += sqrt(pow(coords[a][0]-coords[b][0],2) + pow(coords[a][1]-coords[b][1],2));
    }
    return total;
}

void tsp_neighbor(void *s, void *n, double temp) {
    tsp_state *next = (tsp_state*)n;
    *next = *(tsp_state*)s;
    int i = rand() % next->n;
    int j = rand() % next->n;
    int tmp = next->cities[i];
    next->cities[i] = next->cities[j];
    next->cities[j] = tmp;
}

void test_tsp() {
    tsp_state init = {{0,1,2,3,4}, 5};
    sa_config cfg = {
        .cost = tsp_cost,
        .neighbor = tsp_neighbor,
        .state_size = sizeof(tsp_state),
        .temp_init = 50.0,
        .temp_min = 0.01,
        .cooling_rate = 0.999,
        .max_iterations = 20000,
        .schedule = SA_SCHEDULE_EXPONENTIAL,
        .reheat_interval = 0,
        .reheat_factor = 1.0
    };
    sa_result r = sa_run(&cfg, &init);
    printf("TSP: cost=%.4f (expected <5)\n", r.best_cost);
    assert(r.best_cost < 5.0);
    sa_result_free(&r);
    printf("  PASS: TSP 5-city\n");
}

/* Test 4: Temperature schedules */
void test_schedules() {
    sa_config cfg = {
        .temp_init = 100.0,
        .temp_min = 0.1,
        .cooling_rate = 0.99,
        .max_iterations = 1000,
    };
    
    // Linear should reach exactly temp_min at end
    cfg.schedule = SA_SCHEDULE_LINEAR;
    double t_linear = sa_temperature(&cfg, 1000);
    assert(fabs(t_linear - 0.1) < 1e-10);
    printf("  PASS: linear schedule\n");
    
    // Exponential at step 0 should be temp_init
    cfg.schedule = SA_SCHEDULE_EXPONENTIAL;
    double t_exp = sa_temperature(&cfg, 0);
    assert(fabs(t_exp - 100.0) < 1e-10);
    printf("  PASS: exponential schedule\n");
    
    // Logarithmic should decrease
    cfg.schedule = SA_SCHEDULE_LOGARITHMIC;
    double t_log_0 = sa_temperature(&cfg, 0);
    double t_log_100 = sa_temperature(&cfg, 100);
    assert(t_log_0 > t_log_100);
    printf("  PASS: logarithmic schedule\n");
}

/* Test 5: Reheat mechanism */
void test_reheat() {
    sphere_state init = {8.0};
    sa_config cfg = {
        .cost = sphere_cost,
        .neighbor = sphere_neighbor,
        .state_size = sizeof(sphere_state),
        .temp_init = 10.0,
        .temp_min = 0.001,
        .cooling_rate = 0.999,
        .max_iterations = 5000,
        .schedule = SA_SCHEDULE_EXPONENTIAL,
        .reheat_interval = 500,
        .reheat_factor = 0.5
    };
    sa_result r = sa_run(&cfg, &init);
    printf("Reheat: x=%.6f, reheats=%d\n", ((sphere_state*)r.best_state)->x, r.reheat_count);
    assert(r.iterations == 5000);
    sa_result_free(&r);
    printf("  PASS: reheat mechanism\n");
}

/* Test 6: Cost history tracking */
void test_history() {
    sphere_state init = {3.0};
    sa_config cfg = {
        .cost = sphere_cost,
        .neighbor = sphere_neighbor,
        .state_size = sizeof(sphere_state),
        .temp_init = 10.0,
        .temp_min = 0.001,
        .cooling_rate = 0.999,
        .max_iterations = 100,
        .schedule = SA_SCHEDULE_EXPONENTIAL,
        .reheat_interval = 0,
        .reheat_factor = 1.0
    };
    sa_result r = sa_run(&cfg, &init);
    assert(r.history_len == 101);
    assert(r.cost_history[0] == 9.0); // 3^2
    // Best cost should be non-increasing
    for (int i = 1; i < r.history_len; i++) {
        assert(r.cost_history[i] <= r.cost_history[i-1] + 1e-10);
    }
    sa_result_free(&r);
    printf("  PASS: cost history\n");
}

/* Test 7: Multi-start (manual) */
void test_multistart() {
    double best_cost = 1e9;
    sphere_state best_state = {0};
    
    for (int start = 0; start < 5; start++) {
        sphere_state init = {(double)(start * 3 - 6)};
        sa_config cfg = {
            .cost = sphere_cost,
            .neighbor = sphere_neighbor,
            .state_size = sizeof(sphere_state),
            .temp_init = 50.0,
            .temp_min = 0.001,
            .cooling_rate = 0.999,
            .max_iterations = 5000,
            .schedule = SA_SCHEDULE_EXPONENTIAL,
            .reheat_interval = 0,
            .reheat_factor = 1.0
        };
        sa_result r = sa_run(&cfg, &init);
        if (r.best_cost < best_cost) {
            best_cost = r.best_cost;
            best_state = *(sphere_state*)r.best_state;
        }
        sa_result_free(&r);
    }
    printf("Multi-start: best x=%.6f, cost=%.6f\n", best_state.x, best_cost);
    assert(fabs(best_state.x) < 0.5);
    printf("  PASS: multi-start\n");
}

int main() {
    srand(42);
    printf("=== Simulated Annealing C Tests ===\n\n");
    
    test_sphere();
    test_rastrigin();
    test_tsp();
    test_schedules();
    test_reheat();
    test_history();
    test_multistart();
    
    printf("\nAll 7 test groups passed! ✅\n");
    return 0;
}
