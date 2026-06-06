# simulated-annealing-c

Generic simulated annealing framework in C with multiple cooling schedules, automatic reheating, and multi-start support.

## Features

- **Generic API** — void* state with pluggable cost and neighbor functions
- **Cooling Schedules** — Linear, exponential, and logarithmic
- **Automatic Reheating** — Detect stagnation and reheat to escape local minima
- **Cost History** — Track best cost per iteration for convergence analysis
- **Multi-Start** — Run multiple independent trials, keep best result
- **Tested Problems** — Sphere function, Rastrigin function, TSP (5-city)

## Usage

```c
#include "sa.h"

// Define your state
typedef struct { double x; } MyState;

// Cost function: minimize x^2
double my_cost(const void *state) {
    return pow(((const MyState*)state)->x, 2);
}

// Neighbor function: perturb by random amount scaled by temperature
void my_neighbor(void *state, void *next, double temp) {
    MyState *s = (MyState*)next;
    s->x = ((MyState*)state)->x + (2.0 * rand() / RAND_MAX - 1.0) * sqrt(temp) * 0.5;
}

int main() {
    srand(42);
    MyState init = {10.0};  // Start far from optimum

    sa_config cfg = {
        .cost = my_cost,
        .neighbor = my_neighbor,
        .state_size = sizeof(MyState),
        .temp_init = 100.0,
        .temp_min = 0.001,
        .cooling_rate = 0.999,
        .max_iterations = 10000,
        .schedule = SA_SCHEDULE_EXPONENTIAL,
        .reheat_interval = 0,
        .reheat_factor = 1.0
    };

    sa_result result = sa_run(&cfg, &init);
    printf("Best x = %f, cost = %f\n",
           ((MyState*)result.best_state)->x,
           result.best_cost);

    sa_result_free(&result);
    return 0;
}
```

## Building & Testing

```bash
make test    # Builds and runs 7 test groups
```

### Test Results

| Test | Description | Result |
|------|-------------|--------|
| Sphere | Minimize x² from x=10 | ✅ finds x≈0 |
| Rastrigin | Multi-modal from x=5 | ✅ finds x≈0 |
| TSP 5-city | Traveling salesman | ✅ cost < 5 |
| Schedules | Linear, exponential, logarithmic | ✅ correct |
| Reheat | Auto-reheat on stagnation | ✅ |
| History | Cost tracking, monotonically decreasing | ✅ |
| Multi-start | 5 independent runs | ✅ all converge |

## License

MIT
