# Simulated Annealing (C)

**Simulated annealing** is a probabilistic metaheuristic for global optimization that mimics the physical annealing process — heating a material then slowly cooling it to minimize energy. This header-only C library provides a generic, reheat-capable SA solver with three cooling schedules and support for arbitrary state spaces.

## Why It Matters

Many optimization problems (TSP, graph coloring, job-shop scheduling, neural network architecture search) have rugged fitness landscapes riddled with local minima. Gradient-based methods get trapped; exhaustive search is NP-hard. Simulated annealing escapes local minima by *probabilistically accepting worse solutions* early in the search (at high temperature), then becoming increasingly greedy as temperature drops. The method is guaranteed to converge to the global optimum under sufficient iteration counts (Geman & Geman, 1984). This library makes SA accessible to C projects with zero external dependencies — just include `sa.h`.

## How It Works

### The Metropolis Criterion

At each step, SA generates a neighbor state `s'` from current state `s`. If `cost(s') ≤ cost(s)`, the move is always accepted. If worse, it's accepted with probability:

```
P(accept) = exp(-ΔE / T)
```

where `ΔE = cost(s') - cost(s)` and `T` is the current temperature. At high T, P → 1 (accept everything, explore freely). At low T, P → 0 (reject worse moves, exploit greedily).

### Cooling Schedules

Three schedules are provided:

| Schedule | Formula | Behavior |
|----------|---------|----------|
| Linear | T(i) = T₀ - (T₀ - T_min) · i/N | Steady decrease |
| Exponential | T(i) = T₀ · αⁱ | Geometric decrease (most common) |
| Logarithmic | T(i) = T₀ / (1 + α · i) | Slow, theory-backed convergence |

Exponential cooling with α ≈ 0.999 is the practical default. Logarithmic cooling guarantees asymptotic convergence but is impractically slow for production use.

### Reheat Mechanism

When the solver detects stagnation (no improvement for `reheat_interval` iterations), it raises the temperature:

```
if stale_count ≥ 2:
    T = T_init · reheat_factor
    stale_count = 0
```

This kicks the search out of local minima. The `reheat_count` in the result tracks how many reheats occurred — a high count suggests a rugged landscape.

### Complexity

Each iteration is O(state_copy + cost_evaluation). For the sphere function `f(x) = x²`, this is O(1). For TSP with n cities, `cost()` is O(n) and `neighbor()` is O(1) (swap two cities). Total: O(I · n) where I = max_iterations.

### Convergence Guarantee

The Geman-Geman theorem states that SA with logarithmic cooling T(i) = c/log(i+1) converges to the global optimum if c ≥ the depth of the deepest local minimum. In practice, exponential cooling trades the guarantee for speed and works well empirically.

## Quick Start

```c
#include "sa.h"

// Define your state
typedef struct { double x; } State;

// Cost function: minimize x²
double cost(const void *s) {
    return ((const State*)s)->x * ((const State*)s)->x;
}

// Neighbor: perturb by temperature-scaled noise
void neighbor(void *cur, void *next, double temp) {
    State *n = (State*)next;
    n->x = ((State*)cur)->x + (2.0 * rand() / RAND_MAX - 1.0) * sqrt(temp) * 0.5;
}

int main(void) {
    State initial = { 10.0 }; // Start far from optimum

    sa_config cfg = {
        .cost = cost,
        .neighbor = neighbor,
        .state_size = sizeof(State),
        .temp_init = 100.0,
        .temp_min = 0.001,
        .cooling_rate = 0.999,
        .max_iterations = 10000,
        .schedule = SA_SCHEDULE_EXPONENTIAL,
        .reheat_interval = 500,
        .reheat_factor = 0.5,
    };

    sa_result r = sa_run(&cfg, &initial);
    State *best = (State*)r.best_state;

    printf("Best x = %.6f, cost = %.6f\n", best->x, r.best_cost);
    printf("Iterations: %d, reheats: %d\n", r.iterations, r.reheat_count);

    sa_result_free(&r);
    return 0;
}
```

Compile: `gcc -lm -o demo test.c && ./demo`

## API

| Type | Field/Method | Description |
|------|-------------|-------------|
| `sa_config` | `cost, neighbor` | Function pointers for cost and neighbor generation |
| `sa_config` | `temp_init, temp_min, cooling_rate` | Temperature parameters |
| `sa_config` | `schedule` | LINEAR, EXPONENTIAL, or LOGARITHMIC |
| `sa_config` | `reheat_interval, reheat_factor` | Reheat on stagnation |
| `sa_result` | `best_state, best_cost` | Optimal solution found |
| `sa_result` | `cost_history[]` | Best cost per iteration |
| `sa_run()` | `(cfg, initial_state) → result` | Execute SA |
| `sa_result_free()` | `(result*)` | Free allocated memory |

## Architecture Notes

Simulated Annealing provides the global optimization backbone for η (eta) in SuperInstance — the process of *eliminating* suboptimal configurations. The Metropolis criterion's acceptance of worse solutions mirrors the negative-space principle: sometimes you must move *toward* a worse state to eventually escape it. The cooling schedule controls the γ/η balance — high temperature favors γ (exploration, construction), low temperature favors η (exploitation, elimination). The reheat mechanism is a meta-level η: when the search itself stagnates, the algorithm resets to re-examine its assumptions. See [ARCHITECTURE.md](https://github.com/SuperInstance/SuperInstance/blob/main/ARCHITECTURE.md).

## References

1. Kirkpatrick, S., Gelatt, C. D., & Vecchi, M. P. (1983). "Optimization by Simulated Annealing." *Science*, 220(4598), 671–680.
2. Geman, S., & Geman, D. (1984). "Stochastic Relaxation, Gibbs Distributions, and the Bayesian Restoration of Images." *IEEE TPAMI*, 6(6), 721–741. — Convergence proof.
3. van Laarhoven, P. J. M., & Aarts, E. H. L. (1987). *Simulated Annealing: Theory and Applications*. D. Reidel.

## License

MIT
