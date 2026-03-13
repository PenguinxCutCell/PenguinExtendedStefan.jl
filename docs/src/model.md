# Mathematical Model

Unknowns:
- liquid temperature `T_l`
- solid temperature `T_s`
- liquid concentration `C_l`
- solid concentration `C_s`
- interface speed `V_\Gamma`
- level set `\phi`

Normal convention: `n` points from solid to liquid.

## Bulk Equations

```math
\rho_l c_{p,l}\,\partial_t T_l = \nabla\cdot(\kappa_l\nabla T_l) + s_{T_l}
```

```math
\rho_s c_{p,s}\,\partial_t T_s = \nabla\cdot(\kappa_s\nabla T_s) + s_{T_s}
```

```math
\partial_t C_l = \nabla\cdot(D_l\nabla C_l) + s_{C_l}
```

```math
\partial_t C_s = \nabla\cdot(D_s\nabla C_s) + s_{C_s}
```

## Interface Conditions

```math
T_{l\Gamma} - T_{s\Gamma} = 0
```

```math
T_{l\Gamma} - (T_m + m\,C_{l\Gamma}) = 0
```

```math
C_{s\Gamma} - k\,C_{l\Gamma} = 0
```

```math
(D_s\partial_n C_s - D_l\partial_n C_l) - V_\Gamma(C_{l\Gamma} - C_{s\Gamma}) = 0
```

Thermal Stefan closure used after each monolithic solve:

```math
V_\Gamma = \frac{\kappa_s\partial_n T_s - \kappa_l\partial_n T_l}{\rho L}
```

## Nonlinear Time-Step Iteration

For each time step:
1. Start from lagged speed `V_\Gamma^{(m)}`.
2. Predict temporary interface geometry from `V_\Gamma^{(m)}`.
3. Assemble and solve one monolithic linear system for `(T_l, T_s, C_l, C_s)`.
4. Update `V_\Gamma^{(m+1)}` from thermal Stefan flux balance.
5. Under-relax speed and iterate until residual tolerances are reached.
6. Advance level-set field with converged speed.
