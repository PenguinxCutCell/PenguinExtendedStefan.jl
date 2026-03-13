# Validation

`PenguinExtendedStefan.jl` uses a staged validation ladder.

1. API and construction smoke tests.
2. Monolithic assembly and solve smoke checks.
3. Fixed-interface interface-algebra checks (`TlΓ=TsΓ`, `CsΓ=kClΓ`, `TlΓ=Tm+mClΓ`).
4. Reduction to classical Stefan behavior (`k=1`, `m=0`, decoupled concentration).
5. 1D Martin-Kauffman benchmark (primary binary-alloy reference).
6. 1D manufactured fixed-interface verification.
7. 2D manufactured fixed-interface curved-interface smoke verification.
8. 1D prescribed-motion manufactured verification.
9. 2D prescribed-motion and curved-interface physics smoke checks.
10. Exploratory planar-front perturbation example (not a validated dendrite benchmark).

## Primary Benchmark

Primary 1D binary-alloy benchmark:
- [AFiD-MuRPhFi Stefan examples: 1-D multicomponent melting (Martin and Kauffman, 1977)](https://chowland.github.io/AFiD-MuRPhFi/examples/stefan/)

The package includes:
- reference helper constructors,
- analytical field evaluators,
- interface-position comparison,
- mesh-convergence checks on interface motion.

## Notes On 2D CI Tests

2D cut-cell interface configurations can be sensitive to grid/interface alignment.
Automated CI tests therefore emphasize deterministic bounded-error smoke checks in 2D,
while stricter convergence assertions are enforced in robust 1D validation cases.
