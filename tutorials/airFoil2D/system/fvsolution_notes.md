# fvSolution — how SIMPLE actually works

Notes from working through `system/fvSolution` for the airfoil2D case
(steady, incompressible, Spalart–Allmaras, SIMPLE).

---

## 1. The problem

Unknowns per cell: **U**, **p**, **ν̃** (`nuTilda`).

| Unknown | Equation |
|---|---|
| U | momentum |
| ν̃ | Spalart–Allmaras transport |
| p | **none** |

Continuity (∇·U = 0) has no p in it. It's a constraint on U, not an equation for p.
Pressure is whatever it has to be so that the velocity it produces conserves mass.

**How p affects continuity:** p → ∇p (force in momentum) → U → ∇·U.
Through velocity, *not* density.

> My first guess was p → ρ via the equation of state. That's the **compressible**
> picture. Here ρ is constant, so there's no p–ρ link. Incompressible solvers use
> kinematic pressure p/ρ (units m²/s²). Losing the p–ρ link is exactly why
> incompressible pressure is hard.

The physics isn't iterative — there's one p field (up to a constant) that works.
The *algorithm* is iterative. And it doesn't guess; each pass derives a correction.

---

## 2. Building an equation for p

Discretized momentum for one cell:

    A·U_P = H − ∇p

- A = the cell's own (diagonal) coefficient
- H = neighbour contributions + sources

Solve for velocity:

    U_P = H/A − (1/A)∇p

Plug into continuity:

    ∇·( (1/A)∇p ) = ∇·( H/A )

→ a Poisson equation for p. Momentum gave U in terms of p; continuity then fixes p.

**Catch:** H is built from the current (wrong) velocities, so the p from this is
only approximate → you have to loop.

---

## 3. The SIMPLE loop

1. Momentum predictor → U* (obeys momentum, not continuity)
2. Compute H/A from U*
3. Solve p equation
4. Correct U and face fluxes with new ∇p
5. Solve ν̃, update ν_t
6. Repeat until nothing changes

Steady case → iterations are pseudo-time. Only the converged state means anything.

- `solvers {}` = **inner** loop (how steps 1, 3, 5 are solved)
- `SIMPLE {}` = **outer** loop

**Why it overshoots:** the p equation holds H fixed — it assumes neighbours won't
react to the new pressure. So pressure is made to fix continuity alone, and the
correction comes out **too large**.

---

## 4. Relaxation

Full correction every iteration → back-and-forth oscillation. If it grows → divergence.
Fix: only take part of the correction.

- **Field relaxation** (`fields { p 0.3; }`) — after the solve:
  `p = p_old + 0.3·(p_new − p_old)`. Blunt cap.
- **Equation relaxation** (`equations { U 0.7; }`) — before the solve:
  diagonal A → A/0.7 plus a matching source. Each cell resists change
  (acts like a pseudo-time step). Bonus: more diagonally dominant → easier for Gauss–Seidel.

Rule of thumb: **α_p + α_U ≈ 1**. Too low = slow; too high = oscillates/diverges.

**Does relaxation change the converged answer? No.**
At convergence p_new = p_old, so the bracket is 0 and α multiplies nothing.
Equation relaxation terms cancel when U = U_old.
Relaxation changes the **path**, not the destination.

Caveats:
- Only true if it actually converges.
- Problems with multiple steady solutions (airfoil near stall) can land on
  different solutions via different paths.

---

## 5. Residuals and tolerances

Residual for a linear system Ax = b: **r = b − Ax** (summed, normalized → dimensionless).

Log line:

    GAMG: Solving for p, Initial residual = 3.2e-03, Final residual = 2.9e-04, No Iterations 4

- **Initial** = imbalance before the inner solve (old field in new equation)
- **Final** = imbalance after the inner solve

Inner solver stops at whichever comes first:
- `tolerance` — absolute: final residual < value
- `relTol` — relative: final/initial < value

Early on, `relTol` decides. `tolerance` only matters near the end.
Why `relTol 0.1`: the matrix gets rebuilt next outer iteration, so don't solve it precisely.

`residualControl` checks the **initial** residual — "does last iteration's field
already satisfy this iteration's equations?" Final residual is always small by
construction, so it says nothing about convergence.

**Small residuals ≠ correct answer.** They mean the discrete equations are satisfied
on this mesh with this model. Still need: mesh study, model choice, Cl/Cd levelled off.

### Check: p `tolerance 1e-4` but `residualControl p 1e-5`?

> My answer: "stops early." Wrong — mixed up the inner solver with the whole run.

1. Once p initial residual < 1e-4, it's already below `tolerance`
2. Inner solver does **0 iterations** — p stops being solved
3. p residual stalls around 1e-4
4. Never reaches 1e-5 → `residualControl` never satisfied
5. Run goes to `endTime`, looks unconverged

**Rule:** inner `tolerance` must be tighter than `residualControl`.
"No Iterations 0" in a log = solver thinks the field is already good enough.

---

## 6. Solver choices

- **p → GAMG.** Elliptic: every cell depends on every other. Gauss–Seidel moves
  info ~1 cell per sweep, too slow. Multigrid solves on coarser grids too, so
  domain-wide error gets removed cheaply.
- **GAMG smoother = GaussSeidel.** GAMG uses a smoother on every level: GS removes
  jagged cell-to-cell error; coarse levels handle the smooth domain-wide error.
- **U, ν̃ → smoothSolver + GaussSeidel.** Transport equations are local, and
  relaxation makes them diagonally dominant → GS alone works, 2 sweeps is enough.
- **ν̃ is SA's modified turbulent viscosity** (m²/s), *not* turbulent kinetic
  energy (that's k; SA has no k). ν_t = ν̃·f_v1.

---

## 7. Other

- `nNonOrthogonalCorrectors` — re-solves p to converge the explicit non-orthogonal
  correction. 0 is fine on a good mesh; raise to 1–2 if `checkMesh` max
  non-orthogonality > ~70°.
- No `pFinal`/`UFinal` — those are for PIMPLE (transient) only.
- SIMPLEC (`consistent yes;`) fixes the neighbour approximation properly → can run
  p near 1.0, U ~0.9. Check v14 still accepts the keyword.

---

## Next

- Run with p relaxation 0.3 vs 0.7 → watch the residual history (see the overshoot)
- `relTol 0` on p → compare iteration count vs wall time
- SIMPLEC with p 1.0 / U 0.9 → iterations to Cl levelling off
- Add a `forceCoeffs` function object to `controlDict` to track Cl/Cd