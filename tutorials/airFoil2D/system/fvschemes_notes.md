# Reading `fvSchemes`: The Numerics Behind the Keywords (USING CLAUDE TO TRY AND UNDERSTAND THIS)

Notes on the `airfoil2D` case (steady, incompressible, Spalart–Allmaras), OpenFOAM Foundation v14.

> **Notation.** OpenFOAM calls the face flux field `phi`. To avoid confusion with a generic transported scalar, these notes write the face flux as $F_f$ and the transported quantity as $\phi$.

---

## 1. The basic problem

OpenFOAM stores **one value per cell**, at the cell centre $\mathbf{x}_P$. The governing equations are full of derivatives ($\nabla\cdot$, $\nabla$, $\nabla^2$), and you can't differentiate one number per cell directly.

The finite volume method gets around this by **integrating each term over the cell volume** $V_P$.

---

## 2. "Gauss" = the divergence theorem

Gauss's theorem turns a volume integral of a divergence into a surface integral, which on a polyhedral cell becomes a sum over faces:

$$
\int_{V_P} \nabla\cdot\mathbf{F}\,dV \;=\; \oint_{\partial V_P} \mathbf{F}\cdot d\mathbf{S} \;\approx\; \sum_{f} \mathbf{F}_f\cdot\mathbf{S}_f
$$

where $\mathbf{S}_f$ is the outward face area vector (normal to the face, magnitude equal to its area).

The same idea gives the cell-centred gradient:

$$
(\nabla\phi)_P \;\approx\; \frac{1}{V_P}\sum_f \mathbf{S}_f\,\phi_f
$$

**"Gauss" in the file means: compute this term by summing face contributions.** The only real alternative is `leastSquares`, and it is only offered for gradients.

---

## 3. The second word: how to get the face value $\phi_f$

Gauss's theorem needs $\phi_f$ at each face, but only cell-centre values exist. The word after `Gauss` is the **interpolation scheme**, the rule for estimating $\phi_f$.

### 1D picture (flow left → right, so $F_f > 0$)

```
   W         P         N
   o----|----o----|----o      flow →
                  f
  φ=1       φ=2       φ=4
```

### `linear` (central differencing)

Distance-weighted average of the two neighbouring cells:

$$
\phi_f = w\,\phi_P + (1-w)\,\phi_N, \qquad w = \frac{|\mathbf{x}_N - \mathbf{x}_f|}{|\mathbf{x}_N - \mathbf{x}_P|}
$$

Example: $\phi_f = 3$.

- Second-order accurate.
- **Ignores flow direction.** In convection this can create unphysical oscillations (overshoots and undershoots) near steep gradients.

### `upwind`

Take the value from the cell the flow comes from:

$$
\phi_f =
\begin{cases}
\phi_P, & F_f \ge 0 \\
\phi_N, & F_f < 0
\end{cases}
$$

Example: $\phi_f = 2$.

- Bounded: never creates new extrema. Very robust.
- First order only: it smears gradients, which acts like extra artificial viscosity (**numerical diffusion**). A Taylor expansion shows the leading error looks like $\tfrac{1}{2}|U|\Delta x\,\partial^2\phi/\partial x^2$, which is a diffusion term.

### `linearUpwind` (second-order upwind)

Start from the upwind value and extrapolate to the face using the upwind cell's gradient:

$$
\phi_f = \phi_U + (\nabla\phi)_U\cdot(\mathbf{x}_f - \mathbf{x}_U)
$$

Example: with a central-difference gradient at $P$,

$$
(\nabla\phi)_P \approx \frac{\phi_N - \phi_W}{2\Delta x} = \frac{3}{2\Delta x},
\qquad
\phi_f = 2 + \frac{3}{2\Delta x}\cdot\frac{\Delta x}{2} = 2.75
$$

- Second order **and** respects flow direction.
- Not strictly bounded, so it can overshoot slightly.
- **This is why the keyword says `grad(U)`:** the scheme needs a gradient to extrapolate with, and that gradient is computed using whatever `gradSchemes` specifies.

### Limited schemes (`limitedLinear`, `vanLeer`, `MUSCL`, …)

These blend between a high-order scheme and `upwind` based on how steep the local solution is. They behave like the high-order scheme in smooth regions and fall back toward upwind near sharp gradients to prevent overshoots.

### The core tradeoff

**Godunov's theorem:** no *linear* scheme of order higher than one can guarantee that it creates no new extrema. Every convection scheme is a compromise between accuracy and boundedness. Limiters get around the theorem by being nonlinear, adapting to the solution.

---

## 4. Why convection gets fancy schemes and diffusion doesn't

| Term | Physics | Natural scheme |
|---|---|---|
| Convection | Carries information **downstream**, so direction matters | Upwind-biased (`upwind`, `linearUpwind`, limited) |
| Diffusion | Spreads information **in all directions equally** | `linear` (central) |

That's why `divSchemes` holds the interesting choices and `laplacianSchemes` is plain `linear`.

---

## 5. Convection: `div(phi,U)` and `div(phi,nuTilda)`

The integrated convection term is

$$
\int_{V_P}\nabla\cdot(\mathbf{U}\phi)\,dV \;\approx\; \sum_f (\mathbf{S}_f\cdot\mathbf{U}_f)\,\phi_f \;=\; \sum_f F_f\,\phi_f
$$

The face flux $F_f$ (`phi`) comes from the pressure–velocity coupling step, not from simply interpolating $\mathbf{U}$. The scheme only determines $\phi_f$.

### Deferred correction

`linearUpwind` is not built fully into the matrix. The upwind part is **implicit** and the high-order correction is **explicit**, using values from the previous iteration $n$:

$$
\sum_f F_f\,\phi_f \;\approx\; \underbrace{\sum_f F_f\,\phi_f^{\text{UD}}}_{\text{implicit (iteration } n+1)} \;+\; \underbrace{\sum_f F_f\left(\phi_f^{\text{LUD}} - \phi_f^{\text{UD}}\right)}_{\text{explicit (iteration } n)}
$$

This keeps the matrix diagonally dominant (stable linear solves). At convergence the correction is exact, so you still get second-order accuracy.

### `bounded`

`bounded Gauss <scheme>` actually discretizes

$$
\nabla\cdot(F\,\phi) \;-\; (\nabla\cdot F)\,\phi
$$

Analytically, $\nabla\cdot\mathbf{U} = 0$ for incompressible flow, so the second term is zero. But in a steady SIMPLE run, continuity isn't satisfied until convergence, and $\nabla\cdot F \neq 0$ acts as a spurious source or sink. Subtracting it keeps the iterations stable and the solution bounded. **It is a steady-state device.** It is not appropriate for time-accurate transient runs.

---

## 6. Diffusion: `laplacianSchemes  Gauss linear corrected`

$$
\int_{V_P}\nabla\cdot(\Gamma\nabla\phi)\,dV \;\approx\; \sum_f \Gamma_f\;\mathbf{S}_f\cdot(\nabla\phi)_f
$$

- `Gauss`: sum over faces.
- `linear`: interpolate the diffusivity $\Gamma$ (e.g. $\nu_{\text{eff}}$) to the face.
- `corrected`: the **surface-normal gradient** scheme (see §7).

---

## 7. Surface-normal gradients: `snGradSchemes  corrected`

On a perfectly **orthogonal** mesh, the vector $\mathbf{d} = \mathbf{x}_N - \mathbf{x}_P$ is parallel to $\mathbf{S}_f$, and

$$
\mathbf{S}_f\cdot(\nabla\phi)_f \;=\; |\mathbf{S}_f|\,\frac{\phi_N - \phi_P}{|\mathbf{d}|}
$$

which is compact, second order, and fully implicit.

On a **non-orthogonal** mesh, split $\mathbf{S}_f$ into a part along $\mathbf{d}$ and a remainder:

$$
\mathbf{S}_f = \boldsymbol{\Delta} + \mathbf{k}, \qquad \boldsymbol{\Delta}\parallel\mathbf{d}
$$

$$
\mathbf{S}_f\cdot(\nabla\phi)_f \;\approx\; \underbrace{|\boldsymbol{\Delta}|\,\frac{\phi_N - \phi_P}{|\mathbf{d}|}}_{\text{implicit}} \;+\; \underbrace{\mathbf{k}\cdot(\nabla\phi)_f}_{\text{explicit correction}}
$$

- `corrected` includes the explicit correction.
- `uncorrected` drops it (fine only on near-orthogonal meshes).
- The explicit part is why `fvSolution` has `nNonOrthogonalCorrectors`: extra passes let the correction converge within each iteration.

---

## 8. The viscous stress term: `div((nuEff*dev2(T(grad(U)))))`

The full deviatoric viscous stress divergence is

$$
\nabla\cdot\Big[\nu_{\text{eff}}\big(\nabla\mathbf{U} + (\nabla\mathbf{U})^T\big) - \tfrac{2}{3}\nu_{\text{eff}}(\nabla\cdot\mathbf{U})\,\mathbf{I}\Big]
$$

OpenFOAM splits it into

$$
\underbrace{\nabla\cdot\big(\nu_{\text{eff}}\nabla\mathbf{U}\big)}_{\text{implicit Laplacian}} \;+\; \underbrace{\nabla\cdot\Big(\nu_{\text{eff}}\,\mathrm{dev2}\big((\nabla\mathbf{U})^T\big)\Big)}_{\text{explicit, this keyword}}
$$

with

$$
\mathrm{dev2}(\mathbf{A}) = \mathbf{A} - \tfrac{2}{3}\,\mathrm{tr}(\mathbf{A})\,\mathbf{I}
$$

Notes:

- For constant $\nu$: $\nabla\cdot(\nabla\mathbf{U})^T = \nabla(\nabla\cdot\mathbf{U}) = 0$. The term only matters because $\nu_{\text{eff}} = \nu + \nu_t$ **varies in space**.
- The keyword is the **literal expression string from the C++ source**, so it must match character for character. To find these names, set `default none` and run: the solver errors out naming the missing keyword.

---

## 9. Wall distance: `wallDist  method meshWave`

Spalart–Allmaras has a destruction term

$$
-\,c_{w1}\,f_w\left(\frac{\tilde{\nu}}{d}\right)^2
$$

so the turbulence model depends directly on the wall distance $d$.

- `meshWave`: fast, approximate. It propagates distance across neighbouring faces and can lose accuracy on strongly curved or skewed near-wall meshes.
- `Poisson`: solves a PDE for $d$. Smoother, also approximate.

---

## 10. Decoding the file

| Entry | Meaning |
|---|---|
| `ddtSchemes  steadyState` | Drop the time derivative. Stability comes from under-relaxation in `fvSolution`, not from time stepping |
| `gradSchemes  Gauss linear` | Green–Gauss gradient with linearly interpolated face values |
| `divSchemes  default none` | Force every convection term to be specified explicitly |
| `div(phi,U)  bounded Gauss linearUpwind grad(U)` | Momentum convection: face sum, second-order upwind using ∇U, with the steady-state continuity fix |
| `div(phi,nuTilda)  bounded Gauss linearUpwind grad(nuTilda)` | Same for the SA variable |
| `div((nuEff*dev2(T(grad(U)))))  Gauss linear` | Explicit transpose part of the viscous stress |
| `laplacianSchemes  Gauss linear corrected` | Diffusion: face sum, linear diffusivity, non-orthogonal-corrected face-normal gradient |
| `interpolationSchemes  linear` | Default face interpolation for anything not otherwise specified |
| `snGradSchemes  corrected` | Face-normal gradient with explicit non-orthogonal correction |
| `wallDist  meshWave` | Approximate wall distance for the SA destruction term |

**General grammar:** `[bounded] <integration> <face interpolation> [<extra input the scheme needs>]`

---

## 11. Choices worth questioning

1. **Unlimited `linearUpwind` on $\tilde{\nu}$.** $\tilde{\nu}$ should stay non-negative, but `linearUpwind` can overshoot or undershoot near the boundary-layer edge or in the wake. A more cautious choice is `limitedLinear 1` or `upwind` for turbulence variables. The tutorial picks accuracy over robustness.
2. **"Second order" is only nominal.** Convection is formally second order, but the gradients feeding it are Green–Gauss on high-aspect-ratio boundary-layer cells, where it loses accuracy. `leastSquares` or `cellLimited Gauss linear 1` are the alternatives. Formal order is not the same as actual accuracy on a given mesh.
3. **`fvSchemes` is half the story.** With `steadyState`, the relaxation factors, non-orthogonal correctors and residual controls in `fvSolution` decide whether these schemes converge. Read the two files together.

---

## 12. Questions to work through

1. Why would `bounded` be wrong in a transient run?
2. If `div(phi,U)` were changed to `bounded Gauss linear`, what would you expect in the solution, and where on the airfoil would it show up first?
3. **Experiment:** run `div(phi,U)` with `linearUpwind`, `upwind`, and `limitedLinearV 1`. Compare $C_L$, $C_D$, and residual histories. Predict the ordering before you run.
4. **Source reading:** `src/finiteVolume/interpolation/surfaceInterpolation/schemes/linearUpwind/`. Find where the explicit correction from §5 is built.