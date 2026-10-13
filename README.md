# Induced magnetic field in a hollow sphere

Axisymmetric finite-element solution for a hollow magnetic sphere (shell with relative permeability μr) placed in a uniform applied field **H**s = hₐ **e**z, using the reduced scalar potential formulation.

**Reference:** I. D. Mayergoyz, M. V. K. Chari, J. D'Angelo, *A New Scalar Potential Formulation for Three-Dimensional Magnetostatic Problems*, IEEE Trans. Magn. **23**(6), 3889–3894, 1987.

## Files

| File | Description |
|---|---|
| `hollowsphere.m` | MATLAB implementation (PDE Toolbox for geometry/mesh, hand-written P1 assembly of the r-weighted stiffness and interface load, φ = 0 on the truncation boundary, element-wise and nodal-averaged H). |
| `hollow_sphere_fenicsx.ipynb` | FEniCSx (dolfinx 0.8–0.10) translation of the MATLAB code, plus the paper's two-step "difference potential" method and an exact solution of the truncated problem for error measurement. |

## Problem setup

| Region | Definition | μ |
|---|---|---|
| Cavity | R < a | μ₀ |
| Shell | a < R < b | μr μ₀ |
| Outer air | b < R < c | μ₀ |

Default parameters: a = 0.1, b = 0.2, c = 0.5 (φ = 0 at R = c), μ₀ = 1, μr = 10, hₐ = 100.

## Running

**MATLAB:** requires the Partial Differential Equation Toolbox. Run `hollowsphere.m`.

**FEniCSx:** requires `dolfinx`, `gmsh` (Python API), `matplotlib`. Run the notebook in serial. Outputs: `hollow_sphere_H.png`, `hollow_sphere_axis.png`, `hollow_sphere.xdmf` (for ParaView).

The notebook checks the FE solution against the exact solution of the truncated problem, compares the surface and volume forms of the right-hand side, and compares the reduced-potential and two-step results. It also reports the cavity field against the infinite-domain shielding formula, which shows the error from truncating the domain at R = c.
