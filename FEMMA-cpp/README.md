# FEMMA — C++ port

A C++/Eigen reimplementation of the `FEM/` core of
[FEMMA](https://github.com/kikioh/FEMMA): finite-element time discretisation of
the Liouville–von Neumann equation for quantum optimal control, with the
plain-MMA and GCMMA optimisers.  Three element orders (linear, quadratic,
Hermite), two equations (Liouville–von Neumann and Bloch), two waveform forms
(phase and xy), and offset/power ensembles are supported.  The GRAPE benchmark
in `lib/Spinach_modified/` of the original is out of scope.

This port carries **two licences**: the core is **MIT**, and the optimiser in
`src/mma/` — translated from Krister Svanberg's GPL MMA/GCMMA code — is
**GPL-3.0-or-later**.  A compiled build links the two, so the **distributed
binary is GPL** while the individual core sources stay MIT.  Each file states
its licence with an `SPDX-License-Identifier` line.  See [`LICENSE`](LICENSE),
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md), and the
[Licensing](#licensing) section below.

## Status

| Layer | MATLAB source | C++ status |
|-------|---------------|------------|
| spin operators | `operatorCreate.m` | **done, validated** (`operator_create`) |
| spin states | `stateCreate.m` | **done, validated** (`state_create`) |
| element matrices (linear, forward) | `eleMatrix_linear.m` | **done, validated** (`ele_matrix_linear`) |
| FEM mesh (linear, LvN) | `meshFEM.m` | **done, validated** (`mesh_fem`) |
| global assembly | `globalStiff.m` | **done, cross-validated vs MATLAB** (`global_stiff`) |
| forward solve (Dirichlet BC) | `femSolve.m` | **done, cross-validated vs MATLAB** (`femsolve`) |
| element-matrix gradient (phase) | `eleMatrix_linear.m` | **done, cross-validated vs MATLAB** (`ele_matrix_grad_linear`) |
| adjoint gradient (rectangle) | `femSolve.m` | **done, finite-difference validated** (`femsolve_grad`) |
| Helmholtz smoothing | `helmSolve.m`, `meshFEM.m` | **done, finite-difference validated** (`helm_solve`, `helm_matrix`) |
| ensemble loop | `ensembleSpin.m` | **done, validated** (`ensemble_spin`, OpenMP) |
| chain-rule combine | `combineGrad.m` | **done** (`combine_grad`) |
| objective assembly (phase, least-square) | `spinSolve_phase.m` | **done, cross-validated vs MATLAB** (`spin_solve_phase`) |
| MMA subproblem solver | `subsolv.m` | **done, validated** (`subsolv`, both branches) |
| MMA iteration and KKT check | `mmasub.m`, `kktcheck.m` | **done, toy-problem validated** |
| optimisation driver (plain MMA) | `femMMA.m`, `mmaInit.m` | **done, end-to-end** (`fem_mma`, `mma_init_phase`, `mma_init_xy`) |
| xy-form objective and gradient | `eleMatrix_linear.m`, `spinSolve_xy.m`, `sigmoid.m`, `powerConstr.m` | **done, finite-difference validated** |
| GCMMA path | `gcmmasub.m`, `asymp.m`, `concheck.m`, `raaupdate.m` | **done, toy-problem validated** |
| quadratic element | `eleMatrix_quad.m` | **done** (`ele_matrix_quad`, eq.35 order three, gradient validated) |
| Hermite element | `eleMatrix_hermite.m` | **done, in optimiser** (`ele_matrix_hermite`, eq.35 order four, nodal gradient validated) |
| Bloch equation | `eleMatrix_bloch.m` | **done, in optimiser** (`bloch_solve`, `bloch_solve_phase`, order two, gradient validated) |

## Interface mapping to the MATLAB

The 2025.09 MATLAB reorganised some signatures while leaving the mathematics
unchanged.  The correspondence is:

- `[K,f,Kgrad] = globalStiff(fem,control)` maps to `global_stiff` (which
  returns `K`) plus `ele_matrix_grad_linear` (the per-element `Kgrad` blocks).
  The load `f` is the all-zero homogeneous load for the LvN equation, so it is
  not threaded through; it is reconstructed inside the solve as the Dirichlet
  right-hand side.  For the Bloch equation `f` is non-zero (the relaxation
  source) and is assembled explicitly inside `bloch_solve`.
- `[x,grad] = femSolve(control,K,f,Kgrad_L)` maps to `femsolve` (forward) and
  `femsolve_grad` (forward plus gradient).
- `control.drifts` is now the starting value of the slice Liouvillian and is
  carried in `control_t::drifts`.

## Cross-validation against MATLAB

`tools/dump_case.cpp` reproduces a 1-spin, 10-slice case and writes the
operators, `K`, `x`, and the gradient blocks to `/tmp` for an external diff.
Against the reference `variables.mat` from the `kikioh/FEMMA` repository the
agreement is to machine precision:

| quantity | max&#124;MATLAB-C++&#124; |
|----------|---------------------------|
| `Lx`, `Ly`, `Lz` (Liouville) | 0 |
| `K` (raw global stiffness, 44x44) | 0 |
| `x` (forward solution, 44x1) | 4.4e-16 |
| `Kgrad` (gradient blocks, 80x80) | 0 |
| fidelity | 14 significant figures |

The final adjoint gradient is degenerate (about 1e-17) for the zero-waveform
reference, so it is validated separately against central finite differences of
the forward solve at a random waveform, agreeing to 3e-11.

The full phase-form pipeline (Helmholtz smoothing, the eleven-offset ensemble,
the adjoint gradient, the chain rule, and the least-square assembly) is
cross-validated against `spinSolve_phase` at a random waveform from
`variables.mat`:

| quantity | size | max&#124;MATLAB-C++&#124; |
|----------|------|---------------------------|
| fidelity | scalar | 1e-16 |
| `fval` (constraints) | 23 | 2e-16 |
| `dfdx` (constraint gradient) | 23x10 | 8e-17 |

The MMA solver is validated on Svanberg's convex toy problem, reproducing the
known optimum `x*=[2.0175,1.7800,1.2375]` with both constraints active and a
KKT residual of 1e-6, and on a second problem with more constraints than
variables to exercise the other branch of `subsolv`.  The GCMMA primitives
(`asymp`, `gcmmasub`, `concheck`, `raaupdate`) are validated by replicating the
GCMMA outer loop on the same toy problem and reaching the same optimum.  End to
end, optimising a single-case phase transfer from a random start drives the
fidelity from about 0 to 0.997 (both with plain MMA and with MMA-GCMMA), and
the xy form, with the sigmoid bound and the power constraint active, reaches
0.99998 in nine iterations.

The xy objective and gradient (`spin_solve_xy`, with `sigmoid` and
`power_constr`) are validated by central finite differences of every
constraint, including the power-constraint rows, agreeing to 3e-10.

The quadratic element (`ele_matrix_quad`, three nodes per slice, dispatched
through `global_stiff` and `femsolve_grad` by the `shape` field) reproduces the
coefficient matrices and assembly of `eleMatrix_quad.m` exactly.  Measured with
the paper's metric (eq.35, the trajectory relative error against step-by-step
exact propagation), the linear element converges at order two and the quadratic
element at order three, with the quadratic about an order of magnitude more
accurate at equal slice count, matching the published convergence figure.  The
quadratic control gradient matches finite differences to 9e-11.  An earlier
final-node comparison against a fine FEM reference suggested order two for both;
that was a metric artefact, the eq.35 trajectory metric recovers the expected
orders.  The shape function is threaded through the ensemble and optimiser by
the `shape` field of `ensemble_t`, so a full multi-case pulse optimisation can
run on quadratic elements; only the Liouville-von Neumann solve changes, the
control parametrisation and the MMA driver are untouched.

The Hermite element (`ele_matrix_hermite`, value and derivative degrees of
freedom per node, a node-based piecewise-linear waveform) reproduces the
coefficient matrices of `eleMatrix_hermite.m`.  Its element stiffness uses two
Liouvillians, at the left and right nodes, and its control gradient assembles a
left and a right block per element into the two nodal control columns; the
`femsolve_grad` Hermite branch does this per element directly, avoiding the
strided index bookkeeping of `femSolve.m`.  Measured with the eq.35 metric it
converges at order four, about a thousand times more accurate than the linear
element at the same slice count, and its nodal gradient matches finite
differences to machine precision.  The Hermite path runs through the optimiser
as well: because the objective is parametrised by the waveform width, supplying
a nodal initial guess, a Helmholtz filter and a mesh sized to the node count,
and an MMA dimension of one per node is enough; a multi-offset phase transfer
optimises to 0.9999 on Hermite elements.

The Bloch equation (`bloch_solve`, `bloch_solve_grad`) is a separate real
pipeline: a three-component magnetisation with relaxation, linear elements, and
a nonzero load vector for the longitudinal relaxation source.  The forward
solve converges at order two against step-by-step exact propagation, and the
control gradient matches finite differences to machine precision.  It runs
through the same MMA optimiser via `bloch_solve_phase` and `bloch_solve_xy`,
which return the same objective record as the Liouville-von Neumann objective.
`bloch_ensemble_grad` evaluates an offset and power grid in parallel, so a
robust z-to-x excitation optimises to 0.999 mean fidelity over an offset
ensemble in both the phase and the xy forms (the latter with the sigmoid bound
and an average-power constraint), with no change to the driver.  The xy Bloch
gradient, including the power-constraint row, is checked independently against
finite differences to 2e-10.  Against the MATLAB EX_bloch run (a z-to-x
excitation, 55 offset and power cases, plain MMA), the port reproduces the
iteration-one fidelity 0.03282 and gradient norm 0.5404, the full 23-iteration
fidelity trajectory to all displayed digits including its non-monotonic dips,
the same stopping iteration, and the final 500-point waveform to 4e-10.

As a full-pipeline cross-validation, the `EX_femma2` benchmark (a single spin,
five power levels, eleven offsets, 55 cases, 500 time slices, phase form,
plain MMA, target fidelity 0.995) was run in C++ from the exact random guess
saved from a MATLAB run.  The per-iteration fidelity reproduces the MATLAB
trajectory at every step including its non-monotonic features, the target is
reached at the same iteration, and the final 500-point smoothed waveform agrees
with the MATLAB result to 1.4e-9 (correlation 1.0). The agreement is bounded
near 1e-9 rather than machine precision because the interior-point subproblem
tolerance is 1e-7 and the linear back ends differ; both effects are
deterministic, so the two optimisers track the same path.

## Build

Needs a C++17 compiler and Eigen 3.3+ (with its `unsupported` modules).

```bash
# Eigen via the system package (recommended on Linux)
sudo apt install libeigen3-dev cmake build-essential

cmake -S . -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Eigen is located in three ways, in order: a system install (`Eigen3::Eigen`),
a copy vendored under `third_party/` (so that `third_party/Eigen/Dense`
exists), or, failing both, a network download via CMake `FetchContent`.

On Windows with Visual Studio 2022, either open the folder as a CMake project
or build from the *x64 Native Tools Command Prompt*:

```bat
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release --target ex_femma2
build\Release\ex_femma2.exe
```

To keep build artefacts out of a synced folder (e.g. OneDrive), point the build
directory elsewhere with `-B C:\some\path\outside\onedrive`.

The ensemble loop uses OpenMP when available; CMake links it automatically and
falls back to serial execution otherwise.

### Build options

| Option | Default | Effect |
|--------|---------|--------|
| `FEMMA_BUILD_TESTS` | `ON` | build and register the `ctest` suite |
| `FEMMA_BUILD_EXAMPLES` | `ON` | build the programmes in `examples/` |
| `FEMMA_BUILD_TOOLS` | `OFF` | build the MATLAB cross-validation harnesses in `tools/` (need reference data in `/tmp`) |

```bash
cmake -S . -B build -DFEMMA_BUILD_TESTS=OFF   # library and examples only
```

### Examples

After building with `FEMMA_BUILD_EXAMPLES=ON` (the default):

```bash
./build/ex_femma2           # single-spin Lz->Lx transfer, the EX_femma2.m case
./build/element_orders      # convergence order of the linear/quad/Hermite elements
./build/bloch_excitation    # Bloch z->x excitation, phase and xy forms
./build/lvn_phase           # broadband Liouville-von Neumann phase pulse, full trajectory
```

`ex_femma2.cpp` is a C++ transcription of `EX_femma2.m`: it builds the same
operators, states, control parameters, mesh, and MMA setup, runs the
optimisation, prints the fidelity trajectory, and writes the optimised pulse as
a Bruker JCAMP-DX shape file (`include/femma/shape_export.hpp`) with a header of
pulse parameters.  The sampling grids (offset, power, point count) live only in
its `PARAMETERS` block, so it is the natural entry point for reading the library
against the MATLAB.

### Using the library

Link against the `femma_fem` target and include `femma/…` headers.  The whole
optimiser is driven through one objective handle:

```cpp
femma::objective_fn obj=[&](const Eigen::MatrixXd& x){
    return femma::spin_solve_phase(fem,helm_K,ens,1.0,x);   // or bloch_solve_phase, …
};
femma::optimise_result r=femma::fem_mma(mma,xval0,obj);
```

Anything that returns a `spin_objective` (Liouville–von Neumann phase or xy,
Bloch phase or xy) plugs into the same `fem_mma` driver unchanged.

## Project layout

```
include/femma/   public headers (mma.hpp is the optimiser interface)
src/mma/         MMA/GCMMA optimiser, translated from Svanberg's GPL code (GPL)
src/core/        original FEMMA translation: elements, solves, ensembles, drivers (MIT)
tests/           ctest suite (19 executables)
examples/        standalone demonstration programmes
tools/           MATLAB cross-validation harnesses (FEMMA_BUILD_TOOLS=ON)
third_party/     optional vendored Eigen (not tracked)
```

The build produces two libraries, `femma_mma` (GPL) and `femma_core` (MIT
sources), with `femma_core` linking `femma_mma`; that single edge is the whole
dependency across the GPL boundary.  The `femma_fem` interface target links both
for convenience.  The split is deliberate: replacing `src/mma/` with a
permissive optimiser behind `include/femma/mma.hpp` frees the rest of the
library (see [Licensing](#licensing)).

## Licensing

This C++ port carries two licences, stated per file with an
`SPDX-License-Identifier` line:

- **`src/mma/`** (`subsolv`, `mmasub`, `kktcheck`, `asymp`, `gcmmasub`,
  `concheck`, `raaupdate`) is a translation of Krister Svanberg's MMA/GCMMA
  reference code, which is GPL for academic, non-commercial use, and is
  therefore **GPL-3.0-or-later**.
- **Everything else** — `src/core/`, the headers, examples, and tests — is an
  original translation of the FEMMA MATLAB code and is **MIT**.

The two parts build as separate libraries (`femma_mma` and `femma_core`).  When
they are compiled and linked into one library or executable, the result is a
combined work that includes GPL code, so **the distributed binary is governed
as a whole by GPL-3.0-or-later** — even though most of the source files are MIT.

Because the MIT core does not derive from Svanberg's work, a fully permissive
build is possible: replace `src/mma/` with an independent CCSA/MMA optimiser
(for example NLopt's CCSA) behind the `include/femma/mma.hpp` interface; nothing
in `src/core/` depends on the Svanberg sources beyond that header.  Place the
full GPL-3.0 text in a `COPYING` file before redistributing.  See
[`LICENSE`](LICENSE) for the MIT text and the GPL notice, and
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) for citations.

This port is published as the `femma-cpp/` folder of the FEMMA repository; the
MATLAB implementation at the repository root remains MIT.

## Design decisions that differ from the MATLAB

These follow the Spinach contribution guidelines, adapted to C++:

- **No default arguments.** `operatorCreate`/`stateCreate` used `~exist` to
  default `coherent` and `formalism`; here those are required `enum class`
  arguments (`build_mode`, `formalism`).
- **Fixed output shape.** `stateCreate.m` returned a matrix in Hilbert space
  or a vector in Liouville space. `state_create` always returns the
  Liouville-space column vector, since that is all the pipeline consumes;
  shape-switching is exactly what the guidelines forbid.
- **British spelling, comments above code, no spaces around operators.**
- **Column-major everywhere** (`cmat`, `cvec`, `csparse`) so that storage and
  vectorisation order match MATLAB.
- **No stored assembly index.** `globalStiff.m` scatters element blocks with a
  precomputed `KIndex`; for linear elements the scatter is the closed form
  `dof_node*element+local`, applied directly in `global_stiff`, so no index
  array is stored.
- **Uniform-grid guard.** `eleMatrix_linear.m` silently takes `dt` from the
  first slice; `ele_matrix_linear` keeps that but throws if `pulse_dt` is not
  uniform, turning a silent wrong result into an explicit error.

## Conventions

All public types live in `namespace femma` (`include/femma/spin_types.hpp`).
Spin operators are spin-1/2 with operators equal to the Pauli matrices over
two; the Liouville superoperator is the commutation superoperator
`kron(I,A)-kron(A^T,I)` using the **plain** transpose, not the conjugate
transpose.
