# EECE 306 Electromagnetics Toolbox

**Team:** TODO - add team number and member names.

This MATLAB/GNU Octave package grows into an electrostatic and magnetostatic
field solver over EECE 306. The fixed function contract is in
`EECE306-API-Spec.pdf`; this repository owns the implementations, tests, and
published lab notebooks that demonstrate that contract.

## Run the tests

Start from this directory. Add the toolbox root - the folder containing
`+em` - to the path, never `+em` itself.

```sh
cd EECE306-Toolbox
octave --quiet --eval "runTests"
```

`runTests.m` and `+em/+test/assertClose.m` are provided by the course and must
not be modified. Every `tests/test_labNN.m` file is cumulative: new work must
keep all earlier lab tests passing.

## Lab deliverable pattern

For each lab, commit the implementation, its test file, and its notebook
together.

| Item | Location | Notes |
| --- | --- | --- |
| Toolbox functions | `+em/+module/` | Exact API names, Nx3 shape discipline, SI units |
| Regression tests | `tests/test_labNN.m` | Add checks; do not delete earlier checks |
| Lab notebook | `LabNN_notebook.m` | Run `publish('LabNN_notebook.m')` from this directory |
| Canvas ZIP | `submissions/TeamNN_LabNN.zip` | Local delivery bundle; ignored by Git |

Keep the initial `CONTENTS.m` placeholder until that package gains its first
real function, then replace it with documented `.m` functions. Every function
needs a useful MATLAB/Octave help block because grading includes `help` checks.

## Modules

| Module | Purpose | First lab |
| --- | --- | --- |
| `em.const` | Physical constants | 1 |
| `em.vec` | Vector algebra | 1 |
| `em.coord` | Coordinate and vector-component conversion | 2 |
| `em.src` | Charge and current source construction | 3 |
| `em.field` | Field evaluation and integrals | 3 |
| `em.viz` | Field visualization | 3 |
| `em.quad` | Numerical quadrature | 4 |
| `em.op` | Differential operators | 5 |
| `em.solve` | Boundary-value solvers | 7 |
| `em.derive` | Energy and derived quantities | 8 |
| `em.test` | Verification helpers | 1 onward |

## Lab 1

Implement `em.const.eps0`, `mu0`, `c0`, `eta0` and `em.vec.mag`, `unit`,
`angle`, `fromTo` against `docs/Lab01-3.pdf`. Use `Lab01_notebook.m` to
produce the required report. The provided Lab 1 test currently fails until
those eight functions are in place; that is the expected initial state.


## Lab 8: dielectric interfaces and derived quantities

Run `Lab08_notebook` from this toolbox root. It runs the shared, uncached
experiments in `tests/lab08_studies.m` and the full regression suite. The
completed notebook is the Lab 8 deliverable; no separate formal report is
included. `docs/Lab08.pdf` is the original course handout.

`poisson2d` now discretizes `div(eps grad V)=-rho` as conservative face
fluxes. `opts.faceMean` selects `harmonic` (default) or `arithmetic` for the
interface study. Unfixed exterior faces are insulating. At least one fixed
node anchors the potential. All four solution methods use the same stencil.
The three new derived APIs are `energyE` (J/m), `capacitance` (F/m), and
`resistance` (ohm m). Capacitance uses stored energy and requires harmonic
faces; resistance sums terminal face current. For actual depth d, multiply
capacitance by d and divide resistance by d.

The notebook adjusts the starter's plate geometry to a small gap/width ratio,
places the constant-layer interface between nodes, and explicitly measures
coax refinement. Constant layers then agree with their series formula to
round-off, so fitting that error would be meaningless. The nontrivial order
study uses two *graded* layers with a factor-four jump and an analytic series
integral; it does not pretend to measure second order from round-off.

The microstrip sweep covers w/h=0.3 to 3, one decade. Grounded outer walls
make it a boxed approximation to an open microstrip. Mesh and box changes
are reported separately at the widest strip, not as a universal error bound.
The air-filled solve uses the same box and mask, with quasi-TEM impedance
`Z0=1/(c0*sqrt(C*Cair))`. Staircased conductor and dielectric boundaries
limit geometric accuracy even when the linear solve has a tiny residual.


## Lab 9: current sources and Biot-Savart

`Lab09_notebook.m` runs the shared, uncached experiments in
`tests/lab09_studies.m`, followed by all nine test files through `runTests`.
The new APIs are `em.src.lineCurrent`, `em.src.surfCurrent`, `em.field.H`
(A/m), and `em.field.B` (T). `B` defaults to relative permeability one and
uses `em.const.mu0()`. This gives 47 API functions including the provided
assertion helper.

`em.quad.line` has an optional third output, the directed raw tangent dr/dt.
Its first two outputs are unchanged. `lineCurrent` uses that tangent and
the arc-length weight to form `Idl`; it does not differentiate again.
`surfCurrent` evaluates vector K at the quadrature positions and multiplies
it by area weights. Both constructors store weights inside `Idl` exactly
once. `H` sums `Idl cross (r-source)/R^3`, so it must not multiply by `w`
again. Existing `E` and `merge` are unchanged.

Notebook changes beyond filled explanations remove the obsolete instruction
to fill placeholders, show small loop errors in scientific notation, add the
required aspect sweep, and resolve two starter geometry issues. The finite
wire is plotted using full length L, matching Lab 4; the Lab 9 handout calls
L/2 h. Both threshold conventions are printed. The solenoid uses 200
midpoint turns rather than 25 endpoint turns, and its axial field is checked
against an independent sum of closed-form loop fields and a refined winding.
The reported 99 percent crossings are first sampled brackets, not exact
continuous thresholds. The H timing is the measured review-to-first-passing
kernel-check interval, not an estimate of team labor or total Lab 9 time.

The source sheet test checks a finite square's analytic solid-angle field;
additional checks cover nonuniform K, surface parameter reversal, off-axis
loop refinement, current and curve reversal, B scaling, element-indexed
singularity errors, mixed-source rejection, and 5000/10000 point timings.
