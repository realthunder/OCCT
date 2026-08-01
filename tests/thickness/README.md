# MakeThickSolid / offset regression tests

Regression suite for the fork's `BRepOffset_MakeOffset` ("MakeThickSolid") fix
chain — the code that makes thickness/offset work on shapes with **concave
removed faces**, which upstream OCCT does not support. The chain rewrites
loop building (`BRepAlgo_Loop`), offset-edge intersection
(`BRepOffset_Inter2d/3d`) and face rebuilding (`BRepAlgo_FaceRestrictor`),
and has a history of fixing one model while breaking another; this suite
pins down both the fixed models and the known remaining breakage.

## Running

Any FreeCAD build linked against this OCCT tree:

```
FreeCADCmd tests/thickness/run_tests.py
```

Exit code 0 = no regression (every expected-pass case passed). Known-broken
cases are marked `XFAIL` in the script; when one starts passing the suite
prints `UNEXPECTED-PASS` — promote it to `pass` with its reference volume and
update this file.

Reference volumes were captured on `LinkVibe-801` after the fixes listed
below, cross-checked against the 7.7.2 `LinkVibe` branch (identical to
4 decimals for the document cases), and are stable across runs.

## Document cases (captured user models)

Each `models/*.FCStd` comes from a reported issue. The test forces a full
recompute and requires every feature to produce a valid shape with the
reference volume.

| Case | Model | Issue | What it covers |
|------|-------|-------|----------------|
| `issue1_ellipse_thickness` | PartDesign Pad of an ellipse profile + Thickness | realthunder/OCCT#1 | MakeThickSolid on an elliptic-profile pad: exercises the concave-face loop rebuild (`BRepAlgo_Loop` cut/const edge classification) |
| `issue2_broken_loft` | PartDesign AdditiveLoft (ring sketch → rectangle) | realthunder/OCCT#2 | Thrusection regression guard: the fork reverts upstream 4607bd0747f (`myPercent` 0.01); the loft breaks if that revert is lost |
| `issue3_pad_thickness` | PartDesign Pad + Thickness | realthunder/OCCT#3 | MakeThickSolid producing hollow walls from a padded profile (volume reference detects silently-wrong results, e.g. non-hollowed solid) |
| `issue4_revolution_thickness` | Part Revolution (quarter revolve, two sphere caps) + Thickness, planar faces removed | realthunder/OCCT#4 | Offset of periodic faces: sphere caps with **degenerated pole edges** + seam handling. Historically: 7.7.2 produced an invalid solid; the 8.0.1 port initially threw `StdFail_NotDone` (hash-order nondeterminism); both fixed on `LinkVibe-801` — result must now be a fully valid single-shell solid |

Fix history for these (branch `LinkVibe-801`):
- `1789444318` — port of the 7.7.2 MakeThickSolid fix chain to 8.0.1
- `756262b692` — deterministic wire nesting in `BRepAlgo_FaceRestrictor`
- `933803af41` — hash-order determinism (`myCutEdges`/wire-set iteration) + latent-bug fixes (`AsDes::Replace(x,x)`, vertex chain compression, `Bubble` dedup, `GetCenter`, unguarded `MakeWire`)
- `72aae3644f` — degenerated pole edges re-embedded into the wire crossing the singularity (turns issue4 from invalid to fully valid)

## Programmatic cases (plain solids — the "don't break normal thickness" guard)

The fix chain used to break thickness on ordinary solids. These cases build
primitives in-code (no model files) and hollow them removing **one face at a
time**, both outward (`+1`) and inward (`-1`), with mode=Skin, join=Arc:

- `cyl_*` — `Part.makeCylinder(4, 20)`: Face1 = lateral (has a seam edge),
  Face2 = bottom, Face3 = top.
- `hole_*` — `makeCylinder(5, 10).cut(makeCylinder(2, 10))`: Face1 = outer
  lateral, Face2 = bottom, Face3 = top, Face4 = hole lateral (reversed
  orientation, seam edge).

Status on `LinkVibe-801` (2026-08-01). The same pattern reproduces on the
7.7.2 branch, i.e. these are long-standing fix-chain casualties, not 8.0
regressions:

| Case | Status | Symptom |
|------|--------|---------|
| `cyl_bottom_out`, `cyl_top_out` | PASS | |
| `cyl_top_in` | PASS | |
| `cyl_side_out/in` (remove lateral face) | **XFAIL** | `StdFail_NotDone` — offsetting when the only removed face is the seam-carrying lateral face fails outright |
| `cyl_bottom_in` | **XFAIL** | invalid + open shell — asymmetric with `cyl_top_in` (PASS): the REVERSED bottom plane trips an orientation-dependent path |
| `hole_bottom_out`, `hole_top_out`, `hole_top_in` | PASS | |
| `hole_bottom_in` | **XFAIL** | invalid + open shell (same asymmetry as `cyl_bottom_in`) |
| `hole_outer_out/in` | **XFAIL** | invalid result — outer lateral + hole lateral interaction |
| `hole_inner_out/in` | **XFAIL** | invalid result (7.7.2 additionally threw on `hole_inner_out`) |

Leads for the XFAIL group, from the 2026-08-01 code audit (see the ranked
SUSPECT findings recorded with the fix commits): `TrimEdge`'s end-vertex skip
can leave an extended edge untrimmed; `BRepOffset_Inter3d::ContextIntByArc`'s
extension block mishandles **closed offset edges** (second `UpdateVertex`
clobbers the first trim parameter — the cylinder lateral faces are exactly
this case); the `ContextFaces` pair-skip in `Inter2d::Compute` lacks a
same-face check; `WireInfo` equality is orientation-blind.

## Layout

```
tests/thickness/
  README.md          this file
  run_tests.py       the suite (FreeCADCmd script)
  models/            captured user models, one per reported issue
```
