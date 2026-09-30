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

`issue2_broken_loft` also broke -- volume 2055.8 -- while FreeCAD built
sketch edges on reversed curves (2026-08-23 to 2026-09-30, fixed there,
FreeCAD docs/TransactionLog.md sec 27.87): the ring's radii come from the
ends of an external edge, which came in reversed.

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

The case names are older than a look at the faces: `makeCylinder` puts its
top (z = height) at Face2 and its bottom at Face3, so `*_bottom_*` removes the
top. The names are kept.

Status on `LinkVibe-801` (2026-09-30):

| Case | Status | Symptom |
|------|--------|---------|
| `cyl_*_out/in` (top, bottom), `hole_*` (all eight) | PASS | |
| `cyl_side_out/in` (remove lateral face) | **XFAIL** | `StdFail_NotDone` -- the only removed face is the seam-carrying lateral face. Upstream 8.0.1 throws the same. |

Six of these were XFAIL from 2026-08-01 to 2026-09-30 and are fixed (FreeCAD
docs/TransactionLog.md sec 27.89). All six were casualties of the fix chain,
not limits of OCCT: upstream's `BRepAlgo_Loop`/`BRepOffset` (at 91be8c4c71,
built into scratch libraries and preloaded) passes them with the volumes above.
Three causes, all in the fork's `BRepAlgo_Loop`:

- `cyl_bottom_in`, `hole_bottom_in`: on a periodic face the seam wire replaced
  the closed-edge wires found before it, but the search found the other
  closed edge again afterwards and kept it as a wire of its own -- the inner
  wall came out with wires [1, 3, 1] and the floor went missing. A wire made
  only of a seam wire's edges is now dropped.
- `hole_outer_*`, `hole_inner_*`: a removed cylinder leaves a wall at each
  end, each closed by its own piece of the seam; `FindLoop` allowed one seam
  wire per face, so the second wall was built from two bare circles. Two
  seam wires may now share a face if they share no edge.
- (found by the wider sweep) the concave-face pass in `Perform` rebuilt the
  const edges as FORWARD copies whenever their count changed, so a closed edge
  lost its orientation and the seam wire closed on two circles running the
  same way (`pocket_bottom_in`).

The 2026-08-01 leads (`TrimEdge`'s end-vertex skip, `ContextIntByArc`'s
second `UpdateVertex`, the `Inter2d::Compute` pair skip, `WireInfo`'s
orientation-blind equality) were not the cause of any of these and are
unverified.

## Regression cases from the upstream comparison (2026-09-30)

A sweep of 712 thickness runs -- 16 solids, every face, +-1, intersection
on/off, Arc/Intersection join -- against upstream's libraries found the fork
worse than upstream in 144 runs and better in 51, counting as right only a
valid closed solid whose volume is plausible for a skin (upstream returns the
input unhollowed in places, a valid solid that is wrong). After the fixes of
sec 27.89: worse in 26, better in 53, and in 36 more the fork throws where
upstream returns a wrong solid without a word. Cases added for each cause:

| Case | Model | What it covers |
|------|-------|----------------|
| `ellipse_*` | a single closed ellipse edge, extruded | the offset of an ellipse is a closed B-spline that is not periodic; `ExtentEdge` stretched it 100 lengths past its ends (points at 1e34): the top came back unhollowed, the bottom as two shells |
| `pocket_bottom_*` | cylinder with a blind pocket in its top | the const edge orientation (inward), and the re-found closed edge (outward) |
| `boxhole_hole_*` | box with a through hole, the hole removed | two seam wires on one face |
| `pocket_inter_join_no_crash` | box with a pocket, inward, intersection on, Intersection join | `BuildSplitsOfTrimmedFaces` never sets the trimmed-to-infinite edge map and `UpdateIntersectedEdges` dereferenced it (upstream code, reached only with the fork's edges); still invalid, as upstream's result is |

## Determinism cases (intersection mode, Arc join)

`arc_inter_boss_same_every_run`, `arc_inter_lbox_same_every_run`,
`arc_inter_boxhole_same_every_run` run one thickness (intersection on,
join=Arc, +1) eight times on freshly built shapes and require one result.
`BuildOffsetByArc` iterated its offsets (`MapSF`, a DataMap) in hash order,
and a shape's hash is its TShape's address, so the order -- and with it the
result -- changed from run to run: up to four results in eight runs of the
same input. The offsets are now taken in the shape's topological order
(FreeCAD docs/TransactionLog.md sec 27.88). The cases check that the result is
the same every run, not that it is right. The boss and the box with a blind
hole settled on an invalid solid until sec 27.89's loop fixes; all three are
valid now, the boss and the box with the hole at upstream's volumes.

## Layout

```
tests/thickness/
  README.md          this file
  run_tests.py       the suite (FreeCADCmd script)
  models/            captured user models, one per reported issue
```
