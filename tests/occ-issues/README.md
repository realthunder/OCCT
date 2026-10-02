# OCC-labeled issues from realthunder/FreeCAD — model corpus + scan

All issues carrying the `occ` label on github.com/realthunder/FreeCAD
(26 as of 2026-08-01), with every recoverable attachment stored under
`models/`. This is raw material for growing the regression suite
(see `../thickness/`); none of these are fixed yet.

**Scan method:** each model opened in its own `FreeCADCmd` process
(crash-isolated), all objects touched, full recompute, then per-object
error states and `isValid()` checked. Scan stack: `LinkVibe-801` OCCT
(debug) + fcad `build/conda-debug-occt801`. "recomputes clean" means the
stored failure needs interactive steps (an edit, an export, a command) —
those need a scripted repro before they can join the suite.

## Reproduces on recompute (ready to become suite cases)

| # | Model(s) | Problem | Scan result |
|---|----------|---------|-------------|
| 64 | `issue064_boolean_split_seam.FCStd` | Part Split through a cylinder **seam edge** produces a wrong boolean; RT: split incorrect near seam, workaround = split at a different angle | invalid: Body, Split, Split_i0, Boolean, Reference |
| 172 | `issue172_sketch_unclosed.FCStd` | Edges imported into sketch as external geometry become unconnected: curve `First/LastParameter` evaluation is ~0.003 mm off the end vertices | errors: Pad006, Pad008; invalid: Sketch008 |
| 273 | `issue273_heater_handle*.FCStd` | Fillet whose surface touches a **seam edge** destroys the shape; chamfer on the same edge segfaults (upstream bug 28354; `ShapeUpgrade_ShapeDivideClosed` suggested as seam removal) | invalid: Pad (both variants; crash itself needs the interactive chamfer) |
| 280 | `issue280_pocket_sphere_seam.FCStd` | Pocket through a mirrored sphere fails ("resulting shape is not a solid") once wide enough — cutting with **overlapping seam edges**; tolerance workaround exists | error: Sketch002 |
| 333 | `issue333_draft_pad_crooked.FCStd` | Draft feature produces a broken model; Pad on top comes out crooked past depth 44 | **crash on recompute** |
| 334 | `issue334_draft_artifact.FCStd` | Same model family: artifact after Draft003; RT: side face broken from Draft001 on | **crash on recompute** |
| 346 | `issue346_offset_one_side.FCStd` | Part Offset works from one side, throws "Unknown OCC exception" from the other | error: Offset |
| 360 | `issue360_fillet_spike.FCStd` | Fillet produces a spike in the model (7 MB model) | **crash on recompute** (closed as OCC bug, still crashes) |
| 363 | `issue363_thickness_multi.FCStd` | "Unexpected result with thickness and other issues" — richest thickness repro in the corpus | errors: Chamfer, Thickness, Thickness001, Thickness002; invalid: Body, Thickness001 |
| 474 | `issue474_fillet_edit_crash.FCStd` | Editing Fillet003 crashes; root cause per RT: sketch misalignment OCC fillet can't absorb; release crash needs `BUILD_RELEASE_DISABLE_EXCEPTIONS=Off` builds | **crash on recompute** |
| 521 | `issue521_thickness_intersection.FCStd` | Thickness **Intersection option** not working | error: Thickness |
| 523 | `issue523_fillet_explodes.FCStd` | Face "explodes" when applying a fillet radius | invalid: Fillet |
| 580 | `issue580_mirror_crash.FCStd` | Mirror produces an invalid shape; the following fuse crashes (whole-shape mirror via MultiTransform) | error: Chamfer; invalid: Body |
| 613 | `issue613_thickness_two_faces.FCStd` | Thickness with two faces selected opens only one of them | error: Fillet; invalid: Body, Thickness |
| 617 | `issue617_intersection_pad.FCStd` | Pad from faces built on a revolution/sketch intersection: extremely slow, then broken geometry (`_simple` variant currently recomputes clean) | error: Pad; invalid: Body, Pad |
| 876 | `issue876_fillet_pocket_flip.FCStd` | Fillet produces invalid geometry; downstream Pocket **adds** material instead of cutting | invalid: Fillet001 |
| 962 | `issue962_fillet_artifact.FCStd` | Strange artifact | invalid: Body001, Fillet, Fillet003, Fillet004 |

## Recomputes clean — needs a scripted repro (the bug is in a step, not the stored state)

| # | Model(s) | Problem / needed script |
|---|----------|------------------------|
| 309 | `issue309_fillet_crash_shape.brp` (+ `notes/issue309_*`) | Fillet on an arc edge crashes in `ChFi3d_Builder::PerformIntersectionAtEnd` (upstream #32929). Script: load brp, `BRepFilletAPI_MakeFillet` on the reported edge. Upstream workaround patch + crash stack kept in `notes/` |
| 310 | `issue310_step_color_export.FCStd` + 2 `.step` refs | STEP export loses the **top face color** of an untransformed cylinder — top/bottom faces share a TShape (partner shapes) and the color goes to the wrong one. Script: export STEP, parse colors |
| 423 | `issue423_fill_crash.FCStd`, `issue423_distance_wrong.FCStd` | `BRepOffsetAPI_MakeFilling` crash invoking Part Fill (known OCC bug, see comment at fcad `TopoShapeEx.cpp` filling site); second file: wrong distance measurement |
| 631 | `issue631_fillet_loft_joint*.FCStd` | Fillet at a loft edge joint gives a wrong **but valid** shape (radius > 3 mm); needs a geometric assertion (volume/area), not validity |
| 652 | `issue652_pad_top_uncovered.FCStd` | Pad on a face whose sketch arc coincides with an existing arc: top surface not covered — tolerance-coincidence case; may need the exact re-pad step |
| 662 | `issue662_thickness_pads_holes.FCStd` | Thickness produced an invalid shape making later pads become holes (validate-shape option added in FreeCAD because of this one); stored doc now recomputes clean |
| 672 | `issue672_recompute_ghost.FCStd` | "Ghost" state: adding a circle to Sketch007 corrupts geometry unless `Support Ring` is recomputed first — recompute-order dependence, possibly not OCC |
| 985 | `issue985_uptoface.FCStd` | Extruding circles "up to face" fails (`Bnd_Box is void`); needs the up-to-face pad edit scripted |

## Found in FreeCAD work (no tracker issue)

Problems met while working on FreeCAD itself, recorded here before any fix.
Named `localNN_*`.

| # | Model | Problem | Reproduces |
|---|-------|---------|-----------|
| local01 | `local01_thickness_open_top_cup.FCStd` | PartDesign Thickness with its defaults (Skin, **Arc** join, Reversed = inward) on a plain `Part.makeCylinder(20, 10)` opened at the **top** (Face2): invalid solid, volume 16585.68 where the wall is 4423.36. The `ThicknessIntersection` twin in the same file (Intersection join) is right. Found 2026-09-27 porting upstream's recto-verso thickness (fcad `b498307318`, whose Arc test is an expected failure because of it). The same defect as `../thickness/` case `cyl_bottom_in` (which opens Face2, the **top** -- that README explains why the names are kept): every radius and height tried (r 4-40, h 5-20, offset -1 and -2) fails with the top opened and passes with the bottom opened. **Fixed** in `BRepAlgo_Loop` (a closed edge reached after its seam wire was built kept a wire of its own) -- twice, independently: here on 2026-09-28 and on the peer box as `b93149464b` (2026-09-30); the merge of 2026-10-01 kept `b93149464b`. Now 4423.36, valid. | fixed |
| local02 | `local02_fuse_shared_face_{fail,ok}.brep` + `local02_fuse_shared_face.py` | Booleans of two valid, **overlapping** solids that **share a face TShape** (same orientation in both: the second solid runs from that face back into the first). Each `.brep` is a compound of the two (base, tool) so the sharing survives. In `_fail` every op is wrong: fuse and common come back **empty, with no error**, base-tool 7457.000 (want 8218.270), tool-base 1003.506 (the whole tool; want 501.752). In `_ok` only tool-base is wrong. The two differ by ~1e-14 in their vertices, nothing else. Copying either solid (no sharing) or a fuzzy value of 1e-7 gives the right answers for both. Found 2026-10-01 as a flaky fcad test, `TestLoft.testTwoFacesAdditiveLoftCase` (upstream's issue 19183 model): a PartDesign loft from a pad's curved flank shares that flank face with the pad, and the sketch solver's last-bit noise across processes picks which side of the failure a run lands on (4-6 of 10 runs fail). Same results on upstream FreeCAD 1.1.4's OCCT 7.8.1. A plain box with a prism of its own face, inward or outward, is fine -- the curved face and the B-spline loft walls matter. Only in **non-destructive** mode (which FreeCAD always sets): `UpdateBlocksWithSharedVertices` (nd only) replaces an old vertex both faces share by a new one (`UpdateVertex` copies it, SD 35 -> 96) after it has already filled the faces' `VerticesOn`; nothing fills them again before `MakeBlocks`, which put the old vertex on the section curves, so the section edges ended on vertices the split edges no longer used and the faces they should split stayed whole (`BOPAlgo_AlertSolidBuilderUnusedFaces`). Whether the old vertex is replaced there turns on its on-curve test, hence the 1e-14 sensitivity. **Fixed 2026-10-01** in `BOPAlgo_PaveFiller::UpdateBlocksWithSharedVertices`: the faces' ON vertices are remapped to the replacements. All eight answers right; `TestLoft` 10/10. **Second symptom, same cause:** the ops also **modified their inputs** despite non-destructive mode -- the script now checks that too: FreeCAD 1.1.4 / OCCT 7.8.1 touches 6 input vertices/edges in every op of both pairs (the `_ok` pair included), the fork none. Instrumented, the pre-fix kernel widened the old shared vertices (35, 37, 61, 66) in `PutPaveOnCurve` (`BRep_Builder::UpdateVertex` on the vertex it puts, 2.12e-7 -> 2.15e-7), and with a fuzzy value of 1e-5 in `BOPTools_AlgoTools::MakeEdge`, which widens both ends of every section edge to the curve tolerance (2.12e-7 -> 7.21e-6). The rule both break: in non-destructive mode a section edge must not end on an input vertex, since `PutPaveOnCurve`, `MakeEdge`, `UpdateVertices` (MakePCurves), the micro-edge merge in `PostTreatFF` and `CorrectToleranceOfSE` all write the tolerance of section-edge vertices in place. **Open, not reproduced:** `EstimatePaveOnCurve` (which decides the replacement) tests against the curve tolerance, `PutPaveOnCurve` against it plus the fuzzy value -- never below 1e-7 -- so a shared vertex just inside the wider test could still be put on a curve unreplaced. | fixed |
| local03 | `../thickness/models/issue3_pad_thickness.FCStd` (the thickness suite's `issue3_pad_thickness`) | On the merged kernel of 2026-10-01 (`3b219f881c`, Windows, fcad `PartDesignPort` -- which sets no Immutable flags) the document's **Pad** -- a prism of a one-wire sketch (two arcs and two lines), no boolean -- came out **inside out**: volume -189356.2997, `Unorientable shape`, its Face3 a cylinder of area -9128.6. **Not the Pad: the Thickness on it edited the Pad's shape in place.** The Pad recomputed alone is a valid 9503.3915; a fresh prism given to `Part`'s `makeThickness` (faces 5, 3 and 6 removed) comes back inside out too. The arcs' parameters run past 2 pi (5.99226-6.57411); the generatrix at the arc's far end ended at u = 0.29093 on the removed cylinder, 6.57411 folded by 2 pi. Cause: `BRepAlgo_Loop::FindLoop`'s pruning builds a test face of each wire on the face's own surface and location, from the loop's own edges -- some of them the input's -- and ran `ShapeFix_Shape` on it when it was invalid; ShapeFix shifted a shared edge's pcurve by a period on that surface, i.e. the input face's own pcurve. Reached since `36d34e9d05` (sec 27.95: a periodic face whose edges lie within one period takes the plane path) -- bisected on Linux with scratch TKBool/TKOffset. It reproduces on Linux exactly, but only with `ImmutableShapeValues` off: FreeCAD's Transaction branch freezes shape values when it finds the fork, and copy-on-write protects a frozen input, so the thickness suite (which set nothing) never saw it. **Fixed** in `1e5f89c06f` (ShapeFix works on a copy of the test face; thickness results unchanged, the 712-run sweep identical) and guarded in `e22e4be177`: the thickness suite now runs unfrozen and every case checks its input is left as it was (`sector_outer_arc_input`, `sector_outer_arc_past_period_input`: a padded ring sector straddling angle 0). FreeCAD docs/TransactionLog.md sec 27.104. Windows confirmed 2026-10-01 on `6280dc5215` (TKBool/TKOffset rebuilt, fcad PartDesignPort, unfrozen): thickness suite PASS 94, issue3's Pad valid 9503.3915, Thickness valid 1241.0718. | fixed |
| local04 | `local04_fuzzy_shared_vertex.brep` (made by `local04_fuzzy_shared_vertex.py`) + `local04_fuzzy_shared_vertex.cpp`; picture `local04_fuzzy_shared_vertex.png` (made in the GUI by `local04_fuzzy_shared_vertex_picture.py`): the model, and the two radii at V to scale -- the result shapes are the same before and after the fix, so there is no before/after pair to show | A **non-destructive** Boolean with a **fuzzy value** writes an **input** vertex's tolerance. `UpdateBlocksWithSharedVertices` asks `EstimatePaveOnCurve` whether an old vertex shared by the two faces of an FF pair lies on their section curve, with the curve's tolerance only (`IsVertexOnLine`: distance <= max(2 (tolV + tolC), 1e-6)); a vertex that does is replaced by a new one before anything writes to it. `PutPaveOnCurve` later asks the same question with tolC + fuzzy. A vertex whose distance falls between the two is never replaced, is put on the curve, and `PutPaveOnCurve` raises its tolerance in place (`BRep_Builder::UpdateVertex`) -- the input's own vertex. The model: two prisms whose faces share ONE vertex V (one TShape) at a reflex corner of both, V on A's plane z=0 and 1.74e-6 off B's plane (tilted 10 deg), tolerance 1.94e-6; the section line passes V at 1e-5. Fuzzy 3.162e-7: the estimate's threshold is 9.8e-6 (off the curve), the put's 1.04e-5 (on it), and V goes 1.94e-6 -> 1e-5 in fuse, cut, common and the reversed cut; at 2.512e-7 both agree and nothing changes. A sweep over tilts 5/10/20 deg and distances 1e-5/1e-4 hits 6 of 132 (pair, fuzzy) points, each in the narrow band 2 (tolV + tolC) < d <= 2 (tolV + tolC + fuzzy), tolC there being the curve's tangential tolerance, which grows with the fuzzy value as the angle shrinks. Why only this shape: in non-destructive mode any VF/EF update of an input vertex replaces it first, so the only candidates are vertices shared ON a face of each argument; and at a convex corner the section line crosses V's own edge at the same distance, which splits or snaps first. **FreeCAD's Part booleans do not reach it:** given a fuzzy value they copy the tools (`BRepBuilderAPI_Copy`, an old workaround), which also breaks the sharing -- hence the C++ driver. **Fixed** in `2f91019c73`: `EstimatePaveOnCurve` adds `myFuzzyValue` to the tolerance, as `PutPaveOnCurve` does (it is never below 1e-7, so non-fuzzy runs move by that much too). The same 528 runs: no input touched, every result's volume and validity identical to before; thickness suite PASS 94, local02 0 wrong; FreeCAD on it (Windows, PartDesignPort): Python 3195 OK, ctest 750/750. | fixed |

## No usable repro

| # | Problem |
|---|---------|
| 937 | Fillet crash: `BRep_Tool::CurveOnSurface` null deref under `ChFi3d_Builder::StartSol` (OCC 7.7.2). Stack trace only, reporter never posted the file |

## Themes

- **Seam edges** are the dominant killer: split through seam (#64), fillet
  touching seam (#273, upstream 28354), pocket with overlapping seams
  (#280). Same family as the thickness-chain seam handling in
  `../thickness/`.
- **Fillet/chamfer** is the biggest cluster (273, 309, 360, 474, 523, 631,
  876, 937, 962) — `ChFi3d` crashes and invalid outputs.
- **Thickness/offset** (346, 363, 521, 613, 662) — same `BRepOffset` code
  paths as `../thickness/`; #363 is the best next target.
- Three models **crash a debug build on plain recompute** (333/334, 360,
  474) — good candidates for crash-isolation work regardless of fix.

Source data (issue JSON + all raw attachments) fetched 2026-08-01;
re-fetch with the GitHub API: `issues?labels=occ&state=all`.
