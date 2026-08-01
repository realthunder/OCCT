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
