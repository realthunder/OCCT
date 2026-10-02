# Fillet regression tests

Regression suite for the fork's fixes to fillets (`TKFillet`: `ChFi3d`,
`BRepFilletAPI_MakeFillet`), the cases taken from the fillet issues of
`../occ-issues/`. Same shape as `../thickness/`: every fix gets cases here,
and the fixes are shown before and after in `models/Fillet.md`.

## Running

Any FreeCAD build linked against this OCCT tree:

```
FreeCADCmd tests/fillet/run_tests.py
```

The suite runs with FreeCAD's shape values unfrozen (`ImmutableShapeValues`
off), as any build that does not freeze them runs; `FILLET_FREEZE=1` runs it
frozen. Every case also checks that its input comes back as it went in
(valid, the same volume).

Exit code 0 = no regression (every expected-pass case passed). Known-broken
cases are marked `XFAIL` in the script; when one starts passing the suite
prints `UNEXPECTED-PASS` -- promote it to `pass` with its reference volume and
update this file. A crash ends the run before its summary line: that is a
failure too.

## A seam left on a piece of a closed face (realthunder/FreeCAD#523, 2026-10-02)

The model: a 10 box with a three-quarter cylinder of radius 2 at the corner
on the z axis, joined by Part's Connect (`models/issue523_fillet_explodes.FCStd`;
`models/issue523_box_cylinder.brep` is its Connect, the fillet's input). The
cylinder's seam is the line x=2, y=0 where it meets the box's side, and the
boolean kept it as an edge with **both** pcurves on the cylinder (u=0 and
u=2 pi) although the three-quarter face holds it once. A radius-1 fillet on
the box edge along y=0 that ends at that line came out invalid
("Unorientable shape"): the face's wire did not close in 2D, the end cap's
pcurve a period (u 0..0.52) away from the face (u 1.57..6.81).

Cause: `ChFi3d_Builder::PerformOneCorner` read the pcurve of the common
point's arc on the face at the fillet's end with `BRep_Tool::CurveOnSurface`,
which, for an edge closed on the face, picks a pcurve by the orientation of
the edge it is handed -- here the orientation the fillet data carries, not
the edge's in the face. `IntersectMoreCorner` repeats the code. Fix
`16df68d224`: `PCurveInFace` takes the edge as it occurs in the face when
the edge is closed on the face but not really closed
(`BRepTools::IsReallyClosed`), for every pcurve those corners read of an
edge not taken from the face itself. Upstream master has the same code.

| Case | What it covers |
|------|----------------|
| `issue523_document` | the captured model recomputes, its Fillet valid at 1092.5263 |
| `seam_end_top_r*`, `seam_end_bottom_r*` | the fillet on the box edge ending on the seam, top and bottom, radius 0.5 and 1 (and 2 on top) |
| `mirror_*` | the same fillet on the edge across the box's diagonal, which does not end on a seam: the seam cases must give its volumes |
| `every_edge_r1_e*` | every edge of the shape at radius 1: nothing else moved |

Checked against the fix: every one of the shape's 15 edges at radius 0.5, 1
and 2 gives the same result as before except the two edges ending on the
seam; the other fillet models of `../occ-issues/` recompute as before, but
for #876's first fillets, whose vertex blend now starts on the right period
(volume -0.0173; it meets its neighbour cylinders within 1.51 deg, 1.91
before); FreeCAD's `TestPartApp` (139) and `TestPartDesignApp` (77) pass
before and after.

## A fillet line collapsed to a point (found on #523's shape, 2026-10-02)

A fillet of radius 2 on the cylinder's top or bottom arc -- the radius of the
cylinder -- has no line on the top face, only a point, and the stripe's curve
there is null. `PerformOneCorner` (through `IntersUpdateOnSame`) and
`IntersectMoreCorner` dereferenced it: SIGSEGV, the process gone. Fix
`214d858dc3`: no intersection is reported, and the fillet fails
(`StdFail_NotDone`). Upstream master has the same code.

| Case | What it covers |
|------|----------------|
| `arc_top_r2_no_crash`, `arc_bottom_r2_no_crash` | the fillet is refused, the process lives, the input is as it was |

## Known broken (XFAIL)

| Case | Symptom |
|------|---------|
| `seam_end_bottom_r2` | invalid, where its mirror image across z=5 (`seam_end_top_r2`) is valid: at radius 2 the fillet's end reaches the far end of the cylinder's face too (u = pi/2, the vertex at (0, 2)) |
| `mirror_top_r2` | `StdFail_NotDone`, the same reach on the side without the seam |

## Pictures

`models/pictures/<case>.png` shows, for each case a fix turned from failing
to passing, three results side by side -- upstream's `TKFillet` files at the
fork's base (`91be8c4c71`), the fork just before the fix, the fork now --
each whole and zoomed on the fillet's end. `models/Fillet.md` walks
through them. `pictures/make_pictures.sh` makes them again, on Linux or
macOS; add a case to `pictures/cases.py`, with its stage, to picture it.

## Layout

```
tests/fillet/
  README.md          this file
  run_tests.py       the suite (FreeCADCmd script)
  models/            captured models and shapes
    Fillet.md        the fixes, before and after, in pictures
    pictures/        before and after, one PNG per fixed case
  pictures/          the tools that make them (make_pictures.sh)
```
