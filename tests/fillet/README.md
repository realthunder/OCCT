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

## A corner whose faces are split in coplanar pieces (realthunder/FreeCAD#962, 2026-10-02)

The model is a PartDesign body without Refine: its walls are kept as several
coplanar faces, split by the booleans that made them. Its Fillet (radius 0.8,
twelve edges) was invalid, and five of the edges alone were too
(`models/issue962_pocket002.brep` is the Fillet's input, the body's
Pocket002). Four are the outer edges of slots cut down to a floor -- the
edge ends on the floor, which the fillet extends (`PerformOneCorner`'s
OnSame case) -- where the slot's wall is two faces, its piece at the outer
face only 0.692 or 0.787 wide. The fifth, edge 50, ends under a chamfer, and
the wall beside it (y=10.2) is two faces split at that end.

`PerformOneCorner` builds such a corner from the face Fad holding the
fillet's common point on the end face Fv, the face Fop on the other side, and
the edge Arcprol between Fv and Fop, which it extends to the fillet. Three
assumptions about those broke on split faces. Fix `cf96757c36`:

- **The fillet wider than Fad's piece at the vertex**: its common point lies
  on the next piece of the wall, so the arc it is on is not an edge of the
  vertex, and the first other edge of Fv at the vertex -- the narrow piece's
  own bottom edge -- was taken as Arcprol. Volumes off by +11 to -24 (the
  slots), invalid from radius 0.7 on. Arcprol is now the edge between Fv and
  Fop (or a face tangent to Fop); and the edges of Fv from the vertex to the
  common point's arc, which nothing else cut away (the narrow piece survived
  whole, its bottom edge dangling in the floor), go the way of the spine's
  edge.
- **Arcprol not in Fop but in Fop's tangent neighbour** (edge 50): its
  orientation in Fop defaulted to FORWARD, and Fop lost the fillet's line --
  an open wire. It is now taken from the neighbour, in the shell, and carried
  into Fop's frame.
- **The extension along the split itself**: a fin padded flush with the
  block's side (`fin_on_block_*`) -- the most ordinary of the three. The
  fillet's line on the fin's face ends on the edge between the fin's and the
  block's face, `IntersUpdateOnSame` reset the point off that edge, and the
  extension, which runs along it, was given to the fin's face: a
  self-intersecting wire at every radius. The point stays on the edge, and
  the extension bounds the block's face.

Upstream's files at the fork's base give the same results as the fork before
the fix.

| Case | What it covers |
|------|----------------|
| `slot_split_{wall,outer,both}_r*` | a block with a slot down to a floor, its wall split 0.7 from the outer face, its outer face split below the floor, or both: each must give the unsplit slot's volume, radius 0.5 to 2 |
| `issue962_e10{1,2,4,5}_r0.8` | the model's slots: all four at 11580.7739 (identical slots) |
| `issue962_e50_r0.3`, `issue962_e50_r0.8` | the edge under the chamfer |
| `fin_on_block_r*`, `fin_arc_on_block_r*` | a fin flush with the block's side, flat or with a top arc tangent to its outer face; the flat one must remove what the slot does |
| `fin_on_block_wall_r0.3` | the same fin, its wall split as the slot's, narrower than the piece |

Checked against the fix: every edge of 22 shapes -- the inputs of the fillet
and chamfer features of #273, #474, #631, #876 and #962, #309's and #523's
shapes, and small models of each case -- at radii 0.3, 0.8 and 2: 112
results go from invalid to valid (#962's Pocket002 and Pad004 among them),
the other 3590 are the same as before. FreeCAD's `TestPartApp` (139) and
`TestPartDesignApp` (77) pass.

## A seam cut at both ends by two corners (realthunder/FreeCAD#962, 2026-10-02)

With the five edges above fixed, #962's Fillet was still invalid with all
twelve: edges 36 and 49, the foot of the block where the arm meets it
(x=38.5, z -9.75..-3.75), were each valid alone, and the Fillet without
either one was valid, but the two together came out invalid, 7225 too large.
The bottom is two coplanar faces, the arm's and the block's, split along
x=38.5 from y=27.2 to y=10.2. `PerformIntersectionAtEnd` builds both corners
across that seam: e36's convex fillet cuts it at y=27.186, e49's concave one
meets its line at y=10.136 -- past the seam's end. Each corner gave the seam
a point interference at its parameter, the second one outside the edge's
range; alone, the face rebuild turns that into a lengthened edge, but with
the other corner's cut on the same edge it built both -- the cut edge
(27.186..10.2) and the lengthened one (10.136..27.2) -- in both faces. It
turns on the seam's direction: with the seam running the other way the
point falls before the edge's start, and that composes.

Fix `b38f0d910c`: across a seam of tangent faces, a point past the seam's end at
the corner's vertex is not given to the seam; the piece of the seam's line
from the vertex to the point becomes a curve of its own on both faces (the
way the corner extends a face's edge when it has to extend the face), and
the seam is left whole for the other corner. Only tangent seams: on a sharp
edge (e49's top, where the arm's top meets the block's side) the curve did
not close the faces, and the lengthening works there.

| Case | What it covers |
|------|----------------|
| `issue962_e36_e49_r0.8` | the pair: the sum of the two alone |
| `issue962_fillet_r0.8` | the Fillet, all twelve edges |
| `arm_on_block_foot_r*` | an arm fused to a block, the seam between their bottoms from y=27.2 to 10.2: the pair at radius 0.5, 0.8 and 1.5 gives what the two give alone |

Checked against the fix: the every-edge sweep (22 shapes, plus this slab)
is unchanged, 3762 results -- one edge at a time never met this; FreeCAD's `TestPartApp` (139) and `TestPartDesignApp` (77)
pass. The fork before the fix fails the five new cases.

An aside found on the way: edge 56 alone reports a volume 5.04 smaller than
its input where its mirror image, edge 50, reports 2.35 -- the faces match
change for change, and a mesh of the result gives 2.2. The corner's blend is
a B-spline far larger than its trimmed piece, and the volume integration
over it is off; the shape is right.

## A fillet ending on a wall kept in coplanar pieces (2026-10-02)

The last of the corner cases above, left open by the fix there: the fin
padded flush with the block's side, its wall split 0.7 from the outer face
as the slot's (`fin_on_block_wall_*`). Each of the two alone works -- the
slot with the split wall, the fin with the whole one -- and so did the two
together up to radius 0.5; from 0.7 on, `StdFail_NotDone`, a walking
failure.

With the wall split, the fillet's support at the spine is the wall's
narrow piece, and its line on the wall lies on the wide piece, whose edge
on the floor stops at x=9.3, short of the spine's end vertex (10,3,5). The
walk along the spine stops at z=5, where the fin's face ends; on the slot
the outer face goes on below the floor, and one walk takes the fillet past
it into the extension. `StartSol`, restarting the walk there, takes an
edge of the end that does not touch the spine's end vertex for an
obstacle: the floor became a face to roll along, and the walk failed.
Behind it, `PerformOneCorner` counted only the fin's common point as on
the vertex, took the block's side for the end face instead of the floor,
and threw "bouchon non ecrit" (the cap is not written).

Fix `2ac4fee0c6`: `ChFi3d_EdgeOnSplitToVertex` follows such an edge along
the end face, through faces tangent to its own across edges of the way, to
the vertex. `StartSol` takes it for the end, as the whole wall's edge, and
the corner counts the common point as on the vertex, classifying the
vertex with the edge the arc continues as there. From there on the corner
is the slot's: the narrow piece's edge on the floor goes with the spine's
edge.

| Case | What it covers |
|------|----------------|
| `fin_on_block_wall_r{0.3,0.8,2}` | the fin flush with the block, its wall split as the slot's: what the unsplit fin takes |

Checked against the fix: every edge of 25 shapes (the sweep's 22, the
slab, the split fin and the split slot) at radii 0.3, 0.8 and 2: 2 results
go from failing to valid -- the split fin's own edge at 0.8 and 2 -- and
the other 3919 are the same. FreeCAD's `TestPartApp` (139) and
`TestPartDesignApp` (77) pass.

Still open (XFAIL), on the split slot as on the split fin: a radius
within about 1e-4 of the piece's width, 0.7. The fillet's line on the wall
then runs along the split itself, its whole length; the start of the walk
wants a point inside both faces (`BRepBlend_Walking::PerformFirstSection`),
and it is on the boundary of both pieces everywhere. 0.69999 and 0.701
work; 0.7 and 0.70001 fail to start, 0.7001 comes out invalid. It wants the
walk to take the tangent pieces as one face, a larger change.

## A corner plate whose boundary curve missed its first surface (realthunder/FreeCAD#962, 2026-10-02)

#962's edge 36 alone -- the arm's side against the block's side at its
foot, 21 degrees apart -- failed below radius 0.42 with a faulty vertex at
its top (38.5,27.2,-3.75). Above the arm the block goes on, its side
y=27.2 two coplanar faces split at the arm's top, and five faces meet at
the vertex. `PerformIntersectionAtEnd` hands such a corner to
`PerformMoreThreeCorner`, which fills it with a `GeomPlate` patch. With the
fillet that narrow (its lines 0.19 r off the edge), one of the patch's
boundary curves had no projection at all on the plate's first surface, and
`GeomPlate_BuildPlateSurface::ProjectCurve` read the bounds of the first
projected piece without asking whether there was one:
`Standard_NoSuchObject`. A curve without a continuous projection is what
`Perform` falls back on other first surfaces for -- the plate's
approximation, then the plane -- and the plane takes all four curves.

Fix `2b9df66c48` (TKGeomAlgo): no projected piece is a null curve, as no
continuous one is. The patch fits within 1.4e-3, the fillet is valid at
every radius, and the volume it takes grows with r^2 as the sliver's does
-- measured on the result cut down to a box around the edge, since the
whole shape's volume (11584) is too large for GProp to see 0.002 in.

| Case | What it covers |
|------|----------------|
| `issue962_e36_r{0.1,0.3}` | the edge in #962's Pocket002: valid (no volume: the change is below GProp's precision on the whole shape) |
| `arm_on_tall_block_r{0.1,0.3,2}` | the same corner made small: the arm fused to a block that goes on above it, the block's side split at the arm's top |

Checked against the fix: the every-edge sweep, 25 shapes: 2 results go from
an exception to valid (the edge, and the same corner in #962's Chamfer's
input), 3914 are the same, and 5 go from an exception to an invalid shape
-- #474 Fillet003's edge 52 at every radius, and #876's front edge at
radius 2 in both of its fillets. Those corners now get past the plate and
fail further on (in #474 a second plate fits only to 0.058); the sweep
holds 140 other invalid results, which is how this API fails as often as
it throws. PartDesign's Fillet reports either as an error; `Part::Fillet`
and `makeFillet` hand the invalid shape back. FreeCAD's `TestPartApp`
(139), `TestPartDesignApp` (77) and `TestSurfaceApp` pass.

Found on the way, open (XFAIL): the same corner with the block's side one
face, `arm_on_tall_block_whole_r0.5`, is invalid at every radius, with the
fork before the fix as after -- the arm's top is extended past the block's
wall by the fillet's reach (0.1875 r) and its wire crosses itself.

## A corner plate folded to stay tangent (realthunder/FreeCAD#962, 2026-10-02)

#962's edge 56 alone took 5.04 of volume where its mirror image, edge 50,
took 2.35 -- and a mesh of the same result said 3.5, GProp's adaptive mode
2.0. The shape was not right. Its foot is the five-face corner at
(38.5,27.2,-3.75), which `PerformMoreThreeCorner` fills with a `GeomPlate`
patch held G1 -- tangent -- to the stripe it closes. The plate missed its
own boundary by 0.092 on a fillet of 0.8, and folded over 1/30 of its area
to meet its constraints; the approximation, allowed ten times the plate's
miss, strayed 0.48, so the corner's edges were stored with tolerances of
0.48 to 1.18 (BRepCheck passes them for it), and a pcurve of the patch
looped over itself. Edge 50's plate missed by 0.011 and strayed 0.19: no
fold, but tolerances of 0.195 on a fillet of 0.8.

It is no rare plate. Over the every-edge sweep, of the 495 corner plates
built, 413 missed their boundaries by more than 1e-3 and 87 by more than
0.1, a few by thousands; the code's own comment says G1 constraints "very
often cause unpredictable undulations". A G0 plate -- positions alone --
fitted edge 56's boundary within 3.6e-4.

Fix (`PerformMoreThreeCorner`): a plate that misses its boundary by more
than 1e-2 is built again from the same boundaries at G0, and taken if it
fits better; its edges are then not marked tangent to the stripes. A
crease along the stripe instead of a fold. Edge 56 takes 2.433 by every
measure now (GProp's default and adaptive modes agree to 1e-4), edge 50
2.420, tolerances 0.014 and 0.009.

| Case | What it covers |
|------|----------------|
| `issue962_e56_r0.8` | edge 56: e50's volume, tolerances under 0.05 |
| `issue962_e50_r0.8` | edge 50 now also held to tolerances under 0.05 (it had 0.195) |
| `issue474_f003_e6_r{0.3,0.8,2}` | #474's Fillet003 input, edge 6: invalid at every radius (tolerances to 36, 244 too much volume at r 2), valid now |

Checked against the fix: the every-edge sweep, 25 shapes, with each
result's largest edge or vertex tolerance: 17 results go from invalid to
valid, 107 valid ones end with less than half their tolerance (#474's edge
17 at 0.8: 41.6 to 0.002; #273's edge 16 at 2: 9.6e5 to 0.18), none gets a
larger one, and 2 go from valid to invalid: #962's edge 33 at 0.3, in the
Fillet's input and the Chamfer's. That corner's G1 plate was no better --
its outline crossed itself ten times and its surface folded, and the
"valid" fillet of radius 0.3 removed 9.4 of volume where 0.6 is due --
and the G0 one, which fits ten times better, has two of its boundary
curves meeting head to head, which BRepCheck now rejects. With the rebuild
at 1e-3 instead of 1e-2, 20 results became valid and 156 tighter, at the
price of rebuilding most corners; 1e-2 leaves the merely mediocre G1
plates (edge 36's, 3.3e-3) as they were. FreeCAD's `TestPartApp` (139),
`TestPartDesignApp` (77) and `TestSurfaceApp` pass.

## Known broken (XFAIL)

| Case | Symptom |
|------|---------|
| `slot_split_wall_r0.7`, `fin_on_block_wall_r0.7` | `StdFail_NotDone` (no start for the walk): the radius equal to the wall's piece, the fillet's line on the wall running along the split |
| `arm_on_tall_block_whole_r0.5` | invalid at every radius: #962's edge 36 corner with the block's side one face; the arm's top, extended under the block's wall, crosses itself |
| `seam_end_bottom_r2` | invalid, where its mirror image across z=5 (`seam_end_top_r2`) is valid: at radius 2 the fillet's end reaches the far end of the cylinder's face too (u = pi/2, the vertex at (0, 2)) |
| `mirror_top_r2` | `StdFail_NotDone`, the same reach on the side without the seam |

## Pictures

`models/pictures/<case>.png` shows, for each case a fix turned from failing
to passing, three results side by side -- upstream's `TKFillet` files at the
fork's base (`91be8c4c71`), with the files of other toolkits a fix changed
(`UPSTREAM_EXTRA` in `pictures/cases.py`), the fork just before the fix, the
fork now --
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
