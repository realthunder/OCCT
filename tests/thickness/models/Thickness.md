# Thickness in the OCCT fork, before and after

FreeCAD's Thickness -- `Part::Thickness`, `PartDesign::Thickness`,
`Shape.makeThickness()` -- hollows a solid by removing faces and giving the
rest a skin. Underneath it is OCCT's `BRepOffsetAPI_MakeThickSolid`, and the
fork (this repository, branch `LinkVibe-801`) carries a chain of changes to
it so that a *concave* face can be removed, which upstream does not support
(realthunder/OCCT#1-#4). From 2026-09-30 to 2026-10-01 that chain was
measured against upstream for the first time and the thickness failures were
chased down, one cause at a time. This page shows what each fix did, in
pictures. The build log -- every cause, what was tried, the hand-worked
volumes, the gates -- is FreeCAD's docs/TransactionLog.md sec 27.88 to
27.96 and sec 27.100 to 27.102; "sec" below means a section there.

Where things are:

- The fork's suite: `tests/thickness/run_tests.py` (`FreeCADCmd
  tests/thickness/run_tests.py`; PASS 82) and its `README.md`, the
  case-by-case reference.
- The pictures: `pictures/<case>.png` beside this page, one per case of the
  suite that a fix turned from failing to passing.
- The tools that made them: `tests/thickness/pictures/` (`make_pictures.sh`,
  see "Making the pictures" below).
- FreeCAD's own cases: `parttests.regression_tests` in FreeCAD's `src/Mod/Part`, one
  test per step (`test_thickness_*`).

## Reading a picture

Each picture is one case: a solid, the face removed (`FaceN`, 1-based as in
`shape.Faces`), the direction (`+1` outward, `-1` inward, thickness 1), the
join (Arc or Intersection) and whether intersection mode is on; the expected
volume, or "a valid solid" where no reference exists; and the section of
FreeCAD's TransactionLog.md that fixed it.

Three columns:

- **upstream OCCT 8.0.1** -- upstream's eleven files of the fix chain at the
  fork's base (`91be8c4c71`), compiled into scratch `TKBool`/`TKOffset`
  libraries and preloaded, the rest of the fork as it is.
- **fork before** -- the fork's `TKBool`/`TKOffset` sources at the commit just
  before the fix (named in the heading), the same way.
- **fork after** -- the fork now.

Two rows: the result seen from the removed face's side, and the result cut
open -- the half towards the camera cut away by a plane through the removed
face's centre, the cut faces orange, so the walls show their thickness.

Under each column heading: "valid, volume V" in green, or what was wrong in
red. "Unhollowed" is a result whose volume is the input's: the opening was
never cut, though the solid may be valid. "Inside out" is a negative volume;
such a solid is drawn black (its back faces), and its cut is done by a clip
plane with no caps, as is any cut the boolean refuses ("clipped"). Where the
thickness threw, the input is drawn greyed. Red edges are edges of the result
not used exactly once each way round by its faces (free, doubled, or turned);
red faces fail `isValid()`; faces of no area are left out, and the picture
says how many.

Upstream's column is one run of many: upstream visits offsets in hash order
and its wrong results vary from run to run (sec 27.88). Of the first 46, two came
out different on a second run (`tshape_bar_top_out`, `filletbox_end_out`),
wrong both times. The fork's columns are the same every run.

A result can be right only when it is a valid solid with one closed shell and
the expected volume -- or, where no cavity can reach the removed face (sec
27.100), two closed shells, the skin and a void; the volumes are upstream's
where upstream is right, and otherwise worked out by hand (FreeCAD's
TransactionLog.md has each derivation).

## How it stood, and how it stands

The sweep behind all of this (sec 27.89): 712 thickness runs
-- 16 solids, every face, +/-1, intersection off and on, the Arc and the
Intersection join -- against both builds. A run counts as right only as
above, with the volume plausible for a skin where no reference exists.

| After | Fork worse than upstream | Fork better | Fork right, Arc join (int. off / on) | Fork right, Intersection join (off / on) | Upstream right (four modes) |
|---|---|---|---|---|---|
| start (sec 27.89) | 144 | 51 | | | |
| sec 27.89 | 26 | 53 | | | |
| sec 27.90 | 9 | 59 | 123 / - | | 101 (Arc, off) |
| sec 27.91 | 0 | 78 | 135 / 135 | 129 / 113 | 106 / 107 / 110 / 111 |
| sec 27.92 | 0 | 110 | 135 / 135 | 137 / 137 | 106 / 107 / 110 / 111 |
| sec 27.93 | 0 | 122 | 141 / 141 | 137 / 137 | 106 / 107 / 110 / 111 |
| sec 27.95 | 0 | 138 | 149 / 149 | 137 / 137 | 106 / 107 / 110 / 111 |
| sec 27.96 | 0 | 162 | 149 / 149 | 149 / 149 | 106 / 107 / 110 / 111 |
| sec 27.100 | 0 | 164 | 150 / 150 | 150 / 150 | 106 / 108 / 110 / 112 |
| sec 27.101 | 0 | 164 | 150 / 150 | 150 / 150 | 106 / 108 / 110 / 112 |
| sec 27.102 | 0 | 164 | 150 / 150 | 150 / 150 | 106 / 108 / 110 / 112 |
| sec 27.103 | 0 | 164 | 150 / 150 | 150 / 150 | 106 / 108 / 110 / 112 |
| sec 27.104 | 0 | 164 | 150 / 150 | 150 / 150 | 106 / 108 / 110 / 112 |
| sec 27.105 | 0 | 164 | 150 / 150 | 150 / 150 | 106 / 108 / 110 / 112 |

The "right" counts are of the 150 runs a mode that are in scope; sec 27.90
set the scope (a face whose removal leaves the shell in pieces is out) and
sec 27.91 rebuilt the pocketed box (walls 2.5, not twice the thickness), so
counts before sec 27.91 do not compare with those after it. The last row
counts the holed cone's sealed void right (sec 27.100), which upstream's
intersection mode gives -- by that rule the fork was worse than upstream in
those two runs before it.

Of the 69 pictured cases, 27 are ones upstream gets right and the fork had
broken -- the chain's casualties (sec 27.89, 27.90, part of 27.91, the holed
cone's top with intersection on, sec 27.100, a short box's and a box's bottom
alone, sec 27.101 and 27.102, a pocket's open shell, sec 27.103, and the
input left inside out, sec 27.104). The other 42 fail upstream too: the fork
now does better than upstream there.

Nothing in the suite fails now, and every run of the sweep in scope is
right. Out of scope, sec 27.101 answers the faces left in pieces where each
piece is a plain plate or disc, with either join (sec 27.102), and sec 27.103
the pieces that are pockets and bosses: every one checks by hand. The one
refusal left in the sweep is the torus's face, and it is right: with its one
face removed no face stays (sec 27.105, which refuses the sphere's too).

## The fixes

### Sec 27.88: the same result every run

Not pictured: no single picture shows it. With intersection on and the Arc
join, `BuildOffsetByArc` visited its offsets in hash order -- a shape's hash
is its TShape's address -- and the result depended on the order: up to four
different results in eight runs of one input. The offsets are now taken in
the shape's topological order. The suite's `arc_inter_*_same_every_run` cases
guard it.

### Sec 27.89: the chain's own regressions

Upstream passes every case here; the fork, before, broke all of them.

*A closed edge found twice* (`BRepAlgo_Loop::FindLoop`). On a periodic face
the seam wire replaced the closed-edge wires found before it, but the search
found the other closed edge again afterwards and kept it as a wire of its
own: the inner wall had wires [1, 3, 1] and the floor went missing -- the
open shell, seen as the black disc where the floor should be.

![cyl_bottom_in](pictures/cyl_bottom_in.png)
![hole_bottom_in](pictures/hole_bottom_in.png)

*One seam wire per face* (`FindLoop`). A removed cylinder leaves a wall at
each end, each closed by its own piece of the seam; the loop allowed one seam
wire per face, so the other wall was built from two bare circles -- the face
drawn red, which fails `isValid()`.

![hole_outer_out](pictures/hole_outer_out.png)
![hole_outer_in](pictures/hole_outer_in.png)
![hole_inner_out](pictures/hole_inner_out.png)
![hole_inner_in](pictures/hole_inner_in.png)
![boxhole_hole_out](pictures/boxhole_hole_out.png)
![boxhole_hole_in](pictures/boxhole_hole_in.png)

*A const edge losing its orientation* (`BRepAlgo_Loop::Perform`). The
concave-face pass replaced the const edges with FORWARD copies, and a seam
wire closed on two circles running the same way (inward); outward, the
re-found closed edge of the first cause.

![pocket_bottom_out](pictures/pocket_bottom_out.png)
![pocket_bottom_in](pictures/pocket_bottom_in.png)

*The offset of an ellipse stretched past its ends*
(`BRepOffset_Inter3d.cxx`, `ExtentEdge`). The offset of an ellipse is a
closed B-spline that is not periodic, and the context extension stretched it
100 lengths past its ends, evaluated there at 1e34: the bottom came back
unhollowed, the top as two shells with a face of no area (outward, of
negative volume). A closed edge on a
curve that is not periodic is no longer extended.

![ellipse_bottom_out](pictures/ellipse_bottom_out.png)
![ellipse_bottom_in](pictures/ellipse_bottom_in.png)
![ellipse_top_out](pictures/ellipse_top_out.png)
![ellipse_top_in](pictures/ellipse_top_in.png)

Not pictured: a crash in `BuildSplitsOfTrimmedFaces` (upstream code, reached
only with the fork's edges), now guarded; the suite's
`pocket_inter_join_no_crash`.

### Sec 27.90: past concave corners, and a concave face removed

*The corner piece of an arc face* (`BRepAlgo_Loop`). Outward with the Arc
join past a concave corner, the arc face along one edge is cut by the arc of
the next; the fork's loop keeps every piece of a cut edge (it needs them for
a concave removed face), and the piece beyond the cut, closed by the arc's
own end, came out as a face of its own -- the red sliver at the corner. A
wire through a piece beyond the edge's own span now loses to one that stays
inside. Upstream is right on the L-box and the T; on the pocket box it
returns the input, invalid.

![lbox_arm_end_out](pictures/lbox_arm_end_out.png)
![tshape_arm_end_out](pictures/tshape_arm_end_out.png)
![pocketbox_wall_out](pictures/pocketbox_wall_out.png)

*A neighbour skipped before its edge was renewed*
(`BRepOffset_MakeLoops::BuildFaces`) and *a stretched edge lying on another*
(`FindLoop`). The T's bar top beside the post -- a concave face -- removed,
outward: a corner sphere kept an old edge, and a piece of the removed face's
stretched edge lay on an arc face's tangent line. Upstream returns it inside
out, with two shells; the fork threw.

![tshape_bar_top_out](pictures/tshape_bar_top_out.png)

### Sec 27.91: inward past a concave top, and intersection mode

*Where the stretched removed face crosses an edge* (`BRepAlgo_Loop::Perform`).
The T's bar top removed, inward: the stretched removed face crossed the far
end wall's inner edge above the bar's inner top, and that crossing was taken
for the edge's own end; the corner above the bar's inner arc came out as a
face of its own and the build failed. Upstream returns the input unhollowed.

![tshape_bar_top_in](pictures/tshape_bar_top_in.png)

*An edge two faces share, trimmed twice* (`BRepOffset_MakeOffset.cxx`,
`TrimEdges`). The guard meant to trim each new edge once, an indexed map's
`Add()`, is never 0; with intersection on and the Intersection join the
second pass cut a section short, and the fork threw where upstream is right.

![lbox_top_inter_join_out](pictures/lbox_top_inter_join_out.png)
![lbox_top_inter_join_in](pictures/lbox_top_inter_join_in.png)
![tshape_back_inter_join_out](pictures/tshape_back_inter_join_out.png)
![tshape_back_inter_join_in](pictures/tshape_back_inter_join_in.png)

*A wire running once round a periodic face* (`FindLoop`). A cone with a
through hole, its bottom removed, inward, intersection on: a circle beyond
the seam's span came out as a face of no area, a second shell (the circle
drawn above the hole in the cut). Upstream is right.

![conehole_bottom_inter_in](pictures/conehole_bottom_inter_in.png)

### Sec 27.92: the Intersection join -- inside out, a blind floor, concave faces

Upstream fails every case here.

*A thick solid inside out* (`MakeThickSolid`). The quilt orients a shell by
the faces it meets first, so which way the result faced depended on the order
the shells met: valid solids of negative volume, drawn black. The solid is
now oriented by classification (`BRepLib::OrientClosedSolid`).

![tshape_bar_top_right_join_in](pictures/tshape_bar_top_right_join_in.png)
![pocket_floor_join_out](pictures/pocket_floor_join_out.png)
![pocket_floor_join_in](pictures/pocket_floor_join_in.png)

*A blind floor's section the wrong way round* (`ContextIntByInt`). A blind
hole's floor meets the wall at a concave edge, and the section on the wall's
offset came out reversed: a band closed on two circles running the same way.

![blindhole_floor_join_out](pictures/blindhole_floor_join_out.png)
![blindhole_floor_join_in](pictures/blindhole_floor_join_in.png)

*Concave removed faces with intersection on* (`MakeOffsetShape`). The splits
of the offset faces do not cut a neighbour's section where the rim needs it;
a shape whose removed face meets a neighbour at a concave edge is now built
as with intersection off, which gets all of these right.

![lbox_notch_wall_inter_join_out](pictures/lbox_notch_wall_inter_join_out.png)
![tshape_post_wall_inter_join_in](pictures/tshape_post_wall_inter_join_in.png)
![pocket_wall_inter_join_in](pictures/pocket_wall_inter_join_in.png)

### Sec 27.93: a removed face tangent to its neighbours

A box with its vertical edges filleted, an end or a side face removed. The
fillets' offsets run parallel to the removed face: inward they never meet it
(unhollowed), outward they meet it far round the cylinder -- a lip in the
opening (the side face: a valid solid, but 354.34 where 388.82 is right) or a
failure. A tube round the tangent edge now closes the gap, turning into the
removed face, with an eighth of a sphere at each corner outward. Upstream
fails all four.

![filletbox_end_in](pictures/filletbox_end_in.png)
![filletbox_side_in](pictures/filletbox_side_in.png)
![filletbox_end_out](pictures/filletbox_end_out.png)
![filletbox_side_out](pictures/filletbox_side_out.png)

### Sec 27.94: the loop's minimal wires by angle

Not pictured: FreeCAD's WireJoiner angle walk ported into `BRepAlgo_Loop`.
It changed no result (the walk and the old search were checked against each
other on every face of the suite and the sweep); it is what sec 27.95 uses on
a curved face.

### Sec 27.95: a fillet removed

The same box with a fillet removed -- a curved face tangent to the planes
beside it. Outward it threw; inward it came back valid at 647.56, more than
the input (459.40). The tubes
and corners are built on the fillet's cylinder; the loop on that periodic
surface, whose edges all lie within one quarter turn, walks the angles as on
a plane; and a piece of the floor's offset left hanging on the shell is
dropped. Upstream fails both. The cut is in plan, to show the opening.

![filletbox_fillet_in](pictures/filletbox_fillet_in.png)
![filletbox_fillet_out](pictures/filletbox_fillet_out.png)

### Sec 27.96: the Intersection join at a tangent edge

The Intersection join builds no tubes: the neighbour's offset ran round its
cylinder into a lip (outward) or never met the removed face (inward,
unhollowed; a fillet threw). The gap is now closed by the tube's sharp
counterpart, the user's choice of a square corner: a strip of the
neighbour's tangent plane a thickness into the removed face, offset with it,
and a wall square to the removed face at its far edge, cubes at the corners.
Upstream fails all eight.

![filletbox_end_join_in](pictures/filletbox_end_join_in.png)
![filletbox_end_join_out](pictures/filletbox_end_join_out.png)
![filletbox_side_join_in](pictures/filletbox_side_join_in.png)
![filletbox_side_join_out](pictures/filletbox_side_join_out.png)
![filletbox_fillet_join_in](pictures/filletbox_fillet_join_in.png)
![filletbox_fillet_join_out](pictures/filletbox_fillet_join_out.png)
![filletbox25_fillet_join_in](pictures/filletbox25_fillet_join_in.png)
![filletbox25_fillet_join_out](pictures/filletbox25_fillet_join_out.png)

### Sec 27.100: a cavity sealed below the removed face

A cone with a through hole, its top removed, inward: the wall is 1.4 thick
at the top, thinner than twice the thickness, and the cone's and the hole's
inner offsets cross at z=6.485. No cavity reaches the removed face. The
material is every point within the thickness of a face that stays, so the
cavity -- farther than that from all of them -- is closed, and the top stays
as skin over it: the user's choice of answer, two shells, the input's skin
and a void of 112.92. The fork gave a "valid" solid larger than its input
with intersection off (a shell crossing itself) and invalid ones with it on.

Three causes. With intersection on, the loop on the cone's periodic offset
let the band from the top circle down to the crossing, found first, take
the crossing circle from the band below it -- the cavity -- which was never
built: a band beyond the seam's span now gives way to one on it, and the
band is made straight from the seam's two pcurves (`ShapeFix_Wire` turned
the seam's pcurves round for one band and left the next running the wrong
way). With every offset face then in a closed shell, the thick solid is the
skin, the removed face kept, with those shells as voids. With intersection
off the two offsets were never intersected -- they are not neighbours; walls
beside a removed face too thin for it are now built with intersection on.
That also fixed the cone's bottom removed with intersection off, whose
cavity ran on past the crossing, "valid" and 0.761 too much (upstream says
invalid). Upstream's intersection mode gives the sealed void; its default
mode is invalid.

![conehole_top_in](pictures/conehole_top_in.png)
![conehole_top_inter_in](pictures/conehole_top_inter_in.png)
![conehole_bottom_in](pictures/conehole_bottom_in.png)

### Sec 27.101: faces left in pieces, and down to one face of a short shape

A cylinder with its side removed keeps its two caps, and nothing joins them:
the faces that stay fall apart into pieces sharing no edge. OCCT refuses
that before it starts (`BRepOffset_NotConnectedShell`), upstream too -- the
suite had it as `cyl_side_out/in`, XFAIL, put down to the seam, which has
nothing to do with it. Each piece is now a thick solid of its own, the shape
with the other pieces removed as well, and the result is their union: a
compound of solids (the user's choice), fused where they overlap, and one
solid where one is left. A cylinder 1.5 high, inward by 1, gives two discs
that overlap, and their union is the whole cylinder. A piece that does not
come out as one valid closed shell refuses the whole as before: a pocket's
walls and floor, the outside of the shape removed around them, are built
wrong still, and their union would pass for an answer.

The short cylinder showed a second fault, in the loop. Keep only the bottom
of a box 1.5 high, inward by 1: each removed side is split where the
bottom's offset crosses it, and the loop kept the piece of the side's edge
nearer one of its ends -- above the offset, since the band below it is the
longer one. It now keeps the piece that runs from the end on an edge of a
face that stays. The fork had given the box with the plate's complement as
a void (sec 27.100's sealed cavity took that shell for one), and before that
an empty shell; upstream is right.

![cyl_side_out](pictures/cyl_side_out.png)
![cyl_side_in](pictures/cyl_side_in.png)
![short_cyl_side_in](pictures/short_cyl_side_in.png)
![short_box_bottom_in](pictures/short_box_bottom_in.png)

### Sec 27.102: the Intersection join with one face left

A cylinder down to its top, a box down to its bottom, with the Intersection
join: the fork threw (`BRepAlgo_Image::Bind`) where upstream answers. Its
face list for the intersections carries the removed faces, and with one face
left the removed faces' enlarged copies had been split; binding the splits
as offsets met the image a removed face already has. A removed face now goes
the way it goes when not split. Underneath was a second fault: the removed
side's seam reached the side's loop once, not once each way -- the
intersection took the side's edges from a map, which holds a seam once --
and the loop, with no band to build, kept the whole side: the cap outward
came out as the cylinder. The seam is recorded both ways now. Upstream's cap
outward is inside out. The cylinder's side, in pieces (sec 27.101), is
answered with this join too.

![cyl_cap_alone_join_out](pictures/cyl_cap_alone_join_out.png)
![box_bottom_alone_join_in](pictures/box_bottom_alone_join_in.png)
![cyl_side_join_out](pictures/cyl_side_join_out.png)

### Sec 27.103: pockets left in pieces

A box with a blind hole, its top removed: the hole's wall and floor are a
piece of their own beside the outside's (sec 27.101), and the fork built it
wrong, so the whole was refused. Upstream refuses it as well.

Two causes. A piece was the shape with the other pieces removed as well, and
the removed faces that do not touch it came back as a shell of their own
(upstream does the same); a piece is now an open shell, its faces and the
removed faces beside them. And on that open shell the fork was wrong where
upstream is right: an edge of a removed face that no other face shares -- the
top's square rim here -- was stretched and put in the top's loop, the loop
closed it round the hole's rim and the offset rim, and the wires nested a
step off, so the ring between the two rims was never built. Free edges stay
out of the loop now; upstream drops them only because they happen not to
chain. Every pocket and boss of the sweep is answered, in every mode, and
each was worked out by hand: here, outward, the outside with arcs (300 +
15 pi + 2 pi / 3) and the hole's pot (10 pi); inward the pot sits on the
floor's plate and they fuse. A square pocket's piece inward floats in the
cavity, apart from the walls.

![blind_top_out](pictures/blind_top_out.png)
![blind_top_in](pictures/blind_top_in.png)
![pocket_top_in](pictures/pocket_top_in.png)
![cylpocket_top_join_out](pictures/cylpocket_top_join_out.png)
![blind_top_in_inter_join](pictures/blind_top_in_inter_join.png)
![blind_wall_out_join](pictures/blind_wall_out_join.png)
![pocket_top_out_inter_join](pictures/pocket_top_out_inter_join.png)
![cylboss_shoulder_in](pictures/cylboss_shoulder_in.png)

The free rim alone: the hole's wall and floor and the top beside them as an
open shell, the top removed. Upstream is right; the fork gave the pot
without its ring outward and the hole itself inward.

![blind_pot_shell_out](pictures/blind_pot_shell_out.png)
![blind_pot_shell_in](pictures/blind_pot_shell_in.png)

### Sec 27.104: the input left as it was

The result was right; the input was not. A PartDesign Pad of an arc whose
parameters run past 2 pi, under a Thickness removing the arc's face and the
caps, came back inside out on a box that does not freeze shape values
(FreeCAD's copy-on-write protects a frozen one): the thickness had edited
its input. The loop's pruning builds a test face of each wire on the face's
own surface from the loop's own edges, some of them the input's, and ran
ShapeFix on it when it was invalid; ShapeFix shifted a shared edge's pcurve
by a period, on that surface -- the input face's own pcurve. It works on a
copy now. Here a ring sector straddling angle 0, padded 27, its outer arc's
face and caps removed; the panels show the input after the call. Upstream
throws, and leaves the input alone. The suite now runs unfrozen and checks
every input.

![sector_outer_arc_input](pictures/sector_outer_arc_input.png)

### Sec 27.105: every face removed

A sphere has one face. Remove it and no face stays to be thickened: there is
no answer, and the fork -- like upstream -- gave one, the sphere itself,
"valid" and unhollowed. A box with all six faces removed came back the same
way. The torus, one face too, was refused already. A thickness with every
face removed is refused now (the user's choice); the panels show the input
where it is.

![sphere_face_refused_in](pictures/sphere_face_refused_in.png)
![box_all_faces_refused_in](pictures/box_all_faces_refused_in.png)

## The captured models

The suite's four document cases (realthunder/OCCT#1-#4: an elliptic pad, a
loft, a pad, a revolution) are older fixes, made to the chain on 7.7.2 and
ported. Upstream 8.0.1 gets the three thickness models right too, so they
have no before-and-after here.

## Making the pictures

`tests/thickness/pictures/make_pictures.sh` does it all; about
ten minutes on the dev box, most of it building libraries:

1. `mkold.sh` builds scratch `TKBool`/`TKOffset` libraries -- upstream's
   eleven chain files at `91be8c4c71`, and the fork's sources at each
   stage's "before" commit (`STAGES` in `cases.py`) -- from the build tree's
   own compile commands (`oldbuild.py`).
2. `compute.py` runs every case of `cases.py` under `FreeCADCmd`, once per
   library set (preloaded with `LD_PRELOAD`; the FreeCAD modules carry an
   RPATH, so `LD_LIBRARY_PATH` does not reach them) and once on the fork as
   built, keeping each result as a `.brep` and a `results.json`.
3. `mkjobs.py` and `render.py` draw the panels in the FreeCAD GUI under
   `xvfb-run`; `compose.py` lays them out (Pillow, from the FreeCAD conda
   env).

A case added to `cases.py` with its stage gets its picture on the next run.
A case listed in `INPUT` is pictured for what the call leaves of its input:
the panels show the input after the thickness, judged against its own volume.
A case in `REFUSED` is right when the call throws; its panels show the input.
The cases run unfrozen, as the suite does.
The renderer turns the transaction log off; the pictures need no history.
(While drawing, the log's worker once crashed writing a shape the viewer
was meshing: FreeCAD's docs/TransactionLog.md sec 27.97, fixed in sec 27.98.)
