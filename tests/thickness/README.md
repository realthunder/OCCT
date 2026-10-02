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

The suite runs with FreeCAD's shape values unfrozen (`ImmutableShapeValues`
off), as any build that does not freeze them runs: a frozen input is
protected by copy-on-write from an algorithm that edits it, and the suite
must see such an edit. `THICK_FREEZE=1` runs it frozen. Every case also
checks that its input comes back as it went in (valid, the same volume).

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

Status on `LinkVibe-801` (2026-10-01):

| Case | Status | Symptom |
|------|--------|---------|
| `cyl_*_out/in` (top, bottom), `hole_*` (all eight) | PASS | |
| `cyl_side_out/in` (remove lateral face) | PASS since sec 27.101 | was `StdFail_NotDone`, upstream too: not the seam but the caps left in two pieces (`BRepOffset_NotConnectedShell`). Now a compound of two discs, see "Faces left in pieces" below |

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

## Concave corners and removed faces, Arc join (2026-09-30, sec 27.90)

| Case | Model | What it covers |
|------|-------|----------------|
| `lbox_arm_end_out`, `tshape_arm_end_out` | L-box, T: the end face of an arm, outward | past a concave corner the arc face along one edge is cut by the arc of the next; the corner piece beyond the cut, closed by the arc's own end, came out as a face of its own (the loop keeps every piece of an edge). A wire through a piece beyond the edge's own span now loses to a rival that stays inside |
| `pocketbox_wall_out` | box with a pocket, a pocket wall, outward | the same; no reference volume (upstream's result is invalid) |
| `tshape_bar_top_out` | T: the bar top beside the post (a concave removed face), outward | the removed face's stretched edge ran on along the back of the T and a piece of it lay on an arc face's tangent line (dropped now), and a corner sphere met before the face that renewed its edge kept the old one (faces are gone over again). Upstream fails it; the reference is the fork's own |

## Inward, intersection mode, the Intersection join (2026-09-30, sec 27.91-27.92)

| Case | Model | What it covers |
|------|-------|----------------|
| `tshape_bar_top_in` | T: the concave bar top, inward | the stretched removed face crossed the far end wall's inner edge above the bar's inner top, and the crossing counted as that edge's own end; the corner above the bar's inner arc came out as a face of its own. Worked out by hand (upstream fails it) |
| `lbox_top_inter_join_*`, `tshape_back_inter_join_*` | L-box top, T back; intersection on, Intersection join | an edge two faces share was trimmed twice (an indexed map's `Add()` is never 0) and the second pass cut the section short. Upstream's volumes |
| `conehole_bottom_inter_in` | cone with a through hole, bottom removed, inward, intersection on | a circle running round the hole's offset beyond the seam's span came out as a face of no area. Worked out by hand, = upstream |
| `tshape_bar_top_right_join_in`, `pocket_floor_join_*` | Intersection join, intersection off | the thick solid came back inside out (the quilt's shells met in an order that turned them); now oriented by classification |
| `blindhole_floor_join_*` | box with a blind hole, the floor removed, Intersection join | the floor meets the wall at a concave edge; the section on the wall's offset was the wrong way round. Worked out by hand |
| `lbox_notch_wall_inter_join_out`, `tshape_post_wall_inter_join_in`, `pocket_wall_inter_join_in` | concave removed faces, intersection on, Intersection join | the splits of the offset faces left the rim uncut; such shapes are built as with intersection off |
| `filletbox_end_*`, `filletbox_side_*` | box with its vertical edges filleted, an end or side face removed (Arc join) | the removed face is tangent to the fillets, whose offsets never meet it (inward) or meet it far round the cylinder (outward, a lip). A tube round the tangent edge closes the gap, turning into the removed face, and outward an eighth of a sphere closes each corner with the top and bottom edges' tubes (sec 27.93). Worked out by hand (Steiner for outward) |
| `filletbox_fillet_*` | the same box, a fillet removed (Arc join) | the removed face is curved: the tubes and corners are built on its cylinder; the loop on that periodic surface kept a dozen wires and now walks the angles, its edges lying within one period; inward a piece of the floor's offset cut off by the cylinder's circle hung on the shell and is dropped (sec 27.95). Worked out by hand |
| `filletbox_*_join_*`, `filletbox25_fillet_join_*` | the same box, an end, side or fillet face removed, Intersection join (a fillet of 2.5 too) | that join builds no tubes: the neighbour's offset ran round its cylinder into a lip (outward) or never met the removed face (inward, unhollowed). The gap is closed with the tube's sharp counterpart -- a strip of the neighbour's tangent plane a thickness into the removed face, offset with it, and a wall square to the removed face at its far edge; cubes at the corners (sec 27.96). Worked out by hand |

## A cavity sealed below the removed face (2026-10-01, sec 27.100)

| Case | Model | What it covers |
|------|-------|----------------|
| `conehole_top_in_sealed_*` | cone with a through hole, top removed, inward, all four modes | the wall is 1.4 thick at the top and the inner offsets cross at z=6.485: no cavity reaches the removed face. The right result (the user's choice) is two shells: the input's skin, the top kept, and a closed void. The loop let the band above the crossing take the crossing circle from the cavity's band below it; the cavity's band is built now, from the seam's pcurves, and a thick solid whose offset faces all close up is the skin and its voids. With intersection off the two offsets were never intersected; walls beside a removed face too thin for it are built with intersection on. Hand values: skin 471.2389, void 112.9243. Upstream's intersection mode gives this; its default mode is invalid |
| `conehole_bottom_in`, `conehole_bottom_join_in` | the same cone, bottom removed, inward, intersection off | the same crossing near the top face, never intersected: the cavity ran on to z=7, a sliver inside out, "valid" and 0.761 too much. Built with intersection on now. Worked out by hand |

The 712-run sweep against upstream after these: the fork right in all 150
runs of every mode (upstream 106-112, counting the sealed void right),
worse than upstream in no run.

## Faces left in pieces, and a short wall (2026-10-01, sec 27.101)

`BRepOffset_MakeOffset::CheckInputData` refuses a shape whose faces that stay
fall apart into pieces sharing no edge once the removed faces are gone
(`BRepOffset_NotConnectedShell`), upstream too. `MakeThickSolidByPieces` now
makes each piece a thick solid of its own -- the shape with the other pieces
removed as well -- and returns their union: a compound of solids, fused where
they overlap; one solid left is returned as it is. The pieces' history is
merged, so `Generated`/`Modified` answer as for one shape (FreeCAD's element
names carry `THK`). A piece that is not one valid closed shell of positive
volume refuses the whole, as before: a pocket's walls and floor, with the
outside of the shape removed around them (a blind hole's top removed), did
not come out right, and their union would pass for an answer (sec 27.103
below).

Making the short cases right needed a fix of its own, in `BRepAlgo_Loop`:
keeping only the bottom of a box shorter than twice the thickness, inward,
the loop split each removed side at the bottom's offset and kept the piece
nearer an end of the side's edge -- above the offset there. Now the piece
from the end on a const edge (the face that stays) is kept when only one end
is; otherwise the nearer end's, as before. The thick solid had come out as
the box with the plate's complement as a void (sec 27.100's sealed cavity)
and, before that, with an empty shell; upstream is right.

| Case | Model | What it covers |
|------|-------|----------------|
| `cyl_side_*` | cylinder, side removed (Arc join, intersection off and on) | two discs of 16 pi, a compound |
| `box_sides_*` | box 10 x 10 x 20, four sides removed | two plates of 100 |
| `ring_walls_in` | cylinder with a hole, both walls removed | two rings of 21 pi |
| `short_cyl_side_in` | cylinder 1.5 high, side removed, inward | the discs overlap and fuse: the whole cylinder, 24 pi |
| `short_box_bottom_in*`, `short_cyl_bottom_in` | box and cylinder 1.5 high, everything but the bottom removed, inward | the loop's piece above the offset; a plate, one shell |

The sweep: the 150 runs of every mode in scope are unchanged. Of the runs out
of scope -- a face whose removal leaves pieces -- the cylinder's, the cone's and
the elliptic pad's side, Arc join, are answered now, by hand: 2 x 16 pi; the
frustum slabs 84.58 + 10.36 = 94.94 outward and 72.80 + 15.07 = 87.87 inward;
2 x 50 pi. The pocket shapes (`boxhole2` Face3 and Face7, `pocket` Face3,
`cylpocket` Face2, `cylboss` Face2) stay refused.

## The Intersection join with one face left (2026-10-01, sec 27.102)

Where one face stays beside the removed ones -- a cylinder down to its top, a
box down to its bottom -- the Intersection join threw `BRepAlgo_Image::Bind`
in `BuildOffsetByInter`: the fork puts the removed faces in its face list
for the intersections, a removed face already has its own image, and its
enlarged face had been split. A removed face now goes the way it goes when
not split -- its edges recorded, nothing bound. That gave the box's plate,
and the cylinder's cap outward as the whole cylinder: `ContextIntByInt`
takes a removed face's edges from a map, which holds a seam once, where
`ContextIntByArc`'s explorer meets it once each way; given once, the loop on
the side built no band and kept all of it. The stretched seam is recorded
both ways now. Upstream answers these (the cap outward inside out).

| Case | Model | What it covers |
|------|-------|----------------|
| `cyl_cap_alone_join_*` | cylinder, side and bottom removed, Intersection join | the throw, and the seam once: a disc of 16 pi |
| `box_bottom_alone_join_in` | box, all but the bottom removed, Intersection join, inward | the throw: a plate of 100 |
| `cyl_side_join_*` | cylinder, side removed, Intersection join | the pieces of sec 27.101 with this join, refused until now |

The sweep: in scope unchanged; the cylinder's, cone's and elliptic pad's
side are answered with the Intersection join too, at the Arc join's volumes
(one face to a piece, no edge between faces that stay for the joins to
differ on). The pocket shapes and the torus's face stay refused, NotDone
now instead of `Bind` or `NoSuchObject` thrown.

## Pockets left in pieces (2026-10-01, sec 27.103)

Remove the face round a pocket's rim and the pocket's walls and floor are a
piece of their own. Two faults kept such pieces refused. A piece was the
shape with the other pieces removed as well, and the removed faces that do
not touch it came back as a shell of their own -- the outside of the box
round the pocket (upstream does the same). A piece is now an open shell: its
faces and the removed faces beside it. Then the removed top's free rim: a
removed face's edges that no kept face shares are stretched and put in its
loop, where an edge between two removed faces belongs; a free one bounds no
material, and the loop closed the stretched square round the hole's rim and
the offset rim, so the wires nested a step off and the ring between them was
never built. Free edges are left out of the loop now (upstream drops them
only because they happen not to chain). And with the Intersection join and
intersection on, the outside's piece threw (`NoSuchObject`) asking the
analysis for the ancestors of the top's free rim, which it does not know.

| Case | Model | What it covers |
|------|-------|----------------|
| `blind_top_*` | box 10 x 10 x 5, a blind hole of radius 2, 3 deep, top removed | outward two solids, 300 + 15 pi + 2 pi / 3 and the pot, 10 pi; inward the pot sits on the floor's plate and they fuse (Arc join, and Intersection join with intersection on). By hand |
| `blind_wall_out_join` | the same, the hole's wall removed | the outside with the hole cut through its top plate, and the floor's plate |
| `blind_pot_shell_*` | the hole's wall and floor and the top as an open shell, the top removed | the free rim alone: the fork gave the pot without its ring outward (invalid) and the hole itself inward; upstream is right |
| `pocket_top_*` | box with a square pocket, top removed | inward the pocket's piece floats in the cavity, two solids |
| `cylpocket_top_out_join`, `cylboss_shoulder_in` | round pocket in a cylinder, boss on one, the face round it removed | 127 pi + 19 pi; 69 pi + 24 pi |

The sweep: in scope unchanged. Every pocket and boss shape is answered now,
in every mode -- `boxhole2` Face3 and Face7, `pocket` Face3, `cylpocket`
Faces 1, 2 and 4, `cylboss` Faces 1, 2 and 4 -- and all 68 runs check by
hand. The torus's face is the one refusal left.

## The input left as it was (2026-10-01, sec 27.104)

A PartDesign Pad of an arc whose parameters run past 2 pi, under a Thickness
removing the arc's face and the caps, came back inside out on the Windows
box (`issue3_pad_thickness`: Pad -189356.2997, occ-issues `local03`): the
thickness had edited the Pad's shape. The loop's pruning builds a test face
of each wire on the face's own surface from the loop's own edges, some of
them the input's, and ran `ShapeFix_Shape` on it when it was invalid;
ShapeFix shifted the pcurve of a shared edge by a period, on that surface --
the input face's own pcurve. Reached since a periodic face whose edges lie
within one period takes the plane path (sec 27.95). ShapeFix now works on a
copy. It showed only unfrozen: FreeCAD freezes the values it finds this fork
for, and this suite ran frozen -- it runs unfrozen now, and every case
checks its input.

| Case | Model | What it covers |
|------|-------|----------------|
| `issue3_pad_thickness` | the captured model, unfrozen | the Pad stays valid under the Thickness |
| `sector_outer_arc_input`, `sector_outer_arc_past_period_input` | a ring sector straddling angle 0 padded 27, its outer arc's face and caps removed | the input as it was (inside out before, frozen or not: shapes made in Python are not frozen) |

The sweep, with an input check added: no run alters its input, before the
fix or after, and every line is unchanged -- its 16 solids never cut a
periodic face across its seam.

## Every face removed (2026-10-02, sec 27.105)

The sweep's last refusal was the torus, its one face removed -- rightly: no
face stays to be thickened. The sphere, a one-face shape too, came back as
the sphere itself, "valid" and unhollowed, and so did a box with all six
faces removed: an answer that is the input. `MakeThickSolid` now refuses
when no face stays (the user's choice).

| Case | Model | What it covers |
|------|-------|----------------|
| `sphere_face_refused_*`, `sphere_face_join_refused_in` | sphere, its face removed | refused; was the sphere back, 523.60 |
| `torus_face_refused_*` | torus, its face removed | refused, as before |
| `box_all_faces_refused_in`, `box_all_faces_inter_refused_out` | box, all six faces removed | refused; was the box back, 600 |

## A face closed at a pole (2026-10-02, sec 27.106)

On a periodic face the loop builds a seam wire from a closed edge, a piece
of the seam and the closed edge at its far end. Where the seam ends at a
pole the closed edge there is the pole's degenerated edge, which the loop
keeps out of its vertex map: no seam wire was built, and the one circle was
left as a wire of its own, bounding nothing -- a face of no area and an
invalid result, since the chain was ported (upstream is right). The band is
now closed by the pole's edge (`BRepAlgo_Loop::FindLoop`). The rim is closed
in the removed face's plane, so the volumes are differences of sphere caps,
pi h^2 (3R - h) / 3; the cone's is worked out in FreeCAD's
docs/TransactionLog.md sec 27.106.

| Case | Model | What it covers |
|------|-------|----------------|
| `dome_flat_out/in`, `dome_flat_join_out`, `dome_flat_inter_in` | half a sphere of radius 5, flat face removed, 0.5 | 86.6556 out, 70.9476 in; was invalid, 216.5455 and 231.5055 |
| `bowl_flat_out` | the lower half | the pole below: 86.6556 |
| `cap_flat_in` | the sphere above latitude 30 | 33.6412; was invalid, 71.3403 |
| `segment_flat_out` | the sphere below latitude 30 | 127.8890; was invalid, 387.9306 |
| `cone_base_out`, `cone_up_base_out` | cone with its apex, radius 4, height 6, base removed | 52.4095, the apex rounded; was invalid, 62.5293 |

Still wrong, upstream too: the cone inward with the Arc join throws (apex
down `BRep_Builder::Infinite parameter` from `BRepOffset_Inter2d`, apex up
`NCollection_DataMap::Find` from `BRepOffset_Tool::ExtentFace`); apex down
with the Intersection join the cone comes back unhollowed (100.5310); half a
dome (180 degrees) fails in most modes. The sweep has no such solid.

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

## Pictures

`models/pictures/<case>.png` shows, for each case a fix turned from failing
to passing, three results side by side -- upstream's chain files at
`91be8c4c71`, the fork just before the fix, the fork now -- each whole and
cut open. `models/Thickness.md` walks through them fix by fix and says
how to read them. `pictures/make_pictures.sh` makes them again (scratch
libraries for upstream and for each stage, the cases run on each, rendered
under Xvfb); add a case to `pictures/cases.py`, with its stage, to picture it.

## Layout

```
tests/thickness/
  README.md          this file
  run_tests.py       the suite (FreeCADCmd script)
  models/            captured user models, one per reported issue
    Thickness.md     the fixes, before and after, in pictures
    pictures/        before and after, one PNG per fixed case
  pictures/          the tools that make them (make_pictures.sh)
```
