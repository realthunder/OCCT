# Regression tests for the MakeThickSolid / offset fix chain.
#
# Run with any FreeCAD build linked against this OCCT:
#
#     FreeCADCmd tests/thickness/run_tests.py
#
# Exit code is 0 when no expected-pass case fails.  Known-broken cases are
# declared XFAIL below; when one starts passing the suite prints
# UNEXPECTED-PASS so its expectation (and README.md) can be updated.
#
# See README.md for what each case covers and the history of each problem.

import math
import os
import sys
import traceback

import FreeCAD as App
import Part

HERE = os.path.dirname(os.path.abspath(__file__))


def emit(line):
    # FreeCAD redirects sys.stdout into its own console, which swallows
    # script output under FreeCADCmd - write straight to fd 1.
    os.write(1, (line + "\n").encode())
MODELS = os.path.join(HERE, "models")
RELTOL = 1e-3

# FreeCAD freezes shape values (Immutable TShapes, copy-on-write) when it finds
# this fork; a frozen input is protected from an algorithm that edits it. Run
# unfrozen by default, the way any other build of FreeCAD runs it, so an edit
# of the input shows (THICK_FREEZE=1 for the frozen path).
App.ParamGet("User parameter:BaseApp/Preferences/Mod/Part").SetBool(
    "ImmutableShapeValues", os.environ.get("THICK_FREEZE", "0") != "0")


def signature(shape):
    """What a thickness must leave of its input: the shape valid, the volume,
    and no tolerance grown -- a vertex of a moved cylinder came back with a
    tolerance as large as the move, valid and of the same volume, and the
    next thickness of that shape was wrong."""
    tol = max([0.0] + [x.Tolerance for x in shape.Vertexes + shape.Edges + shape.Faces])
    return (shape.isValid(), round(shape.Volume, 6), round(tol, 9))


def input_problems(shape, before):
    after = signature(shape)
    if after == before:
        return []
    return ["input changed: valid=%s vol=%.4f tol=%g" % after]

results = []  # (name, verdict, detail); verdict in PASS/FAIL/XFAIL/UNEXPECTED-PASS


def report(name, ok, expect_fail, detail):
    if ok and not expect_fail:
        verdict = "PASS"
    elif ok and expect_fail:
        verdict = "UNEXPECTED-PASS"
    elif not ok and expect_fail:
        verdict = "XFAIL"
    else:
        verdict = "FAIL"
    results.append((name, verdict, detail))
    emit("%-32s %-16s %s" % (name, verdict, detail))


# ---------------------------------------------------------------------------
# Document cases: recompute the captured model, everything must stay valid.
# ---------------------------------------------------------------------------

def document_case(name, filename, volumes=None):
    """Open models/<filename>, force a full recompute and require every
    feature to recompute without error into a valid shape.  `volumes` maps
    object Name -> expected volume of its shape."""
    try:
        doc = App.openDocument(os.path.join(MODELS, filename))
        try:
            for obj in doc.Objects:
                obj.touch()
            doc.recompute()
            bad = []
            for o in doc.Objects:
                states = [str(s) for s in o.State]
                if "Invalid" in states or "Error" in states:
                    bad.append("%s:%s" % (o.Name, states))
                if hasattr(o, "Shape") and not o.Shape.isNull():
                    try:
                        if not o.Shape.isValid():
                            bad.append("%s:invalid-shape" % o.Name)
                    except Exception:
                        bad.append("%s:check-threw" % o.Name)
            if not bad and volumes:
                for oname, vol in volumes.items():
                    got = getattr(doc, oname).Shape.Volume
                    if abs(got - vol) > RELTOL * abs(vol):
                        bad.append("%s:volume %.4f != %.4f" % (oname, got, vol))
            report(name, not bad, False, "; ".join(bad) if bad else "ok")
        finally:
            App.closeDocument(doc.Name)
    except Exception:
        report(name, False, False, traceback.format_exc().splitlines()[-1])


def placed(shape):
    """The shape under a location: an object's placement."""
    moved = shape.copy()
    moved.Placement = App.Placement(App.Vector(3, 4, 5), App.Rotation(App.Vector(1, 2, 3), 40))
    return moved


# ---------------------------------------------------------------------------
# Programmatic cases: thickness of simple solids, removing one face at a time.
# ---------------------------------------------------------------------------

def thickness_case(name, shape, face_index, value, expect, ref_volume=None, inter=False, join=0):
    """makeThickness(mode=Skin, join=Arc unless given) removing 1-based face `face_index`
    (or the faces of a list of them).

    expect='pass':  result must be a valid solid, one closed shell, and match
                    ref_volume (captured from a verified-good run).
    expect='xfail': known broken (see README.md); any exception, invalid
                    shape, or open shell counts as the expected failure.
    """
    detail = ""
    try:
        indices = face_index if isinstance(face_index, list) else [face_index]
        faces = [shape.Faces[i - 1] for i in indices]
        before = signature(shape)
        r = shape.makeThickness(faces, value, 1e-7, inter, False, 0, join)
        problems = input_problems(shape, before)
        if r.ShapeType != "Solid":
            problems.append("type=%s" % r.ShapeType)
        if not r.isValid():
            problems.append("invalid")
        if len(r.Shells) != 1:
            problems.append("shells=%d" % len(r.Shells))
        elif not r.Shells[0].isClosed():
            problems.append("open-shell")
        if not problems and ref_volume is not None:
            if abs(r.Volume - ref_volume) > RELTOL * abs(ref_volume):
                problems.append("volume %.4f != %.4f" % (r.Volume, ref_volume))
        ok = not problems
        detail = "; ".join(problems) if problems else "vol=%.4f" % r.Volume
    except Exception as e:
        ok = False
        detail = "EXCEPTION " + str(e).strip().splitlines()[-1]
    report(name, ok, expect == "xfail", detail)


def pieces_case(name, shape, face_indices, value, solids, inter=False, join=0):
    """makeThickness(mode=Skin, join=Arc unless given) removing the 1-based faces in
    `face_indices`. `solids` lists the volume of each solid expected, in any
    order: one is a Solid, several a Compound of them. Each must be valid
    with one closed shell -- a void is a wrong answer here."""
    detail = ""
    try:
        before = signature(shape)
        r = shape.makeThickness([shape.Faces[i - 1] for i in face_indices], value, 1e-7,
                                inter, False, 0, join)
        problems = input_problems(shape, before)
        want = "Solid" if len(solids) == 1 else "Compound"
        if r.ShapeType != want:
            problems.append("type=%s" % r.ShapeType)
        if not r.isValid():
            problems.append("invalid")
        if len(r.Solids) != len(solids):
            problems.append("solids=%d" % len(r.Solids))
        for s in r.Solids:
            if len(s.Shells) != 1 or not s.Shells[0].isClosed():
                problems.append("shells=%d" % len(s.Shells))
                break
        if not problems:
            got = sorted(s.Volume for s in r.Solids)
            for g, w in zip(got, sorted(solids)):
                if abs(g - w) > RELTOL * abs(w):
                    problems.append("volumes %s != %s" % (
                        ["%.4f" % v for v in got], ["%.4f" % v for v in sorted(solids)]))
                    break
        ok = not problems
        detail = "; ".join(problems) if problems else "vols=%s" % ", ".join(
            "%.4f" % s.Volume for s in r.Solids)
    except Exception as e:
        ok = False
        detail = "EXCEPTION " + str(e).strip().splitlines()[-1]
    report(name, ok, False, detail)


document_case("issue1_ellipse_thickness", "issue1_ellipse_thickness.FCStd",
              volumes={"Thickness": 3598.2930})
document_case("issue2_broken_loft", "issue2_broken_loft.FCStd",
              volumes={"AdditiveLoft": 6143.8106})
document_case("issue3_pad_thickness", "issue3_pad_thickness.FCStd",
              volumes={"Thickness": 1241.0718})
document_case("issue4_revolution_thickness", "issue4_revolution_thickness.FCStd",
              volumes={"Thickness": 431.4454})

# Plain cylinder: Face1 = lateral (seam), Face2 = bottom, Face3 = top.
# Removing the lateral face leaves the two caps, which share no edge: each is
# a thick solid of its own, a disc of radius 4 and the thickness, 16 pi.
cyl = Part.makeCylinder(4, 20)
pieces_case("cyl_side_out",       cyl, [1], +1.0, [16 * math.pi] * 2)
pieces_case("cyl_side_in",        cyl, [1], -1.0, [16 * math.pi] * 2)
pieces_case("cyl_side_in_inter",  cyl, [1], -1.0, [16 * math.pi] * 2, True)
thickness_case("cyl_bottom_out",  cyl, 2, +1.0, "pass", 637.5858)
thickness_case("cyl_bottom_in",   cyl, 2, -1.0, "pass", 468.0973)
thickness_case("cyl_top_out",     cyl, 3, +1.0, "pass", 637.5858)
thickness_case("cyl_top_in",      cyl, 3, -1.0, "pass", 468.0973)

# The faces that stay fall apart into pieces (README.md, "Pieces"). Each
# piece is a thick solid of its own; pieces that overlap are fused. Box Faces
# 1-4 are its sides, Face5 the bottom, Face6 the top.
box = Part.makeBox(10, 10, 20)
pieces_case("box_sides_out",      box, [1, 2, 3, 4], +1.0, [100.0, 100.0])
pieces_case("box_sides_in",       box, [1, 2, 3, 4], -1.0, [100.0, 100.0])
pieces_case("box_sides_in_inter", box, [1, 2, 3, 4], -1.0, [100.0, 100.0], True)
ring = Part.makeCylinder(5, 10).cut(Part.makeCylinder(2, 10))
pieces_case("ring_walls_in",      ring, [1, 4], -1.0, [21 * math.pi] * 2)
# Shorter than twice the thickness: inward, the two discs overlap, and their
# union is the whole cylinder.
pieces_case("short_cyl_side_in",  Part.makeCylinder(4, 1.5), [1], -1.0, [24 * math.pi])

# Only the bottom of a box shorter than twice the thickness stays: a plate of
# the thickness. The loop split each removed side at the bottom's offset and
# kept the piece nearer an end of the side's edge, the one beyond the offset
# here; the thick solid came out as the box with the plate's complement as a
# void, or, before the sealed cavity, with an empty shell (upstream is right).
pieces_case("short_box_bottom_in", Part.makeBox(10, 10, 1.5), [1, 2, 3, 4, 6], -1.0, [100.0])
pieces_case("short_box_bottom_in_inter", Part.makeBox(10, 10, 1.5), [1, 2, 3, 4, 6], -1.0,
            [100.0], True)
pieces_case("short_cyl_bottom_in", Part.makeCylinder(4, 1.5), [1, 2], -1.0, [16 * math.pi])

# The Intersection join with one face left (sec 12). The removed faces'
# enlarged faces were split, and binding the splits as offsets threw
# (BRepAlgo_Image::Bind); and the removed side's stretched seam was recorded
# once, so its loop built no band: the cap outward came out as the whole
# cylinder. Upstream answers these, the cap outward inside out.
pieces_case("cyl_cap_alone_join_out", cyl, [1, 3], +1.0, [16 * math.pi], False, 2)
pieces_case("cyl_cap_alone_join_in",  cyl, [1, 3], -1.0, [16 * math.pi], False, 2)
pieces_case("box_bottom_alone_join_in", Part.makeBox(10, 10, 6), [1, 2, 3, 4, 6], -1.0,
            [100.0], False, 2)
pieces_case("cyl_side_join_out",      cyl, [1], +1.0, [16 * math.pi] * 2, False, 2)
pieces_case("cyl_side_join_in_inter", cyl, [1], -1.0, [16 * math.pi] * 2, True, 2)

# A pocket's walls and floor are a piece of their own when the face round its
# rim is removed (sec 13). A piece is now an open shell, its faces and the
# removed faces beside it: the shape with the other pieces removed brought the
# rest of the removed faces back as a shell of their own (upstream too). And
# the top's loop closed its free rim, stretched, into a wire round the rest,
# so the piece between the hole's rim and the offset rim was never built.
# Each value is worked out by hand: box 10 x 10 x 5 with a blind hole of
# radius 2, 3 deep (Face3 the top, Face7 the hole's wall); outward, the outside
# with arcs (300 + 15 pi + 2 pi / 3) and the hole's pot (10 pi); inward, the
# pot sits on the floor's plate and they fuse.
blind = Part.makeBox(10, 10, 5).cut(Part.makeCylinder(2, 3, App.Vector(5, 5, 2)))
pieces_case("blind_top_out", blind, [3], +1.0,
            [300 + 15 * math.pi + 2 * math.pi / 3, 10 * math.pi])
pieces_case("blind_top_in", blind, [3], -1.0,
            [244 + 19 * math.pi + math.pi ** 2 / 2 * (2 + 4 / (3 * math.pi))])
pieces_case("blind_top_in_inter_join", blind, [3], -1.0, [244 + 24 * math.pi], True, 2)
pieces_case("blind_wall_out_join", blind, [7], +1.0, [508 - 4 * math.pi, 4 * math.pi], False, 2)
# The hole's wall and floor and the top beside them as an open shell, the top
# removed: its square rim is free. The fork gave the pot without its rim's
# ring outward (invalid, two shells) and the hole itself inward; upstream is
# right.
pot = Part.Shell([blind.Faces[i - 1] for i in (7, 8, 3)])
pieces_case("blind_pot_shell_out", pot, [3], +1.0, [10 * math.pi])
pieces_case("blind_pot_shell_in", pot, [3], -1.0,
            [19 * math.pi + math.pi ** 2 / 2 * (2 + 4 / (3 * math.pi))])
# Box 10 x 10 x 6 with a square pocket 5 x 5, 3 deep: inward the pocket's
# piece floats in the cavity, apart from the walls.
pocket = Part.makeBox(10, 10, 6).cut(Part.makeBox(5, 5, 3, App.Vector(2.5, 2.5, 3)))
pieces_case("pocket_top_in", pocket, [3], -1.0, [280, 85 + 8 * math.pi + 2 * math.pi / 3])
pieces_case("pocket_top_out_inter_join", pocket, [3], +1.0, [408, 57], True, 2)
# Round pocket in a cylinder 6 x 6, and a boss on one, the shoulder removed.
cylpocket_top = Part.makeCylinder(6, 6).cut(Part.makeCylinder(3, 3, App.Vector(0, 0, 3)))
pieces_case("cylpocket_top_out_join", cylpocket_top, [2], +1.0,
            [127 * math.pi, 19 * math.pi], False, 2)
cylboss = Part.makeCylinder(6, 4).fuse(Part.makeCylinder(3, 4, App.Vector(0, 0, 4)))
pieces_case("cylboss_shoulder_in", cylboss, [2], -1.0, [69 * math.pi, 24 * math.pi])

# Cylinder with a centered hole: Face1 = outer lateral, Face2 = bottom,
# Face3 = top, Face4 = hole lateral.
ann = Part.makeCylinder(5, 10).cut(Part.makeCylinder(2, 10))
thickness_case("hole_outer_out",  ann, 1, +1.0, "pass", 241.7451)
thickness_case("hole_outer_in",   ann, 1, -1.0, "pass", 257.6106)
thickness_case("hole_bottom_out", ann, 2, +1.0, "pass", 540.3400)
thickness_case("hole_bottom_in",  ann, 2, -1.0, "pass", 461.8141)
thickness_case("hole_top_out",    ann, 3, +1.0, "pass", 540.3400)
thickness_case("hole_top_in",     ann, 3, -1.0, "pass", 461.8141)
thickness_case("hole_inner_out",  ann, 4, +1.0, "pass", 531.0589)
thickness_case("hole_inner_in",   ann, 4, -1.0, "pass", 358.1416)

# Elliptic pad (a single closed ellipse edge, extruded): Face2 = bottom,
# Face3 = top. The offset of the ellipse is a closed B-spline that is not
# periodic, which the context extension used to stretch 100 lengths past its
# ends: the top came back unhollowed, the bottom as two shells.
ell = Part.Face(Part.Wire(Part.Ellipse(App.Vector(0, 0, 0), 10, 5).toShape())).extrude(
    App.Vector(0, 0, 8))
thickness_case("ellipse_bottom_out", ell, 2, +1.0, "pass", 608.6622)
thickness_case("ellipse_bottom_in",  ell, 2, -1.0, "pass", 473.2070)
thickness_case("ellipse_top_out",    ell, 3, +1.0, "pass", 608.6631)
thickness_case("ellipse_top_in",     ell, 3, -1.0, "pass", 473.2073)

# Cylinder with a blind pocket in its top: Face3 = bottom. The loop on the
# offset pocket wall lost the orientation of a circle and closed its seam
# wire with two circles running the same way.
cylpocket = Part.makeCylinder(6, 6).cut(Part.makeCylinder(3, 3, App.Vector(0, 0, 3)))
thickness_case("pocket_bottom_out", cylpocket, 3, +1.0, "pass", 433.9707)
thickness_case("pocket_bottom_in",  cylpocket, 3, -1.0, "pass", 346.7660)

# Box with a through hole: Face7 = the hole. The removed cylinder leaves a
# wall at each end, two seam wires on one face, where only one was allowed.
boxhole = Part.makeBox(10, 10, 5).cut(Part.makeCylinder(2, 5, App.Vector(5, 5, 0)))
thickness_case("boxhole_hole_out", boxhole, 7, +1.0, "pass", 457.5959)
thickness_case("boxhole_hole_in",  boxhole, 7, -1.0, "pass", 282.8673)


# Outward with the Arc join past a concave corner: the arc face along one
# edge is cut by the arc of the next, and the loop kept the corner piece
# beyond the cut as a face of its own -- the loop keeps every piece of an
# edge, and the arc's own end closed the corner. L-box Face4 and T Face1 are
# end faces of an arm; the pocketed box's Face7 is a wall of the pocket
# (no reference volume: upstream's result is invalid).
lbox = Part.makeBox(10, 10, 5).cut(Part.makeBox(5, 5, 5, App.Vector(5, 5, 0)))
thickness_case("lbox_arm_end_out", lbox, 4, +1.0, "pass", 388.5671)
tshape = Part.makeBox(12, 4, 4).fuse(Part.makeBox(4, 4, 10, App.Vector(4, 0, 0))).removeSplitter()
thickness_case("tshape_arm_end_out", tshape, 1, +1.0, "pass", 372.9204)
pocketbox = Part.makeBox(10, 10, 6).cut(Part.makeBox(6, 6, 3, App.Vector(2, 2, 3)))
thickness_case("pocketbox_wall_out", pocketbox, 7, +1.0, "pass", None)

# A concave removed face -- the T's bar top beside the post (Face2) --
# outward with the Arc join. The removed face's stretched edge ran on along
# the back of the T and one of its pieces lay on the tangent line of an arc
# face: the plane got the line twice. And a corner sphere, met before the arc
# face that renewed its edge, kept the old edge. Upstream fails this one too;
# the reference is the fork's own.
thickness_case("tshape_bar_top_out", tshape, 2, +1.0, "pass", 382.4425)

# The same face inward. The stretched edge of the removed face crossed the
# far end wall's inner edge above the bar's inner top, and that crossing
# counted as the end wall edge's own end: the corner above the bar's inner
# arc came out as a face of its own. The reference is worked out by hand:
# 288 less the cavity -- bar 40, post 24, under the removed top 8, and the
# corner outside the concave arc, 2 * (1 - pi/4).
thickness_case("tshape_bar_top_in", tshape, 2, -1.0, "pass", 215.5708)

# A box with a pocket, thickened inward with intersection on and the
# Intersection join, crashed: splitting the trimmed faces had no map from
# trimmed to infinite edges and dereferenced it. The result is still invalid
# (as upstream's is); the case is here to run to the end. Since sec 5
# a removed face at a concave edge is built as with intersection off, so
# this input no longer reaches the guarded code; it throws (its walls are
# twice the thickness, the offsets coincide).
def nocrash_case(name, shape, face_index, value, inter, join):
    try:
        r = shape.makeThickness([shape.Faces[face_index - 1]], value, 1e-7, inter, False, 0,
                                join)
        detail = "valid=%s vol=%.4f" % (r.isValid(), r.Volume)
    except Exception as e:
        detail = "EXCEPTION " + type(e).__name__
    report(name, True, False, "no crash; " + detail)


# Intersection on, the Intersection join, a face of the L-box's top or the T's
# back removed. The section of the removed face with an offset face is
# trimmed to its outermost crossings, then trimmed again from the other face
# that shares it: the guard meant to trim it once, an indexed map's Add(),
# returns the key's index and is never 0. The second pass took the edge's
# own end for a crossing and cut the section short -- the rim face had an
# edge its wall did not (models/Thickness.md sec 4).
# Upstream's volumes.
thickness_case("lbox_top_inter_join_out", lbox, 3, +1.0, "pass", 339.0, True, 2)
thickness_case("lbox_top_inter_join_in",  lbox, 3, -1.0, "pass", 219.0, True, 2)
thickness_case("tshape_back_inter_join_out", tshape, 8, +1.0, "pass", 312.0, True, 2)
thickness_case("tshape_back_inter_join_in",  tshape, 8, -1.0, "pass", 192.0, True, 2)

# A cone with a through hole, its bottom removed, inward with intersection
# on. The cone's and the hole's offsets meet at z=6.485, below the top face's
# offset at z=7, which vanishes; each cut the other at z=7 too, and that
# circle, beyond the seam's span, came out as a face of its own with no area
# (models/Thickness.md sec 4). The reference is worked out by
# hand -- the cavity is the band between r=2.5 and the cone's offset up to
# where they meet -- and is upstream's.
conehole = Part.makeCone(6, 3, 8).cut(Part.makeCylinder(1.5, 8))
thickness_case("conehole_bottom_inter_in", conehole, 3, -1.0, "pass", 307.1946, True, 0)


def sealed_case(name, shape, face_index, value, skin_volume, void_volume, inter, join):
    """A thick solid whose cavity reaches no removed face: a valid solid of two
    closed shells, the skin the input with the removed face kept, and a void."""
    try:
        before = signature(shape)
        r = shape.makeThickness([shape.Faces[face_index - 1]], value, 1e-7, inter, False, 0, join)
        problems = input_problems(shape, before)
        if r.ShapeType != "Solid":
            problems.append("type=%s" % r.ShapeType)
        if not r.isValid():
            problems.append("invalid")
        vols = sorted((abs(Part.Solid(sh).Volume) for sh in r.Shells), reverse=True)
        if len(r.Shells) != 2 or not all(sh.isClosed() for sh in r.Shells):
            problems.append("shells=%s" % ["%.4f" % v for v in vols])
        else:
            for got, want in zip(vols, (skin_volume, void_volume)):
                if abs(got - want) > RELTOL * want:
                    problems.append("shell %.4f != %.4f" % (got, want))
            want = skin_volume - void_volume
            if abs(r.Volume - want) > RELTOL * want:
                problems.append("volume %.4f != %.4f" % (r.Volume, want))
        detail = "; ".join(problems) if problems else "vol=%.4f" % r.Volume
        report(name, not problems, False, detail)
    except Exception as e:
        report(name, False, False, "EXCEPTION " + str(e).strip().splitlines()[-1])


# The same cone, its top removed, inward: the wall is 1.4 thick at the top,
# the cone's and the hole's inner offsets cross at z=6.485, and no cavity
# reaches the removed face. The material is every point within the thickness
# of a face that stays, so the cavity -- r > 2.5, inside the cone's offset,
# z > 1 -- is closed, and the top stays as skin over it (the user's choice,
# models/Thickness.md sec 10). The fork gave a "valid" 577.918
# with intersection off, a shell crossing itself, and invalid shapes with it
# on: the loop let the band above the crossing take the crossing circle from
# the cavity's band, and intersection off never meets the two offsets.
# Upstream's intersection mode gives this result; its default mode is invalid.
# Hand values: skin 471.2389, void 112.9243 (the band between r=2.5 and the
# cone's offset, z from 1 to 6.4853).
for inter in (False, True):
    for join in (0, 2):
        sealed_case("conehole_top_in_sealed_%s_j%d" % ("inter" if inter else "nointer", join),
                    conehole, 2, -1.0, 471.2389, 112.9243, inter, join)

# Its bottom removed with intersection off: the two offsets were never
# intersected and the cavity ran on to z=7, a sliver of it inside out --
# "valid", 0.761 too much (upstream: invalid). The walls beside the removed
# face are too thin for intersection off, which builds it with intersection
# on now (sec 10).
thickness_case("conehole_bottom_in", conehole, 3, -1.0, "pass", 307.1946, False, 0)
thickness_case("conehole_bottom_join_in", conehole, 3, -1.0, "pass", 307.1946, False, 2)

# The Intersection join, intersection off: the T's right bar top inward, and
# a box pocketed 5 x 5 with its pocket floor removed, came back as valid
# solids of negative volume -- inside out, the quilt's shells met in an order
# that flipped them (the T's left bar top, its mirror image, was right). The
# thick solid is now oriented by classification (FreeCAD
# models/Thickness.md sec 5). Worked out by hand: 288 - 72, and the
# box less its pocket grown by the floor's removal.
pocket5 = Part.makeBox(10, 10, 6).cut(Part.makeBox(5, 5, 3, App.Vector(2.5, 2.5, 3)))
thickness_case("tshape_bar_top_right_join_in", tshape, 7, -1.0, "pass", 216.0, False, 2)
thickness_case("pocket_floor_join_out", pocket5, 11, +1.0, "pass", 591.0, False, 2)
thickness_case("pocket_floor_join_in",  pocket5, 11, -1.0, "pass", 367.0, False, 2)

# The Intersection join, a blind hole's floor removed. The floor meets the
# hole's wall at a concave edge, and the section of the floor with the wall's
# offset was oriented as for a convex one: the band on the offset cylinder
# closed on two circles running the same way -- a face of negative area
# (models/Thickness.md sec 5). Worked out by hand: inward, 462.3
# less the box shrunk by 1 around the hole's offset; outward, the box grown by
# 1 less the hole's offset below its top, less the input.
blindhole = Part.makeBox(10, 10, 5).cut(Part.makeCylinder(2, 3, App.Vector(5, 5, 2)))
thickness_case("blindhole_floor_join_out", blindhole, 8, +1.0, "pass", 533.1327, False, 2)
thickness_case("blindhole_floor_join_in",  blindhole, 8, -1.0, "pass", 326.8496, False, 2)

# Intersection on and the Intersection join, a removed face meeting a
# neighbour at a concave edge -- the L-box's notch wall, the T's post wall, a
# pocket's wall. The splits of the offset faces did not cut the neighbour's
# section where the rim needs it; such a shape is now built as with
# intersection off (models/Thickness.md sec 5). The volumes are
# intersection-off's, which the sweep checked.
thickness_case("lbox_notch_wall_inter_join_out", lbox, 7, +1.0, "pass", 423.0, True, 2)
thickness_case("tshape_post_wall_inter_join_in", tshape, 3, -1.0, "pass", 212.0, True, 2)
thickness_case("pocket_wall_inter_join_in", pocket5, 7, -1.0, "pass", 395.0, True, 2)

# A box with its vertical edges filleted, an end or a side face removed,
# inward: the removed face is tangent to the fillets, whose offsets run
# parallel to it and never meet it; the rim had no edge but the removed
# face's own and the result came back unhollowed. A tube round the tangent
# edge closes the gap, as the Arc join closes a convex edge (FreeCAD
# models/Thickness.md sec 6). Worked out by hand: the input less the
# shrunk rounded box, less the channel to the opening minus two quarter
# discs of the thickness's radius.
filletbox = Part.makeBox(10, 8, 6)
filletbox = filletbox.makeFillet(2, [filletbox.Edges[i] for i in (0, 2, 4, 6)])
thickness_case("filletbox_end_in",  filletbox, 1, -1.0, "pass", 261.1150)
thickness_case("filletbox_side_in", filletbox, 6, -1.0, "pass", 253.1150)
# Outward: the tube turns into the removed face too (the fillet's own offset
# covers the other side), and where it meets the tube round the convex top
# or bottom edge, the corner is an eighth of a sphere. The side face used to
# come back valid at 354.3386: the fillet's offset stretched round the
# cylinder until it crossed the removed plane, a lip in the opening further
# than the thickness from every face kept. Worked out by hand (Steiner): the
# whole skin 422.72, less what lies in front of the removed face, plus the
# two tubes and four eighths of a sphere.
thickness_case("filletbox_end_out",  filletbox, 1, +1.0, "pass", 403.9604)
thickness_case("filletbox_side_out", filletbox, 6, +1.0, "pass", 388.8188)
# A fillet removed: the removed face is curved, tangent to the planes beside
# it, and the tubes and corners above are built on its cylinder. The loop on
# that periodic surface searched and kept a dozen wires, the fillet's own
# outline among them; its edges all lie within the fillet's quarter turn, so
# it now walks the angles as on a plane. Inward, the floor's and the
# ceiling's offsets cut the cylinder in a circle that crosses the walls'
# offsets before it reaches the tubes, and the piece of floor cut off there
# hangs on the shell; it is dropped (models/Thickness.md sec 8). Before: outward an exception, inward a valid 647.5624 -- more than
# the input. Worked out by hand: inward the input less the cavity (its
# section by Green's theorem), outward (Steiner) the whole skin less the
# fillet's slab and two quarter tori, plus two tubes and four sphere pieces.
thickness_case("filletbox_fillet_in",  filletbox, 3, -1.0, "pass", 267.0193)
thickness_case("filletbox_fillet_out", filletbox, 3, +1.0, "pass", 405.9034)
# The same faces with the Intersection join, which builds no tubes: the
# neighbour's offset ran round its cylinder to the removed face (a lip,
# 382.017 for the end face outward) or never met it (inward: the input came
# back unhollowed; a fillet face threw). The gap is closed by the tube's sharp
# counterpart: a strip of the neighbour's tangent plane a thickness into the
# removed face, offset with it, and a wall square to the removed face at the
# strip's far edge; cubes where the Arc join has sphere eighths (FreeCAD
# models/Thickness.md sec 9). Worked out by hand: the sharp-grown skin
# less the slab in front of the removed face, plus the squares; inward the
# channel narrowed by a square at each side; for a fillet the squares reach
# its cylinder (Green's theorem for the cavity's corner). At r=2 the wall of
# one tangent edge lies in the plane of the other side's offset; the r=2.5
# fillet has no such coincidence.
thickness_case("filletbox_end_join_in",    filletbox, 1, -1.0, "pass", 262.8319, False, 2)
thickness_case("filletbox_end_join_out",   filletbox, 1, +1.0, "pass", 422.7964, False, 2)
thickness_case("filletbox_side_join_in",   filletbox, 6, -1.0, "pass", 254.8319, False, 2)
thickness_case("filletbox_side_join_out",  filletbox, 6, +1.0, "pass", 406.7964, False, 2)
thickness_case("filletbox_fillet_join_in",  filletbox, 3, -1.0, "pass", 268.7129, False, 2)
thickness_case("filletbox_fillet_join_out", filletbox, 3, +1.0, "pass", 424.7690, False, 2)
filletbox25 = Part.makeBox(10, 8, 6)
filletbox25 = filletbox25.makeFillet(2.5, [filletbox25.Edges[i] for i in (0, 2, 4, 6)])
thickness_case("filletbox25_fillet_join_in",  filletbox25, 3, -1.0, "pass", 258.4221, False, 2)
thickness_case("filletbox25_fillet_join_out", filletbox25, 3, +1.0, "pass", 407.4611, False, 2)

nocrash_case("pocket_inter_join_no_crash",
             Part.makeBox(10, 10, 6).cut(Part.makeBox(6, 6, 3, App.Vector(2, 2, 3))),
             7, -1.0, True, 2)

# A face closed at a pole. A dome -- half a sphere, its flat face removed --
# leaves one spherical face, bounded by the seam, the pole's degenerated edge
# and the equator. The loop keeps degenerated edges out of its vertex map, so
# the seam wire (closed edge, seam, closed edge) found no closed edge at the
# pole and was never built: the offset sphere came back with the equator as
# its only wire, a face of no area, the result invalid (216.5455 outward and
# 231.5055 inward by 0.5; upstream is right). The band is now closed by the
# pole's edge. The rim is closed in the removed face's plane, so each volume
# is a difference of sphere caps, pi h^2 (3R - h) / 3; the cone's skin has
# the apex rounded by a ball of the thickness (models/Thickness.md, "Sec
# 16").
dome = Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), 0, 90, 360)
thickness_case("dome_flat_out",       dome, 2, +0.5, "pass", 86.6556)
thickness_case("dome_flat_in",        dome, 2, -0.5, "pass", 70.9476)
thickness_case("dome_flat_join_out",  dome, 2, +0.5, "pass", 86.6556, False, 2)
thickness_case("dome_flat_inter_in",  dome, 2, -0.5, "pass", 70.9476, True)
thickness_case("bowl_flat_out",
               Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), -90, 0, 360),
               2, +0.5, "pass", 86.6556)
thickness_case("cap_flat_in",
               Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), 30, 90, 360),
               2, -0.5, "pass", 33.6412)
thickness_case("segment_flat_out",
               Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), -90, 30, 360),
               2, +0.5, "pass", 127.8890)
thickness_case("cone_base_out", Part.makeCone(0, 4, 6), 2, +0.5, "pass", 52.4095)
thickness_case("cone_up_base_out", Part.makeCone(4, 0, 6), 2, +0.5, "pass", 52.4095)

# The cone with its apex, the other ways. The offset cone is trimmed again at
# its own apex, and the edge re-trimmed there kept an infinite range on its
# 3d line: inward with the Arc join an intersection on it got a parameter of
# 2e100 and threw; with the Intersection join, apex down, the face's bounds
# were infinite, the enlarged cone kept the wrong side of its apex, never met
# the base's plane, and the cone came back unhollowed, 100.5310 (upstream the
# same in all of them). Inward the skin is the cone less the same cone with
# its apex 0.5 / sin a higher: 38.8428; outward and sharp, 52.4563.
for name, cone in (("cone", Part.makeCone(0, 4, 6)), ("cone_up", Part.makeCone(4, 0, 6))):
    thickness_case(name + "_base_in",       cone, 2, -0.5, "pass", 38.8428)
    thickness_case(name + "_base_join_in",  cone, 2, -0.5, "pass", 38.8428, False, 2)
    thickness_case(name + "_base_join_out", cone, 2, +0.5, "pass", 52.4563, False, 2)

# Half a dome: a quarter ball, its sphere bounded by two meridians meeting at
# the pole, its flat side in two coplanar faces (Face3, Face4) and its bottom
# Face2. Bottom removed, outward: the pole is a vertex on a free border with a
# tube on each meridian, two images, and rebinding the first threw
# (BRepAlgo_Image::Bind). A side removed: the other side is its tangent
# neighbour, closed by a tube round the axis whose edge on the removed face
# crosses the face's own outline -- a crossing inside both edges, which got
# no vertex; the tube kept the edge whole where the rim took its pieces; and
# the pole's edge was dropped where the wire reached the pole on a vertex of
# its own. Inward, the offset sphere came with two more wires than it has.
# Hand values, models/Thickness.md "Sec 17".
halfdome = Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), 0, 90, 180)
thickness_case("halfdome_bottom_out",      halfdome, 2, +0.5, "pass", 66.1779)
thickness_case("halfdome_bottom_in",       halfdome, 2, -0.5, "pass", 51.3127)
thickness_case("halfdome_bottom_join_in",  halfdome, 2, -0.5, "pass", 51.3127, False, 2)
thickness_case("halfdome_side_out",        halfdome, 3, +0.5, "pass", 79.7628)
thickness_case("halfdome_side_in",         halfdome, 3, -0.5, "pass", 58.8944)
thickness_case("halfdome_other_side_out",  halfdome, 4, +0.5, "pass", 79.7628)
# The mirror image inward: the removed side's meridian is stretched over the
# pole and down the far side, and the piece past the pole closed a wire round
# what the neighbour's offset cut away. Its area told nothing (its pcurve
# runs past the pole, out of the sphere's parameters): right for one side,
# kept for the other. A piece beyond a stretched edge's own ends is outside.
thickness_case("halfdome_other_side_in",   halfdome, 4, -0.5, "pass", 58.8944)
# A side inward with the Intersection join: the closing wall ends on the
# sphere at the pole, where its end edge -- a line -- has no pcurve to be
# called convex by, and the wall was never cut by the sphere; and of the two
# half circles the removed side's plane cuts the enlarged sphere in, the one
# beyond the pole was kept. As the Arc join's 58.8944, with the square column
# beside the axis, 0.9954, where that has the quarter tube, 0.7827: 59.1071.
thickness_case("halfdome_side_join_in",       halfdome, 3, -0.5, "pass", 59.1071, False, 2)
thickness_case("halfdome_other_side_join_in", halfdome, 4, -0.5, "pass", 59.1071, False, 2)
thickness_case("halfdome_side_inter_join_in", halfdome, 3, -0.5, "pass", 59.1071, True, 2)
# The Intersection join outward: the flat neighbours' offsets cut the offset
# sphere behind its pole, and the sphere face, grown past the meridians that
# bound it, has to run round the pole -- a whole turn of U with a seam, where
# U grows by a tenth of what is left of the turn. The face is put on the same
# sphere with its axis turned, where the region is a plain patch
# (BRepOffset_Tool::EnLargeFace). The sharp skin: the ball of 5.5 above the
# bottom's plane and beside the side's plane a thickness out, less the
# quarter ball, 67.0206; a side removed, the bottom's plane a thickness down
# and the neighbour's a thickness out as far as the wall, 81.7345.
thickness_case("halfdome_bottom_join_out",       halfdome, 2, +0.5, "pass", 67.0206, False, 2)
thickness_case("halfdome_bottom_inter_join_out", halfdome, 2, +0.5, "pass", 67.0206, True, 2)
thickness_case("halfdome_side_join_out",         halfdome, 3, +0.5, "pass", 81.7345, False, 2)
thickness_case("halfdome_other_side_join_out",   halfdome, 4, +0.5, "pass", 81.7345, False, 2)
# The face reaching both poles, a quarter turn of it: the same quarter ball
# lying on its side, each flat side one face.
lune = Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), -90, 90, 90)
thickness_case("lune_side_join_out", lune, 2, +0.5, "pass", 67.0206, False, 2)
thickness_case("lune_side_join_in",  lune, 3, -0.5, "pass", 51.3127, False, 2)
# An eighth of a ball, and a third of a dome (its sides not coplanar).
eighth = Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), 0, 90, 90)
thickness_case("eighth_ball_bottom_join_out", eighth, 2, +0.5, "pass", 46.7279, False, 2)
thickness_case("eighth_ball_side_join_out",   eighth, 3, +0.5, "pass", 46.7279, False, 2)
dome120 = Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), 0, 90, 120)
thickness_case("dome120_bottom_join_out", dome120, 2, +0.5, "pass", 53.3701, False, 2)
thickness_case("dome120_side_join_out",   dome120, 3, +0.5, "pass", 57.4660, False, 2)
thickness_case("dome120_side_join_in",    dome120, 3, -0.5, "pass", 41.2951, False, 2)

# Sec 19. Three quarters of a dome: more than half a turn, and the edge
# at its axis concave. The bottom removed outward came back the input,
# 196.3495; a side removed was refused, or invalid with the Arc join inward
# (the wall runs on past the pole to the other side's offset). Sharp: the
# ball of 5.5 beside the kept planes a thickness out, less the dome; inward
# the dome less the ball of 4.5 a thickness in from the kept faces. At a
# thickness of 1 the side removed outward was a valid solid on the wrong
# side of the kept one, 89.196.
dome270 = Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), 0, 90, 270)
thickness_case("dome270_bottom_join_out",       dome270, 2, +0.5, "pass", 87.3133, False, 2)
thickness_case("dome270_bottom_inter_join_out", dome270, 2, +0.5, "pass", 87.3133, True, 2)
thickness_case("dome270_side_join_out",         dome270, 3, +0.5, "pass", 113.7486, False, 2)
thickness_case("dome270_other_side_join_out",   dome270, 4, +0.5, "pass", 113.7486, False, 2)
thickness_case("dome270_side_join_in",          dome270, 3, -0.5, "pass", 83.7681, False, 2)
thickness_case("dome270_other_side_join_in",    dome270, 4, -0.5, "pass", 83.7681, False, 2)
thickness_case("dome270_side_in",               dome270, 3, -0.5, "pass", 83.7681)
thickness_case("dome270_other_side_in",         dome270, 4, -0.5, "pass", 83.7681)
thickness_case("dome270_side_join_out_thick",   dome270, 3, +1.0, "pass", 260.9367, False, 2)
# The Arc join: three tubes meet at the pole of an eighth of a ball or a
# third of a dome, the bottom removed outward, and the piece of sphere
# between them was never built; a third of a dome, a side removed inward,
# had the other side's offset extended into three quarters of a disc.
thickness_case("eighth_ball_bottom_out", eighth,  2, +0.5, "pass", 45.5612)
thickness_case("dome120_bottom_out",     dome120, 2, +0.5, "pass", 52.4334)
thickness_case("dome120_side_in",        dome120, 3, -0.5, "pass", 41.2951)
thickness_case("dome120_other_side_in",  dome120, 4, -0.5, "pass", 41.2951)
# The sphere itself removed: the wall lies on the ball of 5, between the
# outline and the kept planes' offsets, past the pole. Outward every shape
# was refused or threw; inward two threw with the Arc join.
thickness_case("eighth_ball_sphere_out",      eighth,   1, +0.5, "pass", 32.3576)
thickness_case("eighth_ball_sphere_join_out", eighth,   1, +0.5, "pass", 33.2167, False, 2)
thickness_case("eighth_ball_sphere_in",       eighth,   1, -0.5, "pass", 25.7418)
thickness_case("dome120_sphere_out",          dome120,  1, +0.5, "pass", 35.2709)
thickness_case("dome120_sphere_join_out",     dome120,  1, +0.5, "pass", 35.8993, False, 2)
thickness_case("halfdome_sphere_out",         halfdome, 1, +0.5, "pass", 41.0976)
thickness_case("halfdome_sphere_join_out",    halfdome, 1, +0.5, "pass", 41.6307, False, 2)
thickness_case("halfdome_sphere_in",          halfdome, 1, -0.5, "pass", 36.6474)
thickness_case("lune_sphere_out",             lune,     1, +0.5, "pass", 41.0976)
thickness_case("lune_sphere_join_out",        lune,     1, +0.5, "pass", 41.6307, False, 2)
thickness_case("lune_sphere_in",              lune,     1, -0.5, "pass", 36.6474)
thickness_case("dome270_sphere_out",          dome270,  1, +0.5, "pass", 49.5532)
thickness_case("dome270_sphere_join_out",     dome270,  1, +0.5, "pass", 50.0446, False, 2)
thickness_case("dome270_sphere_in",           dome270,  1, -0.5, "pass", 47.3132)
thickness_case("dome270_sphere_join_in",      dome270,  1, -0.5, "pass", 47.5529, False, 2)
thickness_case("dome270_sphere_out_thick",      dome270, 1, +1.0, "pass", 99.0413)
thickness_case("dome270_sphere_join_out_thick", dome270, 1, +1.0, "pass", 100.7985, False, 2)
# Half a ball cut through both its poles. Its sphere's outline is a whole
# great circle and no axis clears it; grown past that circle it holds its
# poles inside. It is cut in two first (MakeThickSolidOfSplit): along its
# equator where it stays, along a meridian where it is removed. 113.0909 and
# 89.0272 for a flat half removed with the Intersection join, the wall a
# thickness past the axis; 39.1390 for the sphere removed, a slab of the
# ball half a unit thick.
halfball = Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), -90, 90, 180)
thickness_case("halfball_flat_join_out",   halfball, 2, +0.5, "pass", 113.0909, False, 2)
thickness_case("halfball_flat_join_in",    halfball, 2, -0.5, "pass", 89.0272, False, 2)
thickness_case("halfball_flat2_join_out",  halfball, 3, +0.5, "pass", 113.0909, False, 2)
thickness_case("halfball_flat2_join_in",   halfball, 3, -0.5, "pass", 89.0272, False, 2)
thickness_case("halfball_sphere_out",      halfball, 1, +0.5, "pass", 39.1390)
thickness_case("halfball_sphere_join_out", halfball, 1, +0.5, "pass", 39.1390, False, 2)
thickness_case("halfball_sphere_in",       halfball, 1, -0.5, "pass", 39.1390)
thickness_case("halfball_sphere_join_in",  halfball, 1, -0.5, "pass", 39.1390, False, 2)
thickness_case("halfball_flat_out",        halfball, 2, +0.5, "pass", 111.6001)
thickness_case("halfball_flat_in",         halfball, 2, -0.5, "pass", 88.5482)
# The same, placed.
thickness_case("placed_halfball_flat_join_out", placed(halfball), 2, +0.5, "pass", 113.0909, False, 2)
thickness_case("placed_halfball_sphere_out",    placed(halfball), 1, +0.5, "pass", 39.1390)

# Half a ball that comes with its sphere in two faces already -- a fuse not
# refined. In two lunes, both removed: the edge between them is stretched
# past the poles for the wall, where a twin's pcurve, interpolated between
# the edge's ends, was no curve (invalid), and inward the section of the two
# lunes' offsets has a block without an edge (a segmentation fault). In two
# domes, a flat half removed with the Intersection join: the equator between
# them is stretched most of a turn, far past what its pcurve on a turned
# sphere was prolonged to -- valid, 0.23 off under a tolerance to match, and
# 90.2265 for 89.0272.
luneball = Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), -90, 90, 180).generalFuse(
    [Part.Arc(App.Vector(0, 0, -5), App.Vector(0, 5, 0), App.Vector(0, 0, 5)).toShape()])[0].Solids[0]
pieces_case("luneball_spheres_out",      luneball, [1, 2], +0.5, [39.1390])
pieces_case("luneball_spheres_join_out", luneball, [1, 2], +0.5, [39.1390], False, 2)
pieces_case("luneball_spheres_in",       luneball, [1, 2], -0.5, [39.1390])
pieces_case("luneball_spheres_join_in",  luneball, [1, 2], -0.5, [39.1390], False, 2)
eqball = Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), -90, 90, 180).generalFuse(
    [Part.ArcOfCircle(Part.Circle(App.Vector(), App.Vector(0, 0, 1), 5), 0, math.pi).toShape()])[0].Solids[0]
thickness_case("eqball_flat_join_out", eqball, 3, +0.5, "pass", 113.0909, False, 2)
thickness_case("eqball_flat_join_in",  eqball, 3, -0.5, "pass", 89.0272, False, 2)
thickness_case("eqball_flat_out",      eqball, 3, +0.5, "pass", 111.6001)
thickness_case("eqball_flat_in",       eqball, 3, -0.5, "pass", 88.5482)

# A box fused of two and not refined: every face across the joint is in two
# coplanar pieces. One piece of the top removed: its neighbour is tangent to
# it, closed by a tube with the Arc join (417.3038 outward: the rounded skin
# 455.5869 less the slab, three quarter tubes and two ball eighths over the
# removed half, plus the tube and two eighths at the joint; 274.7124 inward:
# 480 - 192 - 6 (3 - pi / 4)).
splitbox = Part.makeBox(4, 8, 6).fuse(Part.makeBox(6, 8, 6, App.Vector(4, 0, 0)))
split_top = [i + 1 for i, f in enumerate(splitbox.Faces)
             if abs(f.BoundBox.ZMin - 6) < 1e-9 and f.BoundBox.XMax < 4.5][0]
thickness_case("splitbox_top_piece_out", splitbox, split_top, +1.0, "pass", 417.3038)
thickness_case("splitbox_top_piece_in",  splitbox, split_top, -1.0, "pass", 274.7124)
# The Intersection join closes the gap with a wall a thickness into the
# removed piece: the sharp skin, 480 either way, less the slab over the
# removed piece short of the wall, 4 x 10 x 1 (440), or less the cavity 192
# and the shaft up to the wall, 2 x 6 x 1 (276). Three faults, models/
# Thickness.md "Sec 18": the wall and the face at the joint's end were
# intersected twice (their section two edges, the wall never built); coplanar
# pieces meeting at a vertex were intersected with each other's neighbours
# (two edges on one line, no piece of a split face built -- an end face
# removed failed too, and upstream gives the box back); and the side beside
# the removed piece had no edge to the other piece of the top, between the
# wall and the joint.
thickness_case("splitbox_top_piece_join_out", splitbox, split_top, +1.0, "pass", 440.0, False, 2)
thickness_case("splitbox_top_piece_join_in",  splitbox, split_top, -1.0, "pass", 276.0, False, 2)
thickness_case("splitbox_top_piece_inter_join_out", splitbox, split_top, +1.0, "pass", 440.0, True, 2)
thickness_case("splitbox_top_piece_inter_join_in",  splitbox, split_top, -1.0, "pass", 276.0, True, 2)
# An end face: no removed face beside a coplanar one, the joint alone.
split_end = [i + 1 for i, f in enumerate(splitbox.Faces) if f.BoundBox.XMax < 1e-9][0]
thickness_case("splitbox_end_join_out", splitbox, split_end, +1.0, "pass", 400.0, False, 2)
thickness_case("splitbox_end_join_in",  splitbox, split_end, -1.0, "pass", 264.0, False, 2)
# A piece of a side, and the larger piece of the top: 480 - 4 x 8 x 1 and
# 480 - 192 - 2 x 4 x 1; 480 - 6 x 10 x 1 and 480 - 192 - 4 x 6 x 1.
split_side = [i + 1 for i, f in enumerate(splitbox.Faces)
              if f.BoundBox.YMax < 1e-9 and f.BoundBox.XMax < 4.5][0]
split_top_2 = [i + 1 for i, f in enumerate(splitbox.Faces)
               if abs(f.BoundBox.ZMin - 6) < 1e-9 and f.BoundBox.XMin > 3.5][0]
thickness_case("splitbox_side_piece_join_out", splitbox, split_side, +1.0, "pass", 448.0, False, 2)
thickness_case("splitbox_side_piece_join_in",  splitbox, split_side, -1.0, "pass", 280.0, False, 2)
thickness_case("splitbox_top_other_piece_join_out", splitbox, split_top_2, +1.0, "pass", 420.0, False, 2)
thickness_case("splitbox_top_other_piece_join_in",  splitbox, split_top_2, -1.0, "pass", 264.0, False, 2)
# The top alone in two pieces, the sides whole: the face at the joint's end
# is one face.
from BOPTools import SplitAPI
topsplit = SplitAPI.slice(Part.makeBox(10, 8, 6),
                          [Part.makeLine(App.Vector(4, 0, 6), App.Vector(4, 8, 6))],
                          "Split").Solids[0]
topsplit_top = [i + 1 for i, f in enumerate(topsplit.Faces)
                if abs(f.BoundBox.ZMin - 6) < 1e-9 and f.BoundBox.XMax < 4.5][0]
thickness_case("topsplit_top_piece_join_out", topsplit, topsplit_top, +1.0, "pass", 440.0, False, 2)
thickness_case("topsplit_top_piece_join_in",  topsplit, topsplit_top, -1.0, "pass", 276.0, False, 2)


# Where the fork has no answer yet it must have none: the input back, "valid"
# and unhollowed, is the one result that is never right. MakeThickSolid
# refuses a result of the shape's own volume (inward) or of its volume with
# the removed faces whole again (outward).
def never_the_input_case(name, shape, face_index, value, inter=False, join=0):
    before = signature(shape)
    try:
        r = shape.makeThickness([shape.Faces[face_index - 1]], value, 1e-7, inter, False, 0, join)
        back = r.isValid() and abs(r.Volume - shape.Volume) <= RELTOL * abs(shape.Volume)
        problems = input_problems(shape, before)
        if back:
            problems.append("the input back, vol=%.4f" % r.Volume)
        report(name, not problems, False, "; ".join(problems) if problems else "vol=%.4f" % r.Volume)
    except Exception:
        problems = input_problems(shape, before)
        report(name, not problems, False, "; ".join(problems) if problems else "refused")


never_the_input_case("splitbox_top_piece_join_in_not_the_box", splitbox, split_top, -1.0, False, 2)
never_the_input_case("halfdome_side_join_in_not_the_input", halfdome, 3, -0.5, False, 2)
never_the_input_case("halfdome_side_join_out_not_the_input", halfdome, 3, +0.5, False, 2)

# Every face removed: no face stays to be thickened, and the call is refused.
# A sphere with its one face removed came back as the sphere itself, "valid"
# and unhollowed; the torus was refused already (FreeCAD
# models/Thickness.md sec 15).
def refused_case(name, shape, face_indices, value, inter=False, join=0):
    before = signature(shape)
    try:
        r = shape.makeThickness([shape.Faces[i - 1] for i in face_indices], value, 1e-7,
                                inter, False, 0, join)
        report(name, False, False, "answered: %s valid=%s vol=%.4f" % (
            r.ShapeType, r.isValid(), r.Volume))
    except Exception:
        problems = input_problems(shape, before)
        report(name, not problems, False, "; ".join(problems) if problems else "refused")


for value, tag in ((+1.0, "out"), (-1.0, "in")):
    refused_case("sphere_face_refused_" + tag, Part.makeSphere(5), [1], value)
    refused_case("torus_face_refused_" + tag, Part.makeTorus(8, 2), [1], value)
refused_case("sphere_face_join_refused_in", Part.makeSphere(5), [1], -1.0, False, 2)
refused_case("box_all_faces_refused_in", Part.makeBox(10, 10, 6), [1, 2, 3, 4, 5, 6], -1.0)
refused_case("box_all_faces_inter_refused_out", Part.makeBox(10, 10, 6), [1, 2, 3, 4, 5, 6], +1.0,
             True)


# The input left as it was. A pad of a ring sector straddling angle 0, its
# outer arc's face and both caps removed: the loop on the removed cylinder
# found its test face invalid and ran ShapeFix on it, which shifted the pcurve
# of an edge the test face shares with the input by a period, on the
# cylinder's own surface -- the input came back inside out (FreeCAD's
# PartDesign Pad under a Thickness, issue3 above, with shape values not
# frozen; occ-issues local03, models/Thickness.md sec 14).
# Only the input and the result's validity are checked here.
def input_case(name, shape, face_indices, value, inter=False, join=0):
    try:
        before = signature(shape)
        r = shape.makeThickness([shape.Faces[i - 1] for i in face_indices], value, 1e-7,
                                inter, False, 0, join)
        problems = input_problems(shape, before)
        if not r.isValid():
            problems.append("result invalid")
        detail = "; ".join(problems) if problems else "input intact, result vol=%.4f" % r.Volume
        report(name, not problems, False, detail)
    except Exception as e:
        report(name, False, False, "EXCEPTION " + str(e).strip().splitlines()[-1])


def ring_sector(past):
    """A ring sector of radii 46.85 and 59.3 straddling angle 0, padded 27:
    its arcs' parameters run from past - a to past + a."""
    V = App.Vector
    outer = Part.ArcOfCircle(Part.Circle(V(0, 0, 0), V(0, 0, 1), 59.3), past - 0.29,
                             past + 0.29).toShape()
    inner = Part.ArcOfCircle(Part.Circle(V(0, 0, 0), V(0, 0, 1), 46.85), past - 0.23,
                             past + 0.23).toShape()
    sides = [Part.makeLine(inner.Vertexes[k].Point, outer.Vertexes[k].Point) for k in (0, 1)]
    wire = Part.Wire(Part.__sortEdges__([outer, sides[0], inner, sides[1]]))
    return Part.Face(wire).extrude(V(0, 0, 27))


# Face3 the outer arc's, Face5 and Face6 the caps.
input_case("sector_outer_arc_input", ring_sector(0.0), [3, 5, 6], 1.0)
input_case("sector_outer_arc_past_period_input", ring_sector(2 * math.pi), [3, 5, 6], 1.0)


# Determinism: thickness with intersection and the Arc join iterated its
# offsets in hash order (a shape's hash is its TShape's address), and the
# result depended on that order -- up to four different results in eight runs
# of one input. Each case must give one result over repeated runs on fresh
# shapes (models/Thickness.md sec 1).
def determinism_case(name, make, face_index, value, runs=8):
    outcomes = set()
    for _ in range(runs):
        shape = make()
        try:
            r = shape.makeThickness([shape.Faces[face_index - 1]], value, 1e-3,
                                    True, False, 0, 0)
            outcomes.add("valid=%s vol=%.4f" % (r.isValid(), r.Volume))
        except Exception as e:
            outcomes.add("EXCEPTION " + type(e).__name__)
    report(name, len(outcomes) == 1, False, "; ".join(sorted(outcomes)))


determinism_case("arc_inter_boss_same_every_run",
                 lambda: Part.makeCylinder(6, 4).fuse(
                     Part.makeCylinder(3, 4, App.Vector(0, 0, 4))), 3, +1.0)
determinism_case("arc_inter_lbox_same_every_run",
                 lambda: Part.makeBox(10, 10, 5).cut(
                     Part.makeBox(5, 5, 5, App.Vector(5, 5, 0))), 7, +1.0)
determinism_case("arc_inter_boxhole_same_every_run",
                 lambda: Part.makeBox(10, 10, 5).cut(
                     Part.makeCylinder(2, 3, App.Vector(5, 5, 2))), 2, +1.0)

# A shape that carries a location -- an object with a placement -- is the
# shape it is where it was made. Four faults, three of them upstream's, each
# met only under a location: a seam's second pcurve checked without the
# edge's location (BRepCheck_Edge::Tolerance), which gave a new seam edge and
# the input vertices on it a tolerance as large as the move -- the thickness
# after it on the same shape was a valid solid of 272.73 for 89.93; a pole
# not known for one where the surface is sampled (CheckInputData), every
# shape with a pole refused; the apex of an offset cone left where the cone
# was made (BRepOffset_Offset); the normal of a whole circle taken from its
# start, middle and end (CorrectConicalFaces), which is rounding, or nothing.
# And the turned sphere of sec 18 was not made for a face that is placed.
# A shape turned in space is the same shape too. The corner arc of a removed
# face tangent to its neighbour was taken from a section whose lines join
# where the intersection cares to cut them, and for a filleted box turned by
# 40 degrees its ends fell on two of them: no arc, no sphere at the corner,
# and a valid solid of 459.398 for 405.903.
def turned(shape):
    """The shape with its geometry turned and moved, no location on it."""
    moved = shape.copy()
    moved.transformShape(
        App.Placement(App.Vector(3, 4, 5), App.Rotation(App.Vector(1, 2, 3), 40)).Matrix, True)
    return moved


for _i in (3, 4, 8, 9):
    thickness_case("turned_filletbox_fillet%d_out" % _i, turned(filletbox), _i, +1.0, "pass", 405.9034)
thickness_case("turned_filletbox_fillet_in", turned(filletbox), 3, -1.0, "pass", 267.0193)
thickness_case("placed_filletbox_fillet_out", placed(filletbox), 4, +1.0, "pass", 405.9034)

# Known broken, each with the volume it should give (models/Thickness.md,
# sec 20, "Found beside these").
# - (answered) The half ball as a cut or a refine leaves it: one disc for its
#   flat. The disc removed with the Intersection join was refused (and
#   searched its loops without end until the search was given a number of
#   steps), the sphere removed was refused outward. The sphere is put on one
#   turned onto its middle: a dome with its rim in two arcs.
# - A ball wedge from pole to pole on more than half a turn, its sphere
#   removed: three lunes outward, refused; on 240 degrees inward with the
#   Arc join a valid solid of 40.4447. On 150 degrees outward, Arc: refused.
# - (answered) The half ball cut the other way round: in two lunes with a
#   flat half removed, in two domes with both removed. Faces of one sphere
#   that have one thing done with them are joined first.
onedisc = Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), -90, 90, 180).removeSplitter()
thickness_case("halfball1_disc_out",         onedisc, 2, +0.5, "pass", 86.6556)
thickness_case("halfball1_disc_in",          onedisc, 2, -0.5, "pass", 70.9476)
thickness_case("halfball1_sphere_join_in",   onedisc, 1, -0.5, "pass", 39.1390, False, 2)
thickness_case("halfball1_disc_join_out",    onedisc, 2, +0.5, "pass", 86.6556, False, 2)
thickness_case("halfball1_disc_join_in",     onedisc, 2, -0.5, "pass", 70.9476, False, 2)
thickness_case("halfball1_sphere_out",       onedisc, 1, +0.5, "pass", 39.1390)
thickness_case("halfball1_sphere_join_out",  onedisc, 1, +0.5, "pass", 39.1390, False, 2)
thickness_case("halfball1_sphere_in",        onedisc, 1, -0.5, "pass", 39.1390)
# The same from a cut, which hands over a compound of one solid; and under a
# location, and turned in space.
_cutball = Part.makeSphere(5).cut(Part.makeBox(20, 20, 20, App.Vector(-10, -20, -10)))
thickness_case("cutball_disc_join_out",      _cutball, 2, +0.5, "pass", 86.6556, False, 2)
thickness_case("cutball_disc_join_in",       _cutball, 2, -0.5, "pass", 70.9476, False, 2)
thickness_case("cutball_sphere_out",         _cutball, 1, +0.5, "pass", 39.1390)
thickness_case("cutball_sphere_join_out",    _cutball, 1, +0.5, "pass", 39.1390, False, 2)
thickness_case("cutball_sphere_join_in",     _cutball, 1, -0.5, "pass", 39.1390, False, 2)
thickness_case("placed_halfball1_disc_join_out",   placed(onedisc), 2, +0.5, "pass", 86.6556, False, 2)
thickness_case("placed_halfball1_sphere_join_out", placed(onedisc), 1, +0.5, "pass", 39.1390, False, 2)
thickness_case("turned_halfball1_disc_join_in",    turned(onedisc), 2, -0.5, "pass", 70.9476, False, 2)
thickness_case("turned_halfball1_sphere_out",      turned(onedisc), 1, +0.5, "pass", 39.1390)
ball270 = Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), -90, 90, 270)
ball240 = Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), -90, 90, 240)
ball150 = Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), -90, 90, 150)
thickness_case("ball270_sphere_in",          ball270, 1, -0.5, "pass", 41.0978)
thickness_case("ball270_sphere_join_in",     ball270, 1, -0.5, "pass", 41.6307, False, 2)
thickness_case("ball270_sphere_out",         ball270, 1, +0.5, "xfail", 36.6474)
thickness_case("ball270_sphere_join_out",    ball270, 1, +0.5, "xfail", 36.6474, False, 2)
thickness_case("ball240_sphere_out",         ball240, 1, +0.5, "xfail", 37.6996)
thickness_case("ball240_sphere_join_out",    ball240, 1, +0.5, "xfail", 37.6996, False, 2)
thickness_case("ball240_sphere_in",          ball240, 1, -0.5, "xfail", 39.9613)
thickness_case("ball240_sphere_join_in",     ball240, 1, -0.5, "pass", 40.5784, False, 2)
thickness_case("ball150_sphere_out",         ball150, 1, +0.5, "xfail", 39.7919)
_meridian = Part.Arc(App.Vector(0, 0, -5), App.Vector(0, 5, 0), App.Vector(0, 0, 5)).toShape()
_equator = Part.ArcOfCircle(Part.Circle(App.Vector(), App.Vector(0, 0, 1), 5), 0, math.pi).toShape()
_halfball = Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), -90, 90, 180)
_lunes = _halfball.generalFuse([_meridian])[0].Solids[0]
_domes = _halfball.generalFuse([_equator])[0].Solids[0]
thickness_case("luneball_flat_join_out",     _lunes, 3, +0.5, "pass", 113.0909, False, 2)
thickness_case("luneball_flat_join_in",      _lunes, 3, -0.5, "pass", 89.0272, False, 2)
thickness_case("eqball_spheres_join_in",     _domes, [1, 2], -0.5, "pass", 39.1390, False, 2)
thickness_case("eqball_spheres_out",         _domes, [1, 2], +0.5, "pass", 39.1390)
thickness_case("eqball_spheres_in",          _domes, [1, 2], -0.5, "pass", 39.1390)
thickness_case("eqball_spheres_join_out",    _domes, [1, 2], +0.5, "pass", 39.1390, False, 2)

# A dome whose rim is in two arcs: the seam's band was only built on a rim
# that is one closed edge, and the offset of the dome came out as a face of no
# area -- an invalid solid of 240.68 for 86.6556 where upstream is right. The
# face is walked in (u, v) with its seam instead. With the sphere removed the
# flat's two arcs were both replaced by the one half of the section circle;
# with the flat in two halves as well, the circle was cut by the line between
# them at the one point it starts on.
_dome_y = Part.makeSphere(5, App.Vector(), App.Vector(0, 1, 0), 0, 90, 360)
_rim = [e for e in _dome_y.Edges if not e.Degenerated and abs(e.Length - 10 * math.pi) < 1e-6][0]
_seam_end = _rim.Vertexes[0].Point
dome2 = _dome_y.generalFuse([Part.Vertex(_seam_end * -1)])[0].Solids[0]
dome2f = _dome_y.generalFuse([Part.makeLine(_seam_end, _seam_end * -1)])[0].Solids[0]
thickness_case("dome2_flat_out",          dome2, 2, +0.5, "pass", 86.6556)
thickness_case("dome2_flat_in",           dome2, 2, -0.5, "pass", 70.9476)
thickness_case("dome2_flat_join_out",     dome2, 2, +0.5, "pass", 86.6556, False, 2)
thickness_case("dome2_flat_join_in",      dome2, 2, -0.5, "pass", 70.9476, False, 2)
thickness_case("dome2_sphere_out",        dome2, 1, +0.5, "pass", 39.1390)
thickness_case("dome2_sphere_in",         dome2, 1, -0.5, "pass", 39.1390)
thickness_case("dome2_sphere_join_out",   dome2, 1, +0.5, "pass", 39.1390, False, 2)
thickness_case("dome2_sphere_join_in",    dome2, 1, -0.5, "pass", 39.1390, False, 2)
thickness_case("dome2f_sphere_out",       dome2f, 1, +0.5, "pass", 39.1390)
thickness_case("dome2f_sphere_in",        dome2f, 1, -0.5, "pass", 39.1390)
thickness_case("dome2f_sphere_join_out",  dome2f, 1, +0.5, "pass", 39.1390, False, 2)
thickness_case("dome2f_sphere_join_in",   dome2f, 1, -0.5, "pass", 39.1390, False, 2)
# Known broken: one half of the flat removed, the Intersection join.
thickness_case("dome2f_flat_join_out",    dome2f, 2, +0.5, "xfail", 113.0909, False, 2)
thickness_case("dome2f_flat_join_in",     dome2f, 2, -0.5, "xfail", 89.0272, False, 2)

# A circle of section is a whole turn, closed on a vertex of its own that the
# intersection puts where it likes -- elsewhere once the shape is turned. Cut
# by the edges that cross it, the piece that vertex lies in was lost
# (TrimEdge, ExtentFace), and it was the piece needed as often as not: half
# of a sphere's cap, a face removed inward, gave a valid solid of 11.598 for
# 22.043 and of 8.126 for 19.232 turned, or was refused. And a cylinder cuts
# the sphere round its end in two circles as far from their edge as each
# other, of which the first found was taken: a dome on a cylinder, the
# cylinder removed, 151.12 for 100.73 turned.
halfcap = Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), 30, 90, 180)
bullet = Part.makeCylinder(5, 4, App.Vector(0, 0, -4)).fuse(
    Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), 0, 90, 360)).removeSplitter()
for _tag, _how in (("", lambda s: s), ("turned_", turned), ("placed_", placed)):
    thickness_case(_tag + "halfcap_sphere_in",          _how(halfcap), 1, -0.5, "pass", 19.2316)
    thickness_case(_tag + "halfcap_sphere_join_in",     _how(halfcap), 1, -0.5, "pass", 19.2316, False, 2)
    thickness_case(_tag + "halfcap_bottom_join_in",     _how(halfcap), 2, -0.5, "pass", 22.0431, False, 2)
    thickness_case(_tag + "halfcap_side_join_in",       _how(halfcap), 3, -0.5, "pass", 28.8642, False, 2)
    thickness_case(_tag + "halfcap_other_side_join_in", _how(halfcap), 4, -0.5, "pass", 28.8642, False, 2)
    pieces_case(_tag + "bullet_side_join_out", _how(bullet), [1], +0.5, [39.2699, 61.4616], False, 2)

placed_cyl = placed(Part.makeCylinder(4, 6))
pieces_case("placed_cyl_side_pieces_out", placed_cyl, [1], +0.5, [25.1327, 25.1327])
thickness_case("placed_cyl_top_in_after",   placed_cyl, 2, -0.5, "pass", 89.9281)
thickness_case("placed_cyl_top_out_after",  placed_cyl, 2, +0.5, "pass", 110.4400)
thickness_case("placed_dome_flat_out",      placed(Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), 0, 90, 360)), 2, +0.5, "pass", 86.6556)
thickness_case("placed_dome_flat_in",       placed(Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), 0, 90, 360)), 2, -0.5, "pass", 70.9476)
thickness_case("placed_cone_base_join_out", placed(Part.makeCone(0, 4, 6)), 2, +0.5, "pass", 52.4563, False, 2)
thickness_case("placed_cone_base_in",       placed(Part.makeCone(0, 4, 6)), 2, -0.5, "pass", 38.8428)
thickness_case("placed_coneup_base_out",    placed(Part.makeCone(4, 0, 6)), 2, +0.5, "pass", 52.4095)
thickness_case("placed_halfdome_bottom_join_out", placed(Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), 0, 90, 180)), 2, +0.5, "pass", 67.0206, False, 2)
thickness_case("placed_halfdome_side_join_in", placed(Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), 0, 90, 180)), 3, -0.5, "pass", 59.1071, False, 2)
thickness_case("placed_dome270_side_join_out", placed(Part.makeSphere(5, App.Vector(), App.Vector(0, 0, 1), 0, 90, 270)), 3, +0.5, "pass", 113.7486, False, 2)

counts = {}
for _, verdict, _ in results:
    counts[verdict] = counts.get(verdict, 0) + 1
emit("---")
emit("summary: " + "  ".join("%s=%d" % kv for kv in sorted(counts.items())))
failed = counts.get("FAIL", 0)
if counts.get("UNEXPECTED-PASS"):
    emit("NOTE: xfail case(s) now pass - promote them to 'pass' with their volume.")
sys.stdout.flush()
os._exit(1 if failed else 0)
