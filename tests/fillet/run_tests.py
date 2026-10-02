# Regression tests for the fork's fillet (TKFillet / ChFi3d) fixes.
#
# Run with any FreeCAD build linked against this OCCT:
#
#     FreeCADCmd tests/fillet/run_tests.py
#
# Exit code is 0 when no expected-pass case fails.  Known-broken cases are
# declared XFAIL below; when one starts passing the suite prints
# UNEXPECTED-PASS so its expectation (and README.md) can be updated.
#
# See README.md for what each case covers and the history of each problem.

import os
import sys
import traceback

import FreeCAD as App
import Part

HERE = os.path.dirname(os.path.abspath(__file__))
MODELS = os.path.join(HERE, "models")
RELTOL = 1e-4


def emit(line):
    # FreeCAD redirects sys.stdout into its own console, which swallows
    # script output under FreeCADCmd - write straight to fd 1.
    os.write(1, (line + "\n").encode())


# FreeCAD freezes shape values when it finds this fork, and a frozen input is
# protected from an algorithm that edits it. Run unfrozen by default, as any
# other build of FreeCAD runs, so such an edit shows (FILLET_FREEZE=1 for the
# frozen path).
App.ParamGet("User parameter:BaseApp/Preferences/Mod/Part").SetBool(
    "ImmutableShapeValues", os.environ.get("FILLET_FREEZE", "0") != "0")


def signature(shape):
    """What a fillet must leave of its input: the shape valid, the volume."""
    return (shape.isValid(), round(shape.Volume, 6))


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


def edge_between(shape, a, b):
    """The 1-based index of the edge of `shape` whose ends are `a` and `b`."""
    want = sorted([tuple(round(c, 6) for c in a), tuple(round(c, 6) for c in b)])
    for i, e in enumerate(shape.Edges, 1):
        if sorted(tuple(round(c, 6) for c in v.Point) for v in e.Vertexes) == want:
            return i
    raise ValueError("no edge between %s and %s" % (a, b))


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
                # the property only: getattr() falls back to Part.getShape()
                if "Shape" not in o.PropertiesList:
                    continue
                s = o.getPropertyByName("Shape")
                if hasattr(s, "isValid") and not s.isNull() and not s.isValid():
                    bad.append("%s:invalid-shape" % o.Name)
            if not bad and volumes:
                for oname, vol in volumes.items():
                    got = doc.getObject(oname).Shape.Volume
                    if abs(got - vol) > RELTOL * abs(vol):
                        bad.append("%s:volume %.4f != %.4f" % (oname, got, vol))
            report(name, not bad, False, "; ".join(bad) if bad else "ok")
        finally:
            App.closeDocument(doc.Name)
    except Exception:
        report(name, False, False, traceback.format_exc().splitlines()[-1])


# ---------------------------------------------------------------------------
# Shape cases: one fillet on a stored or built shape.
# ---------------------------------------------------------------------------

def fillet_case(name, shape, edge_index, radius, expect, ref_volume=None):
    """makeFillet(radius, [Edge<edge_index>]); edge_index may be a list.

    expect='pass':  a valid solid, one closed shell, the input left as it was,
                    and ref_volume when given.
    expect='xfail': known broken (see README.md); an exception, an invalid
                    shape or an open shell counts as the expected failure.
    """
    try:
        before = signature(shape)
        indices = edge_index if isinstance(edge_index, list) else [edge_index]
        r = shape.makeFillet(radius, [shape.Edges[i - 1] for i in indices])
        problems = []
        if signature(shape) != before:
            problems.append("input changed: valid=%s vol=%.4f" % signature(shape))
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


def nocrash_case(name, shape, edge_index, radius):
    """The fillet may fail, but must not take the process down; the input
    must be left as it was either way. A crash ends the run: the missing
    summary line is the failure."""
    before = signature(shape)
    try:
        r = shape.makeFillet(radius, [shape.Edges[edge_index - 1]])
        detail = "valid=%s vol=%.4f" % (r.isValid(), r.Volume)
    except Exception as e:
        detail = "refused: " + str(e).strip().splitlines()[-1]
    ok = signature(shape) == before
    report(name, ok, False, detail if ok else "input changed; " + detail)


# realthunder/FreeCAD#523: a 10 box with a three-quarter cylinder (r 2) at its
# corner on the z axis, made by Part's Connect. The cylinder's seam -- the
# line x=2, y=0 where it meets the box's side -- stayed an edge with both
# pcurves on the cylinder, and a fillet on the box edge ending there read
# the pcurve a period away from the face.
document_case("issue523_document", "issue523_fillet_explodes.FCStd",
              {"Fillet": 1092.5263})

bc = Part.read(os.path.join(MODELS, "issue523_box_cylinder.brep"))
top = edge_between(bc, (2, 0, 10), (10, 0, 10))     # ends on the seam, at the top
bottom = edge_between(bc, (2, 0, 0), (10, 0, 0))    # the same at the bottom
mtop = edge_between(bc, (0, 2, 10), (0, 10, 10))    # its mirror image: no seam
mbottom = edge_between(bc, (0, 2, 0), (0, 10, 0))
arc_top = edge_between(bc, (0, 2, 10), (2, 0, 10))  # the cylinder's top arc
arc_bottom = edge_between(bc, (0, 2, 0), (2, 0, 0))

# The fillets on the edges ending on the seam match their mirror images.
for r, vol in ((0.5, 1093.8183), (1.0, 1092.5263)):
    fillet_case("seam_end_top_r%g" % r, bc, top, r, "pass", vol)
    fillet_case("seam_end_bottom_r%g" % r, bc, bottom, r, "pass", vol)
    fillet_case("mirror_top_r%g" % r, bc, mtop, r, "pass", vol)
    fillet_case("mirror_bottom_r%g" % r, bc, mbottom, r, "pass", vol)
# At the cylinder's radius the fillet's end reaches the far end of the
# cylinder's face as well.
fillet_case("seam_end_top_r2", bc, top, 2.0, "pass", 1087.3009)
fillet_case("seam_end_bottom_r2", bc, bottom, 2.0, "xfail")
fillet_case("mirror_top_r2", bc, mtop, 2.0, "xfail")

# Every edge at radius 1: nothing else moved.
for i in range(1, len(bc.Edges) + 1):
    fillet_case("every_edge_r1_e%02d" % i, bc, i, 1.0, "pass")

# A fillet of the cylinder's own radius on its arcs: the fillet's line on the
# top face is a point. It crashed in IntersUpdateOnSame.
nocrash_case("arc_top_r2_no_crash", bc, arc_top, 2.0)
nocrash_case("arc_bottom_r2_no_crash", bc, arc_bottom, 2.0)

# realthunder/FreeCAD#962: a fillet whose end lies on a face split into
# coplanar pieces (a PartDesign body without Refine). A 10x10x14 block with a
# slot cut down to z=5 (y 3..6); the fillet is on the outer edge of the slot,
# x=10, y=3, which ends on the slot's floor.
V = App.Vector


def slotted(wall=False, outer=False):
    """The slotted block; `wall` splits the slot's wall (y=3) by a vertical
    edge 0.7 from the outer face, `outer` splits the outer face (x=10) by the
    spine's own line below the floor."""
    s = Part.makeBox(10, 10, 14).cut(Part.makeBox(12, 3, 10, V(-1, 3, 5)))
    for on, (a, b) in ((outer, (V(10, 3, 0), V(10, 3, 5))),
                       (wall, (V(9.3, 3, 5), V(9.3, 3, 20)))):
        if on:
            s = s.generalFuse([Part.LineSegment(a, b).toShape()])[0].Solids[0]
    return s


plain = slotted()
spine = edge_between(plain, (10, 3, 5), (10, 3, 14))
for r in (0.5, 0.8, 1.0, 2.0):
    ref = plain.makeFillet(r, [plain.Edges[spine - 1]]).Volume
    for tag, s in (("wall", slotted(wall=True)), ("outer", slotted(outer=True)),
                   ("both", slotted(wall=True, outer=True))):
        # wider than the wall's piece at the vertex from 0.7 on
        fillet_case("slot_split_%s_r%g" % (tag, r), s,
                    edge_between(s, (10, 3, 5), (10, 3, 14)), r, "pass", ref)

# The model's own shape: the Pocket002 the fillet is made on. Edges 101, 102,
# 104 and 105 are the outer edges of four slots like the block's, the wall's
# piece at the floor 0.692 and 0.787 wide; edge 50 ends where the wall
# (y=10.2) continues as another face above the end face, a chamfer.
p2 = Part.read(os.path.join(MODELS, "issue962_pocket002.brep"))
for name, a, b in (("e101", (50, 13.7, 22), (50, 13.7, 14)),
                   ("e102", (50, 17.1, 14), (50, 17.1, 22)),
                   ("e104", (50, 20.3, 22), (50, 20.3, 14)),
                   ("e105", (50, 23.7, 14), (50, 23.7, 22))):
    fillet_case("issue962_%s_r0.8" % name, p2, edge_between(p2, a, b), 0.8, "pass",
                11580.7739)
e50 = edge_between(p2, (38.5, 10.2, -3.75), (38.5, 10.2, 14))
fillet_case("issue962_e50_r0.3", p2, e50, 0.3, "pass", 11583.8075)
fillet_case("issue962_e50_r0.8", p2, e50, 0.8, "pass", 11581.7861)
# Edges 36 and 49, the foot of the block where the arm meets it (x=38.5,
# z -9.75..-3.75): each fine alone, together they gave an invalid solid 7225
# too large. The bottom is two coplanar faces, the arm's and the block's,
# split along x=38.5; one corner cut that seam, the other met its line past
# the seam's end, and both versions of the seam stayed.
e36 = edge_between(p2, (38.5, 27.2, -3.75), (38.5, 27.2, -9.75))
e49 = edge_between(p2, (38.5, 10.2, -3.75), (38.5, 10.2, -9.75))
fillet_case("issue962_e36_e49_r0.8", p2, [e36, e49], 0.8, "pass", 11584.2053)
# and the Fillet itself, all twelve edges
twelve = [edge_between(p2, a, b) for a, b in (
    ((50, 10.2, -3.75), (50, 10.2, 22)), ((50, 13.7, 22), (50, 13.7, 14)),
    ((50, 17.1, 14), (50, 17.1, 22)), ((50, 20.3, 22), (50, 20.3, 14)),
    ((50, 23.7, 14), (50, 23.7, 22)), ((50, 27.2, -3.75), (50, 27.2, 22)),
    ((38.5, 27.2, -3.75), (38.5, 27.2, 14)), ((38.5, 27.2, -3.75), (38.5, 27.2, -9.75)),
    ((38.5, 10.2, -3.75), (38.5, 10.2, 14)), ((38.5, 10.2, -3.75), (38.5, 10.2, -9.75)),
    ((10.898402, 16.470802, -3.75), (10.898402, 16.470802, -9.75)),
    ((17.356038, -9.42499, -3.75), (17.356038, -9.42499, -9.75)))]
fillet_case("issue962_fillet_r0.8", p2, twelve, 0.8, "pass", 11552.7830)


# The same foot made small: an arm (a prism of a quadrilateral) fused to a
# block, the seam between their bottoms running from y=27.2 to y=10.2 -- the
# direction that put the second corner's point past the seam's end.
def arm_on_block():
    pts = [V(17.356, -9.425, -9.75), V(38.5, 10.2, -9.75), V(38.5, 27.2, -9.75),
           V(10.898, 16.471, -9.75)]
    arm = Part.Face(Part.makePolygon(list(reversed(pts)) + [pts[-1]])).extrude(V(0, 0, 6))
    bp = [V(38.5, 27.2, -9.75), V(38.5, 10.2, -9.75), V(50, 10.2, -9.75), V(50, 27.2, -9.75)]
    return Part.Face(Part.makePolygon(bp + [bp[0]])).extrude(V(0, 0, 6)).fuse(arm)


ab = arm_on_block()
foot = [edge_between(ab, (38.5, 27.2, -3.75), (38.5, 27.2, -9.75)),
        edge_between(ab, (38.5, 10.2, -3.75), (38.5, 10.2, -9.75))]
for r in (0.5, 0.8, 1.5):
    alone = sum(ab.Volume - ab.makeFillet(r, [ab.Edges[i - 1]]).Volume for i in foot)
    fillet_case("arm_on_block_foot_r%g" % r, ab, foot, r, "pass", ab.Volume - alone)

# A fin padded flush with the block's side: the outer face is two faces, the
# fin's and the block's, split at the floor's height, and the fillet's line on
# the fin's face ends on that split. It must take what it takes from the slot.
def fin_on_block(wall=False, arc=False):
    if arc:  # the fin's top an arc tangent to its outer face
        prof = Part.Wire([Part.LineSegment(V(0, 0, 5), V(10, 0, 5)).toShape(),
                          Part.LineSegment(V(10, 0, 5), V(10, 0, 8)).toShape(),
                          Part.Arc(V(10, 0, 8), V(4 + 6 * 0.5 ** 0.5, 0, 8 + 6 * 0.5 ** 0.5),
                                   V(4, 0, 14)).toShape(),
                          Part.LineSegment(V(4, 0, 14), V(0, 0, 14)).toShape(),
                          Part.LineSegment(V(0, 0, 14), V(0, 0, 5)).toShape()])
        f = Part.Face(prof).extrude(V(0, 3, 0))
    else:
        f = Part.makeBox(10, 3, 9, V(0, 0, 5))
    s = Part.makeBox(10, 10, 5).fuse(f)
    if wall:
        s = s.generalFuse([Part.LineSegment(V(9.3, 3, 5), V(9.3, 3, 20)).toShape()])[0].Solids[0]
    return s


for r in (0.3, 0.8, 2.0):
    taken = plain.Volume - plain.makeFillet(r, [plain.Edges[spine - 1]]).Volume
    fin = fin_on_block()
    fillet_case("fin_on_block_r%g" % r, fin, edge_between(fin, (10, 3, 5), (10, 3, 14)), r,
                "pass", fin.Volume - taken)
fin = fin_on_block(arc=True)
for r in (0.3, 0.8):
    fillet_case("fin_arc_on_block_r%g" % r, fin, edge_between(fin, (10, 3, 5), (10, 3, 8)), r,
                "pass")
# and its wall split as the slot's: a walking failure from 0.8 on
fin = fin_on_block(wall=True)
fillet_case("fin_on_block_wall_r0.3", fin, edge_between(fin, (10, 3, 5), (10, 3, 14)), 0.3,
            "pass", 770 - (plain.Volume - plain.makeFillet(0.3, [plain.Edges[spine - 1]]).Volume))
fillet_case("fin_on_block_wall_r0.8", fin, edge_between(fin, (10, 3, 5), (10, 3, 14)), 0.8,
            "xfail")

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
