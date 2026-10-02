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
    """makeFillet(radius, [Edge<edge_index>]).

    expect='pass':  a valid solid, one closed shell, the input left as it was,
                    and ref_volume when given.
    expect='xfail': known broken (see README.md); an exception, an invalid
                    shape or an open shell counts as the expected failure.
    """
    try:
        before = signature(shape)
        r = shape.makeFillet(radius, [shape.Edges[edge_index - 1]])
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
