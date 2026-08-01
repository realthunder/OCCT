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


# ---------------------------------------------------------------------------
# Programmatic cases: thickness of simple solids, removing one face at a time.
# ---------------------------------------------------------------------------

def thickness_case(name, shape, face_index, value, expect, ref_volume=None):
    """makeThickness(mode=Skin, join=Arc) removing 1-based face `face_index`.

    expect='pass':  result must be a valid solid, one closed shell, and match
                    ref_volume (captured from a verified-good run).
    expect='xfail': known broken (see README.md); any exception, invalid
                    shape, or open shell counts as the expected failure.
    """
    detail = ""
    try:
        faces = [shape.Faces[face_index - 1]]
        r = shape.makeThickness(faces, value, 1e-7, False, False, 0, 0)
        problems = []
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


document_case("issue1_ellipse_thickness", "issue1_ellipse_thickness.FCStd",
              volumes={"Thickness": 3598.2930})
document_case("issue2_broken_loft", "issue2_broken_loft.FCStd",
              volumes={"AdditiveLoft": 6143.8106})
document_case("issue3_pad_thickness", "issue3_pad_thickness.FCStd",
              volumes={"Thickness": 1241.0718})
document_case("issue4_revolution_thickness", "issue4_revolution_thickness.FCStd",
              volumes={"Thickness": 431.4454})

# Plain cylinder: Face1 = lateral (seam), Face2 = bottom, Face3 = top.
cyl = Part.makeCylinder(4, 20)
thickness_case("cyl_side_out",    cyl, 1, +1.0, "xfail")
thickness_case("cyl_side_in",     cyl, 1, -1.0, "xfail")
thickness_case("cyl_bottom_out",  cyl, 2, +1.0, "pass", 637.5858)
thickness_case("cyl_bottom_in",   cyl, 2, -1.0, "xfail")
thickness_case("cyl_top_out",     cyl, 3, +1.0, "pass", 637.5858)
thickness_case("cyl_top_in",      cyl, 3, -1.0, "pass", 468.0973)

# Cylinder with a centered hole: Face1 = outer lateral, Face2 = bottom,
# Face3 = top, Face4 = hole lateral.
ann = Part.makeCylinder(5, 10).cut(Part.makeCylinder(2, 10))
thickness_case("hole_outer_out",  ann, 1, +1.0, "xfail")
thickness_case("hole_outer_in",   ann, 1, -1.0, "xfail")
thickness_case("hole_bottom_out", ann, 2, +1.0, "pass", 540.3400)
thickness_case("hole_bottom_in",  ann, 2, -1.0, "xfail")
thickness_case("hole_top_out",    ann, 3, +1.0, "pass", 540.3400)
thickness_case("hole_top_in",     ann, 3, -1.0, "pass", 461.8141)
thickness_case("hole_inner_out",  ann, 4, +1.0, "xfail")
thickness_case("hole_inner_in",   ann, 4, -1.0, "xfail")

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
