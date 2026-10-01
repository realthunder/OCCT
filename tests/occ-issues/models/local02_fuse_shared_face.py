# local02: booleans of two solids that share a face TShape, see ../README.md.
# Run: FreeCADCmd local02_fuse_shared_face.py   (prints one line per op)
# FreeCAD's booleans run non-destructive, so each op must leave its inputs
# as they were; the tolerances of every vertex and edge are compared too.
import os
import Part

HERE = os.path.dirname(os.path.abspath(__file__))
# The right answers, from the same pair with the tool copied (nothing shared)
TRUTH = {"fuse": 9221.776, "common": 501.754, "base-tool": 8218.270, "tool-base": 501.752}


def tolerances(shapes):
    # each vertex and edge once, though the two shapes share some
    res = []
    for shape in shapes:
        for s in shape.Vertexes + shape.Edges:
            if not any(s.isSame(x) for x, _ in res):
                res.append((s, s.Tolerance))
    return res


def run(path):
    bad = 0
    for label, op, swap in (("fuse", "fuse", False), ("common", "common", False),
                            ("base-tool", "cut", False), ("tool-base", "cut", True)):
        base, tool = Part.read(path).childShapes()
        b, t = (tool, base) if swap else (base, tool)
        before = tolerances([b, t])
        r = getattr(b, op)(t)
        touched = sum(1 for s, tol in before if s.Tolerance != tol)
        ok = abs(r.Volume - TRUTH[label]) < 1e-2 and not touched
        bad += not ok
        print("%s %-9s solids=%d vol=%9.3f inputs touched=%d %s" % (
            os.path.basename(path), label, len(r.Solids), r.Volume, touched,
            "ok" if ok else "WRONG (want %.3f, inputs untouched)" % TRUTH[label]))
    return bad


bad = sum(run(os.path.join(HERE, "local02_fuse_shared_face_%s.brep" % k)) for k in ("fail", "ok"))
print("local02: %d wrong" % bad)
