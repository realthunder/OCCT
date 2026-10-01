# local02: booleans of two solids that share a face TShape, see ../README.md.
# Run: FreeCADCmd local02_fuse_shared_face.py   (prints one line per op)
import os
import Part

HERE = os.path.dirname(os.path.abspath(__file__))
# The right answers, from the same pair with the tool copied (nothing shared)
TRUTH = {"fuse": 9221.776, "common": 501.754, "base-tool": 8218.270, "tool-base": 501.752}


def run(path):
    base, tool = Part.read(path).childShapes()
    bad = 0
    for label, b, t, op in (("fuse", base, tool, "fuse"), ("common", base, tool, "common"),
                            ("base-tool", base, tool, "cut"), ("tool-base", tool, base, "cut")):
        r = getattr(b, op)(t)
        ok = abs(r.Volume - TRUTH[label]) < 1e-2
        bad += not ok
        print("%s %-9s solids=%d vol=%9.3f %s" % (os.path.basename(path), label, len(r.Solids),
              r.Volume, "ok" if ok else "WRONG (want %.3f)" % TRUTH[label]))
    return bad


bad = sum(run(os.path.join(HERE, "local02_fuse_shared_face_%s.brep" % k)) for k in ("fail", "ok"))
print("local02: %d wrong" % bad)
