# FreeCADCmd compute.py -- run the cases and keep what came out.
#
# env: FILLET_TOOLS = this directory, OUT = where the results go (a .brep per
# case and its input, results.json), STAGES = "s523 ..." or "all".
import json
import os
import sys

sys.path.insert(0, os.environ["FILLET_TOOLS"])
import FreeCAD as App
import Part
from cases import CASES, SHAPES, UVFACE, edge_between

OUT = os.environ["OUT"]
os.makedirs(OUT, exist_ok=True)
# Unfrozen, as the suite runs: a frozen input would hide an edit of it.
App.ParamGet("User parameter:BaseApp/Preferences/Mod/Part").SetBool("ImmutableShapeValues", False)
stages = os.environ.get("STAGES", "all").split()


def emit(s):
    os.write(1, (s + "\n").encode())


def uv_outline(shape, pick):
    """The pcurves of the face `pick` chooses, each edge as it occurs in its
    wire (what BRepCheck reads), sampled: [[[u, v], ...], ...] per wire."""
    out = []
    for f in shape.Faces:
        if not pick(f):
            continue
        for w in f.childShapes():
            wire = []
            for e in w.childShapes():
                if e.ShapeType != "Edge":
                    continue
                c, a, b = f.curveOnSurface(e)[:3]
                pts = [c.value(a + (b - a) * i / 24.0) for i in range(25)]
                if e.Orientation == "Reversed":
                    pts.reverse()
                wire.append([[p.x, p.y] for p in pts])
            out.append(wire)
    return out


res = {}
for name, (st, sk, ends, radius, ref, focus, look) in CASES.items():
    if "all" not in stages and st not in stages:
        continue
    s = SHAPES[sk]()
    s.exportBrep(os.path.join(OUT, name + ".input.brep"))
    before = (s.isValid(), round(s.Volume, 6))
    d = {"input_volume": s.Volume, "ref": ref}
    try:
        r = s.makeFillet(radius, [edge_between(s, *ends)])
        p = []
        if (s.isValid(), round(s.Volume, 6)) != before:
            p.append("input changed")
        if not r.isValid():
            p.append("invalid")
        if len(r.Shells) != 1:
            p.append("shells=%d" % len(r.Shells))
        elif not r.Shells[0].isClosed():
            p.append("open shell")
        vol = r.Volume
        if not p and ref is not None and abs(vol - ref) > 1e-4 * abs(ref):
            p.append("volume %.4f, expected %.4f" % (vol, ref))
        r.exportBrep(os.path.join(OUT, name + ".brep"))
        d["uv"] = uv_outline(r, UVFACE[sk])
        d.update(ok=not p, problems=p, volume=vol)
    except Exception as e:
        d.update(ok=False, problems=["threw " + str(e).strip().splitlines()[-1]], volume=None)
    res[name] = d
    emit("%-28s %s %s" % (name, "OK " if d["ok"] else "BAD", "; ".join(d["problems"])
                          or "vol=%.4f" % d["volume"]))
json.dump(res, open(os.path.join(OUT, "results.json"), "w"), indent=1)
emit("DONE")
os._exit(0)
