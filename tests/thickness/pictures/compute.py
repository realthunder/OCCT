# FreeCADCmd compute.py -- run the cases and keep what came out.
#
# env: THICK_TOOLS = this directory, OUT = where the results go (a .brep per
# case and its input, results.json), STAGES = "s89 s90 ..." or "all", ONLY =
# case names, DOCS = 1 to run the captured models too.
import os, sys, json
sys.path.insert(0, os.environ["THICK_TOOLS"])
import FreeCAD as App, Part
from cases import SHAPES, CASES, DOCS, SHELLS
OUT = os.environ["OUT"]; os.makedirs(OUT, exist_ok=True)
stages = os.environ.get("STAGES", "all").split()
only = os.environ.get("ONLY", "").split()
def emit(s): os.write(1, (s + "\n").encode())
def judge(r, ref, shells=1):
    p = []
    if r.ShapeType != "Solid": p.append("type=%s" % r.ShapeType)
    try:
        if not r.isValid(): p.append("invalid")
    except Exception: p.append("check threw")
    if len(r.Shells) != shells: p.append("shells=%d" % len(r.Shells))
    elif not all(sh.isClosed() for sh in r.Shells): p.append("open shell")
    try: vol = r.Volume
    except Exception: vol = None
    if not p and ref is not None and vol is not None and abs(vol - ref) > 1e-3 * abs(ref):
        p.append("volume %.2f, expected %.2f" % (vol, ref))
    return p, vol
def face_info(f):
    u0, u1, v0, v1 = f.ParameterRange
    c = f.CenterOfMass if f.Area > 1e-9 else f.valueAt((u0 + u1) / 2, (v0 + v1) / 2)
    uv = f.Surface.parameter(c)
    n = f.normalAt(*uv)
    return [c.x, c.y, c.z], [n.x, n.y, n.z]
res = {}
for name, (st, sk, fi, val, inter, join, ref) in CASES.items():
    if "all" not in stages and st not in stages: continue
    if only and name not in only: continue
    s = SHAPES[sk]()
    f = s.Faces[fi - 1]
    c, n = face_info(f)
    s.exportBrep(os.path.join(OUT, name + ".input.brep"))
    d = {"center": c, "normal": n, "value": val, "ref": ref, "input_volume": s.Volume}
    try:
        r = s.makeThickness([f], val, 1e-7, inter, False, 0, join)
        p, vol = judge(r, ref, SHELLS.get(name, 1))
        r.exportBrep(os.path.join(OUT, name + ".brep"))
        d.update(ok=not p, problems=p, volume=vol, shells=len(r.Shells))
    except Exception as e:
        d.update(ok=False, problems=["threw " + (str(e).strip().splitlines() or [type(e).__name__])[-1]], volume=None)
    res[name] = d
    emit("%-32s %s %s" % (name, "OK " if d["ok"] else "BAD", "; ".join(d["problems"]) or "vol=%.4f" % d["volume"]))
if os.environ.get("DOCS"):
    M = os.path.join(os.environ["THICK_TOOLS"], "..", "models")
    for name, (fn, on, ref) in DOCS.items():
        if only and name not in only: continue
        doc = App.openDocument(os.path.join(M, fn))
        for o in doc.Objects: o.touch()
        doc.recompute()
        o = doc.getObject(on)
        if hasattr(o, "Base"):
            base, subs = o.Base
            val = o.Value.Value * (-1 if o.Reversed else 1)
        else:
            base, subs = o.Faces
            val = o.Value.Value
        bs = base.Shape
        bs.exportBrep(os.path.join(OUT, name + ".input.brep"))
        c, n = face_info(bs.getElement(subs[0]))
        d = {"center": c, "normal": n, "value": val, "ref": ref}
        st = [str(x) for x in o.State]
        if "Invalid" in st or "Error" in st or o.Shape.isNull():
            d.update(ok=False, problems=["failed: " + (" ".join(st))], volume=None)
            # PartDesign keeps the previous shape on error: do not save it
        else:
            p, vol = judge(o.Shape.Solids[0] if len(o.Shape.Solids) == 1 else o.Shape, ref)
            o.Shape.exportBrep(os.path.join(OUT, name + ".brep"))
            d.update(ok=not p, problems=p, volume=vol)
        res[name] = d
        emit("%-32s %s %s" % (name, "OK " if d["ok"] else "BAD", "; ".join(d["problems"]) or "vol=%.4f" % d["volume"]))
        App.closeDocument(doc.Name)
with open(os.path.join(OUT, "results.json"), "w") as fp: json.dump(res, fp, indent=1)
emit("DONE")
os._exit(0)
