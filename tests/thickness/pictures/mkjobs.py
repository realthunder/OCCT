# python3 mkjobs.py [names...] -- the render jobs for every case (or those named).
#
# Four columns -- the shape given with its removed faces marked, upstream
# ($THICK_WORK/r/up), the fork before the case's fix (r/<stage>), the fork now
# (r/after) -- and two rows, the whole and cut open. One camera per case,
# framed on the input and every result of a sane size, so the columns compare.
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cases import CASES
W = os.environ["THICK_WORK"]
R = {k: json.load(open(W + "/r/%s/results.json" % k)) for k in os.listdir(W + "/r")
     if os.path.exists(W + "/r/%s/results.json" % k)}
names = sys.argv[1:] or list(CASES)
panels = []
os.makedirs(W + "/png", exist_ok=True)
for name in names:
    stage = CASES[name][0]
    var = [("up", "up"), ("before", stage), ("after", "after")]
    after = R["after"][name]
    inp = W + "/r/after/%s.input.brep" % name
    breps = [inp]
    for _, k in var:
        f = W + "/r/%s/%s.brep" % (k, name)
        if os.path.exists(f) and R[k][name]["volume"] and abs(R[k][name]["volume"]) < 1e5:
            breps.append(f)
    fi = CASES[name][2]
    removed = [i - 1 for i in (fi if isinstance(fi, list) else [fi])]
    for cut in (0, 1):
        panels.append(dict(brep=inp, ghost=False, input=True, removed=removed, mark=False, cut=cut,
                           normal=after["normal"], center=after["center"], bbox_breps=breps,
                           cutaxis=[0, 0, 1] if "_fillet_" in name else None,
                           png=W + "/png/%s.in.%d.png" % (name, cut)))
    for v, k in var:
        d = R[k][name]
        f = W + "/r/%s/%s.brep" % (k, name)
        ghost = not os.path.exists(f) or bool(d["problems"]) and d["problems"][0].startswith("threw")
        for cut in (0, 1):
            panels.append(dict(brep=inp if ghost else f, ghost=ghost, mark=not ghost, cut=cut,
                               normal=after["normal"], center=after["center"], bbox_breps=breps,
                               # a fillet's opening shows best in plan
                               cutaxis=[0, 0, 1] if "_fillet_" in name else None,
                               png=W + "/png/%s.%s.%d.png" % (name, v, cut)))
json.dump(dict(size=[360, 300], panels=panels), open(W + "/jobs.json", "w"), indent=0)
print(len(panels), "panels")
