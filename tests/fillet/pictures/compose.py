# python mkjobs|compose ... -- the render jobs, and one PNG per case from the panels.
#
#   python compose.py jobs            write $FILLET_WORK/jobs.json
#   python compose.py <outdir> [names...]
#
# Three columns -- upstream ($FILLET_WORK/r/up), the fork before the case's fix
# (r/<stage>), the fork now (r/after) -- and three rows: the whole result, the
# result zoomed on the fillet's end, and the (u, v) outline of the face at the
# fillet's end (UVFACE in cases.py), its open joints circled red. Needs Pillow (the FreeCAD conda env
# has it) and the DejaVu fonts (matplotlib's copy on macOS).
import glob
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cases import CASES, EDGES, NAMES, STAGES, UVLABEL, VIEW

P = os.environ["FILLET_WORK"]
R = {k: json.load(open(P + "/r/%s/results.json" % k)) for k in os.listdir(P + "/r")
     if os.path.exists(P + "/r/%s/results.json" % k)}
SIZE = (360, 300)

if sys.argv[1] == "jobs":
    import math
    panels = []
    os.makedirs(P + "/png", exist_ok=True)
    for name, (st, sk, ends, radius, ref, focus, look) in CASES.items():
        inp = P + "/r/after/%s.input.brep" % name
        for v, k in (("up", "up"), ("before", st), ("after", "after")):
            f = P + "/r/%s/%s.brep" % (k, name)
            ghost = not os.path.exists(f)
            # the whole shape, then the fillet's end
            for row, (center, height) in enumerate((VIEW[sk],
                                                    (focus, 4 * radius + 3))):
                panels.append(dict(brep=inp if ghost else f, ghost=ghost, mark=not ghost,
                                   eye=list(look), center=list(center), height=height,
                                   png=P + "/png/%s.%s.%d.png" % (name, v, row)))
    json.dump(dict(size=list(SIZE), panels=panels), open(P + "/jobs.json", "w"), indent=0)
    print(len(panels), "panels")
    sys.exit(0)

from PIL import Image, ImageDraw, ImageFont


def font(bold):
    n = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
    for d in ["/usr/share/fonts/truetype/dejavu/"] + glob.glob(
            sys.prefix + "/lib/python3*/site-packages/matplotlib/mpl-data/fonts/ttf/"):
        if os.path.exists(d + n):
            return d + n
    raise SystemExit("no " + n)


F = ImageFont.truetype(font(False), 15)
FB = ImageFont.truetype(font(True), 16)
FT = ImageFont.truetype(font(True), 19)


def status(d, name):
    if d["ok"]:
        return "valid, volume %.4f" % d["volume"], (20, 120, 40)
    s = "; ".join(x.replace(" BRep_API: command not done", "") for x in d["problems"])
    if d["volume"] is not None:
        s += "; volume %.4f" % d["volume"]
    return s, (190, 30, 30)


OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
for name in sys.argv[2:] or list(CASES):
    st, sk, ends, radius, ref, focus, look = CASES[name]
    before, fix = STAGES[st]
    W, H = SIZE
    top, lab = 64, 64
    UVH = 240
    img = Image.new("RGB", (3 * W, top + lab + 2 * H + UVH), "white")
    d = ImageDraw.Draw(img)
    d.text((10, 8), name, font=FT, fill=(0, 0, 0))
    if isinstance(ends, list):
        on = "%d edges (%s)" % (len(ends), EDGES.get(name, "see cases.py"))
    else:
        on = "the edge %s-%s" % (ends[0], ends[1])
    desc = "%s; fillet r%g on %s; expected %.4f  (fix: %s)" % (NAMES[sk], radius, on, ref, fix)
    words, line, lines = desc.split(" "), "", []
    for wd in words:
        if d.textlength(line + " " + wd, font=F) > 3 * W - 20 and line:
            lines.append(line)
            line = wd
        else:
            line = (line + " " + wd).strip()
    lines.append(line)
    for li, t in enumerate(lines[:2]):
        d.text((10, 30 + li * 17), t, font=F, fill=(60, 60, 60))
    cols = [("up", "up", "upstream OCCT 8.0.1"), ("before", st, "fork before (%s)" % before),
            ("after", "after", "fork after")]
    for ci, (v, k, title) in enumerate(cols):
        x = ci * W
        s, colr = status(R[k][name], name)
        d.text((x + 10, top), title, font=FB, fill=(0, 0, 0))
        words, line, lines = s.split(" "), "", []
        for wd in words:
            if d.textlength(line + " " + wd, font=F) > W - 20 and line:
                lines.append(line)
                line = wd
            else:
                line = (line + " " + wd).strip()
        lines.append(line)
        for li, t in enumerate(lines[:2]):
            d.text((x + 10, top + 21 + li * 19), t, font=F, fill=colr)
        for row in (0, 1):
            pf = P + "/png/%s.%s.%d.png" % (name, v, row)
            img.paste(Image.open(pf).convert("RGB"), (x, top + lab + row * H))
            meta = json.load(open(pf + ".json")) if os.path.exists(pf + ".json") else {}
            if row == 0 and R[k][name]["problems"] and R[k][name]["problems"][0].startswith("threw"):
                d.text((x + 10, top + lab + 8), "threw: the input is shown", font=F, fill=(190, 30, 30))
            if row == 1:
                d.text((x + 10, top + lab + H + 8), "the fillet's end", font=F, fill=(110, 110, 110))
            if meta.get("dropped") and row == 0:
                d.text((x + 10, top + lab + 28), "%d face(s) of no area, not drawn" % meta["dropped"],
                       font=F, fill=(190, 30, 30))
    # the third row: one (u, v) scale for the three columns
    uvs = [R[k][name].get("uv") or [] for _, k, _ in cols]
    pts = [p for uv in uvs for wire in uv for c in wire for p in c]
    if pts:
        u0 = min(p[0] for p in pts); u1 = max(p[0] for p in pts)
        v0 = min(p[1] for p in pts); v1 = max(p[1] for p in pts)
        y0 = top + lab + 2 * H
        mx, my = 34, 40
        sx = (W - 2 * mx) / max(u1 - u0, 1e-9); sy = (UVH - my - 22) / max(v1 - v0, 1e-9)
        for ci, uv in enumerate(uvs):
            x = ci * W
            P2 = lambda p: (x + mx + (p[0] - u0) * sx, y0 + my + (v1 - p[1]) * sy)
            label, grid = UVLABEL[sk]
            d.text((x + 10, y0 + 6), label + " in (u, v)", font=F, fill=(110, 110, 110))
            k = 0
            while grid and k * (math.pi / 2) <= u1 + 1e-9:
                if k * (math.pi / 2) >= u0 - 1e-9:
                    gx = P2((k * math.pi / 2, v0))[0]
                    d.line([(gx, y0 + my - 4), (gx, y0 + UVH - 20)], fill=(225, 225, 225), width=1)
                    lbl = ["0", "pi/2", "pi", "3pi/2", "2pi", "5pi/2", "3pi"][k] if k < 7 else ""
                    d.text((gx - d.textlength(lbl, font=F) / 2, y0 + UVH - 19), lbl, font=F,
                           fill=(140, 140, 140))
                k += 1
            ends = [c[0] for wire in uv for c in wire] + [c[-1] for wire in uv for c in wire]
            for wire in uv:
                for c in wire:
                    d.line([P2(p) for p in c], fill=(40, 40, 60), width=2)
            for e in ends:
                if sum(1 for o in ends if abs(o[0] - e[0]) < 1e-4 and abs(o[1] - e[1]) < 1e-4) < 2:
                    cx, cy = P2(e)
                    d.ellipse([cx - 7, cy - 7, cx + 7, cy + 7], outline=(210, 30, 30), width=3)
    for ci in (1, 2):
        d.line([(ci * W, top), (ci * W, top + lab + 2 * H + UVH)], fill=(200, 200, 200), width=1)
    d.line([(0, top + lab + H), (3 * W, top + lab + H)], fill=(225, 225, 225), width=1)
    d.line([(0, top + lab + 2 * H), (3 * W, top + lab + 2 * H)], fill=(225, 225, 225), width=1)
    img = img.quantize(colors=128, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE)
    img.save(os.path.join(OUT, name + ".png"), optimize=True)
    print(name, os.path.getsize(os.path.join(OUT, name + ".png")))
