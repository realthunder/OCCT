# python compose.py <outdir> [names...] -- one PNG per case from the panels.
#
# Needs Pillow (the FreeCAD conda env has it) and the DejaVu fonts.
import json, os, re, sys
from PIL import Image, ImageDraw, ImageFont
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cases import CASES, INPUT, REFUSED, STAGES as STAGE
P = os.environ["THICK_WORK"]
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
R = {k: json.load(open(P + "/r/%s/results.json" % k)) for k in os.listdir(P + "/r")
     if os.path.exists(P + "/r/%s/results.json" % k)}
FD = "/usr/share/fonts/truetype/dejavu/"
F = ImageFont.truetype(FD + "DejaVuSans.ttf", 15)
FB = ImageFont.truetype(FD + "DejaVuSans-Bold.ttf", 16)
FT = ImageFont.truetype(FD + "DejaVuSans-Bold.ttf", 19)
NAMES = {"cyl": "cylinder", "ann": "cylinder with a through hole", "ell": "elliptic pad",
         "cylpocket": "cylinder with a blind pocket", "boxhole": "box with a through hole", "lbox": "L-box",
         "tshape": "T", "pocketbox": "box with a 6x6 pocket", "conehole": "cone with a through hole",
         "pocket5": "box with a 5x5 pocket", "blindhole": "box with a blind hole",
         "filletbox": "box with filleted vertical edges, r2", "filletbox25": "box with filleted vertical edges, r2.5",
         "shortcyl": "cylinder 1.5 high", "shortbox": "box 10 x 10 x 1.5", "box6": "box 10 x 10 x 6",
         "blindpot": "a blind hole's wall and floor and the top, an open shell",
         "cylboss": "cylinder with a boss", "sector": "ring sector pad across angle 0", "sphere": "sphere"}
JOIN = {0: "Arc join", 2: "Intersection join"}
def status(d):
    if d["ok"] and d.get("refused"):
        return "refused", (20, 120, 40)
    if d["ok"]:
        if d.get("solids", 1) > 1:
            sealed = " (%d solids)" % d["solids"]
        else:
            sealed = " (a skin and a void)" if d.get("shells", 1) == 2 else ""
        note = " (%s)" % d["note"] if d.get("note") else ""
        return "valid, volume %.2f%s%s" % (d["volume"], sealed, note), (20, 120, 40)
    p = []
    for x in d["problems"]:
        x = re.sub(r"threw \d+", "threw ", x)
        x = x.replace(" BRep_API: command not done", "").replace("BRep_Tool:: ", "")
        p.append(x)
    s = "; ".join(p)
    if d["volume"] is not None and d["volume"] < 0:
        p.insert(0, "inside out")
        s = "; ".join(p)
    if d["volume"] is not None and not any(x.startswith("volume") for x in d["problems"]):
        s += "; volume %.2f" % d["volume"]
    iv = R["after"][name].get("input_volume")
    if d["volume"] is not None and iv and abs(d["volume"] - iv) < 1e-3 * iv:
        s += " (unhollowed)"
    return s, (190, 30, 30)
names = sys.argv[2:] or list(CASES)
for name in names:
    st, sk, fi, val, inter, join, ref = CASES[name]
    commit, sec = STAGE[st]
    W, H = 360, 300
    top = 64; lab = 64
    img = Image.new("RGB", (3 * W, top + lab + 2 * H), "white")
    d = ImageDraw.Draw(img)
    d.text((10, 8), name, font=FT, fill=(0, 0, 0))
    nf = 6 if sk.startswith(("box", "short")) and "cyl" not in sk else 3
    if not isinstance(fi, list):
        faces = "Face%d" % fi
    elif len(fi) == nf - 1:
        faces = "all but Face%d" % [i for i in range(1, nf + 1) if i not in fi][0]
    else:
        faces = "Faces " + ", ".join(str(i) for i in fi)
    if name in REFUSED:
        expect = "a refusal, no face stays"
    elif name in INPUT:
        expect = "the input unchanged, %.2f" % ref
    else:
        expect = "%.2f" % ref if ref else "a valid solid"
    desc = "%s, %s removed, %s, %s%s; expected %s  (fix: sec %s)" % (
        NAMES[sk], faces, "outward +1" if val > 0 else "inward -1", JOIN[join], ", intersection on" if inter else "",
        expect, sec)
    d.text((10, 36), desc, font=F, fill=(60, 60, 60))
    cols = [("up", "up", "upstream OCCT 8.0.1"), ("before", st, "fork before (%s)" % commit), ("after", "after", "fork after")]
    for ci, (v, k, title) in enumerate(cols):
        x = ci * W
        s, colr = status(R[k][name])
        d.text((x + 10, top), title, font=FB, fill=(0, 0, 0))
        # wrap the status to the column
        words = s.split(" "); line = ""; lines = []
        for wd in words:
            if d.textlength(line + " " + wd, font=F) > W - 20 and line:
                lines.append(line); line = wd
            else:
                line = (line + " " + wd).strip()
        lines.append(line)
        d.text((x + 10, top + 21), lines[0], font=F, fill=colr)
        for cut in (0, 1):
            pf = P + "/png/%s.%s.%d.png" % (name, v, cut)
            im = Image.open(pf).convert("RGB")
            img.paste(im, (x, top + lab + cut * H))
            if R[k][name]["problems"] and R[k][name]["problems"][0].startswith("threw") and cut == 0:
                d.text((x + 10, top + lab + 8), "threw: the input is shown", font=F, fill=(190, 30, 30))
            if R[k][name].get("refused") and cut == 0:
                d.text((x + 10, top + lab + 8), "refused: the input is shown", font=F, fill=(110, 110, 110))
            meta = json.load(open(pf + ".json")) if os.path.exists(pf + ".json") else {}
            if meta.get("huge") and cut == 0:
                d.text((x + 10, top + lab + 8), "own scale: it spans %.0e mm" % meta["huge"], font=F, fill=(190, 30, 30))
            if meta.get("dropped") and cut == 0:
                d.text((x + 10, top + lab + 8), "%d face(s) of no area, not drawn" % meta["dropped"], font=F, fill=(190, 30, 30))
            if meta.get("clipped") and cut == 1:
                d.text((x + 10, top + lab + H + 8), "clipped (the cut failed or is inside out)", font=F, fill=(110, 110, 110))
        if len(lines) > 1:
            d.text((x + 10, top + 40), " ".join(lines[1:]), font=F, fill=colr)
    for ci in (1, 2):
        d.line([(ci * W, top), (ci * W, top + lab + 2 * H)], fill=(200, 200, 200), width=1)
    d.line([(0, top + lab + H), (3 * W, top + lab + H)], fill=(225, 225, 225), width=1)
    img = img.quantize(colors=128, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE)
    img.save(os.path.join(OUT, name + ".png"), optimize=True)
    print(name, os.path.getsize(os.path.join(OUT, name + ".png")))
