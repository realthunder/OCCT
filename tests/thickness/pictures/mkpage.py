#!/usr/bin/env python3
"""mkpage.py OUTDIR -- the "Thickness Before and After" page from the pictures.

Writes OUTDIR/index.html and copies every picture of models/pictures that
cases.py names into OUTDIR/img/. One section per fix of models/Thickness.md
that has pictures, oldest first, its title taken from the "### Sec N: ..."
heading and its summary from SUMMARY below. With THICK_WORK set (the
make_pictures.sh work directory) the counts of what upstream gets right come
from its results; otherwise they are left out.

Run it after every picture run, then publish OUTDIR (index.html plus img/).

env: SUITE  the suite's standing, e.g. "277 passing, 8 marked known broken"
     SWEEP  the sweep's, e.g. "928 of 928"
"""
import html, json, os, re, shutil, struct, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cases  # noqa: E402

DOC = os.path.join(HERE, "..", "models", "Thickness.md")
PICS = os.path.join(HERE, "..", "models", "pictures")

# One or two sentences per section: what was wrong and what the fix does.
SUMMARY = {
    2: "Upstream gets every case here right; the fork had broken all of them. Four causes in the loop and the context extension: a closed edge found twice on a periodic face (the floor went missing), one seam wire allowed per face (the second wall built from two bare circles), a const edge losing its orientation, and the offset of an ellipse stretched 100 lengths past its ends.",
    3: "Outward with the Arc join past a concave corner, the corner piece of an arc face came out as a face of its own. Plus a neighbour skipped before its edge was renewed, and a stretched edge lying on another (the T's bar top).",
    4: "The stretched removed face's crossing taken for an edge's own end (T bar top inward); an edge two faces share trimmed twice (the fork threw where upstream is right); a wire running once round a periodic face (a face of no area on the holed cone).",
    5: "Upstream fails all of these. Thick solids that came back inside out (drawn black) are now oriented by classification; a blind floor's section was the wrong way round; concave removed faces with intersection on are built as with it off.",
    6: "A box with filleted vertical edges, an end or side face removed: the fillets' offsets never meet the removed face (unhollowed) or meet it far round (a lip). A tube round the tangent edge closes the gap, with sphere eighths at the corners outward.",
    8: "The removed face is curved. Outward it threw, inward it came back valid but larger than the input. The tubes are built on the fillet's cylinder, the loop walks its angles as on a chart, and a hanging piece of the floor's offset is dropped. The cut is in plan.",
    9: "That join builds no tubes. The gap is closed by their sharp counterpart: a strip of the neighbour's tangent plane a thickness into the removed face, and a wall square to it; cubes at the corners.",
    10: "The holed cone's top removed inward: the cavity closes below the removed face, and the right result is a skin and a sealed void, two shells.",
    11: "Removing a face can leave the kept faces in separate pieces: each piece is thickened on its own and the result is a compound. A short box or cylinder taken down to one face is answered too.",
    12: "The same pieces with the Intersection join: a seam edge given once to the loop left the whole cap; one face left beside the removed ones threw.",
    13: "A top face with a pocket or blind hole removed leaves the pocket's wall and floor as a piece of its own; the fork built it wrong even as an open shell, where upstream is right.",
    14: "The thickness was right and left its input inside out: a test face was healed on the input's own edges. The panels show the input after the call.",
    15: "With every face removed no face stays to be thickened: refused. A sphere used to come back as itself, valid and unhollowed. The panels show the input.",
    16: "A dome, a bowl, a cap, a cone with its apex: the seam band's far edge can be the pole's degenerated edge. Upstream is right on these; the fork had been wrong since the chain was ported.",
    17: "An edge re-trimmed at an offset cone's apex had an infinite range; half a dome threw rebinding a border vertex. A result that is the input itself is refused.",
    18: "Half a dome's mirror side; the Intersection join on faces in coplanar pieces (a box fused of two); the closing wall on a sphere; and the offset sphere grown round its pole by turning the sphere's axis off the face.",
    19: "Three quarters of a dome (the axis turn stopped at half a turn; the old pole's edge kept two meridians apart; a meridian's piece past the pole kept a pcurve off the sphere); the piece of sphere where three tubes meet at a pole; a face extended to the far crossing; and the sphere itself removed, its wall built on a turned twin of the removed face. At twice the thickness one case had been a valid solid on the wrong side.",
    20: "Half a ball cut through both its poles: its sphere is cut in two first, along the equator where it stays and along a meridian where it is removed. And shapes under a location or turned in space, which gave other results than plain: a seam checked without its location, a cone's apex left where it was, a corner arc taken from the wrong section, a whole circle of section that lost the piece its own vertex lies in.",
    21: "A dome whose rim is in two arcs was broken since the chain was ported, upstream right: the face is now walked in (u, v) with its seam, and its flat's two arcs no longer both become one half circle. The half ball with one disc, as a refine or a cut leaves it, is put on a sphere turned onto its middle, which makes it that dome; and a half ball with its sphere in two faces has them joined first.",
}

TAGLINE = "Every pictured case of the OCCT fork's thickness suite, each run three ways"


def sections():
    """{number: title} from the write-up's headings."""
    out = {}
    for line in open(DOC, encoding="utf-8"):
        m = re.match(r"### Sec (\d+): (.*)", line)
        if m:
            t = m.group(2).strip()
            out[int(m.group(1))] = t[:1].upper() + t[1:]
    return out


def png_size(path):
    with open(path, "rb") as f:
        head = f.read(24)
    return struct.unpack(">II", head[16:24])


def main(out):
    titles = sections()
    work = os.environ.get("THICK_WORK")
    up = None
    if work and os.path.exists(os.path.join(work, "r", "up", "results.json")):
        up = json.load(open(os.path.join(work, "r", "up", "results.json")))
    # section -> [(stage, [case, ...]), ...] in the order the stages are listed
    by_sec = {}
    for st, (commit, sec) in cases.STAGES.items():
        names = [n for n, c in cases.CASES.items() if c[0] == st]
        if names:
            by_sec.setdefault(int(sec), []).append((st, commit, names))
    order = sorted(by_sec)
    os.makedirs(os.path.join(out, "img"), exist_ok=True)

    total = len(cases.CASES)
    up_right = sum(1 for n in cases.CASES if up and up.get(n, {}).get("ok"))
    newest = order[-1]
    intro = "%s: %d cases across %d fixes, the newest last (sec %d)." % (TAGLINE, total, len(order), newest)
    if up:
        intro += " Upstream gets %d of them right (the fix chain's own regressions) and fails %d." % (
            up_right, total - up_right)
    if os.environ.get("SUITE"):
        intro += " The suite stands at %s." % os.environ["SUITE"]
    if os.environ.get("SWEEP"):
        intro += " The sweep: %s." % os.environ["SWEEP"]
    intro += (" Section numbers are those of the write-up, which lives in the fork with the pictures:"
              " <code>tests/thickness/models/Thickness.md</code> and <code>models/pictures/</code>.")

    e = html.escape
    parts = []
    for sec in order:
        figs, metas = [], []
        for st, commit, names in by_sec[sec]:
            ok = sum(1 for n in names if up and up.get(n, {}).get("ok"))
            meta = "Before = fork at <code>%s</code> &middot; %d case%s" % (
                e(commit), len(names), "" if len(names) == 1 else "s")
            if up:
                meta += ", upstream right on %d" % ok
            meta += ": " + " ".join('<a href="#%s">%s</a>' % (e(n), e(n)) for n in names)
            metas.append('<p class="meta">%s</p>' % meta)
            for n in names:
                src = os.path.join(PICS, n + ".png")
                shutil.copyfile(src, os.path.join(out, "img", n + ".png"))
                w, h = png_size(src)
                figs.append(
                    '<figure id="%s"><img src="img/%s.png" alt="%s: the shape given, upstream, fork before, fork after"'
                    ' loading="lazy" width="%d" height="%d"><figcaption><code>%s</code></figcaption></figure>'
                    % (e(n), e(n), e(n), w, h, e(n)))
        parts.append(
            '<section id="sec%d"><header class="sh"><span class="tag">sec %d</span><h2>%s</h2></header>\n'
            '<p>%s</p>%s\n%s</section>'
            % (sec, sec, e(titles.get(sec, "")), e(SUMMARY.get(sec, "")), "".join(metas), "\n".join(figs)))
    nav = " ".join('<a href="#sec%d">%d</a>' % (s, s) for s in order)
    page = TEMPLATE.replace("@INTRO@", intro).replace("@NAV@", nav).replace("@SECTIONS@", "\n".join(parts))
    with open(os.path.join(out, "index.html"), "w", encoding="utf-8") as f:
        f.write(page)
    print("%d sections, %d cases -> %s" % (len(order), total, out))


TEMPLATE = """<title>Thickness Before and After</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
/* A long scroll of figure plates under one sticky index of fixes. */
:root{--bg:#f3f5f8;--panel:#ffffff;--ink:#1d2430;--muted:#5a6575;--line:#d8dee7;--accent:#3d5f8f;--good:#1f7a3a;--bad:#b3261e;--cut:#c66a2a;--plate:#ffffff;--plate-line:#e3e7ee;--plate-cap:#fafbfc;--plate-ink:#5a6575;color-scheme:light}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#141922;--panel:#1c2330;--ink:#e3e8ef;--muted:#9aa6b6;--line:#2d3746;--accent:#8fb0dd;--good:#6cc58a;--bad:#ef8a80;--cut:#e39a62;--plate:#ffffff;--plate-line:#e3e7ee;--plate-cap:#fafbfc;--plate-ink:#5a6575;color-scheme:dark}}
:root[data-theme="dark"]{--bg:#141922;--panel:#1c2330;--ink:#e3e8ef;--muted:#9aa6b6;--line:#2d3746;--accent:#8fb0dd;--good:#6cc58a;--bad:#ef8a80;--cut:#e39a62;--plate:#ffffff;--plate-line:#e3e7ee;--plate-cap:#fafbfc;--plate-ink:#5a6575;color-scheme:dark}
body{background:var(--bg);color:var(--ink);font:15px/1.55 "IBM Plex Sans",system-ui,sans-serif;padding-inline:16px;padding-block:24px 64px}
main{max-width:1120px;margin:0 auto;display:grid;gap:40px}
h1{font-size:28px;line-height:1.2;margin:0 0 8px;text-wrap:balance}
h2{font-size:20px;margin:0;text-wrap:balance}
p{max-width:68ch;margin:0}
code{font-family:"IBM Plex Mono",ui-monospace,monospace;font-size:.88em}
a{color:var(--accent)}
a:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
.intro{display:grid;gap:12px}
.key{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:8px 24px;padding:14px 16px;border:1px solid var(--line);border-radius:6px;background:var(--panel);font-size:14px}
.key b{font-weight:600}
.sw{display:inline-block;width:.8em;height:.8em;border-radius:2px;vertical-align:-1px;margin-right:6px}
nav{display:flex;flex-wrap:wrap;gap:6px 14px;font-family:"IBM Plex Mono",monospace;font-size:14px;position:sticky;top:env(safe-area-inset-top,0px);background:var(--bg);padding-block:8px;z-index:1;border-bottom:1px solid var(--line)}
section{display:grid;gap:14px;min-width:0}
.sh{display:flex;align-items:baseline;gap:12px;flex-wrap:wrap}
.tag{font-family:"IBM Plex Mono",monospace;font-size:13px;letter-spacing:.04em;color:var(--accent);text-transform:uppercase}
.meta{color:var(--muted);font-size:13.5px;max-width:none}
.meta a{font-family:"IBM Plex Mono",monospace;font-size:12.5px;margin-right:4px}
figure{margin:0;background:var(--plate);border:1px solid var(--line);border-radius:6px;overflow:hidden}
figure img{display:block;width:100%;height:auto}
figcaption{padding:6px 12px;font-size:13px;color:var(--plate-ink);border-top:1px solid var(--plate-line);background:var(--plate-cap)}
@media (prefers-reduced-motion:no-preference){html{scroll-behavior:smooth}}
</style>
<main>
<div class="intro">
<h1>Thickness, before and after</h1>
<p>@INTRO@</p>
<div class="key">
<div><b>Columns</b> &mdash; the shape given, its removed faces magenta &middot; upstream OCCT 8.0.1 chain files at <code>91be8c4c71</code> &middot; the fork just before the fix &middot; the fork now</div>
<div><b>Rows</b> &mdash; the result from the removed face's side &middot; cut open, <span class="sw" style="background:var(--cut)"></span>cut faces orange</div>
<div><b>Labels</b> &mdash; <span style="color:var(--good)">green</span> valid at the expected volume, <span style="color:var(--bad)">red</span> what went wrong; <i>unhollowed</i> = the input's volume</div>
<div><b>Marks</b> &mdash; red edges not used once each way round, red faces fail <code>isValid()</code>, black = inside out, grey = it threw (input shown)</div>
</div>
</div>
<nav aria-label="Fixes">@NAV@</nav>
@SECTIONS@
</main>
"""

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.environ.get("THICK_WORK", "/tmp/thickness-pictures"), "page"))
