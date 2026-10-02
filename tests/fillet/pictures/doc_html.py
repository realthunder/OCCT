# python doc_html.py <out.html> -- models/Fillet.md as one HTML page.
#
# The page references its pictures as pictures/<case>.png, beside it; publish
# models/pictures/*.png at those paths with it. Needs python-markdown (the
# FreeCAD conda env has it).
import html
import os
import re
import sys

import markdown

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "models", "Fillet.md")

text = open(SRC, encoding="utf-8").read()
# the first line is the page's heading, set apart from the body
first, body = text.split("\n", 1)
title = first.lstrip("# ").strip()

md = markdown.Markdown(extensions=["toc", "tables", "fenced_code"],
                       extension_configs={"toc": {"toc_depth": "2-3"}})
out = md.convert(body)


def figures(m):
    figs = []
    for alt, src in re.findall(r'<img alt="([^"]*)" src="([^"]*)"\s*/?>', m.group(0)):
        figs.append('<figure><a href="%s"><img alt="%s" src="%s" loading="lazy"></a>'
                    '<figcaption><code>%s</code></figcaption></figure>' % (src, alt, src, alt))
    return "\n".join(figs)


# a paragraph of pictures (one per line in the source) becomes figures
out = re.sub(r"<p>(?:\s*<img [^>]*>\s*)+</p>", figures, out)

# the fixes, for the contents list: h3 under "The fixes", with the commit
fixes = []
for hid, inner in re.findall(r'<h3 id="([^"]+)">(.*?)</h3>', out):
    commit = re.search(r"<code>([0-9a-f]{10})</code>", inner)
    issue = re.search(r"#(\d+)\)", inner)
    name = re.sub(r"\s*\(.*\)\s*$", "", re.sub(r"<[^>]+>", "", inner)).strip()
    fixes.append((hid, name, commit.group(1) if commit else "", issue.group(1) if issue else ""))
nfig = len(re.findall(r"<figure>", out))

toc = "\n".join(
    '<li><a href="#%s"><span class="fix-name">%s</span>'
    '<span class="fix-meta">%s%s</span></a></li>' % (
        hid, html.escape(name), ("<code>%s</code>" % c) if c else "",
        (' <span class="issue">#%s</span>' % i) if i else "")
    for hid, name, c, i in fixes)

page = """<title>TKFillet Fixes</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:ital,wght@0,400;0,500;0,600;1,400&family=IBM+Plex+Sans+Condensed:wght@500;600&display=swap">
<style>
:root {
  --bg: #f3f5f8;
  --paper: #ffffff;
  --ink: #18202b;
  --ink-2: #4a5566;
  --rule: #d8dee7;
  --accent: #33648f;
  --accent-soft: #e3ecf5;
  --code-bg: #e9edf2;
  --fig-bg: #ffffff;
  --sans: "IBM Plex Sans", "Helvetica Neue", Arial, sans-serif;
  --cond: "IBM Plex Sans Condensed", "IBM Plex Sans", "Arial Narrow", sans-serif;
  --mono: "IBM Plex Mono", ui-monospace, Menlo, Consolas, monospace;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --bg: #12161c;
    --paper: #1a2029;
    --ink: #e2e7ee;
    --ink-2: #a3aebd;
    --rule: #2e3744;
    --accent: #8db7de;
    --accent-soft: #223244;
    --code-bg: #252d39;
    --fig-bg: #f4f6f9;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --bg: #12161c;
  --paper: #1a2029;
  --ink: #e2e7ee;
  --ink-2: #a3aebd;
  --rule: #2e3744;
  --accent: #8db7de;
  --accent-soft: #223244;
  --code-bg: #252d39;
  --fig-bg: #f4f6f9;
}
body {
  background: var(--bg);
  color: var(--ink);
  font: 16px/1.62 var(--sans);
  padding-inline: 16px;
  padding-block: 0 64px;
}
.wrap { max-width: 1112px; margin: 0 auto; }
header.top {
  padding-block: 40px 28px;
  border-bottom: 1px solid var(--rule);
  display: grid;
  gap: 10px;
}
.eyebrow {
  font: 500 12px/1.4 var(--mono);
  letter-spacing: .06em;
  text-transform: uppercase;
  color: var(--ink-2);
}
h1 {
  font: 600 clamp(28px, 4.2vw, 42px)/1.12 var(--cond);
  margin: 0;
  text-wrap: balance;
}
.lede { max-width: 66ch; margin: 0; color: var(--ink-2); }
.layout {
  display: grid;
  grid-template-columns: minmax(0, 1fr);
  gap: 32px;
  margin-top: 28px;
}
@media (min-width: 1000px) {
  .layout { grid-template-columns: 250px minmax(0, 1fr); gap: 44px; }
  nav.fixes { position: sticky; top: calc(env(safe-area-inset-top, 0px) + 20px); align-self: start; }
}
nav.fixes h2 {
  font: 600 13px/1.3 var(--cond);
  letter-spacing: .05em;
  text-transform: uppercase;
  color: var(--ink-2);
  margin: 0 0 10px;
}
nav.fixes ol { list-style: none; margin: 0; padding: 0; display: grid; gap: 4px; }
nav.fixes a {
  display: grid;
  gap: 2px;
  padding: 8px 10px;
  border-radius: 6px;
  color: var(--ink);
  text-decoration: none;
  border-left: 3px solid var(--rule);
}
nav.fixes a:hover, nav.fixes a:focus-visible { background: var(--accent-soft); border-left-color: var(--accent); }
.fix-name { font-size: 14px; line-height: 1.35; }
.fix-meta { font-size: 12px; color: var(--ink-2); }
.fix-meta code { background: none; padding: 0; }
.issue { font-family: var(--mono); }
nav.fixes .count { margin-top: 12px; font-size: 13px; color: var(--ink-2); }
article { min-width: 0; }
article > * { max-width: 70ch; }
article > figure, article > table { max-width: none; }
article h2 {
  font: 600 25px/1.2 var(--cond);
  margin: 44px 0 10px;
  padding-top: 18px;
  border-top: 1px solid var(--rule);
  text-wrap: balance;
}
article h2:first-child { margin-top: 0; padding-top: 0; border-top: 0; }
article h3 {
  font: 600 20px/1.3 var(--cond);
  margin: 40px 0 8px;
  text-wrap: balance;
  scroll-margin-top: 20px;
}
article h3 code { font-size: .78em; font-weight: 500; }
article p, article ul, article ol { margin: 0 0 14px; }
article li + li { margin-top: 6px; }
a { color: var(--accent); }
a:focus-visible { outline: 2px solid var(--accent); outline-offset: 2px; }
code {
  font: .86em/1 var(--mono);
  background: var(--code-bg);
  padding: .16em .36em;
  border-radius: 4px;
}
pre { overflow-x: auto; background: var(--code-bg); padding: 12px 14px; border-radius: 6px; }
pre code { background: none; padding: 0; }
figure {
  margin: 18px 0 26px;
  background: var(--fig-bg);
  border: 1px solid var(--rule);
  border-radius: 6px;
  overflow: hidden;
}
figure a { display: block; }
figure img { display: block; width: 100%; height: auto; max-width: 1080px; margin: 0 auto; }
figcaption {
  padding: 7px 12px;
  border-top: 1px solid var(--rule);
  background: var(--paper);
  font-size: 13px;
  color: var(--ink-2);
}
figcaption code { background: none; padding: 0; }
.source { margin-top: 40px; font-size: 13px; color: var(--ink-2); }
@media (prefers-reduced-motion: no-preference) { html { scroll-behavior: smooth; } }
</style>
<div class="wrap">
<header class="top">
  <div class="eyebrow">realthunder OCCT fork &middot; LinkVibe-801 &middot; TKFillet / ChFi3d</div>
  <h1>%(title)s</h1>
  <p class="lede">Each fix the fork makes to OCCT's fillet builder, with the cases that show it:
  upstream OCCT, the fork before the fix and the fork after, side by side.</p>
</header>
<div class="layout">
<nav class="fixes" aria-label="The fixes">
  <h2>The fixes</h2>
  <ol>
%(toc)s
  </ol>
  <p class="count">%(nfig)d pictures &middot; source <code>tests/fillet/models/Fillet.md</code></p>
</nav>
<article>
%(body)s
<p class="source">Rendered from <code>tests/fillet/models/Fillet.md</code> in the OCCT fork; the suite is
<code>tests/fillet/run_tests.py</code>, the pictures come from <code>tests/fillet/pictures/make_pictures.sh</code>.</p>
</article>
</div>
</div>
"""
for key, value in (("%(title)s", html.escape(title)), ("%(toc)s", toc), ("%(nfig)d", str(nfig)),
                   ("%(body)s", out)):
    page = page.replace(key, value)

open(sys.argv[1], "w", encoding="utf-8").write(page)
print(sys.argv[1], len(page), "bytes,", len(fixes), "fixes,", nfig, "pictures")
