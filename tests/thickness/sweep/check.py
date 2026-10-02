#!/usr/bin/env python3
# check.py <sweep output> [--write-ref <upstream's sweep output>]
#
# Judges a sweep by reference.txt: a run is right when it gives the
# reference's outcome -- the volume to 2e-3, the number of solids and of
# shells, or a refusal -- and leaves its input alone. Exit 1 when any run is
# not. (The sweep used to be judged by a ratio: valid, closed, and a volume
# plausible for a skin. That passed a lip of 382.017 where 422.796 was right.)
#
# Each reference line says where its value comes from:
#   U  upstream's chain gives the same result: two implementations agree
#   H  worked by hand, hand.py
#   S  the reference volume of a case in run_tests.py, worked by hand there
#   R  a ruling: a shape with every face removed is refused; the holed
#      cone's top inward is a skin and a sealed void
# --write-ref builds reference.txt from a sweep of the fork and one of
# upstream; it refuses a run that none of the four confirms.
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from hand import hand  # noqa: E402

TOL = 2e-3
REF = os.path.join(HERE, "reference.txt")


def load(path):
    runs = {}
    for line in open(path):
        p = line.split(None, 3)
        if len(p) == 4 and p[0] == "R":
            runs[p[2]] = p[3].strip()
    return runs


def outcome(text):
    """(status, volume, solids, shells) of a sweep line's result."""
    m = re.match(r"OK vol=(-?[\d.]+) solids=(\d+) shells=(\d+)$", text)
    if m:
        return ("OK", float(m.group(1)), int(m.group(2)), int(m.group(3)))
    if text.startswith("EXC") and "INPUT-CHANGED" not in text:
        return ("REFUSED", None, None, None)
    return ("BAD", text, None, None)


def same(a, b):
    if a[0] != b[0] or a[0] == "BAD":
        return False
    return a[0] == "REFUSED" or (abs(a[1] - b[1]) <= TOL and a[2:] == b[2:])


def suite_volumes():
    text = open(os.path.join(HERE, "..", "run_tests.py")).read()
    return [float(x) for x in re.findall(r"(?<![\w.])\d+\.\d+(?!\d)|(?<![\w.])\d{2,4}(?![\d.])", text)]


def write_ref(fork, upstream):
    suite = suite_volumes()
    lines, unconfirmed = [], []
    for tag, text in fork.items():
        o = outcome(text)
        shape, face, direction = tag.split("_")[:3]
        h = hand(tag)
        if o[0] == "REFUSED" and shape in ("sphere", "torus"):
            by = "R"
        elif o[0] == "OK" and o[3] == 2 and (shape, face, direction) == ("conehole", "f2", "-1"):
            by = "R"
        elif o[0] == "OK" and tag in upstream and same(o, outcome(upstream[tag])):
            by = "U"
        elif o[0] == "OK" and h is not None and abs(o[1] - h) <= TOL:
            by = "H"
        elif o[0] == "OK" and any(abs(o[1] - v) <= TOL for v in suite):
            by = "S"
        else:
            unconfirmed.append("%s %s" % (tag, text))
            continue
        lines.append("%-28s %s %s" % (tag, by, "REFUSED" if o[0] == "REFUSED" else
                                      "vol=%.3f solids=%d shells=%d" % o[1:]))
    if unconfirmed:
        print("not confirmed by upstream, hand.py, the suite or a ruling:")
        print("\n".join("  " + u for u in unconfirmed))
        return 1
    open(REF, "w").write("\n".join(lines) + "\n")
    print("reference.txt: %d runs" % len(lines))
    return 0


def check(runs):
    wrong, by, modes = [], {}, {}
    ref = {}
    for line in open(REF):
        tag, src, text = line.split(None, 2)
        ref[tag] = (src, text.strip())
    for tag, (src, text) in ref.items():
        want = ("REFUSED", None, None, None) if text == "REFUSED" else outcome("OK " + text)
        got = runs.get(tag)
        mode = tag.split("_", 3)[3]
        modes.setdefault(mode, [0, 0])[1] += 1
        if got is not None and same(outcome(got), want):
            by[src] = by.get(src, 0) + 1
            modes[mode][0] += 1
        else:
            wrong.append("%-28s got %s | want %s (%s)" % (tag, got or "nothing", text, src))
    for tag in runs:
        if tag not in ref:
            wrong.append("%-28s not in reference.txt" % tag)
    print("right: %d of %d  (%s)" % (len(ref) - sum(1 for w in wrong if "| want" in w), len(ref),
                                     "  ".join("%s=%d" % kv for kv in sorted(by.items()))))
    print("by mode: " + "  ".join("%s %d/%d" % (m, a, b) for m, (a, b) in sorted(modes.items())))
    for w in wrong:
        print("WRONG " + w)
    return 1 if wrong else 0


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[2] == "--write-ref":
        sys.exit(write_ref(load(sys.argv[1]), load(sys.argv[3])))
    sys.exit(check(load(sys.argv[1])))
