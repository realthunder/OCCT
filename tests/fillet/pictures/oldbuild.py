# Compile the toolkits in $TOOLKITS (default TKFillet) again from the build
# tree's own compile commands, with the sources under $U/src and the headers
# under $U/inc taking the place of the tree's, and link them into $U/lib.
# Called by mkold.sh. Linux (.so) and macOS (.dylib) alike.
import os
import re
import shlex
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

U = os.environ["U"]
SRC = os.environ["OCCT"].rstrip("/") + "/"
B = os.environ["OCCT_BUILD"]
RUN = os.environ.get("RUN", os.path.expanduser("~/works/sw/fcad/.conda/run.sh"))
TKS = os.environ.get("TOOLKITS", "TKFillet").split()
JOBS = int(os.environ.get("JOBS", "4"))
out = []
for tk in TKS:
    out += subprocess.run([RUN, "ninja", "-C", B, "-t", "commands", tk],
                          capture_output=True, text=True).stdout.splitlines()
dirs = tuple("%s.dir" % tk for tk in TKS)
libre = re.compile(r"-o \S*lib(%s)\.(so|[0-9.]*dylib)(\.[0-9.]+)?( |$)" % "|".join(TKS))
comp = sorted({l for l in out if any(d in l for d in dirs) and " -c " in l})
links = sorted({l for l in out if libre.search(l)})
print(len(comp), "objects,", len(links), "libraries")


def fix(l):
    a = shlex.split(l)
    r = []
    i = 0
    while i < len(a):
        x = a[i]
        if x == "-MD":
            i += 1
            continue
        if x in ("-MT", "-MF"):
            i += 2
            continue
        if x == "-o":
            o = U + "/obj/" + a[i + 1]
            os.makedirs(os.path.dirname(o), exist_ok=True)
            r += ["-o", o]
            i += 2
            continue
        if x.startswith(SRC + "src/") and os.path.exists(U + "/" + x[len(SRC):]):
            x = U + "/" + x[len(SRC):]
        r.append(x)
        i += 1
    r.insert(1, "-I" + U + "/inc")
    return r


def run(l):
    p = subprocess.run(fix(l), cwd=B, capture_output=True, text=True)
    if p.returncode:
        print(p.stderr[-3000:])
        return 1
    return 0


with ThreadPoolExecutor(JOBS) as ex:
    bad = sum(ex.map(run, comp))
print("compile failures", bad)
if bad:
    sys.exit(1)
for l in links:
    seg = [s for s in l.split("&&") if " -o " in s][0]
    a = shlex.split(seg)
    r = [U + "/obj/" + x if x.endswith(".o") and any(d in x for d in dirs) else x for x in a]
    j = r.index("-o")
    name = os.path.basename(r[j + 1])
    r[j + 1] = U + "/lib/" + name
    p = subprocess.run(r, cwd=B, capture_output=True, text=True)
    print("link", name, p.returncode, p.stderr[-2000:])
    if p.returncode:
        sys.exit(1)
    # the names the loader asks for: libTKFillet.8.0.dylib, libTKFillet.so.8.0
    m = re.match(r"(lib\w+)\.(\d+)\.(\d+)\.(\d+)\.dylib$", name)
    if m:
        os.symlink(name, U + "/lib/%s.%s.%s.dylib" % m.group(1, 2, 3))
    m = re.match(r"(lib\w+\.so)\.(\d+)\.(\d+)\.(\d+)$", name)
    if m:
        os.symlink(name, U + "/lib/%s.%s.%s" % m.group(1, 2, 3))
