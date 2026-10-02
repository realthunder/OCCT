# FreeCAD (GUI, under xvfb-run) render.py -- one PNG per panel of $JOBS.
#
# A panel is a result (or, when the thickness threw, the input, greyed) seen
# from the removed face's side, z up. An "input" panel is the shape given,
# its removed faces ("removed", indices into its faces) in magenta. With "cut", the half towards the camera
# is cut away by a plane through the removed face's centre (axis-aligned where
# one splits the input evenly, else through a curved face's axis) and the cut
# faces are orange; when the boolean fails or the solid is inside out a clip
# plane does it instead, with no caps. Without "cut", with "mark", the edges
# not used once each way round are drawn red, and faces failing isValid() are
# red. Faces of no area are left out (the mesher draws nothing otherwise).
# Next to each PNG a .json says what was done: clipped, huge, dropped.
import os, json
import FreeCAD as App, FreeCADGui as Gui, Part
V = App.Vector
def emit(s): os.write(1, (s + "\n").encode())
# The transaction log is off: the pictures need no history. (Its worker once
# crashed here writing a shape the viewer was meshing -- FreeCAD
# docs/TransactionLog.md sec 27.97, fixed in sec 27.98.)
App.ParamGet("User parameter:BaseApp/Preferences/Document").SetInt("TransactionLog", int(os.environ.get("TLOG", "0")))
jobs = json.load(open(os.environ["JOBS"]))
W, H = jobs["size"]
Gui.showMainWindow() if not Gui.getMainWindow() else None
def defects(shape):
    uses = {}
    for f in shape.Faces:
        for e in f.Edges:
            if e.Degenerated: continue
            h = e.hashCode()
            l = uses.setdefault(h, (e, []))[1]
            try:
                seam = e.isSeam(f)
            except Exception:
                seam = False
            l.extend(["Forward", "Reversed"] if seam else [e.Orientation])
    bad = []
    for h, (e, l) in uses.items():
        if sorted(l) != ["Forward", "Reversed"]:
            bad.append(e)
    badf = [i for i, f in enumerate(shape.Faces) if not f.isValid()]
    return bad, badf
def camera(eye, up, center, height):
    z = V(eye); z.normalize()
    y = up - z * up.dot(z); y.normalize()
    x = y.cross(z)
    m = App.Matrix(x.x, y.x, z.x, 0, x.y, y.y, z.y, 0, x.z, y.z, z.z, 0)
    r = App.Rotation(m)
    ax = r.Axis; ang = r.Angle
    pos = center + z * 1000
    return ("#Inventor V2.1 ascii\nOrthographicCamera {\n viewportMapping ADJUST_CAMERA\n"
            " position %g %g %g\n orientation %g %g %g %g\n nearDistance 1\n farDistance 3000\n"
            " aspectRatio 1\n focalDistance 1000\n height %g\n}\n"
            % (pos.x, pos.y, pos.z, ax.x, ax.y, ax.z, ang, height))
bbcache = {}
def bbcache_get(f):
    if f not in bbcache:
        t = Part.Shape(); t.read(f); bbcache[f] = t.BoundBox
    return bbcache[f]
for job in jobs["panels"]:
    doc = App.newDocument("R")
    Gui.updateGui()
    n = V(*job["normal"]); n.normalize()
    c = V(*job["center"])
    isos = [V(1, -1, .7), V(-1, -1, .7), V(1, 1, .7), V(-1, 1, .7)]
    for i in isos: i.normalize()
    iso = max(isos, key=lambda i: i.dot(n))
    eye = n + iso * 0.9; eye.normalize()
    ib = bbcache_get(job["bbox_breps"][0])
    axes = [V(1, 0, 0), V(0, 1, 0), V(0, 0, 1)]
    axes = [a for a in axes if abs(a.dot(n)) < 1e-6]
    if abs(n.z) < 0.99:
        nz = n.cross(V(0, 0, 1)); nz.normalize()
        if all(abs(abs(nz.dot(a)) - 1) > 1e-6 for a in axes): axes.append(nz)
    if job.get("cutaxis"):
        u = V(*job["cutaxis"])
        if u.dot(eye) < 0: u = u * -1
    elif axes:
        def split(a):
            ext = max(ib.XLength * abs(a.x) + ib.YLength * abs(a.y) + ib.ZLength * abs(a.z), 1e-9)
            if abs(a.x) > 1e-6 and abs(a.y) > 1e-6:
                return 0.005 - 1e-3 * abs(a.dot(eye))  # through the face's own axis
            return abs((c - ib.Center).dot(a)) / ext + (0.01 if abs(a.z) > 0.9 else 0) - 1e-3 * abs(a.dot(eye))
        u = V(min(axes, key=split))
        if u.dot(eye) < 0: u = u * -1
    else:
        u = eye - n * eye.dot(n)
        if u.Length < 1e-6: u = iso - n * iso.dot(n)
        u.normalize()
    w = n.cross(u)
    bb = App.BoundBox()
    for f in job["bbox_breps"]:
        # A result without bounds does not set the row's scale.
        if bbcache_get(f).DiagonalLength < 1e6:
            bb.add(bbcache_get(f))
    height = bb.DiagonalLength * 1.02
    shape = Part.Shape(); shape.read(job["brep"])
    ghost = job.get("ghost", False)
    meta = {}
    # The removed faces, as they are in the shape read: after a cut, the
    # faces lying on one of them are found again by their edges.
    removed = [shape.Faces[i] for i in job.get("removed", [])] if job.get("input") else []
    def on_removed(g):
        if not removed: return False
        for f in removed:
            if f.Surface.__class__ is not g.Surface.__class__: continue
            pts = [e.valueAt((e.FirstParameter + e.LastParameter) / 2) for e in g.Edges if not e.Degenerated]
            if pts and all(f.distToShape(Part.Vertex(p))[0] < 1e-5 for p in pts): return True
        return False
    # A face without bounds (upstream's inward cone, a volume of 2e100)
    # cannot be meshed -- the mesher crashes on it; it is left out.
    def bounded(f):
        try:
            return f.BoundBox.DiagonalLength < 1e6 and abs(f.Area) < 1e12
        except Exception:
            return False
    if not ghost and not all(bounded(f) for f in shape.Faces):
        keep = [f for f in shape.Faces if bounded(f)]
        meta["unbounded"] = len(shape.Faces) - len(keep)
        shape = Part.makeCompound(keep)
    sb = shape.BoundBox
    if not ghost and sb.DiagonalLength > 20 * bb.DiagonalLength:
        meta["huge"] = sb.DiagonalLength
        bb = sb
        height = sb.DiagonalLength * 1.02
    try:
        inside_out = not ghost and shape.Volume < 0
    except Exception:
        inside_out = False
    cap = None
    clipped = False
    if job.get("cut") and (inside_out or meta.get("huge")):
        clipped = True
        meta["clipped"] = 1
    elif job.get("cut"):
        # keep the half behind the plane through the removed face's centre, normal u
        big = bb.DiagonalLength * 4
        p0 = c - n * (big / 2) + w * (big / 2)
        half = Part.makePlane(big, big, p0, u, n).extrude(u * -big)
        try:
            cutres = shape.common(half)
            if cutres.isNull() or cutres.Volume < 1e-6 or not cutres.Faces:
                raise RuntimeError("empty")
            shape2 = cutres
            capf = [i for i, f in enumerate(shape2.Faces) if f.Surface.__class__.__name__ == "Plane"
                    and abs(abs(f.normalAt(0, 0).dot(u)) - 1) < 1e-6 and abs((f.CenterOfMass - c).dot(u)) < 1e-6]
            shape = shape2; cap = capf
        except Exception:
            clipped = True
            meta["clipped"] = 1
    # A cut that cannot be measured is no cut: the whole shape, clipped.
    try:
        [f.Area for f in shape.Faces]
    except Exception:
        shape = Part.Shape(); shape.read(job["brep"])
        cap = None
        clipped = True
        meta["clipped"] = 1
    if any(f.Area < 1e-9 for f in shape.Faces):
        keep = [f for f in shape.Faces if f.Area >= 1e-9]
        meta["dropped"] = len(shape.Faces) - len(keep)
        shape = Part.makeCompound(keep)
        cap = None
    bad, badf = ([], [])
    if job.get("mark") and not job.get("cut"):
        bad, badf = defects(shape)
    remf = [i for i, f in enumerate(shape.Faces) if on_removed(f)]
    obj = doc.addObject("Part::Feature", "S"); obj.Shape = shape
    if bad:
        eo = doc.addObject("Part::Feature", "Bad"); eo.Shape = Part.makeCompound(bad)
    doc.recompute()
    vo = obj.ViewObject
    base = (0.62, 0.70, 0.80) if not ghost else (0.75, 0.75, 0.75)
    vo.ShapeColor = base
    vo.LineColor = (0.15, 0.15, 0.2)
    vo.LineWidth = 1.5
    vo.Deviation = 0.05
    if ghost:
        vo.Transparency = 60
    if cap or badf or remf:
        cols = [base + (1.0,)] * len(shape.Faces)
        for i in cap or []: cols[i] = (0.93, 0.55, 0.25, 1.0)
        for i in badf: cols[i] = (0.95, 0.25, 0.25, 1.0)
        for i in remf: cols[i] = (0.85, 0.25, 0.65, 1.0)
        try:
            vo.DiffuseColor = cols
        except Exception:
            pass
    if bad:
        eo.ViewObject.LineColor = (0.9, 0.0, 0.0)
        eo.ViewObject.LineWidth = 5
        eo.ViewObject.PointSize = 1
        eo.ViewObject.PointColor = (0.9, 0.0, 0.0)
    view = Gui.ActiveDocument.ActiveView
    if clipped:
        from pivy import coin
        sg = view.getSceneGraph()
        cp = coin.SoClipPlane()
        cp.plane.setValue(coin.SbPlane(coin.SbVec3f(-u.x, -u.y, -u.z), coin.SbVec3f(c.x, c.y, c.z)))
        sg.insertChild(cp, 0)
    up = V(0, 0, 1) if abs(eye.z) < 0.95 else V(0, 1, 0)
    view.setCamera(camera(eye, up, bb.Center, height))
    for _ in range(3): Gui.updateGui()
    view.redraw()
    view.saveImage(job["png"], W, H, "#FFFFFF")
    if clipped:
        sg.removeChild(cp)
    App.closeDocument(doc.Name)
    json.dump(meta, open(job["png"] + ".json", "w"))
    emit("PNG %s %s badedges=%d badfaces=%d" % (job["png"], "clipped" if clipped else "", len(bad), len(badf)))
emit("DONE")
os._exit(0)
