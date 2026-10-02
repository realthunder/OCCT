# FreeCAD (GUI) render.py -- one PNG per panel of $JOBS.
#
# Run as a macro (macOS: copied to render.FCMacro, since a .py argument opens
# in the editor there) or under xvfb-run on Linux. A panel is a result (or,
# when the fillet threw, the input, greyed) seen from the job's direction, on
# the whole shape or zoomed on the fillet's end. With "mark", edges not used
# once each way round are drawn red and faces failing isValid() are red.
# Faces of no area are left out. Next to each PNG a .json says what was done.
import json
import os

import FreeCAD as App
import FreeCADGui as Gui
import Part

V = App.Vector


def emit(s):
    os.write(1, (s + "\n").encode())


# The transaction log is off: the pictures need no history.
App.ParamGet("User parameter:BaseApp/Preferences/Document").SetInt("TransactionLog", 0)
# Coin draws, not the bgfx renderer: saveImage cannot capture a backend frame
# on every platform (macOS Metal falls back to a black GL capture), and the
# pictures are about the geometry. No navigation cube in them either.
App.ParamGet("User parameter:BaseApp/Preferences/View").SetInt("RenderCache", 0)
App.ParamGet("User parameter:BaseApp/Preferences/View").SetBool("ShowNaviCube", False)
jobs = json.load(open(os.environ["JOBS"]))
W, H = jobs["size"]


def defects(shape):
    # Each edge's uses, counted through the faces' wires: a seam occurs in
    # its face twice, once each way. (isSeam() is BRep_Tool::IsClosed, true
    # for an edge with two pcurves even where the face holds it once.)
    uses = {}
    for f in shape.Faces:
        for w in f.childShapes():
            for e in w.childShapes():
                if e.ShapeType != "Edge" or e.Degenerated:
                    continue
                uses.setdefault(e.hashCode(), (e, []))[1].append(e.Orientation)
    bad = [e for e, l in uses.values() if sorted(l) != ["Forward", "Reversed"]]
    badf = [i for i, f in enumerate(shape.Faces) if not f.isValid()]
    return bad, badf


def camera(eye, up, center, height):
    z = V(eye)
    z.normalize()
    y = up - z * up.dot(z)
    y.normalize()
    x = y.cross(z)
    m = App.Matrix(x.x, y.x, z.x, 0, x.y, y.y, z.y, 0, x.z, y.z, z.z, 0)
    r = App.Rotation(m)
    ax, ang = r.Axis, r.Angle
    pos = center + z * 1000
    return ("#Inventor V2.1 ascii\nOrthographicCamera {\n viewportMapping ADJUST_CAMERA\n"
            " position %g %g %g\n orientation %g %g %g %g\n nearDistance 1\n farDistance 3000\n"
            " aspectRatio 1\n focalDistance 1000\n height %g\n}\n"
            % (pos.x, pos.y, pos.z, ax.x, ax.y, ax.z, ang, height))


for job in jobs["panels"]:
    doc = App.newDocument("R")
    Gui.updateGui()
    shape = Part.Shape()
    shape.read(job["brep"])
    ghost = job.get("ghost", False)
    meta = {}
    if any(f.Area < 1e-9 for f in shape.Faces):
        keep = [f for f in shape.Faces if f.Area >= 1e-9]
        meta["dropped"] = len(shape.Faces) - len(keep)
        shape = Part.makeCompound(keep)
    bad, badf = defects(shape) if job.get("mark") else ([], [])
    obj = doc.addObject("Part::Feature", "S")
    obj.Shape = shape
    if bad:
        eo = doc.addObject("Part::Feature", "Bad")
        eo.Shape = Part.makeCompound(bad)
    doc.recompute()
    vo = obj.ViewObject
    base = (0.62, 0.70, 0.80) if not ghost else (0.75, 0.75, 0.75)
    vo.ShapeColor = base
    vo.LineColor = (0.15, 0.15, 0.2)
    vo.LineWidth = 1.5
    vo.Deviation = 0.02
    if ghost:
        vo.Transparency = 60
    if badf:
        cols = [base + (1.0,)] * len(shape.Faces)
        for i in badf:
            cols[i] = (0.95, 0.25, 0.25, 1.0)
        try:
            vo.DiffuseColor = cols
        except Exception:
            pass
    if bad:
        eo.ViewObject.LineColor = (0.9, 0.0, 0.0)
        eo.ViewObject.LineWidth = 5
        eo.ViewObject.PointColor = (0.9, 0.0, 0.0)
    view = Gui.ActiveDocument.ActiveView
    eye = V(*job["eye"])
    eye.normalize()
    up = V(0, 0, 1) if abs(eye.z) < 0.95 else V(0, 1, 0)
    view.setCamera(camera(eye, up, V(*job["center"]), job["height"]))
    for _ in range(3):
        Gui.updateGui()
    view.redraw()
    view.saveImage(job["png"], W, H, "#FFFFFF")
    App.closeDocument(doc.Name)
    meta.update(badedges=len(bad), badfaces=len(badf))
    json.dump(meta, open(job["png"] + ".json", "w"))
    emit("PNG %s badedges=%d badfaces=%d" % (job["png"], len(bad), len(badf)))
emit("DONE")
os._exit(0)
