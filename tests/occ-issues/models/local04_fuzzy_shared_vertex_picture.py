# local04: makes local04_fuzzy_shared_vertex.png, see ../README.md.
# Run inside the FreeCAD GUI (it renders a 3D view), e.g. through the MCP
# console: python scripts/mcp_run.py local04_fuzzy_shared_vertex_picture.py
#
# Left: the model -- A and B, and the vertex V they share. Right: the plane
# x = 0 through V at micrometre scale: FA on z = 0, FB tilted 10 deg and
# passing 1.74 um above V, their section 10 um from V, and the two radii the
# pave filler measures that distance against. The result shapes are the same
# before and after the fix; what changed is the input vertex's tolerance.
import math
import os
import tempfile

import FreeCAD as App
import FreeCADGui as Gui
import Part
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.image import imread
from matplotlib.patches import Circle, Polygon

HERE = os.path.dirname(os.path.abspath(__file__))
BREP = os.path.join(HERE, "local04_fuzzy_shared_vertex.brep")
OUT = os.path.join(HERE, "local04_fuzzy_shared_vertex.png")

THETA = 10.0      # deg, the tilt of FB
D = 10.0          # um, V to the section line
TOLV = 1.94       # um, V's tolerance in the model
TOLC = 2.96       # um, the section's (tangential) tolerance at fuzzy 3.162e-7
FUZZY = 0.3162    # um
EST = 2 * (TOLV + TOLC)            # EstimatePaveOnCurve before the fix
PUT = 2 * (TOLV + TOLC + FUZZY)    # PutPaveOnCurve, and the estimate after it


def render(path):
    doc = App.newDocument("local04")
    a, b = Part.read(BREP).childShapes()
    for name, shape, color in (("A", a, (0.55, 0.65, 0.8)), ("B", b, (0.85, 0.75, 0.5))):
        obj = doc.addObject("Part::Feature", name)
        obj.Shape = shape
        obj.ViewObject.ShapeColor = color
        obj.ViewObject.Transparency = 30 if name == "B" else 0
    mark = doc.addObject("Part::Feature", "V")
    mark.Shape = Part.makeSphere(0.35)
    mark.ViewObject.ShapeColor = (0.85, 0.1, 0.1)
    doc.recompute()
    view = Gui.ActiveDocument.ActiveView
    view.viewIsometric()
    view.setCameraOrientation(App.Rotation(App.Vector(0, 0, 1), -150).multiply(
        App.Rotation(App.Vector(1, 0, 0), 55)).Q)
    # Frame B and V: centre the (orthographic) camera on V, 16 mm across
    view.fitAll()
    from pivy import coin
    cam = view.getCameraNode()
    look = cam.orientation.getValue().multVec(coin.SbVec3f(0, 0, -1))
    focal = cam.position.getValue() + look * cam.focalDistance.getValue()
    cam.position.setValue(cam.position.getValue() - focal)
    if hasattr(cam, "height"):
        cam.height.setValue(16)
    Gui.updateGui()
    view.saveImage(path, 1200, 900, "White")
    App.closeDocument(doc.Name)


def closeup(ax):
    t = math.radians(THETA)
    lo, hi = -16.0, 8.0
    # A below FA (z = 0), B above FB (z = tan(t) (y + D))
    ax.add_patch(Polygon([(lo, 0), (hi, 0), (hi, -8), (lo, -8)], color=(0.55, 0.65, 0.8), alpha=0.5))
    zb = lambda y: math.tan(t) * (y + D)
    ax.add_patch(Polygon([(lo, zb(lo)), (hi, zb(hi)), (hi, 8), (lo, 8)], color=(0.85, 0.75, 0.5), alpha=0.5))
    ax.plot([lo, hi], [0, 0], color=(0.2, 0.3, 0.5), lw=1.5)
    ax.plot([lo, hi], [zb(lo), zb(hi)], color=(0.5, 0.4, 0.1), lw=1.5)
    ax.text(hi - 0.3, -0.9, "FA (A's face)", ha="right", va="top", fontsize=9)
    ax.text(hi - 0.3, zb(hi) + 0.6, "FB (B's face, %g deg)" % THETA, ha="right", va="bottom", fontsize=9)
    # the radii, centred on V
    est = Circle((0, 0), EST, fill=False, ls="--", color="tab:blue", lw=1.3,
                 label="EstimatePaveOnCurve (before the fix): replace V if the section is\n"
                       "within 2(tolV+tolC) = %.1f um -- it is not" % EST)
    put = Circle((0, 0), PUT, fill=False, color="tab:red", lw=1.3,
                 label="PutPaveOnCurve: put V on the section if within\n"
                       "2(tolV+tolC+fuzzy) = %.1f um -- it is" % PUT)
    ax.add_patch(est)
    ax.add_patch(put)
    ax.legend(handles=[est, put], loc="lower right", fontsize=8, framealpha=0.95)
    # V and the section
    ax.plot([0], [0], "o", color=(0.85, 0.1, 0.1), ms=7)
    ax.text(0.5, 0.6, "V (shared)", fontsize=9)
    ax.plot([-D], [0], "kx", ms=9, mew=2)
    ax.annotate("section of FA and FB\n%g um from V" % D, (-D, 0), (-15.5, 5.0), fontsize=9,
                arrowprops=dict(arrowstyle="->", lw=0.8))
    ax.set_xlim(lo, hi)
    ax.set_ylim(-8, 8)
    ax.set_aspect("equal")
    ax.set_xlabel("y (um)")
    ax.set_ylabel("z (um)")
    ax.set_title("At V, plane x = 0 (to scale)", fontsize=11)


tmp = os.path.join(tempfile.gettempdir(), "local04_model.png")
render(tmp)
fig = Figure(figsize=(15, 6.6), dpi=100)
FigureCanvasAgg(fig)
left = fig.add_axes([0.0, 0.14, 0.45, 0.78])
left.imshow(imread(tmp))
left.axis("off")
left.set_title("local04: A (blue) and B (tan) share vertex V (red)", fontsize=11)
right = fig.add_axes([0.52, 0.18, 0.45, 0.74])
closeup(right)
fig.text(0.5, 0.025,
         "Fuzzy 3.162e-7: 9.8 um < 10 um <= 10.4 um, so V is never replaced, goes on the curve, and the "
         "INPUT vertex's tolerance is raised 1.94 um -> 10 um.\nFixed: the estimate uses the put's radius "
         "(red), V is replaced by a new vertex first, the input keeps 1.94 um; results unchanged.",
         ha="center", va="bottom", fontsize=9.5)
fig.savefig(OUT)
os.remove(tmp)
print("wrote", OUT)
