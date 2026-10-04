# local05: makes local05_step_hidden_component.png, see ../README.md.
# Plain matplotlib (Agg canvas, no pyplot, no FreeCAD):
#   python local05_step_hidden_component_picture.py
#
# The three cases local05_step_hidden_component.cpp writes: an assembly of
# two boxes, the second component invisible, the colour on the component
# itself ("own"), on the referred shapes only ("ref") or nowhere ("none").
# Nothing differs on screen between before and after -- before, the writer
# faults; after, it writes -- so the picture shows where each colour lives
# and what the writer made of the hidden instance.
import os

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "local05_step_hidden_component.png")

RED = "#d62728"
GREY = "#f2f2f2"
CASES = [
    ("own", "colour on the hidden component", False, True,
     "before: written, instance red", "after: the same"),
    ("ref", "colour on the shapes only\n(every FreeCAD GUI export)", True, False,
     "before: ACCESS VIOLATION, read at 0x18", "after: written, instance default white"),
    ("none", "no colour anywhere", False, False,
     "before: ACCESS VIOLATION, read at 0x18", "after: written, instance default white"),
]


def node(ax, x, y, text, fill, dashed=False, w=0.38, h=0.12):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle="round,pad=0.01", facecolor=fill,
                                edgecolor="black", linewidth=1.4,
                                linestyle="--" if dashed else "-"))
    ax.text(x, y, text, ha="center", va="center", fontsize=9,
            color="white" if fill == RED else "black")


def arrow(ax, x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                                 mutation_scale=12, linewidth=1.1, color="#444444"))


fig = Figure(figsize=(12, 5.2), dpi=100)
FigureCanvasAgg(fig)
fig.suptitle("local05 -- STEPCAFControl_Writer and an invisible assembly component\n"
             "MakeSTEPStyles: invisible + no own colour -> \"default white\" branch -> "
             "setDefaultInstanceColor(theOverride = null)", fontsize=11)
for i, (name, title, shapeRed, ownRed, before, after) in enumerate(CASES):
    ax = fig.add_axes([0.02 + i * 0.33, 0.04, 0.30, 0.74])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.05)
    ax.axis("off")
    ax.text(0.5, 1.0, '"%s": %s' % (name, title), ha="center", va="top", fontsize=10,
            fontweight="bold")
    node(ax, 0.5, 0.80, "assembly", GREY)
    node(ax, 0.25, 0.56, "component 1", GREY)
    node(ax, 0.75, 0.56, "component 2\nINVISIBLE", RED if ownRed else "white", dashed=True)
    node(ax, 0.25, 0.32, "Box1 shape", RED if shapeRed else GREY)
    node(ax, 0.75, 0.32, "Box2 shape", RED if shapeRed else GREY)
    arrow(ax, 0.45, 0.74, 0.30, 0.62)
    arrow(ax, 0.55, 0.74, 0.70, 0.62)
    arrow(ax, 0.25, 0.50, 0.25, 0.38)
    arrow(ax, 0.75, 0.50, 0.75, 0.38)
    bad = "VIOLATION" in before
    ax.text(0.5, 0.17, before, ha="center", fontsize=9,
            color=RED if bad else "black", fontweight="bold" if bad else "normal")
    ax.text(0.5, 0.08, after, ha="center", fontsize=9)
fig.text(0.5, 0.01, "The writer styles a component from its own label only; FreeCAD's "
         "workaround tested GetInstanceColor(), which falls back to the shape's colour, "
         "so \"ref\" went through unprotected.", ha="center", fontsize=8.5, color="#444444")
fig.savefig(OUT)
print(OUT)
