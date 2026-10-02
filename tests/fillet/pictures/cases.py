# The cases pictured in models/pictures: the suite's (run_tests.py) that a fix
# turned from failing to passing, each with the stage of the fix -- the fork
# commit just before it (STAGES) is the "before" column.
import os

try:
    import FreeCAD as App, Part
except ImportError:  # plain Python reads the tables only
    App = Part = None

MODELS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "models")

# stage -> (the fork commit just before the fix, the fix)
STAGES = {
    "s523": ("f723999a15", "16df68d224"),
}
# Upstream: the toolkit's files that differ from the fork, at the fork's base.
UPSTREAM = "91be8c4c71"


def _brep(name):
    def make():
        s = Part.Shape()
        s.read(os.path.join(MODELS, name))
        return s
    return make


SHAPES = {
    "boxcyl": _brep("issue523_box_cylinder.brep"),
}
# shape -> the face whose (u, v) outline the third row draws, found in a result
UVFACE = {
    # the corner's cylinder, not a fillet of the same radius
    "boxcyl": (lambda f: f.Surface.__class__.__name__ == "Cylinder"
               and abs(f.Surface.Radius - 2) < 1e-6 and abs(abs(f.Surface.Axis.z) - 1) < 1e-9),
}
NAMES = {
    "boxcyl": "10 box with a 3/4 cylinder r2 at a corner (#523, Part Connect)",
}
# name -> (stage, shape, the edge by its ends, radius, expected volume,
#          the point to zoom on, the direction the camera looks from)
CASES = {
    "seam_end_top_r1": ("s523", "boxcyl", ((2, 0, 10), (10, 0, 10)), 1.0, 1092.5263,
                        (2, 0, 10), (0.55, -1, 0.75)),
    "seam_end_bottom_r1": ("s523", "boxcyl", ((2, 0, 0), (10, 0, 0)), 1.0, 1092.5263,
                           (2, 0, 0), (0.55, -1, -0.75)),
    "seam_end_top_r0.5": ("s523", "boxcyl", ((2, 0, 10), (10, 0, 10)), 0.5, 1093.8183,
                          (2, 0, 10), (0.55, -1, 0.75)),
    "seam_end_top_r2": ("s523", "boxcyl", ((2, 0, 10), (10, 0, 10)), 2.0, 1087.3009,
                        (2, 0, 10), (0.55, -1, 0.75)),
}


def edge_between(shape, a, b):
    want = sorted([tuple(round(c, 6) for c in a), tuple(round(c, 6) for c in b)])
    for i, e in enumerate(shape.Edges):
        if sorted(tuple(round(c, 6) for c in v.Point) for v in e.Vertexes) == want:
            return e
    raise ValueError("no edge between %s and %s" % (a, b))
