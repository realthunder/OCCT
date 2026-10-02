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
    "s962": ("d8ba4ad7eb", "cf96757c36"),
    "s962b": ("1482f7d4e2", "b38f0d910c"),
}
# Upstream: the toolkit's files that differ from the fork, at the fork's base.
UPSTREAM = "91be8c4c71"


def _brep(name):
    def make():
        s = Part.Shape()
        s.read(os.path.join(MODELS, name))
        return s
    return make


def _slot_wall():
    # run_tests.py's slotted block, its slot's wall split 0.7 from the outer face
    V = App.Vector
    s = Part.makeBox(10, 10, 14).cut(Part.makeBox(12, 3, 10, V(-1, 3, 5)))
    return s.generalFuse([Part.LineSegment(V(9.3, 3, 5), V(9.3, 3, 20)).toShape()])[0].Solids[0]


def _fin():
    # run_tests.py's fin padded flush with the block's side
    V = App.Vector
    return Part.makeBox(10, 10, 5).fuse(Part.makeBox(10, 3, 9, V(0, 0, 5)))


def _arm_on_block():
    # run_tests.py's arm fused to a block, the seam of their bottoms from y=27.2 to 10.2
    V = App.Vector
    pts = [V(17.356, -9.425, -9.75), V(38.5, 10.2, -9.75), V(38.5, 27.2, -9.75),
           V(10.898, 16.471, -9.75)]
    arm = Part.Face(Part.makePolygon(list(reversed(pts)) + [pts[-1]])).extrude(V(0, 0, 6))
    bp = [V(38.5, 27.2, -9.75), V(38.5, 10.2, -9.75), V(50, 10.2, -9.75), V(50, 27.2, -9.75)]
    return Part.Face(Part.makePolygon(bp + [bp[0]])).extrude(V(0, 0, 6)).fuse(arm)


SHAPES = {
    "boxcyl": _brep("issue523_box_cylinder.brep"),
    "slotwall": _slot_wall,
    "p962": _brep("issue962_pocket002.brep"),
    "fin": _fin,
    "armfoot": _arm_on_block,
}


def _plane_at(axis, value, lo=None, hi=None):
    """A planar face lying in <axis> = value, its bounding box's centre
    between lo and hi on the other axes (dicts axis -> bound)."""
    def pick(f):
        if f.Surface.__class__.__name__ != "Plane":
            return False
        b = f.BoundBox
        mn, mx = getattr(b, axis.upper() + "Min"), getattr(b, axis.upper() + "Max")
        if abs(mn - value) > 1e-6 or abs(mx - value) > 1e-6:
            return False
        c = b.Center
        return all(lo[a] <= getattr(c, a) for a in lo or {}) and all(
            getattr(c, a) <= hi[a] for a in hi or {})
    return pick


# shape -> the face whose (u, v) outline the third row draws, found in a result
UVFACE = {
    # the corner's cylinder, not a fillet of the same radius
    "boxcyl": (lambda f: f.Surface.__class__.__name__ == "Cylinder"
               and abs(f.Surface.Radius - 2) < 1e-6 and abs(abs(f.Surface.Axis.z) - 1) < 1e-9),
    # the slot's floor, where the fillet ends
    "slotwall": _plane_at("z", 5),
    # the slot's floor at e101's foot (y 13.7..17.1); for e50 see UVFACE_CASE
    "p962": _plane_at("z", 14, {"y": 13.0}, {"y": 17.1}),
    # the fin's outer face, which the fillet's line runs on
    "fin": _plane_at("x", 10, {"z": 5.0 + 1e-6}),
    # the arm's bottom, which holds the seam at x=38.5
    "armfoot": _plane_at("z", -9.75, None, {"x": 38.5}),
}
# case -> a third-row face other than its shape's
UVFACE_CASE = {
    # the wall y=10.2 beside e50, its piece from x=38.5 on (the hexagon's)
    "issue962_e50_r0.8": _plane_at("y", 10.2, {"x": 40.0}, {"x": 50.0}),
    # the arm's bottom, which holds the seam at x=38.5
    "issue962_fillet_r0.8": _plane_at("z", -9.75, None, {"x": 38.5}),
}
NAMES = {
    "boxcyl": "10 box with a 3/4 cylinder r2 at a corner (#523, Part Connect)",
    "slotwall": "10x10x14 block, a slot down to z=5, the slot's wall two faces split 0.7 "
                "from the outer face",
    "p962": "#962's Pocket002, the Fillet's input (a PartDesign body without Refine)",
    "fin": "10x10x5 block with a 10x3x9 fin padded flush with its side x=10",
    "armfoot": "an arm fused to an 11.5x17x6 block, their bottoms two faces split along x=38.5",
}
# shape -> the whole shape's view: (center, height)
VIEW = {
    "boxcyl": ((4, 4, 5), 19.0),  # the box and cylinder span -2..10 on each axis
    "slotwall": ((5, 5, 7), 20.0),
    "p962": ((16, 4, 11), 62.0),
    "fin": ((5, 5, 7), 20.0),
    "armfoot": ((30, 9, -6.75), 42.0),
}
# shape -> (what the third row's face is, whether to draw u gridlines at pi/2)
UVLABEL = {
    "boxcyl": ("the cylinder face", True),
    "slotwall": ("the slot's floor", False),
    "p962": ("the face at the fillet's end", False),
    "fin": ("the fin's outer face", False),
    "armfoot": ("the arm's bottom", False),
}
_FOOT = [((38.5, 27.2, -3.75), (38.5, 27.2, -9.75)), ((38.5, 10.2, -3.75), (38.5, 10.2, -9.75))]
_FILLET962 = [((50, 10.2, -3.75), (50, 10.2, 22)), ((50, 13.7, 22), (50, 13.7, 14)),
              ((50, 17.1, 14), (50, 17.1, 22)), ((50, 20.3, 22), (50, 20.3, 14)),
              ((50, 23.7, 14), (50, 23.7, 22)), ((50, 27.2, -3.75), (50, 27.2, 22)),
              ((38.5, 27.2, -3.75), (38.5, 27.2, 14)), ((38.5, 10.2, -3.75), (38.5, 10.2, 14))] + _FOOT + [
              ((10.898402, 16.470802, -3.75), (10.898402, 16.470802, -9.75)),
              ((17.356038, -9.42499, -3.75), (17.356038, -9.42499, -9.75))]
# name -> (stage, shape, the edge by its ends -- or a list of them --, radius, expected volume,
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
    "slot_split_wall_r0.8": ("s962", "slotwall", ((10, 3, 5), (10, 3, 14)), 0.8, 1128.7639,
                             (10, 3, 5), (1, 0.9, 0.8)),
    "issue962_e101_r0.8": ("s962", "p962", ((50, 13.7, 22), (50, 13.7, 14)), 0.8, 11580.7739,
                           (50, 13.7, 14), (1, 0.9, 0.8)),
    "issue962_e50_r0.8": ("s962", "p962", ((38.5, 10.2, -3.75), (38.5, 10.2, 14)), 0.8,
                          11581.7861, (38.5, 10.2, 14), (-1, -1, 0.6)),
    "fin_on_block_r0.8": ("s962", "fin", ((10, 3, 5), (10, 3, 14)), 0.8, 768.7639,
                          (10, 3, 5), (1, 0.9, 0.8)),
    "arm_on_block_foot_r0.8": ("s962b", "armfoot", _FOOT, 0.8, 4603.6145,
                               (38.5, 10.2, -9.75), (-0.3, -1, -0.8)),
    "issue962_fillet_r0.8": ("s962b", "p962", _FILLET962, 0.8, 11552.7830,
                             (38.5, 10.2, -9.75), (-0.3, -1, -0.8)),
}
# multi-edge case -> what its edges are, for the picture's heading
EDGES = {
    "arm_on_block_foot_r0.8": "the two at the foot, x=38.5, y=27.2 and y=10.2",
    "issue962_fillet_r0.8": "the Fillet's twelve; zoomed on edge 49's foot",
}


def edge_between(shape, a, b):
    want = sorted([tuple(round(c, 6) for c in a), tuple(round(c, 6) for c in b)])
    for i, e in enumerate(shape.Edges):
        if sorted(tuple(round(c, 6) for c in v.Point) for v in e.Vertexes) == want:
            return e
    raise ValueError("no edge between %s and %s" % (a, b))
