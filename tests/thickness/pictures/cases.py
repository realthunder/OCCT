# The cases pictured in models/pictures: the suite's (run_tests.py) that a fix
# turned from failing to passing, each with the stage of the fix -- the fork
# commit just before it (STAGES) is the "before" column.
import math

try:
    import FreeCAD as App, Part
except ImportError:  # plain Python reads the tables only
    App = Part = None


def V(*a):
    return App.Vector(*a)


# stage -> (the fork commit just before the fix, FreeCAD docs/TransactionLog.md section)
STAGES = {
    "s89": ("050c66d58e", "27.89"),
    "s90": ("2a15bbb52b", "27.90"),
    "s91": ("cbc47c91a7", "27.91"),
    "s92": ("51a0cc9b38", "27.92"),
    "s93": ("e5e6d02f70", "27.93"),
    "s95": ("9a14fb9db5", "27.95"),
    "s96": ("c4be5548f9", "27.96"),
    "s100": ("457b652f42", "27.100"),
    "s101": ("b21dabe1c4", "27.101"),
    "s102": ("d4d2fe0719", "27.102"),
    "s103": ("d8d480ef65", "27.103"),
    "s104": ("839a3e4606", "27.104"),
    "s105": ("f723999a15", "27.105"),
    "s106": ("dc3a7b85d1", "27.106"),
}
# Upstream: the eleven files the fix chain touches, at the fork's base.
UPSTREAM = ("91be8c4c71", [
    "src/ModelingAlgorithms/TKBool/BRepAlgo/" + f for f in (
        "BRepAlgo_Loop.cxx", "BRepAlgo_Loop.hxx", "BRepAlgo_AsDes.cxx",
        "BRepAlgo_FaceRestrictor.cxx", "BRepAlgo_Image.cxx")] + [
    "src/ModelingAlgorithms/TKOffset/BRepOffset/" + f for f in (
        "BRepOffset_Tool.cxx", "BRepOffset_Inter2d.cxx", "BRepOffset_Inter2d.hxx",
        "BRepOffset_Inter3d.cxx", "BRepOffset_MakeLoops.cxx", "BRepOffset_MakeOffset.cxx")])
def _blindpot():
    # The blind hole's wall and floor and the top beside them, an open shell:
    # Face3 is the top.
    b = SHAPES["blindhole"]()
    return Part.Shell([b.Faces[i - 1] for i in (7, 8, 3)])


def _sector(past):
    # A ring sector straddling angle 0, padded 27, its arcs' parameters
    # running from past - a to past + a: Face3 is the outer arc's.
    o = Part.ArcOfCircle(Part.Circle(V(0, 0, 0), V(0, 0, 1), 59.3), past - 0.29, past + 0.29).toShape()
    i = Part.ArcOfCircle(Part.Circle(V(0, 0, 0), V(0, 0, 1), 46.85), past - 0.23, past + 0.23).toShape()
    sides = [Part.makeLine(i.Vertexes[k].Point, o.Vertexes[k].Point) for k in (0, 1)]
    return Part.Face(Part.Wire(Part.__sortEdges__([o, sides[0], i, sides[1]]))).extrude(V(0, 0, 27))


def _fillet(r):
    b = Part.makeBox(10, 8, 6)
    return b.makeFillet(r, [b.Edges[i] for i in (0, 2, 4, 6)])
SHAPES = {
    "cyl": lambda: Part.makeCylinder(4, 20),
    "ann": lambda: Part.makeCylinder(5, 10).cut(Part.makeCylinder(2, 10)),
    "ell": lambda: Part.Face(Part.Wire(Part.Ellipse(V(0, 0, 0), 10, 5).toShape())).extrude(V(0, 0, 8)),
    "cylpocket": lambda: Part.makeCylinder(6, 6).cut(Part.makeCylinder(3, 3, V(0, 0, 3))),
    "boxhole": lambda: Part.makeBox(10, 10, 5).cut(Part.makeCylinder(2, 5, V(5, 5, 0))),
    "lbox": lambda: Part.makeBox(10, 10, 5).cut(Part.makeBox(5, 5, 5, V(5, 5, 0))),
    "tshape": lambda: Part.makeBox(12, 4, 4).fuse(Part.makeBox(4, 4, 10, V(4, 0, 0))).removeSplitter(),
    "pocketbox": lambda: Part.makeBox(10, 10, 6).cut(Part.makeBox(6, 6, 3, V(2, 2, 3))),
    "conehole": lambda: Part.makeCone(6, 3, 8).cut(Part.makeCylinder(1.5, 8)),
    "pocket5": lambda: Part.makeBox(10, 10, 6).cut(Part.makeBox(5, 5, 3, V(2.5, 2.5, 3))),
    "blindhole": lambda: Part.makeBox(10, 10, 5).cut(Part.makeCylinder(2, 3, V(5, 5, 2))),
    "filletbox": lambda: _fillet(2),
    "filletbox25": lambda: _fillet(2.5),
    "shortcyl": lambda: Part.makeCylinder(4, 1.5),
    "shortbox": lambda: Part.makeBox(10, 10, 1.5),
    "box6": lambda: Part.makeBox(10, 10, 6),
    "blindpot": lambda: _blindpot(),
    "cylboss": lambda: Part.makeCylinder(6, 4).fuse(Part.makeCylinder(3, 4, V(0, 0, 4))),
    "sector": lambda: _sector(2 * math.pi),
    "sphere": lambda: Part.makeSphere(5),
    "dome": lambda: Part.makeSphere(5, V(), V(0, 0, 1), 0, 90, 360),
    "cap": lambda: Part.makeSphere(5, V(), V(0, 0, 1), 30, 90, 360),
    "cone": lambda: Part.makeCone(0, 4, 6),
}
# name: (stage, shape, face or list of faces removed, value, inter, join, reference volume)
CASES = {
    "cyl_bottom_in": ("s89", "cyl", 2, -1, False, 0, 468.0973),
    "hole_bottom_in": ("s89", "ann", 2, -1, False, 0, 461.8141),
    "hole_outer_out": ("s89", "ann", 1, +1, False, 0, 241.7451),
    "hole_outer_in": ("s89", "ann", 1, -1, False, 0, 257.6106),
    "hole_inner_out": ("s89", "ann", 4, +1, False, 0, 531.0589),
    "hole_inner_in": ("s89", "ann", 4, -1, False, 0, 358.1416),
    "ellipse_bottom_out": ("s89", "ell", 2, +1, False, 0, 608.6622),
    "ellipse_bottom_in": ("s89", "ell", 2, -1, False, 0, 473.2070),
    "ellipse_top_out": ("s89", "ell", 3, +1, False, 0, 608.6631),
    "ellipse_top_in": ("s89", "ell", 3, -1, False, 0, 473.2073),
    "pocket_bottom_out": ("s89", "cylpocket", 3, +1, False, 0, 433.9707),
    "pocket_bottom_in": ("s89", "cylpocket", 3, -1, False, 0, 346.7660),
    "boxhole_hole_out": ("s89", "boxhole", 7, +1, False, 0, 457.5959),
    "boxhole_hole_in": ("s89", "boxhole", 7, -1, False, 0, 282.8673),
    "lbox_arm_end_out": ("s90", "lbox", 4, +1, False, 0, 388.5671),
    "tshape_arm_end_out": ("s90", "tshape", 1, +1, False, 0, 372.9204),
    "pocketbox_wall_out": ("s90", "pocketbox", 7, +1, False, 0, None),
    "tshape_bar_top_out": ("s90", "tshape", 2, +1, False, 0, 382.4425),
    "tshape_bar_top_in": ("s91", "tshape", 2, -1, False, 0, 215.5708),
    "lbox_top_inter_join_out": ("s91", "lbox", 3, +1, True, 2, 339.0),
    "lbox_top_inter_join_in": ("s91", "lbox", 3, -1, True, 2, 219.0),
    "tshape_back_inter_join_out": ("s91", "tshape", 8, +1, True, 2, 312.0),
    "tshape_back_inter_join_in": ("s91", "tshape", 8, -1, True, 2, 192.0),
    "conehole_bottom_inter_in": ("s91", "conehole", 3, -1, True, 0, 307.1946),
    "tshape_bar_top_right_join_in": ("s92", "tshape", 7, -1, False, 2, 216.0),
    "pocket_floor_join_out": ("s92", "pocket5", 11, +1, False, 2, 591.0),
    "pocket_floor_join_in": ("s92", "pocket5", 11, -1, False, 2, 367.0),
    "blindhole_floor_join_out": ("s92", "blindhole", 8, +1, False, 2, 533.1327),
    "blindhole_floor_join_in": ("s92", "blindhole", 8, -1, False, 2, 326.8496),
    "lbox_notch_wall_inter_join_out": ("s92", "lbox", 7, +1, True, 2, 423.0),
    "tshape_post_wall_inter_join_in": ("s92", "tshape", 3, -1, True, 2, 212.0),
    "pocket_wall_inter_join_in": ("s92", "pocket5", 7, -1, True, 2, 395.0),
    "filletbox_end_in": ("s93", "filletbox", 1, -1, False, 0, 261.1150),
    "filletbox_side_in": ("s93", "filletbox", 6, -1, False, 0, 253.1150),
    "filletbox_end_out": ("s93", "filletbox", 1, +1, False, 0, 403.9604),
    "filletbox_side_out": ("s93", "filletbox", 6, +1, False, 0, 388.8188),
    "filletbox_fillet_in": ("s95", "filletbox", 3, -1, False, 0, 267.0193),
    "filletbox_fillet_out": ("s95", "filletbox", 3, +1, False, 0, 405.9034),
    "filletbox_end_join_in": ("s96", "filletbox", 1, -1, False, 2, 262.8319),
    "filletbox_end_join_out": ("s96", "filletbox", 1, +1, False, 2, 422.7964),
    "filletbox_side_join_in": ("s96", "filletbox", 6, -1, False, 2, 254.8319),
    "filletbox_side_join_out": ("s96", "filletbox", 6, +1, False, 2, 406.7964),
    "filletbox_fillet_join_in": ("s96", "filletbox", 3, -1, False, 2, 268.7129),
    "filletbox_fillet_join_out": ("s96", "filletbox", 3, +1, False, 2, 424.7690),
    "filletbox25_fillet_join_in": ("s96", "filletbox25", 3, -1, False, 2, 258.4221),
    "filletbox25_fillet_join_out": ("s96", "filletbox25", 3, +1, False, 2, 407.4611),
    "conehole_top_in": ("s100", "conehole", 2, -1, False, 0, 358.3146),
    "conehole_top_inter_in": ("s100", "conehole", 2, -1, True, 0, 358.3146),
    "conehole_bottom_in": ("s100", "conehole", 3, -1, False, 0, 307.1946),
    "cyl_side_out": ("s101", "cyl", 1, +1, False, 0, 100.5310),
    "cyl_side_in": ("s101", "cyl", 1, -1, False, 0, 100.5310),
    "short_cyl_side_in": ("s101", "shortcyl", 1, -1, False, 0, 75.3982),
    "short_box_bottom_in": ("s101", "shortbox", [1, 2, 3, 4, 6], -1, False, 0, 100.0),
    "cyl_cap_alone_join_out": ("s102", "cyl", [1, 3], +1, False, 2, 50.2655),
    "box_bottom_alone_join_in": ("s102", "box6", [1, 2, 3, 4, 6], -1, False, 2, 100.0),
    "cyl_side_join_out": ("s102", "cyl", 1, +1, False, 2, 100.5310),
    "blind_top_out": ("s103", "blindhole", 3, +1, False, 0, 380.6342),
    "blind_top_in": ("s103", "blindhole", 3, -1, False, 0, 315.6543),
    "pocket_top_in": ("s103", "pocket5", 3, -1, False, 0, 392.2271),
    "cylpocket_top_join_out": ("s103", "cylpocket", 2, +1, False, 2, 458.6725),
    "blind_top_in_inter_join": ("s103", "blindhole", 3, -1, True, 2, 319.3982),
    "blind_wall_out_join": ("s103", "blindhole", 7, +1, False, 2, 508.0),
    "blind_pot_shell_out": ("s103", "blindpot", 3, +1, False, 0, 31.4159),
    "blind_pot_shell_in": ("s103", "blindpot", 3, -1, False, 0, 71.6543),
    "pocket_top_out_inter_join": ("s103", "pocket5", 3, +1, True, 2, 465.0),
    "cylboss_shoulder_in": ("s103", "cylboss", 2, -1, False, 0, 292.1681),
    "sector_outer_arc_input": ("s104", "sector", [3, 5, 6], +1, False, 0, 9405.6558),
    "sphere_face_refused_in": ("s105", "sphere", 1, -1, False, 0, None),
    "box_all_faces_refused_in": ("s105", "box6", [1, 2, 3, 4, 5, 6], -1, False, 0, None),
    "dome_flat_out": ("s106", "dome", 2, +0.5, False, 0, 86.6556),
    "dome_flat_in": ("s106", "dome", 2, -0.5, False, 0, 70.9476),
    "cap_flat_in": ("s106", "cap", 2, -0.5, False, 0, 33.6412),
    "cone_base_out": ("s106", "cone", 2, +0.5, False, 0, 52.4095),
}
# Cases whose right result is more than one shell: the holed cone's top, its
# cavity sealed below the removed face (a skin and a void).
SHELLS = {
    "conehole_top_in": 2,
    "conehole_top_inter_in": 2,
}
# Cases whose right result is several solids, a compound: the faces that
# stay fall apart into pieces once the removed ones are gone (sec 27.101).
SOLIDS = {
    "cyl_side_out": 2,
    "cyl_side_in": 2,
    "cyl_side_join_out": 2,
    "blind_top_out": 2,
    "pocket_top_in": 2,
    "cylpocket_top_join_out": 2,
    "blind_wall_out_join": 2,
    "pocket_top_out_inter_join": 2,
    "cylboss_shoulder_in": 2,
}
# Cases pictured for what the thickness leaves of its input: the panels show
# the input after the call, and the reference volume is the input's own
# (sec 27.104: the thickness was right and left its input inside out).
INPUT = {
    "sector_outer_arc_input",
}
# Cases whose right outcome is a refusal: every face is removed, no face
# stays to be thickened (sec 27.105). The panels show the input.
REFUSED = {
    "sphere_face_refused_in",
    "box_all_faces_refused_in",
}
# The captured user models; upstream 8.0.1 passes them all, so not pictured.
DOCS = {
    "issue1_ellipse_thickness": ("issue1_ellipse_thickness.FCStd", "Thickness", 3598.2930),
    "issue3_pad_thickness": ("issue3_pad_thickness.FCStd", "Thickness", 1241.0718),
    "issue4_revolution_thickness": ("issue4_revolution_thickness.FCStd", "Thickness", 431.4454),
}
