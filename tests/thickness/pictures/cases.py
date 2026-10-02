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


# stage -> (the fork commit just before the fix, the section of models/Thickness.md)
STAGES = {
    "s2": ("050c66d58e", "2"),
    "s3": ("2a15bbb52b", "3"),
    "s4": ("cbc47c91a7", "4"),
    "s5": ("51a0cc9b38", "5"),
    "s6": ("e5e6d02f70", "6"),
    "s8": ("9a14fb9db5", "8"),
    "s9": ("c4be5548f9", "9"),
    "s10": ("457b652f42", "10"),
    "s11": ("b21dabe1c4", "11"),
    "s12": ("d4d2fe0719", "12"),
    "s13": ("d8d480ef65", "13"),
    "s14": ("839a3e4606", "14"),
    "s15": ("f723999a15", "15"),
    "s16": ("dc3a7b85d1", "16"),
    "s17": ("b6404e079a", "17"),
    "s18": ("705644654e", "18"),
    "s19": ("e63b3f86cc", "19"),
    "s20": ("ffc01775a7", "20"),
    "s21a": ("76691d6f87", "21"),
    "s21b": ("b45cff22c5", "21"),
    "s21c": ("20372de2b3", "21"),
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
def _placement():
    return App.Placement(V(3, 4, 5), App.Rotation(V(1, 2, 3), 40))


def _placed(shape):
    # Under a location: an object's placement.
    shape.Placement = _placement()
    return shape


def _turned(shape):
    # The geometry itself turned and moved, no location on it.
    shape.transformShape(_placement().Matrix, True)
    return shape


def _halfball(cut=None):
    # Half a ball cut through both its poles; with cut, its sphere in two
    # faces already, by a general fuse with a meridian or the equator.
    b = Part.makeSphere(5, V(), V(0, 0, 1), -90, 90, 180)
    if cut == "meridian":
        e = Part.Arc(V(0, 0, -5), V(0, 5, 0), V(0, 0, 5)).toShape()
    elif cut == "equator":
        e = Part.ArcOfCircle(Part.Circle(V(), V(0, 0, 1), 5), 0, math.pi).toShape()
    else:
        return b
    return b.generalFuse([e])[0].Solids[0]


def _dome2(line=False):
    # A dome on the y axis, its rim in two arcs: fused with the vertex across
    # from its seam's end, or -- line -- with the line between them, which
    # cuts its flat in two halves as well.
    d = Part.makeSphere(5, V(), V(0, 1, 0), 0, 90, 360)
    rim = [e for e in d.Edges if not e.Degenerated and abs(e.Length - 10 * math.pi) < 1e-6][0]
    p = rim.Vertexes[0].Point
    tool = Part.makeLine(p, p * -1) if line else Part.Vertex(p * -1)
    return d.generalFuse([tool])[0].Solids[0]


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
    "halfdome": lambda: Part.makeSphere(5, V(), V(0, 0, 1), 0, 90, 180),
    "splitbox": lambda: Part.makeBox(4, 8, 6).fuse(Part.makeBox(6, 8, 6, V(4, 0, 0))),
    "eighth": lambda: Part.makeSphere(5, V(), V(0, 0, 1), 0, 90, 90),
    "dome120": lambda: Part.makeSphere(5, V(), V(0, 0, 1), 0, 90, 120),
    "dome270": lambda: Part.makeSphere(5, V(), V(0, 0, 1), 0, 90, 270),
    "lune": lambda: Part.makeSphere(5, V(), V(0, 0, 1), -90, 90, 90),
    "placed_dome": lambda: _placed(Part.makeSphere(5, V(), V(0, 0, 1), 0, 90, 360)),
    "placed_cone": lambda: _placed(Part.makeCone(0, 4, 6)),
    "placed_coneup": lambda: _placed(Part.makeCone(4, 0, 6)),
    "placed_halfdome": lambda: _placed(Part.makeSphere(5, V(), V(0, 0, 1), 0, 90, 180)),
    "turned_filletbox": lambda: _turned(_fillet(2)),
    "turned_halfcap": lambda: _turned(Part.makeSphere(5, V(), V(0, 0, 1), 30, 90, 180)),
    "turned_bullet": lambda: _turned(Part.makeCylinder(5, 4, V(0, 0, -4)).fuse(
        Part.makeSphere(5, V(), V(0, 0, 1), 0, 90, 360)).removeSplitter()),
    "halfball": lambda: _halfball(),
    "luneball": lambda: _halfball("meridian"),
    "eqball": lambda: _halfball("equator"),
    "dome2": lambda: _dome2(),
    "dome2f": lambda: _dome2(True),
    "halfball1": lambda: _halfball().removeSplitter(),
    "cutball": lambda: Part.makeSphere(5).cut(Part.makeBox(20, 20, 20, V(-10, -20, -10))),
}
# name: (stage, shape, face or list of faces removed, value, inter, join, reference volume)
CASES = {
    "cyl_bottom_in": ("s2", "cyl", 2, -1, False, 0, 468.0973),
    "hole_bottom_in": ("s2", "ann", 2, -1, False, 0, 461.8141),
    "hole_outer_out": ("s2", "ann", 1, +1, False, 0, 241.7451),
    "hole_outer_in": ("s2", "ann", 1, -1, False, 0, 257.6106),
    "hole_inner_out": ("s2", "ann", 4, +1, False, 0, 531.0589),
    "hole_inner_in": ("s2", "ann", 4, -1, False, 0, 358.1416),
    "ellipse_bottom_out": ("s2", "ell", 2, +1, False, 0, 608.6622),
    "ellipse_bottom_in": ("s2", "ell", 2, -1, False, 0, 473.2070),
    "ellipse_top_out": ("s2", "ell", 3, +1, False, 0, 608.6631),
    "ellipse_top_in": ("s2", "ell", 3, -1, False, 0, 473.2073),
    "pocket_bottom_out": ("s2", "cylpocket", 3, +1, False, 0, 433.9707),
    "pocket_bottom_in": ("s2", "cylpocket", 3, -1, False, 0, 346.7660),
    "boxhole_hole_out": ("s2", "boxhole", 7, +1, False, 0, 457.5959),
    "boxhole_hole_in": ("s2", "boxhole", 7, -1, False, 0, 282.8673),
    "lbox_arm_end_out": ("s3", "lbox", 4, +1, False, 0, 388.5671),
    "tshape_arm_end_out": ("s3", "tshape", 1, +1, False, 0, 372.9204),
    "pocketbox_wall_out": ("s3", "pocketbox", 7, +1, False, 0, None),
    "tshape_bar_top_out": ("s3", "tshape", 2, +1, False, 0, 382.4425),
    "tshape_bar_top_in": ("s4", "tshape", 2, -1, False, 0, 215.5708),
    "lbox_top_inter_join_out": ("s4", "lbox", 3, +1, True, 2, 339.0),
    "lbox_top_inter_join_in": ("s4", "lbox", 3, -1, True, 2, 219.0),
    "tshape_back_inter_join_out": ("s4", "tshape", 8, +1, True, 2, 312.0),
    "tshape_back_inter_join_in": ("s4", "tshape", 8, -1, True, 2, 192.0),
    "conehole_bottom_inter_in": ("s4", "conehole", 3, -1, True, 0, 307.1946),
    "tshape_bar_top_right_join_in": ("s5", "tshape", 7, -1, False, 2, 216.0),
    "pocket_floor_join_out": ("s5", "pocket5", 11, +1, False, 2, 591.0),
    "pocket_floor_join_in": ("s5", "pocket5", 11, -1, False, 2, 367.0),
    "blindhole_floor_join_out": ("s5", "blindhole", 8, +1, False, 2, 533.1327),
    "blindhole_floor_join_in": ("s5", "blindhole", 8, -1, False, 2, 326.8496),
    "lbox_notch_wall_inter_join_out": ("s5", "lbox", 7, +1, True, 2, 423.0),
    "tshape_post_wall_inter_join_in": ("s5", "tshape", 3, -1, True, 2, 212.0),
    "pocket_wall_inter_join_in": ("s5", "pocket5", 7, -1, True, 2, 395.0),
    "filletbox_end_in": ("s6", "filletbox", 1, -1, False, 0, 261.1150),
    "filletbox_side_in": ("s6", "filletbox", 6, -1, False, 0, 253.1150),
    "filletbox_end_out": ("s6", "filletbox", 1, +1, False, 0, 403.9604),
    "filletbox_side_out": ("s6", "filletbox", 6, +1, False, 0, 388.8188),
    "filletbox_fillet_in": ("s8", "filletbox", 3, -1, False, 0, 267.0193),
    "filletbox_fillet_out": ("s8", "filletbox", 3, +1, False, 0, 405.9034),
    "filletbox_end_join_in": ("s9", "filletbox", 1, -1, False, 2, 262.8319),
    "filletbox_end_join_out": ("s9", "filletbox", 1, +1, False, 2, 422.7964),
    "filletbox_side_join_in": ("s9", "filletbox", 6, -1, False, 2, 254.8319),
    "filletbox_side_join_out": ("s9", "filletbox", 6, +1, False, 2, 406.7964),
    "filletbox_fillet_join_in": ("s9", "filletbox", 3, -1, False, 2, 268.7129),
    "filletbox_fillet_join_out": ("s9", "filletbox", 3, +1, False, 2, 424.7690),
    "filletbox25_fillet_join_in": ("s9", "filletbox25", 3, -1, False, 2, 258.4221),
    "filletbox25_fillet_join_out": ("s9", "filletbox25", 3, +1, False, 2, 407.4611),
    "conehole_top_in": ("s10", "conehole", 2, -1, False, 0, 358.3146),
    "conehole_top_inter_in": ("s10", "conehole", 2, -1, True, 0, 358.3146),
    "conehole_bottom_in": ("s10", "conehole", 3, -1, False, 0, 307.1946),
    "cyl_side_out": ("s11", "cyl", 1, +1, False, 0, 100.5310),
    "cyl_side_in": ("s11", "cyl", 1, -1, False, 0, 100.5310),
    "short_cyl_side_in": ("s11", "shortcyl", 1, -1, False, 0, 75.3982),
    "short_box_bottom_in": ("s11", "shortbox", [1, 2, 3, 4, 6], -1, False, 0, 100.0),
    "cyl_cap_alone_join_out": ("s12", "cyl", [1, 3], +1, False, 2, 50.2655),
    "box_bottom_alone_join_in": ("s12", "box6", [1, 2, 3, 4, 6], -1, False, 2, 100.0),
    "cyl_side_join_out": ("s12", "cyl", 1, +1, False, 2, 100.5310),
    "blind_top_out": ("s13", "blindhole", 3, +1, False, 0, 380.6342),
    "blind_top_in": ("s13", "blindhole", 3, -1, False, 0, 315.6543),
    "pocket_top_in": ("s13", "pocket5", 3, -1, False, 0, 392.2271),
    "cylpocket_top_join_out": ("s13", "cylpocket", 2, +1, False, 2, 458.6725),
    "blind_top_in_inter_join": ("s13", "blindhole", 3, -1, True, 2, 319.3982),
    "blind_wall_out_join": ("s13", "blindhole", 7, +1, False, 2, 508.0),
    "blind_pot_shell_out": ("s13", "blindpot", 3, +1, False, 0, 31.4159),
    "blind_pot_shell_in": ("s13", "blindpot", 3, -1, False, 0, 71.6543),
    "pocket_top_out_inter_join": ("s13", "pocket5", 3, +1, True, 2, 465.0),
    "cylboss_shoulder_in": ("s13", "cylboss", 2, -1, False, 0, 292.1681),
    "sector_outer_arc_input": ("s14", "sector", [3, 5, 6], +1, False, 0, 9405.6558),
    "sphere_face_refused_in": ("s15", "sphere", 1, -1, False, 0, None),
    "box_all_faces_refused_in": ("s15", "box6", [1, 2, 3, 4, 5, 6], -1, False, 0, None),
    "dome_flat_out": ("s16", "dome", 2, +0.5, False, 0, 86.6556),
    "dome_flat_in": ("s16", "dome", 2, -0.5, False, 0, 70.9476),
    "cap_flat_in": ("s16", "cap", 2, -0.5, False, 0, 33.6412),
    "cone_base_out": ("s16", "cone", 2, +0.5, False, 0, 52.4095),
    "cone_base_in": ("s17", "cone", 2, -0.5, False, 0, 38.8428),
    "cone_base_join_in": ("s17", "cone", 2, -0.5, False, 2, 38.8428),
    "halfdome_bottom_out": ("s17", "halfdome", 2, +0.5, False, 0, 66.1779),
    "halfdome_side_out": ("s17", "halfdome", 3, +0.5, False, 0, 79.7628),
    "halfdome_side_in": ("s17", "halfdome", 3, -0.5, False, 0, 58.8944),
    "halfdome_other_side_in": ("s18", "halfdome", 4, -0.5, False, 0, 58.8944),
    "splitbox_top_piece_join_out": ("s18", "splitbox", 3, +1, False, 2, 440.0),
    "splitbox_top_piece_join_in": ("s18", "splitbox", 3, -1, False, 2, 276.0),
    "splitbox_end_join_in": ("s18", "splitbox", 1, -1, False, 2, 264.0),
    "halfdome_side_join_in": ("s18", "halfdome", 3, -0.5, False, 2, 59.1071),
    "halfdome_bottom_join_out": ("s18", "halfdome", 2, +0.5, False, 2, 67.0206),
    "halfdome_side_join_out": ("s18", "halfdome", 3, +0.5, False, 2, 81.7345),
    "eighth_ball_bottom_join_out": ("s18", "eighth", 2, +0.5, False, 2, 46.7279),
    "dome270_bottom_join_out": ("s19", "dome270", 2, +0.5, False, 2, 87.3133),
    "dome270_side_join_out": ("s19", "dome270", 3, +0.5, False, 2, 113.7486),
    "dome270_side_join_in": ("s19", "dome270", 3, -0.5, False, 2, 83.7681),
    "dome270_side_in": ("s19", "dome270", 3, -0.5, False, 0, 83.7681),
    "dome270_side_join_out_thick": ("s19", "dome270", 3, +1.0, False, 2, 260.9367),
    "eighth_ball_bottom_out": ("s19", "eighth", 2, +0.5, False, 0, 45.5612),
    "dome120_bottom_out": ("s19", "dome120", 2, +0.5, False, 0, 52.4334),
    "dome120_side_in": ("s19", "dome120", 3, -0.5, False, 0, 41.2951),
    "eighth_ball_sphere_out": ("s19", "eighth", 1, +0.5, False, 0, 32.3576),
    "eighth_ball_sphere_join_out": ("s19", "eighth", 1, +0.5, False, 2, 33.2167),
    "dome120_sphere_join_out": ("s19", "dome120", 1, +0.5, False, 2, 35.8993),
    "halfdome_sphere_out": ("s19", "halfdome", 1, +0.5, False, 0, 41.0976),
    "halfdome_sphere_join_out": ("s19", "halfdome", 1, +0.5, False, 2, 41.6307),
    "lune_sphere_out": ("s19", "lune", 1, +0.5, False, 0, 41.0976),
    "lune_sphere_in": ("s19", "lune", 1, -0.5, False, 0, 36.6474),
    "dome270_sphere_join_out": ("s19", "dome270", 1, +0.5, False, 2, 50.0446),
    "dome270_sphere_out_thick": ("s19", "dome270", 1, +1.0, False, 0, 99.0413),
    "placed_dome_flat_out": ("s20", "placed_dome", 2, +0.5, False, 0, 86.6556),
    "placed_cone_base_join_out": ("s20", "placed_cone", 2, +0.5, False, 2, 52.4563),
    "placed_coneup_base_out": ("s20", "placed_coneup", 2, +0.5, False, 0, 52.4095),
    "placed_halfdome_bottom_join_out": ("s20", "placed_halfdome", 2, +0.5, False, 2, 67.0206),
    "turned_filletbox_out": ("s20", "turned_filletbox", 4, +1, False, 0, 405.9034),
    "turned_halfcap_sphere_in": ("s20", "turned_halfcap", 1, -0.5, False, 0, 19.2316),
    "turned_halfcap_sphere_join_in": ("s20", "turned_halfcap", 1, -0.5, False, 2, 19.2316),
    "turned_bullet_side_join_out": ("s20", "turned_bullet", 1, +0.5, False, 2, 100.7315),
    "halfball_flat_join_out": ("s20", "halfball", 2, +0.5, False, 2, 113.0909),
    "halfball_flat_join_in": ("s20", "halfball", 2, -0.5, False, 2, 89.0272),
    "halfball_sphere_out": ("s20", "halfball", 1, +0.5, False, 0, 39.1390),
    "halfball_sphere_join_out": ("s20", "halfball", 1, +0.5, False, 2, 39.1390),
    "luneball_spheres_out": ("s20", "luneball", [1, 2], +0.5, False, 0, 39.1390),
    "eqball_flat_join_out": ("s20", "eqball", 3, +0.5, False, 2, 113.0909),
    "eqball_flat_join_in": ("s20", "eqball", 3, -0.5, False, 2, 89.0272),
    "dome2_flat_out": ("s21a", "dome2", 2, +0.5, False, 0, 86.6556),
    "dome2_flat_in": ("s21a", "dome2", 2, -0.5, False, 0, 70.9476),
    "dome2_sphere_out": ("s21a", "dome2", 1, +0.5, False, 0, 39.1390),
    "dome2_sphere_join_in": ("s21a", "dome2", 1, -0.5, False, 2, 39.1390),
    "dome2f_sphere_out": ("s21a", "dome2f", 1, +0.5, False, 0, 39.1390),
    "halfball1_sphere_out": ("s21a", "halfball1", 1, +0.5, False, 0, 39.1390),
    "halfball1_sphere_in": ("s21a", "halfball1", 1, -0.5, False, 0, 39.1390),
    "halfball1_disc_join_out": ("s21b", "halfball1", 2, +0.5, False, 2, 86.6556),
    "halfball1_disc_join_in": ("s21b", "halfball1", 2, -0.5, False, 2, 70.9476),
    "halfball1_sphere_join_out": ("s21b", "halfball1", 1, +0.5, False, 2, 39.1390),
    "cutball_disc_join_out": ("s21b", "cutball", 2, +0.5, False, 2, 86.6556),
    "cutball_sphere_out": ("s21b", "cutball", 1, +0.5, False, 0, 39.1390),
    "luneball_flat_join_out": ("s21c", "luneball", 3, +0.5, False, 2, 113.0909),
    "luneball_flat_join_in": ("s21c", "luneball", 3, -0.5, False, 2, 89.0272),
    "eqball_spheres_out": ("s21c", "eqball", [1, 2], +0.5, False, 0, 39.1390),
    "eqball_spheres_join_out": ("s21c", "eqball", [1, 2], +0.5, False, 2, 39.1390),
}
# Cases whose right result is more than one shell: the holed cone's top, its
# cavity sealed below the removed face (a skin and a void).
SHELLS = {
    "conehole_top_in": 2,
    "conehole_top_inter_in": 2,
}
# Cases whose right result is several solids, a compound: the faces that
# stay fall apart into pieces once the removed ones are gone (sec 11).
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
    "turned_bullet_side_join_out": 2,
}
# Cases pictured for what the thickness leaves of its input: the panels show
# the input after the call, and the reference volume is the input's own
# (sec 14: the thickness was right and left its input inside out).
INPUT = {
    "sector_outer_arc_input",
}
# Cases whose right outcome is a refusal: every face is removed, no face
# stays to be thickened (sec 15). The panels show the input.
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
