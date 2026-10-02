# Hand values for the sweep's runs that neither upstream nor a suite case
# confirms (README.md, "The sweep"). Thickness 1 throughout.
#
# The rules they are worked from, each settled by an earlier fix or ruling
# (FreeCAD docs/TransactionLog.md sec 27.89 to 27.105):
#  - the skin is closed in the removed face's own surface, extended;
#  - the Intersection join (j2) meets offsets in sharp corners, the Arc join
#    (j0) rounds them where the offsets part -- outward at a convex edge of
#    the solid, inward at a concave one -- and is sharp where they cross;
#  - faces left in pieces give one solid a piece, a pocket's piece floating
#    in the cavity where nothing joins it to the walls.
import os
from math import pi

EDGE = 1 - pi / 4          # a sharp unit corner less a quarter disc, per unit length
CORNER = 1 - pi / 6        # a unit cube less an eighth of a ball
# Where two rounded edges end on a crossing (sharp) one, each stops a
# thickness short and the two quarter tubes overlap: counted whole, the two
# lose 2 EDGE in the unit cube there, where they lose 1 - (pi / 2 - 2 / 3).
MIXED = 2 * EDGE - (1 - (pi / 2 - 2. / 3))


def ring(r, toward):
    """A sharp ring corner on a circle of radius r less the quarter torus
    that rounds it, bulging to larger (+1) or smaller (-1) radius."""
    sharp = pi * abs((r + toward) ** 2 - r ** 2)
    return sharp - (pi / 4) * 2 * pi * (r + toward * 4 / (3 * pi))


def frustum(r1, r2, h=1.0):
    return pi * h * (r1 * r1 + r1 * r2 + r2 * r2) / 3


BOX = 104 * EDGE + 8 * CORNER      # a box 10 x 10 x 6 grown sharp, less grown round
BOX5 = 100 * EDGE + 8 * CORNER     # 10 x 10 x 5

# (shape, face, direction): {join: volume}; both intersection settings.
HAND = {
    # cone(5, 2, 8), its side removed: two discs closed in the cone, extended.
    ("cone", 1, +1): {0: frustum(5.375, 5) + frustum(2, 1.625), 2: frustum(5.375, 5) + frustum(2, 1.625)},
    ("cone", 1, -1): {0: frustum(5, 4.625) + frustum(2.375, 2), 2: frustum(5, 4.625) + frustum(2.375, 2)},
    # The elliptic pad, its side removed: two elliptic discs.
    ("ellipse", 1, +1): {0: 2 * pi * 50, 2: 2 * pi * 50},
    ("ellipse", 1, -1): {0: 2 * pi * 50, 2: 2 * pi * 50},
    # The L box, a notch wall removed, inward: the cavity runs to the removed
    # wall, 3 x (8 x 3 + 4 x 5).
    ("lbox", 7, -1): {0: 375 - 132, 2: 375 - 132},
    ("lbox", 8, -1): {0: 375 - 132, 2: 375 - 132},
    # The T, a post wall removed. Inward the cavity is the bar's 10 x 2 x 2
    # and the post's 3 x 2 x 6; the Arc join rounds the other post wall's
    # foot over the 2 between the side skins. Outward the post grows from the
    # removed wall's plane: 504 + 180 - 288; rounded, 96 of edge, 10 corners,
    # two mixed ends.
    ("tshape", 3, -1): {0: 212 - 2 * EDGE},
    ("tshape", 5, -1): {0: 212 - 2 * EDGE},
    ("tshape", 3, +1): {0: 396 - 96 * EDGE - 10 * CORNER + 2 * MIXED, 2: 396},
    ("tshape", 5, +1): {0: 396 - 96 * EDGE - 10 * CORNER + 2 * MIXED, 2: 396},
    # The T, its front or back removed, outward: the profile grown, 5 deep,
    # 600 - 288; rounded, 6 edges of 4 and the far profile's 44, 6 corners,
    # two mixed ends.
    ("tshape", 8, +1): {0: 312 - 68 * EDGE - 6 * CORNER + 2 * MIXED},
    ("tshape", 9, +1): {0: 312 - 68 * EDGE - 6 * CORNER + 2 * MIXED},
    # The box with a blind hole (radius 2, 3 deep). Wall removed: inward the
    # floor's disc floats in the cavity, 308 - 4 pi + 4 pi; outward 508 sharp.
    # Floor removed, outward: 500 + 12 pi - 4 pi sharp, the hole's rim rounded too.
    ("boxhole2", 7, -1): {0: 308, 2: 308},
    ("boxhole2", 7, +1): {0: 508 - BOX5},
    ("boxhole2", 8, +1): {0: 508 + 8 * pi - BOX5 - ring(2, -1)},
    # The box with a pocket 5 x 5 x 3. Top removed: the cup 408 and the pot
    # 57 outward; 280 and 121 inward with the Intersection join.
    ("pocket", 3, +1): {0: 465 - (64 * EDGE + 4 * CORNER)},
    ("pocket", 3, -1): {2: 401},
    # A pocket wall removed. Inward 395 sharp; rounded, the pocket's two
    # upright edges over 2, three floor edges of 5, two corners. Outward 591
    # sharp; rounded, the box, three rim edges of 5 and two mixed ends.
    ("pocket", 7, -1): {0: 395 - 19 * EDGE - 2 * CORNER},
    ("pocket", 7, +1): {0: 591 - BOX - 15 * EDGE + 2 * MIXED},
    # The floor removed: 367 inward, four upright edges over 2; 591 outward,
    # four rim edges of 5 and four mixed ends.
    ("pocket", 11, -1): {0: 367 - 8 * EDGE},
    ("pocket", 11, +1): {0: 591 - BOX - 20 * EDGE + 4 * MIXED},
    # The cylinder with a boss (6 x 4 under 3 x 4). Side removed: the bottom
    # disc 36 pi, and shoulder, boss wall and top: 64 pi outward, 56 pi inward.
    ("cylboss", 1, +1): {0: 100 * pi - ring(3, +1), 2: 100 * pi},
    ("cylboss", 1, -1): {0: 92 * pi - ring(3, -1), 2: 92 * pi},
    # Shoulder removed, outward: the cup 101 pi and the boss's cap 44 pi.
    ("cylboss", 2, +1): {0: 145 * pi - ring(6, +1) - ring(3, +1), 2: 145 * pi},
    # Boss wall removed, inward: 85 pi and the boss top's disc 9 pi.
    ("cylboss", 4, -1): {0: 94 * pi, 2: 94 * pi},
    # The cylinder with a pocket (6 x 6, pocket 3 x 3). Side removed: the
    # bottom disc 36 pi and top, pocket wall and floor.
    ("cylpocket", 1, +1): {0: 87 * pi - ring(3, -1), 2: 87 * pi},
    ("cylpocket", 1, -1): {0: 93 * pi - ring(3, +1)},
    # Top removed: the cup and the pot.
    ("cylpocket", 2, +1): {0: 146 * pi - ring(6, +1), 2: 146 * pi},
    ("cylpocket", 2, -1): {0: 128 * pi - ring(3, +1), 2: 128 * pi},
    # Pocket wall removed: the floor's disc floats inward, 107 pi + 9 pi;
    # outward 167 pi + 9 pi, both outer rims rounded.
    ("cylpocket", 4, -1): {0: 116 * pi, 2: 116 * pi},
    ("cylpocket", 4, +1): {0: 176 * pi - 2 * ring(6, +1), 2: 176 * pi},
    # Floor removed: 189 pi - 68 pi inward; outward 392 pi - 189 pi - 16 pi,
    # both outer rims and the pocket's rounded.
    ("cylpocket", 5, -1): {0: 121 * pi, 2: 121 * pi},
    ("cylpocket", 5, +1): {0: 187 * pi - 2 * ring(6, +1) - ring(3, -1), 2: 187 * pi},
}
# Worked in the suite's own cases too (README.md), where the reference is
# written another way -- a list of pieces, or the other intersection setting.
HAND.update({
    # The cylinder, its side removed: two discs of 16 pi.
    ("cyl", 1, +1): {0: 32 * pi, 2: 32 * pi},
    ("cyl", 1, -1): {0: 32 * pi, 2: 32 * pi},
    # The boss's shoulder removed, inward: the cup 69 pi, the cap 24 pi.
    ("cylboss", 2, -1): {0: 93 * pi, 2: 93 * pi},
    # The boss's wall removed, outward: 150 pi sharp, two rims rounded.
    ("cylboss", 4, +1): {0: 150 * pi - 2 * ring(6, +1)},
    # The blind hole's top removed (sec 27.103).
    ("boxhole2", 3, +1): {0: 300 + 15 * pi + 2 * pi / 3 + 10 * pi, 2: 364 + 10 * pi},
    ("boxhole2", 3, -1): {0: 244 + 19 * pi + (pi * pi / 2) * (2 + 4 / (3 * pi)), 2: 244 + 24 * pi},
    # The L box, a notch wall removed, outward: 423 sharp; 90 of edge, 8 corners.
    ("lbox", 7, +1): {0: 423 - 90 * EDGE - 8 * CORNER},
    ("lbox", 8, +1): {0: 423 - 90 * EDGE - 8 * CORNER},
})
HAND[("cylpocket", 1, -1)][2] = 93 * pi
HAND[("pocket", 3, +1)][2] = 465
# The pocket's top removed, inward, rounded: the pot's four upright edges of
# 3, four floor edges of 5, four corners.
HAND[("pocket", 3, -1)][0] = 280 + 121 - 32 * EDGE - 4 * CORNER
for a, b in ((7, 8), (7, 9), (7, 10)):
    for d in (+1, -1):
        HAND[("pocket", b, d)] = HAND[("pocket", a, d)]


# The pole solids -- a ball cut by its equator and by planes through its axis
# -- are worked numerically, by polehand.py, which says how; polehand.txt is
# its table.
POLES = {}
for _line in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "polehand.txt")):
    _tag, _vol = _line.split()
    POLES[_tag] = float(_vol)


def hand(tag):
    """The hand value of a run's tag (shape_fN_+1_n_j0), or None."""
    shape, face, direction, _, join = tag.split("_")
    v = HAND.get((shape, int(face[1:]), int(direction)))
    if v is None:
        return POLES.get("%s_%s_%s_%s" % (shape, face, direction, join))
    return v.get(int(join[1:]))
