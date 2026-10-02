# polehand.py [thickness]  -- writes the table polehand.txt holds (thickness 1)
#
# Hand values for the sweep's pole solids: a ball of R cut by its equator
# (z > 0, the domes) or left whole from pole to pole, and by two planes
# through its axis that leave a wedge of angle A in (x, y), from the +x half
# plane round to angle A. Numeric: a
# midpoint rule over (x, y), the extent in z worked exactly at each point.
#
# hand.py reads the table, so that check.py needs no numpy; this script is how
# it is made:  python polehand.py > polehand.txt  (FreeCAD's Python has numpy).
# At a thickness of 0.5 it gives the 76 volumes the suite and the probes of
# FreeCAD docs/TransactionLog.md sec 27.106 to 27.109 had settled one by one.
#
# The rules (tests/thickness/models/Thickness.md):
#  - the skin is closed in the removed face's own surface, extended;
#  - j2 (Intersection) meets the kept faces' offsets in sharp corners;
#  - j0 (Arc) rounds them where they part: outward at a convex edge, inward
#    at a concave one (the axis of the three-quarter dome);
#  - a removed side in the plane of the kept one (half a dome) is closed by a
#    wall a thickness past the edge between them, and with the Arc join by a
#    tube round that edge instead.
import math
import sys

import numpy as np

R = 5.0
N = 4800
SPAN = 6.0                      # h = 0.0025: the offsets and the radii fall on cell edges
h = 2 * SPAN / N
c = -SPAN + (np.arange(N) + .5) * h
X, Y = np.meshgrid(c, c, indexing='ij')
RHO = np.hypot(X, Y)
INF = 1e9


def integ(ztop, zbot, cond=True):
    return (np.clip(ztop - zbot, 0, None) * cond).sum() * h * h


def cap(r, rho2=None):
    rho2 = RHO * RHO if rho2 is None else rho2
    return np.sqrt(np.clip(r * r - rho2, 0, None))


def half(ang):
    """inward distance from the half plane through the axis at angle ang, the solid on its left"""
    return -math.sin(ang) * X + math.cos(ang) * Y


class Wedge:
    def __init__(self, A):
        self.A = A
        self.a = half(0.0)                 # side at angle 0: inside where > 0
        self.b = -half(A)                  # side at angle A: inside where > 0
        self.reflex = A > math.pi + 1e-9
        self.flat = abs(A - math.pi) < 1e-9
        self.whole = A > 2 * math.pi - 1e-9

    def inside(self, ga=0., gb=0.):
        """the wedge with its sides moved out by ga, gb (in by negative), sharp"""
        if self.whole:
            return np.ones_like(X, bool)
        if self.reflex:
            return (self.a > -ga) | (self.b > -gb)
        return (self.a > -ga) & (self.b > -gb)

    def dist(self):
        """distance from the wedge, 0 inside"""
        if self.whole:
            return np.zeros_like(X)
        da, db = np.clip(-self.a, 0, None), np.clip(-self.b, 0, None)
        if self.reflex:
            return np.minimum(da, db)
        # convex: nearest a side, or the axis where both feet fall off their sides
        ua = X                                            # along side a
        ub = math.cos(self.A) * X + math.sin(self.A) * Y  # along side b
        d = np.where((da > 0) & (ua >= 0) & ((db == 0) | (da >= db) | (ub < 0)), da, 0.)
        d = np.where((db > 0) & (ub >= 0) & (d == 0), db, d)
        d = np.where((da > 0) & (db > 0) & (ua >= 0) & (ub >= 0), np.maximum(da, db), d)
        d = np.where(((da > 0) | (db > 0)) & (ua < 0) & (ub < 0), RHO, d)
        d = np.where((da > 0) & (ua < 0) & (db == 0), RHO, d)
        d = np.where((db > 0) & (ub < 0) & (da == 0), RHO, d)
        return d

    def depth(self):
        """distance to the complement, 0 outside (for the erosion)"""
        if self.whole:
            return np.full_like(X, INF)
        if self.reflex:                                   # complement convex: quadrant between the sides
            ca, cb = np.clip(self.a, 0, None), np.clip(self.b, 0, None)
            both = (self.a > 0) & (self.b > 0)
            return np.where(both, np.hypot(ca, cb), np.maximum(ca, cb))
        return np.clip(np.minimum(self.a, self.b), 0, None)


def dil_ball_part(t, sdist, w):
    """z extent (top, bottom) of the points within t of ball & {z > 0 or all z} & one half space,
    at horizontal distance sdist >= 0 outside the half space, w the coordinate along its edge;
    returns (top, bottom_dome) where bottom_dome is the lower end when z > 0 is a face"""
    tau = np.sqrt(np.clip(t * t - sdist * sdist, 0, None))
    ok = sdist < t
    top = np.where(ok, cap(1, 1 - (R + tau) ** 2 + w * w), 0)
    aw = np.abs(w)
    bot = np.where(aw <= R, -tau, -np.sqrt(np.clip(tau * tau - (aw - R) ** 2, 0, None)))
    bot = np.where(ok & (aw < R + tau), bot, 0)
    top = np.where(aw < R + tau, top, 0)
    return top, bot


def values(A, dome, t):
    W = Wedge(A)
    V0 = integ(cap(R), 0. if dome else -cap(R), W.inside())
    zlo = (lambda r: 0 * X) if dome else (lambda r: -cap(r))
    res = {}
    # --- the sphere removed: the wall on the ball of R
    k = lambda g: W.inside(g, g)
    if W.flat and not W.whole:
        k = lambda g: Y > -g
    res['sphere', +1, 2] = integ(cap(R), np.maximum(-cap(R), -t) if dome else -cap(R), k(t)) - V0
    d = np.clip(-Y, 0, None) if W.flat else W.dist()
    low = -np.sqrt(np.clip(t * t - d * d, 0, None)) if dome else -cap(R)
    res['sphere', +1, 0] = integ(cap(R), np.maximum(-cap(R), low), d < t) - V0
    kin = (Y > t) if W.flat else W.inside(-t, -t)
    res['sphere', -1, 2] = V0 - integ(cap(R), t if dome else -cap(R), kin)
    dep = np.clip(Y, 0, None) if W.flat else W.depth()
    res['sphere', -1, 0] = V0 - integ(cap(R), t if dome else -cap(R), dep > t)
    # --- the bottom removed
    if dome:
        res['bottom', +1, 2] = integ(cap(R + t), 0, k(t)) - V0
        res['bottom', -1, 2] = V0 - integ(cap(R - t), 0, kin)
        res['bottom', -1, 0] = V0 - integ(cap(R - t), 0, dep > t)
        if W.whole:
            res['bottom', +1, 0] = res['bottom', +1, 2]
        else:
            # within t of ball & wedge: inside the wedge the ball of R + t, off it the
            # sides' slabs and the tubes round their arcs (the axis is no edge of a
            # reflex wedge, and of a convex one its tube lies off both feet)
            top = np.where(W.inside(), cap(R + t), 0)
            for s, u in ((W.a, X), (W.b, math.cos(A) * X + math.sin(A) * Y)):
                sd = np.clip(-s, 0, None)
                tp, _ = dil_ball_part(t, sd, u)
                # foot on the side's half disc (u >= 0), else on the axis
                tau = np.sqrt(np.clip(t * t - sd * sd, 0, None))
                onaxis = R + np.sqrt(np.clip(tau * tau - u * u, 0, None))
                tp = np.where(u >= 0, tp, np.where((sd < t) & (np.abs(u) < tau), onaxis, 0))
                top = np.maximum(top, np.where(s <= 0, tp, 0))
            if W.flat:
                top = np.where(Y > 0, cap(R + t), dil_ball_part(t, np.clip(-Y, 0, None), X)[0])
            res['bottom', +1, 0] = integ(top, 0) - V0
    # --- a side removed (the one at angle 0; the other is its mirror image)
    if not W.whole:
        zb = np.maximum(-cap(R + t), -t) if dome else None
        zi = t if dome else None
        if W.flat:
            reg = np.where(X < t, Y > -t, Y > 0)
            res['side', +1, 2] = integ(cap(R + t), zb if dome else -cap(R + t), reg) - V0
            reg = np.where(X < t, Y > t, Y > 0)
            res['side', -1, 2] = V0 - integ(cap(R - t), zi if dome else -cap(R - t), reg)
            # Arc: the tube round the edge between the two sides. Inward the
            # cavity keeps t from the kept side, a half plane ending at the axis.
            ray = np.where(X < 0, np.abs(Y), RHO)
            res['side', -1, 0] = V0 - integ(cap(R - t), zi if dome else -cap(R - t), (Y > 0) & (ray > t))
            # Outward: beside the removed face within t of the ball & bottom, across
            # it within t of the kept side itself, a quarter (half) disc.
            topH, botH = dil_ball_part(t, 0 * X, RHO)
            tau = np.sqrt(np.clip(t * t - Y * Y, 0, None))
            ax = np.abs(X)
            topQ = np.where(X <= 0, np.where(ax < R + tau, cap(1, 1 - (R + tau) ** 2 + X * X), 0),
                            R + np.sqrt(np.clip(tau * tau - X * X, 0, None)))
            botQ = np.where(X <= 0, np.where(ax <= R, -tau, -np.sqrt(np.clip(tau * tau - (ax - R) ** 2, 0, None))),
                            -np.sqrt(np.clip(tau * tau - X * X, 0, None)))
            okQ = (Y <= 0) & (np.abs(Y) < t) & np.where(X <= 0, ax < R + tau, ax < tau)
            if not dome:
                botH, botQ = -topH, -topQ
            top = np.where(Y > 0, topH, np.where(okQ, topQ, 0))
            bot = np.where(Y > 0, botH, np.where(okQ, botQ, 0))
            res['side', +1, 0] = integ(top, bot) - V0
        else:
            reg = W.inside(0, t)
            res['side', +1, 2] = integ(cap(R + t), zb if dome else -cap(R + t), reg) - V0
            reg = W.inside(0, -t)
            res['side', -1, 2] = V0 - integ(cap(R - t), zi if dome else -cap(R - t), reg)
            # inward the offsets cross and the Arc join is sharp too; a reflex wedge has
            # no kept edge at the axis once a side is gone
            res['side', -1, 0] = res['side', -1, 2]
            # outward, Arc: beside the removed face (a > 0 side of its plane, and beyond
            # it where the wedge is reflex) within t of the ball & bottom alone; across
            # it, within t of the ball & bottom & the kept side
            sd = np.clip(-W.b, 0, None)
            ub = math.cos(A) * X + math.sin(A) * Y
            topB, botB = dil_ball_part(t, sd, np.where(sd > 0, ub, RHO))
            topH, botH = dil_ball_part(t, 0 * X, RHO)
            if not dome:
                botB, botH = -topB, -topH
            if W.reflex:
                top = np.where(W.a > 0, topH, topB)
                bot = np.where(W.a > 0, botH, botB)
            else:
                top = np.where(W.a > 0, topB, 0)
                bot = np.where(W.a > 0, botB, 0)
            res['side', +1, 0] = integ(top, bot) - V0
    return V0, res


SHAPES = {"dome": (2 * math.pi, True, ("sphere", "bottom")),
          "halfdome": (math.pi, True, ("sphere", "bottom", "side", "side")),
          "quartdome": (math.pi / 2, True, ("sphere", "bottom", "side", "side")),
          "dome120": (2 * math.pi / 3, True, ("sphere", "bottom", "side", "side")),
          "dome270": (1.5 * math.pi, True, ("sphere", "bottom", "side", "side")),
          "lune90": (math.pi / 2, False, ("sphere", "side", "side")),
          "ball120": (2 * math.pi / 3, False, ("sphere", "side", "side")),
          "halfball": (math.pi, False, ("sphere", "side", "side"))}

if __name__ == "__main__":
    t = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    for name, (A, dome, faces) in SHAPES.items():
        V0, res = values(A, dome, t)
        for i, role in enumerate(faces):
            for d in (+1, -1):
                for j in (0, 2):
                    if (role, d, j) in res:
                        print("%s_f%d_%+d_j%d %.4f" % (name, i + 1, d, j, res[role, d, j]))
