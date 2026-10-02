# local04: makes local04_fuzzy_shared_vertex.brep, see ../README.md.
# Run: FreeCADCmd local04_fuzzy_shared_vertex.py [theta_deg d]   (default 10 1e-5)
#
# Two prisms whose faces share ONE vertex V (one TShape) at the origin:
# A's top face FA lies on z=0 and has a reflex corner at V (a 60 deg notch
# toward +y); B's bottom face FB has the same kind of corner at V, its other
# corners on the plane PB: z = tan(theta) (y + d). V lies d sin(theta) off PB,
# and V and its two FB edges get a tolerance covering that. FA and FB then
# meet on the line y=-d, z=0, which passes V at distance d through the solid
# side of both corners and crosses neither of V's edges.
import math
import os
import sys

import Part
from FreeCAD import Vector

S3 = math.sqrt(3.0)
NOTCH_A = [(0, 0), (5, 5 * S3), (10, 5 * S3), (10, -10), (-10, -10), (-10, 5 * S3), (-5, 5 * S3)]
NOTCH_B = [(0, 0), (2, 2 * S3), (4, 2 * S3), (4, -4), (-4, -4), (-4, 2 * S3), (-2, 2 * S3)]


def build(theta, d):
    t = math.radians(theta)
    tolV = d * math.sin(t) + 2e-7
    ptsA = [Vector(x, y, 0) for x, y in NOTCH_A]
    A = Part.Face(Part.makePolygon(ptsA + [ptsA[0]])).extrude(Vector(0, 0, -5))
    V = [v for v in A.Vertexes if v.Point.Length < 1e-12][0]
    V.Tolerance = tolV
    verts = [V] + [Part.Vertex(Vector(x, y, math.tan(t) * (y + d))) for x, y in NOTCH_B[1:]]
    edges = [Part.Edge(verts[i], verts[(i + 1) % len(verts)]) for i in range(len(verts))]
    for e in (edges[0], edges[-1]):
        e.Tolerance = tolV
    normal = Vector(0, -math.sin(t), math.cos(t))
    FB = Part.Face(Part.Plane(Vector(0, -d, 0), normal), Part.Wire(edges))
    B = FB.extrude(normal * 2)
    assert any(v.isSame(V) for v in B.Vertexes), "V is not shared"
    return A, B


args = [a for a in sys.argv[1:] if not a.endswith(".py")]
theta, d = (float(args[0]), float(args[1])) if len(args) >= 2 else (10.0, 1e-5)
A, B = build(theta, d)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "local04_fuzzy_shared_vertex.brep")
Part.Compound([A, B]).exportBrep(out)
print("wrote", out, "valid", A.isValid(), B.isValid())
