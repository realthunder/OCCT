// Copyright (c) 2026 OPEN CASCADE SAS
//
// This file is part of Open CASCADE Technology software library.
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.
//
// Alternatively, this file may be used under the terms of Open CASCADE
// commercial license or contractual agreement.

// BRepOffset_Tool::AnalyticTwin: a solid whose B-spline faces and edges are
// planes, cylinders, cones, lines and circles under another name, made again
// of those.

#include <BRepOffset_Tool.hxx>

#include <BRepAdaptor_Curve.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepLib.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <ElSLib.hxx>
#include <GProp_GProps.hxx>
#include <Geom2d_Curve.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <GeomConvert_CurveToAnaCurve.hxx>
#include <GeomConvert_SurfToAnaSurf.hxx>
#include <GeomLib_IsPlanarSurface.hxx>
#include <GeomLib.hxx>
#include <GeomProjLib.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_BezierSurface.hxx>
#include <Geom_Circle.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Line.hxx>
#include <Geom_ElementarySurface.hxx>
#include <Geom_Plane.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom2d_Line.hxx>
#include <gp_Lin.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <NCollection_IndexedMap.hxx>
#include <NCollection_List.hxx>
#include <NCollection_Map.hxx>
#include <Precision.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax3.hxx>
#include <gp_Cone.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Pln.hxx>

#include <algorithm>
#include <cmath>
#include <vector>

namespace
{
typedef NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher> ShapeToShape;

//! The twin holds what the solid holds: nothing further apart than this is
//! taken for the same.
constexpr double THE_TWIN_TOL = Precision::Confusion();

//! Whether the surface is a plane, a cylinder, a cone, a torus or a sphere.
bool isTwinSurface(const occ::handle<Geom_Surface>& theS)
{
  return theS->IsKind(STANDARD_TYPE(Geom_Plane))
         || theS->IsKind(STANDARD_TYPE(Geom_CylindricalSurface))
         || theS->IsKind(STANDARD_TYPE(Geom_ConicalSurface))
         || theS->IsKind(STANDARD_TYPE(Geom_ToroidalSurface))
         || theS->IsKind(STANDARD_TYPE(Geom_SphericalSurface));
}

//! The middle of the widest gap in a set of angles of [0, 2 pi), and the
//! gap; with theIsHeld, only a gap that has the origin in it or at its end
//! (an angle within 1e-6 of the origin is at either end of the turn), none
//! else: the angles are then held in one turn from the origin, and
//! theMiddle is the middle of what they span there.
bool widestGap(std::vector<double> theAngles,
               const bool          theIsHeld,
               double&             theMiddle,
               double&             theGap)
{
  if (theAngles.empty())
  {
    return false;
  }
  std::sort(theAngles.begin(), theAngles.end());
  theGap    = theAngles.front() + 2. * M_PI - theAngles.back();
  theMiddle = theAngles.back() + theGap / 2.;
  double aLow = theAngles.back(), aHigh = theAngles.front() + 2. * M_PI;
  for (size_t i = 1; i < theAngles.size(); ++i)
  {
    if (theAngles[i] - theAngles[i - 1] > theGap)
    {
      theGap    = theAngles[i] - theAngles[i - 1];
      theMiddle = (theAngles[i] + theAngles[i - 1]) / 2.;
      aLow      = theAngles[i - 1];
      aHigh     = theAngles[i];
    }
  }
  if (!theIsHeld)
  {
    return true;
  }
  const double anEps = 1.e-6;
  if (aHigh >= 2. * M_PI - anEps)
  {
    // the gap runs to the end of the turn: the angles span [front, aLow]
    theMiddle = (theAngles.front() + aLow) / 2.;
    return true;
  }
  if (aLow <= anEps)
  {
    // the angles below the gap are the origin, which is the end of the turn
    theMiddle = (aHigh + 2. * M_PI) / 2.;
    return true;
  }
  return false;
}

//! The plane a B-spline or Bezier surface's poles lie on, each within the
//! twin's tolerance of it; false where they do not. The test the offset has
//! (GeomLib_IsPlanarSurface) holds the root of the points' least moment of
//! inertia against the tolerance, and that root is 1e-7 by rounding alone
//! on a slope 27 long: a wedge's slope with a fillet beside it was no plane
//! by it, its poles 4e-15 off one.
bool planeOfPoles(const occ::handle<Geom_Surface>& theS, gp_Pln& thePln)
{
  std::vector<gp_Pnt> aPoles;
  if (occ::handle<Geom_BSplineSurface> aBS = occ::down_cast<Geom_BSplineSurface>(theS))
  {
    for (int i = 1; i <= aBS->NbUPoles(); ++i)
    {
      for (int j = 1; j <= aBS->NbVPoles(); ++j)
      {
        aPoles.push_back(aBS->Pole(i, j));
      }
    }
  }
  else if (occ::handle<Geom_BezierSurface> aBz = occ::down_cast<Geom_BezierSurface>(theS))
  {
    for (int i = 1; i <= aBz->NbUPoles(); ++i)
    {
      for (int j = 1; j <= aBz->NbVPoles(); ++j)
      {
        aPoles.push_back(aBz->Pole(i, j));
      }
    }
  }
  if (aPoles.size() < 3)
  {
    return false;
  }
  NCollection_Array1<gp_Pnt> anArr(1, static_cast<int>(aPoles.size()));
  for (size_t i = 0; i < aPoles.size(); ++i)
  {
    anArr(static_cast<int>(i) + 1) = aPoles[i];
  }
  gp_Pnt aBary;
  gp_Dir aDX, aDY;
  double aGx = 0., aGy = 0., aGz = 0.;
  GeomLib::Inertia(anArr, aBary, aDX, aDY, aGx, aGy, aGz);
  if (aGy <= Precision::Confusion())
  {
    return false;
  }
  thePln = gp_Pln(gp_Ax3(aBary, aDX ^ aDY, aDX));
  for (const gp_Pnt& aP : aPoles)
  {
    if (thePln.Distance(aP) > THE_TWIN_TOL)
    {
      return false;
    }
  }
  return true;
}

//! The normal of a surface at a point's parameters, as its derivatives have
//! it; false where they have none.
bool normalAt(const occ::handle<Geom_Surface>& theS, const double theU, const double theV,
              gp_Vec& theN)
{
  gp_Pnt aP;
  gp_Vec aDU, aDV;
  theS->D1(theU, theV, aP, aDU, aDV);
  theN = aDU ^ aDV;
  if (theN.Magnitude() < gp::Resolution())
  {
    return false;
  }
  theN.Normalize();
  return true;
}

//! The angle of a point round a cylinder's or a cone's axis, from its X
//! direction, in [0, 2 pi).
double angleRound(const gp_Ax3& theAx, const gp_Pnt& theP)
{
  const gp_Vec aV(theAx.Location(), theP);
  double       anA =
    std::atan2(aV.Dot(gp_Vec(theAx.YDirection())), aV.Dot(gp_Vec(theAx.XDirection())));
  if (anA < 0.)
  {
    anA += 2. * M_PI;
  }
  return anA;
}

//! The plane, cylinder or cone a face's B-spline or Bezier surface is, with
//! the surface's own normal, and for a cylinder or a cone its seam where the
//! face has its seam edge -- or, the face having none, away from the face.
//! Null where it is none of them, where the face has a degenerated edge (an
//! apex), or where it goes all the way round with no seam edge. A torus only
//! where the face lies within one turn round its tube from the outer
//! equator, theVMiddle the middle of what it spans there; and its seam edge,
//! if any, a meridian.
occ::handle<Geom_Surface> analyticSurface(const TopoDS_Face& theF, double& theVMiddle)
{
  theVMiddle = 0.;
  occ::handle<Geom_Surface> aNull;
  TopLoc_Location           aLoc;
  occ::handle<Geom_Surface> aS = BRep_Tool::Surface(theF, aLoc);
  if (aS.IsNull())
  {
    return aNull;
  }
  if (aS->IsKind(STANDARD_TYPE(Geom_RectangularTrimmedSurface)))
  {
    aS = occ::down_cast<Geom_RectangularTrimmedSurface>(aS)->BasisSurface();
  }
  if (!aS->IsKind(STANDARD_TYPE(Geom_BSplineSurface))
      && !aS->IsKind(STANDARD_TYPE(Geom_BezierSurface)))
  {
    return aNull;
  }
  TopoDS_Edge   aSeam;
  TopoDS_Vertex aPole;
  for (TopExp_Explorer anExp(theF, TopAbs_EDGE); anExp.More(); anExp.Next())
  {
    const TopoDS_Edge& anE = TopoDS::Edge(anExp.Current());
    if (BRep_Tool::Degenerated(anE))
    {
      aPole = TopExp::FirstVertex(anE);
      continue;
    }
    if (BRep_Tool::IsClosed(anE, theF))
    {
      aSeam = anE;
    }
  }
  aS = occ::down_cast<Geom_Surface>(aS->Transformed(aLoc.Transformation()));
  GeomConvert_SurfToAnaSurf aConv(aS);
  aConv.SetConvType(GeomConvert_MinGap);
  occ::handle<Geom_Surface> aNew = aConv.ConvertToAnalytical(THE_TWIN_TOL);
  if (aNew.IsNull())
  {
    // a plane the conversion passes over (a wedge's slope), by the test the
    // offset has for it
    GeomLib_IsPlanarSurface aPlanar(aS, THE_TWIN_TOL);
    gp_Pln                  aPln;
    if (aPlanar.IsPlanar())
    {
      aPln = aPlanar.Plan();
    }
    else if (!planeOfPoles(aS, aPln))
    {
      return aNull;
    }
    aNew = new Geom_Plane(aPln);
  }
  else if (aConv.Gap() > THE_TWIN_TOL || !isTwinSurface(aNew))
  {
    return aNull;
  }
  // a degenerated edge is a sphere's pole, and nothing else here
  const bool isSphere = aNew->IsKind(STANDARD_TYPE(Geom_SphericalSurface));
  if (!aPole.IsNull() && !isSphere)
  {
    return aNull;
  }

  // a point of the face, on both surfaces
  TopExp_Explorer anExpE(theF, TopAbs_EDGE);
  if (!anExpE.More())
  {
    return aNull;
  }
  double                    aF2, aL2;
  occ::handle<Geom2d_Curve> aPC =
    BRep_Tool::CurveOnSurface(TopoDS::Edge(anExpE.Current()), theF, aF2, aL2);
  if (aPC.IsNull())
  {
    return aNull;
  }
  const gp_Pnt2d aUV = aPC->Value((aF2 + aL2) / 2.);
  gp_Vec         aNOld;
  if (!normalAt(aS, aUV.X(), aUV.Y(), aNOld))
  {
    return aNull;
  }
  const gp_Pnt aP = aS->Value(aUV.X(), aUV.Y());

  if (occ::handle<Geom_Plane> aPlane = occ::down_cast<Geom_Plane>(aNew))
  {
    gp_Ax3 anAx = aPlane->Position();
    gp_Vec aNNew;
    normalAt(aPlane, 0., 0., aNNew);
    if (aNNew.Dot(aNOld) < 0.)
    {
      anAx.ZReverse();
    }
    return new Geom_Plane(anAx);
  }

  // A cylinder, a cone or a torus: where its seam is to be.
  const bool   isCylinder = aNew->IsKind(STANDARD_TYPE(Geom_CylindricalSurface));
  const bool   isTorus    = aNew->IsKind(STANDARD_TYPE(Geom_ToroidalSurface));
  gp_Ax3       anAx       = occ::down_cast<Geom_ElementarySurface>(aNew)->Position();
  if (isSphere && !aPole.IsNull())
  {
    // the axis through the pole the face has
    const gp_Vec aToPole(anAx.Location(), BRep_Tool::Pnt(aPole));
    if (aToPole.Magnitude() < gp::Resolution())
    {
      return aNull;
    }
    anAx = gp_Ax3(gp_Ax2(anAx.Location(), gp_Dir(aToPole)));
  }
  double       aSeamAngle = 0.;
  if (!aSeam.IsNull())
  {
    BRepAdaptor_Curve aC(aSeam);
    aSeamAngle = angleRound(anAx, aC.Value((aC.FirstParameter() + aC.LastParameter()) / 2.));
    for (int i = 0; i <= 4 && isTorus; ++i)
    {
      // a meridian, not a turn round the axis
      double aD = std::abs(
        angleRound(anAx,
                   aC.Value(aC.FirstParameter()
                            + (aC.LastParameter() - aC.FirstParameter()) * (i + 0.5) / 5.))
        - aSeamAngle);
      aD = std::min(aD, 2. * M_PI - aD);
      if (aD > 1.e-6)
      {
        return aNull;
      }
    }
  }
  else
  {
    // the middle of the widest gap between the angles of the face's outline
    std::vector<double> anAngles;
    for (TopExp_Explorer anExp(theF, TopAbs_EDGE); anExp.More(); anExp.Next())
    {
      if (BRep_Tool::Degenerated(TopoDS::Edge(anExp.Current())))
      {
        continue;
      }
      BRepAdaptor_Curve aC(TopoDS::Edge(anExp.Current()));
      for (int i = 0; i <= 16; ++i)
      {
        const gp_Pnt aPi =
          aC.Value(aC.FirstParameter() + (aC.LastParameter() - aC.FirstParameter()) * i / 16.);
        // not a point on the axis: a sphere's pole has no angle round it
        if (gp_Lin(anAx.Axis()).Distance(aPi) > 1.e-6)
        {
          anAngles.push_back(angleRound(anAx, aPi));
        }
      }
    }
    double aGap = 0.;
    widestGap(anAngles, false, aSeamAngle, aGap);
    if (aGap < 0.01)
    {
      return aNull;
    }
  }
  const gp_Dir aX(gp_Vec(anAx.XDirection()) * std::cos(aSeamAngle)
                  + gp_Vec(anAx.YDirection()) * std::sin(aSeamAngle));
  for (int aTry = 0; aTry < 2; ++aTry)
  {
    gp_Ax3 aSeated(anAx.Location(), anAx.Direction(), aX);
    if (aTry == 1)
    {
      aSeated.YReverse();
    }
    occ::handle<Geom_Surface> aSeatedS;
    if (isSphere)
    {
      aSeatedS =
        new Geom_SphericalSurface(aSeated,
                                  occ::down_cast<Geom_SphericalSurface>(aNew)->Radius());
    }
    else if (isTorus)
    {
      const occ::handle<Geom_ToroidalSurface> aTorus = occ::down_cast<Geom_ToroidalSurface>(aNew);
      aSeatedS = new Geom_ToroidalSurface(aSeated, aTorus->MajorRadius(), aTorus->MinorRadius());
    }
    else if (isCylinder)
    {
      aSeatedS = new Geom_CylindricalSurface(
        aSeated,
        occ::down_cast<Geom_CylindricalSurface>(aNew)->Radius());
    }
    else
    {
      const occ::handle<Geom_ConicalSurface> aCone = occ::down_cast<Geom_ConicalSurface>(aNew);
      aSeatedS = new Geom_ConicalSurface(aSeated, aCone->SemiAngle(), aCone->RefRadius());
    }
    GeomAdaptor_Surface anAd(aSeatedS);
    double              aU, aV;
    if (isSphere)
    {
      ElSLib::Parameters(anAd.Sphere(), aP, aU, aV);
    }
    else if (isTorus)
    {
      ElSLib::Parameters(anAd.Torus(), aP, aU, aV);
    }
    else if (isCylinder)
    {
      ElSLib::Parameters(anAd.Cylinder(), aP, aU, aV);
    }
    else
    {
      ElSLib::Parameters(anAd.Cone(), aP, aU, aV);
    }
    gp_Vec aNNew;
    if (aSeatedS->Value(aU, aV).Distance(aP) <= 10. * THE_TWIN_TOL
        && normalAt(aSeatedS, aU, aV, aNNew)
        && aNNew.Dot(aNOld) > 0.)
    {
      if (isTorus)
      {
        // within one turn round the tube
        std::vector<double> aVs;
        for (TopExp_Explorer anExp(theF, TopAbs_EDGE); anExp.More(); anExp.Next())
        {
          BRepAdaptor_Curve aC(TopoDS::Edge(anExp.Current()));
          for (int i = 0; i <= 16; ++i)
          {
            double aUi, aVi;
            ElSLib::Parameters(
              anAd.Torus(),
              aC.Value(aC.FirstParameter() + (aC.LastParameter() - aC.FirstParameter()) * i / 16.),
              aUi,
              aVi);
            aVs.push_back(aVi >= 2. * M_PI ? 0. : aVi);
          }
        }
        double aVGap = 0.;
        if (!widestGap(aVs, true, theVMiddle, aVGap) || aVGap < 0.01)
        {
          return aNull;
        }
      }
      return aSeatedS;
    }
  }
  return aNull;
}

//! The line or circle an edge's B-spline or Bezier curve is, running the same
//! way between the same ends, with its range. Null where it is neither.
occ::handle<Geom_Curve> analyticCurve(const TopoDS_Edge& theE, double& theF, double& theL)
{
  occ::handle<Geom_Curve> aNull;
  double                  aF, aL;
  TopLoc_Location         aLoc;
  occ::handle<Geom_Curve> aC = BRep_Tool::Curve(theE, aLoc, aF, aL);
  if (aC.IsNull())
  {
    return aNull;
  }
  if (aC->IsKind(STANDARD_TYPE(Geom_TrimmedCurve)))
  {
    aC = occ::down_cast<Geom_TrimmedCurve>(aC)->BasisCurve();
  }
  if (!aC->IsKind(STANDARD_TYPE(Geom_BSplineCurve)) && !aC->IsKind(STANDARD_TYPE(Geom_BezierCurve)))
  {
    return aNull;
  }
  aC = occ::down_cast<Geom_Curve>(aC->Transformed(aLoc.Transformation()));
  double                  aGap = 0.;
  occ::handle<Geom_Curve> aNew =
    GeomConvert_CurveToAnaCurve::ComputeLine(aC, THE_TWIN_TOL, aF, aL, theF, theL, aGap);
  if (aNew.IsNull() || aGap > THE_TWIN_TOL)
  {
    aNew = GeomConvert_CurveToAnaCurve::ComputeCircle(aC, THE_TWIN_TOL, aF, aL, theF, theL, aGap);
  }
  if (aNew.IsNull() || aGap > THE_TWIN_TOL
      || (!aNew->IsKind(STANDARD_TYPE(Geom_Line)) && !aNew->IsKind(STANDARD_TYPE(Geom_Circle))))
  {
    return aNull;
  }
  // the same ends, left the same way (a circle from a rational B-spline does
  // not run at the B-spline's pace between them)
  gp_Pnt aP1, aP2;
  gp_Vec aT1, aT2;
  aC->D1(aF, aP1, aT1);
  aNew->D1(theF, aP2, aT2);
  if (aP1.Distance(aP2) > 10. * THE_TWIN_TOL
      || aC->Value(aL).Distance(aNew->Value(theL)) > 10. * THE_TWIN_TOL || aT1.Dot(aT2) <= 0.
      || theL <= theF)
  {
    return aNull;
  }
  return aNew;
}

//! An edge's pcurve on a plane, a cylinder, a cone or a torus, by its 3D
//! curve; for a seam edge the two of them (theVMiddle: analyticSurface). False where there is none.
bool addPCurve(const TopoDS_Edge& theE,
               const TopoDS_Face& theF,
               const bool         theIsSeam,
               const double       theVMiddle)
{
  BRep_Builder                    aB;
  const occ::handle<Geom_Surface> aS = BRep_Tool::Surface(theF);
  if (aS->IsKind(STANDARD_TYPE(Geom_Plane)))
  {
    occ::handle<Geom2d_Curve> aPC;
    bool                      isToUpdate = false;
    BRepLib::BuildPCurveForEdgeOnPlane(theE, theF, aPC, isToUpdate);
    if (aPC.IsNull())
    {
      return false;
    }
    aB.UpdateEdge(theE, aPC, theF, BRep_Tool::Tolerance(theE));
    return true;
  }
  double                  aF, aL;
  occ::handle<Geom_Curve> aC = BRep_Tool::Curve(theE, aF, aL);
  if (aC.IsNull())
  {
    return false;
  }
  double                    aTol = BRep_Tool::Tolerance(theE);
  occ::handle<Geom2d_Curve> aPC  = GeomProjLib::Curve2d(aC, aF, aL, aS, aTol);
  if (aPC.IsNull())
  {
    return false;
  }
  // On a torus, in the turn round the tube the face lies in: its outer
  // equator is the start of that turn and its end.
  if (aS->IsKind(STANDARD_TYPE(Geom_ToroidalSurface)))
  {
    const double aVShift =
      std::round((theVMiddle - aPC->Value((aF + aL) / 2.).Y()) / (2. * M_PI)) * 2. * M_PI;
    if (aVShift != 0.)
    {
      aPC = occ::down_cast<Geom2d_Curve>(aPC->Translated(gp_Vec2d(0., aVShift)));
    }
  }
  // within the surface's own turn
  const gp_Pnt2d aMid = aPC->Value((aF + aL) / 2.);
  if (theIsSeam)
  {
    // The seam runs along u = 0 and u = 2 pi; the face lies on the left of
    // the one its edge has as it is.
    const double aShift = -std::round(aMid.X() / (2. * M_PI)) * 2. * M_PI;
    occ::handle<Geom2d_Curve> aAt0 =
      occ::down_cast<Geom2d_Curve>(aPC->Translated(gp_Vec2d(aShift, 0.)));
    occ::handle<Geom2d_Curve> aAt2Pi =
      occ::down_cast<Geom2d_Curve>(aPC->Translated(gp_Vec2d(aShift + 2. * M_PI, 0.)));
    const bool isUp = aPC->Value(aL).Y() > aPC->Value(aF).Y();
    if (isUp)
    {
      aB.UpdateEdge(theE, aAt2Pi, aAt0, theF, aTol);
    }
    else
    {
      aB.UpdateEdge(theE, aAt0, aAt2Pi, theF, aTol);
    }
    return true;
  }
  if (aMid.X() < 0. || aMid.X() > 2. * M_PI)
  {
    aPC = occ::down_cast<Geom2d_Curve>(
      aPC->Translated(gp_Vec2d(-std::floor(aMid.X() / (2. * M_PI)) * 2. * M_PI, 0.)));
  }
  aB.UpdateEdge(theE, aPC, theF, aTol);
  return true;
}
} // namespace

//=================================================================================================

bool BRepOffset_Tool::AnalyticTwin(
  const TopoDS_Shape&                                                       theSolid,
  TopoDS_Shape&                                                             theTwin,
  NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>& theOrigin)
{
  theTwin.Nullify();
  theOrigin.Clear();
  if (theSolid.IsNull() || theSolid.ShapeType() != TopAbs_SOLID)
  {
    return false;
  }
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aFaces, aVertices;
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    anEdgeFaces;
  TopExp::MapShapes(theSolid, TopAbs_FACE, aFaces);
  TopExp::MapShapes(theSolid, TopAbs_VERTEX, aVertices);
  TopExp::MapShapesAndAncestors(theSolid, TopAbs_EDGE, TopAbs_FACE, anEdgeFaces);

  // The faces: each on its surface as a plane, a cylinder or a cone, or as
  // it is.
  BRep_Builder aB;
  ShapeToShape aNew;
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aConverted;
  bool         isAnyConverted = false;
  NCollection_DataMap<TopoDS_Shape, double, TopTools_ShapeMapHasher> aVMiddles;
  for (int i = 1; i <= aFaces.Extent(); ++i)
  {
    const TopoDS_Face&              aF = TopoDS::Face(aFaces(i));
    double                          aVMiddle = 0.;
    const occ::handle<Geom_Surface> aS       = analyticSurface(aF, aVMiddle);
    TopoDS_Face                     aNF;
    if (aS.IsNull())
    {
      aNF = TopoDS::Face(aF.EmptyCopied().Oriented(TopAbs_FORWARD));
    }
    else
    {
      aB.MakeFace(aNF, aS, BRep_Tool::Tolerance(aF));
      aConverted.Add(aF);
      aVMiddles.Bind(aF, aVMiddle);
      isAnyConverted = true;
    }
    aNew.Bind(aF, aNF);
  }
  for (int i = 1; i <= aVertices.Extent(); ++i)
  {
    const TopoDS_Vertex& aV = TopoDS::Vertex(aVertices(i));
    TopoDS_Vertex        aNV;
    aB.MakeVertex(aNV, BRep_Tool::Pnt(aV), BRep_Tool::Tolerance(aV));
    aNew.Bind(aV, aNV);
  }
  // No face to convert, no twin: lines and circles alone, between faces
  // that are what they were, change nothing an offset is refused for.
  if (!isAnyConverted)
  {
    return false;
  }

  // The edges: each on its curve as a line or a circle where every face it
  // bounds is a plane, a cylinder or a cone in the twin -- on another
  // surface it keeps its pcurve, which goes with the curve it has.
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aRecurved, aTouched;
  for (int i = 1; i <= anEdgeFaces.Extent(); ++i)
  {
    const TopoDS_Edge& anE = TopoDS::Edge(anEdgeFaces.FindKey(i));
    for (TopoDS_Iterator anIt(anE.Oriented(TopAbs_FORWARD)); anIt.More(); anIt.Next())
    {
      if (anIt.Value().Orientation() != TopAbs_FORWARD
          && anIt.Value().Orientation() != TopAbs_REVERSED)
      {
        return false;
      }
    }
    bool isOnAnalytic = !BRep_Tool::Degenerated(anE);
    for (NCollection_List<TopoDS_Shape>::Iterator anIt(anEdgeFaces(i)); anIt.More() && isOnAnalytic;
         anIt.Next())
    {
      isOnAnalytic = isTwinSurface(BRep_Tool::Surface(TopoDS::Face(aNew(anIt.Value()))));
    }
    double                  aF = 0., aL = 0.;
    occ::handle<Geom_Curve> aC;
    if (isOnAnalytic)
    {
      aC = analyticCurve(anE, aF, aL);
    }
    TopoDS_Edge aNE;
    if (aC.IsNull())
    {
      aNE = TopoDS::Edge(anE.EmptyCopied().Oriented(TopAbs_FORWARD));
    }
    else
    {
      aB.MakeEdge(aNE, aC, BRep_Tool::Tolerance(anE));
      aB.Range(aNE, aF, aL);
      aRecurved.Add(anE);
      isAnyConverted = true;
    }
    for (TopoDS_Iterator anIt(anE.Oriented(TopAbs_FORWARD)); anIt.More(); anIt.Next())
    {
      aB.Add(aNE, aNew(anIt.Value()).Oriented(anIt.Value().Orientation()));
    }
    aNew.Bind(anE, aNE);
  }
  if (!isAnyConverted)
  {
    return false;
  }

  // The faces' wires, and the pcurves that are new: of an edge on another
  // curve, and of any edge on another surface.
  try
  {
    for (int i = 1; i <= aFaces.Extent(); ++i)
    {
      const TopoDS_Face& aF  = TopoDS::Face(aFaces(i));
      TopoDS_Face        aNF = TopoDS::Face(aNew(aF));
      const bool         isNewSurface = aConverted.Contains(aF);
      for (TopoDS_Iterator anItW(aF.Oriented(TopAbs_FORWARD)); anItW.More(); anItW.Next())
      {
        if (anItW.Value().ShapeType() != TopAbs_WIRE)
        {
          return false;
        }
        TopoDS_Wire aNW;
        aB.MakeWire(aNW);
        NCollection_List<TopoDS_Shape> aPoles;
        for (TopoDS_Iterator anItE(anItW.Value().Oriented(TopAbs_FORWARD)); anItE.More();
             anItE.Next())
        {
          const TopoDS_Edge& anE = TopoDS::Edge(anItE.Value());
          const TopoDS_Edge  aNE = TopoDS::Edge(aNew(anE));
          double             aPF, aPL;
          if (isNewSurface && BRep_Tool::Degenerated(anE))
          {
            aPoles.Append(anE);
            aB.Add(aNW, aNE.Oriented(anE.Orientation()));
            continue;
          }
          if ((isNewSurface || aRecurved.Contains(anE))
              && BRep_Tool::CurveOnSurface(aNE, aNF, aPF, aPL).IsNull())
          {
            const double* aVMiddle = aVMiddles.Seek(aF);
            if (!addPCurve(aNE, aNF, BRep_Tool::IsClosed(anE, aF), aVMiddle ? *aVMiddle : 0.))
            {
              return false;
            }
            aTouched.Add(aNE);
          }
          aB.Add(aNW, aNE.Oriented(anE.Orientation()));
        }
        // A sphere's pole: its edge runs along the pole's line of the sphere,
        // between the parameters the meridians reach it at, the face on its
        // left.
        for (NCollection_List<TopoDS_Shape>::Iterator anItP(aPoles); anItP.More(); anItP.Next())
        {
          const TopoDS_Edge&   anE    = TopoDS::Edge(anItP.Value());
          const TopoDS_Edge    aNE    = TopoDS::Edge(aNew(anE));
          const TopoDS_Vertex  aPoleV = TopoDS::Vertex(aNew(TopExp::FirstVertex(anE)));
          const gp_Pnt         aPoleP = BRep_Tool::Pnt(aPoleV);
          const occ::handle<Geom_SphericalSurface> aSph =
            occ::down_cast<Geom_SphericalSurface>(BRep_Tool::Surface(aNF));
          if (aSph.IsNull())
          {
            return false;
          }
          const bool isNorth =
            gp_Vec(aSph->Location(), aPoleP).Dot(gp_Vec(aSph->Position().Direction())) > 0.;
          double aU1 = Precision::Infinite(), aU2 = -Precision::Infinite();
          for (TopoDS_Iterator anItE(aNW); anItE.More(); anItE.Next())
          {
            const TopoDS_Edge& aME = TopoDS::Edge(anItE.Value());
            if (BRep_Tool::Degenerated(aME))
            {
              continue;
            }
            for (int aSide = 0; aSide < 2; ++aSide)
            {
              double                    aMF, aML;
              occ::handle<Geom2d_Curve> aMPC = BRep_Tool::CurveOnSurface(
                TopoDS::Edge(aME.Oriented(aSide == 0 ? TopAbs_FORWARD : TopAbs_REVERSED)),
                aNF,
                aMF,
                aML);
              if (aMPC.IsNull())
              {
                continue;
              }
              for (int anEnd = 0; anEnd < 2; ++anEnd)
              {
                const gp_Pnt2d aUV = aMPC->Value(anEnd == 0 ? aMF : aML);
                if (aSph->Value(aUV.X(), aUV.Y()).Distance(aPoleP) <= 10. * THE_TWIN_TOL)
                {
                  aU1 = std::min(aU1, aUV.X());
                  aU2 = std::max(aU2, aUV.X());
                }
              }
            }
          }
          if (aU2 - aU1 < Precision::PConfusion())
          {
            return false;
          }
          // forward in the wire it runs back along u at the north pole and
          // on along u at the south
          const bool isOnward = isNorth == (anE.Orientation() == TopAbs_REVERSED);
          const double aVPole = isNorth ? M_PI / 2. : -M_PI / 2.;
          occ::handle<Geom2d_Line> aLine =
            new Geom2d_Line(gp_Pnt2d(0., aVPole), gp_Dir2d(isOnward ? 1. : -1., 0.));
          aB.UpdateEdge(aNE, aLine, aNF, BRep_Tool::Tolerance(anE));
          aB.Range(aNE, isOnward ? aU1 : -aU2, isOnward ? aU2 : -aU1);
        }
        aNW.Closed(anItW.Value().Closed());
        aB.Add(aNF, aNW.Oriented(anItW.Value().Orientation()));
      }
      aNew.ChangeFind(aF) = aNF;
    }
  }
  catch (Standard_Failure const&)
  {
    return false;
  }

  // The solid, shell by shell.
  TopoDS_Solid aSolid;
  aB.MakeSolid(aSolid);
  for (TopoDS_Iterator anItS(theSolid.Oriented(TopAbs_FORWARD)); anItS.More(); anItS.Next())
  {
    if (anItS.Value().ShapeType() != TopAbs_SHELL)
    {
      return false;
    }
    TopoDS_Shell aShell;
    aB.MakeShell(aShell);
    for (TopoDS_Iterator anItF(anItS.Value().Oriented(TopAbs_FORWARD)); anItF.More(); anItF.Next())
    {
      if (anItF.Value().ShapeType() != TopAbs_FACE)
      {
        return false;
      }
      aB.Add(aShell, aNew(anItF.Value()).Oriented(anItF.Value().Orientation()));
    }
    aShell.Closed(anItS.Value().Closed());
    aB.Add(aSolid, aShell.Oriented(anItS.Value().Orientation()));
  }
  aSolid.Orientation(theSolid.Orientation());
  for (NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>::Iterator anIt(aTouched); anIt.More();
       anIt.Next())
  {
    aB.SameParameter(TopoDS::Edge(anIt.Value()), false);
  }
  try
  {
    BRepLib::SameParameter(aSolid, THE_TWIN_TOL);
  }
  catch (Standard_Failure const&)
  {
    return false;
  }

  // The same solid: valid, of the same volume and the same area.
  if (!BRepCheck_Analyzer(aSolid).IsValid())
  {
    return false;
  }
  GProp_GProps aV1, aV2, aA1, aA2;
  BRepGProp::VolumeProperties(theSolid, aV1);
  BRepGProp::VolumeProperties(aSolid, aV2);
  BRepGProp::SurfaceProperties(theSolid, aA1);
  BRepGProp::SurfaceProperties(aSolid, aA2);
  if (std::abs(aV1.Mass() - aV2.Mass()) > 1.e-6 * std::max(std::abs(aV1.Mass()), 1.)
      || std::abs(aA1.Mass() - aA2.Mass()) > 1.e-6 * std::max(aA1.Mass(), 1.))
  {
    return false;
  }
  theTwin = aSolid;
  for (ShapeToShape::Iterator anIt(aNew); anIt.More(); anIt.Next())
  {
    theOrigin.Bind(anIt.Value(), anIt.Key());
  }
  return true;
}
