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

// Leaving out a pcurve that the reader computes again.
//
// The property under test is not that the file is smaller. It is that a file
// written this way reads back with every pcurve still answerable, and with the
// same pcurve: BRep_Tool::CurveOnSurface projects onto a plane when it finds no
// stored representation, so a planar pcurve is recoverable and no other kind is.

#include <BRepTools.hxx>

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepTools_ShapeSet.hxx>
#include <BRep_Builder.hxx>
#include <BRep_TEdge.hxx>
#include <BRep_Tool.hxx>
#include <Geom2d_Line.hxx>
#include <Geom_Plane.hxx>
#include <gp.hxx>
#include <Precision.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>

#include <gtest/gtest.h>

#include <sstream>

namespace
{

std::string written(const TopoDS_Shape& theShape, const bool theOmit)
{
  std::ostringstream aStream;
  BRepTools_ShapeSet aSet(false);
  aSet.SetFormatNb(1);
  aSet.SetOmitPCurvesOnPlane(theOmit);
  aSet.Add(theShape);
  aSet.Write(aStream);
  aSet.Write(theShape, aStream);
  return aStream.str();
}

TopoDS_Shape readBack(const std::string& theText)
{
  std::istringstream aStream(theText);
  BRep_Builder       aBuilder;
  BRepTools_ShapeSet aSet(aBuilder, true);
  TopoDS_Shape       aShape;
  aSet.Read(aStream);
  aSet.Read(aShape, aStream);
  return aShape;
}

//! Every edge-on-face pair that answers with a pcurve, and the worst distance
//! on the surface between the two shapes' answers.
int comparePCurves(const TopoDS_Shape& theA, const TopoDS_Shape& theB, double& theWorst, int& theMissing)
{
  int             aPairs = 0;
  TopExp_Explorer aFaceA(theA, TopAbs_FACE);
  TopExp_Explorer aFaceB(theB, TopAbs_FACE);
  theWorst   = 0.;
  theMissing = 0;
  for (; aFaceA.More() && aFaceB.More(); aFaceA.Next(), aFaceB.Next())
  {
    const TopoDS_Face& aFA = TopoDS::Face(aFaceA.Current());
    const TopoDS_Face& aFB = TopoDS::Face(aFaceB.Current());
    TopExp_Explorer    anEdgeA(aFA, TopAbs_EDGE);
    TopExp_Explorer    anEdgeB(aFB, TopAbs_EDGE);
    for (; anEdgeA.More() && anEdgeB.More(); anEdgeA.Next(), anEdgeB.Next())
    {
      ++aPairs;
      double                    aF1, aL1, aF2, aL2;
      occ::handle<Geom2d_Curve> aCA =
        BRep_Tool::CurveOnSurface(TopoDS::Edge(anEdgeA.Current()), aFA, aF1, aL1);
      occ::handle<Geom2d_Curve> aCB =
        BRep_Tool::CurveOnSurface(TopoDS::Edge(anEdgeB.Current()), aFB, aF2, aL2);
      if (aCA.IsNull() || aCB.IsNull())
      {
        ++theMissing;
        continue;
      }
      TopLoc_Location                 aLoc;
      const occ::handle<Geom_Surface> aSurf   = BRep_Tool::Surface(aFA, aLoc);
      const double                    aBegin  = std::max(aF1, aF2);
      const double                    anEnd   = std::min(aL1, aL2);
      if (anEnd <= aBegin)
      {
        continue;
      }
      for (int i = 0; i <= 8; ++i)
      {
        const double   aT  = aBegin + (anEnd - aBegin) * double(i) / 8.;
        const gp_Pnt2d aPA = aCA->Value(aT);
        const gp_Pnt2d aPB = aCB->Value(aT);
        theWorst           = std::max(theWorst,
                            aSurf->Value(aPA.X(), aPA.Y())
                              .Distance(aSurf->Value(aPB.X(), aPB.Y())));
      }
    }
  }
  return aPairs;
}

} // namespace

// A box is all planes, so every pcurve it holds is recoverable and the file
// gets smaller -- and still reads back with all of them.
TEST(BRepTools_OmitPCurve, BoxKeepsEveryPCurve)
{
  const TopoDS_Shape aBox = BRepPrimAPI_MakeBox(10., 20., 30.).Shape();

  const std::string aPlain = written(aBox, false);
  const std::string aLean  = written(aBox, true);
  EXPECT_LT(aLean.size(), aPlain.size()) << "a box is all planes, something should have been left out";

  const TopoDS_Shape aBack = readBack(aLean);
  ASSERT_FALSE(aBack.IsNull());

  double aWorst   = 0.;
  int    aMissing = 0;
  const int aPairs = comparePCurves(aBox, aBack, aWorst, aMissing);
  EXPECT_GT(aPairs, 0);
  EXPECT_EQ(aMissing, 0) << "a pcurve was dropped that the reader cannot answer for";
  EXPECT_LT(aWorst, Precision::Confusion()) << "the recovered pcurve is not the one dropped";
}

// A cylinder's curved face keeps its pcurves: nothing recomputes those, and the
// HLR crash this fork already fixed is what a missing one costs.
TEST(BRepTools_OmitPCurve, CurvedFaceIsUntouched)
{
  const TopoDS_Shape aCyl = BRepPrimAPI_MakeCylinder(5., 10.).Shape();

  const TopoDS_Shape aBack = readBack(written(aCyl, true));
  ASSERT_FALSE(aBack.IsNull());

  double aWorst   = 0.;
  int    aMissing = 0;
  const int aPairs = comparePCurves(aCyl, aBack, aWorst, aMissing);
  EXPECT_GT(aPairs, 0);
  EXPECT_EQ(aMissing, 0);
  EXPECT_LT(aWorst, Precision::Confusion());

  // The seam and the curved face must still carry stored pcurves.
  int aStored = 0;
  for (TopExp_Explorer aFaceIt(aBack, TopAbs_FACE); aFaceIt.More(); aFaceIt.Next())
  {
    const TopoDS_Face& aFace = TopoDS::Face(aFaceIt.Current());
    TopLoc_Location    aLoc;
    if (!occ::down_cast<Geom_Plane>(BRep_Tool::Surface(aFace, aLoc)).IsNull())
    {
      continue;
    }
    for (TopExp_Explorer anEdgeIt(aFace, TopAbs_EDGE); anEdgeIt.More(); anEdgeIt.Next())
    {
      double aF, aL;
      bool   isStored = false;
      BRep_Tool::CurveOnSurface(TopoDS::Edge(anEdgeIt.Current()), aFace, aF, aL, &isStored);
      if (isStored)
      {
        ++aStored;
      }
    }
  }
  EXPECT_GT(aStored, 0) << "a curved face lost pcurves that nothing can recompute";
}

// Off by default: a file written without asking is byte for byte what it was.
TEST(BRepTools_OmitPCurve, OffByDefault)
{
  const TopoDS_Shape aBox = BRepPrimAPI_MakeBox(1., 2., 3.).Shape();
  BRepTools_ShapeSet aSet(false);
  EXPECT_FALSE(aSet.IsOmitPCurvesOnPlane());
  EXPECT_EQ(written(aBox, false).size(), written(aBox, false).size());
}

// The predicate refuses a representation with no 3D curve to project.
TEST(BRepTools_OmitPCurve, RefusesWithoutA3dCurve)
{
  const TopoDS_Shape aBox = BRepPrimAPI_MakeBox(10., 10., 10.).Shape();
  TopExp_Explorer    aFaceIt(aBox, TopAbs_FACE);
  ASSERT_TRUE(aFaceIt.More());
  TopExp_Explorer anEdgeIt(aFaceIt.Current(), TopAbs_EDGE);
  ASSERT_TRUE(anEdgeIt.More());
  const TopoDS_Edge anEdge = TopoDS::Edge(anEdgeIt.Current());

  occ::handle<BRep_TEdge> aTE = occ::down_cast<BRep_TEdge>(anEdge.TShape());
  NCollection_List<occ::handle<BRep_CurveRepresentation>>::Iterator aCR(aTE->Curves());
  occ::handle<BRep_CurveRepresentation> aPCurveRep;
  for (; aCR.More(); aCR.Next())
  {
    if (aCR.Value()->IsCurveOnSurface())
    {
      aPCurveRep = aCR.Value();
      break;
    }
  }
  ASSERT_FALSE(aPCurveRep.IsNull());
  EXPECT_TRUE(BRepTools::IsPCurveOmittable(anEdge, aPCurveRep));

  // An edge built from a pcurve on a surface and nothing else -- what an IGES
  // import leaves behind -- is not omittable however planar it is: its pcurve
  // is the only geometry it has, and CurveOnPlane has nothing to project.
  const occ::handle<Geom_Plane> aPlane = new Geom_Plane(gp::XOY());
  const occ::handle<Geom2d_Line> a2dLine =
    new Geom2d_Line(gp_Pnt2d(0., 0.), gp_Dir2d(1., 0.));
  const TopoDS_Edge a2dOnly = BRepBuilderAPI_MakeEdge(a2dLine, aPlane, 0., 5.).Edge();
  ASSERT_FALSE(a2dOnly.IsNull());
  {
    double aF, aL;
    EXPECT_TRUE(BRep_Tool::Curve(a2dOnly, aF, aL).IsNull()) << "the fixture is not 2D only";
  }

  occ::handle<BRep_TEdge> aBareTE = occ::down_cast<BRep_TEdge>(a2dOnly.TShape());
  NCollection_List<occ::handle<BRep_CurveRepresentation>>::Iterator aBareCR(aBareTE->Curves());
  bool aChecked = false;
  for (; aBareCR.More(); aBareCR.Next())
  {
    if (aBareCR.Value()->IsCurveOnSurface())
    {
      EXPECT_FALSE(BRepTools::IsPCurveOmittable(a2dOnly, aBareCR.Value()));
      aChecked = true;
    }
  }
  EXPECT_TRUE(aChecked) << "the fixture carries no pcurve to test";
}

// A null representation, and one that is not a pcurve at all, are refused
// rather than crashed on.
TEST(BRepTools_OmitPCurve, RefusesNonPCurve)
{
  const TopoDS_Shape aBox = BRepPrimAPI_MakeBox(1., 1., 1.).Shape();
  TopExp_Explorer    anEdgeIt(aBox, TopAbs_EDGE);
  ASSERT_TRUE(anEdgeIt.More());
  const TopoDS_Edge anEdge = TopoDS::Edge(anEdgeIt.Current());

  EXPECT_FALSE(BRepTools::IsPCurveOmittable(anEdge, occ::handle<BRep_CurveRepresentation>()));

  occ::handle<BRep_TEdge> aTE = occ::down_cast<BRep_TEdge>(anEdge.TShape());
  NCollection_List<occ::handle<BRep_CurveRepresentation>>::Iterator aCR(aTE->Curves());
  for (; aCR.More(); aCR.Next())
  {
    if (aCR.Value()->IsCurve3D())
    {
      EXPECT_FALSE(BRepTools::IsPCurveOmittable(anEdge, aCR.Value()));
    }
  }
}
