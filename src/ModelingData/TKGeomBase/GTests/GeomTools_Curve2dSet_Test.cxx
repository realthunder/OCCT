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

#include <gtest/gtest.h>

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepTools.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2d_Line.hxx>
#include <GeomTools_Curve2dSet.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <gp_Dir2d.hxx>
#include <gp_Pnt2d.hxx>

#include <sstream>

//=================================================================================================
// A pcurve table holds one entry per distinct written form, and the shape
// records still find their curve through it.
//
// The set is keyed by handle, so independently computed but identical pcurves
// used to be written once each. Merging them is only safe while Index() answers
// for a curve that was merged away -- that is what the records are written
// through, and without it they are written with index 0 and every pcurve in the
// file is lost on reading.
//=================================================================================================

namespace
{
//! Every face of the shape has a pcurve for every one of its edges.
bool allPCurvesPresent(const TopoDS_Shape& theShape)
{
  for (TopExp_Explorer aFaceExp(theShape, TopAbs_FACE); aFaceExp.More(); aFaceExp.Next())
  {
    const TopoDS_Face& aFace = TopoDS::Face(aFaceExp.Current());
    for (TopExp_Explorer anEdgeExp(aFace, TopAbs_EDGE); anEdgeExp.More(); anEdgeExp.Next())
    {
      double aFirst = 0.0, aLast = 0.0;
      const occ::handle<Geom2d_Curve> aPCurve =
        BRep_Tool::CurveOnSurface(TopoDS::Edge(anEdgeExp.Current()), aFace, aFirst, aLast);
      if (aPCurve.IsNull())
      {
        return false;
      }
    }
  }
  return true;
}

//! Write to a string and read back, the way a project file round trips.
TopoDS_Shape roundTrip(const TopoDS_Shape& theShape)
{
  std::ostringstream anOut;
  BRepTools::Write(theShape, anOut);

  std::istringstream anIn(anOut.str());
  BRep_Builder       aBuilder;
  TopoDS_Shape       aBack;
  BRepTools::Read(aBack, anIn, aBuilder);
  return aBack;
}
} // namespace

//=================================================================================================

TEST(GeomTools_Curve2dSet_Test, EqualCurvesShareOneEntry)
{
  GeomTools_Curve2dSet aSet;

  occ::handle<Geom2d_Curve> aFirst = new Geom2d_Line(gp_Pnt2d(0.0, 0.0), gp_Dir2d(1.0, 0.0));
  occ::handle<Geom2d_Curve> aSame  = new Geom2d_Line(gp_Pnt2d(0.0, 0.0), gp_Dir2d(1.0, 0.0));
  occ::handle<Geom2d_Curve> aOther = new Geom2d_Line(gp_Pnt2d(0.0, 3.0), gp_Dir2d(0.0, 1.0));

  const int aFirstIndex = aSet.Add(aFirst);
  const int aSameIndex  = aSet.Add(aSame);
  const int aOtherIndex = aSet.Add(aOther);

  EXPECT_EQ(aSameIndex, aFirstIndex) << "a curve written identically is not a second entry";
  EXPECT_NE(aOtherIndex, aFirstIndex);
  EXPECT_EQ(aSet.Extent(), 2);
}

//=================================================================================================

TEST(GeomTools_Curve2dSet_Test, MergedCurveStillHasAnIndex)
{
  GeomTools_Curve2dSet aSet;

  occ::handle<Geom2d_Curve> aFirst = new Geom2d_Line(gp_Pnt2d(0.0, 0.0), gp_Dir2d(1.0, 0.0));
  occ::handle<Geom2d_Curve> aSame  = new Geom2d_Line(gp_Pnt2d(0.0, 0.0), gp_Dir2d(1.0, 0.0));

  const int aIndex = aSet.Add(aFirst);
  aSet.Add(aSame);

  // The shape records are written through Index(). Answering 0 here is what
  // silently empties every pcurve of a written file.
  EXPECT_EQ(aSet.Index(aSame), aIndex);
  EXPECT_EQ(aSet.Index(aFirst), aIndex);
}

//=================================================================================================

TEST(GeomTools_Curve2dSet_Test, ClearForgetsTheMergedOnesToo)
{
  GeomTools_Curve2dSet aSet;

  occ::handle<Geom2d_Curve> aCurve = new Geom2d_Line(gp_Pnt2d(0.0, 0.0), gp_Dir2d(1.0, 0.0));
  aSet.Add(aCurve);
  aSet.Add(aCurve);
  aSet.Clear();

  EXPECT_EQ(aSet.Extent(), 0);
  EXPECT_EQ(aSet.Index(aCurve), 0);
}

//=================================================================================================

TEST(GeomTools_Curve2dSet_Test, BoxKeepsItsPCurvesThroughAWrite)
{
  // A box is the case that matters: its faces are planes and their edges give
  // pcurves that are written identically over and over.
  const TopoDS_Shape aBox = BRepPrimAPI_MakeBox(10.0, 20.0, 30.0).Shape();
  ASSERT_TRUE(allPCurvesPresent(aBox));

  const TopoDS_Shape aBack = roundTrip(aBox);
  ASSERT_FALSE(aBack.IsNull());
  EXPECT_TRUE(allPCurvesPresent(aBack)) << "pcurves did not survive the round trip";
}

//=================================================================================================

TEST(GeomTools_Curve2dSet_Test, RepeatedCurvedFacesKeepTheirPCurves)
{
  // The case that actually exercises the table. A pcurve lives in its surface's
  // own parameter space, so identical holes have identical pcurves however far
  // apart they are drilled -- several entries the set will merge into one.
  //
  // The faces have to be curved. A planar face carries no stored pcurve at all
  // and BRep_Tool::CurveOnSurface computes one on demand, which hides a table
  // that came back empty.
  TopoDS_Shape aShape = BRepPrimAPI_MakeBox(40.0, 40.0, 10.0).Shape();
  const gp_Pnt aCentres[4] = {gp_Pnt(10.0, 10.0, -5.0),
                              gp_Pnt(30.0, 10.0, -5.0),
                              gp_Pnt(10.0, 30.0, -5.0),
                              gp_Pnt(30.0, 30.0, -5.0)};
  for (const gp_Pnt& aCentre : aCentres)
  {
    const TopoDS_Shape aDrill =
      BRepPrimAPI_MakeCylinder(gp_Ax2(aCentre, gp_Dir(0.0, 0.0, 1.0)), 3.0, 20.0).Shape();
    aShape = BRepAlgoAPI_Cut(aShape, aDrill).Shape();
  }

  int aCurvedFaces = 0;
  for (TopExp_Explorer anExp(aShape, TopAbs_FACE); anExp.More(); anExp.Next())
  {
    if (BRepAdaptor_Surface(TopoDS::Face(anExp.Current())).GetType() != GeomAbs_Plane)
    {
      ++aCurvedFaces;
    }
  }
  ASSERT_GE(aCurvedFaces, 4) << "the fixture lost its holes";
  ASSERT_TRUE(allPCurvesPresent(aShape));

  const TopoDS_Shape aBack = roundTrip(aShape);
  ASSERT_FALSE(aBack.IsNull());
  EXPECT_TRUE(allPCurvesPresent(aBack)) << "pcurves did not survive the round trip";
}
