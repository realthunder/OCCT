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
#include <BRepAdaptor_Curve2d.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepTools_ReShape.hxx>
#include <Geom_Curve.hxx>
#include <Geom_SphericalSurface.hxx>
#include <HLRAlgo_Projector.hxx>
#include <HLRBRep_Algo.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <Standard_NullObject.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>

//=================================================================================================
// Hidden line removal over a defective shape: a face on a curved surface holding an
// edge that carries no pcurve. HLR looks for the silhouette of such a face with
// Contap, which walks the face's edges as 2D arcs; BRepAdaptor_Curve2d leaves its
// adaptor empty for a pcurve-less edge, and Contap_ArcFunction::Value then evaluates
// it. That used to dereference a null curve inside Geom2dAdaptor_Curve::EvalD0.
//
// A planar face would hide the defect -- OCCT stores no pcurve for one and computes
// it on demand -- so every fixture here sits on a sphere.
//=================================================================================================

namespace
{
//! A spherical face whose wire holds a single edge built from a 3D curve alone.
TopoDS_Face faceWithPCurvelessEdge(TopoDS_Edge& theEdge)
{
  occ::handle<Geom_SphericalSurface> aSurface = new Geom_SphericalSurface(gp_Ax3(), 10.0);

  BRep_Builder aBuilder;
  TopoDS_Face  aFace;
  aBuilder.MakeFace(aFace, aSurface, Precision::Confusion());

  theEdge           = BRepBuilderAPI_MakeEdge(gp_Pnt(10.0, 0.0, 0.0), gp_Pnt(0.0, 10.0, 0.0));
  TopoDS_Wire aWire = BRepBuilderAPI_MakeWire(theEdge);
  aBuilder.Add(aFace, aWire);
  return aFace;
}

//! A whole sphere with every non-degenerate edge replaced by a 3D-only copy, so no
//! edge keeps a pcurve on the face that holds it.
TopoDS_Shape sphereWithoutPCurves()
{
  TopoDS_Shape aSphere = BRepPrimAPI_MakeSphere(10.0).Shape();

  BRepTools_ReShape aReShape;
  for (TopExp_Explorer anExp(aSphere, TopAbs_EDGE); anExp.More(); anExp.Next())
  {
    const TopoDS_Edge& anOldEdge = TopoDS::Edge(anExp.Current());
    if (BRep_Tool::Degenerated(anOldEdge))
    {
      continue;
    }
    double                  aFirst = 0.0, aLast = 0.0;
    occ::handle<Geom_Curve> aCurve = BRep_Tool::Curve(anOldEdge, aFirst, aLast);
    if (aCurve.IsNull())
    {
      continue;
    }
    TopoDS_Edge aNewEdge = BRepBuilderAPI_MakeEdge(aCurve, aFirst, aLast);
    aNewEdge.Orientation(anOldEdge.Orientation());
    aReShape.Replace(anOldEdge, aNewEdge);
  }
  return aReShape.Apply(aSphere);
}

//! Projects along Z. Returns false only if the algorithm let something other than a
//! Standard_Failure out; a crash fails the suite outright, which is the point.
bool runHiddenLineRemoval(const TopoDS_Shape& theShape)
{
  try
  {
    occ::handle<HLRBRep_Algo> anAlgo = new HLRBRep_Algo();
    anAlgo->Add(theShape, 0);
    anAlgo->Projector(HLRAlgo_Projector(gp_Ax2(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0))));
    anAlgo->Update();
    anAlgo->Hide();
  }
  catch (const Standard_Failure&)
  {
    // Refusing the defective face is a valid outcome; crashing on it is not.
  }
  catch (...)
  {
    return false;
  }
  return true;
}
} // namespace

//=================================================================================================

TEST(HLRBRep_PCurveLess_Test, EdgeWithoutPCurve_AdaptorStaysEmpty)
{
  TopoDS_Edge anEdge;
  TopoDS_Face aFace = faceWithPCurvelessEdge(anEdge);

  BRepAdaptor_Curve2d aCurve(anEdge, aFace);

  EXPECT_FALSE(aCurve.IsInitialized());
  EXPECT_THROW(aCurve.Value(0.0), Standard_NullObject);
}

//=================================================================================================

TEST(HLRBRep_PCurveLess_Test, HiddenLineRemovalOverBareFace)
{
  TopoDS_Edge anEdge;
  EXPECT_TRUE(runHiddenLineRemoval(faceWithPCurvelessEdge(anEdge)));
}

//=================================================================================================

TEST(HLRBRep_PCurveLess_Test, HiddenLineRemovalOverSphereWithoutPCurves)
{
  EXPECT_TRUE(runHiddenLineRemoval(sphereWithoutPCurves()));
}
