// Copyright (c) 1999-2014 OPEN CASCADE SAS
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

// 25.12.98 pdn: adding empty constructor

#include <BRep_TEdge.hxx>
#include <BRep_TFace.hxx>
#include <BRep_Tool.hxx>
#include <BRep_TVertex.hxx>
#include <BRepLib.hxx>
#include <ShapeFix_ShapeTolerance.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>

#include <map>
#include <utility>

namespace
{
//! The tolerances a call sets on Immutable parts (a fork flag: a value its
//! owner shares, a feature's Shape in FreeCAD), by TShape: each is set on a
//! thawed copy of the part instead, put into the shape the call was given
//! in its place (BRepLib_ThawedCopy, BRepLib_PutThawed) -- as BRepLib's
//! UpdateTolerances does. The last one set wins, as it would in place.
using FrozenTolerances = std::map<const TopoDS_TShape*, std::pair<TopoDS_Shape, double>>;

double toleranceOf(const TopoDS_Shape& theS);

//! Sets the tolerance of the vertex, edge or face <theS>; kept in <theFrozen>
//! for a copy where <theS> is Immutable and the tolerance would change.
void setTolerance(const TopoDS_Shape& theS, const double theTol, FrozenTolerances& theFrozen)
{
  const TopoDS_TShape* aTS = theS.TShape().get();
  if (theS.Immutable())
  {
    const auto aKept = theFrozen.find(aTS);
    if (aKept != theFrozen.end())
    {
      aKept->second.second = theTol;
    }
    else if (theTol != toleranceOf(theS))
    {
      theFrozen[aTS] = {theS, theTol};
    }
    return;
  }
  switch (theS.ShapeType())
  {
    case TopAbs_VERTEX:
      static_cast<BRep_TVertex*>(theS.TShape().get())->Tolerance(theTol);
      break;
    case TopAbs_EDGE:
      static_cast<BRep_TEdge*>(theS.TShape().get())->Tolerance(theTol);
      break;
    case TopAbs_FACE:
      static_cast<BRep_TFace*>(theS.TShape().get())->Tolerance(theTol);
      break;
    default:
      break;
  }
}

//! The tolerance of the vertex, edge or face <theS>.
double toleranceOf(const TopoDS_Shape& theS)
{
  switch (theS.ShapeType())
  {
    case TopAbs_VERTEX:
      return BRep_Tool::Tolerance(TopoDS::Vertex(theS));
    case TopAbs_EDGE:
      return BRep_Tool::Tolerance(TopoDS::Edge(theS));
    default:
      return BRep_Tool::Tolerance(TopoDS::Face(theS));
  }
}

//! The copies of the Immutable parts <theFrozen> keeps, their tolerances set,
//! put into <theRoot> in their places. An Immutable <theRoot> is a refusal, as
//! a change made in place to it was.
void putFrozen(const TopoDS_Shape& theRoot, const FrozenTolerances& theFrozen)
{
  if (theFrozen.empty())
  {
    return;
  }
  NCollection_List<TopoDS_Shape> aCopies;
  for (const auto& aKept : theFrozen)
  {
    const TopoDS_Shape aCopy = BRepLib_ThawedCopy(aKept.second.first);
    switch (aCopy.ShapeType())
    {
      case TopAbs_VERTEX:
        static_cast<BRep_TVertex*>(aCopy.TShape().get())->Tolerance(aKept.second.second);
        break;
      case TopAbs_EDGE:
        static_cast<BRep_TEdge*>(aCopy.TShape().get())->Tolerance(aKept.second.second);
        break;
      default:
        static_cast<BRep_TFace*>(aCopy.TShape().get())->Tolerance(aKept.second.second);
        break;
    }
    aCopies.Append(aCopy);
  }
  BRepLib_PutThawed(theRoot, aCopies);
}

//! ShapeFix_ShapeTolerance::LimitTolerance, the changes to Immutable parts kept
//! in <theFrozen>.
bool limitTolerance(const TopoDS_Shape&    theShape,
                    const double           theMin,
                    const double           theMax,
                    const TopAbs_ShapeEnum theType,
                    FrozenTolerances&      theFrozen)
{
  const bool isMax = theMax >= theMin;
  bool       isDone = false;
  if (theType == TopAbs_VERTEX || theType == TopAbs_EDGE || theType == TopAbs_FACE)
  {
    for (TopExp_Explorer ex(theShape, theType); ex.More(); ex.Next())
    {
      const TopoDS_Shape& aS = ex.Current();
      // a frozen part already changed: as the copy has it
      const auto   aKept = theFrozen.find(aS.TShape().get());
      const double aPrec = aKept != theFrozen.end() ? aKept->second.second : toleranceOf(aS);
      if (isMax && aPrec > theMax)
      {
        setTolerance(aS, theMax, theFrozen);
        isDone = true;
      }
      else if (aPrec < theMin)
      {
        setTolerance(aS, theMin, theFrozen);
        isDone = true;
      }
    }
  }
  else if (theType == TopAbs_WIRE)
  {
    for (TopExp_Explorer ex(theShape, TopAbs_EDGE); ex.More(); ex.Next())
    {
      const TopoDS_Edge E = TopoDS::Edge(ex.Current());
      limitTolerance(E, theMin, theMax, TopAbs_EDGE, theFrozen);
      TopoDS_Vertex V1, V2;
      TopExp::Vertices(E, V1, V2);
      if (!V1.IsNull())
      {
        isDone |= limitTolerance(V1, theMin, theMax, TopAbs_VERTEX, theFrozen);
      }
      if (!V2.IsNull())
      {
        isDone |= limitTolerance(V2, theMin, theMax, TopAbs_VERTEX, theFrozen);
      }
    }
  }
  else
  {
    isDone |= limitTolerance(theShape, theMin, theMax, TopAbs_VERTEX, theFrozen);
    isDone |= limitTolerance(theShape, theMin, theMax, TopAbs_EDGE, theFrozen);
    isDone |= limitTolerance(theShape, theMin, theMax, TopAbs_FACE, theFrozen);
  }
  return isDone;
}

//! ShapeFix_ShapeTolerance::SetTolerance, the changes to Immutable parts kept
//! in <theFrozen>.
void setTolerances(const TopoDS_Shape&    theShape,
                   const double           thePrec,
                   const TopAbs_ShapeEnum theType,
                   FrozenTolerances&      theFrozen)
{
  //   VERTEX ou EDGE ou FACE : ces types seulement
  //   WIRE : EDGE + VERTEX
  //   Autres : TOUT (donc == WIRE + FACE)
  if (theType == TopAbs_VERTEX || theType == TopAbs_EDGE || theType == TopAbs_FACE)
  {
    for (TopExp_Explorer ex(theShape, theType); ex.More(); ex.Next())
    {
      setTolerance(ex.Current(), thePrec, theFrozen);
    }
  }
  else if (theType == TopAbs_WIRE)
  {
    for (TopExp_Explorer ex(theShape, TopAbs_EDGE); ex.More(); ex.Next())
    {
      const TopoDS_Edge E = TopoDS::Edge(ex.Current());
      setTolerance(E, thePrec, theFrozen);
      TopoDS_Vertex V1, V2;
      TopExp::Vertices(E, V1, V2);
      if (!V1.IsNull())
      {
        setTolerance(V1, thePrec, theFrozen);
      }
      if (!V2.IsNull())
      {
        setTolerance(V2, thePrec, theFrozen);
      }
    }
  }
  else
  {
    setTolerances(theShape, thePrec, TopAbs_VERTEX, theFrozen);
    setTolerances(theShape, thePrec, TopAbs_EDGE, theFrozen);
    setTolerances(theShape, thePrec, TopAbs_FACE, theFrozen);
  }
}
} // namespace

//=================================================================================================

ShapeFix_ShapeTolerance::ShapeFix_ShapeTolerance() = default;

//=================================================================================================

bool ShapeFix_ShapeTolerance::LimitTolerance(const TopoDS_Shape&    shape,
                                             const double           tmin,
                                             const double           tmax,
                                             const TopAbs_ShapeEnum styp) const
{
  if (shape.IsNull() || tmin < 0)
  {
    return false;
  }
  FrozenTolerances aFrozen;
  const bool       isDone = limitTolerance(shape, tmin, tmax, styp, aFrozen);
  putFrozen(shape, aFrozen);
  return isDone;
}

//=================================================================================================

void ShapeFix_ShapeTolerance::SetTolerance(const TopoDS_Shape&    shape,
                                           const double           preci,
                                           const TopAbs_ShapeEnum styp) const
{
  if (shape.IsNull() || preci <= 0)
  {
    return;
  }
  FrozenTolerances aFrozen;
  setTolerances(shape, preci, styp, aFrozen);
  putFrozen(shape, aFrozen);
}
