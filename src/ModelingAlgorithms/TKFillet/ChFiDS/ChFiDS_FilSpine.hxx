// Created on: 1995-04-24
// Created by: Modelistation
// Copyright (c) 1995-1999 Matra Datavision
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

#ifndef _ChFiDS_FilSpine_HeaderFile
#define _ChFiDS_FilSpine_HeaderFile

#include <ChFiDS_Spine.hxx>
#include <ChFiDS_ElSpine.hxx>
#include <Law_Function.hxx>
#include <NCollection_List.hxx>
#include <gp_XY.hxx>
#include <NCollection_Sequence.hxx>
#include <NCollection_DataMap.hxx>
#include <TopoDS_Shape.hxx>
#include <TopTools_ShapeMapHasher.hxx>

class TopoDS_Edge;
class TopoDS_Vertex;
class gp_XY;
class Law_Function;
class Law_Composite;

//! Provides data specific to the fillets -
//! vector or rule of evolution (C2).
class ChFiDS_FilSpine : public ChFiDS_Spine
{

public:
  Standard_EXPORT ChFiDS_FilSpine();

  Standard_EXPORT ChFiDS_FilSpine(const double Tol);

  Standard_EXPORT void Reset(const bool AllData = false) override;

  //! initializes the constant vector on edge E.
  Standard_EXPORT void SetRadius(const double Radius, const TopoDS_Edge& E);

  //! resets the constant vector on edge E.
  Standard_EXPORT void UnSetRadius(const TopoDS_Edge& E);

  //! initializes the vector on Vertex V.
  Standard_EXPORT void SetRadius(const double Radius, const TopoDS_Vertex& V);

  //! resets the vector on Vertex V.
  Standard_EXPORT void UnSetRadius(const TopoDS_Vertex& V);

  //! initializes the vector on the point of parameter W.
  Standard_EXPORT void SetRadius(const gp_XY& UandR, const int IinC);

  //! initializes the constant vector on all spine.
  Standard_EXPORT void SetRadius(const double Radius);

  //! initializes the rule of evolution on all spine.
  Standard_EXPORT void SetRadius(const occ::handle<Law_Function>& C, const int IinC);

  //! returns true if the radius is constant
  //! all along the spine.
  Standard_EXPORT bool IsConstant() const;

  //! returns true if the radius is constant
  //! all along the edge E.
  Standard_EXPORT bool IsConstant(const int IE) const;

  //! returns the radius if the fillet is constant
  //! all along the spine.
  Standard_EXPORT double Radius() const;

  //! returns the radius if the fillet is constant
  //! all along the edge E.
  Standard_EXPORT double Radius(const int IE) const;

  //! returns the radius if the fillet is constant
  //! all along the edge E.
  Standard_EXPORT double Radius(const TopoDS_Edge& E) const;

  Standard_EXPORT void AppendElSpine(const occ::handle<ChFiDS_ElSpine>& Els) override;

  Standard_EXPORT occ::handle<Law_Composite> Law(const occ::handle<ChFiDS_ElSpine>& Els) const;

  //! returns the elementary law
  Standard_EXPORT occ::handle<Law_Function>& ChangeLaw(const TopoDS_Edge& E);

  //! returns the maximum radius if the fillet is non-constant
  Standard_EXPORT double MaxRadFromSeqAndLaws() const;

  //! Sets the setback of the contour's first end (<isFirst>) or last end:
  //! the stripe stops <theDist> from the end's vertex, measured along the
  //! spine, and the corner there is filled by one patch over the opening
  //! (ChFi3d_Builder::PerformMoreThreeCorner). 0 asks for the smallest
  //! setback the corner allows -- where the stripes meet today; less than 0
  //! removes the setback. Kept over Reset(), as the radii are.
  Standard_EXPORT void SetSetback(const bool isFirst, const double theDist);

  //! The setback of the first or last end, less than 0 where none is set.
  Standard_EXPORT double Setback(const bool isFirst) const;

  //! Sets the depth of the face <theFace> at the first end's corner
  //! (<isFirst>) or the last end's, for a setback corner: the patch's
  //! boundary on that face bows into it, away from the vertex, <theDepth>
  //! at its middle, measured from the straight line between its ends. 0 or
  //! less removes it; the boundary is then the fairest curve between its
  //! ends. Kept over Reset(), as the setbacks are.
  Standard_EXPORT void SetFaceDepth(const bool          isFirst,
                                    const TopoDS_Shape& theFace,
                                    const double        theDepth);

  //! The depth of <theFace> at the first or last end's corner, 0 where none
  //! is set.
  Standard_EXPORT double FaceDepth(const bool isFirst, const TopoDS_Shape& theFace) const;

  DEFINE_STANDARD_RTTIEXT(ChFiDS_FilSpine, ChFiDS_Spine)

private:
  Standard_EXPORT occ::handle<Law_Composite> ComputeLaw(const occ::handle<ChFiDS_ElSpine>& Els);

  Standard_EXPORT void AppendLaw(const occ::handle<ChFiDS_ElSpine>& Els);

  NCollection_Sequence<gp_XY>                 parandrad;
  NCollection_List<occ::handle<Law_Function>> laws;
  double                                      mySetbackFirst = -1.;
  double                                      mySetbackLast  = -1.;
  NCollection_DataMap<TopoDS_Shape, double, TopTools_ShapeMapHasher> myFaceDepthFirst;
  NCollection_DataMap<TopoDS_Shape, double, TopTools_ShapeMapHasher> myFaceDepthLast;
};

#endif // _ChFiDS_FilSpine_HeaderFile
