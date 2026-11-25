// Created on: 1996-09-05
// Created by: Yves FRICAUD
// Copyright (c) 1996-1999 Matra Datavision
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

#ifndef _BRepOffset_MakeLoops_HeaderFile
#define _BRepOffset_MakeLoops_HeaderFile

#include <Message_ProgressRange.hxx>

#include <Standard.hxx>
#include <Standard_DefineAlloc.hxx>
#include <Standard_Handle.hxx>

#include <TopoDS_Shape.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <NCollection_DataMap.hxx>
#include <NCollection_List.hxx>
class BRepAlgo_AsDes;
class BRepAlgo_Image;
class BRepOffset_Analyse;

class BRepOffset_MakeLoops
{
public:
  DEFINE_STANDARD_ALLOC

  Standard_EXPORT BRepOffset_MakeLoops();

  Standard_EXPORT void Build(const NCollection_List<TopoDS_Shape>& LF,
                             const occ::handle<BRepAlgo_AsDes>&    AsDes,
                             BRepAlgo_Image&                       Image,
                             BRepAlgo_Image&                       theImageVV,
                             const Message_ProgressRange&          theRange);

  Standard_EXPORT void BuildOnContext(const NCollection_List<TopoDS_Shape>& LContext,
                                      const BRepOffset_Analyse&             Analyse,
                                      const occ::handle<BRepAlgo_AsDes>&    AsDes,
                                      BRepAlgo_Image&                       Image,
                                      const bool                            InSide,
                                      const Message_ProgressRange&          theRange);

  Standard_EXPORT void BuildFaces(const NCollection_List<TopoDS_Shape>& LF,
                                  const occ::handle<BRepAlgo_AsDes>&    AsDes,
                                  BRepAlgo_Image&                       Image,
                                  const Message_ProgressRange&          theRange);

private:
  NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher> myVerVerMap;
};

//! (realthunder) The thawed copies of the frozen vertices the loops widened
//! in the computation under way (the TKBRep Immutable flag: a vertex of a
//! shape given that a FreeCAD document froze), for BRepLib_PutThawed to put
//! in the place of those vertices in the result; with <theIsCleared>, none
//! after.
Standard_EXPORT NCollection_List<TopoDS_Shape> BRepOffset_MakeLoops_ThawedVertices(
  const bool theIsCleared);

//! (realthunder) The distance within which the loops of this thread take two
//! vertices for one (BRepAlgo_Loop::SetTolConf), <theTolConf> from now on; 0
//! leaves the loop's own, 0.001. Returns the one it replaces.
Standard_EXPORT double BRepOffset_MakeLoops_TolConf(const double theTolConf);

#endif // _BRepOffset_MakeLoops_HeaderFile
