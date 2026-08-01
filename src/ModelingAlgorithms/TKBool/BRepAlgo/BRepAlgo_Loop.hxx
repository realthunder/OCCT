// Created on: 1995-11-10
// Created by: Yves FRICAUD
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

#ifndef _BRepAlgo_Loop_HeaderFile
#define _BRepAlgo_Loop_HeaderFile

#include <Standard.hxx>
#include <Standard_DefineAlloc.hxx>
#include <Standard_Handle.hxx>

#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <NCollection_List.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <NCollection_DataMap.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <BRepAlgo_Image.hxx>
class BRepAlgo_AsDes;
class TopoDS_Edge;

//! Builds the loops from a set of edges on a face.
class BRepAlgo_Loop
{
public:
  DEFINE_STANDARD_ALLOC

  Standard_EXPORT BRepAlgo_Loop();

  //! Init with <F> the set of edges must have
  //! pcurves on <F>.
  Standard_EXPORT void Init(const TopoDS_Face& F);

  //! Add E with <LV>. <E> will be copied and trim
  //! by vertices in <LV>.
  Standard_EXPORT void AddEdge(TopoDS_Edge& E, const NCollection_List<TopoDS_Shape>& LV);

  //! Add <E> as const edge, E can be in the result.
  Standard_EXPORT void AddConstEdge(const TopoDS_Edge& E);

  //! Add <LE> as a set of const edges.
  Standard_EXPORT void AddConstEdges(const NCollection_List<TopoDS_Shape>& LE);

  //! Sets the Image Vertex - Vertex
  Standard_EXPORT void SetImageVV(const BRepAlgo_Image& theImageVV);

  //! Make loops.
  Standard_EXPORT void Perform();

  //! Make loops, additionally checking for intersecting edges (used when
  //! making thick solid with concave removed faces).
  Standard_EXPORT void Perform(
    const NCollection_List<TopoDS_Shape>* ContextFaces,
    const occ::handle<BRepAlgo_AsDes>&    AsDes = occ::handle<BRepAlgo_AsDes>());

  //! Update VE map according to Image Vertex - Vertex
  Standard_EXPORT void UpdateVEmap(NCollection_IndexedDataMap<TopoDS_Shape,
                                                              NCollection_List<TopoDS_Shape>,
                                                              TopTools_ShapeMapHasher>& theVEmap);

  //! Cut the edge <E> in several edges <NE> on the
  //! vertices<VonE>.
  Standard_EXPORT void CutEdge(const TopoDS_Edge&                    E,
                               const NCollection_List<TopoDS_Shape>& VonE,
                               NCollection_List<TopoDS_Shape>&       NE) const;

  //! Returns the list of wires performed.
  //! can be an empty list.
  Standard_EXPORT const NCollection_List<TopoDS_Shape>& NewWires() const;

  //! Build faces from the wires result.
  Standard_EXPORT void WiresToFaces();

  //! Returns the list of faces.
  //! Warning: The method <WiresToFaces> as to be called before.
  //! can be an empty list.
  Standard_EXPORT const NCollection_List<TopoDS_Shape>& NewFaces() const;

  //! Returns the list of new edges built from an edge <E>
  //! it can be an empty list.
  Standard_EXPORT const NCollection_List<TopoDS_Shape>& NewEdges(const TopoDS_Edge& E) const;

  //! Returns the map edge => list of cut edges.
  //! Indexed so that iteration follows the order in which edges were cut -
  //! consumers must not depend on hash order.
  const NCollection_IndexedDataMap<TopoDS_Shape,
                                   NCollection_List<TopoDS_Shape>,
                                   TopTools_ShapeMapHasher>&
    CutEdges() const
  {
    return myCutEdges;
  }

  //! Returns the datamap of vertices with their substitutes.
  Standard_EXPORT void GetVerticesForSubstitute(
    NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>& VerVerMap) const;

  Standard_EXPORT void VerticesForSubstitute(
    NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>& VerVerMap);

  //! Set maximal tolerance used for comparing distances between vertices.
  void SetTolConf(const double theTolConf) { myTolConf = theTolConf; }

  //! Get maximal tolerance used for comparing distances between vertices.
  double GetTolConf() const { return myTolConf; }

private:
  void CutEdge(const TopoDS_Edge&                    E,
               const NCollection_List<TopoDS_Shape>& VonE,
               NCollection_List<TopoDS_Shape>&       NE,
               bool                                  KeepAll) const;

  void FindLoop();

  TopoDS_Face                    myFace;
  NCollection_List<TopoDS_Shape> myConstEdges;
  NCollection_List<TopoDS_Shape> myEdges;
  NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
                                 myVerOnEdges;
  NCollection_List<TopoDS_Shape> myNewWires;
  NCollection_List<TopoDS_Shape> myNewFaces;
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
                                                                           myCutEdges;
  NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher> myVerticesForSubstitute;
  BRepAlgo_Image                                                           myImageVV;
  double                                                                   myTolConf;
};

//! Helper class to enable collecting intersecting concave edges found when
//! building loop
class BRepAlgo_LoopIntersectingEdgeMap
{
public:
  Standard_EXPORT BRepAlgo_LoopIntersectingEdgeMap();
  Standard_EXPORT ~BRepAlgo_LoopIntersectingEdgeMap();

  Standard_EXPORT static NCollection_DataMap<TopoDS_Shape,
                                             NCollection_List<TopoDS_Shape>,
                                             TopTools_ShapeMapHasher>&
    EdgeMap();
};

#endif // _BRepAlgo_Loop_HeaderFile
