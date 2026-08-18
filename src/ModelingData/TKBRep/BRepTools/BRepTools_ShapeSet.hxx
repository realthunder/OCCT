// Created on: 1993-07-19
// Created by: Remi LEQUETTE
// Copyright (c) 1993-1999 Matra Datavision
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

#ifndef _BRepTools_ShapeSet_HeaderFile
#define _BRepTools_ShapeSet_HeaderFile

#include <Standard.hxx>
#include <Standard_DefineAlloc.hxx>
#include <Standard_Handle.hxx>

#include <BRep_Builder.hxx>
#include <GeomTools_SurfaceSet.hxx>
#include <GeomTools_CurveSet.hxx>
#include <GeomTools_Curve2dSet.hxx>
#include <Standard_Transient.hxx>
#include <NCollection_IndexedMap.hxx>
#include <TopTools_ShapeSet.hxx>
#include <Standard_OStream.hxx>
#include <Standard_IStream.hxx>
#include <TopAbs_ShapeEnum.hxx>

class TopoDS_Shape;

//! Contains a Shape and all its subshapes, locations
//! and geometries.
//!
//! The topology is inherited from TopTools.
class BRepTools_ShapeSet : public TopTools_ShapeSet
{
public:
  DEFINE_STANDARD_ALLOC

  //! Builds an empty ShapeSet.
  //! @param theWithTriangles flag to write triangulation data
  Standard_EXPORT BRepTools_ShapeSet(const bool theWithTriangles = true,
                                     const bool theWithNormals   = false);

  //! Builds an empty ShapeSet.
  //! @param theWithTriangles flag to write triangulation data
  Standard_EXPORT BRepTools_ShapeSet(const BRep_Builder& theBuilder,
                                     const bool          theWithTriangles = true,
                                     const bool          theWithNormals   = false);

  Standard_EXPORT ~BRepTools_ShapeSet() override;

  //! Return true if shape should be stored with triangles.
  bool IsWithTriangles() const { return myWithTriangles; }

  //! Return true if shape should be stored triangulation with normals.
  bool IsWithNormals() const { return myWithNormals; }

  //! Define if shape will be stored with triangles.
  //! Ignored (always written) if face defines only triangulation (no surface).
  void SetWithTriangles(const bool theWithTriangles) { myWithTriangles = theWithTriangles; }

  //! Define if shape will be stored triangulation with normals.
  //! Ignored (always written) if face defines only triangulation (no surface).
  void SetWithNormals(const bool theWithNormals) { myWithNormals = theWithNormals; }

  //! Return true if pcurves a plane can rebuild are left out of the file.
  bool IsOmitPCurvesOnPlane() const { return myOmitPCurvesOnPlane; }

  //! Define whether to leave out the pcurves that reading the file back
  //! computes again anyway, which is those on a plane -- see
  //! BRepTools::IsPCurveOmittable for the exact condition and how it is
  //! checked. Off by default, so a file written without asking is unchanged.
  //!
  //! Nothing is needed on the reading side: a representation that is not in
  //! the file is one BRep_Tool::CurveOnSurface answers for by projecting the
  //! edge's 3D curve onto the plane.
  void SetOmitPCurvesOnPlane(const bool theOmit) { myOmitPCurvesOnPlane = theOmit; }

  //! Clears the content of the set.
  Standard_EXPORT void Clear() override;

  //! Stores the geometry of <S>.
  Standard_EXPORT void AddGeometry(const TopoDS_Shape& S) override;

  //! Dumps the geometry of me on the stream <OS>.
  Standard_EXPORT void DumpGeometry(Standard_OStream& OS) const override;

  //! Writes the geometry of me on the stream <OS> in a
  //! format that can be read back by Read.
  Standard_EXPORT void WriteGeometry(
    Standard_OStream&            OS,
    const Message_ProgressRange& theProgress = Message_ProgressRange()) override;

  //! Reads the geometry of me from the stream <IS>.
  Standard_EXPORT void ReadGeometry(
    Standard_IStream&            IS,
    const Message_ProgressRange& theProgress = Message_ProgressRange()) override;

  //! Dumps the geometry of <S> on the stream <OS>.
  Standard_EXPORT void DumpGeometry(const TopoDS_Shape& S, Standard_OStream& OS) const override;

  //! Writes the geometry of <S> on the stream <OS> in a
  //! format that can be read back by Read.
  Standard_EXPORT void WriteGeometry(const TopoDS_Shape& S, Standard_OStream& OS) const override;

  //! Reads the geometry of a shape of type <T> from the
  //! stream <IS> and returns it in <S>.
  Standard_EXPORT void ReadGeometry(const TopAbs_ShapeEnum T,
                                    Standard_IStream&      IS,
                                    TopoDS_Shape&          S) override;

  //! Inserts the shape <S2> in the shape <S1>. This
  //! method must be redefined to use the correct
  //! builder.
  Standard_EXPORT void AddShapes(TopoDS_Shape& S1, const TopoDS_Shape& S2) override;

  Standard_EXPORT void Check(const TopAbs_ShapeEnum T, TopoDS_Shape& S) override;

  //! Reads the 3d polygons of me
  //! from the stream <IS>.
  Standard_EXPORT void ReadPolygon3D(
    Standard_IStream&            IS,
    const Message_ProgressRange& theProgress = Message_ProgressRange());

  //! Writes the 3d polygons
  //! on the stream <OS> in a format that can
  //! be read back by Read.
  Standard_EXPORT void WritePolygon3D(
    Standard_OStream&            OS,
    const bool                   Compact     = true,
    const Message_ProgressRange& theProgress = Message_ProgressRange()) const;

  //! Dumps the 3d polygons
  //! on the stream <OS>.
  Standard_EXPORT void DumpPolygon3D(Standard_OStream& OS) const;

  //! Reads the triangulation of me
  //! from the stream <IS>.
  Standard_EXPORT void ReadTriangulation(
    Standard_IStream&            IS,
    const Message_ProgressRange& theProgress = Message_ProgressRange());

  //! Writes the triangulation
  //! on the stream <OS> in a format that can
  //! be read back by Read.
  Standard_EXPORT void WriteTriangulation(
    Standard_OStream&            OS,
    const bool                   Compact     = true,
    const Message_ProgressRange& theProgress = Message_ProgressRange()) const;

  //! Dumps the triangulation
  //! on the stream <OS>.
  Standard_EXPORT void DumpTriangulation(Standard_OStream& OS) const;

  //! Reads the polygons on triangulation of me
  //! from the stream <IS>.
  Standard_EXPORT void ReadPolygonOnTriangulation(
    Standard_IStream&            IS,
    const Message_ProgressRange& theProgress = Message_ProgressRange());

  //! Writes the polygons on triangulation
  //! on the stream <OS> in a format that can
  //! be read back by Read.
  Standard_EXPORT void WritePolygonOnTriangulation(
    Standard_OStream&            OS,
    const bool                   Compact     = true,
    const Message_ProgressRange& theProgress = Message_ProgressRange()) const;

  //! Dumps the polygons on triangulation
  //! on the stream <OS>.
  Standard_EXPORT void DumpPolygonOnTriangulation(Standard_OStream& OS) const;

  //! The geometry tables the shape records index into.
  //!
  //! A subclass that writes or reads the tables itself needs them: the indices
  //! WriteGeometry() emits per shape are positions in these sets, so a format
  //! that stores a table entry differently -- naming geometry held elsewhere
  //! rather than writing it out -- has to see and seed them. Reading them is
  //! also the only way to ask which curve or surface an index stands for.
  //! @{
  const GeomTools_SurfaceSet& Surfaces() const { return mySurfaces; }

  GeomTools_SurfaceSet& ChangeSurfaces() { return mySurfaces; }

  const GeomTools_CurveSet& Curves() const { return myCurves; }

  GeomTools_CurveSet& ChangeCurves() { return myCurves; }

  const GeomTools_Curve2dSet& Curves2d() const { return myCurves2d; }

  GeomTools_Curve2dSet& ChangeCurves2d() { return myCurves2d; }

  //! @}

private:
  BRep_Builder                                            myBuilder;
  GeomTools_SurfaceSet                                    mySurfaces;
  GeomTools_CurveSet                                      myCurves;
  GeomTools_Curve2dSet                                    myCurves2d;
  NCollection_IndexedMap<occ::handle<Standard_Transient>> myPolygons2D;
  NCollection_IndexedMap<occ::handle<Standard_Transient>> myPolygons3D;
  NCollection_IndexedDataMap<occ::handle<Poly_Triangulation>,
                             // clang-format off
                             bool> myTriangulations; //!< Contains a boolean flag with information
                                                                 //!  to save normals for triangulation
  // clang-format on
  NCollection_IndexedMap<occ::handle<Standard_Transient>> myNodes;
  bool                                                    myWithTriangles;
  bool                                                    myWithNormals;
  bool                                                    myOmitPCurvesOnPlane = false;
};

#endif // _BRepTools_ShapeSet_HeaderFile
