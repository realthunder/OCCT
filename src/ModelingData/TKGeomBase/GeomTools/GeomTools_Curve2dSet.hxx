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

#ifndef _GeomTools_Curve2dSet_HeaderFile
#define _GeomTools_Curve2dSet_HeaderFile

#include <Standard.hxx>
#include <Standard_DefineAlloc.hxx>
#include <Standard_Handle.hxx>

#include <Standard_Transient.hxx>
#include <NCollection_IndexedMap.hxx>
#include <Standard_Integer.hxx>
#include <Standard_OStream.hxx>
#include <Standard_IStream.hxx>
#include <Message_ProgressRange.hxx>

class Geom2d_Curve;

//! Stores a set of Curves from Geom2d.
class GeomTools_Curve2dSet
{
public:
  DEFINE_STANDARD_ALLOC

  //! Returns an empty set of Curves.
  Standard_EXPORT GeomTools_Curve2dSet();

  //! Clears the content of the set.
  Standard_EXPORT void Clear();

  //! Incorporate a new Curve in the set and returns its index.
  //!
  //! A curve that would be written exactly as one already in the set is not
  //! added again: the index of that one is returned instead. The set is keyed
  //! by handle otherwise, so a file used to carry the same pcurve as often as
  //! independent computation happened to produce it -- on a real model the
  //! largest single duplication in a shape file.
  //!
  //! This is sound for 2D curves specifically. A pcurve is reached through the
  //! surface an edge's representation names, never by the identity of the curve
  //! object itself, so two representations sharing one Geom2d_Curve are
  //! indistinguishable from two holding equal copies. The 3D curve and surface
  //! sets are deliberately NOT deduplicated this way: a vertex's parameter is
  //! keyed on its edge's curve, and an edge's pcurve on its face's surface, so
  //! merging equal-but-distinct objects there would make those lookups
  //! ambiguous.
  //!
  //! Turned off with SetMerging(false), which restores keying by handle alone.
  Standard_EXPORT int Add(const occ::handle<Geom2d_Curve>& C);

  //! True if Add() merges a curve into an entry that would be written
  //! identically. On by default.
  Standard_EXPORT bool IsMerging() const;

  //! Define whether Add() merges equal curves. Set it before anything is added:
  //! turning it off part way leaves already merged entries merged, and turning
  //! it on part way only merges what comes after.
  Standard_EXPORT void SetMerging(const bool theMerging);

  //! Returns the Curve of index <I>.
  Standard_EXPORT occ::handle<Geom2d_Curve> Curve2d(const int I) const;

  //! Returns the index of <L>.
  Standard_EXPORT int Index(const occ::handle<Geom2d_Curve>& C) const;

  //! The number of entries. A subclass of BRepTools_ShapeSet that writes the
  //! tables itself needs it: the shape records index into this set, so the
  //! count is part of what has to be emitted and read back.
  int Extent() const { return myMap.Extent(); }

  //! Dumps the content of me on the stream <OS>.
  Standard_EXPORT void Dump(Standard_OStream& OS) const;

  //! Writes the content of me on the stream <OS> in a
  //! format that can be read back by Read.
  Standard_EXPORT void Write(
    Standard_OStream&            OS,
    const Message_ProgressRange& theProgress = Message_ProgressRange()) const;

  //! Reads the content of me from the stream <IS>.
  //! me is first cleared.
  Standard_EXPORT void Read(Standard_IStream&            IS,
                            const Message_ProgressRange& theProgress = Message_ProgressRange());

  //! Dumps the curve on the stream, if compact is True
  //! use the compact format that can be read back.
  Standard_EXPORT static void PrintCurve2d(const occ::handle<Geom2d_Curve>& C,
                                           Standard_OStream&                OS,
                                           const bool                       compact = false);

  //! Reads the curve from the stream. The curve is
  //! assumed to have been written with the Print
  //! method (compact = True).
  Standard_EXPORT static occ::handle<Geom2d_Curve> ReadCurve2d(Standard_IStream& IS);

private:
  //! The one member upstream's class has, and the class is to stay that
  //! size: a client built against upstream's library holds a set, or a
  //! BRepTools_ShapeSet with one in it, at upstream's size, and the
  //! constructor here must not write past it (conda-forge's ifcopenshell
  //! does, in ifcopenshell.geom.serialise; it crashed). What merging needs
  //! beyond the map is kept behind the map's allocator, see the source.
  NCollection_IndexedMap<occ::handle<Standard_Transient>> myMap;
};

//! GeomTools_Curve2dSet::SetMerging as a C function, for a client to look up
//! by name at run time (dlsym, GetProcAddress) and do without where the
//! library is upstream's, which merges nothing: a client that calls a
//! function of the fork by its C++ name does not load there.
extern "C" Standard_EXPORT void GeomTools_Curve2dSet_SetMerging(GeomTools_Curve2dSet& theSet,
                                                                const bool theMerging);

#endif // _GeomTools_Curve2dSet_HeaderFile
