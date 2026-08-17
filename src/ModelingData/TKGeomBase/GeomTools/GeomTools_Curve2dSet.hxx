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

#include <cstddef>
#include <unordered_map>
#include <utility>
#include <vector>

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
  Standard_EXPORT int Add(const occ::handle<Geom2d_Curve>& C);

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
  NCollection_IndexedMap<occ::handle<Standard_Transient>> myMap;
  //! Hash of a curve's written form -> the indices in myMap that hash to it.
  //! Only Add() maintains this; Read() fills myMap directly, so a file whose
  //! table does hold equal entries still reads back with its indices intact.
  std::unordered_map<std::size_t, std::vector<int>> myByValue;
  //! A curve Add() merged into an equal entry -> that entry's index. Index() is
  //! how the shape records are written, and it looks up by handle, so a merged
  //! curve has to be answerable even though the map holds its twin instead. The
  //! handle is kept so the address remains this curve's for the set's lifetime.
  std::unordered_map<const Standard_Transient*, std::pair<occ::handle<Standard_Transient>, int>>
    myAlias;
};

#endif // _GeomTools_Curve2dSet_HeaderFile
