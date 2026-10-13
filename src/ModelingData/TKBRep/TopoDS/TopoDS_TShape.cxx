// Created on: 1991-03-20
// Created by: Remi Lequette
// Copyright (c) 1991-1999 Matra Datavision
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

#include <TopoDS_TShape.hxx>
#include <TopoDS_LockedShape.hxx>
#include <TopoDS_Shape.hxx>

#include <Standard_Dump.hxx>

#include <NCollection_BaseAllocator.hxx>

IMPLEMENT_STANDARD_RTTIEXT(TopoDS_TShape, Standard_Transient)

//=================================================================================================

unsigned int TopoDS_LockedShape_Raised()
{
  return TopoDS_LockedShape::Raised();
}

namespace
{
// A thawed copy's original (TopoDS_TShape::Thaw), held by the copy: the
// allocator of its list of children. The flag and this are all a TShape
// carries, and the layout stays upstream's. The list lets go of the allocator
// with itself, whoever compiled the TShape's destructor -- upstream leaves it
// implicit, so with MSVC a client that makes a wire has one of its own, and
// the library is not told when the wire dies. A table beside the TShapes
// needed a destructor here to forget a copy, which is a symbol upstream's
// library does not have.
class ThawedOriginal : public NCollection_BaseAllocator
{
public:
  explicit ThawedOriginal(const occ::handle<TopoDS_TShape>& theOriginal)
      : myOriginal(theOriginal)
  {
  }

  occ::handle<TopoDS_TShape> myOriginal;
};
} // namespace

//=================================================================================================

void TopoDS_TShape::Thaw(const occ::handle<TopoDS_TShape>& theCopy,
                         const occ::handle<TopoDS_TShape>& theOriginal)
{
  if (theCopy.IsNull() || theOriginal.IsNull() || theCopy == theOriginal)
  {
    return;
  }
  // A list's nodes are freed by the allocator that made them: the children
  // the copy has already are taken out and put back under the new one. The
  // copy is the caller's alone until it hands it on, as any TShape it builds.
  NCollection_List<TopoDS_Shape> aChildren;
  for (NCollection_List<TopoDS_Shape>::Iterator anIt(theCopy->myShapes); anIt.More(); anIt.Next())
  {
    aChildren.Append(anIt.Value());
  }
  theCopy->myShapes.Clear(new ThawedOriginal(theOriginal));
  for (NCollection_List<TopoDS_Shape>::Iterator anIt(aChildren); anIt.More(); anIt.Next())
  {
    theCopy->myShapes.Append(anIt.Value());
  }
  theCopy->setBit(Bit_Thawed, true);
}

//=================================================================================================

occ::handle<TopoDS_TShape> TopoDS_TShape::ThawedFrom(const TopoDS_TShape* theCopy)
{
  if (theCopy == nullptr || !theCopy->Thawed())
  {
    return occ::handle<TopoDS_TShape>();
  }
  // The flag may be on a TShape that took it with another's state and was
  // never thawed itself: no original then.
  const ThawedOriginal* aState =
    dynamic_cast<const ThawedOriginal*>(theCopy->myShapes.Allocator().get());
  return aState != nullptr ? aState->myOriginal : occ::handle<TopoDS_TShape>();
}

//=================================================================================================

void TopoDS_TShape_ThawedFrom(const TopoDS_TShape* theCopy, occ::handle<TopoDS_TShape>& theOriginal)
{
  theOriginal = TopoDS_TShape::ThawedFrom(theCopy);
}

//=================================================================================================

void TopoDS_TShape::DumpJson(Standard_OStream& theOStream, int) const
{
  OCCT_DUMP_TRANSIENT_CLASS_BEGIN(theOStream)

  OCCT_DUMP_FIELD_VALUE_POINTER(theOStream, this)

  OCCT_DUMP_FIELD_VALUE_NUMERICAL(theOStream, ShapeType())
  OCCT_DUMP_FIELD_VALUE_NUMERICAL(theOStream, NbChildren())

  OCCT_DUMP_FIELD_VALUE_NUMERICAL(theOStream, myState)

  OCCT_DUMP_FIELD_VALUE_NUMERICAL(theOStream, Free())
  OCCT_DUMP_FIELD_VALUE_NUMERICAL(theOStream, Locked())
  OCCT_DUMP_FIELD_VALUE_NUMERICAL(theOStream, Modified())
  OCCT_DUMP_FIELD_VALUE_NUMERICAL(theOStream, Checked())

  OCCT_DUMP_FIELD_VALUE_NUMERICAL(theOStream, Orientable())
  OCCT_DUMP_FIELD_VALUE_NUMERICAL(theOStream, Closed())
  OCCT_DUMP_FIELD_VALUE_NUMERICAL(theOStream, Infinite())
  OCCT_DUMP_FIELD_VALUE_NUMERICAL(theOStream, Convex())
}
