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
#include <TopoDS_Shape.hxx>

#include <Standard_Dump.hxx>

#include <mutex>
#include <unordered_map>

IMPLEMENT_STANDARD_RTTIEXT(TopoDS_TShape, Standard_Transient)

namespace
{
// Thawed copy -> its original (TopoDS_TShape::Thaw). A side table, so the flag
// is all a TShape carries: copies are rare, and the layout stays as it is.
std::mutex& thawedMutex()
{
  static std::mutex aMutex;
  return aMutex;
}

std::unordered_map<const TopoDS_TShape*, occ::handle<TopoDS_TShape>>& thawedOriginals()
{
  static std::unordered_map<const TopoDS_TShape*, occ::handle<TopoDS_TShape>> aMap;
  return aMap;
}
} // namespace

//=================================================================================================

void TopoDS_TShape::Thaw(const occ::handle<TopoDS_TShape>& theCopy,
                         const occ::handle<TopoDS_TShape>& theOriginal)
{
  if (theCopy.IsNull() || theOriginal.IsNull() || theCopy == theOriginal)
  {
    return;
  }
  std::lock_guard<std::mutex> aLock(thawedMutex());
  thawedOriginals()[theCopy.get()] = theOriginal;
  theCopy->setBit(Bit_Thawed, true);
}

//=================================================================================================

occ::handle<TopoDS_TShape> TopoDS_TShape::ThawedFrom(const TopoDS_TShape* theCopy)
{
  if (theCopy == nullptr || !theCopy->Thawed())
  {
    return occ::handle<TopoDS_TShape>();
  }
  std::lock_guard<std::mutex> aLock(thawedMutex());
  auto anIt = thawedOriginals().find(theCopy);
  return anIt == thawedOriginals().end() ? occ::handle<TopoDS_TShape>() : anIt->second;
}

//=================================================================================================

TopoDS_TShape::~TopoDS_TShape()
{
  if (Thawed())
  {
    occ::handle<TopoDS_TShape> anOriginal; // released after the lock
    std::lock_guard<std::mutex> aLock(thawedMutex());
    auto anIt = thawedOriginals().find(this);
    if (anIt != thawedOriginals().end())
    {
      anOriginal = anIt->second;
      thawedOriginals().erase(anIt);
    }
  }
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
