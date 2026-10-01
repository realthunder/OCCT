// Copyright (c) 2026 realthunder
//
// This file is part of Open CASCADE Technology software library.
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.

#include <BRep_RepresentationLock.hxx>

#include <TopoDS_Shape.hxx>
#include <TopoDS_TShape.hxx>

#include <cstdint>
#include <mutex>

namespace
{
constexpr std::size_t THE_NB_STRIPES = 64;

//! A table per kind, so that an edge's lock and a vertex's are never the
//! same stripe: the edge-then-vertex order then holds between the tables.
std::recursive_mutex& stripe(const TopoDS_TShape* theTShape)
{
  static std::recursive_mutex THE_VERTICES[THE_NB_STRIPES];
  static std::recursive_mutex THE_EDGES[THE_NB_STRIPES];
  static std::recursive_mutex THE_OTHERS[THE_NB_STRIPES];
  // TShapes are heap blocks: the low bits carry no information.
  const std::uintptr_t anAddr  = reinterpret_cast<std::uintptr_t>(theTShape);
  const std::size_t    anIndex = static_cast<std::size_t>((anAddr >> 4) ^ (anAddr >> 12)) % THE_NB_STRIPES;
  switch (theTShape->ShapeType())
  {
    case TopAbs_VERTEX:
      return THE_VERTICES[anIndex];
    case TopAbs_EDGE:
      return THE_EDGES[anIndex];
    default:
      return THE_OTHERS[anIndex];
  }
}
} // namespace

//=================================================================================================

BRep_RepresentationLock::BRep_RepresentationLock(const TopoDS_TShape* theTShape)
    : myMutex(nullptr)
{
  if (theTShape == nullptr || !theTShape->Immutable())
  {
    return;
  }
  std::recursive_mutex& aMutex = stripe(theTShape);
  aMutex.lock();
  myMutex = &aMutex;
}

//=================================================================================================

BRep_RepresentationLock::~BRep_RepresentationLock()
{
  if (myMutex != nullptr)
  {
    static_cast<std::recursive_mutex*>(myMutex)->unlock();
  }
}
