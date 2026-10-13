// Created on: 1990-12-11
// Created by: Remi Lequette
// Copyright (c) 1990-1999 Matra Datavision
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

#ifndef _TopoDS_LockedShape_HeaderFile
#define _TopoDS_LockedShape_HeaderFile

#include <NCollection_ForkState.hxx>
#include <Standard_DefineException.hxx>
#include <Standard_SStream.hxx>
#include <Standard_DomainError.hxx>

//! An attempt was made to modify a geometry of Shape already
//! shared or protected.

#if !defined No_Exception && !defined No_TopoDS_LockedShape
  #define TopoDS_LockedShape_Raise_if(CONDITION, MESSAGE)                                          \
    if (CONDITION)                                                                                 \
      throw TopoDS_LockedShape(MESSAGE);
#else
  #define TopoDS_LockedShape_Raise_if(CONDITION, MESSAGE)
#endif

//! Each one raised on this thread is counted (Raised()): a caller that has
//! the refusal turned into an error message by the algorithm it called -- as
//! FreeCAD's features do -- can still tell a frozen input was refused, and
//! run again on a copy.
//!
//! The count is the fork's library's, reached through NCollection_ForkState
//! and not through a function of this class: the tolerance setters of
//! BRep_TVertex, BRep_TEdge and BRep_TFace are inline and raise one, so a
//! client that calls a setter constructs it, and a counter this class
//! exported would be a symbol upstream's library does not have in every
//! such client. With upstream's library nothing is counted, and nothing is
//! frozen to be refused.
class TopoDS_LockedShape : public Standard_DomainError
{
public:
  TopoDS_LockedShape(const char* theMessage = "")
      : Standard_DomainError(theMessage)
  {
    count();
  }

  TopoDS_LockedShape(const char* theMessage, const char* theStackTrace)
      : Standard_DomainError(theMessage, theStackTrace)
  {
    count();
  }

  const char* ExceptionType() const noexcept override { return "TopoDS_LockedShape"; }

  //! How many have been raised on this thread; 0 where the library is
  //! upstream's.
  static unsigned int Raised()
  {
    NCollection_ForkState* aFork = NCollection_ForkState::Get();
    return aFork != nullptr ? aFork->LockedShapesRaised() : 0;
  }

private:
  static void count()
  {
    if (NCollection_ForkState* aFork = NCollection_ForkState::Get())
    {
      ++aFork->LockedShapesRaised();
    }
  }
};

//! TopoDS_LockedShape::Raised as a C function, for a client to look up by
//! name at run time (dlsym, GetProcAddress): found, the library is the
//! fork's and counts.
extern "C" Standard_EXPORT unsigned int TopoDS_LockedShape_Raised();

#endif // _TopoDS_LockedShape_HeaderFile
