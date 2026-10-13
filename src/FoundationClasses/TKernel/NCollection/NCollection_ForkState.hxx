// Copyright (c) 2026 realthunder
//
// This file is part of Open CASCADE Technology software library.
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.

#ifndef _NCollection_ForkState_HeaderFile
#define _NCollection_ForkState_HeaderFile

#include <NCollection_BaseAllocator.hxx>

//! What the fork's library keeps for the whole process, behind its common
//! allocator (realthunder): NCollection_BaseAllocator::CommonBaseAllocator()
//! is of this class in the fork's TKernel, and of NCollection_BaseAllocator
//! itself in upstream's.
//!
//! It is how an inline function of a header reaches state of the fork's
//! library. Such a function is compiled into the clients as well, and a
//! client built against the fork's headers is to load where the libraries are
//! upstream's: a function of the fork's it calls by name is a symbol
//! upstream's library does not have, for which the loader refuses the client
//! whole -- whether or not the call is ever made. Get() calls
//! CommonBaseAllocator(), which upstream's library has too, and what it gives
//! is asked through virtual functions, which are no symbols of the client's:
//! null with upstream's library, where there is nothing to keep.
//!
//! Which module's copy of the inline function runs makes no difference -- a
//! client's, a toolkit's, or one the loader chose for all of them -- as it
//! would to a body that differs between the library and a client.
//!
//! State of an object goes behind a collection of that object instead
//! (tests/fork/abi/README.md); this is for what belongs to no object.
class NCollection_ForkState : public NCollection_BaseAllocator
{
public:
  //! The fork's state, or null where the library is upstream's.
  static NCollection_ForkState* Get()
  {
    return dynamic_cast<NCollection_ForkState*>(
      NCollection_BaseAllocator::CommonBaseAllocator().get());
  }

  //! How many TopoDS_LockedShape have been raised on the calling thread
  //! (TopoDS_LockedShape::Raised()), to be read and incremented.
  virtual unsigned int& LockedShapesRaised() = 0;
};

#endif // _NCollection_ForkState_HeaderFile
