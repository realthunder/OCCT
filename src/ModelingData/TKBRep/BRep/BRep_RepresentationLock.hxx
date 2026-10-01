// Copyright (c) 2026 realthunder
//
// This file is part of Open CASCADE Technology software library.
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.

#ifndef _BRep_RepresentationLock_HeaderFile
#define _BRep_RepresentationLock_HeaderFile

#include <Standard.hxx>

class TopoDS_TShape;

//! Holds the lock on an Immutable (a fork flag) vertex's, edge's or face's
//! representations while it is in scope; for any other TShape it does
//! nothing.
//!
//! A frozen shape is a value: another thread may read it -- FreeCAD's
//! transaction log writes values on a worker thread while the document goes
//! on. But a frozen TShape still takes caches (BRep_Builder: the edge
//! polygons and face triangulations a mesher writes, a pcurve or a vertex
//! parameter for a new face, the sweep of the caches of faces that are gone,
//! BRepTools::Clean), and each of those edits the very lists a writer walks
//! (BRep_TEdge::Curves(), BRep_TVertex::Points(), the face's
//! triangulations): a node freed under the walker. Every edit of an
//! Immutable TShape's representations and every walk of them by a writer
//! (BRepTools_ShapeSet, BinTools_ShapeSet, BinTools_ShapeWriter) holds this
//! lock.
//!
//! The locks are striped by TShape address, a table per kind, and
//! recursive: the builder's methods call one another on the same shape. A
//! thread may take an edge's lock and then one of its vertices' (the sweep
//! does), never the other way round, and never two of one kind; the writers
//! take one at a time. A thread that keeps to that cannot deadlock.
//!
//! No class changes size for this: the table lives here, in the library.
class BRep_RepresentationLock
{
public:
  //! Locks <theTShape>'s representations if it is Immutable.
  Standard_EXPORT explicit BRep_RepresentationLock(const TopoDS_TShape* theTShape);

  Standard_EXPORT ~BRep_RepresentationLock();

  BRep_RepresentationLock(const BRep_RepresentationLock&)            = delete;
  BRep_RepresentationLock& operator=(const BRep_RepresentationLock&) = delete;

private:
  void* myMutex; //!< the stripe held, or null
};

#endif // _BRep_RepresentationLock_HeaderFile
