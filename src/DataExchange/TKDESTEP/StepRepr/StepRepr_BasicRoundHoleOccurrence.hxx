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

#ifndef _StepRepr_BasicRoundHoleOccurrence_HeaderFile
#define _StepRepr_BasicRoundHoleOccurrence_HeaderFile

#include <Standard.hxx>

#include <StepRepr_ShapeAspect.hxx>

//! An occurrence of a round hole (basic_round_hole_occurrence, AP242 ed3: a shape aspect, and
//! the basic_round_hole it is one of): NIST STC-09's 25 holes, each identified with its faces
//! by an item_identified_representation_usage. Known as the shape aspect it is, a composite
//! shape aspect may have it as a member (STC-09's positions of a pattern of holes); its
//! definition, an entity the toolkit does not know, is kept as read.
class StepRepr_BasicRoundHoleOccurrence : public StepRepr_ShapeAspect
{

public:
  Standard_EXPORT StepRepr_BasicRoundHoleOccurrence();

  //! The basic_round_hole it is an occurrence of
  const occ::handle<Standard_Transient>& Definition() const { return myDefinition; }

  void SetDefinition(const occ::handle<Standard_Transient>& theDefinition)
  {
    myDefinition = theDefinition;
  }

  DEFINE_STANDARD_RTTIEXT(StepRepr_BasicRoundHoleOccurrence, StepRepr_ShapeAspect)

private:
  occ::handle<Standard_Transient> myDefinition;
};
#endif // _StepRepr_BasicRoundHoleOccurrence_HeaderFile
