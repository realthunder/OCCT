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

#include <MoniTool_Macros.hxx>
#include <Standard_Transient.hxx>
#include <StepData_UndefinedEntity.hxx>
#include <StepVisual_DraughtingCallout.hxx>
#include <StepVisual_InvisibleItem.hxx>
#include <StepVisual_PresentationLayerAssignment.hxx>
#include <StepVisual_PresentationRepresentation.hxx>
#include <StepVisual_StyledItem.hxx>
#include <TCollection_AsciiString.hxx>

StepVisual_InvisibleItem::StepVisual_InvisibleItem() = default;

int StepVisual_InvisibleItem::CaseNum(const occ::handle<Standard_Transient>& ent) const
{
  if (ent.IsNull())
  {
    return 0;
  }
  if (ent->IsKind(STANDARD_TYPE(StepVisual_StyledItem)))
  {
    return 1;
  }
  if (ent->IsKind(STANDARD_TYPE(StepVisual_PresentationLayerAssignment)))
  {
    return 2;
  }
  if (ent->IsKind(STANDARD_TYPE(StepVisual_PresentationRepresentation)))
  {
    return 3;
  }
  // A callout (AP242's invisible_item has draughting_callout: NIST FTC-09 hides 44)
  if (ent->IsKind(STANDARD_TYPE(StepVisual_DraughtingCallout)))
  {
    return 4;
  }
  // An annotation placeholder (AP242 ed2, a styled item no class is for: NIST STC-06 hides 19)
  occ::handle<StepData_UndefinedEntity> anUndefined = occ::down_cast<StepData_UndefinedEntity>(ent);
  if (!anUndefined.IsNull()
      && TCollection_AsciiString(anUndefined->StepType()).Search("PLACEHOLDER_OCCURRENCE") > 0)
  {
    return 5;
  }
  return 0;
}

occ::handle<StepVisual_StyledItem> StepVisual_InvisibleItem::StyledItem() const
{
  return GetCasted(StepVisual_StyledItem, Value());
}

occ::handle<StepVisual_PresentationLayerAssignment> StepVisual_InvisibleItem::
  PresentationLayerAssignment() const
{
  return GetCasted(StepVisual_PresentationLayerAssignment, Value());
}

occ::handle<StepVisual_PresentationRepresentation> StepVisual_InvisibleItem::
  PresentationRepresentation() const
{
  return GetCasted(StepVisual_PresentationRepresentation, Value());
}
