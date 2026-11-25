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

#include "RWStepRepr_RWBasicRoundHoleOccurrence.pxx"

#include <Interface_Check.hxx>
#include <Interface_EntityIterator.hxx>
#include <StepData_StepReaderData.hxx>
#include <StepData_StepWriter.hxx>
#include <StepRepr_BasicRoundHoleOccurrence.hxx>
#include <StepRepr_ProductDefinitionShape.hxx>
#include <TCollection_HAsciiString.hxx>

//=================================================================================================

RWStepRepr_RWBasicRoundHoleOccurrence::RWStepRepr_RWBasicRoundHoleOccurrence() = default;

//=================================================================================================

void RWStepRepr_RWBasicRoundHoleOccurrence::ReadStep(
  const occ::handle<StepData_StepReaderData>&           data,
  const int                                             num,
  occ::handle<Interface_Check>&                         ach,
  const occ::handle<StepRepr_BasicRoundHoleOccurrence>& ent) const
{
  if (!data->CheckNbParams(num, 5, ach, "basic_round_hole_occurrence"))
  {
    return;
  }
  // Of shape_aspect
  occ::handle<TCollection_HAsciiString> aName;
  data->ReadString(num, 1, "shape_aspect.name", ach, aName);
  occ::handle<TCollection_HAsciiString> aDescription;
  if (data->IsParamDefined(num, 2))
  {
    data->ReadString(num, 2, "shape_aspect.description", ach, aDescription);
  }
  occ::handle<StepRepr_ProductDefinitionShape> anOfShape;
  data->ReadEntity(num,
                   3,
                   "shape_aspect.of_shape",
                   ach,
                   STANDARD_TYPE(StepRepr_ProductDefinitionShape),
                   anOfShape);
  StepData_Logical aProductDefinitional;
  data->ReadLogical(num, 4, "shape_aspect.product_definitional", ach, aProductDefinitional);
  ent->Init(aName, aDescription, anOfShape, aProductDefinitional);

  // Its definition, the basic_round_hole: an entity of any type, as read
  occ::handle<Standard_Transient> aDefinition;
  data->ReadEntity(num, 5, "definition", ach, STANDARD_TYPE(Standard_Transient), aDefinition);
  ent->SetDefinition(aDefinition);
}

//=================================================================================================

void RWStepRepr_RWBasicRoundHoleOccurrence::WriteStep(
  StepData_StepWriter&                                  SW,
  const occ::handle<StepRepr_BasicRoundHoleOccurrence>& ent) const
{
  SW.Send(ent->Name());
  // A description is optional here, as NIST writes it ($)
  if (ent->Description().IsNull())
  {
    SW.SendUndef();
  }
  else
  {
    SW.Send(ent->Description());
  }
  SW.Send(ent->OfShape());
  SW.SendLogical(ent->ProductDefinitional());
  SW.Send(ent->Definition());
}

//=================================================================================================

void RWStepRepr_RWBasicRoundHoleOccurrence::Share(
  const occ::handle<StepRepr_BasicRoundHoleOccurrence>& ent,
  Interface_EntityIterator&                             iter) const
{
  iter.AddItem(ent->OfShape());
  iter.AddItem(ent->Definition());
}
