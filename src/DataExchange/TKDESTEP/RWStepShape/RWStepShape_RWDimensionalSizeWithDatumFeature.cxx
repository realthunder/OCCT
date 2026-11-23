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

#include "RWStepShape_RWDimensionalSizeWithDatumFeature.pxx"

#include <Interface_Check.hxx>
#include <Interface_EntityIterator.hxx>
#include <StepData_StepReaderData.hxx>
#include <StepData_StepWriter.hxx>
#include <StepRepr_ProductDefinitionShape.hxx>
#include <StepShape_DimensionalSizeWithDatumFeature.hxx>
#include <TCollection_HAsciiString.hxx>

//=================================================================================================

RWStepShape_RWDimensionalSizeWithDatumFeature::RWStepShape_RWDimensionalSizeWithDatumFeature() =
  default;

//=================================================================================================

void RWStepShape_RWDimensionalSizeWithDatumFeature::ReadStep(
  const occ::handle<StepData_StepReaderData>&                   data,
  const int                                                     num,
  occ::handle<Interface_Check>&                                 ach,
  const occ::handle<StepShape_DimensionalSizeWithDatumFeature>& ent) const
{
  if (!data->CheckNbParams(num, 6, ach, "dimensional_size_with_datum_feature"))
  {
    return;
  }
  // Of shape_aspect, as datum_feature has them
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

  // Of dimensional_size: applies_to, which is the entity itself, and name
  occ::handle<Standard_Transient> anAppliesTo;
  data->ReadEntity(num,
                   5,
                   "dimensional_size.applies_to",
                   ach,
                   STANDARD_TYPE(Standard_Transient),
                   anAppliesTo);
  if (!anAppliesTo.IsNull() && anAppliesTo != ent)
  {
    ach->AddWarning("Parameter #5 (applies_to) is not the dimensional size itself");
  }
  occ::handle<TCollection_HAsciiString> aSizeName;
  data->ReadString(num, 6, "dimensional_size.name", ach, aSizeName);
  ent->SetSizeName(aSizeName);
}

//=================================================================================================

void RWStepShape_RWDimensionalSizeWithDatumFeature::WriteStep(
  StepData_StepWriter&                                          SW,
  const occ::handle<StepShape_DimensionalSizeWithDatumFeature>& ent) const
{
  SW.Send(ent->Name());
  SW.Send(ent->Description());
  SW.Send(ent->OfShape());
  SW.SendLogical(ent->ProductDefinitional());
  SW.Send(ent);
  SW.Send(ent->SizeName());
}

//=================================================================================================

void RWStepShape_RWDimensionalSizeWithDatumFeature::Share(
  const occ::handle<StepShape_DimensionalSizeWithDatumFeature>& ent,
  Interface_EntityIterator&                                     iter) const
{
  iter.AddItem(ent->OfShape());
}
