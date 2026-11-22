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

#include "RWStepAP242_RWDraughtingModelItemAssociationWithPlaceholder.pxx"
#include "RWStepAP242_RWDraughtingModelItemAssociation.pxx"
#include <Interface_EntityIterator.hxx>
#include <StepData_StepReaderData.hxx>
#include <StepData_StepWriter.hxx>
#include <StepAP242_DraughtingModelItemAssociationWithPlaceholder.hxx>
#include <StepRepr_RepresentationItem.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_HArray1.hxx>
#include <StepRepr_Representation.hxx>

//=================================================================================================

RWStepAP242_RWDraughtingModelItemAssociationWithPlaceholder::
  RWStepAP242_RWDraughtingModelItemAssociationWithPlaceholder() = default;

//=================================================================================================

void RWStepAP242_RWDraughtingModelItemAssociationWithPlaceholder::ReadStep(
  const occ::handle<StepData_StepReaderData>&                                 theData,
  const int                                                                   theNum,
  occ::handle<Interface_Check>&                                               theCheck,
  const occ::handle<StepAP242_DraughtingModelItemAssociationWithPlaceholder>& theEnt) const
{
  if (!theData->CheckNbParams(theNum,
                              6,
                              theCheck,
                              "draughting_model_item_association_with_placeholder"))
  {
    return;
  }

  // Inherited fields of ItemIdentifiedRepresentationUsage
  occ::handle<TCollection_HAsciiString> aName;
  theData->ReadString(theNum, 1, "item_identified_representation_usage.name", theCheck, aName);

  occ::handle<TCollection_HAsciiString> aDescription;
  if (theData->IsParamDefined(theNum, 2))
  {
    theData->ReadString(theNum,
                        2,
                        "item_identified_representation_usage.description",
                        theCheck,
                        aDescription);
  }

  StepAP242_ItemIdentifiedRepresentationUsageDefinition aDefinition;
  theData->ReadEntity(theNum,
                      3,
                      "item_identified_representation_usage.definition",
                      theCheck,
                      aDefinition);

  occ::handle<StepRepr_Representation> aRepresentation;
  theData->ReadEntity(theNum,
                      4,
                      "item_identified_representation_usage.used_representation",
                      theCheck,
                      STANDARD_TYPE(StepRepr_Representation),
                      aRepresentation);

  occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> anItems;
  occ::handle<StepRepr_RepresentationItem>                                   anItem;
  int                                                                        aSub = 0;
  if (theData->ParamType(theNum, 5) == Interface_ParamIdent)
  {
    theData->ReadEntity(theNum,
                        5,
                        "item_identified_representation_usage.identified_item",
                        theCheck,
                        STANDARD_TYPE(StepRepr_RepresentationItem),
                        anItem);
    anItems = new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, 1);
    anItems->SetValue(1, anItem);
  }
  else if (theData->ReadSubList(theNum,
                                5,
                                "item_identified_representation_usage.identified_item",
                                theCheck,
                                aSub))
  {
    const int aNb = theData->NbParams(aSub);
    anItems       = new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, aNb);
    for (int i = 1; i <= aNb; i++)
    {
      if (theData->ReadEntity(aSub,
                              i,
                              "representation_item",
                              theCheck,
                              STANDARD_TYPE(StepRepr_RepresentationItem),
                              anItem))
      {
        anItems->SetValue(i, anItem);
      }
    }
  }

  // Own field: the placeholder, an entity of a type this toolkit may not know
  occ::handle<Standard_Transient> aPlaceholder;
  theData->ReadEntity(theNum,
                      6,
                      "draughting_model_item_association_with_placeholder.annotation_placeholder",
                      theCheck,
                      STANDARD_TYPE(Standard_Transient),
                      aPlaceholder);

  theEnt->Init(aName, aDescription, aDefinition, aRepresentation, anItems);
  theEnt->SetAnnotationPlaceholder(aPlaceholder);
}

//=================================================================================================

void RWStepAP242_RWDraughtingModelItemAssociationWithPlaceholder::WriteStep(
  StepData_StepWriter&                                                        theSW,
  const occ::handle<StepAP242_DraughtingModelItemAssociationWithPlaceholder>& theEnt) const
{
  RWStepAP242_RWDraughtingModelItemAssociation().WriteStep(theSW, theEnt);
  theSW.Send(theEnt->AnnotationPlaceholder());
}

//=================================================================================================

void RWStepAP242_RWDraughtingModelItemAssociationWithPlaceholder::Share(
  const occ::handle<StepAP242_DraughtingModelItemAssociationWithPlaceholder>& theEnt,
  Interface_EntityIterator&                                                   theIter) const
{
  RWStepAP242_RWDraughtingModelItemAssociation().Share(theEnt, theIter);
  if (!theEnt->AnnotationPlaceholder().IsNull())
  {
    theIter.AddItem(theEnt->AnnotationPlaceholder());
  }
}
