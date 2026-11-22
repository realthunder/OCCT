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

#ifndef _StepAP242_DraughtingModelItemAssociationWithPlaceholder_HeaderFile
#define _StepAP242_DraughtingModelItemAssociationWithPlaceholder_HeaderFile

#include <Standard.hxx>

#include <StepAP242_DraughtingModelItemAssociation.hxx>

//! A draughting model item association that names, besides the graphic, the
//! placeholder in it where a reader may write the PMI's own text (AP242 ed2,
//! draughting_model_item_association_with_placeholder). Its definition is the
//! semantic PMI the graphic presents, whatever its name: NIST's STC-07 names
//! them 'PPH - ...', CTC-02 ''. The placeholder, an
//! annotation_placeholder_occurrence this toolkit does not model, is kept as
//! read.
class StepAP242_DraughtingModelItemAssociationWithPlaceholder
    : public StepAP242_DraughtingModelItemAssociation
{

public:
  Standard_EXPORT StepAP242_DraughtingModelItemAssociationWithPlaceholder();

  //! The annotation placeholder occurrence (or the entity read for it)
  const occ::handle<Standard_Transient>& AnnotationPlaceholder() const
  {
    return myAnnotationPlaceholder;
  }

  void SetAnnotationPlaceholder(const occ::handle<Standard_Transient>& thePlaceholder)
  {
    myAnnotationPlaceholder = thePlaceholder;
  }

  DEFINE_STANDARD_RTTIEXT(StepAP242_DraughtingModelItemAssociationWithPlaceholder,
                          StepAP242_DraughtingModelItemAssociation)

private:
  occ::handle<Standard_Transient> myAnnotationPlaceholder;
};
#endif // _StepAP242_DraughtingModelItemAssociationWithPlaceholder_HeaderFile
