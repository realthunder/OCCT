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

#ifndef _StepShape_DimensionalSizeWithDatumFeature_HeaderFile
#define _StepShape_DimensionalSizeWithDatumFeature_HeaderFile

#include <Standard.hxx>

#include <StepDimTol_DatumFeature.hxx>

class TCollection_HAsciiString;

//! A feature that is both a datum feature and the size of itself
//! (dimensional_size_with_datum_feature, AP242: a subtype of datum_feature
//! and dimensional_size): NIST STC-07's holes, datums B and C, each with its
//! diameter. Read as the datum feature, the shape aspect a datum and a shape
//! aspect relationship name; its size's applies_to is itself, and its name the
//! size's ('diameter').
class StepShape_DimensionalSizeWithDatumFeature : public StepDimTol_DatumFeature
{

public:
  Standard_EXPORT StepShape_DimensionalSizeWithDatumFeature();

  //! The dimensional size's name ('diameter', 'radius' ...)
  const occ::handle<TCollection_HAsciiString>& SizeName() const { return mySizeName; }

  void SetSizeName(const occ::handle<TCollection_HAsciiString>& theName) { mySizeName = theName; }

  DEFINE_STANDARD_RTTIEXT(StepShape_DimensionalSizeWithDatumFeature, StepDimTol_DatumFeature)

private:
  occ::handle<TCollection_HAsciiString> mySizeName;
};
#endif // _StepShape_DimensionalSizeWithDatumFeature_HeaderFile
