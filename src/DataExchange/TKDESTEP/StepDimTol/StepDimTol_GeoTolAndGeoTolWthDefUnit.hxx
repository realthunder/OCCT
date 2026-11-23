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

#ifndef _StepDimTol_GeoTolAndGeoTolWthDefUnit_HeaderFile
#define _StepDimTol_GeoTolAndGeoTolWthDefUnit_HeaderFile

#include <Standard.hxx>

#include <StepDimTol_AreaUnitType.hxx>
#include <StepDimTol_GeometricTolerance.hxx>
#include <StepDimTol_GeometricToleranceType.hxx>

//! A tolerance per unit length or area: the complex of a tolerance's type,
//! geometric_tolerance, geometric_tolerance_with_defined_unit and, per area,
//! geometric_tolerance_with_defined_area_unit -- a flatness 0.01 on any
//! 25 x 25 (NIST CTC-03, FTC-07, STC-07), a straightness 0.02 on any 15
//! (FTC-10, STC-10). A unit size is kept as read, as the magnitude is: a
//! length_measure_with_unit, or the complex with a measure_representation_item
//! NIST writes.
class StepDimTol_GeoTolAndGeoTolWthDefUnit : public StepDimTol_GeometricTolerance
{

public:
  Standard_EXPORT StepDimTol_GeoTolAndGeoTolWthDefUnit();

  StepDimTol_GeometricToleranceType GetToleranceType() const { return myToleranceType; }

  void SetToleranceType(const StepDimTol_GeometricToleranceType theType)
  {
    myToleranceType = theType;
  }

  //! The unit's length; per area, its first side or its diameter
  const occ::handle<Standard_Transient>& UnitSize() const { return myUnitSize; }

  void SetUnitSize(const occ::handle<Standard_Transient>& theSize) { myUnitSize = theSize; }

  //! Per area (geometric_tolerance_with_defined_area_unit), not per length
  bool IsPerArea() const { return myIsPerArea; }

  void SetPerArea(const bool theIsPerArea) { myIsPerArea = theIsPerArea; }

  StepDimTol_AreaUnitType AreaType() const { return myAreaType; }

  void SetAreaType(const StepDimTol_AreaUnitType theType) { myAreaType = theType; }

  //! A rectangle's second side; null for none
  const occ::handle<Standard_Transient>& SecondUnitSize() const { return mySecondUnitSize; }

  void SetSecondUnitSize(const occ::handle<Standard_Transient>& theSize)
  {
    mySecondUnitSize = theSize;
  }

  DEFINE_STANDARD_RTTIEXT(StepDimTol_GeoTolAndGeoTolWthDefUnit, StepDimTol_GeometricTolerance)

private:
  StepDimTol_GeometricToleranceType myToleranceType = StepDimTol_GTTFlatnessTolerance;
  occ::handle<Standard_Transient>   myUnitSize;
  bool                              myIsPerArea = false;
  StepDimTol_AreaUnitType           myAreaType  = StepDimTol_Square;
  occ::handle<Standard_Transient>   mySecondUnitSize;
};
#endif // _StepDimTol_GeoTolAndGeoTolWthDefUnit_HeaderFile
