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


#include "RWStepDimTol_RWGeoTolAndGeoTolWthDefUnit.pxx"

#include <Interface_Check.hxx>
#include <Interface_EntityIterator.hxx>
#include <StepData_StepReaderData.hxx>
#include <StepData_StepWriter.hxx>
#include <StepDimTol_GeoTolAndGeoTolWthDefUnit.hxx>
#include <StepDimTol_GeometricToleranceTarget.hxx>
#include <TCollection_HAsciiString.hxx>

#include <cstring>

namespace
{
// The tolerances' own entities, in StepDimTol_GeometricToleranceType's order
const char* const THE_TYPE_NAMES[] = {"ANGULARITY_TOLERANCE",
                                      "CIRCULAR_RUNOUT_TOLERANCE",
                                      "COAXIALITY_TOLERANCE",
                                      "CONCENTRICITY_TOLERANCE",
                                      "CYLINDRICITY_TOLERANCE",
                                      "FLATNESS_TOLERANCE",
                                      "LINE_PROFILE_TOLERANCE",
                                      "PARALLELISM_TOLERANCE",
                                      "PERPENDICULARITY_TOLERANCE",
                                      "POSITION_TOLERANCE",
                                      "ROUNDNESS_TOLERANCE",
                                      "STRAIGHTNESS_TOLERANCE",
                                      "SURFACE_PROFILE_TOLERANCE",
                                      "SYMMETRY_TOLERANCE",
                                      "TOTAL_RUNOUT_TOLERANCE"};
constexpr int THE_NB_TYPES = int(sizeof(THE_TYPE_NAMES) / sizeof(THE_TYPE_NAMES[0]));

// Before GEOMETRIC_TOLERANCE in a complex, its components sorted by name
bool isBeforeGeometricTolerance(const StepDimTol_GeometricToleranceType theType)
{
  return std::strcmp(THE_TYPE_NAMES[theType], "GEOMETRIC_TOLERANCE") < 0;
}
} // namespace

//=================================================================================================

RWStepDimTol_RWGeoTolAndGeoTolWthDefUnit::RWStepDimTol_RWGeoTolAndGeoTolWthDefUnit() = default;

//=================================================================================================

void RWStepDimTol_RWGeoTolAndGeoTolWthDefUnit::ReadStep(
  const occ::handle<StepData_StepReaderData>&              data,
  const int                                                num0,
  occ::handle<Interface_Check>&                            ach,
  const occ::handle<StepDimTol_GeoTolAndGeoTolWthDefUnit>& ent) const
{
  int num = 0;
  data->NamedForComplex("GEOMETRIC_TOLERANCE", "GMTTLR", num0, num, ach);
  if (!data->CheckNbParams(num, 4, ach, "geometric_tolerance"))
  {
    return;
  }
  occ::handle<TCollection_HAsciiString> aName;
  data->ReadString(num, 1, "name", ach, aName);
  occ::handle<TCollection_HAsciiString> aDescription;
  data->ReadString(num, 2, "description", ach, aDescription);
  occ::handle<Standard_Transient> aMagnitude;
  data->ReadEntity(num, 3, "magnitude", ach, STANDARD_TYPE(Standard_Transient), aMagnitude);
  StepDimTol_GeometricToleranceTarget aTarget;
  data->ReadEntity(num, 4, "toleranced_shape_aspect", ach, aTarget);
  ent->Init(aName, aDescription, aMagnitude, aTarget);

  NCollection_Sequence<TCollection_AsciiString> aTypes;
  data->ComplexType(num0, aTypes);
  bool isType = false;
  for (const TCollection_AsciiString& aType : aTypes)
  {
    if (aType.IsEqual("GEOMETRIC_TOLERANCE_WITH_DEFINED_AREA_UNIT"))
    {
      ent->SetPerArea(true);
    }
    for (int i = 0; i < THE_NB_TYPES; i++)
    {
      if (aType.IsEqual(THE_TYPE_NAMES[i]))
      {
        ent->SetToleranceType(StepDimTol_GeometricToleranceType(i));
        isType = true;
      }
    }
  }
  if (!isType)
  {
    ach->AddFail("The type of geometric tolerance is not supported");
  }

  if (ent->IsPerArea())
  {
    data->NamedForComplex("GEOMETRIC_TOLERANCE_WITH_DEFINED_AREA_UNIT", num0, num, ach);
    if (data->CheckNbParams(num, 2, ach, "geometric_tolerance_with_defined_area_unit"))
    {
      if (data->ParamType(num, 1) == Interface_ParamEnum)
      {
        const char* aText = data->ParamCValue(num, 1);
        if (std::strcmp(aText, ".CIRCULAR.") == 0)
        {
          ent->SetAreaType(StepDimTol_Circular);
        }
        else if (std::strcmp(aText, ".RECTANGULAR.") == 0)
        {
          ent->SetAreaType(StepDimTol_Rectangular);
        }
        else if (std::strcmp(aText, ".SQUARE.") == 0)
        {
          ent->SetAreaType(StepDimTol_Square);
        }
        else
        {
          ach->AddFail("Parameter #1 (area_type) has not allowed value");
        }
      }
      else
      {
        ach->AddFail("Parameter #1 (area_type) is not enumeration");
      }
      occ::handle<Standard_Transient> aSecond;
      if (data->IsParamDefined(num, 2))
      {
        data->ReadEntity(num,
                         2,
                         "second_unit_size",
                         ach,
                         STANDARD_TYPE(Standard_Transient),
                         aSecond);
      }
      ent->SetSecondUnitSize(aSecond);
    }
  }

  data->NamedForComplex("GEOMETRIC_TOLERANCE_WITH_DEFINED_UNIT", num0, num, ach);
  if (data->CheckNbParams(num, 1, ach, "geometric_tolerance_with_defined_unit"))
  {
    occ::handle<Standard_Transient> aSize;
    data->ReadEntity(num, 1, "unit_size", ach, STANDARD_TYPE(Standard_Transient), aSize);
    ent->SetUnitSize(aSize);
  }
}

//=================================================================================================

void RWStepDimTol_RWGeoTolAndGeoTolWthDefUnit::WriteStep(
  StepData_StepWriter&                                     SW,
  const occ::handle<StepDimTol_GeoTolAndGeoTolWthDefUnit>& ent) const
{
  const StepDimTol_GeometricToleranceType aType = ent->GetToleranceType();
  if (isBeforeGeometricTolerance(aType))
  {
    SW.StartEntity(THE_TYPE_NAMES[aType]);
  }
  SW.StartEntity("GEOMETRIC_TOLERANCE");
  SW.Send(ent->Name());
  SW.Send(ent->Description());
  SW.Send(ent->Magnitude());
  SW.Send(ent->TolerancedShapeAspect().Value());
  if (ent->IsPerArea())
  {
    SW.StartEntity("GEOMETRIC_TOLERANCE_WITH_DEFINED_AREA_UNIT");
    switch (ent->AreaType())
    {
      case StepDimTol_Circular:
        SW.SendEnum(".CIRCULAR.");
        break;
      case StepDimTol_Rectangular:
        SW.SendEnum(".RECTANGULAR.");
        break;
      case StepDimTol_Square:
        SW.SendEnum(".SQUARE.");
        break;
    }
    if (ent->SecondUnitSize().IsNull())
    {
      SW.SendUndef();
    }
    else
    {
      SW.Send(ent->SecondUnitSize());
    }
  }
  SW.StartEntity("GEOMETRIC_TOLERANCE_WITH_DEFINED_UNIT");
  SW.Send(ent->UnitSize());
  if (!isBeforeGeometricTolerance(aType))
  {
    SW.StartEntity(THE_TYPE_NAMES[aType]);
  }
}

//=================================================================================================

void RWStepDimTol_RWGeoTolAndGeoTolWthDefUnit::Share(
  const occ::handle<StepDimTol_GeoTolAndGeoTolWthDefUnit>& ent,
  Interface_EntityIterator&                                iter) const
{
  iter.AddItem(ent->Magnitude());
  iter.AddItem(ent->TolerancedShapeAspect().Value());
  iter.AddItem(ent->UnitSize());
  if (!ent->SecondUnitSize().IsNull())
  {
    iter.AddItem(ent->SecondUnitSize());
  }
}
