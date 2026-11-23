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

#ifndef _XCAFDimTolObjects_GeomToleranceUnitBasis_HeaderFile
#define _XCAFDimTolObjects_GeomToleranceUnitBasis_HeaderFile

//! What a tolerance's value holds on: its whole feature (None), or any part of it of a given
//! size -- per unit length ("0.02/15"), per unit area, a square, a rectangle or a circle
//! ("0.01/25x25", "0.01/25x30", "0.01/DIA 25"). AP242's geometric_tolerance_with_defined_unit
//! and geometric_tolerance_with_defined_area_unit.
enum XCAFDimTolObjects_GeomToleranceUnitBasis
{
  XCAFDimTolObjects_GeomToleranceUnitBasis_None,
  XCAFDimTolObjects_GeomToleranceUnitBasis_Length,
  XCAFDimTolObjects_GeomToleranceUnitBasis_Square,
  XCAFDimTolObjects_GeomToleranceUnitBasis_Rectangle,
  XCAFDimTolObjects_GeomToleranceUnitBasis_Circle
};

#endif // _XCAFDimTolObjects_GeomToleranceUnitBasis_HeaderFile
