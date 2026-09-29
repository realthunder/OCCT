// Created on: 1993-01-21
// Created by: Remi LEQUETTE
// Copyright (c) 1993-1999 Matra Datavision
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

#include <GeomTools.hxx>

#include <cctype>
#include <cstring>
#include <string>

#include <Geom2d_Curve.hxx>
#include <Geom_Surface.hxx>
#include <GeomTools_Curve2dSet.hxx>
#include <GeomTools_CurveSet.hxx>
#include <GeomTools_SurfaceSet.hxx>
#include <GeomTools_UndefinedTypeHandler.hxx>

static occ::handle<GeomTools_UndefinedTypeHandler> theActiveHandler =
  new GeomTools_UndefinedTypeHandler;

void GeomTools::Dump(const occ::handle<Geom_Surface>& S, Standard_OStream& OS)
{
  GeomTools_SurfaceSet::PrintSurface(S, OS);
}

void GeomTools::Write(const occ::handle<Geom_Surface>& S, Standard_OStream& OS)
{
  GeomTools_SurfaceSet::PrintSurface(S, OS, true);
}

void GeomTools::Read(occ::handle<Geom_Surface>& S, Standard_IStream& IS)
{
  S = GeomTools_SurfaceSet::ReadSurface(IS);
}

void GeomTools::Dump(const occ::handle<Geom_Curve>& C, Standard_OStream& OS)
{
  GeomTools_CurveSet::PrintCurve(C, OS);
}

void GeomTools::Write(const occ::handle<Geom_Curve>& C, Standard_OStream& OS)
{
  GeomTools_CurveSet::PrintCurve(C, OS, true);
}

void GeomTools::Read(occ::handle<Geom_Curve>& C, Standard_IStream& IS)
{
  C = GeomTools_CurveSet::ReadCurve(IS);
}

void GeomTools::Dump(const occ::handle<Geom2d_Curve>& C, Standard_OStream& OS)
{
  GeomTools_Curve2dSet::PrintCurve2d(C, OS);
}

void GeomTools::Write(const occ::handle<Geom2d_Curve>& C, Standard_OStream& OS)
{
  GeomTools_Curve2dSet::PrintCurve2d(C, OS, true);
}

void GeomTools::Read(occ::handle<Geom2d_Curve>& C, Standard_IStream& IS)
{
  C = GeomTools_Curve2dSet::ReadCurve2d(IS);
}

//=================================================================================================

void GeomTools::SetUndefinedTypeHandler(const occ::handle<GeomTools_UndefinedTypeHandler>& aHandler)
{
  if (!aHandler.IsNull())
  {
    theActiveHandler = aHandler;
  }
}

//=================================================================================================

occ::handle<GeomTools_UndefinedTypeHandler> GeomTools::GetUndefinedTypeHandler()
{
  return theActiveHandler;
}

//=================================================================================================

void GeomTools::GetReal(Standard_IStream& IS, double& theValue)
{
  theValue = 0.;
  if (IS.eof())
  {
    return;
  }
  // According IEEE-754 Specification and standard stream parameters
  // the most optimal buffer length not less then 25
  constexpr size_t THE_BUFFER_SIZE = 32;
  char             aBuffer[THE_BUFFER_SIZE];

  aBuffer[0]                = '\0';
  std::streamsize anOldWide = IS.width(THE_BUFFER_SIZE - 1);
  IS >> aBuffer;
  IS.width(anOldWide);
  // A real longer than the buffer is the rest of the same token, not the
  // next field. Writers with fixed notation produce such tokens -- an
  // unbounded edge's +-2e100 range as a 101-digit integer part -- and the
  // 256-byte buffer this one replaced read them whole; split, every field
  // after it came from the wrong token and the reader spun on a failed
  // stream. The common short token pays one strlen.
  // (width(N) extracts at most N - 1 characters: a full buffer holds
  // THE_BUFFER_SIZE - 2 of them.)
  if (std::strlen(aBuffer) >= THE_BUFFER_SIZE - 2 && IS.good())
  {
    const int aNext = IS.peek();
    if (aNext != std::char_traits<char>::eof() && !std::isspace(aNext))
    {
      std::string aToken(aBuffer);
      std::string aRest;
      IS >> aRest;
      aToken += aRest;
      theValue = Strtod(aToken.c_str(), nullptr);
      return;
    }
  }
  theValue = Strtod(aBuffer, nullptr);
}
