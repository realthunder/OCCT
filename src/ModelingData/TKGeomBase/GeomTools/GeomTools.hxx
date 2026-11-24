// Created on: 1992-08-28
// Created by: Remi LEQUETTE
// Copyright (c) 1992-1999 Matra Datavision
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

#ifndef _GeomTools_HeaderFile
#define _GeomTools_HeaderFile

#include <Standard.hxx>
#include <Standard_DefineAlloc.hxx>
#include <Standard_Handle.hxx>

#include <Standard_OStream.hxx>
#include <Standard_IStream.hxx>
#include <Standard_Real.hxx>
#include <gp_Dir.hxx>
#include <gp_Dir2d.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax22d.hxx>
#include <gp_Ax3.hxx>

#include <cmath>
#include <cstring>
#include <limits>
#include <type_traits>
class Geom_Surface;
class Geom_Curve;
class Geom2d_Curve;
class GeomTools_UndefinedTypeHandler;

//! The GeomTools package provides utilities for Geometry.
//!
//! * SurfaceSet, CurveSet, Curve2dSet : Tools used
//! for dumping, writing and reading.
//!
//! * Methods to dump, write, read curves and surfaces.
class GeomTools
{
public:
  DEFINE_STANDARD_ALLOC

  //! A set of Curves from Geom2d.
  //! Dumps the surface on the stream.
  Standard_EXPORT static void Dump(const occ::handle<Geom_Surface>& S, Standard_OStream& OS);

  //! Writes the surface on the stream.
  Standard_EXPORT static void Write(const occ::handle<Geom_Surface>& S, Standard_OStream& OS);

  //! Reads the surface from the stream.
  Standard_EXPORT static void Read(occ::handle<Geom_Surface>& S, Standard_IStream& IS);

  //! Dumps the Curve on the stream.
  Standard_EXPORT static void Dump(const occ::handle<Geom_Curve>& C, Standard_OStream& OS);

  //! Writes the Curve on the stream.
  Standard_EXPORT static void Write(const occ::handle<Geom_Curve>& C, Standard_OStream& OS);

  //! Reads the Curve from the stream.
  Standard_EXPORT static void Read(occ::handle<Geom_Curve>& C, Standard_IStream& IS);

  //! Dumps the Curve on the stream.
  Standard_EXPORT static void Dump(const occ::handle<Geom2d_Curve>& C, Standard_OStream& OS);

  //! Writes the Curve on the stream.
  Standard_EXPORT static void Write(const occ::handle<Geom2d_Curve>& C, Standard_OStream& OS);

  //! Reads the Curve from the stream.
  Standard_EXPORT static void Read(occ::handle<Geom2d_Curve>& C, Standard_IStream& IS);

  Standard_EXPORT static void SetUndefinedTypeHandler(
    const occ::handle<GeomTools_UndefinedTypeHandler>& aHandler);

  Standard_EXPORT static occ::handle<GeomTools_UndefinedTypeHandler> GetUndefinedTypeHandler();

  //! Reads the double value from the stream. Zero is read
  //! in case of error
  Standard_EXPORT static void GetReal(Standard_IStream& IS, double& theValue);

  //! Fork: the direction of coordinates read back. One written unit to its
  //! last bits -- its square norm 1 within 4 epsilon, as a transformation
  //! can leave a direction -- keeps them as written; any other is normalized,
  //! as gp_Dir does. gp_Dir keeps a vector as given only when its square norm
  //! is 1 within gp::Resolution(), so it normalized such a direction again,
  //! which moved it by an ulp: the shape read back was not the shape written,
  //! and what was made of it differed. Inline, so no symbol of the fork's.
  static gp_Dir DirAsWritten(const double theX, const double theY, const double theZ)
  {
    static_assert(std::is_trivially_copyable<gp_Dir>::value
                    && sizeof(gp_Dir) == 3 * sizeof(double),
                  "gp_Dir is three doubles");
    const double aSqMod = theX * theX + theY * theY + theZ * theZ;
    if (std::abs(aSqMod - 1.0) <= 4.0 * std::numeric_limits<double>::epsilon())
    {
      const double aXYZ[3] = {theX, theY, theZ};
      gp_Dir       aDir;
      std::memcpy(static_cast<void*>(&aDir), aXYZ, sizeof(aXYZ));
      if (aDir.X() == theX && aDir.Y() == theY && aDir.Z() == theZ)
      {
        return aDir;
      }
    }
    return gp_Dir(theX, theY, theZ);
  }

  //! Fork: the frame of a conic read back, made as gp_Ax2(theP, theA, theAX)
  //! makes it, and kept as written -- theAY too -- where what that made is
  //! the written frame to rounding (FrameNear()). gp_Ax2 makes its X and Y
  //! directions again from the main and X directions, which moved a frame a
  //! rotation left by an ulp, as DirAsWritten() is for a direction.
  static gp_Ax2 Ax2AsWritten(const gp_Pnt& theP,
                             const gp_Dir& theA,
                             const gp_Dir& theAX,
                             const gp_Dir& theAY)
  {
    return frameAsWritten(gp_Ax2(theP, theA, theAX), theP, theA, theAX, theAY);
  }

  //! Fork: Ax2AsWritten() for an elementary surface's frame, which may be
  //! left-handed: theAY's side of the main and X directions says which
  static gp_Ax3 Ax3AsWritten(const gp_Pnt& theP,
                             const gp_Dir& theA,
                             const gp_Dir& theAX,
                             const gp_Dir& theAY)
  {
    gp_Ax3 aMade(theP, theA, theAX);
    if (theAY.DotCross(theA, theAX) < 0)
    {
      aMade.YReverse();
    }
    return frameAsWritten(aMade, theP, theA, theAX, theAY);
  }

  //! Fork: Ax2AsWritten() in the plane, made as gp_Ax22d(theP, theAX, theAY)
  static gp_Ax22d Ax22dAsWritten(const gp_Pnt2d& theP, const gp_Dir2d& theAX, const gp_Dir2d& theAY)
  {
    static_assert(std::is_trivially_copyable<gp_Ax22d>::value
                    && sizeof(gp_Ax22d) == 6 * sizeof(double),
                  "gp_Ax22d is a point and two directions");
    const gp_Ax22d aMade(theP, theAX, theAY);
    if (!frameNear(aMade.XDirection().XY(), theAX.XY())
        || !frameNear(aMade.YDirection().XY(), theAY.XY()))
    {
      return aMade;
    }
    // point, vydir, vxdir
    const double aV[6] = {theP.X(), theP.Y(), theAY.X(), theAY.Y(), theAX.X(), theAX.Y()};
    gp_Ax22d     aRes;
    std::memcpy(static_cast<void*>(&aRes), aV, sizeof(aV));
    if (aRes.Location().X() == theP.X() && aRes.Location().Y() == theP.Y()
        && aRes.XDirection().X() == theAX.X() && aRes.XDirection().Y() == theAX.Y()
        && aRes.YDirection().X() == theAY.X() && aRes.YDirection().Y() == theAY.Y())
    {
      return aRes;
    }
    return aMade;
  }

  //! Fork: DirAsWritten() in the plane
  static gp_Dir2d Dir2dAsWritten(const double theX, const double theY)
  {
    static_assert(std::is_trivially_copyable<gp_Dir2d>::value
                    && sizeof(gp_Dir2d) == 2 * sizeof(double),
                  "gp_Dir2d is two doubles");
    const double aSqMod = theX * theX + theY * theY;
    if (std::abs(aSqMod - 1.0) <= 4.0 * std::numeric_limits<double>::epsilon())
    {
      const double aXY[2] = {theX, theY};
      gp_Dir2d     aDir;
      std::memcpy(static_cast<void*>(&aDir), aXY, sizeof(aXY));
      if (aDir.X() == theX && aDir.Y() == theY)
      {
        return aDir;
      }
    }
    return gp_Dir2d(theX, theY);
  }

private:
  //! Whether a direction made again is the written one to rounding: each
  //! coordinate within 16 epsilon
  static bool frameNear(const gp_XYZ& theMade, const gp_XYZ& theWritten)
  {
    const double aTol = 16.0 * std::numeric_limits<double>::epsilon();
    return std::abs(theMade.X() - theWritten.X()) <= aTol
           && std::abs(theMade.Y() - theWritten.Y()) <= aTol
           && std::abs(theMade.Z() - theWritten.Z()) <= aTol;
  }

  static bool frameNear(const gp_XY& theMade, const gp_XY& theWritten)
  {
    const double aTol = 16.0 * std::numeric_limits<double>::epsilon();
    return std::abs(theMade.X() - theWritten.X()) <= aTol
           && std::abs(theMade.Y() - theWritten.Y()) <= aTol;
  }

  static bool sameBits(const gp_XYZ& theA, const gp_XYZ& theB)
  {
    return theA.X() == theB.X() && theA.Y() == theB.Y() && theA.Z() == theB.Z();
  }

  //! theMade with the written frame put in as it was, where it is theMade to
  //! rounding: a gp_Ax2 or gp_Ax3, a point and three directions
  template <class Frame>
  static Frame frameAsWritten(const Frame&  theMade,
                              const gp_Pnt& theP,
                              const gp_Dir& theA,
                              const gp_Dir& theAX,
                              const gp_Dir& theAY)
  {
    static_assert(std::is_trivially_copyable<Frame>::value
                    && sizeof(Frame) == 12 * sizeof(double),
                  "a frame is a point and three directions");
    if (!frameNear(theMade.Direction().XYZ(), theA.XYZ())
        || !frameNear(theMade.XDirection().XYZ(), theAX.XYZ())
        || !frameNear(theMade.YDirection().XYZ(), theAY.XYZ()))
    {
      return theMade;
    }
    // axis (location, direction), vydir, vxdir
    const double aV[12] = {theP.X(),  theP.Y(),  theP.Z(),  theA.X(),  theA.Y(),  theA.Z(),
                           theAY.X(), theAY.Y(), theAY.Z(), theAX.X(), theAX.Y(), theAX.Z()};
    Frame        aRes;
    std::memcpy(static_cast<void*>(&aRes), aV, sizeof(aV));
    if (sameBits(aRes.Location().XYZ(), theP.XYZ()) && sameBits(aRes.Direction().XYZ(), theA.XYZ())
        && sameBits(aRes.XDirection().XYZ(), theAX.XYZ())
        && sameBits(aRes.YDirection().XYZ(), theAY.XYZ()))
    {
      return aRes;
    }
    return theMade;
  }
};

#endif // _GeomTools_HeaderFile
