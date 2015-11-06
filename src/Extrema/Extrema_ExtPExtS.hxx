// Created on: 1999-09-16
// Created by: Edward AGAPOV
// Copyright (c) 1999 Matra Datavision
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

#ifndef _Extrema_ExtPExtS_HeaderFile
#define _Extrema_ExtPExtS_HeaderFile

#include <Standard.hxx>
#include <Standard_Type.hxx>

#include <Standard_Real.hxx>
#include <Extrema_FuncExtPS.hxx>
#include <gp_Vec.hxx>
#include <gp_Ax2.hxx>
#include <Extrema_GenExtPS.hxx>
#include <Standard_Boolean.hxx>
#include <Standard_Integer.hxx>
#include <Extrema_POnSurf.hxx>
#include <Standard_Transient.hxx>
class Adaptor3d_HCurve;
class GeomAdaptor_HSurfaceOfLinearExtrusion;
class StdFail_NotDone;
class Standard_OutOfRange;
class gp_Pnt;
class Extrema_POnSurf;
class gp_Ax2;


class Extrema_ExtPExtS;
DEFINE_STANDARD_HANDLE(Extrema_ExtPExtS, Standard_Transient)

//! It calculates all the extremum (minimum and
//! maximum) distances between a point and a linear
//! extrusion surface.
class Extrema_ExtPExtS : public Standard_Transient
{

public:

  
  Standard_EXPORT Extrema_ExtPExtS();
  
  //! It calculates all the distances between a point
  //! from gp and a Surface.
  Standard_EXPORT Extrema_ExtPExtS(const gp_Pnt& P, const Handle(GeomAdaptor_HSurfaceOfLinearExtrusion)& S, const Standard_Real Umin, const Standard_Real Usup, const Standard_Real Vmin, const Standard_Real Vsup, const Standard_Real TolU, const Standard_Real TolV);
  
  //! It calculates all the distances between a point
  //! from gp and a Surface.
  Standard_EXPORT Extrema_ExtPExtS(const gp_Pnt& P, const Handle(GeomAdaptor_HSurfaceOfLinearExtrusion)& S, const Standard_Real TolU, const Standard_Real TolV);
  
  //! Initializes the fields of the algorithm.
  Standard_EXPORT void Initialize (const Handle(GeomAdaptor_HSurfaceOfLinearExtrusion)& S, const Standard_Real Uinf, const Standard_Real Usup, const Standard_Real Vinf, const Standard_Real Vsup, const Standard_Real TolU, const Standard_Real TolV);
  
  Standard_EXPORT void Perform (const gp_Pnt& P);
  
  //! Returns True if the distances are found.
  Standard_EXPORT Standard_Boolean IsDone() const;
  
  //! Returns the number of extremum distances.
  Standard_EXPORT Standard_Integer NbExt() const;
  
  //! Returns the value of the Nth resulting square distance.
  Standard_EXPORT Standard_Real SquareDistance (const Standard_Integer N) const;
  
  //! Returns the point of the Nth resulting distance.
  Standard_EXPORT const Extrema_POnSurf& Point (const Standard_Integer N) const;




  DEFINE_STANDARD_RTTI(Extrema_ExtPExtS,Standard_Transient)

protected:




private:

  
  Standard_EXPORT void MakePreciser (Standard_Real& U, const gp_Pnt& P, const Standard_Boolean isMin, const gp_Ax2& OrtogSection) const;

  Standard_Real myuinf;
  Standard_Real myusup;
  Standard_Real mytolu;
  Standard_Real myvinf;
  Standard_Real myvsup;
  Standard_Real mytolv;
  Extrema_FuncExtPS myF;
  Handle(Adaptor3d_HCurve) myC;
  Handle(GeomAdaptor_HSurfaceOfLinearExtrusion) myS;
  gp_Vec myDirection;
  gp_Ax2 myPosition;
  Extrema_GenExtPS myExtPS;
  Standard_Boolean myIsAnalyticallyComputable;
  Standard_Boolean myDone;
  Standard_Integer myNbExt;
  Standard_Real mySqDist[4];
  Extrema_POnSurf myPoint[4];


};







#endif // _Extrema_ExtPExtS_HeaderFile