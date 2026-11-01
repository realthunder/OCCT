// Copyright (c) 2008-2015 OPEN CASCADE SAS
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

#include <math.hxx>
#include <Standard_Assert.hxx>
#include <BRepGProp_Face.hxx>
#include <BRepGProp_Domain.hxx>
#include <BRepGProp_Gauss.hxx>
#include <BRep_Tool.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_OffsetSurface.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_SurfaceOfLinearExtrusion.hxx>
#include <Geom_SurfaceOfRevolution.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <NCollection_Vector.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Dir.hxx>

#include <algorithm>
#include <cmath>

// If the following is defined the error of algorithm is calculated by static moments
#define IS_MIN_DIM

namespace
{
// Minimal value of interval's range for computation | minimal value of "dim" | ...
static const double EPS_PARAM          = 1.e-12;
static const double EPS_DIM            = 1.e-30;
static const double ERROR_ALGEBR_RATIO = 2.0 / 3.0;

// Maximum of GaussPoints on a subinterval and maximum of subintervals
static const int GPM        = math::GaussPointsMax();
static const int SUBS_POWER = 32;
static const int SM         = SUBS_POWER * GPM + 1;

// Auxiliary inner functions to perform arithmetic operations.
static double Add(const double theA, const double theB)
{
  return theA + theB;
}

static double AddInf(const double theA, const double theB)
{
  if (Precision::IsPositiveInfinite(theA))
  {
    if (Precision::IsNegativeInfinite(theB))
    {
      return 0.0;
    }
    else
    {
      return Precision::Infinite();
    }
  }

  if (Precision::IsPositiveInfinite(theB))
  {
    if (Precision::IsNegativeInfinite(theA))
    {
      return 0.0;
    }
    else
    {
      return Precision::Infinite();
    }
  }

  if (Precision::IsNegativeInfinite(theA))
  {
    if (Precision::IsPositiveInfinite(theB))
    {
      return 0.0;
    }
    else
    {
      return -Precision::Infinite();
    }
  }

  if (Precision::IsNegativeInfinite(theB))
  {
    if (Precision::IsPositiveInfinite(theA))
    {
      return 0.0;
    }
    else
    {
      return -Precision::Infinite();
    }
  }

  return theA + theB;
}

static double Mult(const double theA, const double theB)
{
  return theA * theB;
}

static double MultInf(const double theA, const double theB)
{
  if ((theA == 0.0) || (theB == 0.0))
  { // strictly zerro (without any tolerances)
    return 0.0;
  }

  if (Precision::IsPositiveInfinite(theA))
  {
    if (theB < 0.0)
    {
      return -Precision::Infinite();
    }
    else
    {
      return Precision::Infinite();
    }
  }

  if (Precision::IsPositiveInfinite(theB))
  {
    if (theA < 0.0)
    {
      return -Precision::Infinite();
    }
    else
    {
      return Precision::Infinite();
    }
  }

  if (Precision::IsNegativeInfinite(theA))
  {
    if (theB < 0.0)
    {
      return +Precision::Infinite();
    }
    else
    {
      return -Precision::Infinite();
    }
  }

  if (Precision::IsNegativeInfinite(theB))
  {
    if (theA < 0.0)
    {
      return +Precision::Infinite();
    }
    else
    {
      return -Precision::Infinite();
    }
  }

  return theA * theB;
}

//! The knots of a B-spline -- a face's surface in u or v, or an edge's curve
//! on it -- strictly inside [theA, theB], when a rule of GPM points over the
//! range cannot follow it: more spans of its degree than such a rule
//! integrates. Nothing otherwise, unless theIsAll: every knot inside. Returns
//! the order a span needs, 0 for none.
static int spanKnots(const NCollection_Array1<double>& theKnots,
                     const int                         theDegree,
                     const double                      theA,
                     const double                      theB,
                     NCollection_Vector<double>&       theOut,
                     const bool                        theIsAll = false)
{
  if (!theIsAll && 2 * (theDegree + 1) * (theKnots.Length() - 1) <= GPM)
  {
    return 0;
  }
  const int    aNb   = theOut.Length();
  const double anEps = EPS_PARAM * std::max(1., std::abs(theB - theA));
  for (int i = theKnots.Lower(); i <= theKnots.Upper(); ++i)
  {
    if (theKnots(i) > theA + anEps && theKnots(i) < theB - anEps)
    {
      theOut.Append(theKnots(i));
    }
  }
  return theOut.Length() == aNb ? 0 : 2 * (theDegree + 1);
}

//! Where to cut the range of an offset surface in one parameter -- u (theIsU)
//! from theA to theB, on iso lines across theW1 to theW2, or v likewise --
//! when its normal turns unevenly along it. An offset carries its basis's
//! unit normal, and the integrands with it the square root of the basis's
//! speed: its zeros lie off the real axis by about the inverse of the
//! turning where it peaks, and further off as one moves away from the peak,
//! about as far as one has moved (an ellipse's turning is d / (t^2 + d^2)
//! about a peak of 1 / d). A rule over the whole range cannot follow them
//! (an offset ellipse of 10 x 2, a rule of 18 points over the turn, came 5%
//! off); pieces graded toward the peaks, each about as long as the zeros
//! are far, each with a rule of OFFSET_ORDER points, do. They are spaced by
//! the integral of sqrt(turning x its maximum) -- 1 / sqrt(t^2 + d^2) about
//! such a peak, a few pieces each time the distance to it grows tenfold. A
//! normal turning evenly, a circle's, needs none: there is no zero.
//! Returns the order a piece needs, 0 for none.
static const int OFFSET_ORDER = 20;

static int turningBreaks(const occ::handle<Geom_Surface>& theS,
                         const bool                       theIsU,
                         const double                     theA,
                         const double                     theB,
                         const double                     theW1,
                         const double                     theW2,
                         NCollection_Vector<double>&      theOut)
{
  if (!(theB - theA > Precision::PConfusion()) || Precision::IsInfinite(theB - theA)
      || Precision::IsInfinite(theW2 - theW1))
  {
    return 0;
  }
  // the turning rate on each step, the fastest of three iso lines
  const int aNb = 128;
  double    aRate[aNb];
  std::fill(aRate, aRate + aNb, 0.);
  double aMin = Precision::Infinite(), aMax = 0.;
  for (int w = 1; w <= 3; ++w)
  {
    const double aW = theW1 + (theW2 - theW1) * w / 4.;
    gp_Dir       aPrev;
    bool         hasPrev = false;
    for (int i = 0; i <= aNb; ++i)
    {
      const double aT = theA + (theB - theA) * i / aNb;
      gp_Pnt       aP;
      gp_Vec       aDU, aDV;
      try
      {
        theS->D1(theIsU ? aT : aW, theIsU ? aW : aT, aP, aDU, aDV);
      }
      catch (const Standard_Failure&)
      {
        hasPrev = false;
        continue;
      }
      const gp_Vec aN = aDU ^ aDV;
      if (aN.Magnitude() < gp::Resolution())
      {
        hasPrev = false;
        continue;
      }
      const gp_Dir aDir(aN);
      if (hasPrev)
      {
        const double aR = aPrev.Angle(aDir) * aNb / (theB - theA);
        aRate[i - 1]    = std::max(aRate[i - 1], aR);
        aMin            = std::min(aMin, aR);
        aMax            = std::max(aMax, aR);
      }
      aPrev   = aDir;
      hasPrev = true;
    }
  }
  if (aMax <= 0. || aMax < 1.2 * aMin)
  {
    return 0;
  }
  double aSum[aNb + 1];
  aSum[0] = 0.;
  for (int i = 0; i < aNb; ++i)
  {
    aSum[i + 1] = aSum[i] + std::sqrt(aRate[i] * aMax) * (theB - theA) / aNb;
  }
  const int aNbPieces = std::min(64, static_cast<int>(std::ceil(aSum[aNb])));
  for (int k = 1, i = 0; k < aNbPieces; ++k)
  {
    const double aS = aSum[aNb] * k / aNbPieces;
    while (aSum[i + 1] < aS)
    {
      ++i;
    }
    const double aF = (aS - aSum[i]) / (aSum[i + 1] - aSum[i]);
    theOut.Append(theA + (theB - theA) * (i + aF) / aNb);
  }
  return aNbPieces > 1 ? OFFSET_ORDER : 0;
}

//! The knots a face's surface has in u and in v that a rule of GPM points
//! cannot follow (spanKnots): a B-spline surface's, or those of the B-spline
//! curve an extrusion (u) or a revolution (v) is made of. Returns the order
//! a span needs, 0 for none.
static int surfaceSpanKnots(const TopoDS_Face&          theFace,
                            const double                theU1,
                            const double                theU2,
                            const double                theV1,
                            const double                theV2,
                            NCollection_Vector<double>& theU,
                            NCollection_Vector<double>& theV)
{
  occ::handle<Geom_Surface> aS = BRep_Tool::Surface(theFace);
  occ::handle<Geom_Surface> anOffset;
  // an offset is parametrized as its basis, its spans the basis's
  while (!aS.IsNull()
         && (aS->IsKind(STANDARD_TYPE(Geom_RectangularTrimmedSurface))
             || aS->IsKind(STANDARD_TYPE(Geom_OffsetSurface))))
  {
    if (aS->IsKind(STANDARD_TYPE(Geom_OffsetSurface)) && anOffset.IsNull())
    {
      anOffset = aS;
    }
    aS = aS->IsKind(STANDARD_TYPE(Geom_OffsetSurface))
           ? occ::down_cast<Geom_OffsetSurface>(aS)->BasisSurface()
           : occ::down_cast<Geom_RectangularTrimmedSurface>(aS)->BasisSurface();
  }
  if (aS.IsNull())
  {
    return 0;
  }
  auto aBasis = [](occ::handle<Geom_Curve> theC) {
    while (!theC.IsNull() && theC->IsKind(STANDARD_TYPE(Geom_TrimmedCurve)))
    {
      theC = occ::down_cast<Geom_TrimmedCurve>(theC)->BasisCurve();
    }
    return occ::down_cast<Geom_BSplineCurve>(theC);
  };
  // An offset is cut at every knot of its basis: the offset's derivatives
  // take in the basis's next ones, which jump at a knot -- a degree 2 basis
  // curved on one side of a knot and flat on the other has an offset whose
  // area element doubles there, and a rule across it follows nothing.
  const bool isAll = !anOffset.IsNull();
  int        anOrder = 0;
  if (occ::handle<Geom_BSplineSurface> aBS = occ::down_cast<Geom_BSplineSurface>(aS))
  {
    anOrder = std::max(spanKnots(aBS->UKnots(), aBS->UDegree(), theU1, theU2, theU, isAll),
                       spanKnots(aBS->VKnots(), aBS->VDegree(), theV1, theV2, theV, isAll));
  }
  else if (occ::handle<Geom_SurfaceOfLinearExtrusion> anExt =
             occ::down_cast<Geom_SurfaceOfLinearExtrusion>(aS))
  {
    if (occ::handle<Geom_BSplineCurve> aC = aBasis(anExt->BasisCurve()))
    {
      anOrder = spanKnots(aC->Knots(), aC->Degree(), theU1, theU2, theU, isAll);
    }
  }
  else if (occ::handle<Geom_SurfaceOfRevolution> aRev =
             occ::down_cast<Geom_SurfaceOfRevolution>(aS))
  {
    if (occ::handle<Geom_BSplineCurve> aC = aBasis(aRev->BasisCurve()))
    {
      anOrder = spanKnots(aC->Knots(), aC->Degree(), theV1, theV2, theV, isAll);
    }
  }
  if (isAll && anOrder > 0)
  {
    anOrder = std::max(anOrder, OFFSET_ORDER);
  }
  if (!anOffset.IsNull())
  {
    // the basis's normal, the one the offset carries: the offset's own turns
    // back where the offset folds (beyond a radius of curvature), the
    // integrands smooth through it
    const occ::handle<Geom_Surface> aBasis =
      occ::down_cast<Geom_OffsetSurface>(anOffset)->BasisSurface();
    const int aU = turningBreaks(aBasis, true, theU1, theU2, theV1, theV2, theU);
    const int aV = turningBreaks(aBasis, false, theV1, theV2, theU1, theU2, theV);
    if (aU > 0 || aV > 0)
    {
      anOrder = std::max({anOrder, aU, aV});
      for (NCollection_Vector<double>* aKnots : {&theU, &theV})
      {
        std::sort(aKnots->begin(), aKnots->end());
      }
    }
  }
  return anOrder;
}

//! Where to break the integral along the edge loaded in theSurface, from
//! theL1 to theL2: at its curve's own knots when a rule of GPM points cannot
//! follow it, and where the curve crosses one of the surface's knots kept
//! (surfaceSpanKnots) -- between them the integrand is a span's polynomial.
//! theBreaks gets theL1, the breaks in order, and theL2. Returns the order a
//! span of the edge's curve needs, 0 when it needs none.
static int edgeBreaks(const BRepGProp_Face&             theSurface,
                      const double                      theL1,
                      const double                      theL2,
                      const NCollection_Vector<double>& theU,
                      const NCollection_Vector<double>& theV,
                      NCollection_Vector<double>&       theBreaks)
{
  NCollection_Vector<double> aBreaks;
  int                        anOrder = 0;
  const int                  aNbSubs = theSurface.LIntSubs();
  if (theSurface.IntegrationOrder() > GPM && aNbSubs > 1)
  {
    NCollection_Array1<double> aKnots(1, aNbSubs + 1);
    theSurface.LKnots(aKnots);
    // IntegrationOrder() is twice the degree plus one, times the spans
    anOrder = spanKnots(aKnots,
                        theSurface.IntegrationOrder() / (2 * aNbSubs) - 1,
                        theL1,
                        theL2,
                        aBreaks);
  }
  // where the curve crosses a surface knot: found on a sampling of the
  // curve, then by halving
  const int aNbSamples = 64 + 8 * (theU.Length() + theV.Length() + aNbSubs);
  for (int aDir = 0; aDir < 2; ++aDir)
  {
    const NCollection_Vector<double>& aKnots = aDir == 0 ? theU : theV;
    if (aKnots.IsEmpty())
    {
      continue;
    }
    auto aCoord = [&](const double theL) {
      const gp_Pnt2d aP = theSurface.Value2d(theL);
      return aDir == 0 ? aP.X() : aP.Y();
    };
    double aLPrev = theL1, aPrev = aCoord(theL1);
    for (int i = 1; i <= aNbSamples; ++i)
    {
      const double aL   = theL1 + (theL2 - theL1) * i / aNbSamples;
      const double aCur = aCoord(aL);
      for (int k = 0; k < aKnots.Length(); ++k)
      {
        const double aK = aKnots(k);
        if ((aPrev - aK) * (aCur - aK) >= 0.)
        {
          continue;
        }
        double aA = aLPrev, aB = aL, aFA = aPrev - aK;
        for (int it = 0; it < 60 && aB - aA > EPS_PARAM * std::abs(theL2 - theL1); ++it)
        {
          const double aM = 0.5 * (aA + aB), aFM = aCoord(aM) - aK;
          if (aFA * aFM <= 0.)
          {
            aB = aM;
          }
          else
          {
            aA  = aM;
            aFA = aFM;
          }
        }
        aBreaks.Append(0.5 * (aA + aB));
      }
      aLPrev = aL;
      aPrev  = aCur;
    }
  }
  std::sort(aBreaks.begin(), aBreaks.end());
  theBreaks.Clear();
  theBreaks.Append(theL1);
  const double anEps = EPS_PARAM * std::max(1., std::abs(theL2 - theL1));
  for (int i = 0; i < aBreaks.Length(); ++i)
  {
    if (aBreaks(i) > theBreaks.Last() + anEps && aBreaks(i) < theL2 - anEps)
    {
      theBreaks.Append(aBreaks(i));
    }
  }
  theBreaks.Append(theL2);
  return anOrder;
}

//! The pieces [theA, theB] is cut into at the knots kept (surfaceSpanKnots),
//! from theA to theB: theA, the knots between in that direction, theB.
static void pieces(const double                      theA,
                   const double                      theB,
                   const NCollection_Vector<double>& theKnots,
                   NCollection_Vector<double>&       theOut)
{
  theOut.Clear();
  theOut.Append(theA);
  const double aLo = std::min(theA, theB), aHi = std::max(theA, theB);
  const double anEps = EPS_PARAM * std::max(1., aHi - aLo);
  const int    aNb   = theKnots.Length();
  for (int i = 0; i < aNb; ++i)
  {
    const double aK = theKnots(theA <= theB ? i : aNb - 1 - i);
    if (aK > aLo + anEps && aK < aHi - anEps)
    {
      theOut.Append(aK);
    }
  }
  theOut.Append(theB);
}
} // namespace

//=================================================================================================

BRepGProp_Gauss::Inertia::Inertia()
    : Mass(0.0),
      Ix(0.0),
      Iy(0.0),
      Iz(0.0),
      Ixx(0.0),
      Iyy(0.0),
      Izz(0.0),
      Ixy(0.0),
      Ixz(0.0),
      Iyz(0.0)
{
}

//=================================================================================================

void BRepGProp_Gauss::Inertia::Reset()
{
  memset(reinterpret_cast<void*>(this), 0, sizeof(BRepGProp_Gauss::Inertia));
}

//=================================================================================================

BRepGProp_Gauss::BRepGProp_Gauss(const BRepGProp_GaussType theType)
    : myType(theType)
{
  add  = (::Add);
  mult = (::Mult);
}

//=================================================================================================

int BRepGProp_Gauss::MaxSubs(const int theN, const int theCoeff)
{
  return IntegerLast() / theCoeff < theN ? IntegerLast() : theN * theCoeff + 1;
}

//=================================================================================================

void BRepGProp_Gauss::Init(NCollection_Handle<math_Vector>& theOutVec,
                           const double                     theValue,
                           const int                        theFirst,
                           const int                        theLast)
{
  if (theLast - theFirst == 0)
  {
    theOutVec->Init(theValue);
  }
  else
  {
    for (int i = theFirst; i <= theLast; ++i)
    {
      theOutVec->Value(i) = theValue;
    }
  }
}

//=================================================================================================

void BRepGProp_Gauss::InitMass(const double  theValue,
                               const int     theFirst,
                               const int     theLast,
                               InertiaArray& theArray)
{
  if (theArray.IsNull())
  {
    return;
  }

  int aFirst = theFirst;
  int aLast  = theLast;

  if (theLast - theFirst == 0)
  {
    aFirst = theArray->Lower();
    aLast  = theArray->Upper();
  }

  for (int i = aFirst; i <= aLast; ++i)
  {
    theArray->ChangeValue(i).Mass = theValue;
  }
}

//=================================================================================================

int BRepGProp_Gauss::FillIntervalBounds(const double                      theA,
                                        const double                      theB,
                                        const NCollection_Array1<double>& theKnots,
                                        const int                         theNumSubs,
                                        InertiaArray&                     theInerts,
                                        NCollection_Handle<math_Vector>&  theParam1,
                                        NCollection_Handle<math_Vector>&  theParam2,
                                        NCollection_Handle<math_Vector>&  theError,
                                        NCollection_Handle<math_Vector>&  theCommonError)
{
  const int aSize = std::max(theKnots.Upper(), MaxSubs(theKnots.Upper() - 1, theNumSubs));

  if (aSize - 1 > theParam1->Upper())
  {
    theInerts = new NCollection_Array1<Inertia>(1, aSize);
    theParam1 = new math_Vector(1, aSize);
    theParam2 = new math_Vector(1, aSize);
    theError  = new math_Vector(1, aSize, 0.0);

    if (!theCommonError.IsNull())
    {
      theCommonError = new math_Vector(1, aSize, 0.0);
    }
  }

  int j = 1, k = 1;
  theParam1->Value(j++) = theA;

  const int aLength = theKnots.Upper();
  for (int i = 1; i <= aLength; ++i)
  {
    const double kn = theKnots(i);
    if (theA < kn)
    {
      if (kn < theB)
      {
        theParam1->Value(j++) = kn;
        theParam2->Value(k++) = kn;
      }
      else
      {
        break;
      }
    }
  }

  theParam2->Value(k) = theB;
  return k;
}

//=================================================================================================

void BRepGProp_Gauss::computeVInertiaOfElementaryPart(const gp_Pnt&             thePoint,
                                                      const gp_Vec&             theNormal,
                                                      const gp_Pnt&             theLocation,
                                                      const double              theWeight,
                                                      const double              theCoeff[],
                                                      const bool                theIsByPoint,
                                                      BRepGProp_Gauss::Inertia& theOutInertia)
{
  double x = thePoint.X() - theLocation.X();
  double y = thePoint.Y() - theLocation.Y();
  double z = thePoint.Z() - theLocation.Z();

  const double xn = theNormal.X() * theWeight;
  const double yn = theNormal.Y() * theWeight;
  const double zn = theNormal.Z() * theWeight;

  if (theIsByPoint)
  {
    /////////////////////                        ///////////////////////
    //    OFV code     //                        //    Initial code   //
    /////////////////////                        ///////////////////////
    // modified by APO

    double dv = x * xn + y * yn + z * zn; // xyz  = x * y * z;
    theOutInertia.Mass += dv / 3.0;       // Ixyi += zn * xyz;
    theOutInertia.Ix += 0.25 * x * dv;    // Iyzi += xn * xyz;
    theOutInertia.Iy += 0.25 * y * dv;    // Ixzi += yn * xyz;
    theOutInertia.Iz += 0.25 * z * dv;    // xi = x * x * x * xn / 3.0;
    x -= theCoeff[0];                     // yi = y * y * y * yn / 3.0;
    y -= theCoeff[1];                     // zi = z * z * z * zn / 3.0;
    z -= theCoeff[2];                     // Ixxi += (yi + zi);
    dv *= 0.2;                            // Iyyi += (xi + zi);
    theOutInertia.Ixy -= x * y * dv;      // Izzi += (xi + yi);
    theOutInertia.Iyz -= y * z * dv;      // x -= Coeff[0];
    theOutInertia.Ixz -= x * z * dv;      // y -= Coeff[1];
    x *= x;                               // z -= Coeff[2];
    y *= y;                               // dv = x * xn + y * yn + z * zn;
    z *= z;                               // dvi +=  dv;
    theOutInertia.Ixx += (y + z) * dv;    // Ixi += x * dv;
    theOutInertia.Iyy += (x + z) * dv;    // Iyi += y * dv;
    theOutInertia.Izz += (x + y) * dv;    // Izi += z * dv;
  }
  else
  { // By plane
    const double s = xn * theCoeff[0] + yn * theCoeff[1] + zn * theCoeff[2];

    double d1 = theCoeff[0] * x + theCoeff[1] * y + theCoeff[2] * z - theCoeff[3];
    double d2 = d1 * d1;
    double d3 = d1 * d2 / 3.0;
    double dv = s * d1;

    theOutInertia.Mass += dv;
    theOutInertia.Ix += (x - (theCoeff[0] * d1 * 0.5)) * dv;
    theOutInertia.Iy += (y - (theCoeff[1] * d1 * 0.5)) * dv;
    theOutInertia.Iz += (z - (theCoeff[2] * d1 * 0.5)) * dv;

    const double px = x - theCoeff[0] * d1;
    const double py = y - theCoeff[1] * d1;
    const double pz = z - theCoeff[2] * d1;

    x = px * px * d1 + px * theCoeff[0] * d2 + theCoeff[0] * theCoeff[0] * d3;
    y = py * py * d1 + py * theCoeff[1] * d2 + theCoeff[1] * theCoeff[1] * d3;
    z = pz * pz * d1 + pz * theCoeff[2] * d2 + theCoeff[2] * theCoeff[2] * d3;

    theOutInertia.Ixx += (y + z) * s;
    theOutInertia.Iyy += (x + z) * s;
    theOutInertia.Izz += (x + y) * s;

    d2 *= 0.5;
    x = (py * pz * d1) + (py * theCoeff[2] * d2) + (pz * theCoeff[1] * d2)
        + (theCoeff[1] * theCoeff[2] * d3);
    y = (px * pz * d1) + (pz * theCoeff[0] * d2) + (px * theCoeff[2] * d2)
        + (theCoeff[0] * theCoeff[2] * d3);
    z = (px * py * d1) + (px * theCoeff[1] * d2) + (py * theCoeff[0] * d2)
        + (theCoeff[0] * theCoeff[1] * d3);

    theOutInertia.Ixy -= z * s;
    theOutInertia.Iyz -= x * s;
    theOutInertia.Ixz -= y * s;
  }
}

//=================================================================================================

void BRepGProp_Gauss::computeSInertiaOfElementaryPart(const gp_Pnt&             thePoint,
                                                      const gp_Vec&             theNormal,
                                                      const gp_Pnt&             theLocation,
                                                      const double              theWeight,
                                                      BRepGProp_Gauss::Inertia& theOutInertia)
{
  // ds - Jacobien (x, y, z) -> (u, v) = ||n||
  const double ds = mult(theNormal.Magnitude(), theWeight);
  const double x  = add(thePoint.X(), -theLocation.X());
  const double y  = add(thePoint.Y(), -theLocation.Y());
  const double z  = add(thePoint.Z(), -theLocation.Z());

  theOutInertia.Mass = add(theOutInertia.Mass, ds);

  const double XdS = mult(x, ds);
  const double YdS = mult(y, ds);
  const double ZdS = mult(z, ds);

  theOutInertia.Ix  = add(theOutInertia.Ix, XdS);
  theOutInertia.Iy  = add(theOutInertia.Iy, YdS);
  theOutInertia.Iz  = add(theOutInertia.Iz, ZdS);
  theOutInertia.Ixy = add(theOutInertia.Ixy, mult(x, YdS));
  theOutInertia.Iyz = add(theOutInertia.Iyz, mult(y, ZdS));
  theOutInertia.Ixz = add(theOutInertia.Ixz, mult(x, ZdS));

  const double XXdS = mult(x, XdS);
  const double YYdS = mult(y, YdS);
  const double ZZdS = mult(z, ZdS);

  theOutInertia.Ixx = add(theOutInertia.Ixx, add(YYdS, ZZdS));
  theOutInertia.Iyy = add(theOutInertia.Iyy, add(XXdS, ZZdS));
  theOutInertia.Izz = add(theOutInertia.Izz, add(XXdS, YYdS));
}

//=================================================================================================

void BRepGProp_Gauss::checkBounds(const double theU1,
                                  const double theU2,
                                  const double theV1,
                                  const double theV2)
{
  if (Precision::IsInfinite(theU1) || Precision::IsInfinite(theU2) || Precision::IsInfinite(theV1)
      || Precision::IsInfinite(theV2))
  {
    add  = (::AddInf);
    mult = (::MultInf);
  }
}

//=================================================================================================

void BRepGProp_Gauss::addAndRestoreInertia(const BRepGProp_Gauss::Inertia& theInInertia,
                                           BRepGProp_Gauss::Inertia&       theOutInertia)
{
  theOutInertia.Mass = add(theOutInertia.Mass, theInInertia.Mass);
  theOutInertia.Ix   = add(theOutInertia.Ix, theInInertia.Ix);
  theOutInertia.Iy   = add(theOutInertia.Iy, theInInertia.Iy);
  theOutInertia.Iz   = add(theOutInertia.Iz, theInInertia.Iz);
  theOutInertia.Ixx  = add(theOutInertia.Ixx, theInInertia.Ixx);
  theOutInertia.Iyy  = add(theOutInertia.Iyy, theInInertia.Iyy);
  theOutInertia.Izz  = add(theOutInertia.Izz, theInInertia.Izz);
  theOutInertia.Ixy  = add(theOutInertia.Ixy, theInInertia.Ixy);
  theOutInertia.Ixz  = add(theOutInertia.Ixz, theInInertia.Ixz);
  theOutInertia.Iyz  = add(theOutInertia.Iyz, theInInertia.Iyz);
}

//=================================================================================================

void BRepGProp_Gauss::multAndRestoreInertia(const double              theValue,
                                            BRepGProp_Gauss::Inertia& theInOutInertia)
{
  theInOutInertia.Mass = mult(theInOutInertia.Mass, theValue);
  theInOutInertia.Ix   = mult(theInOutInertia.Ix, theValue);
  theInOutInertia.Iy   = mult(theInOutInertia.Iy, theValue);
  theInOutInertia.Iz   = mult(theInOutInertia.Iz, theValue);
  theInOutInertia.Ixx  = mult(theInOutInertia.Ixx, theValue);
  theInOutInertia.Iyy  = mult(theInOutInertia.Iyy, theValue);
  theInOutInertia.Izz  = mult(theInOutInertia.Izz, theValue);
  theInOutInertia.Ixy  = mult(theInOutInertia.Ixy, theValue);
  theInOutInertia.Ixz  = mult(theInOutInertia.Ixz, theValue);
  theInOutInertia.Iyz  = mult(theInOutInertia.Iyz, theValue);
}

//=================================================================================================

void BRepGProp_Gauss::convert(const BRepGProp_Gauss::Inertia& theInertia,
                              gp_Pnt&                         theOutGravityCenter,
                              gp_Mat&                         theOutMatrixOfInertia,
                              double&                         theOutMass)
{
  if (std::abs(theInertia.Mass) >= EPS_DIM)
  {
    const double anInvMass = 1.0 / theInertia.Mass;
    theOutGravityCenter.SetX(theInertia.Ix * anInvMass);
    theOutGravityCenter.SetY(theInertia.Iy * anInvMass);
    theOutGravityCenter.SetZ(theInertia.Iz * anInvMass);

    theOutMass = theInertia.Mass;
  }
  else
  {
    theOutMass = 0.0;
    theOutGravityCenter.SetCoord(0.0, 0.0, 0.0);
  }

  theOutMatrixOfInertia = gp_Mat(gp_XYZ(theInertia.Ixx, -theInertia.Ixy, -theInertia.Ixz),
                                 gp_XYZ(-theInertia.Ixy, theInertia.Iyy, -theInertia.Iyz),
                                 gp_XYZ(-theInertia.Ixz, -theInertia.Iyz, theInertia.Izz));
}

//=================================================================================================

void BRepGProp_Gauss::convert(const BRepGProp_Gauss::Inertia& theInertia,
                              const double                    theCoeff[],
                              const bool                      theIsByPoint,
                              gp_Pnt&                         theOutGravityCenter,
                              gp_Mat&                         theOutMatrixOfInertia,
                              double&                         theOutMass)
{
  convert(theInertia, theOutGravityCenter, theOutMatrixOfInertia, theOutMass);
  if (std::abs(theInertia.Mass) >= EPS_DIM && theIsByPoint)
  {
    const double anInvMass = 1.0 / theInertia.Mass;
    if (theIsByPoint)
    {
      theOutGravityCenter.SetX(theCoeff[0] + theInertia.Ix * anInvMass);
      theOutGravityCenter.SetY(theCoeff[1] + theInertia.Iy * anInvMass);
      theOutGravityCenter.SetZ(theCoeff[2] + theInertia.Iz * anInvMass);
    }
    else
    {
      theOutGravityCenter.SetX(theInertia.Ix * anInvMass);
      theOutGravityCenter.SetY(theInertia.Iy * anInvMass);
      theOutGravityCenter.SetZ(theInertia.Iz * anInvMass);
    }

    theOutMass = theInertia.Mass;
  }
  else
  {
    theOutMass = 0.0;
    theOutGravityCenter.SetCoord(0.0, 0.0, 0.0);
  }

  theOutMatrixOfInertia = gp_Mat(gp_XYZ(theInertia.Ixx, theInertia.Ixy, theInertia.Ixz),
                                 gp_XYZ(theInertia.Ixy, theInertia.Iyy, theInertia.Iyz),
                                 gp_XYZ(theInertia.Ixz, theInertia.Iyz, theInertia.Izz));
}

//=================================================================================================

double BRepGProp_Gauss::Compute(BRepGProp_Face&   theSurface,
                                BRepGProp_Domain& theDomain,
                                const gp_Pnt&     theLocation,
                                const double      theEps,
                                const double      theCoeff[],
                                const bool        theIsByPoint,
                                double&           theOutMass,
                                gp_Pnt&           theOutGravityCenter,
                                gp_Mat&           theOutInertia)
{
  const bool isErrorCalculation  = 0.0 > theEps || theEps < 0.001;
  const bool isVerifyComputation = 0.0 < theEps && theEps < 0.001;

  double anEpsilon = std::abs(theEps);

  BRepGProp_Gauss::Inertia anInertia;
  InertiaArray             anInertiaL = new NCollection_Array1<Inertia>(1, SM);
  InertiaArray             anInertiaU = new NCollection_Array1<Inertia>(1, SM);

  // Prepare Gauss points and weights
  NCollection_Handle<math_Vector> LGaussP[2];
  NCollection_Handle<math_Vector> LGaussW[2];
  NCollection_Handle<math_Vector> UGaussP[2];
  NCollection_Handle<math_Vector> UGaussW[2];

  const int aNbGaussPoint = RealToInt(std::ceil(ERROR_ALGEBR_RATIO * GPM));

  LGaussP[0] = new math_Vector(1, GPM);
  LGaussP[1] = new math_Vector(1, aNbGaussPoint);
  LGaussW[0] = new math_Vector(1, GPM);
  LGaussW[1] = new math_Vector(1, aNbGaussPoint);

  UGaussP[0] = new math_Vector(1, GPM);
  UGaussP[1] = new math_Vector(1, aNbGaussPoint);
  UGaussW[0] = new math_Vector(1, GPM);
  UGaussW[1] = new math_Vector(1, aNbGaussPoint);

  NCollection_Handle<math_Vector> L1 = new math_Vector(1, SM);
  NCollection_Handle<math_Vector> L2 = new math_Vector(1, SM);
  NCollection_Handle<math_Vector> U1 = new math_Vector(1, SM);
  NCollection_Handle<math_Vector> U2 = new math_Vector(1, SM);

  NCollection_Handle<math_Vector> ErrL  = new math_Vector(1, SM, 0.0);
  NCollection_Handle<math_Vector> ErrU  = new math_Vector(1, SM, 0.0);
  NCollection_Handle<math_Vector> ErrUL = new math_Vector(1, SM, 0.0);

  // Face parametrization in U and V direction
  double BV1, BV2, BU1, BU2;
  theSurface.Bounds(BU1, BU2, BV1, BV2);
  checkBounds(BU1, BU2, BV1, BV2);

  //
  const int          NumSubs = SUBS_POWER;
  const TopoDS_Face& aF      = theSurface.GetFace();
  // clang-format off
  const bool isNaturalRestriction = (aF.NbChildren () == 0); //theSurface.NaturalRestriction();
  // clang-format on

  double CIx, CIy, CIz, CIxy, CIxz, CIyz;
  double CDim[2], CIxx[2], CIyy[2], CIzz[2];

  // Boundary curve parametrization
  double u1 = BU1, u2, l1, l2, lm, lr, l, v;

  // On the boundary curve u-v
  gp_Pnt2d Puv;
  gp_Vec2d Vuv;
  double   Dul; // Dul = Du / Dl

  int iLS, iLSubEnd, iGL, iGLEnd, NbLGaussP[2], LRange[2], iL, kL, kLEnd, IL, JL;
  int i, iUSubEnd, NbUGaussP[2], URange[2], kU, kUEnd, IU, JU;
  int UMaxSubs, LMaxSubs;

  double ErrorU, ErrorL, ErrorLMax = 0.0, Eps = 0.0, EpsL = 0.0, EpsU = 0.0;
  iGLEnd = isErrorCalculation ? 2 : 1;

  NbUGaussP[0] = theSurface.SIntOrder(anEpsilon);
  NbUGaussP[1] = RealToInt(std::ceil(ERROR_ALGEBR_RATIO * NbUGaussP[0]));

  math::GaussPoints(NbUGaussP[0], *UGaussP[0]);
  math::GaussWeights(NbUGaussP[0], *UGaussW[0]);
  math::GaussPoints(NbUGaussP[1], *UGaussP[1]);
  math::GaussWeights(NbUGaussP[1], *UGaussW[1]);

  const int                  aNbUSubs = theSurface.SUIntSubs();
  NCollection_Array1<double> UKnots(1, aNbUSubs + 1);
  theSurface.UKnots(UKnots);

  while (isNaturalRestriction || theDomain.More())
  {
    if (isNaturalRestriction)
    {
      NbLGaussP[0] = std::min(2 * NbUGaussP[0], math::GaussPointsMax());
    }
    else
    {
      if (!theSurface.Load(theDomain.Value()))
      {
        return Precision::Infinite();
      }
      NbLGaussP[0] = theSurface.LIntOrder(anEpsilon);
    }

    NbLGaussP[1] = RealToInt(std::ceil(ERROR_ALGEBR_RATIO * NbLGaussP[0]));

    math::GaussPoints(NbLGaussP[0], *LGaussP[0]);
    math::GaussWeights(NbLGaussP[0], *LGaussW[0]);
    math::GaussPoints(NbLGaussP[1], *LGaussP[1]);
    math::GaussWeights(NbLGaussP[1], *LGaussW[1]);

    const int aNbLSubs = isNaturalRestriction ? theSurface.SVIntSubs() : theSurface.LIntSubs();
    NCollection_Array1<double> LKnots(1, aNbLSubs + 1);

    if (isNaturalRestriction)
    {
      theSurface.VKnots(LKnots);
      l1 = BV1;
      l2 = BV2;
    }
    else
    {
      theSurface.LKnots(LKnots);
      l1 = theSurface.FirstParameter();
      l2 = theSurface.LastParameter();
    }
    ErrorL = 0.0;
    kLEnd  = 1;
    JL     = 0;

    if (std::abs(l2 - l1) > EPS_PARAM)
    {
      iLSubEnd = FillIntervalBounds(l1, l2, LKnots, NumSubs, anInertiaL, L1, L2, ErrL, ErrUL);
      LMaxSubs = BRepGProp_Gauss::MaxSubs(iLSubEnd);

      if (LMaxSubs > SM)
      {
        LMaxSubs = SM;
      }

      BRepGProp_Gauss::InitMass(0.0, 1, LMaxSubs, anInertiaL);
      BRepGProp_Gauss::Init(ErrL, 0.0, 1, LMaxSubs);
      BRepGProp_Gauss::Init(ErrUL, 0.0, 1, LMaxSubs);

      do // while: L
      {
        if (++JL > iLSubEnd)
        {
          LRange[0] = IL = ErrL->Max();
          LRange[1]      = JL;
          L1->Value(JL)  = (L1->Value(IL) + L2->Value(IL)) * 0.5;
          L2->Value(JL)  = L2->Value(IL);
          L2->Value(IL)  = L1->Value(JL);
        }
        else
        {
          LRange[0] = IL = JL;
        }

        if (JL == LMaxSubs || std::abs(L2->Value(JL) - L1->Value(JL)) < EPS_PARAM)
        {
          if (kLEnd == 1)
          {
            anInertiaL->ChangeValue(JL).Reset();
            ErrL->Value(JL) = 0.0;
          }
          else
          {
            --JL;
            EpsL = ErrorL;
            Eps  = EpsL / 0.9;
            break;
          }
        }
        else
        {
          for (kL = 0; kL < kLEnd; kL++)
          {
            iLS = LRange[kL];
            lm  = 0.5 * (L2->Value(iLS) + L1->Value(iLS));
            lr  = 0.5 * (L2->Value(iLS) - L1->Value(iLS));

            CIx = CIy = CIz = CIxy = CIxz = CIyz = 0.0;

            for (iGL = 0; iGL < iGLEnd; ++iGL)
            {
              CDim[iGL] = CIxx[iGL] = CIyy[iGL] = CIzz[iGL] = 0.0;

              for (iL = 1; iL <= NbLGaussP[iGL]; iL++)
              {
                l = lm + lr * LGaussP[iGL]->Value(iL);
                if (isNaturalRestriction)
                {
                  v   = l;
                  u2  = BU2;
                  Dul = LGaussW[iGL]->Value(iL);
                }
                else
                {
                  theSurface.D12d(l, Puv, Vuv);
                  Dul = Vuv.Y() * LGaussW[iGL]->Value(iL); // Dul = Du / Dl

                  if (std::abs(Dul) < EPS_PARAM)
                  {
                    continue;
                  }

                  v  = Puv.Y();
                  u2 = Puv.X();

                  // Check on cause out off bounds of value current parameter
                  if (v < BV1)
                  {
                    v = BV1;
                  }
                  else if (v > BV2)
                  {
                    v = BV2;
                  }

                  if (u2 < BU1)
                  {
                    u2 = BU1;
                  }
                  else if (u2 > BU2)
                  {
                    u2 = BU2;
                  }
                }

                ErrUL->Value(iLS) = 0.0;
                kUEnd             = 1;
                JU                = 0;

                if (std::abs(u2 - u1) < EPS_PARAM)
                {
                  continue;
                }

                NCollection_Handle<math_Vector> aDummy;
                iUSubEnd =
                  FillIntervalBounds(u1, u2, UKnots, NumSubs, anInertiaU, U1, U2, ErrU, aDummy);
                UMaxSubs = BRepGProp_Gauss::MaxSubs(iUSubEnd);

                if (UMaxSubs > SM)
                {
                  UMaxSubs = SM;
                }

                BRepGProp_Gauss::InitMass(0.0, 1, UMaxSubs, anInertiaU);
                BRepGProp_Gauss::Init(ErrU, 0.0, 1, UMaxSubs);
                ErrorU = 0.0;

                do
                { // while: U
                  if (++JU > iUSubEnd)
                  {
                    URange[0] = IU = ErrU->Max();
                    URange[1]      = JU;

                    U1->Value(JU) = (U1->Value(IU) + U2->Value(IU)) * 0.5;
                    U2->Value(JU) = U2->Value(IU);
                    U2->Value(IU) = U1->Value(JU);
                  }
                  else
                  {
                    URange[0] = IU = JU;
                  }

                  if (JU == UMaxSubs || std::abs(U2->Value(JU) - U1->Value(JU)) < EPS_PARAM)
                  {
                    if (kUEnd == 1)
                    {
                      ErrU->Value(JU) = 0.0;
                      anInertiaU->ChangeValue(JU).Reset();
                    }
                    else
                    {
                      --JU;
                      EpsU = ErrorU;
                      Eps  = 10. * EpsU * std::abs((u2 - u1) * Dul);
                      EpsL = 0.9 * Eps;
                      break;
                    }
                  }
                  else
                  {
                    gp_Pnt aPoint;
                    gp_Vec aNormal;

                    for (kU = 0; kU < kUEnd; ++kU)
                    {
                      BRepGProp_Gauss::Inertia aLocal[2];

                      int       iUS     = URange[kU];
                      const int aLength = iGLEnd - iGL;

                      const double um = 0.5 * (U2->Value(iUS) + U1->Value(iUS));
                      const double ur = 0.5 * (U2->Value(iUS) - U1->Value(iUS));

                      for (int iGU = 0; iGU < aLength; ++iGU)
                      {
                        for (int iU = 1; iU <= NbUGaussP[iGU]; ++iU)
                        {
                          double       w = UGaussW[iGU]->Value(iU);
                          const double u = um + ur * UGaussP[iGU]->Value(iU);

                          theSurface.Normal(u, v, aPoint, aNormal);

                          if (myType == Vinert)
                          {
                            computeVInertiaOfElementaryPart(aPoint,
                                                            aNormal,
                                                            theLocation,
                                                            w,
                                                            theCoeff,
                                                            theIsByPoint,
                                                            aLocal[iGU]);
                          }
                          else
                          {
                            if (iGU > 0)
                            {
                              aLocal[iGU].Mass += (w * aNormal.Magnitude());
                            }
                            else
                            {
                              computeSInertiaOfElementaryPart(aPoint,
                                                              aNormal,
                                                              theLocation,
                                                              w,
                                                              aLocal[iGU]);
                            }
                          }
                        }
                      }

                      BRepGProp_Gauss::Inertia& anUI = anInertiaU->ChangeValue(iUS);

                      anUI.Mass = mult(aLocal[0].Mass, ur);

                      if (myType == Vinert)
                      {
                        anUI.Ixx = mult(aLocal[0].Ixx, ur);
                        anUI.Iyy = mult(aLocal[0].Iyy, ur);
                        anUI.Izz = mult(aLocal[0].Izz, ur);
                      }

                      if (iGL > 0)
                      {
                        continue;
                      }

                      double aDMass = std::abs(aLocal[1].Mass - aLocal[0].Mass);

                      if (myType == Vinert)
                      {
                        aLocal[1].Ixx = std::abs(aLocal[1].Ixx - aLocal[0].Ixx);
                        aLocal[1].Iyy = std::abs(aLocal[1].Iyy - aLocal[0].Iyy);
                        aLocal[1].Izz = std::abs(aLocal[1].Izz - aLocal[0].Izz);

                        anUI.Ix = mult(aLocal[0].Ix, ur);
                        anUI.Iy = mult(aLocal[0].Iy, ur);
                        anUI.Iz = mult(aLocal[0].Iz, ur);

                        anUI.Ixy = mult(aLocal[0].Ixy, ur);
                        anUI.Ixz = mult(aLocal[0].Ixz, ur);
                        anUI.Iyz = mult(aLocal[0].Iyz, ur);

#ifndef IS_MIN_DIM
                        aDMass = aLocal[1].Ixx + aLocal[1].Iyy + aLocal[1].Izz;
#endif

                        ErrU->Value(iUS) = mult(aDMass, ur);
                      }
                      else
                      {
                        anUI.Ix  = mult(aLocal[0].Ix, ur);
                        anUI.Iy  = mult(aLocal[0].Iy, ur);
                        anUI.Iz  = mult(aLocal[0].Iz, ur);
                        anUI.Ixx = mult(aLocal[0].Ixx, ur);
                        anUI.Iyy = mult(aLocal[0].Iyy, ur);
                        anUI.Izz = mult(aLocal[0].Izz, ur);
                        anUI.Ixy = mult(aLocal[0].Ixy, ur);
                        anUI.Ixz = mult(aLocal[0].Ixz, ur);
                        anUI.Iyz = mult(aLocal[0].Iyz, ur);

                        ErrU->Value(iUS) = mult(aDMass, ur);
                      }
                    }
                  }

                  if (JU == iUSubEnd)
                  {
                    kUEnd  = 2;
                    ErrorU = ErrU->Value(ErrU->Max());
                  }
                } while ((ErrorU - EpsU > 0.0 && EpsU != 0.0) || kUEnd == 1);

                for (i = 1; i <= JU; ++i)
                {
                  const BRepGProp_Gauss::Inertia& anIU = anInertiaU->Value(i);

                  CDim[iGL] = add(CDim[iGL], mult(anIU.Mass, Dul));
                  CIxx[iGL] = add(CIxx[iGL], mult(anIU.Ixx, Dul));
                  CIyy[iGL] = add(CIyy[iGL], mult(anIU.Iyy, Dul));
                  CIzz[iGL] = add(CIzz[iGL], mult(anIU.Izz, Dul));
                }

                if (iGL > 0)
                {
                  continue;
                }

                ErrUL->Value(iLS) = ErrorU * std::abs((u2 - u1) * Dul);

                for (i = 1; i <= JU; ++i)
                {
                  const BRepGProp_Gauss::Inertia& anIU = anInertiaU->Value(i);

                  CIx = add(CIx, mult(anIU.Ix, Dul));
                  CIy = add(CIy, mult(anIU.Iy, Dul));
                  CIz = add(CIz, mult(anIU.Iz, Dul));

                  CIxy = add(CIxy, mult(anIU.Ixy, Dul));
                  CIxz = add(CIxz, mult(anIU.Ixz, Dul));
                  CIyz = add(CIyz, mult(anIU.Iyz, Dul));
                }
              } // for: iL
            } // for: iGL

            BRepGProp_Gauss::Inertia& aLI = anInertiaL->ChangeValue(iLS);

            aLI.Mass = mult(CDim[0], lr);
            aLI.Ixx  = mult(CIxx[0], lr);
            aLI.Iyy  = mult(CIyy[0], lr);
            aLI.Izz  = mult(CIzz[0], lr);

            if (iGLEnd == 2)
            {
              double aSubDim = std::abs(CDim[1] - CDim[0]);

              if (myType == Vinert)
              {
                ErrorU = ErrUL->Value(iLS);

                CIxx[1] = std::abs(CIxx[1] - CIxx[0]);
                CIyy[1] = std::abs(CIyy[1] - CIyy[0]);
                CIzz[1] = std::abs(CIzz[1] - CIzz[0]);

#ifndef IS_MIN_DIM
                aSubDim = CIxx[1] + CIyy[1] + CIzz[1];
#endif

                ErrL->Value(iLS) = add(mult(aSubDim, lr), ErrorU);
              }
              else
              {
                ErrL->Value(iLS) = add(mult(aSubDim, lr), ErrUL->Value(iLS));
              }
            }

            aLI.Ix = mult(CIx, lr);
            aLI.Iy = mult(CIy, lr);
            aLI.Iz = mult(CIz, lr);

            aLI.Ixy = mult(CIxy, lr);
            aLI.Ixz = mult(CIxz, lr);
            aLI.Iyz = mult(CIyz, lr);
          } // for: (kL)iLS
        }

        // Calculate/correct epsilon of computation by current value of dim
        // That is need for not spend time for
        if (JL == iLSubEnd)
        {
          kLEnd = 2;

          double DDim = 0.0;
          for (i = 1; i <= JL; ++i)
          {
            DDim += anInertiaL->Value(i).Mass;
          }

#ifndef IS_MIN_DIM
          {
            if (myType == Vinert)
            {
              double DIxx = 0.0, DIyy = 0.0, DIzz = 0.0;
              for (i = 1; i <= JL; ++i)
              {
                const BRepGProp_Gauss::Inertia& aLocalL = anInertiaL->Value(i);

                DIxx += aLocalL.Ixx;
                DIyy += aLocalL.Iyy;
                DIzz += aLocalL.Izz;
              }

              DDim = std::abs(DIxx) + std::abs(DIyy) + std::abs(DIzz);
            }
          }
#endif

          DDim = std::abs(DDim * anEpsilon);

          if (DDim > Eps)
          {
            Eps  = DDim;
            EpsL = 0.9 * Eps;
          }
        }
        if (kLEnd == 2)
        {
          ErrorL = ErrL->Value(ErrL->Max());
        }
      } while ((ErrorL - EpsL > 0.0 && isVerifyComputation) || kLEnd == 1);

      for (i = 1; i <= JL; i++)
      {
        addAndRestoreInertia(anInertiaL->Value(i), anInertia);
      }

      ErrorLMax = std::max(ErrorLMax, ErrorL);
    }

    if (isNaturalRestriction)
    {
      break;
    }

    theDomain.Next();
  }

  if (myType == Vinert)
  {
    convert(anInertia, theCoeff, theIsByPoint, theOutGravityCenter, theOutInertia, theOutMass);
  }
  else
  {
    convert(anInertia, theOutGravityCenter, theOutInertia, theOutMass);
  }

  if (iGLEnd == 2)
  {
    if (theOutMass != 0.0)
    {
      Eps = ErrorLMax / std::abs(theOutMass);

#ifndef IS_MIN_DIM
      {
        if (myType == Vinert)
          Eps = ErrorLMax
                / (std::abs(anInertia.Ixx) + std::abs(anInertia.Iyy) + std::abs(anInertia.Izz));
      }
#endif
    }
    else
    {
      Eps = 0.0;
    }
  }
  else
  {
    Eps = anEpsilon;
  }

  return Eps;
}

//=================================================================================================

double BRepGProp_Gauss::Compute(BRepGProp_Face&   theSurface,
                                BRepGProp_Domain& theDomain,
                                const gp_Pnt&     theLocation,
                                const double      theEps,
                                double&           theOutMass,
                                gp_Pnt&           theOutGravityCenter,
                                gp_Mat&           theOutInertia)
{
  Standard_ASSERT_RAISE(myType == Sinert, "BRepGProp_Gauss: Incorrect type");

  return Compute(theSurface,
                 theDomain,
                 theLocation,
                 theEps,
                 nullptr,
                 true,
                 theOutMass,
                 theOutGravityCenter,
                 theOutInertia);
}

//=================================================================================================

void BRepGProp_Gauss::Compute(BRepGProp_Face&   theSurface,
                              BRepGProp_Domain& theDomain,
                              const gp_Pnt&     theLocation,
                              double&           theOutMass,
                              gp_Pnt&           theOutGravityCenter,
                              gp_Mat&           theOutInertia)
{
  Standard_ASSERT_RAISE(myType == Sinert, "BRepGProp_Gauss: Incorrect type");

  double u1, u2, v1, v2;
  theSurface.Bounds(u1, u2, v1, v2);
  checkBounds(u1, u2, v1, v2);

  const int NbUGaussgp_Pnts = std::min(theSurface.UIntegrationOrder(), math::GaussPointsMax());

  const int NbVGaussgp_Pnts = std::min(theSurface.VIntegrationOrder(), math::GaussPointsMax());

  const int NbGaussgp_Pnts = std::max(NbUGaussgp_Pnts, NbVGaussgp_Pnts);

  // Number of Gauss points for the integration on the face
  math_Vector GaussSPV(1, NbGaussgp_Pnts);
  math_Vector GaussSWV(1, NbGaussgp_Pnts);
  math::GaussPoints(NbGaussgp_Pnts, GaussSPV);
  math::GaussWeights(NbGaussgp_Pnts, GaussSWV);

  // A surface or an edge of more spans than a rule of GaussPointsMax()
  // points follows is integrated span by span, each with a rule of its own
  // (#617: an extrusion of a degree 8 B-spline of 29 spans, its knots
  // bunched, came 30% off). Anything else is integrated as it was.
  NCollection_Vector<double> aUKnots, aVKnots, aBreaks, aPieces;
  const int aSurfOrder = surfaceSpanKnots(theSurface.GetFace(), u1, u2, v1, v2, aUKnots, aVKnots);

  // the integral along the edge from theA to theB, of the integral across
  // the face from u1, in pieces at theKnots
  BRepGProp_Gauss::Inertia anInertia;
  auto anAlong = [&](const double                      theA,
                     const double                      theB,
                     const math_Vector&                theCP,
                     const math_Vector&                theCW,
                     const math_Vector&                theSP,
                     const math_Vector&                theSW,
                     const NCollection_Vector<double>& theKnots) {
    const double lm = 0.5 * (theB + theA);
    const double lr = 0.5 * (theB - theA);

    BRepGProp_Gauss::Inertia aCInertia;
    for (int i = theCP.Lower(); i <= theCP.Upper(); ++i)
    {
      const double l = lm + lr * theCP(i);

      gp_Pnt2d Puv;
      gp_Vec2d Vuv;
      theSurface.D12d(l, Puv, Vuv);

      const double v   = Puv.Y();
      const double Dul = Vuv.Y() * theCW(i);

      pieces(u1, Puv.X(), theKnots, aPieces);
      for (int p = 1; p < aPieces.Length(); ++p)
      {
        const double um = 0.5 * (aPieces(p) + aPieces(p - 1));
        const double ur = 0.5 * (aPieces(p) - aPieces(p - 1));

        BRepGProp_Gauss::Inertia aLocalInertia;
        for (int j = theSP.Lower(); j <= theSP.Upper(); ++j)
        {
          const double u       = add(um, mult(ur, theSP(j)));
          const double aWeight = Dul * theSW(j);

          gp_Pnt aPoint;
          gp_Vec aNormal;
          theSurface.Normal(u, v, aPoint, aNormal);

          computeSInertiaOfElementaryPart(aPoint, aNormal, theLocation, aWeight, aLocalInertia);
        }

        multAndRestoreInertia(ur, aLocalInertia);
        addAndRestoreInertia(aLocalInertia, aCInertia);
      }
    }

    multAndRestoreInertia(lr, aCInertia);
    addAndRestoreInertia(aCInertia, anInertia);
  };

  const NCollection_Vector<double> aNoKnots;
  for (; theDomain.More(); theDomain.Next())
  {
    if (!theSurface.Load(theDomain.Value()))
    {
      return;
    }

    const double l1 = theSurface.FirstParameter();
    const double l2 = theSurface.LastParameter();

    int NbCGaussgp_Pnts = std::min(theSurface.IntegrationOrder(), math::GaussPointsMax());

    NbCGaussgp_Pnts = std::max(NbCGaussgp_Pnts, NbGaussgp_Pnts);

    math_Vector GaussCP(1, NbCGaussgp_Pnts);
    math_Vector GaussCW(1, NbCGaussgp_Pnts);
    math::GaussPoints(NbCGaussgp_Pnts, GaussCP);
    math::GaussWeights(NbCGaussgp_Pnts, GaussCW);

    const int anEdgeOrder = edgeBreaks(theSurface, l1, l2, aUKnots, aVKnots, aBreaks);
    if (aBreaks.Length() > 2 || !aUKnots.IsEmpty())
    {
      const int   aSpanOrder = std::min(std::max({8, aSurfOrder, anEdgeOrder}), GPM);
      math_Vector aSP(1, aSpanOrder), aSW(1, aSpanOrder);
      math::GaussPoints(aSpanOrder, aSP);
      math::GaussWeights(aSpanOrder, aSW);
      const bool isAcross = !aUKnots.IsEmpty();
      // an edge not broken at its own knots keeps, along each piece, the
      // rule it is integrated with whole, when that is the finer
      const bool isEdgeRule = anEdgeOrder == 0 && NbCGaussgp_Pnts > aSpanOrder;
      for (int b = 1; b < aBreaks.Length(); ++b)
      {
        anAlong(aBreaks(b - 1),
                aBreaks(b),
                isEdgeRule ? GaussCP : aSP,
                isEdgeRule ? GaussCW : aSW,
                isAcross ? aSP : GaussSPV,
                isAcross ? aSW : GaussSWV,
                aUKnots);
      }
      continue;
    }

    anAlong(l1, l2, GaussCP, GaussCW, GaussSPV, GaussSWV, aNoKnots);
  }

  convert(anInertia, theOutGravityCenter, theOutInertia, theOutMass);
}

//=================================================================================================

void BRepGProp_Gauss::Compute(BRepGProp_Face&   theSurface,
                              BRepGProp_Domain& theDomain,
                              const gp_Pnt&     theLocation,
                              const double      theCoeff[],
                              const bool        theIsByPoint,
                              double&           theOutMass,
                              gp_Pnt&           theOutGravityCenter,
                              gp_Mat&           theOutInertia)
{
  Standard_ASSERT_RAISE(myType == Vinert, "BRepGProp_Gauss: Incorrect type");

  double u1, v1, u2, v2;
  theSurface.Bounds(u1, u2, v1, v2);
  checkBounds(u1, u2, v1, v2);

  double _u2 = u2; // OCC104

  // In spans where a rule of GaussPointsMax() points cannot follow the
  // surface or an edge, as for the area above (#617).
  NCollection_Vector<double> aUKnots, aVKnots, aBreaks, aPieces;
  const int aSurfOrder = surfaceSpanKnots(theSurface.GetFace(), u1, u2, v1, v2, aUKnots, aVKnots);

  // the integral along the edge from theA to theB, of the integral across
  // the face from u1, in pieces at theKnots
  BRepGProp_Gauss::Inertia anInertia;
  auto anAlong = [&](const double                      theA,
                     const double                      theB,
                     const math_Vector&                theCP,
                     const math_Vector&                theCW,
                     const math_Vector&                theSP,
                     const math_Vector&                theSW,
                     const NCollection_Vector<double>& theKnots) {
    const double lm = 0.5 * (theB + theA);
    const double lr = 0.5 * (theB - theA);

    BRepGProp_Gauss::Inertia aCInertia;
    for (int i = theCP.Lower(); i <= theCP.Upper(); ++i)
    {
      const double l = lm + lr * theCP(i);

      gp_Pnt2d Puv;
      gp_Vec2d Vuv;

      theSurface.D12d(l, Puv, Vuv);

      const double aU2 = std::min(std::max(u1, Puv.X()), _u2); // OCC104
      const double v   = std::min(std::max(Puv.Y(), v1), v2);
      const double Dul = Vuv.Y() * theCW(i);

      pieces(u1, aU2, theKnots, aPieces);
      for (int p = 1; p < aPieces.Length(); ++p)
      {
        const double um = 0.5 * (aPieces(p) + aPieces(p - 1));
        const double ur = 0.5 * (aPieces(p) - aPieces(p - 1));

        BRepGProp_Gauss::Inertia aLocalInertia;
        for (int j = theSP.Lower(); j <= theSP.Upper(); ++j)
        {
          const double u       = um + ur * theSP(j);
          const double aWeight = Dul * theSW(j);

          gp_Pnt aPoint;
          gp_Vec aNormal;

          theSurface.Normal(u, v, aPoint, aNormal);

          computeVInertiaOfElementaryPart(aPoint,
                                          aNormal,
                                          theLocation,
                                          aWeight,
                                          theCoeff,
                                          theIsByPoint,
                                          aLocalInertia);
        }

        multAndRestoreInertia(ur, aLocalInertia);
        addAndRestoreInertia(aLocalInertia, aCInertia);
      }
    }

    multAndRestoreInertia(lr, aCInertia);
    addAndRestoreInertia(aCInertia, anInertia);
  };

  const NCollection_Vector<double> aNoKnots;
  for (; theDomain.More(); theDomain.Next())
  {
    if (!theSurface.Load(theDomain.Value()))
    {
      return;
    }

    const double l1 = theSurface.FirstParameter();
    const double l2 = theSurface.LastParameter();

    const int aVNbCGaussgp_Pnts = theSurface.VIntegrationOrder();

    const int aNbGaussgp_Pnts =
      std::min(std::max(theSurface.IntegrationOrder(), aVNbCGaussgp_Pnts), math::GaussPointsMax());

    math_Vector GaussP(1, aNbGaussgp_Pnts);
    math_Vector GaussW(1, aNbGaussgp_Pnts);
    math::GaussPoints(aNbGaussgp_Pnts, GaussP);
    math::GaussWeights(aNbGaussgp_Pnts, GaussW);

    const int anEdgeOrder = edgeBreaks(theSurface, l1, l2, aUKnots, aVKnots, aBreaks);
    if (aBreaks.Length() > 2 || !aUKnots.IsEmpty())
    {
      const int   aSpanOrder = std::min(std::max({8, aSurfOrder, anEdgeOrder}), GPM);
      math_Vector aSP(1, aSpanOrder), aSW(1, aSpanOrder);
      math::GaussPoints(aSpanOrder, aSP);
      math::GaussWeights(aSpanOrder, aSW);
      const bool isAcross = !aUKnots.IsEmpty();
      // an edge not broken at its own knots keeps, along each piece, the
      // rule it is integrated with whole, when that is the finer
      const bool isEdgeRule = anEdgeOrder == 0 && aNbGaussgp_Pnts > aSpanOrder;
      for (int b = 1; b < aBreaks.Length(); ++b)
      {
        anAlong(aBreaks(b - 1),
                aBreaks(b),
                isEdgeRule ? GaussP : aSP,
                isEdgeRule ? GaussW : aSW,
                isAcross ? aSP : GaussP,
                isAcross ? aSW : GaussW,
                aUKnots);
      }
      continue;
    }

    anAlong(l1, l2, GaussP, GaussW, GaussP, GaussW, aNoKnots);
  }

  convert(anInertia, theCoeff, theIsByPoint, theOutGravityCenter, theOutInertia, theOutMass);
}

//=================================================================================================

void BRepGProp_Gauss::Compute(const BRepGProp_Face& theSurface,
                              const gp_Pnt&         theLocation,
                              const double          theCoeff[],
                              const bool            theIsByPoint,
                              double&               theOutMass,
                              gp_Pnt&               theOutGravityCenter,
                              gp_Mat&               theOutInertia)
{
  double LowerU, UpperU, LowerV, UpperV;
  theSurface.Bounds(LowerU, UpperU, LowerV, UpperV);
  checkBounds(LowerU, UpperU, LowerV, UpperV);

  const int UOrder = std::min(theSurface.UIntegrationOrder(), math::GaussPointsMax());
  const int VOrder = std::min(theSurface.VIntegrationOrder(), math::GaussPointsMax());

  // Gauss points and weights
  math_Vector GaussPU(1, UOrder);
  math_Vector GaussWU(1, UOrder);
  math_Vector GaussPV(1, VOrder);
  math_Vector GaussWV(1, VOrder);

  math::GaussPoints(UOrder, GaussPU);
  math::GaussWeights(UOrder, GaussWU);
  math::GaussPoints(VOrder, GaussPV);
  math::GaussWeights(VOrder, GaussWV);

  const double um = 0.5 * add(UpperU, LowerU);
  const double vm = 0.5 * add(UpperV, LowerV);
  double       ur = 0.5 * add(UpperU, -LowerU);
  double       vr = 0.5 * add(UpperV, -LowerV);

  gp_Pnt aPoint;
  gp_Vec aNormal;

  BRepGProp_Gauss::Inertia anInertia;
  for (int j = 1; j <= VOrder; ++j)
  {
    BRepGProp_Gauss::Inertia anInertiaOfElementaryPart;
    const double             v = add(vm, mult(vr, GaussPV(j)));

    for (int i = 1; i <= UOrder; ++i)
    {
      const double aWeight = GaussWU(i);
      const double u       = add(um, mult(ur, GaussPU(i)));
      theSurface.Normal(u, v, aPoint, aNormal);

      if (myType == Vinert)
      {
        computeVInertiaOfElementaryPart(aPoint,
                                        aNormal,
                                        theLocation,
                                        aWeight,
                                        theCoeff,
                                        theIsByPoint,
                                        anInertiaOfElementaryPart);
      }
      else // Sinert
      {
        computeSInertiaOfElementaryPart(aPoint,
                                        aNormal,
                                        theLocation,
                                        aWeight,
                                        anInertiaOfElementaryPart);
      }
    }

    multAndRestoreInertia(GaussWV(j), anInertiaOfElementaryPart);
    addAndRestoreInertia(anInertiaOfElementaryPart, anInertia);
  }
  vr            = mult(vr, ur);
  anInertia.Ixx = mult(vr, anInertia.Ixx);
  anInertia.Iyy = mult(vr, anInertia.Iyy);
  anInertia.Izz = mult(vr, anInertia.Izz);
  anInertia.Ixy = mult(vr, anInertia.Ixy);
  anInertia.Ixz = mult(vr, anInertia.Ixz);
  anInertia.Iyz = mult(vr, anInertia.Iyz);

  if (myType == Vinert)
  {
    convert(anInertia, theCoeff, theIsByPoint, theOutGravityCenter, theOutInertia, theOutMass);
  }
  else // Sinert
  {
    convert(anInertia, theOutGravityCenter, theOutInertia, theOutMass);
  }

  theOutMass *= vr;
}

//=================================================================================================

void BRepGProp_Gauss::Compute(const BRepGProp_Face& theSurface,
                              const gp_Pnt&         theLocation,
                              double&               theOutMass,
                              gp_Pnt&               theOutGravityCenter,
                              gp_Mat&               theOutInertia)
{
  Standard_ASSERT_RAISE(myType == Sinert, "BRepGProp_Gauss: Incorrect type");

  Compute(theSurface, theLocation, nullptr, true, theOutMass, theOutGravityCenter, theOutInertia);
}
