// Created on: 1993-11-18
// Created by: Isabelle GRIGNON
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

#include <Adaptor2d_Curve2d.hxx>
#include <Blend_FuncInv.hxx>
#include <BRepBlend_Line.hxx>
#include <BRepLib.hxx>
#include <BRepLib_CheckCurveOnSurface.hxx>
#include <BRepTopAdaptor_TopolTool.hxx>
#include <ChFi3d_Builder.hxx>
#include <ChFi3d_Builder_0.hxx>
#include <ChFiDS_CommonPoint.hxx>
#include <ChFiDS_SurfData.hxx>
#include <NCollection_Sequence.hxx>
#include <NCollection_HSequence.hxx>
#include <ChFiDS_Stripe.hxx>
#include <NCollection_List.hxx>
#include <ChFiDS_FilSpine.hxx>
#include <ChFiDS_Spine.hxx>
#include <Geom2d_Curve.hxx>
#include <gp_Pnt2d.hxx>
#include <Precision.hxx>
#include <ShapeFix.hxx>
#include <Standard_ErrorHandler.hxx>
#include <Standard_Failure.hxx>
#include <Standard_NotImplemented.hxx>
#include <Standard_Integer.hxx>
#include <NCollection_Map.hxx>
#include <TopAbs.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopOpeBRepBuild_HBuilder.hxx>
#include <TopOpeBRepBuild_PaveSet.hxx>
#include <TopOpeBRepDS_Curve.hxx>
#include <TopOpeBRepDS_CurveExplorer.hxx>
#include <TopOpeBRepDS_CurvePointInterference.hxx>
#include <TopOpeBRepDS_DataStructure.hxx>
#include <TopOpeBRepDS_HDataStructure.hxx>
#include <TopOpeBRepDS_Interference.hxx>
#include <TopOpeBRepDS_PointIterator.hxx>

#include <atomic>
#include <exception>
#include <sstream>
#include <vector>
#include <BRep_Builder.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepLProp_SLProps.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepTools.hxx>
#include <BRep_TEdge.hxx>
#include <BRep_TVertex.hxx>
#include <BRep_Tool.hxx>
#include <Extrema_ExtPC.hxx>
#include <GeomProjLib.hxx>
#include <TopOpeBRepDS_Point.hxx>
#include <NCollection_Array1.hxx>
#include <TopOpeBRepDS_Transition.hxx>
#include <TopOpeBRepDS_Surface.hxx>
#include <TopOpeBRepDS_SolidSurfaceInterference.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Conic.hxx>
#include <Geom_Line.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Geom2d_Conic.hxx>
#include <Geom2d_Line.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <BRepTools_History.hxx>
#include <BRepTools_ReShape.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <NCollection_DataMap.hxx>
#include <NCollection_Sequence.hxx>
#include <TopExp.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Wire.hxx>
#include <NCollection_IndexedMap.hxx>
#include <BOPTools_AlgoTools.hxx>
#include <IntTools_Context.hxx>
#include <IntTools_Curve.hxx>
#include <IntTools_FaceFace.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <TopLoc_Location.hxx>
#include <ElCLib.hxx>
#include <ElSLib.hxx>
#include <Extrema_ExtPS.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <NCollection_HArray1.hxx>

#ifdef OCCT_DEBUG
  #include <OSD_Chronometer.hxx>

// variables for performances

OSD_Chronometer cl_total, cl_extent, cl_perfsetofsurf, cl_perffilletonvertex, cl_filds,
  cl_reconstruction, cl_setregul, cl_perform1corner, cl_perform2corner, cl_performatend,
  cl_perform3corner, cl_performmore3corner;

Standard_EXPORT double t_total, t_extent, t_perfsetofsurf, t_perffilletonvertex, t_filds,
  t_reconstruction, t_setregul, t_perfsetofkgen, t_perfsetofkpart, t_makextremities, t_performatend,
  t_startsol, t_performsurf, t_perform1corner, t_perform2corner, t_perform3corner,
  t_performmore3corner, t_batten, t_inter, t_sameinter, t_same, t_plate, t_approxplate,
  t_t2cornerinit, t_perf2cornerbyinter, t_chfikpartcompdata, t_cheminement, t_remplissage,
  t_t3cornerinit, t_spherique, t_torique, t_notfilling, t_filling, t_sameparam, t_computedata,
  t_completedata, t_t2cornerDS, t_t3cornerDS;

extern void ChFi3d_InitChron(OSD_Chronometer& ch);
extern void ChFi3d_ResultChron(OSD_Chronometer& ch, double& time);
extern bool ChFi3d_GettraceCHRON();
#endif

//=================================================================================================

namespace
{
std::atomic<double> THE_PLATE_G0_FALLBACK(Precision::Infinite());
std::atomic<double> THE_PLATE_G0_FALLBACK_RATIO(0.01);
std::atomic<double> THE_CORNER_SETBACK_FALLBACK(2.);
// set while Compute() runs inside the setback fallback, or the plain
// computation the fallback starts from: those do not fall back again
thread_local bool THE_IN_SETBACK_COMPUTE = false;
// SetFlatEndAtVertex(): the calling thread's
thread_local bool THE_FLAT_END_AT_VERTEX = false;
// The vertices where two fillets on no common face meet (OppositeEnds) that
// the computation running is to build with the opposite corner
// (PerformOppositeCorner), and those it met and built with the plate; both
// the computation's own (ComputeOppositeCorners)
thread_local NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> THE_OPPOSITE_CORNERS;
thread_local NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> THE_OPPOSITE_SEEN;
// The tolerances of the input's vertices and edges as the computation found
// them: a computation raises those it builds on in place (the result shares
// them), and one that failed would hand its raises on to the fallbacks' --
// #962's V36 at 0.3, whose plate set a vertex 3.5 away to 3.82, a tolerance
// the opposite corner then kept. Reset() puts them back.
thread_local NCollection_DataMap<TopoDS_Shape, double, TopTools_ShapeMapHasher>
  THE_INPUT_TOLERANCES;
} // namespace

void ChFi3d_Builder::SetPlateG0FallbackRatio(const double theRatio)
{
  THE_PLATE_G0_FALLBACK_RATIO.store(theRatio, std::memory_order_relaxed);
}

double ChFi3d_Builder::PlateG0FallbackRatio()
{
  return THE_PLATE_G0_FALLBACK_RATIO.load(std::memory_order_relaxed);
}

double ChFi3d_SetPlateG0FallbackRatio(const double theRatio)
{
  return THE_PLATE_G0_FALLBACK_RATIO.exchange(theRatio, std::memory_order_relaxed);
}

void ChFi3d_Builder::SetPlateG0Fallback(const double theDistance)
{
  THE_PLATE_G0_FALLBACK.store(theDistance, std::memory_order_relaxed);
}

double ChFi3d_Builder::PlateG0Fallback()
{
  return THE_PLATE_G0_FALLBACK.load(std::memory_order_relaxed);
}

double ChFi3d_SetPlateG0Fallback(const double theDistance)
{
  return THE_PLATE_G0_FALLBACK.exchange(theDistance, std::memory_order_relaxed);
}

void ChFi3d_Builder::SetCornerSetbackFallback(const double theMultiple)
{
  THE_CORNER_SETBACK_FALLBACK.store(theMultiple, std::memory_order_relaxed);
}

double ChFi3d_Builder::CornerSetbackFallback()
{
  return THE_CORNER_SETBACK_FALLBACK.load(std::memory_order_relaxed);
}

double ChFi3d_SetCornerSetbackFallback(const double theMultiple)
{
  return THE_CORNER_SETBACK_FALLBACK.exchange(theMultiple, std::memory_order_relaxed);
}

void ChFi3d_Builder::SetFlatEndAtVertex(const bool theToUse)
{
  THE_FLAT_END_AT_VERTEX = theToUse;
}

bool ChFi3d_Builder::FlatEndAtVertex()
{
  return THE_FLAT_END_AT_VERTEX;
}

int ChFi3d_SetFlatEndAtVertex(const int theToUse)
{
  const bool wasSet      = THE_FLAT_END_AT_VERTEX;
  THE_FLAT_END_AT_VERTEX = theToUse != 0;
  return wasSet ? 1 : 0;
}

//=================================================================================================

static void CompleteDS(TopOpeBRepDS_DataStructure& DStr, const TopoDS_Shape& S)
{
  ChFiDS_Map MapEW, MapFS;
  MapEW.Fill(S, TopAbs_EDGE, TopAbs_WIRE);
  MapFS.Fill(S, TopAbs_FACE, TopAbs_SHELL);

  TopExp_Explorer ExpE;
  for (ExpE.Init(S, TopAbs_EDGE); ExpE.More(); ExpE.Next())
  {
    const TopoDS_Edge& E       = TopoDS::Edge(ExpE.Current());
    bool               hasgeom = DStr.HasGeometry(E);
    if (hasgeom)
    {
      const NCollection_List<TopoDS_Shape>&    WireListAnc = MapEW(E);
      NCollection_List<TopoDS_Shape>::Iterator itaW(WireListAnc);
      while (itaW.More())
      {
        const TopoDS_Shape& WireAnc = itaW.Value();
        DStr.AddShape(WireAnc);
        itaW.Next();
      }
    }
  }

  TopExp_Explorer ExpF;
  for (ExpF.Init(S, TopAbs_FACE); ExpF.More(); ExpF.Next())
  {
    const TopoDS_Face& F       = TopoDS::Face(ExpF.Current());
    bool               hasgeom = DStr.HasGeometry(F);
    if (hasgeom)
    {
      const NCollection_List<TopoDS_Shape>&    ShellListAnc = MapFS(F);
      NCollection_List<TopoDS_Shape>::Iterator itaS(ShellListAnc);
      while (itaS.More())
      {
        const TopoDS_Shape& ShellAnc = itaS.Value();
        DStr.AddShape(ShellAnc);
        itaS.Next();
      }
    }
  }

  // set the range on the DS Curves
  for (int ic = 1; ic <= DStr.NbCurves(); ic++)
  {
    double parmin = RealLast(), parmax = RealFirst();
    const NCollection_List<occ::handle<TopOpeBRepDS_Interference>>& LI =
      DStr.CurveInterferences(ic);
    for (TopOpeBRepDS_PointIterator it(LI); it.More(); it.Next())
    {
      double par = it.Parameter();
      parmin     = std::min(parmin, par);
      parmax     = std::max(parmax, par);
    }
    DStr.ChangeCurve(ic).SetRange(parmin, parmax);
  }
}

//=================================================================================================

ChFi3d_Builder::~ChFi3d_Builder() = default;

//=================================================================================================

void ChFi3d_Builder::ExtentAnalyse()
{
  int nbedges, nbs;
  for (int iv = 1; iv <= myVDataMap.Extent(); iv++)
  {
    nbs                      = myVDataMap(iv).Extent();
    const TopoDS_Vertex& Vtx = myVDataMap.FindKey(iv);
    // nbedges = ChFi3d_NumberOfEdges(Vtx, myVEMap);
    nbedges = ChFi3d_NumberOfSharpEdges(Vtx, myVEMap, myEFMap);
    switch (nbs)
    {
      case 1:
        ExtentOneCorner(Vtx, myVDataMap.FindFromIndex(iv).First());
        break;
      case 2:
        if (nbedges <= 3)
        {
          ExtentTwoCorner(Vtx, myVDataMap.FindFromIndex(iv));
        }
        break;
      case 3:
        if (nbedges <= 3)
        {
          ExtentThreeCorner(Vtx, myVDataMap.FindFromIndex(iv));
        }
        break;
      default:
        break;
    }
  }
}

//=================================================================================================
// Adds to <theOld>, the input's shapes of type <theType>, each one of that type
// in <theResult> that is a thawed copy of one of them (TopoDS_TShape::Thawed):
// a frozen input's shape the fillet copied for a vertex it changed stands for
// the shape it was copied from, as that shape would be, changed in place.

template <class TheMap>
static void ChFi3d_AddThawed(TheMap&                theOld,
                             const TopoDS_Shape&    theResult,
                             const TopAbs_ShapeEnum theType)
{
  for (TopExp_Explorer ex(theResult, theType); ex.More(); ex.Next())
  {
    const TopoDS_Shape& aS = ex.Current();
    if (!aS.TShape()->Thawed())
    {
      continue;
    }
    TopoDS_Shape anOriginal = aS;
    for (occ::handle<TopoDS_TShape> aFrom = TopoDS_TShape::ThawedFrom(aS.TShape().get());
         !aFrom.IsNull();
         aFrom = TopoDS_TShape::ThawedFrom(aFrom.get()))
    {
      anOriginal.TShape(aFrom);
      if (theOld.Contains(anOriginal))
      {
        theOld.Add(aS);
        break;
      }
    }
  }
}

//=================================================================================================
// An edge the fillet made whose curve ends off its vertex by more than the
// edge's tolerance, but within the vertex's, takes that distance: at a corner
// the extension of an edge ends at a point computed on the fillet's line, not
// on the extension, and the cut there leaves it tangent -- BRepCheck's wire
// check excuses a crossing near the vertex only while both edges keep within
// twice their tolerance of the chord from the vertex, which the vertex's
// offset alone exceeded (#631's Fillet002: 3.7e-7 off, tolerance 1e-7).
// Edges of the input are left alone: the input shares them.

//! The thawed copies <thePart> of <theResult> holds that <theSpread> has not
//! seen, put in their originals' place in the whole of <theResult>. BRepLib
//! raises a frozen vertex or edge of the part it is given on a copy, and puts
//! the copy in that part alone (BRepLib::UpdateTolerances' copy-on-write):
//! another face of the result that shares the original -- one of the input's
//! kept as it was -- would keep it, and the wire the two faces share would not
//! close at it (#273's flat end, 0.3 on its edge 39, frozen). Unfrozen, the
//! vertex is raised in place, in every face that holds it.
static void ChFi3d_SpreadThawed(const TopoDS_Shape&                    theResult,
                                const TopoDS_Shape&                    thePart,
                                NCollection_Map<const TopoDS_TShape*>& theSpread)
{
  NCollection_List<TopoDS_Shape> aCopies;
  for (const TopAbs_ShapeEnum aType : {TopAbs_VERTEX, TopAbs_EDGE})
  {
    for (TopExp_Explorer ex(thePart, aType); ex.More(); ex.Next())
    {
      if (ex.Current().TShape()->Thawed() && theSpread.Add(ex.Current().TShape().get()))
      {
        aCopies.Append(ex.Current());
      }
    }
  }
  if (!aCopies.IsEmpty())
  {
    BRepLib_PutThawed(theResult, aCopies);
  }
}

//=================================================================================================

static void ChFi3d_EdgesCoverTheirEnds(const TopoDS_Shape& theResult, const TopoDS_Shape& theInput)
{
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> anOld;
  for (TopExp_Explorer ex(theInput, TopAbs_EDGE); ex.More(); ex.Next())
  {
    anOld.Add(ex.Current());
  }
  ChFi3d_AddThawed(anOld, theResult, TopAbs_EDGE);
  BRep_Builder                                           aB;
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aDone;
  for (TopExp_Explorer ex(theResult, TopAbs_EDGE); ex.More(); ex.Next())
  {
    const TopoDS_Edge& anE = TopoDS::Edge(ex.Current());
    if (anOld.Contains(anE) || !aDone.Add(anE) || BRep_Tool::Degenerated(anE))
    {
      continue;
    }
    double                         aF, aL;
    const occ::handle<Geom_Curve>& aC = BRep_Tool::Curve(anE, aF, aL);
    if (aC.IsNull())
    {
      continue;
    }
    double aTol = BRep_Tool::Tolerance(anE);
    for (TopExp_Explorer exv(anE, TopAbs_VERTEX); exv.More(); exv.Next())
    {
      const TopoDS_Vertex& aV = TopoDS::Vertex(exv.Current());
      // the curve as BRep_Tool::Curve(E, f, l) gives it is located already
      const gp_Pnt aP   = aC->Value(BRep_Tool::Parameter(aV, anE));
      const double aGap = aP.Distance(BRep_Tool::Pnt(aV));
      if (aGap > aTol && aGap <= BRep_Tool::Tolerance(aV))
      {
        aTol = aGap * (1. + 1.e-9);
      }
    }
    if (aTol > BRep_Tool::Tolerance(anE))
    {
      aB.UpdateEdge(anE, aTol);
    }
  }
}

//=================================================================================================
// An edge the fillet made covers its curves on the faces of the result. The
// tolerance the last pass leaves (BRepLib::SameParameter) is the largest
// distance at its samples and a twentieth more; a curve between two fillets
// has its ends put on the points known on each, and is farthest off its
// surface a little inside the end, between two samples (#962's pocket,
// edges 10 and 14 at radius 1: 2.51e-6 at 0.003 of the curve for 2.18e-6 at
// its start, tolerance 2.29e-6 -- BOP's check, which looks for the largest
// distance, calls the curve invalid). Here the distance is taken at ten
// times the samples and followed to its peak from the largest; only an edge
// that comes out within a fiftieth of its tolerance, or over it, is looked
// at the way the check does (BRepLib_CheckCurveOnSurface, a search too
// long for every edge). For a curve or a curve on the surface that is not
// a line or a conic: those are as far off between the samples as at them.
// Edges of the input are left alone, and an edge at a vertex that cannot
// grow.

static double ChFi3d_DistanceToCurveOnSurface(const TopoDS_Edge& theE, const TopoDS_Face& theF)
{
  const BRepAdaptor_Curve aC(theE), aCS(theE, theF);
  const double            aF = aC.FirstParameter(), aL = aC.LastParameter();
  auto                    aDist2 = [&](const double theT) {
    return aC.Value(theT).SquareDistance(aCS.Value(theT));
  };
  const int aNb   = 220;
  const double aStep = (aL - aF) / aNb;
  int       iMax  = 0;
  double    aMax  = -1.;
  for (int i = 0; i <= aNb; i++)
  {
    const double aD = aDist2(aF + i * aStep);
    if (aD > aMax)
    {
      aMax = aD;
      iMax = i;
    }
  }
  // to the peak beside the largest sample
  double a = aF + std::max(0, iMax - 1) * aStep, b = aF + std::min(aNb, iMax + 1) * aStep;
  for (int it = 0; it < 30; it++)
  {
    const double m1 = a + (b - a) / 3., m2 = b - (b - a) / 3.;
    const double d1 = aDist2(m1), d2 = aDist2(m2);
    aMax = std::max(aMax, std::max(d1, d2));
    if (d1 < d2)
    {
      a = m1;
    }
    else
    {
      b = m2;
    }
  }
  return sqrt(aMax);
}

static void ChFi3d_EdgesCoverTheirPCurves(const TopoDS_Shape& theResult,
                                          const TopoDS_Shape& theInput)
{
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> anOld;
  for (TopExp_Explorer ex(theInput, TopAbs_EDGE); ex.More(); ex.Next())
  {
    anOld.Add(ex.Current());
  }
  ChFi3d_AddThawed(anOld, theResult, TopAbs_EDGE);
  BRep_Builder aB;
  for (TopExp_Explorer exf(theResult, TopAbs_FACE); exf.More(); exf.Next())
  {
    const TopoDS_Face& aFace = TopoDS::Face(exf.Current());
    for (TopExp_Explorer ex(aFace, TopAbs_EDGE); ex.More(); ex.Next())
    {
      const TopoDS_Edge& anE = TopoDS::Edge(ex.Current());
      if (anOld.Contains(anE) || BRep_Tool::Degenerated(anE) || anE.Immutable())
      {
        continue;
      }
      double                      aF, aL;
      occ::handle<Geom_Curve>     aC  = BRep_Tool::Curve(anE, aF, aL);
      occ::handle<Geom2d_Curve>   aPC = BRep_Tool::CurveOnSurface(anE, aFace, aF, aL);
      if (aC.IsNull() || aPC.IsNull())
      {
        continue;
      }
      if (aC->IsKind(STANDARD_TYPE(Geom_TrimmedCurve)))
      {
        aC = occ::down_cast<Geom_TrimmedCurve>(aC)->BasisCurve();
      }
      if (aPC->IsKind(STANDARD_TYPE(Geom2d_TrimmedCurve)))
      {
        aPC = occ::down_cast<Geom2d_TrimmedCurve>(aPC)->BasisCurve();
      }
      const bool isPlain3d =
        aC->IsKind(STANDARD_TYPE(Geom_Line)) || aC->IsKind(STANDARD_TYPE(Geom_Conic));
      const bool isPlain2d =
        aPC->IsKind(STANDARD_TYPE(Geom2d_Line)) || aPC->IsKind(STANDARD_TYPE(Geom2d_Conic));
      if (isPlain3d && isPlain2d)
      {
        continue;
      }
      double aDist = 0.;
      try
      {
        OCC_CATCH_SIGNALS
        if (!BRep_Tool::SameParameter(anE)
            || 1.02 * ChFi3d_DistanceToCurveOnSurface(anE, aFace) <= BRep_Tool::Tolerance(anE))
        {
          continue;
        }
        BRepLib_CheckCurveOnSurface aCheck(anE, aFace);
        aCheck.Perform();
        if (!aCheck.IsDone())
        {
          continue;
        }
        aDist = aCheck.MaxDistance() * (1. + 1.e-6);
      }
      catch (const Standard_Failure&)
      {
        continue;
      }
      if (aDist <= BRep_Tool::Tolerance(anE))
      {
        continue;
      }
      bool isFree = true;
      for (TopExp_Explorer exv(anE, TopAbs_VERTEX); exv.More() && isFree; exv.Next())
      {
        const TopoDS_Vertex& aV = TopoDS::Vertex(exv.Current());
        isFree                  = !aV.Immutable() || BRep_Tool::Tolerance(aV) >= aDist;
      }
      if (!isFree)
      {
        continue;
      }
      aB.UpdateEdge(anE, aDist);
      for (TopExp_Explorer exv(anE, TopAbs_VERTEX); exv.More(); exv.Next())
      {
        const TopoDS_Vertex& aV = TopoDS::Vertex(exv.Current());
        if (BRep_Tool::Tolerance(aV) < aDist)
        {
          aB.UpdateVertex(aV, aDist);
        }
      }
    }
  }
}

//=================================================================================================
// A face of the result whose one wire runs through a point twice, at two
// vertices there, the corner's new one and one of the input (#523 at the
// cylinder's radius: the fillet's line on the top is tangent to the
// cylinder's circle at its vertex, and the top is the box's rectangle and
// the cylinder's whole disc, touching there). BRepCheck lets the wire be,
// its vertices 4e-15 apart; BOP's check calls the vertex self-intersecting,
// and made one vertex the wire is unorientable. The face is two regions
// meeting at a point: the corner's vertex is replaced by the input's, and
// the face split there in two, each bounded by one loop of its wire -- only
// where both loops bound a region, not where a loop round a hole touches
// the outer one. The same with two vertices the fillet made, and with one
// of the loops a single edge from one to the other, and with a wire through
// one vertex twice, the face's other wires to the piece they lie in (#876's
// post, the fillet's end flat on the plate's rim, MakeRimFlatEnd)
// (FindPinch). Returns the history of the change, null when there was none.

namespace
{
// a face pinched at a point, and where its wire passes the point
struct PinchedFace
{
  TopoDS_Face                        Face;    // forward
  TopoDS_Vertex                      Keep;    // the input's vertex there
  TopoDS_Vertex                      Drop;    // the other
  NCollection_Sequence<TopoDS_Shape> Ordered; // the wire's edges in order
  int                                K1 = 0;  // the edges ending at the point
  int                                K2 = 0;
  NCollection_List<TopoDS_Shape>     Holes; // the face's other wires
};

// The wire of <theP> from the edge after <theFrom>, <theNb> edges, with the
// edges as <theReShape> has them.
TopoDS_Wire PinchLoop(const PinchedFace&                    theP,
                      const int                             theFrom,
                      const int                             theNb,
                      const occ::handle<BRepTools_ReShape>& theReShape)
{
  BRep_Builder aB;
  TopoDS_Wire  aW;
  aB.MakeWire(aW);
  for (int n = 0; n < theNb; n++)
  {
    const TopoDS_Shape& anE = theP.Ordered((theFrom + n) % theP.Ordered.Length() + 1);
    aB.Add(aW, theReShape.IsNull() ? anE : theReShape->Value(anE));
  }
  aW.Closed(true);
  return aW;
}

// <theF>, one wire, pinched: its wire passes twice through a point, at two
// vertices there, not both the input's, with no edge between them but one
// that is itself a loop of the wire -- longer than ten times what the
// vertices cover: the whole circle a fillet's line makes round a disc that
// touches the rest of the face (#523's bottom, the cylinder's rim and the
// box's edge filleted at radius 1) --, at one point of the face's
// parameters, and each loop of the wire from there bounds a region of its
// own. Kept: the input's vertex, or of two the fillet made the looser.
bool FindPinch(const TopoDS_Face&                                                   theF,
               const NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>& theInputV,
               PinchedFace&                                                         theP)
{
  TopoDS_Wire aW;
  int         aNbW = 0;
  theP.Holes.Clear();
  for (TopExp_Explorer exW(theF, TopAbs_WIRE); exW.More(); exW.Next(), aNbW++)
  {
    aW = TopoDS::Wire(exW.Current());
  }
  if (aNbW != 1)
  {
    // with holes, only a wire through one vertex twice (below): that wire,
    // the others holes in one piece or the other
    TopoDS_Wire aPinched;
    for (TopExp_Explorer exW(theF, TopAbs_WIRE); exW.More(); exW.Next())
    {
      NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aSeen;
      bool                                                   isTwice = false;
      for (TopoDS_Iterator itE(exW.Current()); itE.More() && !isTwice; itE.Next())
      {
        isTwice = !aSeen.Add(TopExp::LastVertex(TopoDS::Edge(itE.Value()), true));
      }
      if (isTwice && !aPinched.IsNull())
      {
        return false;
      }
      if (isTwice)
      {
        aPinched = TopoDS::Wire(exW.Current());
      }
      else
      {
        theP.Holes.Append(exW.Current());
      }
    }
    if (aPinched.IsNull())
    {
      return false;
    }
    aW = aPinched;
  }
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aWV;
  TopExp::MapShapes(aW, TopAbs_VERTEX, aWV);
  int           aNbPinch = 0;
  TopoDS_Vertex aV1, aV2;
  for (int i = 1; i <= aWV.Extent(); i++)
  {
    const TopoDS_Vertex& aVi = TopoDS::Vertex(aWV(i));
    for (int j = i + 1; j <= aWV.Extent(); j++)
    {
      const TopoDS_Vertex& aVj = TopoDS::Vertex(aWV(j));
      if (BRep_Tool::Pnt(aVi).Distance(BRep_Tool::Pnt(aVj))
          > std::max(BRep_Tool::Tolerance(aVi), BRep_Tool::Tolerance(aVj)))
      {
        continue;
      }
      bool isJoined = false;
      for (TopExp_Explorer exE(aW, TopAbs_EDGE); exE.More() && !isJoined; exE.Next())
      {
        const TopoDS_Edge& anE = TopoDS::Edge(exE.Current());
        TopoDS_Vertex      aE1, aE2;
        TopExp::Vertices(anE, aE1, aE2);
        isJoined = ((aE1.IsSame(aVi) && aE2.IsSame(aVj)) || (aE1.IsSame(aVj) && aE2.IsSame(aVi)))
                   && (!BRep_Tool::IsGeometric(anE)
                       || GCPnts_AbscissaPoint::Length(BRepAdaptor_Curve(anE))
                            <= 10. * (BRep_Tool::Tolerance(aVi) + BRep_Tool::Tolerance(aVj)));
      }
      if (!isJoined)
      {
        aNbPinch++;
        aV1 = aVi;
        aV2 = aVj;
      }
    }
  }
  if (aNbPinch > 1 || (aNbPinch == 1 && theInputV.Contains(aV1) && theInputV.Contains(aV2))
      || (aNbPinch == 1 && aNbW != 1))
  {
    return false;
  }
  theP.Ordered.Clear();
  for (BRepTools_WireExplorer wex(aW, theF); wex.More(); wex.Next())
  {
    theP.Ordered.Append(wex.Current());
  }
  if (aNbPinch == 0)
  {
    // Or one vertex the wire passes twice, a fillet's end flat where its
    // line runs off the face over an edge, short of the vertex the face
    // narrows to (MakeRimFlatEnd, #876's post): the line and the edge meet
    // there, and so do the flat face's side and the edge's piece beyond
    NCollection_Sequence<TopoDS_Shape> anAll;
    for (TopoDS_Iterator itE(aW); itE.More(); itE.Next())
    {
      anAll.Append(itE.Value());
    }
    NCollection_IndexedDataMap<TopoDS_Shape, int, TopTools_ShapeMapHasher> aNbEnd;
    for (int k = 1; k <= anAll.Length(); k++)
    {
      const TopoDS_Vertex aV = TopExp::LastVertex(TopoDS::Edge(anAll(k)), true);
      if (aNbEnd.Contains(aV))
      {
        aNbEnd.ChangeFromKey(aV)++;
      }
      else
      {
        aNbEnd.Add(aV, 1);
      }
    }
    // the first such vertex; the pieces are looked at again for the next
    // (a fillet made flat so at both its ends on one face)
    for (int i = 1; i <= aNbEnd.Extent(); i++)
    {
      if (aNbEnd(i) > 2)
      {
        return false;
      }
      if (aNbEnd(i) == 2 && aNbPinch == 0 && !theInputV.Contains(aNbEnd.FindKey(i)))
      {
        aNbPinch++;
        aV1 = aV2 = TopoDS::Vertex(aNbEnd.FindKey(i));
      }
    }
    if (aNbPinch != 1)
    {
      return false;
    }
    // The wire explorer stops at the vertex, its first loop closed: the two
    // loops from it, each from one of its two edges leaving the vertex back
    // to it, one after the other
    std::vector<bool> isUsed(anAll.Length() + 1, false);
    theP.Ordered.Clear();
    for (int aLoop = 0; aLoop < 2; aLoop++)
    {
      int aCur = 0;
      for (int k = 1; k <= anAll.Length() && aCur == 0; k++)
      {
        if (!isUsed[k] && TopExp::FirstVertex(TopoDS::Edge(anAll(k)), true).IsSame(aV1))
        {
          aCur = k;
        }
      }
      while (aCur != 0)
      {
        isUsed[aCur] = true;
        theP.Ordered.Append(anAll(aCur));
        const TopoDS_Vertex aV = TopExp::LastVertex(TopoDS::Edge(anAll(aCur)), true);
        if (aV.IsSame(aV1))
        {
          break;
        }
        // at another such vertex either way on (the loop left is put in
        // below)
        int aNext = 0;
        for (int k = 1; k <= anAll.Length() && aNext == 0; k++)
        {
          if (!isUsed[k] && TopExp::FirstVertex(TopoDS::Edge(anAll(k)), true).IsSame(aV))
          {
            aNext = k;
          }
        }
        if (aNext == 0)
        {
          return false;
        }
        aCur = aNext;
      }
    }
    // a loop left from another vertex the wire passes twice, put in where
    // the wire passes that vertex
    for (int aPass = 0; aPass < anAll.Length() && theP.Ordered.Length() < anAll.Length(); aPass++)
    {
      int aStart = 0, anAt = 0;
      for (int k = 1; k <= anAll.Length() && anAt == 0; k++)
      {
        const TopoDS_Vertex aV = TopExp::FirstVertex(TopoDS::Edge(anAll(k)), true);
        for (int j = 1; j <= theP.Ordered.Length() && anAt == 0 && !isUsed[k]; j++)
        {
          anAt = TopExp::LastVertex(TopoDS::Edge(theP.Ordered(j)), true).IsSame(aV) ? j : 0;
        }
        aStart = anAt != 0 ? k : 0;
      }
      if (anAt == 0)
      {
        return false;
      }
      const TopoDS_Vertex aW0 = TopExp::FirstVertex(TopoDS::Edge(anAll(aStart)), true);
      NCollection_Sequence<TopoDS_Shape> aSub;
      for (int aCur = aStart; aCur != 0;)
      {
        isUsed[aCur] = true;
        aSub.Append(anAll(aCur));
        const TopoDS_Vertex aV = TopExp::LastVertex(TopoDS::Edge(anAll(aCur)), true);
        aCur                   = 0;
        for (int k = 1; k <= anAll.Length() && aCur == 0 && !aV.IsSame(aW0); k++)
        {
          if (!isUsed[k] && TopExp::FirstVertex(TopoDS::Edge(anAll(k)), true).IsSame(aV))
          {
            aCur = k;
          }
        }
        if (aCur == 0 && !aV.IsSame(aW0))
        {
          return false;
        }
      }
      theP.Ordered.InsertAfter(anAt, aSub);
    }
  }
  const bool isFirstKept =
    theInputV.Contains(aV1)
    || (!theInputV.Contains(aV2) && BRep_Tool::Tolerance(aV1) >= BRep_Tool::Tolerance(aV2));
  theP.Face = theF;
  theP.Keep = isFirstKept ? aV1 : aV2;
  theP.Drop = isFirstKept ? aV2 : aV1;
  int aNbE = 0;
  for (TopoDS_Iterator itE(aW); itE.More(); itE.Next())
  {
    aNbE++;
  }
  if (theP.Ordered.Length() != aNbE)
  {
    return false;
  }
  int aNbAt = 0;
  theP.K1 = theP.K2 = 0;
  for (int k = 1; k <= theP.Ordered.Length(); k++)
  {
    const TopoDS_Vertex aV = TopExp::LastVertex(TopoDS::Edge(theP.Ordered(k)), true);
    if (aV.IsSame(theP.Keep) || aV.IsSame(theP.Drop))
    {
      aNbAt++;
      if (theP.K1 == 0)
      {
        theP.K1 = k;
      }
      else
      {
        theP.K2 = k;
      }
    }
  }
  if (aNbAt != 2)
  {
    return false;
  }
  // one point of the face's parameters: on a closed surface the wire can
  // pass a point twice a period apart, as at a seam
  auto anEndUV = [&theF](const TopoDS_Shape& theE) {
    gp_Pnt2d aFirst, aLast;
    BRep_Tool::UVPoints(TopoDS::Edge(theE), theF, aFirst, aLast);
    return theE.Orientation() == TopAbs_REVERSED ? aFirst : aLast;
  };
  const gp_Pnt2d            aUV1 = anEndUV(theP.Ordered(theP.K1));
  const gp_Pnt2d            aUV2 = anEndUV(theP.Ordered(theP.K2));
  const BRepAdaptor_Surface aS(theF, false);
  const double              aTol = 10. * std::max(BRep_Tool::Tolerance(theP.Keep),
                                                  BRep_Tool::Tolerance(theP.Drop));
  if (std::abs(aUV1.X() - aUV2.X()) > aS.UResolution(aTol)
      || std::abs(aUV1.Y() - aUV2.Y()) > aS.VResolution(aTol))
  {
    return false;
  }
  // each loop a region, not a hole touching the other, and longer than ten
  // times what the vertices cover: not a loop of one edge 3.2e-4 long within
  // a vertex of tolerance 0.011, a loose fillet's (#273's handle at 0.8)
  const int aNb = theP.Ordered.Length();
  for (int aLoop = 0; aLoop < 2; aLoop++)
  {
    const int aFrom   = aLoop == 0 ? theP.K1 : theP.K2;
    const int aLen    = aLoop == 0 ? theP.K2 - theP.K1 : aNb - theP.K2 + theP.K1;
    double    aLength = 0.;
    for (int n = 0; n < aLen; n++)
    {
      const TopoDS_Edge& anE = TopoDS::Edge(theP.Ordered((aFrom + n) % aNb + 1));
      if (BRep_Tool::IsGeometric(anE))
      {
        aLength += GCPnts_AbscissaPoint::Length(BRepAdaptor_Curve(anE));
      }
    }
    if (aLength <= 10. * (BRep_Tool::Tolerance(theP.Keep) + BRep_Tool::Tolerance(theP.Drop)))
    {
      return false;
    }
    TopoDS_Face aPiece = TopoDS::Face(theF.EmptyCopied());
    BRep_Builder().Add(aPiece, PinchLoop(theP, aFrom, aLen, nullptr));
    if (BRepTopAdaptor_FClass2d(aPiece, Precision::PConfusion()).PerformInfinitePoint()
        != TopAbs_OUT)
    {
      return false;
    }
  }
  return true;
}
} // namespace

static occ::handle<BRepTools_History> ChFi3d_SplitPinchedFaces(TopoDS_Shape&       theResult,
                                                              const TopoDS_Shape& theInput)
{
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> anInputV;
  TopExp::MapShapes(theInput, TopAbs_VERTEX, anInputV);
  ChFi3d_AddThawed(anInputV, theResult, TopAbs_VERTEX);
  NCollection_List<PinchedFace>                                             aPinches;
  NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher> aDropToKeep;
  for (TopExp_Explorer exF(theResult, TopAbs_FACE); exF.More(); exF.Next())
  {
    PinchedFace aP;
    if (!FindPinch(TopoDS::Face(exF.Current().Oriented(TopAbs_FORWARD)), anInputV, aP)
        || (aDropToKeep.IsBound(aP.Drop) && !aDropToKeep(aP.Drop).IsSame(aP.Keep)))
    {
      continue;
    }
    if (!aP.Drop.IsSame(aP.Keep))
    {
      aDropToKeep.Bind(aP.Drop, aP.Keep);
    }
    aPinches.Append(aP);
  }
  if (aPinches.IsEmpty())
  {
    return occ::handle<BRepTools_History>();
  }

  // the vertices made one, everywhere; the edges at a vertex dropped take
  // the vertex kept at the same parameter
  BRep_Builder                   aB;
  occ::handle<BRepTools_ReShape> aMerge = new BRepTools_ReShape;
  // a frozen vertex kept: its copy raised, in the result once it is made
  NCollection_DataMap<TopoDS_Shape, double, TopTools_ShapeMapHasher> aRaised;
  for (NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>::Iterator it(
         aDropToKeep);
       it.More();
       it.Next())
  {
    const TopoDS_Vertex& aDrop = TopoDS::Vertex(it.Key());
    const TopoDS_Vertex& aKeep = TopoDS::Vertex(it.Value());
    const double         aTol =
      BRep_Tool::Pnt(aDrop).Distance(BRep_Tool::Pnt(aKeep)) + BRep_Tool::Tolerance(aDrop);
    ChFi3d_RaiseKept(aKeep, aTol, aRaised);
    aMerge->Replace(aDrop.Oriented(TopAbs_FORWARD), aKeep.Oriented(TopAbs_FORWARD));
  }
  NCollection_List<TopoDS_Shape> anEdges, aKeeps;
  NCollection_List<double>       aParams;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aResultEdges;
  TopExp::MapShapes(theResult, TopAbs_EDGE, aResultEdges);
  for (int i = 1; i <= aResultEdges.Extent(); i++)
  {
    const TopoDS_Edge& anE = TopoDS::Edge(aResultEdges(i));
    // a loop from the vertex dropped to the one kept is closed on the one
    // kept, its range as it is
    bool isLoop = false;
    for (TopoDS_Iterator itV(anE); itV.More() && !isLoop; itV.Next())
    {
      if (aDropToKeep.IsBound(itV.Value()))
      {
        for (TopoDS_Iterator itO(anE); itO.More() && !isLoop; itO.Next())
        {
          isLoop = itO.Value().IsSame(aDropToKeep(itV.Value()));
        }
      }
    }
    for (TopoDS_Iterator itV(anE); itV.More() && !isLoop; itV.Next())
    {
      if (aDropToKeep.IsBound(itV.Value()))
      {
        anEdges.Append(anE);
        aKeeps.Append(aDropToKeep(itV.Value()));
        aParams.Append(BRep_Tool::Parameter(TopoDS::Vertex(itV.Value()), anE));
      }
    }
  }
  const TopoDS_Shape                       aMerged = aMerge->Apply(theResult);
  NCollection_List<TopoDS_Shape>::Iterator itK(aKeeps);
  NCollection_List<double>::Iterator       itP(aParams);
  for (NCollection_List<TopoDS_Shape>::Iterator itE(anEdges); itE.More();
       itE.Next(), itK.Next(), itP.Next())
  {
    const TopoDS_Vertex& aKeep = TopoDS::Vertex(itK.Value());
    aB.UpdateVertex(aKeep,
                    itP.Value(),
                    TopoDS::Edge(aMerge->Value(itE.Value())),
                    BRep_Tool::Tolerance(aKeep));
  }

  // each face split there in two, a face for each loop of its wire
  occ::handle<BRepTools_ReShape> aSplit = new BRepTools_ReShape;
  for (NCollection_List<PinchedFace>::Iterator itF(aPinches); itF.More(); itF.Next())
  {
    const PinchedFace& aP  = itF.Value();
    const int          aNb = aP.Ordered.Length();
    const TopoDS_Face  aMF = TopoDS::Face(aMerge->Value(aP.Face).Oriented(TopAbs_FORWARD));
    TopoDS_Compound    aPieces;
    aB.MakeCompound(aPieces);
    NCollection_List<TopoDS_Shape> aHoles = aP.Holes;
    for (int aLoop = 0; aLoop < 2; aLoop++)
    {
      TopoDS_Face aPiece = TopoDS::Face(aMF.EmptyCopied());
      aB.Add(aPiece,
             aLoop == 0 ? PinchLoop(aP, aP.K1, aP.K2 - aP.K1, aMerge)
                        : PinchLoop(aP, aP.K2, aNb - aP.K2 + aP.K1, aMerge));
      // the holes within this loop, the rest the other's
      const BRepTopAdaptor_FClass2d aCl(aPiece, Precision::PConfusion());
      for (NCollection_List<TopoDS_Shape>::Iterator itH(aHoles); itH.More();)
      {
        TopoDS_Iterator itE(itH.Value());
        gp_Pnt2d        aUV1, aUV2;
        BRep_Tool::UVPoints(TopoDS::Edge(itE.Value()), aMF, aUV1, aUV2);
        if (aLoop == 1 || aCl.Perform(aUV1) == TopAbs_IN)
        {
          aB.Add(aPiece, aMerge->Value(itH.Value()));
          aHoles.Remove(itH);
        }
        else
        {
          itH.Next();
        }
      }
      aB.Add(aPieces, aPiece);
    }
    aSplit->Replace(aMF, aPieces);
  }
  theResult = aSplit->Apply(aMerged);
  ChFi3d_ThawInto(theResult, aRaised);
  occ::handle<BRepTools_History> aHistory = aMerge->History();
  aHistory->Merge(aSplit->History());
  return aHistory;
}

//=================================================================================================
// A face and the one MakeRimRunOn made on its surface carried up beside it
// (#876's plate side and its rise to the fillet run on to the vertex), two
// faces of one surface across the edge between them, the rim's piece: made
// one, nothing else -- the one's wire with that edge replaced by the
// other's edges, which take their curves on the one's surface (the same
// surface, a copy). Not ShapeUpgrade_UnifySameDomain, which also lowered
// the tolerance of edges it kept below what BOP's check takes (#876's model
// at 2). The edges are found by their middles (ChFi3d_RunOnSeams), which
// are cleared. Returns the history of the change, null when there was none.

static occ::handle<BRepTools_History> ChFi3d_MergeRunOnSeams(TopoDS_Shape& theResult,
                                                             const double  theTol)
{
  NCollection_List<gp_Pnt> aSeams = ChFi3d_RunOnSeams();
  ChFi3d_RunOnSeams().Clear();
  if (aSeams.IsEmpty())
  {
    return occ::handle<BRepTools_History>();
  }
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    anEF;
  TopExp::MapShapesAndAncestors(theResult, TopAbs_EDGE, TopAbs_FACE, anEF);
  occ::handle<BRepTools_ReShape>                         aReShape = new BRepTools_ReShape;
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aDone;
  BRep_Builder                                           aB;
  for (int i = 1; i <= anEF.Extent(); i++)
  {
    const TopoDS_Edge&            anE = TopoDS::Edge(anEF.FindKey(i));
    double                        aF, aL;
    const occ::handle<Geom_Curve> aC = BRep_Tool::Curve(anE, aF, aL);
    if (aC.IsNull() || anEF(i).Extent() != 2 || BRep_Tool::Degenerated(anE))
    {
      continue;
    }
    const gp_Pnt aMid   = aC->Value(0.5 * (aF + aL));
    bool         isSeam = false;
    for (NCollection_List<gp_Pnt>::Iterator it(aSeams); it.More() && !isSeam; it.Next())
    {
      isSeam = aMid.Distance(it.Value()) <= std::max(theTol, BRep_Tool::Tolerance(anE));
    }
    const TopoDS_Face aF1 = TopoDS::Face(anEF(i).First());
    const TopoDS_Face aF2 = TopoDS::Face(anEF(i).Last());
    if (!isSeam || aF1.IsSame(aF2) || aDone.Contains(aF1) || aDone.Contains(aF2)
        || aF1.Orientation() != aF2.Orientation())
    {
      continue;
    }
    // the two on one surface: a point of the other's on the one's
    {
      const BRepAdaptor_Surface aS2(aF2);
      const gp_Pnt aP2 = aS2.Value(0.5 * (aS2.FirstUParameter() + aS2.LastUParameter()),
                                   0.5 * (aS2.FirstVParameter() + aS2.LastVParameter()));
      GeomAPI_ProjectPointOnSurf aPr(aP2, BRep_Tool::Surface(aF1));
      if (!aPr.IsDone() || aPr.NbPoints() == 0 || aPr.LowerDistance() > theTol)
      {
        continue;
      }
    }
    // each face's wire through the edge, in order from the edge on
    auto aLoopAfter = [&anE](const TopoDS_Face&                  theF,
                             TopoDS_Wire&                        theW,
                             NCollection_Sequence<TopoDS_Shape>& theEdges) {
      for (TopExp_Explorer exW(theF, TopAbs_WIRE); exW.More(); exW.Next())
      {
        NCollection_Sequence<TopoDS_Shape> anAll;
        int                                anAt = 0;
        for (BRepTools_WireExplorer wex(TopoDS::Wire(exW.Current()), theF); wex.More(); wex.Next())
        {
          anAll.Append(wex.Current());
          if (wex.Current().IsSame(anE))
          {
            anAt = anAll.Length();
          }
        }
        if (anAt != 0)
        {
          theW = TopoDS::Wire(exW.Current());
          for (int k = 1; k < anAll.Length(); k++)
          {
            theEdges.Append(anAll((anAt + k - 1) % anAll.Length() + 1));
          }
          return true;
        }
      }
      return false;
    };
    TopoDS_Wire                        aW1, aW2;
    NCollection_Sequence<TopoDS_Shape> anE1, anE2;
    if (!aLoopAfter(aF1, aW1, anE1) || !aLoopAfter(aF2, aW2, anE2))
    {
      continue;
    }
    TopoDS_Face aNew = TopoDS::Face(aF1.EmptyCopied());
    TopoDS_Wire aW;
    aB.MakeWire(aW);
    for (const NCollection_Sequence<TopoDS_Shape>* aSeq : {&anE1, &anE2})
    {
      for (NCollection_Sequence<TopoDS_Shape>::Iterator it(*aSeq); it.More(); it.Next())
      {
        aB.Add(aW, it.Value());
      }
    }
    aW.Closed(true);
    aB.Add(aNew, aW);
    for (const TopoDS_Face* aFace : {&aF1, &aF2})
    {
      for (TopExp_Explorer exW(*aFace, TopAbs_WIRE); exW.More(); exW.Next())
      {
        if (!exW.Current().IsSame(aFace == &aF1 ? aW1 : aW2))
        {
          aB.Add(aNew, exW.Current());
        }
      }
    }
    // the other's edges on the one's surface: its curves there
    bool isOk = true;
    for (TopExp_Explorer exE(aF2, TopAbs_EDGE); exE.More() && isOk; exE.Next())
    {
      const TopoDS_Edge& anE2 = TopoDS::Edge(exE.Current());
      if (anE2.IsSame(anE) || BRep_Tool::IsClosed(anE2, aF2))
      {
        isOk = !BRep_Tool::IsClosed(anE2, aF2);
        continue;
      }
      double                          aPF, aPL;
      const occ::handle<Geom2d_Curve> aPc = BRep_Tool::CurveOnSurface(anE2, aF2, aPF, aPL);
      if (aPc.IsNull())
      {
        isOk = false;
        continue;
      }
      aB.UpdateEdge(anE2, aPc, aNew, BRep_Tool::Tolerance(anE2));
    }
    if (!isOk || !BRepCheck_Analyzer(aNew).IsValid())
    {
      continue;
    }
    aReShape->Replace(aF1, aNew);
    aReShape->Remove(aF2);
    aDone.Add(aF1);
    aDone.Add(aF2);
  }
  if (aDone.IsEmpty())
  {
    return occ::handle<BRepTools_History>();
  }
  const TopoDS_Shape aMerged = aReShape->Apply(theResult);
  if (aMerged.IsNull() || !BRepCheck_Analyzer(aMerged).IsValid())
  {
    return occ::handle<BRepTools_History>();
  }
  theResult = aMerged;
  return aReShape->History();
}

//=================================================================================================
// A closed side whose seam ends at the vertex where a fillet was pinched
// (PerformPinchEnds: a round plate's rim, its seam through the point where
// a boss inside touches it) and the side carried up over the gap, a face of
// its own on a copy of the side's surface between the cuts and the side's
// edge either side of the vertex: made one, the seam carried up from the
// vertex to the point where the cuts meet, those two pieces of the edge and
// the vertex gone. The pair of points (vertex, point) of each are in
// ChFi3d_SideSeams, which is cleared. Returns the history of the change,
// null when there was none.

static occ::handle<BRepTools_History> ChFi3d_MergeSideSeams(TopoDS_Shape& theResult,
                                                            const double  theTol)
{
  NCollection_List<std::pair<gp_Pnt, gp_Pnt>> aPairs = ChFi3d_SideSeams();
  ChFi3d_SideSeams().Clear();
  if (aPairs.IsEmpty())
  {
    return occ::handle<BRepTools_History>();
  }
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    anEF, aVE;
  TopExp::MapShapesAndAncestors(theResult, TopAbs_EDGE, TopAbs_FACE, anEF);
  TopExp::MapShapesAndAncestors(theResult, TopAbs_VERTEX, TopAbs_EDGE, aVE);
  auto aVertexAt = [&aVE, theTol](const gp_Pnt& theP) {
    for (int i = 1; i <= aVE.Extent(); i++)
    {
      const TopoDS_Vertex& aV = TopoDS::Vertex(aVE.FindKey(i));
      if (BRep_Tool::Pnt(aV).Distance(theP) <= std::max(theTol, BRep_Tool::Tolerance(aV)))
      {
        return aV;
      }
    }
    return TopoDS_Vertex();
  };
  // the oriented edges of the wire of <theF> holding <theE>, in order
  auto aLoop = [](const TopoDS_Face&                  theF,
                  const TopoDS_Shape&                 theE,
                  TopoDS_Wire&                        theW,
                  NCollection_Sequence<TopoDS_Shape>& theEdges) {
    for (TopExp_Explorer exW(theF, TopAbs_WIRE); exW.More(); exW.Next())
    {
      NCollection_Sequence<TopoDS_Shape> anAll;
      bool                               isIn = false;
      for (BRepTools_WireExplorer wex(TopoDS::Wire(exW.Current()), theF); wex.More(); wex.Next())
      {
        anAll.Append(wex.Current());
        isIn = isIn || wex.Current().IsSame(theE);
      }
      if (isIn)
      {
        theW     = TopoDS::Wire(exW.Current());
        theEdges = anAll;
        return true;
      }
    }
    return false;
  };
  // the uv of <theV> on the edge <theE> of <theF>
  auto anUV = [](const TopoDS_Shape& theE, const TopoDS_Vertex& theV, const TopoDS_Face& theF) {
    double                          aF = 0., aL = 0.;
    const occ::handle<Geom2d_Curve> aPc =
      BRep_Tool::CurveOnSurface(TopoDS::Edge(theE), theF, aF, aL);
    return aPc.IsNull() ? gp_Pnt2d(Precision::Infinite(), Precision::Infinite())
                        : aPc->Value(BRep_Tool::Parameter(theV, TopoDS::Edge(theE)));
  };
  occ::handle<BRepTools_ReShape>                         aReShape = new BRepTools_ReShape;
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aDone;
  BRep_Builder                                           aB;
  for (NCollection_List<std::pair<gp_Pnt, gp_Pnt>>::Iterator itP(aPairs); itP.More(); itP.Next())
  {
    const TopoDS_Vertex aV  = aVertexAt(itP.Value().first);
    const TopoDS_Vertex aVP = aVertexAt(itP.Value().second);
    if (aV.IsNull() || aVP.IsNull() || aV.IsSame(aVP))
    {
      continue;
    }
    // the seam ending at the vertex, the side, and the face carried up:
    // the one beside the side across the two edges of the vertex
    TopoDS_Edge aSeam;
    TopoDS_Face aF1, aFn;
    for (NCollection_List<TopoDS_Shape>::Iterator itE(aVE.FindFromKey(aV)); itE.More(); itE.Next())
    {
      const TopoDS_Edge& anE = TopoDS::Edge(itE.Value());
      if (!anEF.Contains(anE))
      {
        continue;
      }
      for (NCollection_List<TopoDS_Shape>::Iterator itF(anEF.FindFromKey(anE)); itF.More();
           itF.Next())
      {
        if (BRep_Tool::IsClosed(anE, TopoDS::Face(itF.Value())))
        {
          aSeam = anE;
          aF1   = TopoDS::Face(itF.Value());
        }
      }
    }
    if (aSeam.IsNull() || aDone.Contains(aF1))
    {
      continue;
    }
    NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aShared;
    for (NCollection_List<TopoDS_Shape>::Iterator itE(aVE.FindFromKey(aV)); itE.More(); itE.Next())
    {
      if (itE.Value().IsSame(aSeam) || !anEF.Contains(itE.Value()))
      {
        continue;
      }
      const NCollection_List<TopoDS_Shape>& aLF = anEF.FindFromKey(itE.Value());
      if (aLF.Extent() != 2 || (!aLF.First().IsSame(aF1) && !aLF.Last().IsSame(aF1)))
      {
        aShared.Clear();
        break;
      }
      const TopoDS_Face aFo = TopoDS::Face(aLF.First().IsSame(aF1) ? aLF.Last() : aLF.First());
      if (!aFn.IsNull() && !aFn.IsSame(aFo))
      {
        aShared.Clear();
        break;
      }
      aFn = aFo;
      aShared.Add(itE.Value());
    }
    if (aShared.Extent() != 2 || aFn.IsNull() || aDone.Contains(aFn)
        || aFn.Orientation() != aF1.Orientation())
    {
      continue;
    }
    const occ::handle<Geom_Surface> aS1 = BRep_Tool::Surface(aF1);
    if (!aS1->IsUPeriodic())
    {
      continue;
    }
    const double                       aPer = aS1->UPeriod();
    TopoDS_Wire                        aW1, aWn;
    NCollection_Sequence<TopoDS_Shape> aL1, aLn;
    if (!aLoop(aF1, aSeam, aW1, aL1)
        || !aLoop(aFn,
                  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>::Iterator(aShared).Value(),
                  aWn,
                  aLn))
    {
      continue;
    }
    // the face carried up's other edges, from one shared edge round to the
    // other, cut in two at the point: <aQ[0]> to it, <aQ[1]> from it
    NCollection_Sequence<TopoDS_Shape> aQ[2];
    {
      const int aNb = aLn.Length();
      int       i0  = 0;
      for (int i = 1; i <= aNb && i0 == 0; i++)
      {
        if (aShared.Contains(aLn(i)) && !aShared.Contains(aLn(i % aNb + 1)))
        {
          i0 = i % aNb + 1;
        }
      }
      int k = 0;
      for (int i = 0; i0 != 0 && i < aNb; i++)
      {
        const TopoDS_Shape& anE = aLn((i0 - 1 + i) % aNb + 1);
        if (aShared.Contains(anE))
        {
          break;
        }
        aQ[k].Append(anE);
        if (k == 0 && TopExp::LastVertex(TopoDS::Edge(anE), true).IsSame(aVP))
        {
          k = 1;
        }
      }
      if (aQ[0].IsEmpty() || aQ[1].IsEmpty() || aQ[0].Length() + aQ[1].Length() + 2 != aNb)
      {
        continue;
      }
    }
    // each half moved by periods onto the side's: where it meets the shared
    // edge it takes the place of
    double aShift[2] = {0., 0.};
    bool   isOk      = true;
    for (int k = 0; k < 2 && isOk; k++)
    {
      const TopoDS_Edge   aQE = TopoDS::Edge(k == 0 ? aQ[0].First() : aQ[1].Last());
      const TopoDS_Vertex aA =
        k == 0 ? TopExp::FirstVertex(aQE, true) : TopExp::LastVertex(aQE, true);
      TopoDS_Shape aSh;
      for (NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>::Iterator itS(aShared);
           itS.More();
           itS.Next())
      {
        TopoDS_Vertex aS0, aS1;
        TopExp::Vertices(TopoDS::Edge(itS.Value()), aS0, aS1);
        if (aS0.IsSame(aA) || aS1.IsSame(aA))
        {
          aSh = itS.Value();
        }
      }
      if (aSh.IsNull())
      {
        isOk = false;
        break;
      }
      const gp_Pnt2d aUV1 = anUV(aSh, aA, aF1), aUVn = anUV(aQE, aA, aFn);
      if (Precision::IsInfinite(aUV1.X()) || Precision::IsInfinite(aUVn.X()))
      {
        isOk = false;
        break;
      }
      aShift[k] = aPer * std::round((aUV1.X() - aUVn.X()) / aPer);
      isOk      = aS1->Value(aUVn.X() + aShift[k], aUVn.Y())
               .Transformed(aF1.Location().Transformation())
               .Distance(BRep_Tool::Pnt(aA))
             <= std::max(theTol, BRep_Tool::Tolerance(aA));
    }
    if (!isOk)
    {
      continue;
    }
    // the seam carried up from the vertex to the point, along its line
    TopoDS_Vertex aSV0, aSV1;
    TopExp::Vertices(aSeam, aSV0, aSV1);
    double                        aF = 0., aL = 0.;
    const occ::handle<Geom_Curve> aCS = BRep_Tool::Curve(aSeam, aF, aL);
    if (aCS.IsNull() || aSV0.IsSame(aSV1))
    {
      continue;
    }
    GeomAPI_ProjectPointOnCurve aPr(BRep_Tool::Pnt(aVP), aCS);
    if (aPr.NbPoints() == 0 || aPr.LowerDistance() > std::max(theTol, BRep_Tool::Tolerance(aVP)))
    {
      continue;
    }
    const bool   isVLast = aSV1.IsSame(aV);
    const double aTP     = aPr.LowerDistanceParameter();
    if (isVLast ? aTP <= aL : aTP >= aF)
    {
      continue;
    }
    TopoDS_Edge aNewSeam = TopoDS::Edge(aSeam.Oriented(TopAbs_FORWARD).EmptyCopied());
    aB.Add(aNewSeam, (isVLast ? aSV0 : aVP).Oriented(TopAbs_FORWARD));
    aB.Add(aNewSeam, (isVLast ? aVP : aSV1).Oriented(TopAbs_REVERSED));
    aB.Range(aNewSeam, isVLast ? aF : aTP, isVLast ? aTP : aL);
    // the point off the line by as much as the cuts' meeting allows
    aB.UpdateVertex(aVP, std::max(BRep_Tool::Tolerance(aVP), 1.001 * aPr.LowerDistance()));
    // the side's wire with the seam carried up and the two shared edges
    // replaced by the halves: the one from the vertex by the half from the
    // point, the one to it by the half to it
    TopoDS_Wire aW;
    aB.MakeWire(aW);
    for (NCollection_Sequence<TopoDS_Shape>::Iterator it(aL1); it.More(); it.Next())
    {
      const TopoDS_Shape& anE = it.Value();
      if (anE.IsSame(aSeam))
      {
        aB.Add(aW, aNewSeam.Oriented(anE.Orientation()));
      }
      else if (aShared.Contains(anE))
      {
        const int k = TopExp::FirstVertex(TopoDS::Edge(anE), true).IsSame(aV) ? 1 : 0;
        for (NCollection_Sequence<TopoDS_Shape>::Iterator itQ(aQ[k]); itQ.More(); itQ.Next())
        {
          aB.Add(aW, itQ.Value());
        }
      }
      else
      {
        aB.Add(aW, anE);
      }
    }
    aW.Closed(true);
    TopoDS_Face aNew = TopoDS::Face(aF1.EmptyCopied());
    aB.Add(aNew, aW);
    for (TopExp_Explorer exW(aF1, TopAbs_WIRE); exW.More(); exW.Next())
    {
      if (!exW.Current().IsSame(aW1))
      {
        aB.Add(aNew, exW.Current());
      }
    }
    for (TopExp_Explorer exW(aFn, TopAbs_WIRE); exW.More(); exW.Next())
    {
      if (!exW.Current().IsSame(aWn))
      {
        isOk = false;
      }
    }
    // the halves' edges on the side's surface: their curves there
    for (int k = 0; k < 2 && isOk; k++)
    {
      for (NCollection_Sequence<TopoDS_Shape>::Iterator itQ(aQ[k]); itQ.More() && isOk; itQ.Next())
      {
        const TopoDS_Edge&              anE = TopoDS::Edge(itQ.Value());
        double                          aPF = 0., aPL = 0.;
        const occ::handle<Geom2d_Curve> aPc = BRep_Tool::CurveOnSurface(anE, aFn, aPF, aPL);
        if (aPc.IsNull() || BRep_Tool::IsClosed(anE, aFn))
        {
          isOk = false;
          break;
        }
        aB.UpdateEdge(anE,
                      occ::down_cast<Geom2d_Curve>(aPc->Translated(gp_Vec2d(aShift[k], 0.))),
                      aNew,
                      BRep_Tool::Tolerance(anE));
      }
    }
    if (!isOk || !BRepCheck_Analyzer(aNew).IsValid())
    {
      continue;
    }
    aReShape->Replace(aF1, aNew);
    aReShape->Remove(aFn);
    aDone.Add(aF1);
    aDone.Add(aFn);
  }
  if (aDone.IsEmpty())
  {
    return occ::handle<BRepTools_History>();
  }
  const TopoDS_Shape aMerged = aReShape->Apply(theResult);
  if (aMerged.IsNull() || !BRepCheck_Analyzer(aMerged).IsValid())
  {
    return occ::handle<BRepTools_History>();
  }
  theResult = aMerged;
  return aReShape->History();
}

//=================================================================================================
// A face no wider than <theTol> beside a face the fillet made: two long
// edges, within <theTol> of each other, and between their ends short ones
// or none -- what a fillet as wide as the face leaves of it with its radius
// a hair smaller (ChFi3d_SnapSlivers; a wall's top 2 wide, its outer edge at
// 2; a full round, the two edges at half the width). The face goes, its
// short edges with it, one long edge takes the other's place in the face
// beyond, a curve on that face's surface, and the ends are made one. The
// edge kept is the input's where one is, else the one whose face beyond is
// on a surface of the input's: the fillet's face reaches over to it, so a
// run of slivers (#273's ledge, a fillet as wide as a step 1 deep below the
// top: three faces 1e-5 tall around it) all go the same way, to the top.
// Nothing changes when the volume would by more than 1e-6 of itself.
// Returns the history of the change, null when there was none.

static occ::handle<BRepTools_History> ChFi3d_CollapseSlivers(TopoDS_Shape&       theResult,
                                                             const TopoDS_Shape& theInput,
                                                             const double        theTol)
{
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> anInputE;
  TopExp::MapShapes(theInput, TopAbs_EDGE, anInputE);
  ChFi3d_AddThawed(anInputE, theResult, TopAbs_EDGE);
  NCollection_Map<occ::handle<Standard_Transient>> anInputS;
  for (TopExp_Explorer exF(theInput, TopAbs_FACE); exF.More(); exF.Next())
  {
    TopLoc_Location aLoc;
    anInputS.Add(BRep_Tool::Surface(TopoDS::Face(exF.Current()), aLoc));
  }
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    anEF;
  TopExp::MapShapesAndAncestors(theResult, TopAbs_EDGE, TopAbs_FACE, anEF);
  occ::handle<BRepTools_ReShape>                         aReShape = new BRepTools_ReShape;
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aDone, aDropped;
  BRep_Builder                                           aB;
  // the vertices made one, and where on which edges the one dropped was
  NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher> aVDropToKeep;
  // a frozen vertex kept: its copy raised, in the result once it is made
  NCollection_DataMap<TopoDS_Shape, double, TopTools_ShapeMapHasher> aRaised;
  auto aLength = [](const TopoDS_Edge& theE) {
    return BRep_Tool::IsGeometric(theE) ? GCPnts_AbscissaPoint::Length(BRepAdaptor_Curve(theE))
                                        : 0.;
  };
  // the face beyond <theE> from <theS> on a surface of the input's
  auto isBeyondInput = [&](const TopoDS_Shape& theE, const TopoDS_Shape& theS) {
    for (NCollection_List<TopoDS_Shape>::Iterator it(anEF.FindFromKey(theE)); it.More(); it.Next())
    {
      TopLoc_Location aLoc;
      if (!it.Value().IsSame(theS)
          && anInputS.Contains(BRep_Tool::Surface(TopoDS::Face(it.Value()), aLoc)))
      {
        return true;
      }
    }
    return false;
  };
  auto aVolume = [](const TopoDS_Shape& theShape) {
    GProp_GProps aProps;
    BRepGProp::VolumeProperties(theShape, aProps);
    return aProps.Mass();
  };
  for (TopExp_Explorer exF(theResult, TopAbs_FACE); exF.More(); exF.Next())
  {
    const TopoDS_Face aS = TopoDS::Face(exF.Current());
    if (aDone.Contains(aS))
    {
      continue;
    }
    int aNbW = 0;
    for (TopExp_Explorer exW(aS, TopAbs_WIRE); exW.More(); exW.Next())
    {
      aNbW++;
    }
    if (aNbW != 1)
    {
      continue;
    }
    // two long edges, the rest short
    NCollection_List<TopoDS_Shape> aLong, aShort;
    bool                           isOk = true;
    for (TopExp_Explorer exE(aS, TopAbs_EDGE); exE.More() && isOk; exE.Next())
    {
      const TopoDS_Edge& anE = TopoDS::Edge(exE.Current());
      if (BRep_Tool::Degenerated(anE) || BRep_Tool::IsClosed(anE, aS))
      {
        isOk = false;
      }
      else if (aLength(anE) <= theTol)
      {
        aShort.Append(anE);
      }
      else
      {
        aLong.Append(anE);
      }
    }
    if (!isOk || aLong.Size() != 2 || aShort.Size() > 2)
    {
      continue;
    }
    TopoDS_Edge aKeep = TopoDS::Edge(aLong.First()), aDrop = TopoDS::Edge(aLong.Last());
    if (anInputE.Contains(aDrop) && !anInputE.Contains(aKeep))
    {
      std::swap(aKeep, aDrop);
    }
    else if (!anInputE.Contains(aKeep) && isBeyondInput(aDrop, aS) && !isBeyondInput(aKeep, aS))
    {
      std::swap(aKeep, aDrop);
    }
    if (anInputE.Contains(aDrop) || !anEF.Contains(aDrop) || anEF.FindFromKey(aDrop).Extent() != 2)
    {
      continue;
    }
    // within the tolerance of each other all along
    double                        aKF, aKL, aDF, aDL;
    const occ::handle<Geom_Curve> aCK = BRep_Tool::Curve(aKeep, aKF, aKL);
    const occ::handle<Geom_Curve> aCD = BRep_Tool::Curve(aDrop, aDF, aDL);
    if (aCK.IsNull() || aCD.IsNull())
    {
      continue;
    }
    double aGap = 0.;
    for (int i = 0; i <= 10 && aGap <= theTol; i++)
    {
      GeomAPI_ProjectPointOnCurve aPr(aCD->Value(aDF + (aDL - aDF) * i / 10.), aCK, aKF, aKL);
      aGap = aPr.NbPoints() > 0 ? std::max(aGap, aPr.LowerDistance()) : RealLast();
    }
    if (aGap > theTol)
    {
      continue;
    }
    // the face beyond the edge dropped, the keep's curve on its surface
    TopoDS_Face aT;
    for (NCollection_List<TopoDS_Shape>::Iterator it(anEF.FindFromKey(aDrop)); it.More(); it.Next())
    {
      if (!it.Value().IsSame(aS))
      {
        aT = TopoDS::Face(it.Value());
      }
    }
    if (aT.IsNull() || aDone.Contains(aT))
    {
      continue;
    }
    TopLoc_Location                 aLoc;
    const occ::handle<Geom_Surface> aST    = BRep_Tool::Surface(aT, aLoc);
    double                          aTolPj = theTol;
    const occ::handle<Geom_Curve>   aCKL =
      occ::down_cast<Geom_Curve>(aCK->Transformed(aLoc.Transformation().Inverted()));
    occ::handle<Geom2d_Curve> aPc = GeomProjLib::Curve2d(aCKL, aKF, aKL, aST, aTolPj);
    if (aPc.IsNull())
    {
      continue;
    }
    // on a periodic surface, the period the dropped edge's curve is in
    {
      double                          aPF, aPL;
      const occ::handle<Geom2d_Curve> aPcD = BRep_Tool::CurveOnSurface(aDrop, aT, aPF, aPL);
      if (aPcD.IsNull())
      {
        continue;
      }
      const gp_Pnt2d aMidD = aPcD->Value(0.5 * (aPF + aPL));
      const gp_Pnt2d aMidK = aPc->Value(0.5 * (aKF + aKL));
      gp_Vec2d       aShift(0., 0.);
      if (aST->IsUPeriodic())
      {
        const double aT0 = aST->UPeriod();
        aShift.SetX(aT0 * std::round((aMidD.X() - aMidK.X()) / aT0));
      }
      if (aST->IsVPeriodic())
      {
        const double aT0 = aST->VPeriod();
        aShift.SetY(aT0 * std::round((aMidD.Y() - aMidK.Y()) / aT0));
      }
      if (aShift.Magnitude() > 0.)
      {
        aPc = occ::down_cast<Geom2d_Curve>(aPc->Translated(aShift));
      }
    }
    double aDev = 0.;
    for (int i = 0; i <= 20; i++)
    {
      const double   t   = aKF + (aKL - aKF) * i / 20.;
      const gp_Pnt2d aUV = aPc->Value(t);
      aDev               = std::max(aDev, aST->Value(aUV.X(), aUV.Y()).Distance(aCKL->Value(t)));
    }
    if (aDev > 2. * theTol)
    {
      continue;
    }
    // the ends: each of the dropped edge's to the keep's nearest
    TopoDS_Vertex aK1, aK2, aD1, aD2;
    TopExp::Vertices(aKeep, aK1, aK2);
    TopExp::Vertices(aDrop, aD1, aD2);
    const gp_Pnt aPK1 = BRep_Tool::Pnt(aK1), aPD1 = BRep_Tool::Pnt(aD1), aPD2 = BRep_Tool::Pnt(aD2);
    const bool   isSameWay     = aPD1.Distance(aPK1) <= aPD2.Distance(aPK1);
    const TopoDS_Vertex aKfor1 = isSameWay ? aK1 : aK2, aKfor2 = isSameWay ? aK2 : aK1;
    if (BRep_Tool::Pnt(aD1).Distance(BRep_Tool::Pnt(aKfor1)) > 2. * theTol
        || BRep_Tool::Pnt(aD2).Distance(BRep_Tool::Pnt(aKfor2)) > 2. * theTol)
    {
      continue;
    }
    aB.UpdateEdge(aKeep, aPc, aT, std::max(BRep_Tool::Tolerance(aKeep), aDev));
    for (int k = 0; k < 2; k++)
    {
      const TopoDS_Vertex& aV  = k == 0 ? aD1 : aD2;
      const TopoDS_Vertex& aVK = k == 0 ? aKfor1 : aKfor2;
      if (aV.IsSame(aVK))
      {
        continue;
      }
      const double aD = BRep_Tool::Pnt(aV).Distance(BRep_Tool::Pnt(aVK)) + BRep_Tool::Tolerance(aV);
      ChFi3d_RaiseKept(aVK, aD, aRaised);
      aVDropToKeep.Bind(aV, aVK);
      aReShape->Replace(aV.Oriented(TopAbs_FORWARD), aVK.Oriented(TopAbs_FORWARD));
    }
    // the keep in the dropped edge's place, the same way along
    TopoDS_Edge aKeepAs = aKeep;
    aKeepAs.Orientation(isSameWay ? TopAbs_FORWARD : TopAbs_REVERSED);
    TopoDS_Edge aDropF = aDrop;
    aDropF.Orientation(TopAbs_FORWARD);
    aReShape->Replace(aDropF, aKeepAs);
    for (NCollection_List<TopoDS_Shape>::Iterator it(aShort); it.More(); it.Next())
    {
      aReShape->Remove(it.Value());
    }
    aReShape->Remove(aS);
    aDone.Add(aS);
    aDone.Add(aT);
    aDropped.Add(aDrop);
  }
  if (aDone.IsEmpty())
  {
    return occ::handle<BRepTools_History>();
  }
  // the edges at a vertex made one take the one kept at the same parameter
  NCollection_List<TopoDS_Shape>                                anEdges, aKeeps;
  NCollection_List<double>                                      aParams;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aResultE;
  TopExp::MapShapes(theResult, TopAbs_EDGE, aResultE);
  for (int i = 1; i <= aResultE.Extent(); i++)
  {
    const TopoDS_Edge& anE = TopoDS::Edge(aResultE(i));
    // an edge dropped: the one kept in its place has its own parameters
    if (aDropped.Contains(anE))
    {
      continue;
    }
    for (TopoDS_Iterator itV(anE); itV.More(); itV.Next())
    {
      if (aVDropToKeep.IsBound(itV.Value()))
      {
        anEdges.Append(anE);
        aKeeps.Append(aVDropToKeep(itV.Value()));
        aParams.Append(BRep_Tool::Parameter(TopoDS::Vertex(itV.Value()), anE));
      }
    }
  }
  const TopoDS_Shape                       aCollapsed = aReShape->Apply(theResult);
  NCollection_List<TopoDS_Shape>::Iterator itK(aKeeps);
  NCollection_List<double>::Iterator       itP(aParams);
  for (NCollection_List<TopoDS_Shape>::Iterator itE(anEdges); itE.More();
       itE.Next(), itK.Next(), itP.Next())
  {
    const TopoDS_Shape aNewE = aReShape->Value(itE.Value());
    if (aNewE.IsNull() || aNewE.ShapeType() != TopAbs_EDGE)
    {
      continue;
    }
    const TopoDS_Vertex& aVK = TopoDS::Vertex(itK.Value());
    aB.UpdateVertex(aVK, itP.Value(), TopoDS::Edge(aNewE), BRep_Tool::Tolerance(aVK));
  }
  // the frozen vertices kept, their copies raised in their places; refused,
  // in the result as it was all the same, as a vertex raised in place stays
  // raised -- the same copies, which the result shares containers with
  const NCollection_List<TopoDS_Shape> aKeptCopies = ChFi3d_ThawedCopies(aRaised);
  BRepLib_PutThawed(aCollapsed, aKeptCopies);
  if (aCollapsed.IsNull() || !BRepCheck_Analyzer(aCollapsed).IsValid())
  {
    BRepLib_PutThawed(theResult, aKeptCopies);
    return occ::handle<BRepTools_History>();
  }
  const double aVolBefore = aVolume(theResult);
  if (std::abs(aVolume(aCollapsed) - aVolBefore) > 1e-6 * std::abs(aVolBefore))
  {
    BRepLib_PutThawed(theResult, aKeptCopies);
    return occ::handle<BRepTools_History>();
  }
  theResult = aCollapsed;
  return aReShape->History();
}

//=================================================================================================
// Whether <theResult> has a face no wider than <theTol> (twice its area over
// its outline) that is not one of <theInput>'s, or a vertex not <theInput>'s
// that touches another vertex or an edge it does not bound, within their
// tolerances (as BOP's check finds them): a sliver the radii a hair smaller
// left and ChFi3d_CollapseSlivers did not take out. #474 `Fillet003`'s E52
// at 2, where the coil's band is 2 thick, kept a face of area 9e-12; a
// slot's wall 4 wide at the corner's edge filleted at 4, the fillet's line
// 1e-5 short of the wall's far edge in a face that runs on below the slot,
// a vertex of tolerance 1e-4 that far from that edge.

static bool ChFi3d_HasSliverLeft(const TopoDS_Shape& theResult,
                                 const TopoDS_Shape& theInput,
                                 const double        theTol)
{
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> anInputF;
  TopExp::MapShapes(theInput, TopAbs_FACE, anInputF);
  ChFi3d_AddThawed(anInputF, theResult, TopAbs_FACE);
  for (TopExp_Explorer exF(theResult, TopAbs_FACE); exF.More(); exF.Next())
  {
    if (anInputF.Contains(exF.Current()))
    {
      continue;
    }
    GProp_GProps aProps;
    BRepGProp::SurfaceProperties(exF.Current(), aProps);
    double anOutline = 0.;
    for (TopExp_Explorer exE(exF.Current(), TopAbs_EDGE); exE.More(); exE.Next())
    {
      const TopoDS_Edge& anE = TopoDS::Edge(exE.Current());
      if (BRep_Tool::IsGeometric(anE))
      {
        anOutline += GCPnts_AbscissaPoint::Length(BRepAdaptor_Curve(anE));
      }
    }
    if (anOutline > 0. && 2. * std::abs(aProps.Mass()) <= theTol * anOutline)
    {
      return true;
    }
  }
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> anInputV, aV, anE;
  TopExp::MapShapes(theInput, TopAbs_VERTEX, anInputV);
  ChFi3d_AddThawed(anInputV, theResult, TopAbs_VERTEX);
  TopExp::MapShapes(theResult, TopAbs_VERTEX, aV);
  TopExp::MapShapes(theResult, TopAbs_EDGE, anE);
  NCollection_Array1<Bnd_Box> anEBox(1, std::max(1, anE.Extent()));
  double                      aMaxVTol = 0.;
  for (int i = 1; i <= aV.Extent(); i++)
  {
    aMaxVTol = std::max(aMaxVTol, BRep_Tool::Tolerance(TopoDS::Vertex(aV(i))));
  }
  for (int i = 1; i <= anE.Extent(); i++)
  {
    BRepBndLib::Add(anE(i), anEBox(i));
    anEBox(i).Enlarge(aMaxVTol);
  }
  for (int i = 1; i <= aV.Extent(); i++)
  {
    if (anInputV.Contains(aV(i)))
    {
      continue;
    }
    const TopoDS_Vertex& aVi  = TopoDS::Vertex(aV(i));
    const gp_Pnt         aP   = BRep_Tool::Pnt(aVi);
    const double         aTol = BRep_Tool::Tolerance(aVi);
    for (int j = 1; j <= aV.Extent(); j++)
    {
      const TopoDS_Vertex& aVj = TopoDS::Vertex(aV(j));
      if (j != i && aP.Distance(BRep_Tool::Pnt(aVj)) <= aTol + BRep_Tool::Tolerance(aVj))
      {
        return true;
      }
    }
    for (int k = 1; k <= anE.Extent(); k++)
    {
      if (anEBox(k).IsOut(aP))
      {
        continue;
      }
      bool isOn = false;
      for (TopoDS_Iterator itV(anE(k)); itV.More() && !isOn; itV.Next())
      {
        isOn = itV.Value().IsSame(aVi);
      }
      if (isOn)
      {
        continue;
      }
      BRepExtrema_DistShapeShape aDist(aVi, anE(k));
      if (aDist.IsDone() && aDist.NbSolution() > 0
          && aDist.Value() <= aTol + BRep_Tool::Tolerance(TopoDS::Edge(anE(k))))
      {
        return true;
      }
    }
  }
  return false;
}

//=================================================================================================
// An edge the fillet made that is no longer than twice <theTol>, the
// tolerance the fillet's surfaces are made to -- its ends, each known to
// that, are one point -- and that its two vertices cover: nothing of it
// lies between them (BOPTools_AlgoTools::IsMicroEdge, the question BOP's
// argument check asks). A cut across a face of no width makes one -- #474's
// Fillet003 keeps a face 1.4e-7 wide beside a corner of its hexagon, and
// every fillet ending there crossed it with an edge of that length -- and
// so does a fillet thinner at its end than its own tolerance (#876's post
// drafted by a degree, the edge at its foot: 1.1e-6 long, and at a larger
// radius 1.0e-4 between ends 5.1e-5 apart; #273's handle, two degrees from
// the wall it stands on: 1.5e-5, closed on one vertex). BRepCheck lets such
// an edge be; BOP's check reports it too small and its vertices
// self-intersecting, and the next boolean is refused its argument. The edge
// is dropped and its vertices made one, the input's when one of them is. An
// edge between two vertices of the input is left alone, and nothing is done
// where a wire would be left without an edge. A longer edge that its
// vertices cover all the same is the mark of a loose result, not of an edge
// of no length: left as it is (#273's fillets at a tolerance of 0.045, an
// edge of 0.029 -- dropping it moved the volume by 1.9). Returns the
// history of the change, null when there was none.

static occ::handle<BRepTools_History> ChFi3d_DropMicroEdges(TopoDS_Shape&       theResult,
                                                           const TopoDS_Shape& theInput,
                                                           const double        theTol)
{
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> anInputE, anInputV, anEdges;
  TopExp::MapShapes(theInput, TopAbs_EDGE, anInputE);
  TopExp::MapShapes(theInput, TopAbs_VERTEX, anInputV);
  ChFi3d_AddThawed(anInputE, theResult, TopAbs_EDGE);
  ChFi3d_AddThawed(anInputV, theResult, TopAbs_VERTEX);
  TopExp::MapShapes(theResult, TopAbs_EDGE, anEdges);
  NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher> aDropToKeep;
  // the vertex a vertex has become
  auto aKept = [&aDropToKeep](TopoDS_Shape theV) {
    while (aDropToKeep.IsBound(theV))
    {
      theV = aDropToKeep(theV);
    }
    return theV;
  };
  occ::handle<IntTools_Context>                          aCtx = new IntTools_Context;
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aMicro;
  for (int i = 1; i <= anEdges.Extent(); i++)
  {
    const TopoDS_Edge& anE = TopoDS::Edge(anEdges(i));
    if (anInputE.Contains(anE) || BRep_Tool::Degenerated(anE))
    {
      continue;
    }
    TopoDS_Vertex aV1, aV2;
    TopExp::Vertices(anE, aV1, aV2);
    if (aV1.IsNull() || aV2.IsNull() || !BRep_Tool::IsGeometric(anE)
        || BRep_Tool::Pnt(aV1).Distance(BRep_Tool::Pnt(aV2)) > 2. * theTol
        || GCPnts_AbscissaPoint::Length(BRepAdaptor_Curve(anE)) > 2. * theTol
        || !BOPTools_AlgoTools::IsMicroEdge(anE, aCtx))
    {
      continue;
    }
    // closed on one vertex already: the edge alone goes
    TopoDS_Shape aKeep = aKept(aV1.Oriented(TopAbs_FORWARD));
    TopoDS_Shape aDrop = aKept(aV2.Oriented(TopAbs_FORWARD));
    if (!aKeep.IsSame(aDrop))
    {
      if (anInputV.Contains(aDrop))
      {
        if (anInputV.Contains(aKeep))
        {
          continue;
        }
        std::swap(aKeep, aDrop);
      }
      aDropToKeep.Bind(aDrop, aKeep);
    }
    aMicro.Add(anE);
  }
  if (aMicro.IsEmpty())
  {
    return occ::handle<BRepTools_History>();
  }
  for (TopExp_Explorer exW(theResult, TopAbs_WIRE); exW.More(); exW.Next())
  {
    bool isLeft = false;
    for (TopoDS_Iterator itE(exW.Current()); itE.More() && !isLeft; itE.Next())
    {
      isLeft = !aMicro.Contains(itE.Value());
    }
    if (!isLeft)
    {
      return occ::handle<BRepTools_History>();
    }
  }

  occ::handle<BRepTools_ReShape> aReShape = new BRepTools_ReShape;
  // a frozen vertex kept: its copy raised, in the result once it is made
  NCollection_DataMap<TopoDS_Shape, double, TopTools_ShapeMapHasher> aRaised;
  for (NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>::Iterator it(
         aDropToKeep);
       it.More();
       it.Next())
  {
    const TopoDS_Vertex& aDrop = TopoDS::Vertex(it.Key());
    const TopoDS_Vertex  aKeep = TopoDS::Vertex(aKept(it.Value()));
    const double         aTol =
      BRep_Tool::Pnt(aDrop).Distance(BRep_Tool::Pnt(aKeep)) + BRep_Tool::Tolerance(aDrop);
    ChFi3d_RaiseKept(aKeep, aTol, aRaised);
    aReShape->Replace(aDrop.Oriented(TopAbs_FORWARD), aKeep.Oriented(TopAbs_FORWARD));
  }
  for (NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>::Iterator it(aMicro); it.More();
       it.Next())
  {
    aReShape->Remove(it.Key());
  }
  theResult = aReShape->Apply(theResult);
  ChFi3d_ThawInto(theResult, aRaised);
  return aReShape->History();
}

//=================================================================================================
// The shapes of <theShapes> as <theHistory> has made them, in place.

static void ChFi3d_Remap(NCollection_List<TopoDS_Shape>&       theShapes,
                         const occ::handle<BRepTools_History>& theHistory)
{
  NCollection_List<TopoDS_Shape> anImages;
  for (NCollection_List<TopoDS_Shape>::Iterator it(theShapes); it.More(); it.Next())
  {
    if (theHistory->IsRemoved(it.Value()))
    {
      continue;
    }
    const NCollection_List<TopoDS_Shape>& aModified = theHistory->Modified(it.Value());
    if (aModified.IsEmpty())
    {
      anImages.Append(it.Value());
    }
    for (NCollection_List<TopoDS_Shape>::Iterator itM(aModified); itM.More(); itM.Next())
    {
      // a face split: the compound of its pieces
      if (itM.Value().ShapeType() != TopAbs_COMPOUND)
      {
        anImages.Append(itM.Value());
      }
      for (TopoDS_Iterator itC(itM.Value());
           itM.Value().ShapeType() == TopAbs_COMPOUND && itC.More();
           itC.Next())
      {
        anImages.Append(itC.Value());
      }
    }
  }
  theShapes = anImages;
}

//=================================================================================================

namespace
{
// Sets THE_IN_SETBACK_COMPUTE for as long as it lives.
class SetbackComputeGuard
{
public:
  SetbackComputeGuard() { THE_IN_SETBACK_COMPUTE = true; }

  ~SetbackComputeGuard() { THE_IN_SETBACK_COMPUTE = false; }
};

// False where a stripe ending at <theV> has no point at that end, on either
// face: the corner there has not made its end.
bool StripeEndsHavePoints(const NCollection_List<occ::handle<ChFiDS_Stripe>>& theStripes,
                          const TopoDS_Vertex&                                 theV)
{
  for (NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator it(theStripes); it.More();
       it.Next())
  {
    int aSens = 0;
    ChFi3d_IndexOfSurfData(theV, it.Value(), aSens);
    const bool                       isFirst = aSens == 1;
    const occ::handle<ChFiDS_Spine>& aSp     = it.Value()->Spine();
    if (!aSp.IsNull() && aSp->Status(isFirst) == ChFiDS_FreeBoundary)
    {
      continue;
    }
    if (it.Value()->IndexPoint(isFirst, 1) == 0 || it.Value()->IndexPoint(isFirst, 2) == 0)
    {
      return false;
    }
  }
  return true;
}

// The largest tolerance of an edge of <theShape>.
double MaxEdgeTolerance(const TopoDS_Shape& theShape)
{
  double aTol = 0.;
  for (TopExp_Explorer anExp(theShape, TopAbs_EDGE); anExp.More(); anExp.Next())
  {
    aTol = std::max(aTol, BRep_Tool::Tolerance(TopoDS::Edge(anExp.Current())));
  }
  return aTol;
}

// Whether <theShape> has a face of no area.
bool HasFaceWithoutArea(const TopoDS_Shape& theShape)
{
  for (TopExp_Explorer anExp(theShape, TopAbs_FACE); anExp.More(); anExp.Next())
  {
    GProp_GProps aProps;
    BRepGProp::SurfaceProperties(anExp.Current(), aProps);
    if (std::abs(aProps.Mass()) < Precision::SquareConfusion())
    {
      return true;
    }
  }
  return false;
}

// Whether <theShape> is a result the setback fallback may hand back: no edge
// looser than <theMaxTol>, no face of no area, valid as BRepCheck finds it,
// and valid as it reads back from its own text. A corner's plate whose
// boundary ran over several faces was valid in memory, and not once
// written; a set back corner on #523's cylinder was valid with an edge of
// 0.39 at radius 3; corners of #876's Fillet002 were valid with a face of
// area 2e-17 between two edges, which no mesh covers.
bool IsValidResult(const TopoDS_Shape& theShape, const double theMaxTol)
{
  if (theShape.IsNull() || MaxEdgeTolerance(theShape) > theMaxTol
      || HasFaceWithoutArea(theShape) || !BRepCheck_Analyzer(theShape).IsValid())
  {
    return false;
  }
  std::stringstream aStream;
  BRepTools::Write(theShape, aStream);
  TopoDS_Shape aRead;
  BRep_Builder aBuilder;
  BRepTools::Read(aRead, aStream, aBuilder);
  return !aRead.IsNull() && BRepCheck_Analyzer(aRead).IsValid();
}

// The normal of <theF>, outward, at its point nearest <theP>.
gp_Dir OutwardNormalAt(const TopoDS_Face& theF, const gp_Pnt& theP)
{
  BRepAdaptor_Surface aS(theF);
  GeomAPI_ProjectPointOnSurf aProj(theP, BRep_Tool::Surface(theF));
  double aU = 0., aV = 0.;
  if (aProj.NbPoints() > 0)
  {
    aProj.LowerDistanceParameters(aU, aV);
  }
  BRepLProp_SLProps aProps(aS, aU, aV, 1, Precision::Confusion());
  gp_Dir aN = aProps.IsNormalDefined() ? aProps.Normal() : gp::DZ();
  return theF.Orientation() == TopAbs_REVERSED ? aN.Reversed() : aN;
}

// Whether every fillet of <theStripes> is in <theResult>: the midpoint of
// each edge filleted lies off the result's faces by half the distance a
// fillet of its radius moves it, r (1 / cos(a / 2) - 1) for the angle a
// between its faces' normals there, and by 1e-6 at least. A set back
// corner's computation could come out valid and be the input itself, the
// fillet asked for dropped (radius 0.8 on a right angle, 0.33 when it is
// made); on the near flat edges of #876's drafted walls a fillet moves it
// by thousandths of a thousandth. An edge of a spine of several that lies
// whole within the setback of the spine's end it is at -- or within the
// radius there, which is the least a set back corner takes -- is not
// looked at: it is under the corner's patch, which passes nearer to it than
// a fillet would (a fillet along a wall kept in pieces, the piece at the
// corner shorter than the radius: refused for it, each time, valid). Half
// the radius is the most asked: a fillet that fills a narrow corner moves
// the edge's middle further from its faces than from the others round it --
// #273's pocket tip, 22 deg and 1 deep, filled at 0.3, 1.28 from its fillet
// and 0.5 from the top and the floor.
bool AreFilletsMade(const NCollection_List<occ::handle<ChFiDS_Stripe>>& theStripes,
                    const ChFiDS_Map&                                   theEFMap,
                    const TopoDS_Shape&                                 theResult)
{
  TopoDS_Compound aFaces;
  BRep_Builder    aBuilder;
  aBuilder.MakeCompound(aFaces);
  for (TopExp_Explorer anExp(theResult, TopAbs_FACE); anExp.More(); anExp.Next())
  {
    aBuilder.Add(aFaces, anExp.Current());
  }
  for (NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator it(theStripes); it.More();
       it.Next())
  {
    const occ::handle<ChFiDS_FilSpine> aSp = occ::down_cast<ChFiDS_FilSpine>(it.Value()->Spine());
    if (aSp.IsNull())
    {
      continue;
    }
    // the edges under a set back corner, from either end
    const int                aNbE = aSp->NbEdges();
    NCollection_Array1<bool> isUnder(1, std::max(aNbE, 1));
    isUnder.Init(false);
    for (int k = 0; k < 2 && aNbE > 1 && !aSp->IsClosed(); k++)
    {
      const bool isFirst = k == 0;
      if (aSp->Setback(isFirst) < 0.)
      {
        continue;
      }
      const int    aJ0   = isFirst ? 1 : aNbE;
      const double aTake = std::max(aSp->Setback(isFirst),
                                    aSp->IsConstant(aJ0) ? aSp->Radius(aJ0)
                                                         : aSp->MaxRadFromSeqAndLaws());
      double aLen = 0.;
      for (int n = 0; n < aNbE - 1; n++)
      {
        const int aJ = isFirst ? 1 + n : aNbE - n;
        aLen += GCPnts_AbscissaPoint::Length(BRepAdaptor_Curve(aSp->Edges(aJ)));
        if (aLen > aTake)
        {
          break;
        }
        isUnder.SetValue(aJ, true);
      }
    }
    for (int j = 1; j <= aNbE; j++)
    {
      if (isUnder.Value(j))
      {
        continue;
      }
      const TopoDS_Edge& anE = aSp->Edges(j);
      BRepAdaptor_Curve  aC(anE);
      const gp_Pnt       aMid = aC.Value(0.5 * (aC.FirstParameter() + aC.LastParameter()));
      const double       aR   = aSp->IsConstant(j) ? aSp->Radius(j) : aSp->MaxRadFromSeqAndLaws();
      double             aMove = 0.;
      if (theEFMap.Contains(anE))
      {
        const NCollection_List<TopoDS_Shape>& aFs = theEFMap.FindFromKey(anE);
        if (aFs.Extent() >= 2)
        {
          const double anA = OutwardNormalAt(TopoDS::Face(aFs.First()), aMid)
                               .Angle(OutwardNormalAt(TopoDS::Face(aFs.Last()), aMid));
          aMove = aR * (1. / std::cos(0.5 * anA) - 1.);
        }
      }
      BRepExtrema_DistShapeShape aDist(BRepBuilderAPI_MakeVertex(aMid).Vertex(), aFaces);
      if (!aDist.IsDone() || aDist.Value() < std::max(0.5 * std::min(aMove, aR), 1.e-6))
      {
        return false;
      }
    }
  }
  return true;
}

// Whether a face of <theResult> made by the computation -- not on a surface
// of <theInput> -- and not elementary folds over itself within <theReach>
// of one of <theVertices>, or anywhere when there are none: a line of its
// intersection with itself, as BOP's check computes it, longer than 1e-4
// and on two sheets -- its middle at two places in (u, v). A set back
// corner's patch across the 22 deg apex of #273's pocket, from the top's
// edge down to the apex at the wall's line, was valid and folded (#273 E23
// and E24 at 0.8, lines 0.3 long; on Fillet001, whose loose input leaves
// the patch at 0.038, 0.13 long). The line BOP finds on #962's three-sided
// plate at its pole is a point, of no length, at one place in (u, v): no
// fold.
bool HasFoldNear(const TopoDS_Shape&                                    theResult,
                 const TopoDS_Shape&                                    theInput,
                 const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& theVertices,
                 const double                                           theReach)
{
  NCollection_Map<occ::handle<Standard_Transient>> aSurfIn;
  for (TopExp_Explorer anExp(theInput, TopAbs_FACE); anExp.More(); anExp.Next())
  {
    TopLoc_Location aL;
    aSurfIn.Add(BRep_Tool::Surface(TopoDS::Face(anExp.Current()), aL));
  }
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aSeen;
  for (TopExp_Explorer anExp(theResult, TopAbs_FACE); anExp.More(); anExp.Next())
  {
    const TopoDS_Face& aF = TopoDS::Face(anExp.Current());
    TopLoc_Location    aL;
    if (!aSeen.Add(aF.Located(TopLoc_Location())) || aSurfIn.Contains(BRep_Tool::Surface(aF, aL)))
    {
      continue;
    }
    BRepAdaptor_Surface       aS(aF, false);
    const GeomAbs_SurfaceType aType = aS.GetType();
    if (aType == GeomAbs_Plane || aType == GeomAbs_Cylinder || aType == GeomAbs_Cone
        || aType == GeomAbs_Sphere
        || (aType == GeomAbs_Torus
            && aS.Torus().MajorRadius() > aS.Torus().MinorRadius() + Precision::Confusion()))
    {
      continue;
    }
    Bnd_Box aBox;
    BRepBndLib::Add(aF, aBox);
    aBox.Enlarge(theReach);
    bool isNear = theVertices.IsEmpty();
    for (NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>::Iterator itV(theVertices);
         itV.More() && !isNear;
         itV.Next())
    {
      isNear = !aBox.IsOut(BRep_Tool::Pnt(TopoDS::Vertex(itV.Key())));
    }
    if (!isNear)
    {
      continue;
    }
    IntTools_FaceFace anFF;
    anFF.Perform(aF, aF, false);
    if (!anFF.IsDone())
    {
      continue;
    }
    for (int i = 1; i <= anFF.Lines().Length(); i++)
    {
      const IntTools_Curve&            aLine = anFF.Lines()(i);
      const occ::handle<Geom_Curve>&   aC    = aLine.Curve();
      const occ::handle<Geom2d_Curve>& aC1   = aLine.FirstCurve2d();
      const occ::handle<Geom2d_Curve>& aC2   = aLine.SecondCurve2d();
      if (aC.IsNull() || aC1.IsNull() || aC2.IsNull())
      {
        continue;
      }
      const double a = aC->FirstParameter(), b = aC->LastParameter(), m = 0.5 * (a + b);
      const gp_Pnt aP0 = aC->Value(a), aPm = aC->Value(m), aP1 = aC->Value(b);
      if (aP0.Distance(aPm) + aPm.Distance(aP1) > 1.e-4
          && aC1->Value(m).Distance(aC2->Value(m)) > 100. * Precision::PConfusion())
      {
        return true;
      }
    }
  }
  return false;
}

// The largest radius at the end of <theSp> numbered <theIE> in the spine.
double LargestRadiusAt(const occ::handle<ChFiDS_FilSpine>& theSp, const int theIE)
{
  return theSp->IsConstant(theIE) ? theSp->Radius(theIE) : theSp->MaxRadFromSeqAndLaws();
}

// Whether a fillet stripe of <theStripes> ends at one of <theVertices>, as
// the setback fallback sets ends back: a chamfer's does not.
bool HasFilletEndAt(const NCollection_List<occ::handle<ChFiDS_Stripe>>&          theStripes,
                    const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& theVertices)
{
  for (NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator itS(theStripes); itS.More();
       itS.Next())
  {
    const occ::handle<ChFiDS_FilSpine> aSp = occ::down_cast<ChFiDS_FilSpine>(itS.Value()->Spine());
    if (aSp.IsNull())
    {
      continue;
    }
    for (int k = 0; k < 2; k++)
    {
      const bool isFirst = k == 0;
      if (theVertices.Contains(isFirst ? aSp->FirstVertex() : aSp->LastVertex())
          && aSp->Status(isFirst) != ChFiDS_FreeBoundary && (isFirst || !aSp->IsPeriodic()))
      {
        return true;
      }
    }
  }
  return false;
}
} // namespace

//=================================================================================================

bool ChFi3d_Builder::HasSetbackAt(const int Index) const
{
  const TopoDS_Vertex& aV = myVDataMap.FindKey(Index);
  for (NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator it(myVDataMap(Index)); it.More();
       it.Next())
  {
    const occ::handle<ChFiDS_FilSpine> aSp = occ::down_cast<ChFiDS_FilSpine>(it.Value()->Spine());
    if (aSp.IsNull())
    {
      continue;
    }
    int aSens = 0;
    ChFi3d_IndexOfSurfData(aV, it.Value(), aSens);
    const bool isFirst = aSens == 1;
    if (aSp->Setback(isFirst) >= 0. && aSp->Status(isFirst) != ChFiDS_FreeBoundary)
    {
      return true;
    }
  }
  return false;
}

//=================================================================================================

void ChFi3d_Builder::ComputeLongExtension()
{
  // a fillet stripe on edges short for its radius at an end it is extended
  // at (ExtentOneCorner): half its length there, under 1.5 radius
  bool   isShort = false;
  double aRadius = Precision::Infinite();
  for (NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator itS(myListStripe); itS.More();
       itS.Next())
  {
    const occ::handle<ChFiDS_FilSpine> aSp = occ::down_cast<ChFiDS_FilSpine>(itS.Value()->Spine());
    if (aSp.IsNull() || aSp->IsPeriodic())
    {
      continue;
    }
    const double aLength = aSp->LastParameter(aSp->NbEdges());
    for (int k = 0; k < 2; k++)
    {
      const double aRad  = LargestRadiusAt(aSp, k == 0 ? 1 : aSp->NbEdges());
      const bool   isEnd = !aSp->IsTangencyExtremity(k == 0);
      aRadius            = std::min(aRadius, aRad);
      isShort            = isShort || (isEnd && 0.5 * aLength < 1.5 * aRad);
    }
  }
  if (!isShort)
  {
    return;
  }
  // as the setback fallback takes a result: valid, every fillet made, no
  // edge looser than the input's or a twentieth of the radius, and no face
  // it made folded over itself
  const double aMaxTol = std::max(MaxEdgeTolerance(myShape), 0.05 * aRadius);
  auto         FreshComputation = [&]() {
    Reset();
    myCoup = new TopOpeBRepBuild_HBuilder(myCoup->BuildTool());
  };
  // on for the thread while this runs, whatever way it leaves
  struct LongExtensionGuard
  {
    LongExtensionGuard() { ChFi3d_LongSpineExtension() = true; }

    ~LongExtensionGuard() { ChFi3d_LongSpineExtension() = false; }
  } aLongExtension;
  try
  {
    OCC_CATCH_SIGNALS
    FreshComputation();
    Compute();
    if (done && IsValidResult(myShapeResult, aMaxTol)
        && AreFilletsMade(myListStripe, myEFMap, myShapeResult)
        && !HasFoldNear(myShapeResult, myShape, {}, 0.))
    {
      return;
    }
  }
  catch (Standard_Failure const&)
  {
    done = false;
  }
  // and failing at a vertex, set back there as well
  if (!done)
  {
    try
    {
      OCC_CATCH_SIGNALS
      ComputeSetbackFallback();
    }
    catch (Standard_Failure const&)
    {
      done = false;
    }
    if (done)
    {
      return;
    }
  }
  // nothing better: the computation as it was, its failure with it
  done = false;
  try
  {
    OCC_CATCH_SIGNALS
    FreshComputation();
    Compute();
  }
  catch (Standard_Failure const&)
  {
    done = false;
  }
}

//=================================================================================================

void ChFi3d_Builder::ComputeSetbackFallback()
{
  const double aMultiple = CornerSetbackFallback();
  if (done || aMultiple <= 0. || badvertices.IsEmpty())
  {
    return;
  }
  // the fillet stripes' ends at the vertices the computation failed at, the
  // setbacks they had, and the largest radius at each of those vertices
  struct SetbackEnd
  {
    occ::handle<ChFiDS_FilSpine> Spine;
    bool                         IsFirst;
    double                       Before;
    double                       Radius;
  };
  NCollection_Sequence<SetbackEnd>                anEnds;
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aVertices;
  // the vertices failed at that are not set back yet; true where one is added
  auto addFailed = [&]() {
    bool isAdded = false;
    for (NCollection_List<TopoDS_Shape>::Iterator itV(badvertices); itV.More(); itV.Next())
    {
      const TopoDS_Vertex& aV = TopoDS::Vertex(itV.Value());
      if (!aVertices.Add(aV))
      {
        continue;
      }
      const int aFirstE = anEnds.Length() + 1;
      double    aRadius = 0.;
      for (NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator itS(myListStripe); itS.More();
           itS.Next())
      {
        const occ::handle<ChFiDS_FilSpine> aSp =
          occ::down_cast<ChFiDS_FilSpine>(itS.Value()->Spine());
        if (aSp.IsNull())
        {
          continue;
        }
        for (int k = 0; k < 2; k++)
        {
          const bool isFirst = k == 0;
          if (!(isFirst ? aSp->FirstVertex() : aSp->LastVertex()).IsSame(aV)
              || aSp->Status(isFirst) == ChFiDS_FreeBoundary || (!isFirst && aSp->IsPeriodic()))
          {
            continue;
          }
          aRadius = std::max(aRadius, LargestRadiusAt(aSp, isFirst ? 1 : aSp->NbEdges()));
          anEnds.Append({aSp, isFirst, aSp->Setback(isFirst), 0.});
        }
      }
      for (int i = aFirstE; i <= anEnds.Length(); i++)
      {
        anEnds.ChangeValue(i).Radius = aRadius;
      }
      isAdded = isAdded || aFirstE <= anEnds.Length();
    }
    return isAdded;
  };
  if (!addFailed())
  {
    return;
  }
  // how loose an edge of a result may be: as the input's loosest, or a
  // twentieth of the smallest radius set back -- today's fillets keep edges
  // of up to about a thirtieth of theirs (the plate fallback's measure)
  const double aMaxTolIn = MaxEdgeTolerance(myShape);
  // what a computation that failed leaves behind, cleared for the next: the
  // stripes it left without a spine, which Compute reads before its own
  // Reset, and the topological builder, which a Perform broken off by an
  // exception left half cleared (a crash in the next Perform; the same
  // builder computed twice crashes upstream as well)
  auto FreshComputation = [&]() {
    Reset();
    myCoup = new TopOpeBRepBuild_HBuilder(myCoup->BuildTool());
  };
  auto aMaxTol = [&]() {
    double aRadius = Precision::Infinite();
    for (int i = 1; i <= anEnds.Length(); i++)
    {
      if (anEnds.Value(i).Radius > 0.)
      {
        aRadius = std::min(aRadius, anEnds.Value(i).Radius);
      }
    }
    return std::max(aMaxTolIn, aRadius < Precision::Infinite() ? 0.05 * aRadius : 0.);
  };
  auto aRadiusMax = [&]() {
    double aRadius = 0.;
    for (int i = 1; i <= anEnds.Length(); i++)
    {
      aRadius = std::max(aRadius, anEnds.Value(i).Radius);
    }
    return aRadius;
  };

  // The setbacks the fallback tries are consumed by Compute(); the spines keep
  // the caller's own afterwards, kept or not, as a Build() after this one
  // starts from what the caller set.
  auto aRestore = [&]() {
    for (int i = 1; i <= anEnds.Length(); i++)
    {
      const SetbackEnd& anEnd = anEnds.Value(i);
      anEnd.Spine->SetSetback(anEnd.IsFirst, anEnd.Before);
    }
  };

  // where the stripes meet, then 1, 1.5, 2... times the radius; a step that
  // fails at a vertex not set back yet is taken again with that one set back
  // too, as a corner set back can move the failure to the next
  for (double aStep = 0.; aStep <= aMultiple + Precision::Confusion();)
  {
    for (int i = 1; i <= anEnds.Length(); i++)
    {
      const SetbackEnd& anEnd = anEnds.Value(i);
      anEnd.Spine->SetSetback(anEnd.IsFirst, std::max(anEnd.Before, aStep * anEnd.Radius));
    }
    try
    {
      OCC_CATCH_SIGNALS
      FreshComputation();
      Compute();
      // and no patch of the corners set back folded over itself, within the
      // setback and twice the radius of their vertices
      if (done && IsValidResult(myShapeResult, aMaxTol())
          && AreFilletsMade(myListStripe, myEFMap, myShapeResult)
          && !HasFoldNear(myShapeResult,
                          myShape,
                          aVertices,
                          (std::max(aStep, 1.) + 2.) * aRadiusMax()))
      {
        aRestore();
        return;
      }
    }
    catch (Standard_Failure const&)
    {
    }
    if (!addFailed())
    {
      aStep = aStep < 1. ? 1. : aStep + 0.5;
    }
  }

  // nothing valid: the setbacks as they were, and the failure as it was
  aRestore();
  FreshComputation();
  Compute();
}

//=======================================================================
// function : PerformOppositeCorner
// purpose  : Two fillets ending at a vertex of four edges on no common face
//           (OppositeEnds). Each is ended as it would be alone
//           (PerformIntersectionAtEnd): closed by the faces beyond the
//           vertex carried on, so that its cap lies on the other fillet's
//           two faces, and the other's on its own. The DS is then told how
//           the two ends share the vertex, for the result the union of the
//           two made alone would give. Of the four edges, each spine's keeps
//           the stub from the vertex to the other fillet's cap (its own side
//           goes under its fillet); each sharp edge keeps the piece from the
//           vertex to its nearer contact -- concave between the two caps --
//           loses the piece between its two contacts, where the caps meet on
//           one plane, and keeps what lies past the farther. The piece to
//           the nearer contact bounds its faces from the other side, so it
//           is a curve of its own in the DS, reversed on both. The ends must
//           come out as that pattern, one new point on each spine edge beyond
//           the vertex and two on each sharp edge, or the corner is refused.
//=======================================================================

void ChFi3d_Builder::PerformOppositeCorner(const int Index)
{
  const TopoDS_Vertex&        Vtx  = myVDataMap.FindKey(Index);
  TopOpeBRepDS_DataStructure& DStr = myDS->ChangeDS();
  NCollection_List<occ::handle<ChFiDS_Stripe>>& aStripes =
    const_cast<NCollection_List<occ::handle<ChFiDS_Stripe>>&>(myVDataMap(Index));

  // the four edges at the vertex and their interferences so far
  NCollection_Sequence<TopoDS_Edge> anEdges;
  for (NCollection_List<TopoDS_Shape>::Iterator itE(myVEMap(Vtx)); itE.More(); itE.Next())
  {
    const TopoDS_Edge& anE  = TopoDS::Edge(itE.Value());
    bool               isIn = BRep_Tool::Degenerated(anE);
    for (int j = 1; j <= anEdges.Length() && !isIn; j++)
    {
      isIn = anEdges(j).IsSame(anE);
    }
    if (!isIn)
    {
      anEdges.Append(anE);
    }
  }
  if (anEdges.Length() != 4)
  {
    throw Standard_Failure("OppositeCorner : not four edges at the vertex");
  }
  NCollection_Sequence<int> aNbBefore;
  for (int j = 1; j <= anEdges.Length(); j++)
  {
    aNbBefore.Append(DStr.HasShape(anEdges(j)) ? DStr.ShapeInterferences(anEdges(j)).Extent() : 0);
  }
  // the spines' edges at the vertex
  NCollection_Sequence<TopoDS_Edge> aSpineEdges;
  for (NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator It(aStripes); It.More(); It.Next())
  {
    const occ::handle<ChFiDS_Spine>& aS = It.Value()->Spine();
    aSpineEdges.Append(aS->FirstVertex().IsSame(Vtx) ? aS->Edges(1) : aS->Edges(aS->NbEdges()));
  }

  // each end as alone: PerformIntersectionAtEnd takes the list's first
  for (int k = 0; k < 2; k++)
  {
    PerformIntersectionAtEnd(Index);
    aStripes.Append(aStripes.First());
    aStripes.RemoveFirst();
  }

  const int iV = DStr.AddShape(Vtx);
  for (int j = 1; j <= anEdges.Length(); j++)
  {
    const TopoDS_Edge& anE     = anEdges(j);
    const bool         isSpine = anE.IsSame(aSpineEdges(1)) || anE.IsSame(aSpineEdges(2));
    // a spine's edge with no stub: gone from the vertex on, told so, or the
    // state of the edge would be taken from the vertex, where the curves
    // bounding a cap are kept
    auto NoStub = [&]() {
      const double pV0 = BRep_Tool::Parameter(Vtx, anE);
      double       f0, l0;
      BRep_Tool::Range(anE, f0, l0);
      // in the DS first: the interferences of an edge not in it go nowhere
      const int iE = DStr.AddShape(anE);
      DStr.ChangeShapeInterferences(iE).Append(
        ChFi3d_FilPointInDS(std::abs(pV0 - f0) <= std::abs(pV0 - l0) ? TopAbs_REVERSED
                                                                      : TopAbs_FORWARD,
                            iE,
                            iV,
                            pV0,
                            true));
    };
    if (!DStr.HasShape(anE))
    {
      if (isSpine)
      {
        // the other's cap passes the spine's edge by, in the plane of the
        // cap: no stub
        NoStub();
        continue;
      }
      throw Standard_Failure("OppositeCorner : an edge at the vertex got no end");
    }
    NCollection_List<occ::handle<TopOpeBRepDS_Interference>>& aLI =
      DStr.ChangeShapeInterferences(anE);
    const double pV    = BRep_Tool::Parameter(Vtx, anE);
    const double aResV = BRepAdaptor_Curve(anE).Resolution(BRep_Tool::Tolerance(Vtx));
    NCollection_Sequence<occ::handle<TopOpeBRepDS_CurvePointInterference>> aNew;
    int                                                                     n = 0;
    for (NCollection_List<occ::handle<TopOpeBRepDS_Interference>>::Iterator itI(aLI); itI.More();
         itI.Next(), n++)
    {
      occ::handle<TopOpeBRepDS_CurvePointInterference> CP =
        occ::down_cast<TopOpeBRepDS_CurvePointInterference>(itI.Value());
      if (n >= aNbBefore(j) && !CP.IsNull())
      {
        aNew.Append(CP);
      }
    }
    auto RemoveNew = [&](const occ::handle<TopOpeBRepDS_Interference>& theI) {
      for (NCollection_List<occ::handle<TopOpeBRepDS_Interference>>::Iterator itI(aLI); itI.More();)
      {
        if (itI.Value() == theI)
        {
          aLI.Remove(itI);
        }
        else
        {
          itI.Next();
        }
      }
    };
    if (isSpine)
    {
      if (aNew.IsEmpty())
      {
        NoStub();
        continue;
      }
      if (aNew.Length() != 1)
      {
        throw Standard_Failure("OppositeCorner : the ends are not two ends alone");
      }
      // the stub to the other's cap kept, the spine's side gone
      const double pA = aNew(1)->Parameter();
      if (std::abs(pA - pV) <= aResV)
      {
        // the other's cap touches the spine's edge at the vertex only: no
        // stub, the edge goes under its fillet whole
        RemoveNew(aNew(1));
        NoStub();
        continue;
      }
      aLI.Append(ChFi3d_FilPointInDS(pA < pV ? TopAbs_REVERSED : TopAbs_FORWARD,
                                     DStr.AddShape(anE),
                                     iV,
                                     pV,
                                     true));
      continue;
    }
    // A sharp edge's contacts: where the two ends cross it. One of them can
    // be a fillet's own line ending on the edge, which FilDS puts on it
    // after the corners (its point already in the DS): a line that ends
    // beside the vertex, its cap on the other's face past the edge.
    struct Contact
    {
      double                                           Par = 0.;
      int                                              Point = 0;
      occ::handle<TopOpeBRepDS_CurvePointInterference> CP;      // on the edge now
      ChFiDS_CommonPoint*                              Pending = nullptr; // FilDS's to put
    };
    NCollection_Sequence<Contact> aCt;
    for (int k = 1; k <= aNew.Length(); k++)
    {
      Contact c;
      c.Par   = aNew(k)->Parameter();
      c.Point = aNew(k)->Geometry();
      c.CP    = aNew(k);
      aCt.Append(c);
    }
    for (NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator It(aStripes);
         It.More() && aCt.Length() < 2;
         It.Next())
    {
      int       aSens = 0;
      const int aNum  = ChFi3d_IndexOfSurfData(Vtx, It.Value(), aSens);
      const bool isf  = aSens == 1;
      const occ::handle<ChFiDS_SurfData>& aSD = It.Value()->SetOfSurfData()->Value(aNum);
      for (int k = 1; k <= 2; k++)
      {
        ChFiDS_CommonPoint& aCP = aSD->ChangeVertex(isf, k);
        const int           iP  = It.Value()->IndexPoint(isf, k);
        bool                isIn = !aCP.IsOnArc() || !aCP.Arc().IsSame(anE) || iP <= 0
                    || std::abs(aCP.ParameterOnArc() - pV) <= aResV;
        for (int m = 1; m <= aCt.Length() && !isIn; m++)
        {
          isIn = aCt(m).Point == iP;
        }
        if (!isIn)
        {
          Contact c;
          c.Par     = aCP.ParameterOnArc();
          c.Point   = iP;
          c.Pending = &aCP;
          aCt.Append(c);
        }
      }
    }
    if (aCt.Length() == 1)
    {
      // the nearer contact the vertex itself: gone to the one contact, kept
      // past it
      if (aCt(1).CP.IsNull() || std::abs(aCt(1).Par - pV) <= aResV)
      {
        throw Standard_Failure("OppositeCorner : a sharp edge's one contact");
      }
      aCt(1).CP->Transition(
        TopOpeBRepDS_Transition(aCt(1).Par > pV ? TopAbs_FORWARD : TopAbs_REVERSED));
      continue;
    }
    if (aCt.Length() != 2)
    {
      throw Standard_Failure("OppositeCorner : the ends are not two ends alone");
    }
    // the vertex to the nearer contact kept, between the two gone, past the
    // farther kept
    const int    iN = std::abs(aCt(1).Par - pV) <= std::abs(aCt(2).Par - pV) ? 1 : 2;
    const double pN = aCt(iN).Par, pF = aCt(3 - iN).Par;
    const TopAbs_Orientation aFarOr = pF > pV ? TopAbs_FORWARD : TopAbs_REVERSED;
    if (!aCt(3 - iN).CP.IsNull())
    {
      aCt(3 - iN).CP->Transition(TopOpeBRepDS_Transition(aFarOr));
    }
    else
    {
      ChFiDS_CommonPoint& aCP = *aCt(3 - iN).Pending;
      aCP.SetArc(aCP.Tolerance(), anE, aCP.ParameterOnArc(), aFarOr);
    }
    // the piece to the nearer contact: a curve of its own, the faces on its
    // other side; the contact not on the edge
    const int iPN = aCt(iN).Point;
    if (!aCt(iN).CP.IsNull())
    {
      RemoveNew(aCt(iN).CP);
    }
    else
    {
      ChFiDS_CommonPoint& aCP  = *aCt(iN).Pending;
      const double        aTol = aCP.Tolerance();
      const gp_Pnt        aP   = aCP.Point();
      const bool          hasV = aCP.HasVector();
      const gp_Vec        aVec = hasV ? aCP.Vector() : gp_Vec();
      aCP.Reset();
      aCP.SetPoint(aP);
      aCP.SetTolerance(aTol);
      if (hasV)
      {
        aCP.SetVector(aVec);
      }
    }
    double                         f, l;
    occ::handle<Geom_Curve>        aC = BRep_Tool::Curve(anE, f, l);
    const double                   p1 = std::min(pV, pN), p2 = std::max(pV, pN);
    occ::handle<Geom_TrimmedCurve> aTC = new Geom_TrimmedCurve(aC, p1, p2);
    const int iC = DStr.AddCurve(TopOpeBRepDS_Curve(aTC, BRep_Tool::Tolerance(anE)));
    DStr.ChangeCurveInterferences(iC).Append(
      ChFi3d_FilPointInDS(TopAbs_FORWARD, iC, pV < pN ? iV : iPN, p1, pV < pN));
    DStr.ChangeCurveInterferences(iC).Append(
      ChFi3d_FilPointInDS(TopAbs_REVERSED, iC, pV < pN ? iPN : iV, p2, pV > pN));
    for (NCollection_List<TopoDS_Shape>::Iterator itF(myEFMap(anE)); itF.More(); itF.Next())
    {
      const TopoDS_Face& aF  = TopoDS::Face(itF.Value());
      TopAbs_Orientation anO = TopAbs_FORWARD;
      for (TopExp_Explorer exE(aF.Oriented(TopAbs_FORWARD), TopAbs_EDGE); exE.More(); exE.Next())
      {
        if (exE.Current().IsSame(anE))
        {
          anO = exE.Current().Orientation();
        }
      }
      double                    f2, l2;
      occ::handle<Geom2d_Curve> aPC = BRep_Tool::CurveOnSurface(anE, aF, f2, l2);
      const int                 iF  = DStr.AddShape(aF);
      DStr.ChangeShapeInterferences(iF).Append(
        ChFi3d_FilCurveInDS(iC, iF, aPC, TopAbs::Reverse(anO)));
    }
  }
}

//=================================================================================================

void ChFi3d_Builder::ComputeOppositeCorners()
{
  if (THE_OPPOSITE_SEEN.IsEmpty())
  {
    return;
  }
  // The plate's result is kept where it is valid: of the vertex sweep's 392
  // opposite pairs the plate had right or loose, the corner below broke 61
  // and mended 7 of the 21 it had wrong (allopp, 2026-10-21), so it only
  // stands in for a plate that failed. Valid on the setback fallback's
  // terms: no edge looser than the input's or a twentieth of the smallest
  // radius at those vertices.
  double aRadius = Precision::Infinite();
  for (NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator itS(myListStripe); itS.More();
       itS.Next())
  {
    const occ::handle<ChFiDS_FilSpine> aSp = occ::down_cast<ChFiDS_FilSpine>(itS.Value()->Spine());
    if (aSp.IsNull())
    {
      continue;
    }
    for (int k = 0; k < 2; k++)
    {
      const bool isFirst = k == 0;
      if (THE_OPPOSITE_SEEN.Contains(isFirst ? aSp->FirstVertex() : aSp->LastVertex()))
      {
        aRadius = std::min(aRadius, LargestRadiusAt(aSp, isFirst ? 1 : aSp->NbEdges()));
      }
    }
  }
  if (Precision::IsInfinite(aRadius))
  {
    aRadius = 0.;
  }
  const double aMaxTol = std::max(MaxEdgeTolerance(myShape), 0.05 * aRadius);
  if (done && IsValidResult(myShapeResult, aMaxTol))
  {
    return;
  }
  auto FreshComputation = [&]() {
    Reset();
    myCoup = new TopOpeBRepBuild_HBuilder(myCoup->BuildTool());
  };
  // the corners at every such vertex; the pieces of the edges carried past
  // the vertex ask TopOpeBRepBuild to keep the vertex's own pave beside the
  // one at their end (TopOpeBRepBuild_PaveSet_PastEndPaves), for this
  // computation only
  THE_OPPOSITE_CORNERS.Assign(THE_OPPOSITE_SEEN);
  try
  {
    OCC_CATCH_SIGNALS
    struct PastEndGuard
    {
      PastEndGuard() { TopOpeBRepBuild_PaveSet_PastEndPaves() = true; }

      ~PastEndGuard() { TopOpeBRepBuild_PaveSet_PastEndPaves() = false; }
    } aPastEnd;
    FreshComputation();
    Compute();
    if (done && IsValidResult(myShapeResult, aMaxTol)
        && AreFilletsMade(myListStripe, myEFMap, myShapeResult)
        && !HasFoldNear(myShapeResult, myShape, {}, 0.))
    {
      THE_OPPOSITE_CORNERS.Clear();
      return;
    }
  }
  catch (Standard_Failure const&)
  {
  }
  // nothing better: the computation as it was, its failure with it
  THE_OPPOSITE_CORNERS.Clear();
  done = false;
  try
  {
    OCC_CATCH_SIGNALS
    FreshComputation();
    Compute();
  }
  catch (Standard_Failure const&)
  {
    done = false;
  }
}

//=================================================================================================

void ChFi3d_Builder::Compute()
{
  if (!THE_IN_SETBACK_COMPUTE)
  {
    // the computation; where it fails at a vertex, the setback fallback;
    // where the result is still invalid or missing and it met two fillets on
    // no common face, the opposite corner there (ComputeOppositeCorners) --
    // after the setback fallback, which makes such corners where the plate
    // failed (the floor and the inline edge of a slot by a slope), and would
    // start from the opposite corner's failure instead of the plate's; where
    // the result is still invalid, the corners the plate made set back, and
    // where it is invalid after that, the plates made again stepped back on
    // their sharp edges only to points on them; then, on edges short for their radius, the spines extended further
    // (ComputeLongExtension). A failure the fallbacks do not mend is the
    // computation's own
    SetbackComputeGuard aGuard;
    struct InputTolerancesGuard
    {
      InputTolerancesGuard(const TopoDS_Shape& theShape)
      {
        THE_INPUT_TOLERANCES.Clear();
        for (TopExp_Explorer exV(theShape, TopAbs_VERTEX); exV.More(); exV.Next())
        {
          THE_INPUT_TOLERANCES.Bind(exV.Current(),
                                    BRep_Tool::Tolerance(TopoDS::Vertex(exV.Current())));
        }
        for (TopExp_Explorer exE(theShape, TopAbs_EDGE); exE.More(); exE.Next())
        {
          THE_INPUT_TOLERANCES.Bind(exE.Current(),
                                    BRep_Tool::Tolerance(TopoDS::Edge(exE.Current())));
        }
      }

      ~InputTolerancesGuard() { THE_INPUT_TOLERANCES.Clear(); }
    } anInputTolerances(myShape);
    struct OppositeGuard
    {
      OppositeGuard() { Clear(); }

      ~OppositeGuard() { Clear(); }

      static void Clear()
      {
        THE_OPPOSITE_CORNERS.Clear();
        THE_OPPOSITE_SEEN.Clear();
        ChFi3d_PlateCorners().Clear();
        ChFi3d_RunOnSeams().Clear();
        ChFi3d_SideSeams().Clear();
      }
    } anOpposite;
    // A fresh builder of the result, as each fallback takes one: Reset()
    // leaves myCoup holding the pieces of the computation before, and a
    // Build again after a result was refused -- BRepBuilderAPI_MakeShape::
    // Shape() runs one -- made the input from them, unchanged, and done
    myCoup = new TopOpeBRepBuild_HBuilder(myCoup->BuildTool());
    std::exception_ptr aFailure;
    try
    {
      Compute();
    }
    catch (Standard_Failure const&)
    {
      aFailure = std::current_exception();
      done     = false;
    }
    // the corners the first computation built with the plate
    const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aPlates = ChFi3d_PlateCorners();
    ComputeSetbackFallback();
    ComputeOppositeCorners();
    // A result built, but invalid, where a corner was the plate: those
    // corners set back, as where the computation fails at them (#962's
    // groove, its floor's edge at the mouth and its wall's foot together:
    // a corner of opposite concavity at four edges, the plate over a fillet
    // whose line on the outer face ran on under the tooth). Nothing valid,
    // the fallback computes it as it was.
    if (done && !aPlates.IsEmpty() && CornerSetbackFallback() > 0.
        && HasFilletEndAt(myListStripe, aPlates) && !BRepCheck_Analyzer(myShapeResult).IsValid())
    {
      badvertices.Clear();
      for (NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>::Iterator itV(aPlates);
           itV.More();
           itV.Next())
      {
        badvertices.Append(itV.Key());
      }
      done = false;
      ComputeSetbackFallback();
    }
    // Invalid still, where a corner was the plate: the plates made again,
    // their sharp edges stepped back only to the stripes' common points on
    // them (ChFi3d_PlateStepBackOnEdge; #474's coil leaving its block, at
    // 1). Kept where valid, as the setback fallback takes a result; else
    // the result as it was.
    if (done && !aPlates.IsEmpty() && !BRepCheck_Analyzer(myShapeResult).IsValid())
    {
      // the smallest and largest radius of the fillets ending at those
      // corners, for the bounds the setback fallback takes a result by
      double aRadMin = Precision::Infinite(), aRadMax = 0.;
      for (NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator itS(myListStripe); itS.More();
           itS.Next())
      {
        const occ::handle<ChFiDS_FilSpine> aSp =
          occ::down_cast<ChFiDS_FilSpine>(itS.Value()->Spine());
        for (int k = 0; k < 2 && !aSp.IsNull(); k++)
        {
          if (aPlates.Contains(k == 0 ? aSp->FirstVertex() : aSp->LastVertex()))
          {
            const double aRad = LargestRadiusAt(aSp, k == 0 ? 1 : aSp->NbEdges());
            aRadMin           = std::min(aRadMin, aRad);
            aRadMax           = std::max(aRadMax, aRad);
          }
        }
      }
      const double aMaxTol =
        std::max(MaxEdgeTolerance(myShape), aRadMin < Precision::Infinite() ? 0.05 * aRadMin : 0.);
      auto FreshComputation = [&]() {
        Reset();
        myCoup = new TopOpeBRepBuild_HBuilder(myCoup->BuildTool());
      };
      // on for the thread while this runs, whatever way it leaves
      struct StepBackGuard
      {
        StepBackGuard() { ChFi3d_PlateStepBackOnEdge() = true; }

        ~StepBackGuard() { ChFi3d_PlateStepBackOnEdge() = false; }
      };
      bool isKept = false;
      try
      {
        OCC_CATCH_SIGNALS
        StepBackGuard aStepBack;
        FreshComputation();
        Compute();
        isKept = done && IsValidResult(myShapeResult, aMaxTol)
                 && AreFilletsMade(myListStripe, myEFMap, myShapeResult)
                 && !HasFoldNear(myShapeResult, myShape, aPlates, 3. * aRadMax);
      }
      catch (Standard_Failure const&)
      {
        isKept = false;
      }
      if (!isKept)
      {
        done = false;
        try
        {
          OCC_CATCH_SIGNALS
          FreshComputation();
          Compute();
        }
        catch (Standard_Failure const&)
        {
          done = false;
        }
      }
    }
    // A result BRepCheck finds invalid is no result: refused, the computation
    // not done, the shape kept for BadShape(). Before the long extension,
    // which may make it valid, and after it, which may give it back as it was.
    // So is one with the radii a hair smaller that kept a sliver
    const auto RefuseInvalid = [this]() {
      if (done
          && (!BRepCheck_Analyzer(myShapeResult).IsValid()
              || (ChFi3d_SnapSlivers() && ChFi3d_HasSliverLeft(myShapeResult, myShape, tolapp3d))))
      {
        badShape  = myShapeResult;
        hasresult = true;
        done      = false;
      }
    };
    RefuseInvalid();
    if (!done)
    {
      ComputeLongExtension();
    }
    RefuseInvalid();
    if (aFailure && !done)
    {
      std::rethrow_exception(aFailure);
    }
    return;
  }

#ifdef OCCT_DEBUG // perf
  t_total              = 0;
  t_extent             = 0;
  t_perfsetofsurf      = 0;
  t_perffilletonvertex = 0;
  t_filds              = 0;
  t_reconstruction     = 0;
  t_setregul           = 0;
  t_perfsetofkpart     = 0;
  t_perfsetofkgen      = 0;
  t_makextremities     = 0;
  t_performsurf        = 0;
  t_startsol           = 0;
  t_perform1corner     = 0;
  t_perform2corner     = 0;
  t_perform3corner     = 0;
  t_performmore3corner = 0;
  t_inter              = 0;
  t_same               = 0;
  t_sameinter          = 0;
  t_plate              = 0;
  t_approxplate        = 0;
  t_batten             = 0;
  t_remplissage        = 0;
  t_t3cornerinit       = 0;
  t_spherique          = 0;
  t_torique            = 0;
  t_notfilling         = 0;
  t_filling            = 0;
  t_performatend       = 0;
  t_t2cornerinit       = 0;
  t_perf2cornerbyinter = 0;
  t_chfikpartcompdata  = 0;
  t_cheminement        = 0;
  t_sameparam          = 0;
  t_computedata        = 0;
  t_completedata       = 0;
  t_t2cornerDS         = 0;
  t_t3cornerDS         = 0;
  ChFi3d_InitChron(cl_total);
  ChFi3d_InitChron(cl_extent);
#endif
  // the changes the computation makes to the input's frozen vertices, kept for
  // copies of them in the result (ChFi3d_UpdateVertex)
  ChFi3d_ThawScope aThaw;
  UpdateTolesp();

  if (myListStripe.IsEmpty())
  {
    throw Standard_Failure("There are no suitable edges for chamfer or fillet");
  }

  Reset();
  myDS                             = new TopOpeBRepDS_HDataStructure();
  TopOpeBRepDS_DataStructure& DStr = myDS->ChangeDS();
  done                             = true;
  hasresult                        = false;

  // filling of myVDatatMap
  NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator itel;

  for (itel.Initialize(myListStripe); itel.More(); itel.Next())
  {
    if ((itel.Value()->Spine()->FirstStatus() <= ChFiDS_BreakPoint))
    {
      myVDataMap.Add(itel.Value()->Spine()->FirstVertex(), itel.Value());
    }
    else if (itel.Value()->Spine()->FirstStatus() == ChFiDS_FreeBoundary)
    {
      ExtentOneCorner(itel.Value()->Spine()->FirstVertex(), itel.Value());
    }
    if ((itel.Value()->Spine()->LastStatus() <= ChFiDS_BreakPoint))
    {
      myVDataMap.Add(itel.Value()->Spine()->LastVertex(), itel.Value());
    }
    else if (itel.Value()->Spine()->LastStatus() == ChFiDS_FreeBoundary)
    {
      ExtentOneCorner(itel.Value()->Spine()->LastVertex(), itel.Value());
    }
  }
  // preanalysis to evaluate the extensions.
  ExtentAnalyse();

#ifdef OCCT_DEBUG // perf
  ChFi3d_ResultChron(cl_extent, t_extent);
  ChFi3d_InitChron(cl_perfsetofsurf);
#endif

  // Construction of the stripe of fillet on each stripe.
  for (itel.Initialize(myListStripe); itel.More(); itel.Next())
  {
    itel.Value()->Spine()->SetErrorStatus(ChFiDS_Ok);
    try
    {
      OCC_CATCH_SIGNALS
      PerformSetOfSurf(itel.ChangeValue());
    }
    catch (Standard_Failure const& anException)
    {
#ifdef OCCT_DEBUG
      std::cout << "EXCEPTION Stripe compute " << anException << std::endl;
#endif
      (void)anException;
      badstripes.Append(itel.Value());
      done = true;
      if (itel.Value()->Spine()->ErrorStatus() == ChFiDS_Ok)
      {
        itel.Value()->Spine()->SetErrorStatus(ChFiDS_Error);
      }
    }
    if (!done)
    {
      badstripes.Append(itel.Value());
    }
    done = true;
  }
  done = (badstripes.IsEmpty());

#ifdef OCCT_DEBUG // perf
  ChFi3d_ResultChron(cl_perfsetofsurf, t_perfsetofsurf);
  ChFi3d_InitChron(cl_perffilletonvertex);
#endif

  // construct fillets on each vertex + feed the Ds
  if (done)
  {
    int j;
    for (j = 1; j <= myVDataMap.Extent(); j++)
    {
      bool isPartial = false;
      try
      {
        OCC_CATCH_SIGNALS
        const bool hadPartial = hasresult;
        const int  aNbShapes  = DStr.NbShapes();
        PerformFilletOnVertex(j);
        // a corner that keeps only a partial result has failed too: the
        // computation ends without a shape, and this is the vertex
        isPartial = hasresult && !hadPartial;
        // a corner that returns with a stripe's end at it given no point has
        // failed: the DS would be read at point 0 below, out of its sight
        if (!StripeEndsHavePoints(myVDataMap(j), myVDataMap.FindKey(j)))
        {
          throw Standard_Failure("A corner left a stripe's end without its points");
        }
        // nor one that put a null shape in the DS: the topological build
        // dereferences every shape there (two stripes' corner of
        // issue273_Fillet001, edges 38 and 39 at 0.3; a crash as soon as
        // another corner of the fillet was mended)
        for (int k = aNbShapes + 1; k <= DStr.NbShapes(); k++)
        {
          if (DStr.Shape(k).IsNull())
          {
            throw Standard_Failure("A corner put a null shape in the DS");
          }
        }
      }
      catch (Standard_Failure const& anException)
      {
#ifdef OCCT_DEBUG
        std::cout << "EXCEPTION Corner compute " << anException << std::endl;
#endif
        (void)anException;
        badvertices.Append(myVDataMap.FindKey(j));
        hasresult = false;
        done      = true;
        isPartial = false;
      }
      if (!done || isPartial)
      {
        badvertices.Append(myVDataMap.FindKey(j));
      }
      done = true;
    }
    if (!hasresult)
    {
      done = badvertices.IsEmpty();
    }
  }

  // An end OnSame at four sharp edges is one PerformExtremity found a corner
  // of three across a tangent split. Whether the corner can be made of it
  // depends on the radius; where it failed, the walk runs on past the end as
  // at a break point, as it did before, and the end's state is put back for
  // the next computation.
  if (!done && !badvertices.IsEmpty())
  {
    NCollection_List<occ::handle<ChFiDS_Spine>> aFirst, aLast;
    for (itel.Initialize(myListStripe); itel.More(); itel.Next())
    {
      const occ::handle<ChFiDS_Spine>& aSp = itel.Value()->Spine();
      if (aSp.IsNull())
      {
        continue;
      }
      for (NCollection_List<TopoDS_Shape>::Iterator itV(badvertices); itV.More(); itV.Next())
      {
        const TopoDS_Vertex& aV = TopoDS::Vertex(itV.Value());
        const bool isFirst = aSp->FirstStatus() == ChFiDS_OnSame && aSp->FirstVertex().IsSame(aV);
        const bool isLast  = aSp->LastStatus() == ChFiDS_OnSame && aSp->LastVertex().IsSame(aV);
        if (!isFirst && !isLast)
        {
          continue;
        }
        // the corner's own count, which can fail as the corner did
        int aNbSharp = 0;
        try
        {
          OCC_CATCH_SIGNALS
          aNbSharp = ChFi3d_NumberOfSharpEdges(aV, myVEMap, myEFMap);
        }
        catch (Standard_Failure const&)
        {
          continue;
        }
        if (aNbSharp != 4)
        {
          continue;
        }
        if (isFirst)
        {
          aSp->SetFirstStatus(ChFiDS_BreakPoint);
          aFirst.Append(aSp);
        }
        if (isLast)
        {
          aSp->SetLastStatus(ChFiDS_BreakPoint);
          aLast.Append(aSp);
        }
      }
    }
    if (!aFirst.IsEmpty() || !aLast.IsEmpty())
    {
      // the stripes this computation left without a spine go first:
      // Compute reads every spine before its own Reset
      Reset();
      Compute();
      for (NCollection_List<occ::handle<ChFiDS_Spine>>::Iterator itS(aFirst); itS.More();
           itS.Next())
      {
        itS.Value()->SetFirstStatus(ChFiDS_OnSame);
      }
      for (NCollection_List<occ::handle<ChFiDS_Spine>>::Iterator itS(aLast); itS.More();
           itS.Next())
      {
        itS.Value()->SetLastStatus(ChFiDS_OnSame);
      }
      return;
    }
  }

#ifdef OCCT_DEBUG // perf
  ChFi3d_ResultChron(cl_perffilletonvertex, t_perffilletonvertex);
  ChFi3d_InitChron(cl_filds);
#endif

  NCollection_Map<int> MapIndSo;
  TopExp_Explorer      expso(myShape, TopAbs_SOLID);
  for (; expso.More(); expso.Next())
  {
    const TopoDS_Shape& cursol    = expso.Current();
    int                 indcursol = DStr.AddShape(cursol);
    MapIndSo.Add(indcursol);
  }
  TopExp_Explorer expsh(myShape, TopAbs_SHELL, TopAbs_SOLID);
  for (; expsh.More(); expsh.Next())
  {
    const TopoDS_Shape& cursh    = expsh.Current();
    int                 indcursh = DStr.AddShape(cursh);
    MapIndSo.Add(indcursh);
  }
  if (done)
  {
    int i1;
    for (itel.Initialize(myListStripe), i1 = 0; itel.More(); itel.Next(), i1++)
    {
      const occ::handle<ChFiDS_Stripe>& st = itel.Value();
      // 05/02/02 akm vvv : (OCC119) First we'll check ain't there
      //                    intersections between fillets
      NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator itel1;
      int                                                    i2;
      for (itel1.Initialize(myListStripe), i2 = 0; itel1.More(); itel1.Next(), i2++)
      {
        if (i2 <= i1)
        {
          // Do not twice intersect the stripes
          continue;
        }
        occ::handle<ChFiDS_Stripe> aCheckStripe = itel1.Value();
        try
        {
          OCC_CATCH_SIGNALS
          ChFi3d_StripeEdgeInter(st, aCheckStripe, DStr, tol2d);
        }
        catch (Standard_Failure const& anException)
        {
#ifdef OCCT_DEBUG
          std::cout << "EXCEPTION Fillets compute " << anException << std::endl;
#endif
          (void)anException;
          badstripes.Append(itel.Value());
          hasresult = false;
          done      = false;
          break;
        }
      }
      // 05/02/02 akm ^^^
      int solidindex = st->SolidIndex();
      ChFi3d_FilDS(solidindex, st, DStr, myRegul, tolapp3d, tol2d);
      if (!done)
      {
        break;
      }
    }

#ifdef OCCT_DEBUG // perf
    ChFi3d_ResultChron(cl_filds, t_filds);
    ChFi3d_InitChron(cl_reconstruction);
#endif

    if (done)
    {
      BRep_Builder B1;
      CompleteDS(DStr, myShape);
      // Update tolerances on vertex to max adjacent edges or
      // Update tolerances on degenerated edge to max of adjacent vertexes.
      TopOpeBRepDS_CurveExplorer cex(DStr);
      for (; cex.More(); cex.Next())
      {
        TopOpeBRepDS_Curve& c     = *((TopOpeBRepDS_Curve*)(void*)&(cex.Curve()));
        double              tolc  = 0.;
        bool                degen = c.Curve().IsNull();
        if (!degen)
        {
          tolc = c.Tolerance();
        }
        int                        ic = cex.Index();
        TopOpeBRepDS_PointIterator It(myDS->CurvePoints(ic));
        for (; It.More(); It.Next())
        {
          occ::handle<TopOpeBRepDS_CurvePointInterference> II;
          II = occ::down_cast<TopOpeBRepDS_CurvePointInterference>(It.Value());
          if (II.IsNull())
          {
            continue;
          }
          TopOpeBRepDS_Kind gk = II->GeometryType();
          int               gi = II->Geometry();
          if (gk == TopOpeBRepDS_VERTEX)
          {
            const TopoDS_Vertex& v    = TopoDS::Vertex(myDS->Shape(gi));
            double               tolv = ChFi3d_VertexTolerance(v);
            if (tolv > 0.0001)
            {
              tolv += 0.0003;
              if (tolc < tolv)
              {
                tolc = tolv + 0.00001;
              }
            }
            // the vertex covers the curve's end, which may lie a hair farther
            // from it than the curve's tolerance says (a corner's extension
            // starting off the vertex): BRepCheck would have it outside. Not
            // more than a hair: an end farther off is a fault to show.
            double tolcv = tolc;
            if (!degen)
            {
              const double aGap =
                c.Curve()->Value(II->Parameter()).Distance(BRep_Tool::Pnt(v));
              if (aGap > tolcv && aGap <= tolcv + Precision::Confusion())
              {
                tolcv = aGap * (1. + 1.e-9); // and against rounding
              }
            }
            if (degen && tolc < tolv)
            {
              tolc = tolv;
            }
            else if (tolcv > tolv)
            {
              ChFi3d_UpdateVertex(v, tolcv);
            }
          }
          else if (gk == TopOpeBRepDS_POINT)
          {
            TopOpeBRepDS_Point& p    = DStr.ChangePoint(gi);
            double              tolp = p.Tolerance();
            if (degen && tolc < tolp)
            {
              tolc = tolp;
            }
            else if (tolc > tolp)
            {
              p.Tolerance(tolc);
            }
          }
        }
        if (degen)
        {
          c.Tolerance(tolc);
        }
      }
      myCoup->Perform(myDS);
      NCollection_Map<int>::Iterator It(MapIndSo);
      for (; It.More(); It.Next())
      {
        int                 indsol   = It.Key();
        const TopoDS_Shape& curshape = DStr.Shape(indsol);
        myCoup->MergeSolid(curshape, TopAbs_IN);
      }

      int i = 1, n = DStr.NbShapes();
      for (; i <= n; i++)
      {
        const TopoDS_Shape S = DStr.Shape(i);
        if (S.ShapeType() != TopAbs_EDGE)
        {
          continue;
        }
        bool issplitIN = myCoup->IsSplit(S, TopAbs_IN);
        if (!issplitIN)
        {
          continue;
        }
        NCollection_List<TopoDS_Shape>::Iterator it(myCoup->Splits(S, TopAbs_IN));
        for (; it.More(); it.Next())
        {
          const TopoDS_Edge& newE = TopoDS::Edge(it.Value());
          double             tole = BRep_Tool::Tolerance(newE);
          TopExp_Explorer    exv(newE, TopAbs_VERTEX);
          for (; exv.More(); exv.Next())
          {
            const TopoDS_Vertex& v    = TopoDS::Vertex(exv.Current());
            double               tolv = ChFi3d_VertexTolerance(v);
            if (tole > tolv)
            {
              ChFi3d_UpdateVertex(v, tole);
            }
          }
        }
      }
      if (!hasresult)
      {
        B1.MakeCompound(TopoDS::Compound(myShapeResult));
        for (It = NCollection_Map<int>::Iterator(MapIndSo); It.More(); It.Next())
        {
          int                                      indsol   = It.Key();
          const TopoDS_Shape&                      curshape = DStr.Shape(indsol);
          NCollection_List<TopoDS_Shape>::Iterator its      = myCoup->Merged(curshape, TopAbs_IN);
          if (!its.More())
          {
            B1.Add(myShapeResult, curshape);
          }
          else
          {
            // If the old type of Shape is Shell, Shell is placed instead of Solid,
            // However there is a problem for compound of open Shell.
            while (its.More())
            {
              const TopAbs_ShapeEnum letype = curshape.ShapeType();
              if (letype == TopAbs_SHELL)
              {
                TopExp_Explorer     expsh2(its.Value(), TopAbs_SHELL);
                const TopoDS_Shape& cursh = expsh2.Current();
                B1.Add(myShapeResult, cursh);
                its.Next();
              }
              else
              {
                B1.Add(myShapeResult, its.Value());
                its.Next();
              }
            }
          }
        }
        // the changes to the input's frozen vertices so far, in the result
        // before the passes below read its tolerances
        aThaw.Apply(myShapeResult);
        ChFi3d_EdgesCoverTheirEnds(myShapeResult, myShape);
        // the builder's history, which Generated and BRepFilletAPI read,
        // as the result has its shapes now
        occ::handle<BRepTools_History> aTouch =
          ChFi3d_DropMicroEdges(myShapeResult, myShape, tolapp3d);
        // until no face is left pinched: a face is split at one point a time
        for (int aPass = 0; aPass < 16; aPass++)
        {
          const occ::handle<BRepTools_History> aSplit =
            ChFi3d_SplitPinchedFaces(myShapeResult, myShape);
          if (aSplit.IsNull())
          {
            break;
          }
          if (aTouch.IsNull())
          {
            aTouch = aSplit;
          }
          else
          {
            aTouch->Merge(aSplit);
          }
        }
        // a face and the one the run-on made on its surface beside it, one;
        // a closed side and the face its surface carried up made at a pinch
        for (int k = 0; k < 2; k++)
        {
          const occ::handle<BRepTools_History> aSeams =
            k == 0 ? ChFi3d_MergeRunOnSeams(myShapeResult, tolapp3d)
                   : ChFi3d_MergeSideSeams(myShapeResult, tolapp3d);
          if (aSeams.IsNull())
          {
            continue;
          }
          if (aTouch.IsNull())
          {
            aTouch = aSeams;
          }
          else
          {
            aTouch->Merge(aSeams);
          }
        }
        // the radii a hair smaller (ChFi3d_SnapSlivers): the slivers of faces
        // as wide as the fillets beside them, gone
        if (ChFi3d_SnapSlivers())
        {
          const occ::handle<BRepTools_History> aSlivers =
            ChFi3d_CollapseSlivers(myShapeResult, myShape, tolapp3d);
          if (!aSlivers.IsNull())
          {
            if (aTouch.IsNull())
            {
              aTouch = aSlivers;
            }
            else
            {
              aTouch->Merge(aSlivers);
            }
          }
        }
        if (!aTouch.IsNull())
        {
          for (int iS = 1; iS <= DStr.NbShapes(); iS++)
          {
            for (TopAbs_State aSt : {TopAbs_IN, TopAbs_OUT, TopAbs_ON})
            {
              if (myCoup->IsSplit(DStr.Shape(iS), aSt))
              {
                ChFi3d_Remap(myCoup->ChangeBuilder().ChangeSplit(DStr.Shape(iS), aSt), aTouch);
              }
            }
          }
          for (NCollection_DataMap<TopoDS_Shape, NCollection_List<int>,
                                   TopTools_ShapeMapHasher>::Iterator itEV(myEVIMap);
               itEV.More();
               itEV.Next())
          {
            for (NCollection_List<int>::Iterator itI(itEV.Value()); itI.More(); itI.Next())
            {
              // NewFaces is the builder's own list, read-only only through
              // the accessor
              ChFi3d_Remap(
                const_cast<NCollection_List<TopoDS_Shape>&>(myCoup->NewFaces(itI.Value())),
                aTouch);
            }
          }
        }
      }
      else
      {
        done = false;
        B1.MakeCompound(TopoDS::Compound(badShape));
        for (It = NCollection_Map<int>::Iterator(MapIndSo); It.More(); It.Next())
        {
          int                                      indsol   = It.Key();
          const TopoDS_Shape&                      curshape = DStr.Shape(indsol);
          NCollection_List<TopoDS_Shape>::Iterator its      = myCoup->Merged(curshape, TopAbs_IN);
          if (!its.More())
          {
            B1.Add(badShape, curshape);
          }
          else
          {
            while (its.More())
            {
              B1.Add(badShape, its.Value());
              its.Next();
            }
          }
        }
      }
#ifdef OCCT_DEBUG // perf
      ChFi3d_ResultChron(cl_reconstruction, t_reconstruction);
      ChFi3d_InitChron(cl_setregul);
#endif

      // Regularities are coded after cutting.
      SetRegul();
      if (done)
      {
        aThaw.Apply(myShapeResult);
      }

#ifdef OCCT_DEBUG // perf
      ChFi3d_ResultChron(cl_setregul, t_setregul);
#endif
    }
  }
#ifdef OCCT_DEBUG // perf
  ChFi3d_ResultChron(cl_total, t_total);
#endif

  // display of time for perfs

#ifdef OCCT_DEBUG
  if (ChFi3d_GettraceCHRON())
  {
    std::cout << std::endl;
    std::cout << "COMPUTE: temps total " << t_total << "s  dont :" << std::endl;
    std::cout << "- Init + ExtentAnalyse " << t_extent << "s" << std::endl;
    std::cout << "- PerformSetOfSurf " << t_perfsetofsurf << "s" << std::endl;
    std::cout << "- PerformFilletOnVertex " << t_perffilletonvertex << "s" << std::endl;
    std::cout << "- FilDS " << t_filds << "s" << std::endl;
    std::cout << "- Reconstruction " << t_reconstruction << "s" << std::endl;
    std::cout << "- SetRegul " << t_setregul << "s" << std::endl << std::endl;

    std::cout << std::endl;
    std::cout << "temps PERFORMSETOFSURF " << t_perfsetofsurf << "s  dont : " << std::endl;
    std::cout << "- SetofKPart " << t_perfsetofkpart << "s" << std::endl;
    std::cout << "- SetofKGen " << t_perfsetofkgen << "s" << std::endl;
    std::cout << "- MakeExtremities " << t_makextremities << "s" << std::endl << std::endl;

    std::cout << "temps SETOFKGEN " << t_perfsetofkgen << "s dont : " << std::endl;
    std::cout << "- PerformSurf " << t_performsurf << "s" << std::endl;
    std::cout << "- starsol " << t_startsol << "s" << std::endl << std::endl;

    std::cout << "temps PERFORMSURF " << t_performsurf << "s  dont : " << std::endl;
    std::cout << "- computedata " << t_computedata << "s" << std::endl;
    std::cout << "- completedata " << t_completedata << "s" << std::endl << std::endl;

    std::cout << "temps PERFORMFILLETVERTEX " << t_perffilletonvertex << "s dont : " << std::endl;
    std::cout << "- PerformOneCorner " << t_perform1corner << "s" << std::endl;
    std::cout << "- PerformIntersectionAtEnd " << t_performatend << "s" << std::endl;
    std::cout << "- PerformTwoCorner " << t_perform2corner << "s" << std::endl;
    std::cout << "- PerformThreeCorner " << t_perform3corner << "s" << std::endl;
    std::cout << "- PerformMoreThreeCorner " << t_performmore3corner << "s" << std::endl
              << std::endl;

    std::cout << "temps PerformOneCorner " << t_perform1corner << "s dont:" << std::endl;
    std::cout << "- temps condition if (same) " << t_same << "s " << std::endl;
    std::cout << "- temps condition if (inter) " << t_inter << "s " << std::endl;
    std::cout << "- temps condition if (same inter) " << t_sameinter << "s " << std::endl
              << std::endl;

    std::cout << "temps PerformTwocorner " << t_perform2corner << "s  dont:" << std::endl;
    std::cout << "- temps initialisation " << t_t2cornerinit << "s" << std::endl;
    std::cout << "- temps PerformTwoCornerbyInter " << t_perf2cornerbyinter << "s" << std::endl;
    std::cout << "- temps ChFiKPart_ComputeData " << t_chfikpartcompdata << "s" << std::endl;
    std::cout << "- temps cheminement " << t_cheminement << "s" << std::endl;
    std::cout << "- temps remplissage " << t_remplissage << "s" << std::endl;
    std::cout << "- temps mise a jour stripes  " << t_t2cornerDS << "s" << std::endl << std::endl;

    std::cout << " temps PerformThreecorner " << t_perform3corner << "s  dont:" << std::endl;
    std::cout << "- temps initialisation " << t_t3cornerinit << "s" << std::endl;
    std::cout << "- temps cas spherique  " << t_spherique << "s" << std::endl;
    std::cout << "- temps cas torique  " << t_torique << "s" << std::endl;
    std::cout << "- temps notfilling " << t_notfilling << "s" << std::endl;
    std::cout << "- temps filling " << t_filling << "s" << std::endl;
    std::cout << "- temps mise a jour stripes  " << t_t3cornerDS << "s" << std::endl << std::endl;

    std::cout << "temps PerformMore3Corner " << t_performmore3corner << "s dont:" << std::endl;
    std::cout << "-temps plate " << t_plate << "s " << std::endl;
    std::cout << "-temps approxplate " << t_approxplate << "s " << std::endl;
    std::cout << "-temps batten " << t_batten << "s " << std::endl << std::endl;

    std::cout << "TEMPS DIVERS " << std::endl;
    std::cout << "-temps ChFi3d_sameparameter " << t_sameparam << "s" << std::endl << std::endl;
  }
#endif
  //
  // Inspect the new faces to provide sameparameter
  // if it is necessary
  if (IsDone())
  {
    double                                   SameParTol = Precision::Confusion();
    int                                      aNbSurfaces, iF;
    NCollection_List<TopoDS_Shape>::Iterator aIt;
    //
    aNbSurfaces = myDS->NbSurfaces();

    NCollection_Map<const TopoDS_TShape*> aSpread;
    for (iF = 1; iF <= aNbSurfaces; ++iF)
    {
      const NCollection_List<TopoDS_Shape>& aLF = myCoup->NewFaces(iF);
      aIt.Initialize(aLF);
      for (; aIt.More(); aIt.Next())
      {
        const TopoDS_Shape& aF = aIt.Value();
        BRepLib::SameParameter(aF, SameParTol, true);
        ChFi3d_SpreadThawed(myShapeResult, aF, aSpread);
        ShapeFix::SameParameter(aF, false, SameParTol);
      }
    }
    // ShapeFix's pass raises a vertex to its edges' curve ends, not the
    // edges with it: #273's fillet of 0.3 on its edge 5 left a circle
    // ending 0.029 off its vertex at tolerance 1e-7 -- covered before
    // only while the corner's face tolerance, 0.045, had raised the vertex
    // past it in time for the pass above
    ChFi3d_EdgesCoverTheirEnds(myShapeResult, myShape);
    ChFi3d_EdgesCoverTheirPCurves(myShapeResult, myShape);
  }
}

//=======================================================================
// function : PerformSingularCorner
// purpose  : Load vertex and degenerated edges.
//=======================================================================

void ChFi3d_Builder::PerformSingularCorner(const int Index)
{
  NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator It;
  occ::handle<ChFiDS_Stripe>                             stripe;
  TopOpeBRepDS_DataStructure&                            DStr = myDS->ChangeDS();
  const TopoDS_Vertex&                                   Vtx  = myVDataMap.FindKey(Index);

  occ::handle<ChFiDS_SurfData> Fd;
  int                          i, Icurv;
  int                          Ivtx = 0;
  for (It.Initialize(myVDataMap(Index)), i = 0; It.More(); It.Next(), i++)
  {
    stripe = It.Value();
    // SurfData concerned and its CommonPoints,
    int  sens                     = 0;
    int  num                      = ChFi3d_IndexOfSurfData(Vtx, stripe, sens);
    bool isfirst                  = (sens == 1);
    Fd                            = stripe->SetOfSurfData()->Sequence().Value(num);
    const ChFiDS_CommonPoint& CV1 = Fd->Vertex(isfirst, 1);
    const ChFiDS_CommonPoint& CV2 = Fd->Vertex(isfirst, 2);
    // Is it always degenerated ?
    if (CV1.Point().IsEqual(CV2.Point(), 0))
    {
      // if yes the vertex is stored in the stripe
      // and the edge at end is created
      if (i == 0)
      {
        Ivtx = ChFi3d_IndexPointInDS(CV1, DStr);
      }
      double                    tolreached;
      double                    Pardeb, Parfin;
      gp_Pnt2d                  VOnS1, VOnS2;
      occ::handle<Geom_Curve>   C3d;
      occ::handle<Geom2d_Curve> PCurv;
      TopOpeBRepDS_Curve        Crv;
      if (isfirst)
      {
        VOnS1 =
          Fd->InterferenceOnS1().PCurveOnSurf()->Value(Fd->InterferenceOnS1().FirstParameter());
        VOnS2 =
          Fd->InterferenceOnS2().PCurveOnSurf()->Value(Fd->InterferenceOnS2().FirstParameter());
      }
      else
      {
        VOnS1 =
          Fd->InterferenceOnS1().PCurveOnSurf()->Value(Fd->InterferenceOnS1().LastParameter());
        VOnS2 =
          Fd->InterferenceOnS2().PCurveOnSurf()->Value(Fd->InterferenceOnS2().LastParameter());
      }

      ChFi3d_ComputeArete(CV1,
                          VOnS1,
                          CV2,
                          VOnS2,
                          DStr.Surface(Fd->Surf()).Surface(),
                          C3d,
                          PCurv,
                          Pardeb,
                          Parfin,
                          tolapp3d,
                          tolapp2d,
                          tolreached,
                          0);
      Crv   = TopOpeBRepDS_Curve(C3d, tolreached);
      Icurv = DStr.AddCurve(Crv);

      stripe->SetCurve(Icurv, isfirst);
      stripe->SetParameters(isfirst, Pardeb, Parfin);
      stripe->ChangePCurve(isfirst) = PCurv;
      stripe->SetIndexPoint(Ivtx, isfirst, 1);
      stripe->SetIndexPoint(Ivtx, isfirst, 2);
    }
  }
}

//=======================================================================
// function : LineOnSplit
// purpose  : A fillet as wide as the piece of a wall its spine is on (a wall
//           kept in coplanar pieces, a body without Refine) has its line
//           on the piece beyond, along the split: the spine's piece goes
//           under the fillet there (ChFi3d_SplitPieceOfSpine). Where the
//           fillet ends at a vertex with other fillets, the corner is made
//           on the faces at the vertex, and the far piece is not one of
//           them. The line is then shown on the spine's piece while the
//           corner is made (Show), and put back on its own after it
//           (Restore), the split cut where the line now ends (Cut).
//=======================================================================

namespace
{
struct LineOnSplit
{
  occ::handle<ChFiDS_Stripe>   Stripe;
  occ::handle<ChFiDS_SurfData> SD;
  bool                         IsFirst = false;
  int                          Ons     = 0;
  int                          IndFar  = 0;
  occ::handle<Geom2d_Curve>    PcFar;
  TopoDS_Edge                  Split;
  bool                         IsTurned = false;

  //! The parameter on the split of the point of the line at <theW>, when
  //! the point is on the split, inside its ends.
  bool OnSplit(const occ::handle<Geom_Curve>& theC,
               const double                   theW,
               const double                   theTol,
               double&                        thePar,
               double&                        theDist) const
  {
    const BRepAdaptor_Curve aSC(Split);
    Extrema_ExtPC           anExt(theC->Value(theW), aSC);
    double                  aD2 = Precision::Infinite();
    for (int j = 1; anExt.IsDone() && j <= anExt.NbExt(); j++)
    {
      if (anExt.SquareDistance(j) < aD2)
      {
        aD2    = anExt.SquareDistance(j);
        thePar = anExt.Point(j).Parameter();
      }
    }
    if (aD2 > theTol * theTol)
    {
      return false;
    }
    theDist           = std::sqrt(aD2);
    const double aRes = aSC.Resolution(theTol);
    return thePar > aSC.FirstParameter() + aRes && thePar < aSC.LastParameter() - aRes;
  }

  void SetFace(const int theInd, const occ::handle<Geom2d_Curve>& thePc) const
  {
    if (Ons == 1)
    {
      SD->ChangeIndexOfS1(theInd);
    }
    else
    {
      SD->ChangeIndexOfS2(theInd);
    }
    SD->ChangeInterference(Ons).ChangePCurveOnFace() = thePc;
    // the two pieces the other way up on their surfaces: the line's
    // transition, told by its face's surface, is the other one (the
    // stripe's own orientations are those of the spine's faces)
    if (IsTurned)
    {
      ChFiDS_FaceInterference& aFi = SD->ChangeInterference(Ons);
      aFi.SetTransition(TopAbs::Reverse(aFi.Transition()));
    }
  }

  //! The end of <theStripe> at <theVtx>, when its line on a side is such a
  //! line: it is put on the spine's piece.
  bool Show(TopOpeBRepDS_DataStructure&       theDS,
            const TopoDS_Vertex&              theVtx,
            const occ::handle<ChFiDS_Stripe>& theStripe,
            const ChFiDS_Map&                 theEFMap,
            const double                      theTol)
  {
    if (theStripe->Spine().IsNull() || theStripe->SetOfSurfData().IsNull()
        || theStripe->SetOfSurfData()->Length() == 0)
    {
      return false;
    }
    int       aSens = 0;
    const int aNum  = ChFi3d_IndexOfSurfData(theVtx, theStripe, aSens);
    Stripe          = theStripe;
    SD              = theStripe->SetOfSurfData()->Value(aNum);
    IsFirst         = aSens == 1;
    const TopoDS_Edge& aSpineE =
      theStripe->Spine()->Edges(IsFirst ? 1 : theStripe->Spine()->NbEdges());
    if (!theEFMap.Contains(aSpineE))
    {
      return false;
    }
    for (Ons = 1; Ons <= 2; Ons++)
    {
      const ChFiDS_FaceInterference& aFi = SD->Interference(Ons);
      IndFar                             = SD->Index(Ons);
      if (IndFar <= 0 || aFi.LineIndex() <= 0 || aFi.PCurveOnFace().IsNull()
          || theDS.Shape(IndFar).ShapeType() != TopAbs_FACE)
      {
        continue;
      }
      const occ::handle<Geom_Curve>& aC = theDS.Curve(aFi.LineIndex()).Curve();
      if (aC.IsNull())
      {
        continue;
      }
      // not a face of the spine, nor one at the vertex
      const TopoDS_Shape& aFar    = theDS.Shape(IndFar);
      bool                isOther = true;
      for (NCollection_List<TopoDS_Shape>::Iterator itF(theEFMap(aSpineE)); itF.More(); itF.Next())
      {
        isOther = isOther && !itF.Value().IsSame(aFar);
      }
      for (TopExp_Explorer exV(aFar, TopAbs_VERTEX); exV.More() && isOther; exV.Next())
      {
        isOther = !exV.Current().IsSame(theVtx);
      }
      if (!isOther)
      {
        continue;
      }
      // the spine's piece across an edge of it, the same way up
      for (TopExp_Explorer exE(aFar, TopAbs_EDGE); exE.More(); exE.Next())
      {
        const TopoDS_Edge& anE = TopoDS::Edge(exE.Current());
        if (BRep_Tool::Degenerated(anE) || !theEFMap.Contains(anE)
            || theEFMap(anE).Extent() != 2)
        {
          continue;
        }
        TopoDS_Face aFarSh, aNear;
        for (NCollection_List<TopoDS_Shape>::Iterator itF(theEFMap(anE)); itF.More(); itF.Next())
        {
          if (itF.Value().IsSame(aFar))
          {
            aFarSh = TopoDS::Face(itF.Value());
            continue;
          }
          for (NCollection_List<TopoDS_Shape>::Iterator itS(theEFMap(aSpineE)); itS.More();
               itS.Next())
          {
            if (itS.Value().IsSame(itF.Value()))
            {
              aNear = TopoDS::Face(itF.Value());
            }
          }
        }
        if (aFarSh.IsNull() || aNear.IsNull() || !ChFi3d_CoplanarPieces(aFarSh, aNear, theTol))
        {
          continue;
        }
        // the line along it, from its end here on
        Split               = anE;
        const double aW     = aFi.Parameter(IsFirst);
        const double aWIn   = aW + 0.25 * (aFi.Parameter(!IsFirst) - aW);
        double       aPar = 0., aDist = 0.;
        if (!OnSplit(aC, aW, theTol, aPar, aDist) || !OnSplit(aC, aWIn, theTol, aPar, aDist))
        {
          continue;
        }
        const occ::handle<Geom2d_Curve> aPc = GeomProjLib::Curve2d(aC, BRep_Tool::Surface(aNear));
        if (aPc.IsNull())
        {
          continue;
        }
        PcFar    = aFi.PCurveOnFace();
        IsTurned = aFarSh.Orientation() != aNear.Orientation();
        SetFace(theDS.AddShape(aNear), aPc);
        return true;
      }
    }
    return false;
  }

  //! The line back on its piece.
  void Restore() const { SetFace(IndFar, PcFar); }

  //! The split cut at the line's end, its part under the line gone.
  void Cut(TopOpeBRepDS_DataStructure& theDS, const double theTol) const
  {
    const ChFiDS_FaceInterference& aFi = SD->Interference(Ons);
    const ChFiDS_CommonPoint&      aCP = SD->Vertex(IsFirst, Ons);
    if (aCP.IsVertex() || aCP.IsOnArc())
    {
      return; // FILDS cuts the arc
    }
    const int                      anIndP = Stripe->IndexPoint(IsFirst, Ons);
    const occ::handle<Geom_Curve>& aC     = theDS.Curve(aFi.LineIndex()).Curve();
    const double                   aW     = aFi.Parameter(IsFirst);
    const double                   aWIn   = aW + 0.25 * (aFi.Parameter(!IsFirst) - aW);
    double                         aPar = 0., aParIn = 0., aDist = 0., aDistIn = 0.;
    if (anIndP <= 0 || aC.IsNull() || !OnSplit(aC, aW, theTol, aPar, aDist)
        || !OnSplit(aC, aWIn, theTol, aParIn, aDistIn))
    {
      throw Standard_Failure("A corner took a line along a split of its wall off the split");
    }
    const int IS = theDS.AddShape(Split);
    for (NCollection_List<occ::handle<TopOpeBRepDS_Interference>>::Iterator it(
           theDS.ShapeInterferences(IS));
         it.More();
         it.Next())
    {
      if (it.Value()->GeometryType() == TopOpeBRepDS_POINT && it.Value()->Geometry() == anIndP)
      {
        return;
      }
    }
    TopOpeBRepDS_Point& aP = theDS.ChangePoint(anIndP);
    aP.Tolerance(std::max(aP.Tolerance(), aDist));
    theDS.ChangeShapeInterferences(IS).Append(
      ChFi3d_FilPointInDS(aParIn > aPar ? TopAbs_REVERSED : TopAbs_FORWARD, IS, anIndP, aPar));
  }
};

//! The lines shown on the spine's piece at a vertex, put back when the
//! corner is made or fails.
class LinesOnSplit
{
public:
  LinesOnSplit(TopOpeBRepDS_DataStructure& theDS, const double theTol)
      : myDS(theDS),
        myTol(theTol)
  {
  }

  ~LinesOnSplit()
  {
    for (NCollection_List<LineOnSplit>::Iterator it(myLines); it.More(); it.Next())
    {
      it.Value().Restore();
    }
  }

  void Show(const TopoDS_Vertex&              theVtx,
            const occ::handle<ChFiDS_Stripe>& theStripe,
            const ChFiDS_Map&                 theEFMap)
  {
    LineOnSplit aLine;
    if (aLine.Show(myDS, theVtx, theStripe, theEFMap, myTol))
    {
      myLines.Append(aLine);
    }
  }

  //! The corner is made.
  void Finish()
  {
    NCollection_List<LineOnSplit> aLines;
    aLines.Assign(myLines);
    myLines.Clear();
    for (NCollection_List<LineOnSplit>::Iterator it(aLines); it.More(); it.Next())
    {
      it.Value().Restore();
    }
    for (NCollection_List<LineOnSplit>::Iterator it(aLines); it.More(); it.Next())
    {
      it.Value().Cut(myDS, myTol);
    }
  }

private:
  TopOpeBRepDS_DataStructure&   myDS;
  double                        myTol;
  NCollection_List<LineOnSplit> myLines;
};

//=======================================================================
// function : PinchEnds
// purpose  : A closed fillet whose two ends stop short of one vertex: the
//           foot of a boss whose wall touches a side of the plate it
//           stands on along a line through the vertex (a cylinder flush
//           with a plane side, or set against a round plate's rim from
//           inside). The foot's fillet runs round the boss until its line
//           on the plate meets the plate's edge along the side, on either
//           side of the vertex; its line on the boss does not get there,
//           the boss meeting the side's edge only at the vertex. Nothing
//           makes a corner of one stripe ending twice at a vertex, and the
//           fillet was refused.
//           The side is carried up as it is: each end of the fillet is cut
//           by the side's surface, from its point on the plate's edge to
//           where its line on the boss, run on, touches the side (on the
//           boss's line of contact with it), and the two cuts meet there.
//           The side grows over the gap, the boss's foot is the fillet's
//           line all round, closed at that point, and the plate's edge
//           between the two ends goes with the vertex -- the fillet made on
//           a plate reaching past the side, then cut by it.
//           Made where the fillet's surface runs on exactly past its end (a
//           periodic one: the torus round a cylinder on a plane) and the
//           side's surface is an elementary one; refused, as before,
//           otherwise.
//=======================================================================

struct PinchEnd
{
  occ::handle<ChFiDS_Stripe>   St;
  occ::handle<ChFiDS_SurfData> SD;
  bool                         IsFirst = false;
  int                          OnsArc  = 0;
  int                          OnsAir  = 0;
  TopoDS_Face                  Fv;
  double                       WEnd = 0., WAir = 0.;
  gp_Pnt                       PAir;
  occ::handle<Geom_Curve>      Cut;
  occ::handle<Geom2d_Curve>    PcS, PcV;
  double                       Tol = 0.;
};

//! The signed distance of <theP> to the surface <theS>, along its normal.
double SignedDistance(const GeomAdaptor_Surface& theS, const gp_Pnt& theP, gp_Pnt2d& theUV)
{
  if (theS.GetType() == GeomAbs_Plane)
  {
    const gp_Pln aPln = theS.Plane();
    double       aU = 0., aV = 0.;
    ElSLib::Parameters(aPln, theP, aU, aV);
    theUV.SetCoord(aU, aV);
    return gp_Vec(aPln.Location(), theP).Dot(gp_Vec(aPln.Axis().Direction()));
  }
  Extrema_ExtPS anExt(theP, theS, Precision::PConfusion(), Precision::PConfusion());
  if (!anExt.IsDone() || anExt.NbExt() == 0)
  {
    return Precision::Infinite();
  }
  int    iMin = 1;
  double aD2  = anExt.SquareDistance(1);
  for (int i = 2; i <= anExt.NbExt(); i++)
  {
    if (anExt.SquareDistance(i) < aD2)
    {
      aD2  = anExt.SquareDistance(i);
      iMin = i;
    }
  }
  double aU = 0., aV = 0.;
  anExt.Point(iMin).Parameter(aU, aV);
  theUV.SetCoord(aU, aV);
  gp_Pnt aPS;
  gp_Vec aDU, aDV;
  theS.D1(aU, aV, aPS, aDU, aDV);
  const gp_Vec aN = aDU.Crossed(aDV);
  const double aD = std::sqrt(aD2);
  return gp_Vec(aPS, theP).Dot(aN) < 0. ? -aD : aD;
}

bool HasVertex(const TopoDS_Edge& theE, const TopoDS_Vertex& theV)
{
  for (TopExp_Explorer ex(theE, TopAbs_VERTEX); ex.More(); ex.Next())
  {
    if (ex.Current().IsSame(theV))
    {
      return true;
    }
  }
  return false;
}

//! The point of the end's line on the boss where, run on past the end, it
//! touches <theSv>, the side's surface; false if it does not within half
//! a turn.
bool LineTouching(const occ::handle<Geom_Surface>& theS,
                         const ChFiDS_FaceInterference&   theFi,
                         const bool                       theFirst,
                         const GeomAdaptor_Surface&       theSv,
                         const double                     theTol,
                         double&                          theW,
                         gp_Pnt&                          theP)
{
  const occ::handle<Geom2d_Curve>& aPc = theFi.PCurveOnSurf();
  if (aPc.IsNull())
  {
    return false;
  }
  const double aW0  = theFi.Parameter(theFirst);
  const double aDir = theFirst ? -1. : 1.;
  // half a turn of the surface, along the line
  gp_Vec2d aT;
  gp_Pnt2d aUV0;
  aPc->D1(aW0, aUV0, aT);
  if (!theS->IsUPeriodic() || std::abs(aT.Y()) > 1.e-9 * aT.Magnitude()
      || aT.Magnitude() <= gp::Resolution())
  {
    return false;
  }
  const double aSpan = 0.5 * theS->UPeriod() / std::abs(aT.X());
  auto         aDist = [&](const double theWt) {
    const gp_Pnt2d aUV = aPc->Value(theWt);
    gp_Pnt2d       aUVv;
    return SignedDistance(theSv, theS->Value(aUV.X(), aUV.Y()), aUVv);
  };
  const double aD0 = aDist(aW0);
  if (std::abs(aD0) <= theTol || Precision::IsInfinite(aD0))
  {
    return false;
  }
  // the first local minimum of the distance on the far side, the line
  // touching the side there (crossing it is not this corner)
  const int NbS  = 400;
  double    aPrv = std::abs(aD0), aWMin = aW0;
  int       iMin = 0;
  for (int i = 1; i <= NbS; i++)
  {
    const double aW = aW0 + aDir * aSpan * i / NbS;
    const double aD = aDist(aW);
    if (Precision::IsInfinite(aD) || aD * aD0 < -theTol * theTol)
    {
      return false;
    }
    if (std::abs(aD) < aPrv)
    {
      aPrv  = std::abs(aD);
      aWMin = aW;
      iMin  = i;
    }
    else if (iMin > 0)
    {
      break;
    }
  }
  if (iMin == 0 || iMin == NbS)
  {
    return false;
  }
  // golden section on the step either side
  const double aStep = aSpan / NbS;
  double       a = aWMin - aStep, b = aWMin + aStep;
  const double aG = 0.5 * (std::sqrt(5.) - 1.);
  double       c = b - aG * (b - a), d = a + aG * (b - a);
  double       fc = std::abs(aDist(c)), fd = std::abs(aDist(d));
  for (int it = 0; it < 100 && b - a > Precision::PConfusion() * 1.e-3; it++)
  {
    if (fc < fd)
    {
      b  = d;
      d  = c;
      fd = fc;
      c  = b - aG * (b - a);
      fc = std::abs(aDist(c));
    }
    else
    {
      a  = c;
      c  = d;
      fc = fd;
      d  = a + aG * (b - a);
      fd = std::abs(aDist(d));
    }
  }
  theW                = 0.5 * (a + b);
  const gp_Pnt2d aUV  = aPc->Value(theW);
  theP                = theS->Value(aUV.X(), aUV.Y());
  return std::abs(aDist(theW)) <= theTol;
}

//! The cut of the end by the side's surface, from the point on the arc to
//! the line's touching point, its pcurves on the fillet's surface and the
//! side's, and the tolerance it reaches. <theToLine>: the cut tangent to the
//! line at that point.
bool CutEnd(const occ::handle<Geom_Surface>& theS,
                   const GeomAdaptor_Surface&       theSv,
                   const gp_Pnt2d&                  theUVvArc,
                   PinchEnd&                        theEnd,
                   const bool                       theToLine,
                   const double                     theTol)
{
  const ChFiDS_FaceInterference& aFiArc = theEnd.SD->Interference(theEnd.OnsArc);
  const ChFiDS_FaceInterference& aFiAir = theEnd.SD->Interference(theEnd.OnsAir);
  const occ::handle<Geom2d_Curve>& aPcArc = aFiArc.PCurveOnSurf();
  const occ::handle<Geom2d_Curve>& aPcAir = aFiAir.PCurveOnSurf();
  if (aPcArc.IsNull() || aPcAir.IsNull())
  {
    return false;
  }
  // the point of the surface at <theW> along the spine, <theT> across
  // from the arc's side (0) to the line's (1)
  auto aUV = [&](const double theW, const double theT) {
    const gp_Pnt2d a1 = aPcArc->Value(theW), a2 = aPcAir->Value(theW);
    return gp_Pnt2d((1. - theT) * a1.X() + theT * a2.X(), (1. - theT) * a1.Y() + theT * a2.Y());
  };
  auto aDist = [&](const double theW, const double theT, gp_Pnt2d& theUVv) {
    const gp_Pnt2d aP = aUV(theW, theT);
    return SignedDistance(theSv, theS->Value(aP.X(), aP.Y()), theUVv);
  };
  const double aW0 = theEnd.WEnd, aW1 = theEnd.WAir;
  gp_Pnt2d     aBid;
  const double aSideIn = aDist(aW0, 1., aBid);
  if (std::abs(aSideIn) <= theTol)
  {
    return false;
  }
  const int                                 NbP = 41;
  occ::handle<NCollection_HArray1<gp_Pnt>>   aPnts = new NCollection_HArray1<gp_Pnt>(1, NbP);
  occ::handle<NCollection_HArray1<gp_Pnt2d>> aUVs  = new NCollection_HArray1<gp_Pnt2d>(1, NbP);
  occ::handle<NCollection_HArray1<gp_Pnt2d>> aUVvs = new NCollection_HArray1<gp_Pnt2d>(1, NbP);
  occ::handle<NCollection_HArray1<double>>   aPars = new NCollection_HArray1<double>(1, NbP);
  gp_Pnt2d                                   aUVvPrev = theUVvArc;
  for (int k = 1; k <= NbP; k++)
  {
    // tangent to the line, the steps across closing in on it as the
    // square, those along the spine even
    const double s = double(k - 1) / (NbP - 1);
    const double t = theToLine ? 1. - (1. - s) * (1. - s) : s;
    double       w = aW0;
    if (k == NbP)
    {
      w = aW1;
    }
    else if (k > 1)
    {
      // the root between the end, inside, and the line's touching point,
      // outside but at the line
      double a = aW0, b = aW1;
      if (aDist(b, t, aBid) * aSideIn > 0. || aDist(a, t, aBid) * aSideIn < 0.)
      {
        return false;
      }
      for (int it = 0; it < 80; it++)
      {
        const double m = 0.5 * (a + b);
        if (aDist(m, t, aBid) * aSideIn > 0.)
        {
          a = m;
        }
        else
        {
          b = m;
        }
      }
      w = 0.5 * (a + b);
    }
    const gp_Pnt2d aP = aUV(w, t);
    aUVs->SetValue(k, aP);
    aPnts->SetValue(k, theS->Value(aP.X(), aP.Y()));
    gp_Pnt2d aUVv;
    SignedDistance(theSv, aPnts->Value(k), aUVv);
    // on the period of the previous point, of the arc's to start
    if (theSv.IsUPeriodic())
    {
      aUVv.SetX(ElCLib::InPeriod(aUVv.X(),
                                 aUVvPrev.X() - 0.5 * theSv.UPeriod(),
                                 aUVvPrev.X() + 0.5 * theSv.UPeriod()));
    }
    if (theSv.IsVPeriodic())
    {
      aUVv.SetY(ElCLib::InPeriod(aUVv.Y(),
                                 aUVvPrev.Y() - 0.5 * theSv.VPeriod(),
                                 aUVvPrev.Y() + 0.5 * theSv.VPeriod()));
    }
    aUVvs->SetValue(k, aUVv);
    aUVvPrev = aUVv;
    aPars->SetValue(k,
                    k == 1 ? 0.
                           : aPars->Value(k - 1) + aPnts->Value(k).Distance(aPnts->Value(k - 1)));
  }
  if (aPars->Value(NbP) <= Precision::Confusion())
  {
    return false;
  }
  GeomAPI_Interpolate   anI3(aPnts, aPars, false, Precision::Confusion());
  Geom2dAPI_Interpolate anIS(aUVs, aPars, false, Precision::PConfusion());
  Geom2dAPI_Interpolate anIV(aUVvs, aPars, false, Precision::PConfusion());
  anI3.Perform();
  anIS.Perform();
  anIV.Perform();
  if (!anI3.IsDone() || !anIS.IsDone() || !anIV.IsDone())
  {
    return false;
  }
  theEnd.Cut = anI3.Curve();
  theEnd.PcS = anIS.Curve();
  theEnd.PcV = anIV.Curve();
  const occ::handle<GeomAdaptor_Surface> aHS  = new GeomAdaptor_Surface(theS);
  const occ::handle<GeomAdaptor_Surface> aHSv = new GeomAdaptor_Surface(theSv);
  theEnd.Tol = std::max(ChFi3d_EvalTolReached(aHS, theEnd.PcS, aHSv, theEnd.PcV, theEnd.Cut),
                        Precision::Confusion());
  return true;
}

//! A curve through <theP> interpolated by chord length, with its pcurves
//! through <theUV1> and <theUV2>; false if any interpolation fails.
bool InterpolateOn(const NCollection_Array1<gp_Pnt>&   theP,
                   const NCollection_Array1<gp_Pnt2d>& theUV1,
                   const NCollection_Array1<gp_Pnt2d>& theUV2,
                   occ::handle<Geom_Curve>&            theC,
                   occ::handle<Geom2d_Curve>&          thePc1,
                   occ::handle<Geom2d_Curve>&          thePc2)
{
  const int                                  aNb   = theP.Length();
  occ::handle<NCollection_HArray1<gp_Pnt>>   aPnts = new NCollection_HArray1<gp_Pnt>(1, aNb);
  occ::handle<NCollection_HArray1<gp_Pnt2d>> aUV1s = new NCollection_HArray1<gp_Pnt2d>(1, aNb);
  occ::handle<NCollection_HArray1<gp_Pnt2d>> aUV2s = new NCollection_HArray1<gp_Pnt2d>(1, aNb);
  occ::handle<NCollection_HArray1<double>>   aPars = new NCollection_HArray1<double>(1, aNb);
  for (int k = 1; k <= aNb; k++)
  {
    aPnts->SetValue(k, theP(theP.Lower() + k - 1));
    aUV1s->SetValue(k, theUV1(theUV1.Lower() + k - 1));
    aUV2s->SetValue(k, theUV2(theUV2.Lower() + k - 1));
    aPars->SetValue(k,
                    k == 1 ? 0.
                           : aPars->Value(k - 1) + aPnts->Value(k).Distance(aPnts->Value(k - 1)));
  }
  if (aPars->Value(aNb) <= Precision::Confusion())
  {
    return false;
  }
  GeomAPI_Interpolate   anI3(aPnts, aPars, false, Precision::Confusion());
  Geom2dAPI_Interpolate anI1(aUV1s, aPars, false, Precision::PConfusion());
  Geom2dAPI_Interpolate anI2(aUV2s, aPars, false, Precision::PConfusion());
  anI3.Perform();
  anI1.Perform();
  anI2.Perform();
  if (!anI3.IsDone() || !anI1.IsDone() || !anI2.IsDone())
  {
    return false;
  }
  theC   = anI3.Curve();
  thePc1 = anI1.Curve();
  thePc2 = anI2.Curve();
  return true;
}

//! <theC> moved by <theDU> in u; <theC> itself when that is none.
occ::handle<Geom2d_Curve> Shifted(const occ::handle<Geom2d_Curve>& theC, const double theDU)
{
  return theDU == 0. ? theC : occ::down_cast<Geom2d_Curve>(theC->Translated(gp_Vec2d(theDU, 0.)));
}

//! <theUV> moved by periods of <theS> to within half a period of <thePrev>.
void InPeriodOf(const GeomAdaptor_Surface& theS, const gp_Pnt2d& thePrev, gp_Pnt2d& theUV)
{
  if (theS.IsUPeriodic())
  {
    theUV.SetX(ElCLib::InPeriod(theUV.X(),
                                thePrev.X() - 0.5 * theS.UPeriod(),
                                thePrev.X() + 0.5 * theS.UPeriod()));
  }
  if (theS.IsVPeriodic())
  {
    theUV.SetY(ElCLib::InPeriod(theUV.Y(),
                                thePrev.Y() - 0.5 * theS.VPeriod(),
                                thePrev.Y() + 0.5 * theS.VPeriod()));
  }
}

//! The side's seam through <theVtx> (a closed side: a round plate's rim
//! with its seam there), null where there is none. False where there is
//! one that the side carried up cannot take: more than one, a side other
//! than a cylinder or a cone, a seam not a line or closed at the vertex.
bool SideSeamAt(const TopoDS_Face&         theFv,
                const TopoDS_Vertex&       theVtx,
                const GeomAdaptor_Surface& theSv,
                TopoDS_Edge&               theSeam)
{
  theSeam.Nullify();
  for (TopExp_Explorer ex(theFv, TopAbs_EDGE); ex.More(); ex.Next())
  {
    if (HasVertex(TopoDS::Edge(ex.Current()), theVtx)
        && BRep_Tool::IsClosed(TopoDS::Edge(ex.Current()), theFv))
    {
      if (!theSeam.IsNull() && !theSeam.IsSame(ex.Current()))
      {
        return false;
      }
      theSeam = TopoDS::Edge(ex.Current());
    }
  }
  if (theSeam.IsNull())
  {
    return true;
  }
  if ((theSv.GetType() != GeomAbs_Cylinder && theSv.GetType() != GeomAbs_Cone)
      || !theSv.IsUPeriodic())
  {
    return false;
  }
  TopoDS_Vertex aV0, aV1;
  TopExp::Vertices(theSeam, aV0, aV1);
  double                  aF = 0., aL = 0.;
  occ::handle<Geom_Curve> aCS = BRep_Tool::Curve(theSeam, aF, aL);
  while (!occ::down_cast<Geom_TrimmedCurve>(aCS).IsNull())
  {
    aCS = occ::down_cast<Geom_TrimmedCurve>(aCS)->BasisCurve();
  }
  return !occ::down_cast<Geom_Line>(aCS).IsNull() && !aV0.IsSame(aV1);
}

//! Whether <theP> lies on the line of the side's seam <theSeam> past
//! <theVtx>, its end, by <theTol>: the seam carried up to it.
bool SeamRunsTo(const TopoDS_Edge&   theSeam,
                const TopoDS_Vertex& theVtx,
                const gp_Pnt&        theP,
                const double         theTol)
{
  TopoDS_Vertex aV0, aV1;
  TopExp::Vertices(theSeam, aV0, aV1);
  double                  aF = 0., aL = 0.;
  occ::handle<Geom_Curve> aCS = BRep_Tool::Curve(theSeam, aF, aL);
  while (!occ::down_cast<Geom_TrimmedCurve>(aCS).IsNull())
  {
    aCS = occ::down_cast<Geom_TrimmedCurve>(aCS)->BasisCurve();
  }
  const occ::handle<Geom_Line> aLS = occ::down_cast<Geom_Line>(aCS);
  const gp_Pnt                 aPV = BRep_Tool::Pnt(theVtx);
  const gp_Pnt                 aPO = BRep_Tool::Pnt(aV0.IsSame(theVtx) ? aV1 : aV0);
  return !aLS.IsNull() && aLS->Lin().Distance(theP) <= theTol
         && gp_Vec(aPV, theP).Dot(gp_Vec(aPO, aPV)) > 0.;
}

//! The side carried up over the gap of a pinch, a face of its own on a copy
//! of the side's surface <theSv>, lying in the solid as the side <theFv>
//! does: its index in the DS.
int CarriedUpSide(
  TopOpeBRepDS_DataStructure&                                                        DStr,
  const occ::handle<Geom_Surface>&                                                   theSv,
  const TopoDS_Face&                                                                 theFv,
  const int                                                                          theSolid,
  const TopoDS_Vertex&                                                               theVtx,
  const double                                                                       theTol,
  NCollection_DataMap<TopoDS_Shape, NCollection_List<int>, TopTools_ShapeMapHasher>& theEVIMap)
{
  const int IsurfV =
    DStr.AddSurface(TopOpeBRepDS_Surface(occ::down_cast<Geom_Surface>(theSv->Copy()), theTol));
  DStr.ChangeShapeInterferences(theSolid).Append(
    new TopOpeBRepDS_SolidSurfaceInterference(TopOpeBRepDS_Transition(theFv.Orientation()),
                                              TopOpeBRepDS_SOLID,
                                              theSolid,
                                              TopOpeBRepDS_SURFACE,
                                              IsurfV));
  if (!theEVIMap.IsBound(theVtx))
  {
    theEVIMap.Bind(theVtx, NCollection_List<int>());
  }
  theEVIMap.ChangeFind(theVtx).Append(IsurfV);
  return IsurfV;
}

//! The arc of <theCPArc> from its point (<theIndArc> in the DS) to <theVtx>
//! (<theIVtx>), in the side <theFv> (FORWARD; <theIShape> in the DS) as the
//! arc is, and on the side carried up (<theIsurfV>) the other way, its
//! pcurve there moved by <theShift> in u.
void AddArcPiece(TopOpeBRepDS_DataStructure& DStr,
                 const TopoDS_Face&          theFv,
                 const int                   theIShape,
                 const int                   theIsurfV,
                 const TopoDS_Vertex&        theVtx,
                 const int                   theIVtx,
                 const ChFiDS_CommonPoint&   theCPArc,
                 const int                   theIndArc,
                 const double                theShift)
{
  const TopoDS_Edge&              anArc = theCPArc.Arc();
  double                          aF = 0., aL = 0., aPF = 0., aPL = 0.;
  const occ::handle<Geom_Curve>   aCA  = BRep_Tool::Curve(anArc, aF, aL);
  const occ::handle<Geom2d_Curve> aPcA = BRep_Tool::CurveOnSurface(anArc, theFv, aPF, aPL);
  TopoDS_Vertex                   aVA0, aVA1;
  TopExp::Vertices(anArc, aVA0, aVA1);
  const double aT = theCPArc.ParameterOnArc();
  const double aTV =
    aVA0.IsSame(aVA1) ? (aT - aF < aL - aT ? aF : aL) : (aVA0.IsSame(theVtx) ? aF : aL);
  const bool isVLow = aTV < aT;
  const int  IC     = DStr.AddCurve(
    TopOpeBRepDS_Curve(new Geom_TrimmedCurve(aCA, std::min(aT, aTV), std::max(aT, aTV)),
                       BRep_Tool::Tolerance(anArc)));
  DStr.ChangeCurveInterferences(IC).Append(ChFi3d_FilPointInDS(TopAbs_FORWARD,
                                                               IC,
                                                               isVLow ? theIVtx : theIndArc,
                                                               std::min(aT, aTV),
                                                               isVLow || theCPArc.IsVertex()));
  DStr.ChangeCurveInterferences(IC).Append(ChFi3d_FilPointInDS(TopAbs_REVERSED,
                                                               IC,
                                                               isVLow ? theIndArc : theIVtx,
                                                               std::max(aT, aTV),
                                                               !isVLow || theCPArc.IsVertex()));
  TopAbs_Orientation anOrA = TopAbs_FORWARD;
  for (TopExp_Explorer ex(theFv, TopAbs_EDGE); ex.More(); ex.Next())
  {
    if (ex.Current().IsSame(anArc))
    {
      anOrA = ex.Current().Orientation();
      break;
    }
  }
  DStr.ChangeShapeInterferences(theIShape).Append(ChFi3d_FilCurveInDS(IC, theIShape, aPcA, anOrA));
  DStr.ChangeSurfaceInterferences(theIsurfV).Append(
    ChFi3d_FilCurveInDS(IC, theIsurfV, Shifted(aPcA, theShift), TopAbs::Reverse(anOrA)));
}

//=======================================================================
// function : PinchGap
// purpose  : The pinch of PinchEnds where the boss's line never touches the
//           side: a drafted boss (a cone narrowing upward) flush with the
//           side at its foot only. The fillet made on a plate reaching past
//           the side and cut by it is then one ring, its line on the boss
//           closed, its lower edge notched where the side cuts it below that
//           line: the side's cut runs across the gap between the two ends,
//           from one end's point on the plate's edge to the other's.
//           The stripe ends where it does, each end its section across the
//           fillet's surface; the gap between them is closed by a face of
//           that surface run on (a periodic one, the torus round a cone on a
//           plane), bounded by the two sections, the boss's line run on
//           across, and the side's cut. False, nothing changed, when the
//           ends are not such.
//=======================================================================

bool PinchGap(const TopoDS_Vertex&                Vtx,
              const occ::handle<ChFiDS_Stripe>&   aSt,
              PinchEnd                            anEnds[2],
              TopOpeBRepDS_DataStructure&         DStr,
              const TopoDS_Face&                  Fv,
              const occ::handle<Geom_Surface>&    aSv,
              const GeomAdaptor_Surface&          aGSv,
              const bool                          isTrimmed,
              const double                        theTol3d,
              NCollection_DataMap<TopoDS_Shape, NCollection_List<int>, TopTools_ShapeMapHasher>&
                theEVIMap)
{
  // Each end a section of its surface at one u, its lines along u.
  struct Section
  {
    occ::handle<Geom_Surface> S;
    double                    U = 0., VArc = 0., VAir = 0.;
  } aSec[2];
  for (int e = 0; e < 2; e++)
  {
    const PinchEnd&                anEnd = anEnds[e];
    const ChFiDS_FaceInterference& aFiA  = anEnd.SD->Interference(anEnd.OnsArc);
    const ChFiDS_FaceInterference& aFiB  = anEnd.SD->Interference(anEnd.OnsAir);
    if (aFiA.PCurveOnSurf().IsNull() || aFiB.PCurveOnSurf().IsNull())
    {
      return false;
    }
    aSec[e].S = DStr.Surface(anEnd.SD->Surf()).Surface();
    if (!aSec[e].S->IsUPeriodic())
    {
      return false;
    }
    gp_Pnt2d aUVA, aUVB;
    gp_Vec2d aTA, aTB;
    aFiA.PCurveOnSurf()->D1(aFiA.Parameter(anEnd.IsFirst), aUVA, aTA);
    aFiB.PCurveOnSurf()->D1(aFiB.Parameter(anEnd.IsFirst), aUVB, aTB);
    if (aTA.Magnitude() <= gp::Resolution() || aTB.Magnitude() <= gp::Resolution()
        || std::abs(aTA.Y()) > 1.e-9 * aTA.Magnitude()
        || std::abs(aTB.Y()) > 1.e-9 * aTB.Magnitude() || aTA.X() * aTB.X() <= 0.
        || std::abs(aUVA.X() - aUVB.X()) > Precision::PConfusion())
    {
      return false;
    }
    aSec[e].U    = aUVB.X();
    aSec[e].VArc = aUVA.Y();
    aSec[e].VAir = aUVB.Y();
  }
  // The gap on the first end's surface: from the last end, along the spine,
  // to the first. The last end's section there, at the same v.
  const occ::handle<Geom_Surface>& aS   = aSec[0].S;
  const ChFiDS_FaceInterference&   aFiB = anEnds[0].SD->Interference(anEnds[0].OnsAir);
  gp_Pnt2d                         aBid;
  gp_Vec2d                         aT0;
  aFiB.PCurveOnSurf()->D1(aFiB.Parameter(true), aBid, aT0);
  const double aSense = aT0.X() > 0. ? 1. : -1.;
  const double U1 = aSec[0].U, VA = aSec[0].VArc, VB = aSec[0].VAir;
  const gp_Pnt P2 = aSec[1].S->Value(aSec[1].U, aSec[1].VArc);
  const gp_Pnt Q2 = aSec[1].S->Value(aSec[1].U, aSec[1].VAir);
  GeomAPI_ProjectPointOnSurf aProj(Q2, aS);
  if (!aProj.IsDone() || aProj.NbPoints() == 0)
  {
    return false;
  }
  double aUq = 0., aVq = 0.;
  aProj.LowerDistanceParameters(aUq, aVq);
  const double aPer  = aS->UPeriod();
  const double aSpan = ElCLib::InPeriod(aSense * (U1 - aUq), 0., aPer);
  const double U2    = U1 - aSense * aSpan;
  if (aSpan <= Precision::PConfusion() || aSpan >= 0.5 * aPer
      || aS->Value(U2, VB).Distance(Q2) > theTol3d || aS->Value(U2, VA).Distance(P2) > theTol3d
      || std::abs(aSec[1].VArc - VA) > Precision::PConfusion()
      || std::abs(aSec[1].VAir - VB) > Precision::PConfusion())
  {
    return false;
  }

  // The side's cut across the gap, and the boss's line run on: at each u a
  // root across from the line on the plate (outside the side past the end)
  // to the boss's (inside it all the way: the line never reaches the side)
  const TopoDS_Face aFAir = TopoDS::Face(DStr.Shape(anEnds[0].SD->Index(anEnds[0].OnsAir)));
  const GeomAdaptor_Surface aGSair(BRep_Tool::Surface(aFAir));
  auto aDist = [&](const double theU, const double theV, gp_Pnt2d& theUVv) {
    return SignedDistance(aGSv, aS->Value(theU, theV), theUVv);
  };
  const double aIn = aDist(0.5 * (U1 + U2), VB, aBid);
  if (std::abs(aIn) <= 10. * theTol3d || Precision::IsInfinite(aIn))
  {
    return false;
  }
  const int                  NbP = 81;
  NCollection_Array1<gp_Pnt>   aCutP(1, NbP), aAirP(1, NbP);
  NCollection_Array1<gp_Pnt2d> aCutS(1, NbP), aCutV(1, NbP), aAirS(1, NbP), aAirF(1, NbP);
  for (int k = 1; k <= NbP; k++)
  {
    const double u = U2 + (U1 - U2) * (k - 1) / (NbP - 1);
    if (aDist(u, VB, aBid) * aIn <= theTol3d * std::abs(aIn))
    {
      return false;
    }
    double v = VA;
    if (k > 1 && k < NbP)
    {
      double a = VA, b = VB;
      if (aDist(u, a, aBid) * aIn >= 0.)
      {
        return false;
      }
      for (int it = 0; it < 80; it++)
      {
        const double m = 0.5 * (a + b);
        if (aDist(u, m, aBid) * aIn < 0.)
        {
          a = m;
        }
        else
        {
          b = m;
        }
      }
      v = 0.5 * (a + b);
    }
    else if (std::abs(aDist(u, VA, aBid)) > 10. * theTol3d)
    {
      return false;
    }
    aCutS(k) = gp_Pnt2d(u, v);
    aCutP(k) = aS->Value(u, v);
    gp_Pnt2d aUVv;
    SignedDistance(aGSv, aCutP(k), aUVv);
    if (k > 1)
    {
      InPeriodOf(aGSv, aCutV(k - 1), aUVv);
    }
    aCutV(k) = aUVv;
    aAirS(k) = gp_Pnt2d(u, VB);
    aAirP(k) = aS->Value(u, VB);
    gp_Pnt2d aUVf;
    SignedDistance(aGSair, aAirP(k), aUVf);
    if (k > 1)
    {
      InPeriodOf(aGSair, aAirF(k - 1), aUVf);
    }
    aAirF(k) = aUVf;
  }
  occ::handle<Geom_Curve>   aCut, aAir;
  occ::handle<Geom2d_Curve> aCutPcS, aCutPcV, aAirPcS, aAirPcF;
  if (!InterpolateOn(aCutP, aCutS, aCutV, aCut, aCutPcS, aCutPcV)
      || !InterpolateOn(aAirP, aAirS, aAirF, aAir, aAirPcS, aAirPcF))
  {
    return false;
  }
  const occ::handle<GeomAdaptor_Surface> aHS    = new GeomAdaptor_Surface(aS);
  const occ::handle<GeomAdaptor_Surface> aHSv   = new GeomAdaptor_Surface(aSv);
  const occ::handle<GeomAdaptor_Surface> aHSair = new GeomAdaptor_Surface(aGSair.Surface());
  const double aTolCut = std::max(ChFi3d_EvalTolReached(aHS, aCutPcS, aHSv, aCutPcV, aCut),
                                  Precision::Confusion());
  const double aTolAir = std::max(ChFi3d_EvalTolReached(aHS, aAirPcS, aHSair, aAirPcF, aAir),
                                  Precision::Confusion());
  if (aTolCut > 10. * theTol3d || aTolAir > 10. * theTol3d)
  {
    return false;
  }
  // the two sections: each the iso of its surface, its parameter v
  occ::handle<Geom_Curve> aSecC[2];
  for (int e = 0; e < 2; e++)
  {
    const double aLo = std::min(aSec[e].VArc, aSec[e].VAir);
    const double aHi = std::max(aSec[e].VArc, aSec[e].VAir);
    aSecC[e]         = new Geom_TrimmedCurve(aSec[e].S->UIso(aSec[e].U), aLo, aHi);
    if (aSecC[e]->Value(aLo).Distance(aSec[e].S->Value(aSec[e].U, aLo)) > Precision::Confusion()
        || aSecC[e]->Value(aHi).Distance(aSec[e].S->Value(aSec[e].U, aHi))
             > Precision::Confusion())
    {
      return false;
    }
  }
  // The side closed, its seam through the vertex (as for PinchEnds): the
  // cut crosses the seam's line, and is cut in two there; the side carried
  // up over the gap a face of its own on a copy of its surface, between the
  // cut and the side's edge either side of the vertex, made one with the
  // side after, the seam carried up to the crossing (ChFi3d_MergeSideSeams).
  // Each arc's pcurve moved by <aShiftArc> onto the cut's period at its end;
  // the cut on that face the way of <aCutOrV>.
  TopoDS_Edge aSideSeam;
  if (!SideSeamAt(Fv, Vtx, aGSv, aSideSeam))
  {
    return false;
  }
  double             aCutW = 0., aShiftArc[2] = {0., 0.};
  gp_Pnt             aCutX;
  TopAbs_Orientation aCutOrV = TopAbs_FORWARD;
  if (!aSideSeam.IsNull())
  {
    const double                    aPerV = aGSv.UPeriod();
    double                          aF = 0., aL = 0.;
    const occ::handle<Geom2d_Curve> aPcSeam = BRep_Tool::CurveOnSurface(aSideSeam, Fv, aF, aL);
    if (aPcSeam.IsNull())
    {
      return false;
    }
    const gp_Pnt2d aUVV = aPcSeam->Value(BRep_Tool::Parameter(Vtx, aSideSeam));
    const double   aW0 = aCut->FirstParameter(), aW1 = aCut->LastParameter();
    const double   aU0 = aCutPcV->Value(aW0).X(), aU1 = aCutPcV->Value(aW1).X();
    const double   aN0 = std::floor((aU0 - aUVV.X()) / aPerV);
    const double   aN1 = std::floor((aU1 - aUVV.X()) / aPerV);
    if (std::abs(aN0 - aN1) != 1.)
    {
      return false;
    }
    const double aUs = aUVV.X() + std::max(aN0, aN1) * aPerV;
    double       a = aW0, b = aW1;
    for (int it = 0; it < 80; it++)
    {
      const double m = 0.5 * (a + b);
      if ((aCutPcV->Value(m).X() - aUs) * (aU0 - aUs) > 0.)
      {
        a = m;
      }
      else
      {
        b = m;
      }
    }
    aCutW = 0.5 * (a + b);
    aCutX = aCut->Value(aCutW);
    if (!SeamRunsTo(aSideSeam, Vtx, aCutX, 10. * theTol3d))
    {
      return false;
    }
    // the cut runs from the last end's point (its first) to the first's
    for (int e = 0; e < 2; e++)
    {
      const ChFiDS_CommonPoint& aCP = anEnds[e].SD->Vertex(anEnds[e].IsFirst, anEnds[e].OnsArc);
      if (!aCP.IsOnArc() || aCP.Arc().IsSame(aSideSeam) || BRep_Tool::IsClosed(aCP.Arc(), Fv))
      {
        return false;
      }
      const occ::handle<Geom2d_Curve> aPcA = BRep_Tool::CurveOnSurface(aCP.Arc(), Fv, aF, aL);
      if (aPcA.IsNull())
      {
        return false;
      }
      const double aUA = aPcA->Value(aCP.ParameterOnArc()).X();
      const double aUC = (e == 1 ? aU0 : aU1);
      aShiftArc[e]     = aPerV * std::round((aUC - aUA) / aPerV);
    }
    // the face over the cut, between it and the vertex: the loop from the
    // cut's first point over it to its last and back by the vertex, turning
    // left or right in the side's (u, v)
    const int NbA   = 40;
    double    aArea = 0.;
    gp_Pnt2d  aPrev = aCutPcV->Value(aW0);
    for (int k = 1; k <= NbA + 1; k++)
    {
      const gp_Pnt2d aNext =
        k <= NbA ? aCutPcV->Value(aW0 + (aW1 - aW0) * k / NbA) : gp_Pnt2d(aUs, aUVV.Y());
      aArea += aPrev.X() * aNext.Y() - aNext.X() * aPrev.Y();
      aPrev = aNext;
    }
    const gp_Pnt2d aFirst = aCutPcV->Value(aW0);
    aArea += aPrev.X() * aFirst.Y() - aFirst.X() * aPrev.Y();
    aCutOrV = aArea > 0. ? TopAbs_FORWARD : TopAbs_REVERSED;
  }

  // The boss's own seam through the vertex (a boss turned so): the line run
  // on crosses it, and is cut in two there, the seam split at the crossing
  TopoDS_Edge aSeam;
  double      aSeamW = 0., aSeamPar = 0., aSeamParV = 0., aSeamDist = 0.;
  double      aUMin = 0., aUMax = 0., aPerAir = 0.;
  for (TopExp_Explorer ex(aFAir, TopAbs_EDGE); ex.More() && aSeam.IsNull(); ex.Next())
  {
    const TopoDS_Edge& aE = TopoDS::Edge(ex.Current());
    if (!HasVertex(aE, Vtx) || !BRep_Tool::IsClosed(aE, aFAir))
    {
      continue;
    }
    if (!aGSair.IsUPeriodic())
    {
      return false;
    }
    aPerAir = aGSair.UPeriod();
    double aVMin = 0., aVMax = 0.;
    BRepTools::UVBounds(aFAir, aUMin, aUMax, aVMin, aVMax);
    // the seam's u, where the line's u crosses it (a whole number of periods off)
    const double aU0 = aAirPcF->Value(aAir->FirstParameter()).X();
    const double aU1 = aAirPcF->Value(aAir->LastParameter()).X();
    const double aN0 = std::floor((aU0 - aUMin) / aPerAir);
    const double aN1 = std::floor((aU1 - aUMin) / aPerAir);
    if (aN0 == aN1)
    {
      break; // the line does not cross it
    }
    const double aUs = aUMin + std::max(aN0, aN1) * aPerAir;
    double       a = aAir->FirstParameter(), b = aAir->LastParameter();
    for (int it = 0; it < 80; it++)
    {
      const double m = 0.5 * (a + b);
      if ((aAirPcF->Value(m).X() - aUs) * (aU0 - aUs) > 0.)
      {
        a = m;
      }
      else
      {
        b = m;
      }
    }
    aSeamW = 0.5 * (a + b);
    double                  aF = 0., aL = 0.;
    occ::handle<Geom_Curve> aCE = BRep_Tool::Curve(aE, aF, aL);
    if (aCE.IsNull())
    {
      return false;
    }
    GeomAdaptor_Curve aGC(aCE, aF, aL);
    Extrema_ExtPC     anExt(aAir->Value(aSeamW), aGC);
    int               iMin = 0;
    for (int i = 1; anExt.IsDone() && i <= anExt.NbExt(); i++)
    {
      if (iMin == 0 || anExt.SquareDistance(i) < anExt.SquareDistance(iMin))
      {
        iMin = i;
      }
    }
    if (iMin == 0 || anExt.SquareDistance(iMin) > Square(10. * theTol3d))
    {
      return false;
    }
    aSeam     = aE;
    aSeamPar  = anExt.Point(iMin).Parameter();
    aSeamParV = BRep_Tool::Parameter(Vtx, aE);
    aSeamDist = std::sqrt(anExt.SquareDistance(iMin));
  }

  // Into the DS. The side's surface untrimmed: the cut lies past the face.
  if (isTrimmed)
  {
    DStr.SetNewSurface(anEnds[0].Fv, aSv);
  }
  const int IShape = DStr.AddShape(anEnds[0].Fv);
  int       IndArc[2], IndAir[2];
  for (int e = 0; e < 2; e++)
  {
    IndArc[e] = ChFi3d_IndexPointInDS(anEnds[e].SD->Vertex(anEnds[e].IsFirst, anEnds[e].OnsArc),
                                      DStr);
    IndAir[e] = ChFi3d_IndexPointInDS(anEnds[e].SD->Vertex(anEnds[e].IsFirst, anEnds[e].OnsAir),
                                      DStr);
    TopOpeBRepDS_Point& aPA = DStr.ChangePoint(IndArc[e]);
    aPA.Tolerance(std::max(aPA.Tolerance(), aTolCut));
    TopOpeBRepDS_Point& aPB = DStr.ChangePoint(IndAir[e]);
    aPB.Tolerance(std::max(aPB.Tolerance(), aTolAir));
  }
  const int IndSeam =
    aSeam.IsNull()
      ? 0
      : DStr.AddPoint(TopOpeBRepDS_Point(aAir->Value(aSeamW),
                                         std::max({aTolAir, aSeamDist, Precision::Confusion()})));
  // the stripe's ends, its sections, from its side 1 to its side 2
  TopAbs_Orientation aSecOr[2];
  int                ISec[2];
  for (int e = 0; e < 2; e++)
  {
    const PinchEnd& anEnd = anEnds[e];
    const double    aLo   = std::min(aSec[e].VArc, aSec[e].VAir);
    const double    aHi   = std::max(aSec[e].VArc, aSec[e].VAir);
    const double    aV1   = anEnd.OnsArc == 1 ? aSec[e].VArc : aSec[e].VAir;
    aSecOr[e]             = aV1 == aLo ? TopAbs_FORWARD : TopAbs_REVERSED;
    ISec[e]               = DStr.AddCurve(TopOpeBRepDS_Curve(aSecC[e], Precision::Confusion()));
    anEnd.St->SetCurve(ISec[e], anEnd.IsFirst);
    anEnd.St->SetParameters(anEnd.IsFirst, aLo, aHi);
    anEnd.St->SetOrientation(aSecOr[e], anEnd.IsFirst);
    anEnd.St->ChangePCurve(anEnd.IsFirst) =
      new Geom2d_Line(gp_Pnt2d(aSec[e].U, 0.), gp_Dir2d(0., 1.));
    anEnd.St->SetIndexPoint(anEnd.OnsArc == 1 ? IndArc[e] : IndAir[e], anEnd.IsFirst, 1);
    anEnd.St->SetIndexPoint(anEnd.OnsArc == 1 ? IndAir[e] : IndArc[e], anEnd.IsFirst, 2);
  }

  // The face over the gap, the first end's surface run on (ChFi3d_FilDS's
  // orientations: a surface data from the last end's section to the
  // first's, its lines those of the first end's, the cut in the place of
  // its line on the plate)
  const occ::handle<ChFiDS_SurfData>& aFd = anEnds[0].SD;
  const int                           Ia  = anEnds[0].OnsArc;
  const ChFiDS_FaceInterference&      aFi1 = aFd->InterferenceOnS1();
  TopAbs_Orientation trafil1 = TopAbs::Compose(DStr.Shape(aFd->IndexOfS1()).Orientation(),
                                               aFd->Orientation());
  trafil1                    = TopAbs::Compose(TopAbs::Reverse(aFi1.Transition()), trafil1);
  const TopAbs_Orientation trafil2  = TopAbs::Reverse(trafil1);
  const TopAbs_Orientation trafilA  = Ia == 1 ? trafil1 : trafil2;
  const TopAbs_Orientation trafilB  = Ia == 1 ? trafil2 : trafil1;
  const TopAbs_Orientation aOrPlate = DStr.Shape(aFd->Index(Ia)).Orientation();
  const TopAbs_Orientation aOrSide  = DStr.Shape(IShape).Orientation();
  const TopAbs_Orientation aTrSide  = TopAbs::Compose(
    TopAbs::Compose(aFd->Interference(Ia).Transition(), aOrPlate), aOrSide);
  const TopAbs_Orientation aTrAir = aFd->Interference(anEnds[0].OnsAir).Transition();

  const int Isurf = DStr.AddSurface(
    TopOpeBRepDS_Surface(aS, std::max({aTolCut, aTolAir, Precision::Confusion()})));
  DStr.ChangeShapeInterferences(aSt->SolidIndex())
    .Append(new TopOpeBRepDS_SolidSurfaceInterference(TopOpeBRepDS_Transition(aFd->Orientation()),
                                                      TopOpeBRepDS_SOLID,
                                                      aSt->SolidIndex(),
                                                      TopOpeBRepDS_SURFACE,
                                                      Isurf));
  if (!theEVIMap.IsBound(Vtx))
  {
    theEVIMap.Bind(Vtx, NCollection_List<int>());
  }
  theEVIMap.ChangeFind(Vtx).Append(Isurf);
  // its first section the stripe's last, its last the stripe's first
  DStr.ChangeSurfaceInterferences(Isurf).Append(
    ChFi3d_FilCurveInDS(ISec[1],
                        Isurf,
                        new Geom2d_Line(gp_Pnt2d(U2, 0.), gp_Dir2d(0., 1.)),
                        TopAbs::Compose(TopAbs::Reverse(trafil1), aSecOr[1])));
  DStr.ChangeSurfaceInterferences(Isurf).Append(
    ChFi3d_FilCurveInDS(ISec[0],
                        Isurf,
                        new Geom2d_Line(gp_Pnt2d(U1, 0.), gp_Dir2d(0., 1.)),
                        TopAbs::Compose(trafil1, aSecOr[0])));
  // the cut, on the gap's face and the side's; on the side carried up in
  // two across the seam's line, with the arc's pieces either side of the
  // vertex
  if (aSideSeam.IsNull())
  {
    const int ICut = DStr.AddCurve(TopOpeBRepDS_Curve(aCut, aTolCut));
    DStr.ChangeCurveInterferences(ICut).Append(
      ChFi3d_FilPointInDS(TopAbs_FORWARD, ICut, IndArc[1], aCut->FirstParameter()));
    DStr.ChangeCurveInterferences(ICut).Append(
      ChFi3d_FilPointInDS(TopAbs_REVERSED, ICut, IndArc[0], aCut->LastParameter()));
    DStr.ChangeSurfaceInterferences(Isurf).Append(
      ChFi3d_FilCurveInDS(ICut, Isurf, aCutPcS, trafilA));
    DStr.ChangeShapeInterferences(IShape).Append(
      ChFi3d_FilCurveInDS(ICut, IShape, aCutPcV, aTrSide));
  }
  else
  {
    const int IsurfV =
      CarriedUpSide(DStr, aSv, anEnds[0].Fv, aSt->SolidIndex(), Vtx, theTol3d, theEVIMap);
    const int IVtx = DStr.AddShape(Vtx);
    const int IndX =
      DStr.AddPoint(TopOpeBRepDS_Point(aCutX, std::max(aTolCut, Precision::Confusion())));
    const double aWc[3] = {aCut->FirstParameter(), aCutW, aCut->LastParameter()};
    const int    aIc[3] = {IndArc[1], IndX, IndArc[0]};
    for (int k = 0; k < 2; k++)
    {
      const int ICut = DStr.AddCurve(TopOpeBRepDS_Curve(aCut, aTolCut));
      DStr.ChangeCurveInterferences(ICut).Append(
        ChFi3d_FilPointInDS(TopAbs_FORWARD, ICut, aIc[k], aWc[k]));
      DStr.ChangeCurveInterferences(ICut).Append(
        ChFi3d_FilPointInDS(TopAbs_REVERSED, ICut, aIc[k + 1], aWc[k + 1]));
      DStr.ChangeSurfaceInterferences(Isurf).Append(
        ChFi3d_FilCurveInDS(ICut, Isurf, aCutPcS, trafilA));
      DStr.ChangeSurfaceInterferences(IsurfV).Append(
        ChFi3d_FilCurveInDS(ICut, IsurfV, aCutPcV, aCutOrV));
    }
    for (int e = 0; e < 2; e++)
    {
      AddArcPiece(DStr,
                  Fv,
                  IShape,
                  IsurfV,
                  Vtx,
                  IVtx,
                  anEnds[e].SD->Vertex(anEnds[e].IsFirst, anEnds[e].OnsArc),
                  IndArc[e],
                  aShiftArc[e]);
    }
    ChFi3d_SideSeams().Append(std::make_pair(BRep_Tool::Pnt(Vtx), aCutX));
  }
  // the boss's line run on, on the gap's face and the boss's; in two
  // across the boss's seam
  const int IFAir = aFd->Index(anEnds[0].OnsAir);
  const double aW[3]  = {aAir->FirstParameter(), aSeamW, aAir->LastParameter()};
  const int    aIP[3] = {IndAir[1], IndSeam, IndAir[0]};
  for (int k = 0; k < (IndSeam > 0 ? 2 : 1); k++)
  {
    const double aW0 = aW[k], aW1 = IndSeam > 0 ? aW[k + 1] : aW[2];
    const int    aI1 = IndSeam > 0 ? aIP[k + 1] : aIP[2];
    occ::handle<Geom2d_Curve> aPcF = aAirPcF;
    if (IndSeam > 0)
    {
      const double aUMid = aAirPcF->Value(0.5 * (aW0 + aW1)).X();
      if (aUMid > aUMax || aUMid < aUMin)
      {
        aPcF = occ::down_cast<Geom2d_Curve>(
          aAirPcF->Translated(gp_Vec2d(aUMid > aUMax ? -aPerAir : aPerAir, 0.)));
      }
    }
    const int IAir = DStr.AddCurve(TopOpeBRepDS_Curve(aAir, aTolAir));
    DStr.ChangeCurveInterferences(IAir).Append(
      ChFi3d_FilPointInDS(TopAbs_FORWARD, IAir, aIP[k], aW0));
    DStr.ChangeCurveInterferences(IAir).Append(
      ChFi3d_FilPointInDS(TopAbs_REVERSED, IAir, aI1, aW1));
    DStr.ChangeSurfaceInterferences(Isurf).Append(
      ChFi3d_FilCurveInDS(IAir, Isurf, aAirPcS, trafilB));
    DStr.ChangeShapeInterferences(IFAir).Append(ChFi3d_FilCurveInDS(IAir, IFAir, aPcF, aTrAir));
  }
  if (IndSeam > 0)
  {
    const int IE = DStr.AddShape(aSeam);
    DStr.ChangeShapeInterferences(IE).Append(
      ChFi3d_FilPointInDS(aSeamParV > aSeamPar ? TopAbs_REVERSED : TopAbs_FORWARD,
                          IE,
                          IndSeam,
                          aSeamPar));
  }
  return true;
}

//=======================================================================
// function : PerformPinchEnds
// purpose  : See PinchEnds above. False, nothing changed, when the vertex
//           is not such a corner.
//=======================================================================

bool PerformPinchEnds(const TopoDS_Vertex&                                Vtx,
                             const NCollection_List<occ::handle<ChFiDS_Stripe>>& theStripes,
                             TopOpeBRepDS_DataStructure&                         DStr,
                             const ChFiDS_Map&                                   theEFMap,
                             const ChFiDS_Map&                                   theVEMap,
                             const double                                        theTol3d,
                             NCollection_DataMap<TopoDS_Shape, NCollection_List<int>,
                                                 TopTools_ShapeMapHasher>&        theEVIMap)
{
  if (theStripes.Extent() != 2)
  {
    return false;
  }
  // One stripe ending twice here, or two ending here once each: the boss
  // flush with two sides at a corner of the plate, the plate's top cut in
  // two at the points its foot touches them, the foot with it
  const bool isOne = theStripes.First() == theStripes.Last();
  PinchEnd   anEnds[2];
  int        IndAir = 0;
  NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator itSt(theStripes);
  for (int e = 0; e < 2; e++, itSt.Next())
  {
    PinchEnd& anEnd = anEnds[e];
    anEnd.St                             = itSt.Value();
    const occ::handle<ChFiDS_Spine>& aSp = anEnd.St->Spine();
    if (aSp.IsNull() || anEnd.St->SetOfSurfData().IsNull()
        || anEnd.St->SetOfSurfData()->Length() == 0)
    {
      return false;
    }
    const bool isFirstHere = aSp->FirstVertex().IsSame(Vtx);
    const bool isLastHere  = aSp->LastVertex().IsSame(Vtx);
    if (isOne ? !(isFirstHere && isLastHere) : isFirstHere == isLastHere)
    {
      return false;
    }
    anEnd.IsFirst = isOne ? e == 0 : isFirstHere;
    anEnd.SD =
      anEnd.St->SetOfSurfData()->Value(anEnd.IsFirst ? 1 : anEnd.St->SetOfSurfData()->Length());
    const ChFiDS_CommonPoint& aCP1 = anEnd.SD->Vertex(anEnd.IsFirst, 1);
    const ChFiDS_CommonPoint& aCP2 = anEnd.SD->Vertex(anEnd.IsFirst, 2);
    if (aCP1.IsOnArc() == aCP2.IsOnArc() || anEnd.SD->IndexOfS1() <= 0
        || anEnd.SD->IndexOfS2() <= 0)
    {
      return false;
    }
    anEnd.OnsArc                     = aCP1.IsOnArc() ? 1 : 2;
    anEnd.OnsAir                     = 3 - anEnd.OnsArc;
    const ChFiDS_CommonPoint& aCPArc = anEnd.SD->Vertex(anEnd.IsFirst, anEnd.OnsArc);
    const TopoDS_Edge&        anArc  = aCPArc.Arc();
    if (!HasVertex(anArc, Vtx) || !theEFMap.Contains(anArc) || theEFMap(anArc).Extent() != 2)
    {
      return false;
    }
    const TopoDS_Shape& aFArc = DStr.Shape(anEnd.SD->Index(anEnd.OnsArc));
    for (NCollection_List<TopoDS_Shape>::Iterator itF(theEFMap(anArc)); itF.More(); itF.Next())
    {
      if (!itF.Value().IsSame(aFArc))
      {
        anEnd.Fv = TopoDS::Face(itF.Value());
      }
    }
    const int anIndAir = anEnd.SD->Index(anEnd.OnsAir);
    if (anEnd.Fv.IsNull() || (e == 1 && (!anEnd.Fv.IsSame(anEnds[0].Fv) || anIndAir != IndAir)))
    {
      return false;
    }
    IndAir                          = anIndAir;
    anEnd.WEnd                      = anEnd.SD->Interference(anEnd.OnsAir).Parameter(anEnd.IsFirst);
  }
  // the boss touches the side at the vertex, no edge between them
  const TopoDS_Shape& aFop = DStr.Shape(IndAir);
  if (aFop.IsSame(anEnds[0].Fv))
  {
    return false;
  }
  for (NCollection_List<TopoDS_Shape>::Iterator itE(theVEMap(Vtx)); itE.More(); itE.Next())
  {
    bool isV = false, isOp = false;
    for (NCollection_List<TopoDS_Shape>::Iterator itF(theEFMap(itE.Value())); itF.More();
         itF.Next())
    {
      isV  = isV || itF.Value().IsSame(anEnds[0].Fv);
      isOp = isOp || itF.Value().IsSame(aFop);
    }
    if (isV && isOp)
    {
      return false;
    }
  }
  // the side's surface, untrimmed
  TopoDS_Face Fv = anEnds[0].Fv;
  Fv.Orientation(TopAbs_FORWARD);
  occ::handle<Geom_Surface> aSv = BRep_Tool::Surface(Fv);
  const occ::handle<Geom_RectangularTrimmedSurface> aTrimmed =
    occ::down_cast<Geom_RectangularTrimmedSurface>(aSv);
  if (!aTrimmed.IsNull())
  {
    aSv = aTrimmed->BasisSurface();
  }
  const GeomAdaptor_Surface aGSv(aSv);
  switch (aGSv.GetType())
  {
    case GeomAbs_Plane:
    case GeomAbs_Cylinder:
    case GeomAbs_Cone:
    case GeomAbs_Sphere:
    case GeomAbs_Torus:
      break;
    default:
      return false;
  }
  // where each end's line touches the side; one point. A line that never
  // does (a drafted boss): the gap closed by the fillet's surface run on,
  // from the end that is a stripe's last to the one that is a stripe's
  // first (two stripes: one of each, or none made; PinchGap goes along the
  // spine from anEnds[1] to anEnds[0])
  for (int e = 0; e < 2; e++)
  {
    PinchEnd&                        anEnd = anEnds[e];
    const occ::handle<Geom_Surface>& aS    = DStr.Surface(anEnd.SD->Surf()).Surface();
    if (!LineTouching(aS,
                      anEnd.SD->Interference(anEnd.OnsAir),
                      anEnd.IsFirst,
                      aGSv,
                      theTol3d,
                      anEnd.WAir,
                      anEnd.PAir))
    {
      if (!isOne && anEnds[0].IsFirst == anEnds[1].IsFirst)
      {
        return false;
      }
      if (!anEnds[0].IsFirst)
      {
        std::swap(anEnds[0], anEnds[1]);
      }
      return PinchGap(Vtx,
                      anEnds[0].St,
                      anEnds,
                      DStr,
                      Fv,
                      aSv,
                      aGSv,
                      !aTrimmed.IsNull(),
                      theTol3d,
                      theEVIMap);
    }
  }
  // One surface data from end to end (the foot one closed edge, the boss's
  // seam through the vertex), its line on the boss run on from both ends to
  // one point: a whole turn of a periodic line, not a hair over -- the
  // points are found to some 1e-8, and BRepCheck refuses a range past the
  // period (the rim's boss turned so, "BRep_API: command not done")
  if (anEnds[0].SD == anEnds[1].SD && anEnds[0].IsFirst != anEnds[1].IsFirst)
  {
    const int                     anIL = anEnds[0].SD->Interference(anEnds[0].OnsAir).LineIndex();
    const occ::handle<Geom_Curve> aCL  = anIL > 0 ? DStr.Curve(anIL).Curve() : nullptr;
    if (!aCL.IsNull() && aCL->IsPeriodic())
    {
      PinchEnd& aF = anEnds[0].IsFirst ? anEnds[0] : anEnds[1];
      PinchEnd& aL = anEnds[0].IsFirst ? anEnds[1] : anEnds[0];
      if (aL.WAir - aF.WAir > aCL->Period())
      {
        aL.WAir = aF.WAir + aCL->Period();
      }
    }
  }
  const double aGap = anEnds[0].PAir.Distance(anEnds[1].PAir);
  if (aGap > 10 * theTol3d)
  {
    return false;
  }
  // The fillet tangent to the side there, as to the boss: the cut then
  // crosses the line at the point, the side's surface meeting the
  // fillet's in two branches. A chamfer is not: its cut is tangent to its
  // line there, a corner of no angle on its face, the cut's position along
  // the spine going as the root of its distance from the line
  bool isToLine = false;
  for (int e = 0; e < 2; e++)
  {
    const PinchEnd&                  anEnd = anEnds[e];
    const occ::handle<Geom_Surface>& aS    = DStr.Surface(anEnd.SD->Surf()).Surface();
    const gp_Pnt2d aUV = anEnd.SD->Interference(anEnd.OnsAir).PCurveOnSurf()->Value(anEnd.WAir);
    gp_Pnt         aPS, aPV;
    gp_Vec         aDU, aDV, aDUv, aDVv;
    aS->D1(aUV.X(), aUV.Y(), aPS, aDU, aDV);
    gp_Pnt2d aUVv;
    SignedDistance(aGSv, aPS, aUVv);
    aGSv.D1(aUVv.X(), aUVv.Y(), aPV, aDUv, aDVv);
    const gp_Vec aN = aDU.Crossed(aDV), aNv = aDUv.Crossed(aDVv);
    if (aN.Magnitude() <= gp::Resolution() || aNv.Magnitude() <= gp::Resolution())
    {
      return false;
    }
    const bool isTangent = aN.Crossed(aNv).Magnitude() <= 1.e-4 * aN.Magnitude() * aNv.Magnitude();
    if (e == 1 && isTangent == isToLine)
    {
      return false;
    }
    isToLine = !isTangent;
  }
  // the cuts
  for (int e = 0; e < 2; e++)
  {
    PinchEnd&                        anEnd  = anEnds[e];
    const ChFiDS_CommonPoint&        aCPArc = anEnd.SD->Vertex(anEnd.IsFirst, anEnd.OnsArc);
    double                           aF = 0., aL = 0.;
    const occ::handle<Geom2d_Curve>  aPcArc = BRep_Tool::CurveOnSurface(aCPArc.Arc(), Fv, aF, aL);
    if (aPcArc.IsNull()
        || !CutEnd(DStr.Surface(anEnd.SD->Surf()).Surface(),
                   aGSv,
                   aPcArc->Value(aCPArc.ParameterOnArc()),
                   anEnd,
                   isToLine,
                   theTol3d))
    {
      return false;
    }
    if (anEnd.Tol > 10 * theTol3d)
    {
      return false;
    }
  }

  // The side closed, its seam through the vertex (a round plate's rim with
  // its seam there): the cuts meet across the seam, and the side carried up
  // over the gap is a face of its own on a copy of its surface, between the
  // cuts and the side's edge either side of the vertex, put in that edge's
  // stead in the side; after the build the two are made one, the seam
  // carried up from the vertex to the point (ChFi3d_MergeSideSeams). A
  // cylinder or a cone, its seam a line, and the point on that line past
  // the vertex.
  TopoDS_Edge aSideSeam;
  if (!SideSeamAt(Fv, Vtx, aGSv, aSideSeam))
  {
    return false;
  }
  // the second cut's pcurve on the side's surface moved by <aShift> in u to
  // meet the first's
  double aShift = 0.;
  if (!aSideSeam.IsNull())
  {
    if (!SeamRunsTo(aSideSeam,
                    Vtx,
                    anEnds[0].PAir.Translated(0.5 * gp_Vec(anEnds[0].PAir, anEnds[1].PAir)),
                    10 * theTol3d))
    {
      return false;
    }
    for (int e = 0; e < 2; e++)
    {
      const TopoDS_Edge& anArc = anEnds[e].SD->Vertex(anEnds[e].IsFirst, anEnds[e].OnsArc).Arc();
      if (anArc.IsSame(aSideSeam) || BRep_Tool::IsClosed(anArc, Fv))
      {
        return false;
      }
    }
    const double aU0 = anEnds[0].PcV->Value(anEnds[0].PcV->LastParameter()).X();
    const double aU1 = anEnds[1].PcV->Value(anEnds[1].PcV->LastParameter()).X();
    aShift           = aGSv.UPeriod() * std::round((aU0 - aU1) / aGSv.UPeriod());
  }
  // The two cuts meet on the side's surface as in space
  {
    const gp_Pnt2d aUV0 = anEnds[0].PcV->Value(anEnds[0].PcV->LastParameter());
    const gp_Pnt2d aUV1 =
      anEnds[1].PcV->Value(anEnds[1].PcV->LastParameter()).Translated(gp_Vec2d(aShift, 0.));
    const double aRes = std::max(aGSv.UResolution(10 * theTol3d), aGSv.VResolution(10 * theTol3d));
    if (aUV0.Distance(aUV1) > aRes)
    {
      return false;
    }
  }

  // The boss's own seam through the vertex (a boss turned so): the point
  // splits it, its stretch from the vertex to the point under the fillet
  TopoDS_Edge aSeam;
  double      aSeamPar = 0., aSeamParV = 0., aSeamDist = 0.;
  {
    const gp_Pnt aPm = anEnds[0].PAir.Translated(0.5 * gp_Vec(anEnds[0].PAir, anEnds[1].PAir));
    for (NCollection_List<TopoDS_Shape>::Iterator itE(theVEMap(Vtx)); itE.More() && aSeam.IsNull();
         itE.Next())
    {
      const TopoDS_Edge& aE = TopoDS::Edge(itE.Value());
      if (!BRep_Tool::IsClosed(aE, TopoDS::Face(aFop)))
      {
        continue;
      }
      double                  aF = 0., aL = 0.;
      occ::handle<Geom_Curve> aCE = BRep_Tool::Curve(aE, aF, aL);
      if (aCE.IsNull())
      {
        return false;
      }
      GeomAdaptor_Curve aGC(aCE, aF, aL);
      Extrema_ExtPC     anExt(aPm, aGC);
      int               iMin = 0;
      for (int i = 1; anExt.IsDone() && i <= anExt.NbExt(); i++)
      {
        if (iMin == 0 || anExt.SquareDistance(i) < anExt.SquareDistance(iMin))
        {
          iMin = i;
        }
      }
      if (iMin == 0 || anExt.SquareDistance(iMin) > Square(10. * theTol3d))
      {
        return false;
      }
      aSeam     = aE;
      aSeamPar  = anExt.Point(iMin).Parameter();
      aSeamParV = BRep_Tool::Parameter(Vtx, aE);
      aSeamDist = std::sqrt(anExt.SquareDistance(iMin));
    }
  }

  // Into the DS. The side's surface untrimmed: the cuts lie past the face.
  if (!aTrimmed.IsNull())
  {
    DStr.SetNewSurface(anEnds[0].Fv, aSv);
  }
  const int IShape = DStr.AddShape(anEnds[0].Fv);
  const gp_Pnt aP  = anEnds[0].PAir.Translated(0.5 * gp_Vec(anEnds[0].PAir, anEnds[1].PAir));
  const int    IndP =
    DStr.AddPoint(TopOpeBRepDS_Point(aP,
                                     std::max({0.5 * aGap + anEnds[0].Tol,
                                               0.5 * aGap + anEnds[1].Tol,
                                               Precision::Confusion()})));
  // the side carried up, a face of its own
  const int IsurfV = aSideSeam.IsNull() ? 0
                                        : CarriedUpSide(DStr,
                                                        aSv,
                                                        anEnds[0].Fv,
                                                        anEnds[0].St->SolidIndex(),
                                                        Vtx,
                                                        theTol3d,
                                                        theEVIMap);
  const int IVtx   = aSideSeam.IsNull() ? 0 : DStr.AddShape(Vtx);
  for (int e = 0; e < 2; e++)
  {
    PinchEnd&                anEnd  = anEnds[e];
    ChFiDS_CommonPoint&      aCPAir = anEnd.SD->ChangeVertex(anEnd.IsFirst, anEnd.OnsAir);
    ChFiDS_FaceInterference& aFiAir = anEnd.SD->ChangeInterference(anEnd.OnsAir);
    const ChFiDS_CommonPoint& aCPArc = anEnd.SD->Vertex(anEnd.IsFirst, anEnd.OnsArc);
    aFiAir.SetParameter(anEnd.WAir, anEnd.IsFirst);
    aCPAir.Reset();
    aCPAir.SetPoint(aP);
    aCPAir.SetTolerance(DStr.Point(IndP).Tolerance());
    // the cut from the stripe's side 1 to its side 2
    occ::handle<Geom_Curve>   aC  = anEnd.Cut;
    occ::handle<Geom2d_Curve> aPs = anEnd.PcS, aPc = anEnd.PcV;
    double                    aU0 = aC->FirstParameter(), aU1 = aC->LastParameter();
    if (anEnd.OnsArc == 2)
    {
      aC  = aC->Reversed();
      aPs = aPs->Reversed();
      aPc = aPc->Reversed();
      const double aU = aU0;
      aU0             = aC->ReversedParameter(aU1);
      aU1             = aC->ReversedParameter(aU);
    }
    // the orientation of the cut on the side, as the corner of one fillet
    // gives it from the arc (PerformOneCorner)
    TopAbs_Orientation Et = TopAbs_FORWARD;
    for (TopExp_Explorer ex(Fv, TopAbs_EDGE); ex.More(); ex.Next())
    {
      if (ex.Current().IsSame(aCPArc.Arc()))
      {
        Et = TopAbs::Compose(ex.Current().Orientation(), aCPArc.TransitionOnArc());
        if (anEnd.OnsArc == 1)
        {
          Et = TopAbs::Reverse(Et);
        }
        break;
      }
    }
    const int ICurve = DStr.AddCurve(TopOpeBRepDS_Curve(aC, anEnd.Tol));
    if (IsurfV == 0)
    {
      DStr.ChangeShapeInterferences(IShape).Append(ChFi3d_FilCurveInDS(ICurve, IShape, aPc, Et));
    }
    else
    {
      DStr.ChangeSurfaceInterferences(IsurfV).Append(
        ChFi3d_FilCurveInDS(ICurve, IsurfV, Shifted(aPc, e == 1 ? aShift : 0.), Et));
    }
    anEnd.St->SetCurve(ICurve, anEnd.IsFirst);
    anEnd.St->SetParameters(anEnd.IsFirst, aU0, aU1);
    anEnd.St->ChangePCurve(anEnd.IsFirst) = aPs;
    const int IndArc                 = ChFi3d_IndexPointInDS(aCPArc, DStr);
    anEnd.St->SetIndexPoint(anEnd.OnsArc == 1 ? IndArc : IndP, anEnd.IsFirst, 1);
    anEnd.St->SetIndexPoint(anEnd.OnsArc == 1 ? IndP : IndArc, anEnd.IsFirst, 2);
    if (IsurfV == 0)
    {
      continue;
    }
    AddArcPiece(DStr, Fv, IShape, IsurfV, Vtx, IVtx, aCPArc, IndArc, e == 1 ? aShift : 0.);
  }
  if (IsurfV != 0)
  {
    ChFi3d_SideSeams().Append(std::make_pair(BRep_Tool::Pnt(Vtx), aP));
  }
  if (!aSeam.IsNull())
  {
    TopOpeBRepDS_Point& aDP = DStr.ChangePoint(IndP);
    aDP.Tolerance(std::max(aDP.Tolerance(), aSeamDist));
    const int IE = DStr.AddShape(aSeam);
    DStr.ChangeShapeInterferences(IE).Append(
      ChFi3d_FilPointInDS(aSeamParV > aSeamPar ? TopAbs_REVERSED : TopAbs_FORWARD,
                          IE,
                          IndP,
                          aSeamPar));
  }
  return true;
}
} // namespace

// True at a vertex where two fillets end whose faces there share nothing:
// four edges, the two filleted ones opposite, the two sharp ones between.
static bool OppositeEnds(const TopoDS_Vertex&                                theV,
                         const NCollection_List<occ::handle<ChFiDS_Stripe>>& theL)
{
  int aF[2][2];
  int k = 0;
  for (NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator It(theL); It.More(); It.Next(), k++)
  {
    if (k > 1)
    {
      return false;
    }
    int       aSens = 0;
    const int aNum  = ChFi3d_IndexOfSurfData(theV, It.Value(), aSens);
    const occ::handle<ChFiDS_SurfData>& aSD = It.Value()->SetOfSurfData()->Value(aNum);
    aF[k][0]                                = aSD->IndexOfS1();
    aF[k][1]                                = aSD->IndexOfS2();
  }
  if (k != 2 || theL.First() == theL.Last())
  {
    return false;
  }
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      if (aF[0][i] <= 0 || aF[0][i] == aF[1][j])
      {
        return false;
      }
    }
  }
  return true;
}

//=================================================================================================

void ChFi3d_Builder::PerformFilletOnVertex(const int Index)
{

  NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator It;
  occ::handle<ChFiDS_Stripe>                             stripe;
  occ::handle<ChFiDS_Spine>                              sp;
  const TopoDS_Vertex&                                   Vtx = myVDataMap.FindKey(Index);

  occ::handle<ChFiDS_SurfData> Fd;
  int                          i;
  bool                         nondegenere      = true;
  bool                         toujoursdegenere = true;
  bool                         isfirst          = false;
  // A fillet's last surface can be one with a face on one side only: what
  // is left of a line that runs on past the spine's end over the next piece
  // of its face, kept in coplanar pieces (a body without Refine), the split
  // ending at the vertex. The corner of one fillet removes that surplus
  // before it starts (PerformOneCorner, PerformIntersectionAtEnd); the
  // corners of several, and the set back one, took it for the fillet's end,
  // a surface with no face in common with the next fillet.
  const bool toRemoveSurplus = myVDataMap(Index).Extent() >= 2 || HasSetbackAt(Index);
  for (It.Initialize(myVDataMap(Index)), i = 0; It.More(); It.Next(), i++)
  {
    stripe = It.Value();
    sp     = stripe->Spine();
    // SurfData and its CommonPoints,
    int sens                      = 0;
    int num                       = ChFi3d_IndexOfSurfData(Vtx, stripe, sens);
    isfirst                       = (sens == 1);
    if (toRemoveSurplus)
    {
      NCollection_Sequence<occ::handle<ChFiDS_SurfData>>& aSeq =
        stripe->ChangeSetOfSurfData()->ChangeSequence();
      auto isSurplus = [&](const int theI) {
        return aSeq.Value(theI)->IndexOfS1() == 0 || aSeq.Value(theI)->IndexOfS2() == 0;
      };
      if (isfirst)
      {
        for (; num < aSeq.Length() && isSurplus(num);)
        {
          aSeq.Remove(num);
        }
      }
      else
      {
        for (; num > 1 && isSurplus(num);)
        {
          aSeq.Remove(num);
          num--;
        }
      }
    }
    Fd                            = stripe->SetOfSurfData()->Sequence().Value(num);
    const ChFiDS_CommonPoint& CV1 = Fd->Vertex(isfirst, 1);
    const ChFiDS_CommonPoint& CV2 = Fd->Vertex(isfirst, 2);
    // Is it always degenerated ?
    if (CV1.Point().IsEqual(CV2.Point(), 0))
    {
      nondegenere = false;
    }
    else
    {
      toujoursdegenere = false;
    }
  }

  // calcul du nombre de faces = nombre d'aretes
  /*  NCollection_List<TopoDS_Shape>::Iterator ItF,JtF,ItE;
    int nbf = 0, jf = 0;
    for (ItF.Initialize(myVFMap(Vtx)); ItF.More(); ItF.Next()){
      jf++;
      int kf = 1;
      const TopoDS_Shape& cur = ItF.Value();
      for (JtF.Initialize(myVFMap(Vtx)); JtF.More() && (kf < jf); JtF.Next(), kf++){
        if(cur.IsSame(JtF.Value())) break;
      }
      if(kf == jf) nbf++;
    }
    int nba=myVEMap(Vtx).Extent();
    for (ItE.Initialize(myVEMap(Vtx)); ItE.More(); ItE.Next()){
      const TopoDS_Edge& cur = TopoDS::Edge(ItE.Value());
      if (BRep_Tool::Degenerated(cur)) nba--;
    }
    nba=nba/2;*/
  int nba = ChFi3d_NumberOfSharpEdges(Vtx, myVEMap, myEFMap);

  // a closed fillet ending twice here, short of the vertex (see PinchEnds)
  if (nondegenere && i == 2
      && PerformPinchEnds(Vtx,
                          myVDataMap(Index),
                          myDS->ChangeDS(),
                          myEFMap,
                          myVEMap,
                          tolapp3d,
                          myEVIMap))
  {
    return;
  }

  // A setback corner: the stripes ending here are cut back and the opening
  // is filled by one patch, whatever the count of stripes and edges.
  if (nondegenere && nba >= 3 && HasSetbackAt(Index))
  {
    PerformMoreThreeCorner(Index, i);
    return;
  }

  if (nondegenere)
  { // Normal processing
    // the lines along a split of the wall, on the far piece (see LineOnSplit)
    LinesOnSplit aLines(myDS->ChangeDS(), 10 * tolapp3d);
    if (i >= 2)
    {
      for (It.Initialize(myVDataMap(Index)); It.More(); It.Next())
      {
        aLines.Show(Vtx, It.Value(), myEFMap);
      }
    }
    switch (i)
    {
      case 1: {
        if (sp->Status(isfirst) == ChFiDS_FreeBoundary)
        {
          return;
        }
        // OnSame at four sharp edges: a corner of three once a face of the
        // spine and the face it runs on into, tangent, are taken as one
        // (PerformExtremity)
        if (nba > 3 && !(nba == 4 && sp->Status(isfirst) == ChFiDS_OnSame))
        {
#ifdef OCCT_DEBUG // perf
          ChFi3d_InitChron(cl_performatend);
#endif
          PerformIntersectionAtEnd(Index);
#ifdef OCCT_DEBUG
          ChFi3d_ResultChron(cl_performatend, t_performatend);
#endif
        }
        else
        {
#ifdef OCCT_DEBUG // perf
          ChFi3d_InitChron(cl_perform1corner);
#endif
          if (ChFi3d_PerformTangentEnd(Vtx,
                                       myVDataMap(Index).First(),
                                       myDS->ChangeDS(),
                                       myEFMap,
                                       tolapp3d))
          {
            // the last section is the cut
          }
          else if (MoreSurfdata(Index))
          {
            PerformMoreSurfdata(Index);
          }
          else
          {
            PerformOneCorner(Index);
          }
#ifdef OCCT_DEBUG // perf
          ChFi3d_ResultChron(cl_perform1corner, t_perform1corner);
#endif
        }
      }
      break;
      case 2: {
        if (nba > 3)
        {
          // An end OnSame at four sharp edges is a corner of three across a
          // tangent split (PerformExtremity), made for one fillet ending
          // there; the plate of two stripes does not take such an end (the
          // result inside out). Refused, Compute runs again with that end a
          // break point, as before.
          for (It.Initialize(myVDataMap(Index)); It.More() && nba == 4; It.Next())
          {
            const occ::handle<ChFiDS_Spine>& aSp = It.Value()->Spine();
            if ((aSp->FirstVertex().IsSame(Vtx) && aSp->FirstStatus() == ChFiDS_OnSame)
                || (aSp->LastVertex().IsSame(Vtx) && aSp->LastStatus() == ChFiDS_OnSame))
            {
              throw Standard_Failure("Corner of two stripes at an end across a tangent split");
            }
          }
          if (nba == 4 && OppositeEnds(Vtx, myVDataMap(Index)))
          {
            // two fillets on no common face: the plate here, the opposite
            // corner where the plate's result was invalid
            // (ComputeOppositeCorners)
            if (THE_OPPOSITE_CORNERS.Contains(Vtx))
            {
              PerformOppositeCorner(Index);
              break;
            }
            THE_OPPOSITE_SEEN.Add(Vtx);
          }
#ifdef OCCT_DEBUG // perf
          ChFi3d_InitChron(cl_performmore3corner);
#endif
          PerformMoreThreeCorner(Index, i);
#ifdef OCCT_DEBUG // perf
          ChFi3d_ResultChron(cl_performmore3corner, t_performmore3corner);
#endif
        }
        else
        {
#ifdef OCCT_DEBUG // perf
          ChFi3d_InitChron(cl_perform2corner);
#endif
          PerformTwoCorner(Index);
#ifdef OCCT_DEBUG // perf
          ChFi3d_ResultChron(cl_perform2corner, t_perform2corner);
#endif
        }
      }
      break;
      case 3: {
        if (nba > 3)
        {
#ifdef OCCT_DEBUG // perf
          ChFi3d_InitChron(cl_performmore3corner);
#endif
          PerformMoreThreeCorner(Index, i);
#ifdef OCCT_DEBUG // perf
          ChFi3d_ResultChron(cl_performmore3corner, t_performmore3corner);
#endif
        }
        else
        {
#ifdef OCCT_DEBUG // perf
          ChFi3d_InitChron(cl_perform3corner);
#endif
          PerformThreeCorner(Index);
#ifdef OCCT_DEBUG // perf
          ChFi3d_ResultChron(cl_perform3corner, t_perform3corner);
#endif
        }
      }
      break;
      default: {
#ifdef OCCT_DEBUG // perf
        ChFi3d_InitChron(cl_performmore3corner);
#endif
        PerformMoreThreeCorner(Index, i);
#ifdef OCCT_DEBUG // perf
        ChFi3d_ResultChron(cl_performmore3corner, t_performmore3corner);
#endif
      }
    }
    aLines.Finish();
  }
  else
  { // Single case processing
    if (toujoursdegenere)
    {
      PerformSingularCorner(Index);
    }
    else
    {
      PerformMoreThreeCorner(Index, i); // Last chance...
    }
  }
}

//=================================================================================================

void ChFi3d_Builder::Reset()
{
  done = false;
  // the input as the computation found it (THE_INPUT_TOLERANCES)
  for (NCollection_DataMap<TopoDS_Shape, double, TopTools_ShapeMapHasher>::Iterator itT(
         THE_INPUT_TOLERANCES);
       itT.More();
       itT.Next())
  {
    const TopoDS_Shape& aS = itT.Key();
    if (aS.ShapeType() == TopAbs_VERTEX)
    {
      occ::handle<BRep_TVertex> aTV = occ::down_cast<BRep_TVertex>(aS.TShape());
      if (!aTV.IsNull() && !aTV->Immutable() && aTV->Tolerance() != itT.Value())
      {
        aTV->Tolerance(itT.Value());
      }
    }
    else
    {
      occ::handle<BRep_TEdge> aTE = occ::down_cast<BRep_TEdge>(aS.TShape());
      if (!aTE.IsNull() && !aTE->Immutable() && aTE->Tolerance() != itT.Value())
      {
        aTE->Tolerance(itT.Value());
      }
    }
  }
  myVDataMap.Clear();
  myRegul.Clear();
  myEVIMap.Clear();
  badstripes.Clear();
  badvertices.Clear();

  NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator itel;
  for (itel.Initialize(myListStripe); itel.More();)
  {
    if (!itel.Value()->Spine().IsNull())
    {
      itel.Value()->Reset();
      itel.Next();
    }
    else
    {
      myListStripe.Remove(itel);
    }
  }
}

//=================================================================================================

const NCollection_List<TopoDS_Shape>& ChFi3d_Builder::Generated(const TopoDS_Shape& EouV)
{
  myGenerated.Clear();
  if (EouV.IsNull())
  {
    return myGenerated;
  }
  if (EouV.ShapeType() != TopAbs_EDGE && EouV.ShapeType() != TopAbs_VERTEX)
  {
    return myGenerated;
  }
  if (myEVIMap.IsBound(EouV))
  {
    const NCollection_List<int>&    L = myEVIMap.Find(EouV);
    NCollection_List<int>::Iterator IL;
    for (IL.Initialize(L); IL.More(); IL.Next())
    {
      int                                      I  = IL.Value();
      const NCollection_List<TopoDS_Shape>&    LS = myCoup->NewFaces(I);
      NCollection_List<TopoDS_Shape>::Iterator ILS;
      for (ILS.Initialize(LS); ILS.More(); ILS.Next())
      {
        myGenerated.Append(ILS.Value());
      }
    }
  }
  return myGenerated;
}


