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
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopOpeBRepBuild_HBuilder.hxx>
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
#include <BRep_Builder.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepLProp_SLProps.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Curve.hxx>
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
// An edge the fillet made whose curve ends off its vertex by more than the
// edge's tolerance, but within the vertex's, takes that distance: at a corner
// the extension of an edge ends at a point computed on the fillet's line, not
// on the extension, and the cut there leaves it tangent -- BRepCheck's wire
// check excuses a crossing near the vertex only while both edges keep within
// twice their tolerance of the chord from the vertex, which the vertex's
// offset alone exceeded (#631's Fillet002: 3.7e-7 off, tolerance 1e-7).
// Edges of the input are left alone: the input shares them.

static void ChFi3d_EdgesCoverTheirEnds(const TopoDS_Shape& theResult, const TopoDS_Shape& theInput)
{
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> anOld;
  for (TopExp_Explorer ex(theInput, TopAbs_EDGE); ex.More(); ex.Next())
  {
    anOld.Add(ex.Current());
  }
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
    const TopLoc_Location& aLoc = anE.Location();
    double                 aTol = BRep_Tool::Tolerance(anE);
    for (TopExp_Explorer exv(anE, TopAbs_VERTEX); exv.More(); exv.Next())
    {
      const TopoDS_Vertex& aV   = TopoDS::Vertex(exv.Current());
      gp_Pnt               aP   = aC->Value(BRep_Tool::Parameter(aV, anE));
      aP.Transform(aLoc.Transformation());
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
// the outer one. Returns the history of the change, null when there was
// none.

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
// vertices there with no edge between them, one of them the input's, at one
// point of the face's parameters, and each loop of the wire from there
// bounds a region of its own.
bool FindPinch(const TopoDS_Face&                                                   theF,
               const NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>& theInputV,
               PinchedFace&                                                         theP)
{
  TopoDS_Wire aW;
  int         aNbW = 0;
  for (TopExp_Explorer exW(theF, TopAbs_WIRE); exW.More(); exW.Next(), aNbW++)
  {
    aW = TopoDS::Wire(exW.Current());
  }
  if (aNbW != 1)
  {
    return false;
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
        TopoDS_Vertex aE1, aE2;
        TopExp::Vertices(TopoDS::Edge(exE.Current()), aE1, aE2);
        isJoined = (aE1.IsSame(aVi) && aE2.IsSame(aVj)) || (aE1.IsSame(aVj) && aE2.IsSame(aVi));
      }
      if (!isJoined)
      {
        aNbPinch++;
        aV1 = aVi;
        aV2 = aVj;
      }
    }
  }
  if (aNbPinch != 1 || theInputV.Contains(aV1) == theInputV.Contains(aV2))
  {
    return false;
  }
  theP.Face = theF;
  theP.Keep = theInputV.Contains(aV1) ? aV1 : aV2;
  theP.Drop = theInputV.Contains(aV1) ? aV2 : aV1;
  theP.Ordered.Clear();
  for (BRepTools_WireExplorer wex(aW, theF); wex.More(); wex.Next())
  {
    theP.Ordered.Append(wex.Current());
  }
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
  // each loop a region, not a hole touching the other
  const int aNb = theP.Ordered.Length();
  for (int aLoop = 0; aLoop < 2; aLoop++)
  {
    const TopoDS_Wire aLW =
      aLoop == 0 ? PinchLoop(theP, theP.K1, theP.K2 - theP.K1, nullptr)
                 : PinchLoop(theP, theP.K2, aNb - theP.K2 + theP.K1, nullptr);
    TopoDS_Face aPiece = TopoDS::Face(theF.EmptyCopied());
    BRep_Builder().Add(aPiece, aLW);
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
    aDropToKeep.Bind(aP.Drop, aP.Keep);
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
  for (NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>::Iterator it(
         aDropToKeep);
       it.More();
       it.Next())
  {
    const TopoDS_Vertex& aDrop = TopoDS::Vertex(it.Key());
    const TopoDS_Vertex& aKeep = TopoDS::Vertex(it.Value());
    const double         aTol =
      BRep_Tool::Pnt(aDrop).Distance(BRep_Tool::Pnt(aKeep)) + BRep_Tool::Tolerance(aDrop);
    if (aTol > BRep_Tool::Tolerance(aKeep))
    {
      aB.UpdateVertex(aKeep, aTol);
    }
    aMerge->Replace(aDrop.Oriented(TopAbs_FORWARD), aKeep.Oriented(TopAbs_FORWARD));
  }
  NCollection_List<TopoDS_Shape> anEdges, aKeeps;
  NCollection_List<double>       aParams;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aResultEdges;
  TopExp::MapShapes(theResult, TopAbs_EDGE, aResultEdges);
  for (int i = 1; i <= aResultEdges.Extent(); i++)
  {
    const TopoDS_Edge& anE = TopoDS::Edge(aResultEdges(i));
    for (TopoDS_Iterator itV(anE); itV.More(); itV.Next())
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
    for (int aLoop = 0; aLoop < 2; aLoop++)
    {
      TopoDS_Face aPiece = TopoDS::Face(aMF.EmptyCopied());
      aB.Add(aPiece,
             aLoop == 0 ? PinchLoop(aP, aP.K1, aP.K2 - aP.K1, aMerge)
                        : PinchLoop(aP, aP.K2, aNb - aP.K2 + aP.K1, aMerge));
      aB.Add(aPieces, aPiece);
    }
    aSplit->Replace(aMF, aPieces);
  }
  theResult                               = aSplit->Apply(aMerged);
  occ::handle<BRepTools_History> aHistory = aMerge->History();
  aHistory->Merge(aSplit->History());
  return aHistory;
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
// by thousandths of a thousandth.
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
    for (int j = 1; j <= aSp->NbEdges(); j++)
    {
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
      if (!aDist.IsDone() || aDist.Value() < std::max(0.5 * aMove, 1.e-6))
      {
        return false;
      }
    }
  }
  return true;
}

// The largest radius at the end of <theSp> numbered <theIE> in the spine.
double LargestRadiusAt(const occ::handle<ChFiDS_FilSpine>& theSp, const int theIE)
{
  return theSp->IsConstant(theIE) ? theSp->Radius(theIE) : theSp->MaxRadFromSeqAndLaws();
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
  // edge looser than the input's or a twentieth of the radius
  const double aMaxTol = std::max(MaxEdgeTolerance(myShape), 0.05 * aRadius);
  auto         FreshComputation = [&]() {
    Reset();
    myCoup = new TopOpeBRepBuild_HBuilder(myCoup->BuildTool());
  };
  ChFi3d_LongSpineExtension() = true;
  try
  {
    OCC_CATCH_SIGNALS
    FreshComputation();
    Compute();
    if (done && IsValidResult(myShapeResult, aMaxTol)
        && AreFilletsMade(myListStripe, myEFMap, myShapeResult))
    {
      ChFi3d_LongSpineExtension() = false;
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
      ChFi3d_LongSpineExtension() = false;
      return;
    }
  }
  ChFi3d_LongSpineExtension() = false;
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
      if (done && IsValidResult(myShapeResult, aMaxTol())
          && AreFilletsMade(myListStripe, myEFMap, myShapeResult))
      {
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
  for (int i = 1; i <= anEnds.Length(); i++)
  {
    const SetbackEnd& anEnd = anEnds.Value(i);
    anEnd.Spine->SetSetback(anEnd.IsFirst, anEnd.Before);
  }
  FreshComputation();
  Compute();
}

//=================================================================================================

void ChFi3d_Builder::Compute()
{
  if (!THE_IN_SETBACK_COMPUTE)
  {
    // the computation, and where it fails at a vertex, the setback fallback;
    // then, on edges short for their radius, the spines extended further
    // (ComputeLongExtension). A failure the fallbacks do not mend is the
    // computation's own
    SetbackComputeGuard aGuard;
    std::exception_ptr  aFailure;
    try
    {
      Compute();
    }
    catch (Standard_Failure const&)
    {
      aFailure = std::current_exception();
      done     = false;
    }
    ComputeSetbackFallback();
    if (!done)
    {
      ComputeLongExtension();
    }
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
            double               tolv = BRep_Tool::Tolerance(v);
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
              B1.UpdateVertex(v, tolcv);
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
            double               tolv = BRep_Tool::Tolerance(v);
            if (tole > tolv)
            {
              B1.UpdateVertex(v, tole);
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
        ChFi3d_EdgesCoverTheirEnds(myShapeResult, myShape);
        // the builder's history, which Generated and BRepFilletAPI read,
        // as the result has its shapes now
        const occ::handle<BRepTools_History> aTouch =
          ChFi3d_SplitPinchedFaces(myShapeResult, myShape);
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

    for (iF = 1; iF <= aNbSurfaces; ++iF)
    {
      const NCollection_List<TopoDS_Shape>& aLF = myCoup->NewFaces(iF);
      aIt.Initialize(aLF);
      for (; aIt.More(); aIt.Next())
      {
        const TopoDS_Shape& aF = aIt.Value();
        BRepLib::SameParameter(aF, SameParTol, true);
        ShapeFix::SameParameter(aF, false, SameParTol);
      }
    }
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
  for (It.Initialize(myVDataMap(Index)), i = 0; It.More(); It.Next(), i++)
  {
    stripe = It.Value();
    sp     = stripe->Spine();
    // SurfData and its CommonPoints,
    int sens                      = 0;
    int num                       = ChFi3d_IndexOfSurfData(Vtx, stripe, sens);
    isfirst                       = (sens == 1);
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

  // A setback corner: the stripes ending here are cut back and the opening
  // is filled by one patch, whatever the count of stripes and edges.
  if (nondegenere && nba >= 3 && HasSetbackAt(Index))
  {
    PerformMoreThreeCorner(Index, i);
    return;
  }

  if (nondegenere)
  { // Normal processing
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
          if (MoreSurfdata(Index))
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

