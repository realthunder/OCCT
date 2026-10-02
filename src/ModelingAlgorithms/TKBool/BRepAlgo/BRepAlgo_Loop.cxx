// Created on: 1995-11-10
// Created by: Yves FRICAUD
// Copyright (c) 1995-1999 Matra Datavision
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

#include <vector>
#include <algorithm>

#include <BRep_Builder.hxx>
#include <BRep_TEdge.hxx>
#include <BRep_Tool.hxx>
#include <BRep_TVertex.hxx>
#include <BRepAlgo_AsDes.hxx>
#include <BRepAlgo_FaceRestrictor.hxx>
#include <BRepAlgo_Loop.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <Standard_ErrorHandler.hxx>
#include <BRepLib_MakeWire.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <IntTools_FClass2d.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomLib.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Ax2.hxx>
#include <Precision.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeFix_Wire.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Shape.hxx>
#include <TopTools.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_DataMap.hxx>
#include <NCollection_List.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <NCollection_IndexedMap.hxx>
#include <NCollection_Map.hxx>
#include <NCollection_Sequence.hxx>

#include <cstdio>
#include <cstdlib>
#include <cstring>
// #define OCCT_DEBUG_ALGO
#ifdef OCCT_DEBUG_ALGO
bool         AffichLoop = true;
int          NbLoops    = 0;
int          NbWires    = 1;
static char* name       = new char[100];
#endif

static thread_local int _CollectingEdges;
static thread_local
  NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    _EdgeMap;

BRepAlgo_LoopIntersectingEdgeMap::BRepAlgo_LoopIntersectingEdgeMap()
{
  ++_CollectingEdges;
}

BRepAlgo_LoopIntersectingEdgeMap::~BRepAlgo_LoopIntersectingEdgeMap()
{
  if (--_CollectingEdges == 0)
  {
    _EdgeMap.Clear();
  }
}

NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>&
  BRepAlgo_LoopIntersectingEdgeMap::EdgeMap()
{
  return _EdgeMap;
}

//=================================================================================================

BRepAlgo_Loop::BRepAlgo_Loop()
    : myTolConf(0.001)
{
}

//=================================================================================================

void BRepAlgo_Loop::Init(const TopoDS_Face& F)
{
  myConstEdges.Clear();
  myEdges.Clear();
  myVerOnEdges.Clear();
  myNewWires.Clear();
  myNewFaces.Clear();
  myCutEdges.Clear();
  myKeptEdges.Clear();
  myFace = F;
}

//=======================================================================
// function : Bubble
// purpose  : Orders the sequence of vertices by increasing parameter.
//=======================================================================

static void Bubble(const TopoDS_Edge&                  E,
                   NCollection_Sequence<TopoDS_Shape>& Seq,
                   NCollection_Sequence<double>&       SeqU)
{
  double        U;
  TopoDS_Vertex V;

  if (Seq.IsEmpty())
  {
    return;
  }

  for (int i = 1; i <= Seq.Length(); i++)
  {
    TopoDS_Shape aLocalV = Seq(i).Oriented(TopAbs_INTERNAL);
    V                    = TopoDS::Vertex(aLocalV);
    U                    = BRep_Tool::Parameter(V, E);
    SeqU.Append(U);
  }

  // Remove duplicates
  for (int i = 1; i < Seq.Length(); i++)
  {
    for (int j = i + 1; j <= Seq.Length(); j++)
    {
      if (Seq(i) == Seq(j) && SeqU(i) == SeqU(j))
      {
        Seq.Remove(j);
        SeqU.Remove(j);
        j--;
      }
    }
  }

  bool Invert   = true;
  int  NbPoints = Seq.Length();

  while (Invert)
  {
    Invert = false;
    for (int i = 1; i < NbPoints; i++)
    {
      if (SeqU(i + 1) < SeqU(i))
      {
        Seq.Exchange(i, i + 1);
        SeqU.Exchange(i, i + 1);
        Invert = true;
      }
    }
  }
}

//=================================================================================================

void BRepAlgo_Loop::AddEdge(TopoDS_Edge& E, const NCollection_List<TopoDS_Shape>& LV)
{
  myEdges.Append(E);
  myVerOnEdges.Bind(E, LV);
  SHOW_TOPO_SHAPE(E, "AddEdge", LV);
}

//=================================================================================================

void BRepAlgo_Loop::KeepPieces(const TopoDS_Edge& E)
{
  myKeptEdges.Add(E);
}

//=================================================================================================

void BRepAlgo_Loop::AddConstEdge(const TopoDS_Edge& E)
{
  myConstEdges.Append(E);
}

//=================================================================================================

void BRepAlgo_Loop::AddConstEdges(const NCollection_List<TopoDS_Shape>& LE)
{
  NCollection_List<TopoDS_Shape>::Iterator itl(LE);
  for (; itl.More(); itl.Next())
  {
    myConstEdges.Append(itl.Value());
  }
}

//=================================================================================================

void BRepAlgo_Loop::SetImageVV(const BRepAlgo_Image& theImageVV)
{
  myImageVV = theImageVV;
}

//=======================================================================
// function : UpdateClosedEdge
// purpose  : If the first or the last vertex of intersection
//           coincides with the closing vertex, it is removed from SV.
//           it will be added at the beginning and the end of SV by the caller.
//=======================================================================

static TopoDS_Vertex UpdateClosedEdge(const TopoDS_Edge&                  E,
                                      NCollection_Sequence<TopoDS_Shape>& SV,
                                      NCollection_Sequence<double>&       SU)
{
  TopoDS_Vertex VB[2], V1, V2, VRes;
  gp_Pnt        P, PC;
  bool          OnStart = false, OnEnd = false;
  //// modified by jgv, 13.04.04 for OCC5634 ////
  TopExp::Vertices(E, V1, V2);
  double Tol = BRep_Tool::Tolerance(V1);
  ///////////////////////////////////////////////

  if (SV.IsEmpty())
  {
    return VRes;
  }

  VB[0] = TopoDS::Vertex(SV.First());
  VB[1] = TopoDS::Vertex(SV.Last());
  PC    = BRep_Tool::Pnt(V1);

  for (int i = 0; i < 2; i++)
  {
    P = BRep_Tool::Pnt(VB[i]);
    if (P.IsEqual(PC, Tol))
    {
      VRes = VB[i];
      if (i == 0)
      {
        OnStart = true;
      }
      else
      {
        OnEnd = true;
      }
    }
  }
  if (OnStart && OnEnd)
  {
    if (!VB[0].IsSame(VB[1]))
    {
#ifdef OCCT_DEBUG_ALGO
      if (AffichLoop)
        std::cout << "Two different vertices on the closing vertex" << std::endl;
#endif
    }
    else
    {
      SV.Remove(1);
      SU.Remove(1);
      if (!SV.IsEmpty())
      {
        SV.Remove(SV.Length());
        SU.Remove(SU.Length());
      }
    }
  }
  else if (OnStart)
  {
    SV.Remove(1);
    SU.Remove(1);
  }
  else if (OnEnd)
  {
    SV.Remove(SV.Length());
    SU.Remove(SU.Length());
  }

  return VRes;
}

//=================================================================================================

static void PurgeNewEdges(
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>&
                                                                NewEdges,
  const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& UsedEdges,
  const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& theKept)
{
  for (int ii = 1; ii <= NewEdges.Extent(); ++ii)
  {
    if (theKept.Contains(NewEdges.FindKey(ii)))
    {
      continue;
    }
    NCollection_List<TopoDS_Shape>&          LNE = NewEdges.ChangeFromIndex(ii);
    NCollection_List<TopoDS_Shape>::Iterator itL(LNE);
    while (itL.More())
    {
      const TopoDS_Shape& NE = itL.Value();
      if (!UsedEdges.Contains(NE))
      {
        LNE.Remove(itL);
      }
      else
      {
        itL.Next();
      }
    }
  }
}

//=================================================================================================

// A band of a periodic face between two closed edges joined by a piece of
// the seam (theSeam, oriented from theNear's vertex to theFar's): the piece
// one way, the far edge, the piece back, the near edge. The seam's two
// pcurves decide which way round the band runs, so each closed edge is taken
// the way that runs on from them in (u, v) -- the bands either side of a
// closed edge take it opposite ways, whatever way it was stored -- and its
// pcurve is moved a whole period onto the band when it lies one over. Null
// when either closed edge cannot be fitted; ShapeFix_Wire is the fallback,
// but it edits the seam's pcurves in place to suit the wire in hand, which
// leaves the band beside it running the wrong way.

static TopoDS_Wire MakeSeamBand(const TopoDS_Edge& theSeam,
                                const TopoDS_Edge& theNear,
                                const TopoDS_Edge& theFar,
                                const TopoDS_Face& theFace)
{
  TopLoc_Location                  aLoc;
  const occ::handle<Geom_Surface>& aSurf = BRep_Tool::Surface(theFace, aLoc);
  const double aPU = aSurf->IsUPeriodic() ? aSurf->UPeriod() : 0.;
  const double aPV = aSurf->IsVPeriodic() ? aSurf->VPeriod() : 0.;
  const double aTol = 1.e-3 * std::max(1., std::max(aPU, aPV));

  auto ends = [&theFace](const TopoDS_Edge& theE, gp_Pnt2d& theFirst, gp_Pnt2d& theLast) {
    double                    aF, aL;
    occ::handle<Geom2d_Curve> aC = BRep_Tool::CurveOnSurface(theE, theFace, aF, aL);
    if (aC.IsNull())
    {
      return false;
    }
    theFirst = aC->Value(aF);
    theLast  = aC->Value(aL);
    if (theE.Orientation() == TopAbs_REVERSED)
    {
      std::swap(theFirst, theLast);
    }
    return true;
  };
  // The closed edge, one way or the other and moved by whole periods, that
  // runs from theStart to theEnd.
  auto fit = [&](const TopoDS_Edge& theE,
                 const gp_Pnt2d&    theStart,
                 const gp_Pnt2d&    theEnd,
                 TopoDS_Edge&       theFitted) {
    for (int i = 0; i < 2; ++i)
    {
      const TopoDS_Edge anE = i == 0 ? theE : TopoDS::Edge(theE.Reversed());
      gp_Pnt2d          aFirst, aLast;
      if (!ends(anE, aFirst, aLast))
      {
        return false;
      }
      const double aKU = aPU > 0. ? std::round((theStart.X() - aFirst.X()) / aPU) : 0.;
      const double aKV = aPV > 0. ? std::round((theStart.Y() - aFirst.Y()) / aPV) : 0.;
      const gp_Vec2d aShift(aKU * aPU, aKV * aPV);
      if (aFirst.Translated(aShift).Distance(theStart) > aTol
          || aLast.Translated(aShift).Distance(theEnd) > aTol)
      {
        continue;
      }
      if (aKU != 0. || aKV != 0.)
      {
        double                    aF, aL;
        occ::handle<Geom2d_Curve> aC = BRep_Tool::CurveOnSurface(anE, theFace, aF, aL);
        occ::handle<Geom2d_Curve> aMoved =
          occ::down_cast<Geom2d_Curve>(aC->Translated(aShift));
        BRep_Builder().UpdateEdge(anE, aMoved, theFace, BRep_Tool::Tolerance(anE));
      }
      theFitted = anE;
      return true;
    }
    return false;
  };

  const TopoDS_Edge aBack = TopoDS::Edge(theSeam.Reversed());
  gp_Pnt2d          aOutFirst, aOutLast, aBackFirst, aBackLast;
  TopoDS_Edge       aFar, aNear;
  if (!ends(theSeam, aOutFirst, aOutLast) || !ends(aBack, aBackFirst, aBackLast)
      || aOutFirst.Distance(aBackLast) <= aTol // not a seam on this face
      || !fit(theFar, aOutLast, aBackFirst, aFar) || !fit(theNear, aBackLast, aOutFirst, aNear))
  {
    return TopoDS_Wire();
  }
  TopoDS_Wire  aWire;
  BRep_Builder aB;
  aB.MakeWire(aWire);
  aB.Add(aWire, theSeam);
  aB.Add(aWire, aFar);
  aB.Add(aWire, aBack);
  aB.Add(aWire, aNear);
  aWire.Closed(true);
  return aWire;
}

//=================================================================================================

// Whether the wire closes in the face's UV space: its edges' pcurves, each
// taken the way the wire runs it, add up to no displacement. A wire running
// once round a periodic surface adds up to a period.
static bool IsClosedInUV(const TopoDS_Wire& theWire, const TopoDS_Face& theFace)
{
  TopLoc_Location                  aLoc;
  const occ::handle<Geom_Surface>& aSurf = BRep_Tool::Surface(theFace, aLoc);
  gp_XY                            aSum(0., 0.);
  for (TopExp_Explorer anExp(theWire, TopAbs_EDGE); anExp.More(); anExp.Next())
  {
    const TopoDS_Edge&        anEdge = TopoDS::Edge(anExp.Current());
    double                    aF, aL;
    occ::handle<Geom2d_Curve> aC2d = BRep_Tool::CurveOnSurface(anEdge, theFace, aF, aL);
    if (aC2d.IsNull())
    {
      return true;
    }
    gp_XY aD = aC2d->Value(aL).XY() - aC2d->Value(aF).XY();
    if (anEdge.Orientation() == TopAbs_REVERSED)
    {
      aD.Reverse();
    }
    aSum += aD;
  }
  if (aSurf->IsUPeriodic() && std::abs(aSum.X()) > aSurf->UPeriod() / 2.)
  {
    return false;
  }
  if (aSurf->IsVPeriodic() && std::abs(aSum.Y()) > aSurf->VPeriod() / 2.)
  {
    return false;
  }
  return true;
}

//=================================================================================================

// Whether the edges of a periodic face all lie within one period, none of
// them a seam: then its (u, v) is a chart of the whole network, as a plane's
// is, and the network's minimal wires are found the same way -- a fillet
// removed from a box, cut by its neighbours' offsets and the tubes beside it
// within its quarter turn.
static bool IsInOnePeriod(
  const NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>&
                                   theMVE,
  const TopoDS_Face&               theFace,
  const occ::handle<Geom_Surface>& theSurf)
{
  double aUMin = RealLast(), aUMax = RealFirst(), aVMin = RealLast(), aVMax = RealFirst();
  for (int i = 1; i <= theMVE.Extent(); ++i)
  {
    for (NCollection_List<TopoDS_Shape>::Iterator it(theMVE(i)); it.More(); it.Next())
    {
      const TopoDS_Edge& anEdge = TopoDS::Edge(it.Value());
      if (BRep_Tool::IsClosed(anEdge, theFace))
      {
        return false;
      }
      double                    aF, aL;
      occ::handle<Geom2d_Curve> aC2d = BRep_Tool::CurveOnSurface(anEdge, theFace, aF, aL);
      if (aC2d.IsNull())
      {
        return false;
      }
      const int aNbS = 16;
      for (int k = 0; k <= aNbS; ++k)
      {
        const gp_Pnt2d aP = aC2d->Value(aF + (aL - aF) * k / aNbS);
        aUMin             = std::min(aUMin, aP.X());
        aUMax             = std::max(aUMax, aP.X());
        aVMin             = std::min(aVMin, aP.Y());
        aVMax             = std::max(aVMax, aP.Y());
      }
    }
  }
  // Short of a whole period by more than the sampling can miss.
  if (theSurf->IsUPeriodic() && aUMax - aUMin > 0.9 * theSurf->UPeriod())
  {
    return false;
  }
  if (theSurf->IsVPeriodic() && aVMax - aVMin > 0.9 * theSurf->VPeriod())
  {
    return false;
  }
  return aUMin <= aUMax;
}

//=================================================================================================

static void StoreInMVE(
  const TopoDS_Face& F,
  TopoDS_Edge&       E,
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>&
                                                                            MVE,
  bool&                                                                     YaCouture,
  NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>& VerticesForSubstitute,
  const double                                                              theTolConf)
{
  TopoDS_Vertex                  V1, V2, V;
  NCollection_List<TopoDS_Shape> Empty;

  gp_Pnt       P1, P;
  BRep_Builder BB;
  for (int iV = 1; iV <= MVE.Extent(); iV++)
  {
    V = TopoDS::Vertex(MVE.FindKey(iV));
    P = BRep_Tool::Pnt(V);
    NCollection_List<TopoDS_Shape> VList;
    TopoDS_Iterator                VerExp(E);
    for (; VerExp.More(); VerExp.Next())
    {
      VList.Append(VerExp.Value());
    }
    NCollection_List<TopoDS_Shape>::Iterator itl(VList);
    for (; itl.More(); itl.Next())
    {
      V1 = TopoDS::Vertex(itl.Value());
      P1 = BRep_Tool::Pnt(V1);
      if (P.IsEqual(P1, theTolConf) && !V.IsSame(V1))
      {
        V.Orientation(V1.Orientation());
        if (VerticesForSubstitute.IsBound(V1))
        {
          TopoDS_Shape OldNewV = VerticesForSubstitute(V1);
          if (!OldNewV.IsSame(V))
          {
            VerticesForSubstitute.Bind(OldNewV, V);
            VerticesForSubstitute(V1) = V;
          }
        }
        else
        {
          if (VerticesForSubstitute.IsBound(V))
          {
            TopoDS_Shape NewNewV = VerticesForSubstitute(V);
            if (!NewNewV.IsSame(V1))
            {
              VerticesForSubstitute.Bind(V1, NewNewV);
            }
          }
          else
          {
            VerticesForSubstitute.Bind(V1, V);
            NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>::Iterator
              mapit(VerticesForSubstitute);
            for (; mapit.More(); mapit.Next())
            {
              if (mapit.Value().IsSame(V1))
              {
                VerticesForSubstitute(mapit.Key()) = V;
              }
            }
          }
        }
        E.Free(true);
        BB.Remove(E, V1);
        BB.Add(E, V);
        SHOW_TOPO_SHAPE(E, "ReplaceMVEE");
        SHOW_TOPO_SHAPE(V, "ReplaceMVEV");
      }
    }
  }

  TopExp::Vertices(E, V1, V2);
  if (V1.IsNull() && V2.IsNull())
  {
    YaCouture = false;
    return;
  }
  if (!MVE.Contains(V1))
  {
    MVE.Add(V1, Empty);
  }
  MVE.ChangeFromKey(V1).Append(E);
  SHOW_TOPO_SHAPE(V1, "MVE_V1_");
  if (!V1.IsSame(V2))
  {
    if (!MVE.Contains(V2))
    {
      MVE.Add(V2, Empty);
    }
    MVE.ChangeFromKey(V2).Append(E);
    SHOW_TOPO_SHAPE(V2, "MVE_V2_");
  }
  TopLoc_Location           L;
  occ::handle<Geom_Surface> S = BRep_Tool::Surface(F, L);
  if (BRep_Tool::IsClosed(E, S, L))
  {
    MVE.ChangeFromKey(V2).Append(E.Reversed());
    SHOW_TOPO_SHAPE(V2, "MVE_ClosedV2_");
    if (!V1.IsSame(V2))
    {
      MVE.ChangeFromKey(V1).Append(E.Reversed());
      SHOW_TOPO_SHAPE(V1, "MVE_ClosedV1_");
    }
    YaCouture = true;
  }
  SHOW_TOPO_SHAPE(E, "MVE");
}

//=================================================================================================

void BRepAlgo_Loop::Perform()
{
  Perform(nullptr);
}

void BRepAlgo_Loop::Perform(const NCollection_List<TopoDS_Shape>* ContextFaces,
                            const occ::handle<BRepAlgo_AsDes>&    AsDes)
{
  NCollection_List<TopoDS_Shape>::Iterator itl, itl1, itl2, itl3;
  TopoDS_Vertex                            V1, V2, OV1, OV2;
  BRep_Builder                             B;

  //------------------------------------------------
  // Check intersection in myConstEdges which is possible when make thick solid
  // with concave removed face
  //------------------------------------------------
  NCollection_List<TopoDS_Shape> theEdges = myEdges;
  NCollection_List<TopoDS_Shape> ConstEdges;
  NCollection_List<TopoDS_Shape> IntersectingEdges;
  // The span each bounded edge's own vertices give it, before the vertices of
  // the other edges are added (see myOutsideEdges).
  NCollection_DataMap<TopoDS_Shape, std::pair<double, double>, TopTools_ShapeMapHasher> aSpans;
  myOutsideEdges.Clear();
  const bool isPlanarFace = BRepAdaptor_Surface(myFace, false).GetType() == GeomAbs_Plane;
  if (_CollectingEdges)
  {
    NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> EMap;
    NCollection_List<TopoDS_Shape>                         theVerts;
    NCollection_List<TopoDS_Shape>                         LV;
    NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> MV;

    // The const edges' vertices, as they came: where a crossing splits an
    // extended edge, the piece to keep runs from the end on one of them.
    NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aConstVertices;
    for (itl.Initialize(myConstEdges); itl.More(); itl.Next())
    {
      theEdges.Append(itl.Value());
      for (TopoDS_Iterator It(itl.Value()); It.More(); It.Next())
      {
        aConstVertices.Add(It.Value());
      }
    }
    auto isOnConstEdge = [&aConstVertices](const TopoDS_Shape& theV) {
      if (aConstVertices.Contains(theV))
      {
        return true;
      }
      const gp_Pnt aP   = BRep_Tool::Pnt(TopoDS::Vertex(theV));
      const double aTol = BRep_Tool::Tolerance(TopoDS::Vertex(theV));
      for (NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>::Iterator anIt(aConstVertices);
           anIt.More();
           anIt.Next())
      {
        const TopoDS_Vertex& aCV = TopoDS::Vertex(anIt.Key());
        if (aP.Distance(BRep_Tool::Pnt(aCV)) <= std::max(aTol, BRep_Tool::Tolerance(aCV)))
        {
          return true;
        }
      }
      return false;
    };

    // Where an extended edge -- the removed face's section, stretched far past
    // the face -- crosses another edge, the crossing does not bound that
    // edge's own span: the stretched line runs on through the neighbours.
    NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aStretchVertices;
    for (itl.Initialize(theEdges); itl.More(); itl.Next())
    {
      const NCollection_List<TopoDS_Shape>* pLV = myVerOnEdges.Seek(itl.Value());
      if (!pLV)
      {
        continue;
      }
      for (TopoDS_Iterator It(itl.Value()); It.More(); It.Next())
      {
        if (It.Value().Orientation() == TopAbs_INTERNAL)
        {
          for (itl1.Initialize(*pLV); itl1.More(); itl1.Next())
          {
            aStretchVertices.Add(itl1.Value());
          }
          break;
        }
      }
    }

    // Where two edges cross inside both, neither holds a vertex there and
    // the projections below find nothing: the rim of a removed face whose
    // tangent neighbour's tube ends on it, the tube's edge on the face
    // running through the face's own outline (half a dome, one of its two
    // coplanar side faces removed). On a plane, such a crossing gets a
    // vertex of its own, offered to both edges with the others'.
    NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
      aCrossings;
    if (isPlanarFace)
    {
      std::vector<TopoDS_Edge>                               aPlain;
      NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aSeen;
      for (itl.Initialize(theEdges); itl.More(); itl.Next())
      {
        const TopoDS_Edge aE = TopoDS::Edge(itl.Value().Oriented(TopAbs_FORWARD));
        double            aF, aL;
        if (aSeen.Add(aE) && !BRep_Tool::Degenerated(aE) && !BRep_Tool::Curve(aE, aF, aL).IsNull()
            && !Precision::IsInfinite(aF) && !Precision::IsInfinite(aL))
        {
          aPlain.push_back(aE);
        }
      }
      const double aTolX = std::max(myTolConf, Precision::Confusion());
      for (size_t i = 0; i < aPlain.size(); ++i)
      {
        for (size_t j = i + 1; j < aPlain.size(); ++j)
        {
          BRepExtrema_DistShapeShape aDist(aPlain[i], aPlain[j]);
          if (!aDist.IsDone() || aDist.Value() > aTolX)
          {
            continue;
          }
          for (int k = 1; k <= aDist.NbSolution(); ++k)
          {
            if (aDist.SupportTypeShape1(k) != BRepExtrema_IsOnEdge
                || aDist.SupportTypeShape2(k) != BRepExtrema_IsOnEdge)
            {
              continue;
            }
            // Clear of every vertex the edges have or were given.
            const gp_Pnt aP     = aDist.PointOnShape1(k);
            bool         isFree = true;
            for (int m = 0; m < 2 && isFree; ++m)
            {
              const TopoDS_Edge&             aE = aPlain[m == 0 ? i : j];
              NCollection_List<TopoDS_Shape> aLV;
              if (const NCollection_List<TopoDS_Shape>* pLV = myVerOnEdges.Seek(aE))
              {
                aLV = *pLV;
              }
              for (TopoDS_Iterator aVIt(aE); aVIt.More(); aVIt.Next())
              {
                aLV.Append(aVIt.Value());
              }
              for (NCollection_List<TopoDS_Shape>::Iterator aVIt(aLV); aVIt.More(); aVIt.Next())
              {
                const TopoDS_Vertex& aV = TopoDS::Vertex(aVIt.Value());
                if (aP.Distance(BRep_Tool::Pnt(aV)) <= std::max(aTolX, BRep_Tool::Tolerance(aV)))
                {
                  isFree = false;
                  break;
                }
              }
            }
            if (!isFree)
            {
              continue;
            }
            TopoDS_Vertex aNewV;
            B.MakeVertex(aNewV, aP, aTolX);
            SHOW_TOPO_SHAPE(aNewV, "CrossingV");
            for (int m = 0; m < 2; ++m)
            {
              const TopoDS_Edge& aE = aPlain[m == 0 ? i : j];
              if (!aCrossings.IsBound(aE))
              {
                aCrossings.Bind(aE, NCollection_List<TopoDS_Shape>());
              }
              aCrossings.ChangeFind(aE).Append(aNewV);
            }
          }
        }
      }
    }

    for (itl.Initialize(theEdges); itl.More(); itl.Next())
    {
      TopoDS_Edge anEdge = TopoDS::Edge(itl.Value().Oriented(TopAbs_FORWARD));
      // Sewn edges can be doubled or not in myConstEdges
      if (!EMap.Add(anEdge))
      {
        continue;
      }

      SHOW_TOPO_SHAPE(anEdge, "InterEdge");

      LV.Clear();
      MV.Clear();

      bool         Bounded = false;
      double       FP = 0.0, LP = 0.0;
      // The span the edge's own vertices keep: from its first FORWARD vertex
      // to its last REVERSED one, open where there is none.
      double aSpanF = -Precision::Infinite(), aSpanL = Precision::Infinite();
      TopoDS_Shape VF, VL;
      if (myVerOnEdges.IsBound(anEdge))
      {
        const NCollection_List<TopoDS_Shape>& LE = myVerOnEdges(anEdge);
        Bounded                                  = true;
        for (itl1.Initialize(LE); itl1.More(); itl1.Next())
        {
          if (!MV.Add(itl1.Value()))
          {
            continue;
          }
          if (itl1.Value().Orientation() == TopAbs_FORWARD
              || itl1.Value().Orientation() == TopAbs_REVERSED)
          {
            const TopoDS_Vertex& aVertex = TopoDS::Vertex(itl1.Value());
            double               P       = BRep_Tool::Parameter(aVertex, anEdge);
            if (aStretchVertices.Contains(aVertex))
            {
              // not the edge's own
            }
            else if (aVertex.Orientation() == TopAbs_FORWARD)
            {
              aSpanF = Precision::IsInfinite(aSpanF) ? P : std::min(aSpanF, P);
            }
            else
            {
              aSpanL = Precision::IsInfinite(aSpanL) ? P : std::max(aSpanL, P);
            }
            if (VF.IsNull())
            {
              FP = LP = P;
              VF      = aVertex.Oriented(TopAbs_FORWARD);
              VL      = aVertex.Oriented(TopAbs_REVERSED);
            }
            else if (FP > P)
            {
              VF = aVertex.Oriented(TopAbs_FORWARD);
              FP = P;
            }
            else if (LP < P)
            {
              VL = aVertex.Oriented(TopAbs_REVERSED);
              LP = P;
            }
          }
        }
      }

      MV.Clear();
      if (!VF.IsNull())
      {
        MV.Add(VF);
        LV.Append(VF);
      }
      if (!VL.IsNull() && !VL.IsSame(VF))
      {
        MV.Add(VL);
        LV.Append(VL);
      }

      bool                          Extended = false;
      double                        aF, aL;
      const occ::handle<Geom_Curve> C = BRep_Tool::Curve(anEdge, aF, aL);
      TopExp::Vertices(anEdge, V1, V2);

      for (TopoDS_Iterator It(anEdge); It.More(); It.Next())
      {
        if (It.Value().Orientation() == TopAbs_INTERNAL)
        {
          SHOW_TOPO_SHAPE(anEdge, "IEdge2");
          SHOW_TOPO_SHAPE(It.Value(), "Internal2");
          if (!VF.IsNull() && !VL.IsNull())
          {
            Extended = true;
          }
          break;
        }
      }

      for (itl1.Initialize(LV); itl1.More(); itl1.Next())
      {
        SHOW_TOPO_SHAPE(itl1.Value(), "InterV1", 1);
      }
      // Open on one side, the span is trusted on arc faces only: the corner
      // piece it is for is closed by an arc face's own end, while on a plane
      // a lone vertex can be where the extended removed face crosses the
      // edge, far from the face itself, and its orientation says nothing.
      const bool isOpenSpan = Precision::IsInfinite(aSpanF) || Precision::IsInfinite(aSpanL);
      if (Bounded && !Extended
          && (!Precision::IsInfinite(aSpanF) || !Precision::IsInfinite(aSpanL))
          && aSpanL > aSpanF && (!isOpenSpan || !isPlanarFace))
      {
        aSpans.Bind(anEdge, std::make_pair(aSpanF, aSpanL));
      }
      else if (Extended && !isPlanarFace && !C.IsNull())
      {
        // A stretched edge keeps its own ends as INTERNAL vertices, and they
        // are its span: a removed face's meridian stretched over the pole and
        // down the far side of the sphere, where the piece past the pole
        // closes a wire round what the neighbour's offset cut away.
        double aIntF = Precision::Infinite(), aIntL = -Precision::Infinite();
        int    aNbInt = 0;
        for (TopoDS_Iterator It(anEdge); It.More(); It.Next())
        {
          if (It.Value().Orientation() == TopAbs_INTERNAL)
          {
            const double aP = BRep_Tool::Parameter(TopoDS::Vertex(It.Value()), anEdge);
            aIntF           = std::min(aIntF, aP);
            aIntL           = std::max(aIntL, aP);
            ++aNbInt;
          }
        }
        if (aNbInt >= 2 && aIntL > aIntF + Precision::PConfusion())
        {
          aSpans.Bind(anEdge, std::make_pair(aIntF, aIntL));
        }
      }

      for (itl1.Initialize(theEdges); itl1.More(); itl1.Next())
      {
        const TopoDS_Edge& otherEdge = TopoDS::Edge(itl1.Value());
        // YES, we will check vertices from the same edge just like any other
        // edges in order to determine the orientation of the internal
        // vertices. So no need to the check IsSame() here.
        //
        // if (otherEdge.IsSame(anEdge))
        //   continue;

        theVerts.Clear();
        const NCollection_List<TopoDS_Shape>* pLV = myVerOnEdges.Seek(otherEdge);
        if (pLV)
        {
          theVerts = *pLV;
        }
        if (const NCollection_List<TopoDS_Shape>* pLX =
              aCrossings.Seek(otherEdge.Oriented(TopAbs_FORWARD)))
        {
          for (itl2.Initialize(*pLX); itl2.More(); itl2.Next())
          {
            theVerts.Append(itl2.Value());
          }
        }
        // An edge can come without vertices (a closed section not yet cut).
        TopExp::Vertices(otherEdge, OV1, OV2);
        if (!OV1.IsNull())
        {
          theVerts.Append(OV1);
        }
        if (!OV2.IsNull())
        {
          theVerts.Append(OV2);
        }
        for (itl2.Initialize(theVerts); itl2.More(); itl2.Next())
        {
          TopoDS_Vertex aVertex = TopoDS::Vertex(itl2.Value());
          if (!MV.Add(aVertex))
          {
            continue;
          }
          double Tol = BRep_Tool::Tolerance(aVertex);
          gp_Pnt OP  = BRep_Tool::Pnt(aVertex);
          if ((!V1.IsNull() && OP.Distance(BRep_Tool::Pnt(V1)) < Tol)
              || (!V2.IsNull() && OP.Distance(BRep_Tool::Pnt(V2)) < Tol))
          {
            continue;
          }
          if (Extended)
          {
            if (OP.Distance(BRep_Tool::Pnt(TopoDS::Vertex(VF))) < Tol
                || OP.Distance(BRep_Tool::Pnt(TopoDS::Vertex(VL))) < Tol)
            {
              continue;
            }
          }
          // Handle null curve
          if (C.IsNull())
          {
            continue;
          }
          SHOW_TOPO_SHAPE(aVertex, "CheckInterV");
          GeomAPI_ProjectPointOnCurve Proj(BRep_Tool::Pnt(aVertex), C);
          if (Proj.NbPoints() > 0)
          {
            double D              = Proj.LowerDistance();
            double P              = Proj.LowerDistanceParameter();
            auto   ReorientVertex = [&](TopoDS_Shape& V, TopAbs_Orientation Ori) {
              if (V.Orientation() == Ori)
              {
                return;
              }
              V.Orientation(Ori);
              for (itl3.Initialize(LV); itl3.More(); itl3.Next())
              {
                if (itl3.Value().IsSame(V))
                {
                  itl3.ChangeValue().Orientation(Ori);
                  SHOW_TOPO_SHAPE(V, "VertexFlip");
                  break;
                }
              }
            };
            if (C->IsPeriodic())
            {
              while (P < aF)
              {
                P += C->Period();
              }
            }
            if (D < Tol && P > aF && P < aL)
            {
              TopoDS_Shape aLocalShape;
              if (Extended)
              {
                if (P < FP)
                {
                  aLocalShape = aVertex.Oriented(TopAbs_FORWARD);
                  SHOW_TOPO_SHAPE(aLocalShape, "VertexOverF");
                  ReorientVertex(VF, TopAbs_REVERSED);
                }
                else if (P > LP)
                {
                  aLocalShape = aVertex.Oriented(TopAbs_REVERSED);
                  SHOW_TOPO_SHAPE(aLocalShape, "VertexOverR");
                  ReorientVertex(VL, TopAbs_FORWARD);
                }
                // A crossing inside the edge keeps the piece from the end on
                // a const edge -- the face that stays -- when only one end
                // is: a removed face's band between the face that stays and
                // its offset, which is the longer piece where the walls are
                // shorter than twice the thickness (a short box hollowed
                // down to its bottom). Otherwise the nearer end's piece.
                else if (isOnConstEdge(VF) != isOnConstEdge(VL) ? isOnConstEdge(VF)
                                                                : P < (FP + LP) / 2)
                {
                  aLocalShape = aVertex.Oriented(TopAbs_REVERSED);
                  ReorientVertex(VF, TopAbs_FORWARD);
                }
                else
                {
                  aLocalShape = aVertex.Oriented(TopAbs_FORWARD);
                  ReorientVertex(VL, TopAbs_REVERSED);
                }
              }
              else if (P < (aF + aL) / 2)
              {
                aLocalShape = aVertex.Oriented(TopAbs_REVERSED);
              }
              else
              {
                aLocalShape = aVertex.Oriented(TopAbs_FORWARD);
              }
              B.UpdateVertex(TopoDS::Vertex(aLocalShape), P, anEdge, Tol);
              LV.Append(aLocalShape);
              SHOW_TOPO_SHAPE(aLocalShape, "InterV", 1);
            }
            else if (D < Tol)
            {
              SHOW_TOPO_SHAPE(aVertex, "InterVSkip", 1);
            }
          }
        }
      }
      if (LV.Extent())
      {
        if (!Bounded)
        {
          IntersectingEdges.Append(anEdge);
          myVerOnEdges.Bind(anEdge, LV);
        }
        else
        {
          *myVerOnEdges.ChangeSeek(anEdge) = LV;
        }
      }
      else
      {
        // As it came, not anEdge: the loop takes a const edge's orientation
        // as the one it has in the face, and a seam wire made of two closed
        // edges running the same way does not close in UV.
        ConstEdges.Append(itl.Value());
        SHOW_TOPO_SHAPE(itl.Value(), "ConstEdge");
      }
    }
    if (ConstEdges.Extent() != myConstEdges.Extent())
    {
      myConstEdges = ConstEdges;
    }
  }

#ifdef OCCT_DEBUG_ALGO
  if (AffichLoop)
  {
    std::cout << "NewLoop" << std::endl;
    NbLoops++;
  }
#endif

  //------------------------------------------------
  // Cut edges
  //------------------------------------------------
  for (itl.Initialize(theEdges); itl.More(); itl.Next())
  {
    const TopoDS_Edge& anEdge = TopoDS::Edge(itl.Value());
    if (myCutEdges.Seek(anEdge))
    {
      continue;
    }
    NCollection_List<TopoDS_Shape>        LCE;
    const NCollection_List<TopoDS_Shape>* pVertices = myVerOnEdges.Seek(anEdge);
    if (pVertices)
    {
      bool KeepAll = true;
      if (ContextFaces && AsDes && AsDes->HasAscendant(anEdge))
      {
        const NCollection_List<TopoDS_Shape> LF    = AsDes->Ascendant(anEdge);
        int                                  Count = 0;
        for (itl1.Initialize(LF); itl1.More(); itl1.Next())
        {
          for (itl2.Initialize(*ContextFaces); itl2.More(); itl2.Next())
          {
            if (itl2.Value().IsSame(itl1.Value()) && ++Count > 1)
            {
              KeepAll = false;
              SHOW_TOPO_SHAPE(anEdge, "NoKeepAll");
              break;
            }
          }
          if (!KeepAll)
          {
            break;
          }
        }
      }
      CutEdge(anEdge, *pVertices, LCE, KeepAll);
      if (const std::pair<double, double>* pSpan = aSpans.Seek(anEdge))
      {
        for (itl1.Initialize(LCE); itl1.More(); itl1.Next())
        {
          double aF, aL;
          BRep_Tool::Range(TopoDS::Edge(itl1.Value()), aF, aL);
          if (aL <= pSpan->first + Precision::PConfusion()
              || aF >= pSpan->second - Precision::PConfusion())
          {
            myOutsideEdges.Add(itl1.Value());
            SHOW_TOPO_SHAPE(itl1.Value(), "OutsidePiece");
          }
        }
      }
      myCutEdges.Add(anEdge, LCE);
    }
  }

  FindLoop();

  if (_CollectingEdges)
  {
#if 1
    NCollection_IndexedDataMap<TopoDS_Shape,
                               NCollection_List<TopoDS_Shape>,
                               TopTools_ShapeMapHasher>::Iterator itM(myCutEdges);
    for (; itM.More(); itM.Next())
    {
      SHOW_TOPO_SHAPE(itM.Key(), "InterNewEdge", itM.Value());
      NCollection_List<TopoDS_Shape>* pLE = _EdgeMap.ChangeSeek(itM.Key());
      if (!pLE)
      {
        pLE = _EdgeMap.Bound(itM.Key(), NCollection_List<TopoDS_Shape>());
      }
      for (itl.Initialize(itM.Value()); itl.More(); itl.Next())
      {
        pLE->Append(itl.Value());
      }
    }
#else
    for (itl.Initialize(IntersectingEdges); itl.More(); itl.Next())
    {
      const NCollection_List<TopoDS_Shape>& aList = NewEdges(TopoDS::Edge(itl.Value()));
      _EdgeMap.Bind(itl.Value(), aList);
      SHOW_TOPO_SHAPE(itl.Value(), "InterNewEdge", aList);
    }
#endif
  }
}

namespace
{

// boost::hash_combine
inline size_t combine(size_t seed, size_t h) noexcept
{
  seed ^= h + 0x9e3779b9 + (seed << 6U) + (seed >> 2U);
  return seed;
}

struct WireInfo
{
  mutable TopoDS_Wire                          aWire;
  std::vector<std::pair<TopoDS_Shape, size_t>> Edges;
  size_t                                       aHashCode;
  bool                                         HasSeam;
  bool                                         Outside = false;
  mutable TopoDS_Face                          aFace;

  WireInfo(const TopoDS_Wire& W, bool theHasSeam = false)
      : aWire(W),
        HasSeam(theHasSeam)
  {
    TopoDS_Iterator aIt(W);
    for (; aIt.More(); aIt.Next())
    {
      Edges.emplace_back(aIt.Value(), std::hash<TopoDS_Shape>{}(aIt.Value()));
    }
    std::sort(Edges.begin(),
              Edges.end(),
              [](const std::pair<TopoDS_Shape, size_t>& a, const std::pair<TopoDS_Shape, size_t>& b) {
                if (a.first.TShape().get() < b.first.TShape().get())
                {
                  return true;
                }
                if (b.first.TShape().get() < a.first.TShape().get())
                {
                  return false;
                }
                return std::hash<TopLoc_Location>{}(a.first.Location())
                       < std::hash<TopLoc_Location>{}(b.first.Location());
              });
    aHashCode = 0;
    for (const auto& s : Edges)
    {
      aHashCode = combine(aHashCode, s.second);
    }
  }

  WireInfo(const WireInfo& other)
      : aWire(other.aWire),
        Edges(other.Edges),
        aHashCode(other.aHashCode),
        HasSeam(other.HasSeam),
        Outside(other.Outside)
  {
  }

  bool Contains(const TopoDS_Edge& E) const
  {
    for (const auto& s : Edges)
    {
      if (s.first.IsSame(E))
      {
        return true;
      }
    }
    return false;
  }
};

struct WireInfoHasher
{
  size_t operator()(const WireInfo& theInfo) const noexcept { return theInfo.aHashCode; }

  bool operator()(const WireInfo& a, const WireInfo& b) const noexcept
  {
    if (a.Edges.size() != b.Edges.size())
    {
      return false;
    }
    for (std::size_t i = 0; i < a.Edges.size(); ++i)
    {
      if (!a.Edges[i].first.IsSame(b.Edges[i].first))
      {
        return false;
      }
    }
    return true;
  }
};

// Indexed so iteration follows discovery order: wire emission order feeds
// FaceRestrictor input order and must not depend on hash/address order.
typedef NCollection_IndexedMap<WireInfo, WireInfoHasher>           MapOfWire;
typedef NCollection_IndexedMap<WireInfo, WireInfoHasher>::Iterator MapIteratorOfMapOfWire;

void FindAllLoops(const TopoDS_Vertex&                                                  CV,
                  const TopoDS_Edge&                                                    CE,
                  NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>& CurrentVEMap,
                  NCollection_List<TopoDS_Shape>&                                       CurrentEdgeList,
                  const NCollection_IndexedDataMap<TopoDS_Shape,
                                                   NCollection_List<TopoDS_Shape>,
                                                   TopTools_ShapeMapHasher>&            MVE,
                  MapOfWire&                                                            NewWires,
                  const TopoDS_Face&                                                    aFace,
                  const double&                                                         Tol,
                  const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>&         theOutside)
{
  NCollection_List<TopoDS_Shape>::Iterator itl;
  TopoDS_Vertex                            V1, V2, NV;

  SHOW_TOPO_SHAPE(CV, "CV");
  SHOW_TOPO_SHAPE(CE, "CE");

  CurrentVEMap.Bind(CV, CE);
  CurrentEdgeList.Prepend(CE);
  TopExp::Vertices(CE, V1, V2);
  if (CV.IsSame(V1))
  {
    NV = V2;
  }
  else
  {
    NV = V1;
  }

  TopoDS_Edge         EF;
  const TopoDS_Shape* pE = CurrentVEMap.Seek(NV);
  if (pE)
  {
    EF = TopoDS::Edge(*pE);
  }
  if (pE && (!EF.IsSame(CE) || V1.IsSame(V2)))
  {
    TopoDS_Edge      E = EF;
    TopoDS_Vertex    VF, VL;
    BRepLib_MakeWire aMakeWire;
    for (;;)
    {
      TopExp::Vertices(E, V1, V2, false);
      if (VF.IsNull())
      {
        VF = NV;
        VL = NV.IsSame(V1) ? V2 : V1;
        SHOW_TOPO_SHAPE(VF, "NWVF");
      }
      else if (VL.IsSame(V1))
      {
        VL = V2;
      }
      else
      {
        VL = V1;
      }
      SHOW_TOPO_SHAPE(VL, "NWV");
      SHOW_TOPO_SHAPE(E, "NWE");
      aMakeWire.Add(E);
      if (VL.IsSame(VF))
      {
        break;
      }
      const TopoDS_Shape* pNext = CurrentVEMap.Seek(VL);
      if (!pNext)
      {
        // The walk left the current DFS path (possible after tolerance-based
        // vertex substitutions) - this candidate cannot form a loop.
        break;
      }
      E = TopoDS::Edge(*pNext);
    }
    if (aMakeWire.IsDone())
    {
      TopoDS_Wire NW = aMakeWire.Wire();
      WireInfo    anInfo(NW);
      for (const auto& v : anInfo.Edges)
      {
        anInfo.Outside = anInfo.Outside || theOutside.Contains(v.first);
      }
      if (NW.Closed() && !NewWires.Contains(anInfo) && NewWires.Add(anInfo) > 0)
      {
        SHOW_TOPO_SHAPE(NW, "NewWire");
      }
      else
      {
        SHOW_TOPO_SHAPE(NW, "DiscardWire");
      }
    }
    else
    {
      SHOW_TOPO_SHAPE(EF, "DiscardOpenChain");
    }
  }
  else
  {
    for (itl.Initialize(MVE.FindFromKey(NV)); itl.More(); itl.Next())
    {
      const TopoDS_Edge& NE = TopoDS::Edge(itl.Value());
      if (!NE.IsSame(CE) && !EF.IsSame(NE))
      {
        FindAllLoops(NV, NE, CurrentVEMap, CurrentEdgeList, MVE, NewWires, aFace, Tol, theOutside);
      }
    }
  }

  CurrentVEMap.UnBind(CV);
  CurrentEdgeList.RemoveFirst();
  // SHOW_TOPO_SHAPE(CE, "Pop");
}

void SplitWires(NCollection_List<TopoDS_Shape>& OutputWires,
                MapOfWire&                      InputWires,
                NCollection_IndexedDataMap<TopoDS_Shape,
                                           NCollection_List<TopoDS_Shape>,
                                           TopTools_ShapeMapHasher>& MVE,
                const TopoDS_Face&              aFace)
{
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> UsedEdges;
  NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
                                                            aEdgeWireMap;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aWireEdgeMap;
  NCollection_List<TopoDS_Shape>::Iterator                  itl;
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>    aPrunedMap;
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>    aCheckMap;
  TopExp_Explorer                                           aExp;
  BRep_Builder                                              B;
  TopLoc_Location                                           L;
  double                                                    Tol = BRep_Tool::Tolerance(aFace);

  for (MapIteratorOfMapOfWire itM(InputWires); itM.More(); itM.Next())
  {
    for (const auto& v : itM.Value().Edges)
    {
      UsedEdges.Add(v.first);
    }
  }

  for (MapIteratorOfMapOfWire itM(InputWires); itM.More(); itM.Next())
  {
    IntTools_FClass2d aClassifier;
    const WireInfo&   Info  = itM.Value();
    TopoDS_Wire&      aWire = TopoDS::Wire(Info.aWire);
    if (aPrunedMap.Contains(aWire))
    {
      continue;
    }

    TopoDS_Face& aNF = Info.aFace;
    if (aNF.IsNull())
    {
      occ::handle<Geom_Surface> S = BRep_Tool::Surface(aFace, L);
      SHOW_TOPO_SHAPE(aFace, "MakeTestSurface");
      B.MakeFace(aNF, S, L, Tol);

      BRepTopAdaptor_FClass2d FClass2d(aNF, Precision::PConfusion());
      if (FClass2d.PerformInfinitePoint() == TopAbs_OUT)
      {
        B.Add(aNF, aWire);
        SHOW_TOPO_SHAPE(aWire, "MakeTestWire");
      }
      else
      {
        aWire.Reverse();
        B.Add(aNF, aWire);
        SHOW_TOPO_SHAPE(aWire, "MakeTestReverseWire");
      }
      B.NaturalRestriction(aNF, false);

      BRepCheck_Analyzer anAnalyzer(aNF);
      if (!anAnalyzer.IsValid())
      {
        SHOW_TOPO_SHAPE(aNF, "MakeTestBeforeFix");
        // The test face only classifies. It is made of the loop's own edges
        // on the face's own surface, and ShapeFix edits the edges it is given
        // in place: shifting a pcurve by a period rewrote the input face's
        // pcurve of an edge it shares with the shape, and the caller's shape
        // came back inside out (a pad's arc face under a thickness). Fix a
        // copy.
        ShapeFix_Shape aFix(BRepBuilderAPI_Copy(aNF, true).Shape());
        aFix.Perform();
        aNF = TopoDS::Face(aFix.Shape());
      }
      aNF = TopoDS::Face(aNF.Oriented(TopAbs_FORWARD));
      SHOW_TOPO_SHAPE(aNF, "MakeTestFace");

      // WARNING! There must be some bug in IntTools_FClass2d::Init() which
      // makes this class instance not reusable, i.e. Init() with other face
      // will yeild weirdly incorrect result. So we have to use a new
      // instance for each new face.
      aClassifier.Init(aNF, BRep_Tool::Tolerance(aNF));
    }

    for (aExp.Init(aWire, TopAbs_VERTEX); aExp.More(); aExp.Next())
    {
      //-----------------------------------------------
      // Prune the wire if it can be split by other wire. For a given vertex of
      // the wire, we check if there is any edge sharing this vertex has its
      // middle point inside the wire.
      //----------------------------------------------
      const TopoDS_Vertex&                  aVertex(TopoDS::Vertex(aExp.Current()));
      const NCollection_List<TopoDS_Shape>* pEL = MVE.Seek(aVertex);
      if (!pEL || pEL->Extent() <= 2)
      {
        continue;
      }

      aCheckMap.Clear();

      aWireEdgeMap.Clear();
      TopExp::MapShapes(aWire, TopAbs_EDGE, aWireEdgeMap);

      for (itl.Initialize(*pEL); itl.More(); itl.Next())
      {
        const TopoDS_Edge& aEdge = TopoDS::Edge(itl.Value());

        if (aWireEdgeMap.Contains(aEdge))
        {
          continue;
        }

        if (!UsedEdges.Contains(aEdge))
        {
          continue;
        }

        if (!aCheckMap.Add(aEdge))
        {
          continue;
        }

        SHOW_TOPO_SHAPE(aEdge, "CheckSplit");

        // Get 2d curve of the edge on the face
        double                           aT1, aT2;
        const occ::handle<Geom2d_Curve>& aC2D = BRep_Tool::CurveOnSurface(aEdge, aNF, aT1, aT2);
        if (aC2D.IsNull())
        {
          SHOW_TOPO_SHAPE(aEdge, "Prune_NoCurve_");
          continue;
        }

        // Get middle point on the curve
        gp_Pnt2d aP2D = aC2D->Value((aT1 + aT2) / 2.);

        // Classify the point
        TopAbs_State aState = aClassifier.Perform(aP2D);

        if (aClassifier.IsHole() && aState == TopAbs_OUT)
        {
          SHOW_TOPO_SHAPE(aNF, "Prune2_");
          aPrunedMap.Add(aWire);
          break;
        }
        else if (!aClassifier.IsHole() && aState == TopAbs_IN)
        {
          SHOW_TOPO_SHAPE(aNF, "Prune3_");
          aPrunedMap.Add(aWire);
          break;
        }
      }
      if (itl.More())
      {
        break;
      }
    }
  }

  for (MapIteratorOfMapOfWire itM(InputWires); itM.More(); itM.Next())
  {
    if (!aPrunedMap.Contains(itM.Value().aWire))
    {
      for (const auto& v : itM.Value().Edges)
      {
        NCollection_List<TopoDS_Shape>* pLE = aEdgeWireMap.ChangeSeek(v.first);
        if (!pLE)
        {
          pLE = aEdgeWireMap.Bound(v.first, NCollection_List<TopoDS_Shape>());
        }
        pLE->Append(itM.Value().aWire);
      }
    }
  }

  for (MapIteratorOfMapOfWire itM(InputWires); itM.More(); itM.Next())
  {
    const TopoDS_Wire& aWire = itM.Value().aWire;
    if (aPrunedMap.Contains(aWire))
    {
      continue;
    }

    //-----------------------------------------------
    // Prune the wire if all of its edges are shared by some other wire, in
    // which case the wire can be interpreted as inner hole)
    //----------------------------------------------

    aCheckMap.Clear();
    bool Pruned = true;
    for (const auto& v : itM.Value().Edges)
    {
      const NCollection_List<TopoDS_Shape>& LE = aEdgeWireMap.Find(v.first);
      if (LE.Extent() == 1)
      {
        OutputWires.Append(aWire);
        Pruned = false;
        break;
      }
      for (itl.Initialize(LE); itl.More(); itl.Next())
      {
        if (!itl.Value().IsSame(aWire))
        {
          aCheckMap.Add(itl.Value());
        }
      }
    }
    if (Pruned)
    {
      SHOW_TOPO_SHAPE(aWire, "Prune1_");
      NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>::Iterator itM1(aCheckMap);
      for (; itM1.More(); itM1.Next())
      {
        SHOW_TOPO_SHAPE(itM1.Value(), "Prune1W_");
      }
      // aPrunedMap.Add(aWire);
    }
  }
}

} // Anonymous namespace


//=================================================================================================
// The minimal wires by angle, not by search (after FreeCAD's WireJoiner).
//
// The faces of a planar edge network are the orbits of a permutation on its
// darts, a dart being one end of an edge taken as that edge travelled away
// from that end:
//
//     rev(d)   the same edge travelled the other way
//     next(d)  the first dart clockwise from rev(d) in the angular order of
//              the darts leaving the vertex d arrives at
//
// Walking next() from any dart returns to it, every dart lies in one orbit,
// and each orbit is a face boundary: the bounded faces come out with a
// positive signed area, the unbounded face of each connected component with
// a negative one. One pass, an angular sort per vertex, nothing searched and
// nothing undone. It is OCCT's own rule at a branching vertex
// (BOPAlgo_WireSplitter::ClockWiseAngle), measured here in the face's
// parametric space, which keeps the cyclic order of the directions leaving a
// vertex on a surface that is not periodic -- periodic faces keep the search
// and its seam wires.
//=================================================================================================

namespace
{
// Angles at or below this are the same direction.
constexpr double THE_ANGLE_TIE = 1.e-8;
// Where the second order tie break samples: two edges leaving a vertex the
// same way are told apart by the chord to a point this far along each.
constexpr double THE_ANGLE_REFINE = 0.1;

struct LoopDart
{
  int  Edge;    // index in the edge map
  bool AtStart; // leaves the edge's first vertex, runs along it
};

// The turn from the back of the edge we came in on to a way out, measured one
// way round; going straight back is a full turn, the last resort.
double ClockWiseAngle(const double theIn, const double theOut)
{
  double a = theIn - theOut;
  while (a < 0.)
  {
    a += 2. * M_PI;
  }
  while (a >= 2. * M_PI)
  {
    a -= 2. * M_PI;
  }
  return a <= THE_ANGLE_TIE ? 2. * M_PI : a;
}
} // namespace

// Direction a dart leaves its vertex in, as an angle in the face's (u, v); a
// fraction above zero takes the chord to a point that far along instead.
static bool DartAngle(const occ::handle<Geom2d_Curve>& theC2d,
                      const double                     theF,
                      const double                     theL,
                      const bool                       theAtStart,
                      double                           theFraction,
                      double&                          theAngle)
{
  gp_Vec2d aDir;
  if (theFraction <= 0.)
  {
    gp_Pnt2d aP;
    gp_Vec2d aD1;
    theC2d->D1(theAtStart ? theF : theL, aP, aD1);
    if (aD1.SquareMagnitude() > Precision::SquarePConfusion())
    {
      aDir = theAtStart ? aD1 : aD1.Reversed();
    }
    else
    {
      theFraction = THE_ANGLE_REFINE;
    }
  }
  if (theFraction > 0.)
  {
    const double p0 = theAtStart ? theF : theL;
    const double p1 = theAtStart ? theF + (theL - theF) * theFraction
                                 : theL - (theL - theF) * theFraction;
    aDir            = gp_Vec2d(theC2d->Value(p0), theC2d->Value(p1));
  }
  if (aDir.SquareMagnitude() <= Precision::SquarePConfusion())
  {
    return false;
  }
  theAngle = std::atan2(aDir.Y(), aDir.X());
  return true;
}

static bool FindLoopsByAngle(
  const NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>&
                                                                MVE,
  const TopoDS_Face&                                            theFace,
  MapOfWire&                                                    theWires,
  const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& theOutside)
{
  // The edges, and for each its pcurve and its two ends among the vertices.
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> anEdges;
  for (int i = 1; i <= MVE.Extent(); ++i)
  {
    for (NCollection_List<TopoDS_Shape>::Iterator it(MVE(i)); it.More(); it.Next())
    {
      anEdges.Add(it.Value().Oriented(TopAbs_FORWARD));
    }
  }
  const int                              aNbE = anEdges.Extent();
  std::vector<occ::handle<Geom2d_Curve>> aC2d(aNbE + 1);
  std::vector<double>                    aF(aNbE + 1), aL(aNbE + 1);
  std::vector<int>                       aV1(aNbE + 1), aV2(aNbE + 1);
  for (int e = 1; e <= aNbE; ++e)
  {
    const TopoDS_Edge& E = TopoDS::Edge(anEdges(e));
    aC2d[e]              = BRep_Tool::CurveOnSurface(E, theFace, aF[e], aL[e]);
    TopoDS_Vertex V1, V2;
    TopExp::Vertices(E, V1, V2);
    aV1[e] = V1.IsNull() ? 0 : MVE.FindIndex(V1);
    aV2[e] = V2.IsNull() ? 0 : MVE.FindIndex(V2);
    if (aC2d[e].IsNull() || aV1[e] == 0 || aV2[e] == 0)
    {
      return false;
    }
  }
  // Darts: 2e leaves the first vertex, 2e+1 the last; rev(d) = d ^ 1.
  const int           aNbD = 2 * (aNbE + 1);
  std::vector<double> anAngle(aNbD, 0.);
  std::vector<std::vector<int>> aOut(MVE.Extent() + 1);
  for (int e = 1; e <= aNbE; ++e)
  {
    for (int k = 0; k < 2; ++k)
    {
      const int d = 2 * e + k;
      if (!DartAngle(aC2d[e], aF[e], aL[e], k == 0, 0., anAngle[d]))
      {
        return false;
      }
      aOut[k == 0 ? aV1[e] : aV2[e]].push_back(d);
    }
  }
  auto aDartVertex  = [&](int d) { return (d & 1) ? aV2[d / 2] : aV1[d / 2]; };
  auto aArrivalVert = [&](int d) { return (d & 1) ? aV1[d / 2] : aV2[d / 2]; };
  auto aRefined     = [&](int d, double& theA) {
    return DartAngle(aC2d[d / 2], aF[d / 2], aL[d / 2], (d & 1) == 0, THE_ANGLE_REFINE, theA);
  };
  auto aNext = [&](int d) {
    const int         aBack = d ^ 1;
    const auto&       aCand = aOut[aArrivalVert(d)];
    int               aBest = -1;
    double            aBestA = 0.;
    for (int c : aCand)
    {
      const double a = ClockWiseAngle(anAngle[aBack], anAngle[c]);
      if (aBest < 0 || a < aBestA)
      {
        aBest  = c;
        aBestA = a;
      }
    }
    // Ties merge two faces into one: a second look, along the edges.
    std::vector<int> aTies;
    for (int c : aCand)
    {
      if (c != aBest && ClockWiseAngle(anAngle[aBack], anAngle[c]) <= aBestA + THE_ANGLE_TIE)
      {
        aTies.push_back(c);
      }
    }
    double anIn = 0.;
    if (!aTies.empty() && aRefined(aBack, anIn))
    {
      aTies.push_back(aBest);
      int    aRef  = -1;
      double aRefA = 0.;
      for (int c : aTies)
      {
        double a = 0.;
        if (!aRefined(c, a))
        {
          continue;
        }
        const double aTurn = ClockWiseAngle(anIn, a);
        if (aRef < 0 || aTurn < aRefA)
        {
          aRef  = c;
          aRefA = aTurn;
        }
      }
      if (aRef >= 0)
      {
        aBest = aRef;
      }
    }
    return aBest;
  };
  (void)aDartVertex;

  std::vector<int> anOrbit(aNbD, 0);
  int              aNbOrbits = 0;
  std::vector<int> aLoop;
  for (int d0 = 2; d0 < aNbD; ++d0)
  {
    if (anOrbit[d0])
    {
      continue;
    }
    ++aNbOrbits;
    aLoop.clear();
    int d = d0;
    for (;;)
    {
      anOrbit[d] = aNbOrbits;
      aLoop.push_back(d);
      const int n = aNext(d);
      if (n < 0)
      {
        return false;
      }
      if (n == d0)
      {
        break;
      }
      if (anOrbit[n])
      {
        // next() is a permutation: coming back anywhere but to the start
        // means the angles do not describe this network.
        return false;
      }
      d = n;
    }
    // Drop the out-and-back excursions: a bridge or a tail bounds nothing.
    std::vector<int> aPruned;
    for (int x : aLoop)
    {
      if (!aPruned.empty() && aPruned.back() == (x ^ 1))
      {
        aPruned.pop_back();
      }
      else
      {
        aPruned.push_back(x);
      }
    }
    while (aPruned.size() >= 2 && aPruned.front() == (aPruned.back() ^ 1))
    {
      aPruned.pop_back();
      aPruned.erase(aPruned.begin());
    }
    if (aPruned.empty() || (aPruned.size() == 1 && aV1[aPruned[0] / 2] != aV2[aPruned[0] / 2]))
    {
      continue;
    }
    // Twice the signed area in (u, v), sampled along each dart the way it
    // runs -- four points an edge, so a lone circle is not read as empty:
    // bounded faces positive.
    double   anArea = 0.;
    gp_Pnt2d aPrev;
    bool     isFirst = true;
    gp_Pnt2d aFirst;
    for (int x : aPruned)
    {
      const int e = x / 2;
      gp_Pnt2d  aPs[4];
      for (int k = 0; k < 4; ++k)
      {
        const double t = (x & 1) ? 1. - k / 4. : k / 4.;
        aPs[k]         = aC2d[e]->Value(aF[e] + (aL[e] - aF[e]) * t);
      }
      for (const gp_Pnt2d& aP : aPs)
      {
        if (isFirst)
        {
          aFirst  = aP;
          isFirst = false;
        }
        else
        {
          anArea += aPrev.X() * aP.Y() - aP.X() * aPrev.Y();
        }
        aPrev = aP;
      }
    }
    anArea += aPrev.X() * aFirst.Y() - aFirst.X() * aPrev.Y();
    if (anArea <= 0.)
    {
      continue;
    }
    BRepLib_MakeWire aMW;
    for (int x : aPruned)
    {
      aMW.Add(TopoDS::Edge(anEdges(x / 2).Oriented((x & 1) ? TopAbs_REVERSED : TopAbs_FORWARD)));
    }
    if (!aMW.IsDone() || !aMW.Wire().Closed())
    {
      return false;
    }
    WireInfo anInfo(aMW.Wire());
    for (const auto& v : anInfo.Edges)
    {
      anInfo.Outside = anInfo.Outside || theOutside.Contains(v.first);
    }
    if (!theWires.Contains(anInfo))
    {
      theWires.Add(anInfo);
      SHOW_TOPO_SHAPE(anInfo.aWire, "AngleWire");
    }
  }
  return true;
}

// Whether FindLoop walks the angles (1, the default), searches (0), or does
// both and reports where the minimal wires differ (2):
// BREPALGO_LOOP_WALK=0/1/check.
static int LoopWalkMode()
{
  const char* aVal = std::getenv("BREPALGO_LOOP_WALK");
  if (aVal == nullptr)
  {
    return 1;
  }
  if (std::strcmp(aVal, "check") == 0)
  {
    return 2;
  }
  return std::atoi(aVal) != 0 ? 1 : 0;
}

// A wire through a piece of an edge beyond the span its intersections gave it
// (the loop keeps every piece, to find the edges of a concave removed face)
// loses to a wire sharing an edge with it that stays inside: the corner of an
// arc face cut off by the arc of the next face, closed by the arc's own end.
// Alone, it is kept.
static void DropOutsideWires(MapOfWire& theWires, const TopoDS_Face& theFace)
{
  for (int iw = theWires.Extent(); iw >= 1; --iw)
  {
    const WireInfo& anInfo = theWires(iw);
    if (!anInfo.Outside || anInfo.HasSeam)
    {
      continue;
    }
    bool isRivalled = false;
    for (MapIteratorOfMapOfWire itW(theWires); itW.More() && !isRivalled; itW.Next())
    {
      const WireInfo& anOther = itW.Value();
      if (anOther.Outside)
      {
        continue;
      }
      for (const auto& v : anInfo.Edges)
      {
        if (!BRep_Tool::IsClosed(TopoDS::Edge(v.first), theFace)
            && anOther.Contains(TopoDS::Edge(v.first)))
        {
          isRivalled = true;
          break;
        }
      }
    }
    if (isRivalled)
    {
      SHOW_TOPO_SHAPE(anInfo.aWire, "OutsideWire");
      theWires.RemoveFromIndex(iw);
    }
  }
}

//=================================================================================================

void BRepAlgo_Loop::FindLoop()
{
  NCollection_List<TopoDS_Shape>::Iterator itl, itl1, itl2;
  bool                                     YaCouture = false;

  myNewWires.Clear();
  myNewFaces.Clear();

  //-----------------------------------
  // Construction map vertex => edges
  //-----------------------------------
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    MVE;

  // Degenerated edges (pole edges of spheres/cones) must not take part in the
  // vertex-topology loop search: being closed on their single vertex they come
  // out as bogus standalone one-edge wires, while the wire actually crossing
  // the singularity is left open in UV space. Collect them and re-insert each
  // into the wire passing through its vertex once the wires are built.
  NCollection_List<TopoDS_Shape> DegenEdges;

  // A piece of an extended edge (one carrying INTERNAL vertices: the context
  // extension stretches the edges of a removed face far past their ends) can
  // lie exactly on an edge the face has from elsewhere -- the tangent line of
  // an arc face, where the stretched edge of the removed face runs on through
  // the rest of the face. Both in the loop, the face gets the line twice and
  // the shell a free edge. Such a piece is left out; the other edge stays.
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aDuplicates;
  {
    auto isExtended = [](const TopoDS_Shape& theE) {
      for (TopoDS_Iterator aIt(theE); aIt.More(); aIt.Next())
      {
        if (aIt.Value().Orientation() == TopAbs_INTERNAL)
        {
          return true;
        }
      }
      return false;
    };
    NCollection_List<TopoDS_Shape> aFromExtended, aOthers;
    NCollection_IndexedDataMap<TopoDS_Shape,
                               NCollection_List<TopoDS_Shape>,
                               TopTools_ShapeMapHasher>::Iterator itC(myCutEdges);
    for (; itC.More(); itC.Next())
    {
      NCollection_List<TopoDS_Shape>& aTarget = isExtended(itC.Key()) ? aFromExtended : aOthers;
      for (itl1.Initialize(itC.Value()); itl1.More(); itl1.Next())
      {
        aTarget.Append(itl1.Value());
      }
    }
    for (itl1.Initialize(myConstEdges); itl1.More(); itl1.Next())
    {
      (isExtended(itl1.Value()) ? aFromExtended : aOthers).Append(itl1.Value());
    }
    const double aTol = std::max(myTolConf, Precision::Confusion());
    auto         aKey = [](const TopoDS_Edge& theE, gp_Pnt& theP1, gp_Pnt& theP2, gp_Pnt& theM) {
      TopoDS_Vertex aV1, aV2;
      TopExp::Vertices(theE, aV1, aV2);
      if (aV1.IsNull() || aV2.IsNull() || aV1.IsSame(aV2) || BRep_Tool::Degenerated(theE))
      {
        return false;
      }
      theP1 = BRep_Tool::Pnt(aV1);
      theP2 = BRep_Tool::Pnt(aV2);
      BRepAdaptor_Curve aC(theE);
      theM = aC.Value((aC.FirstParameter() + aC.LastParameter()) / 2.);
      return true;
    };
    for (itl1.Initialize(aFromExtended); itl1.More(); itl1.Next())
    {
      gp_Pnt aP1, aP2, aM;
      if (!aKey(TopoDS::Edge(itl1.Value()), aP1, aP2, aM))
      {
        continue;
      }
      for (itl2.Initialize(aOthers); itl2.More(); itl2.Next())
      {
        gp_Pnt aQ1, aQ2, aN;
        if (itl2.Value().IsSame(itl1.Value())
            || !aKey(TopoDS::Edge(itl2.Value()), aQ1, aQ2, aN) || aM.Distance(aN) > aTol)
        {
          continue;
        }
        if ((aP1.Distance(aQ1) <= aTol && aP2.Distance(aQ2) <= aTol)
            || (aP1.Distance(aQ2) <= aTol && aP2.Distance(aQ1) <= aTol))
        {
          aDuplicates.Add(itl1.Value());
          SHOW_TOPO_SHAPE(itl1.Value(), "DuplicateOfExtended");
          break;
        }
      }
    }
  }

  // A degenerated edge whose pcurve stays on one point is the pole of a
  // sphere the face is no longer parametrized by (BRepOffset_Tool::EnLargeFace
  // turns a sphere's axis off a face that has to grow round its pole): the
  // point is an ordinary one there, and the edge has no place in a wire.
  auto hasPCurve = [this](const TopoDS_Edge& theE) {
    double                          aF, aL;
    const occ::handle<Geom2d_Curve> aC2d = BRep_Tool::CurveOnSurface(theE, myFace, aF, aL);
    return !aC2d.IsNull() && aC2d->Value(aF).Distance(aC2d->Value(aL)) > Precision::PConfusion();
  };

  // add cut edges (in the order the edges were cut - hash order here would
  // make vertex canonicalization and loop discovery nondeterministic).
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> Emap;
  NCollection_IndexedDataMap<TopoDS_Shape,
                             NCollection_List<TopoDS_Shape>,
                             TopTools_ShapeMapHasher>::Iterator itM(myCutEdges);
  for (; itM.More(); itM.Next())
  {
    for (itl1.Initialize(itM.Value()); itl1.More(); itl1.Next())
    {
      TopoDS_Edge& E = TopoDS::Edge(itl1.ChangeValue());
      if (!aDuplicates.Contains(E) && Emap.Add(E))
      {
        if (BRep_Tool::Degenerated(E))
        {
          if (hasPCurve(E))
          {
            DegenEdges.Append(E);
          }
          continue;
        }
        StoreInMVE(myFace, E, MVE, YaCouture, myVerticesForSubstitute, myTolConf);
      }
    }
  }

  // add const edges
  // Sewn edges can be doubled or not in myConstEdges
  // => call only once StoreInMVE which should double them
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> DejaVu;
  for (itl.Initialize(myConstEdges); itl.More(); itl.Next())
  {
    TopoDS_Edge& E = TopoDS::Edge(itl.ChangeValue());
    if (!aDuplicates.Contains(E) && DejaVu.Add(E))
    {
      if (BRep_Tool::Degenerated(E))
      {
        if (hasPCurve(E))
        {
          DegenEdges.Append(E);
        }
        continue;
      }
      StoreInMVE(myFace, E, MVE, YaCouture, myVerticesForSubstitute, myTolConf);
    }
  }

  UpdateVEmap(MVE);
  if (MVE.IsEmpty())
  {
    return;
  }

  NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher> CurrentVEMap;
  NCollection_List<TopoDS_Shape>                                           CurrentEdgeList;
  MapOfWire                                                                NewWires;
  TopoDS_Vertex                                                            V1, V2, NV, NNV;

  //-----------------------------------------------
  // Find all possible closed wires
  //----------------------------------------------

  MapIteratorOfMapOfWire    itW;
  TopLoc_Location           L;
  occ::handle<Geom_Surface> S          = BRep_Tool::Surface(myFace, L);
  bool                      IsPeriodic = S->IsUPeriodic() || S->IsVPeriodic();
  if (IsPeriodic && DegenEdges.IsEmpty() && IsInOnePeriod(MVE, myFace, S))
  {
    IsPeriodic = false;
    SHOW_TOPO_SHAPE(myFace, "LoopInOnePeriod");
  }
  DejaVu.Clear();

  // By angle where the face is not periodic, or its edges lie within one
  // period; the search otherwise, or when the walk cannot describe the network.
  const int aWalkMode = IsPeriodic ? 0 : LoopWalkMode();
  MapOfWire aWalked;
  const bool isWalked =
    aWalkMode != 0 && FindLoopsByAngle(MVE, myFace, aWalked, myOutsideEdges);
  if (isWalked && aWalkMode == 1)
  {
    NewWires = aWalked;
  }

  for (int ii = 1; ii <= MVE.Extent() && !(isWalked && aWalkMode == 1); ++ii)
  {
    const TopoDS_Vertex& VF = TopoDS::Vertex(MVE.FindKey(ii));

    for (itl.Initialize(MVE(ii)); itl.More(); itl.Next())
    {
      TopoDS_Edge CE = TopoDS::Edge(itl.Value());
      TopExp::Vertices(CE, V1, V2, true);

      if (!DejaVu.Add(CE))
      {
        continue;
      }

      FindAllLoops(VF,
                   CE,
                   CurrentVEMap,
                   CurrentEdgeList,
                   MVE,
                   NewWires,
                   myFace,
                   myTolConf,
                   myOutsideEdges);

      // Perioidc surface needs a wire with seam edge. Look for wires consists
      // of a wire with two closed edge joined by a seam edge.
      //
      // Note: sphere have one and only one edge (actually two identical edge
      // but opposite orientation) that is the seam. But we can't have a
      // continuous sphere face when dealing with thick solid. Or can we???

      // Start with closed periodic edge first.
      if (!IsPeriodic || !V1.IsSame(V2))
      {
        continue;
      }

      NV = V2;
      SHOW_TOPO_SHAPE(CE, "SeamEdge1", true);
      SHOW_TOPO_SHAPE(NV, "SeamEdgeV1", true);
      for (itl1.Initialize(MVE.FindFromKey(NV)); itl1.More(); itl1.Next())
      {
        TopoDS_Edge NE = TopoDS::Edge(itl1.Value());
        if (NE.IsSame(CE))
        {
          continue;
        }
        TopExp::Vertices(NE, V1, V2, true);
        if (V1.IsSame(V2))
        {
          continue;
        }

        if (NV.IsSame(V1))
        {
          NNV = V2;
        }
        else
        {
          NNV = V1;
          V1  = V2;
          V2  = NNV;
          NE.Reverse();
        }

        // The far end of the seam piece can be a pole: the band is then closed
        // by the pole's degenerated edge, which is kept out of MVE -- a dome
        // whose flat face is removed, its offset bounded by the seam, the
        // pole and one circle. Without it no seam wire is built, and the
        // circle is left as a wire on its own, bounding nothing.
        NCollection_List<TopoDS_Shape> aFarEdges = MVE.FindFromKey(NNV);
        for (itl2.Initialize(DegenEdges); itl2.More(); itl2.Next())
        {
          if (TopExp::FirstVertex(TopoDS::Edge(itl2.Value().Oriented(TopAbs_FORWARD))).IsSame(NNV))
          {
            aFarEdges.Append(itl2.Value());
          }
        }
        for (itl2.Initialize(aFarEdges); itl2.More(); itl2.Next())
        {
          TopoDS_Edge NNE = TopoDS::Edge(itl2.Value());
          if (NNE.IsSame(CE))
          {
            continue;
          }
          TopExp::Vertices(NNE, V1, V2, true);
          if (V1.IsSame(V2))
          {
            // A face can have more than one seam wire -- a removed cylinder
            // leaves a wall at each end, each closed by its own piece of the
            // seam -- but not two sharing an edge: a seam wire already built
            // on one of these edges is this one, or one it cannot sit beside.
            // Except that a band on the seam's own span beats one on a piece
            // beyond it (myOutsideEdges): the holed cone's inner offsets cross
            // below a removed top, and the band from the top circle down to
            // the crossing, found first, took the crossing circle from the
            // band below it -- the cavity -- which was never built.
            const bool isOutside = myOutsideEdges.Contains(NE);
            bool       isTaken   = false;
            for (itW = MapIteratorOfMapOfWire(NewWires); itW.More() && !isTaken; itW.Next())
            {
              const WireInfo& info = itW.Value();
              isTaken = info.HasSeam
                        && (info.Contains(NE)
                            || ((info.Contains(CE) || info.Contains(NNE))
                                && (isOutside || !info.Outside)));
            }
            if (isTaken)
            {
              break;
            }

            SHOW_TOPO_SHAPE(NE, "SeamEdge2", true);
            SHOW_TOPO_SHAPE(NNV, "SeamEdgeV2", true);
            SHOW_TOPO_SHAPE(NNE, "SeamEdge3", true);

            TopoDS_Wire NW = MakeSeamBand(NE, CE, NNE, myFace);
            if (!NW.IsNull())
            {
              SHOW_TOPO_SHAPE(NW, "SeamBand");
            }
            else
            {
              BRepLib_MakeWire aMakeWire;
              aMakeWire.Add(TopoDS::Edge(NE.Reversed()));
              aMakeWire.Add(CE);
              aMakeWire.Add(NE);
              aMakeWire.Add(NNE);
              if (!aMakeWire.IsDone())
              {
                SHOW_TOPO_SHAPE(CE, "DiscardSeamChain");
                continue;
              }
              NW = aMakeWire.Wire();

              ShapeFix_Wire aFixer;
              aFixer.Load(NW);
              aFixer.SetFace(myFace);
              aFixer.FixReorder();
              aFixer.FixConnected();
              aFixer.FixSeam(0);
              aFixer.FixEdgeCurves();
              aFixer.FixDegenerated();
              NW = aFixer.Wire();
            }

            WireInfo aSeamInfo(NW, true);
            aSeamInfo.Outside = isOutside;
            if (!NW.Closed() || NewWires.Contains(aSeamInfo))
            {
              SHOW_TOPO_SHAPE(NW, "DiscardWire2");
              continue;
            }

            // The seam wire replaces the plain wires made of its edges, and
            // the band beyond the span it beat -- only now that it exists, or
            // a failed build would lose them for nothing.
            for (int iw = NewWires.Extent(); iw >= 1; --iw)
            {
              const WireInfo& info = NewWires(iw);
              if (info.Contains(CE) || info.Contains(NE) || info.Contains(NNE))
              {
                SHOW_TOPO_SHAPE(info.aWire, info.HasSeam ? "SeamOutsideRemove" : "SeamRemove");
                NewWires.RemoveFromIndex(iw);
              }
            }
            DejaVu.Add(NE);
            DejaVu.Add(NNE);
            NewWires.Add(aSeamInfo);
            SHOW_TOPO_SHAPE(NW, "NewWire2");
          }
        }
      }
    }
  }

  DropOutsideWires(NewWires, myFace);

  if (!IsPeriodic)
  {
    //-----------------------------------------------
    // Split wires
    //----------------------------------------------
    SplitWires(myNewWires, NewWires, MVE, myFace);
    if (isWalked && aWalkMode == 2)
    {
      // Both ways, and the minimal wires compared.
      DropOutsideWires(aWalked, myFace);
      NCollection_List<TopoDS_Shape> aWalkedOut;
      SplitWires(aWalkedOut, aWalked, MVE, myFace);
      MapOfWire aA, aB;
      for (itl.Initialize(myNewWires); itl.More(); itl.Next())
      {
        aA.Add(WireInfo(TopoDS::Wire(itl.Value())));
      }
      for (itl.Initialize(aWalkedOut); itl.More(); itl.Next())
      {
        aB.Add(WireInfo(TopoDS::Wire(itl.Value())));
      }
      bool isSame = aA.Extent() == aB.Extent();
      for (MapIteratorOfMapOfWire itA(aA); itA.More() && isSame; itA.Next())
      {
        isSame = aB.Contains(itA.Value());
      }
      if (!isSame)
      {
        std::fprintf(stderr, "LOOPWALK MISMATCH search %d walk %d\n", aA.Extent(), aB.Extent());
        SHOW_TOPO_SHAPE(myFace, "LoopWalkMismatchSearch", myNewWires);
        SHOW_TOPO_SHAPE(myFace, "LoopWalkMismatchWalk", aWalkedOut);
      }
      else
      {
        std::fprintf(stderr, "LOOPWALK SAME %d\n", aA.Extent());
      }
    }
  }
  else
  {
    // The seam wire replaces the closed-edge wires found before it, but the
    // search goes on from the other vertices and can find one of its closed
    // edges again on its own afterwards -- the floor circle of a cylinder
    // whose seam wire was built from the top circle. An edge is in one wire
    // of a face (a seam twice, in one wire), so a wire made only of the seam
    // wire's edges is that edge found again, and the face must not get it.
    NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aSeamWireEdges;
    for (itW = MapIteratorOfMapOfWire(NewWires); itW.More(); itW.Next())
    {
      if (itW.Value().HasSeam)
      {
        for (const auto& v : itW.Value().Edges)
        {
          aSeamWireEdges.Add(v.first);
        }
      }
    }
    // The search returns every closed wire, and on a periodic face nothing
    // sorts them afterwards. Where the face's own edges run on past the cut
    // -- an offset sphere trimmed short of its pole -- it finds the wire
    // round the face and, through the pieces beyond, the wire round what was
    // cut away, run the other way, and one of no area along the pieces
    // themselves. An open edge is in one wire of a face: a wire bounding
    // nothing (no area, or a negative one) loses to a wire that shares an
    // open edge with it and bounds something.
    NCollection_Map<int> aLosers;
    {
      std::vector<double> anAreas(NewWires.Extent() + 1, 0.);
      bool                hasBad = false;
      for (int iw = 1; iw <= NewWires.Extent(); ++iw)
      {
        const WireInfo& anInfo = NewWires(iw);
        if (anInfo.HasSeam || !anInfo.aWire.Closed())
        {
          anAreas[iw] = 1.;
          continue;
        }
        try
        {
          OCC_CATCH_SIGNALS
          BRep_Builder aTB;
          TopoDS_Face  aTF = TopoDS::Face(myFace.EmptyCopied().Oriented(TopAbs_FORWARD));
          // A wire of its own: one added to a face is free no longer, and
          // the pole's edge is still to be put into the real one.
          TopoDS_Wire aTW;
          aTB.MakeWire(aTW);
          for (TopoDS_Iterator aEIt(anInfo.aWire); aEIt.More(); aEIt.Next())
          {
            aTB.Add(aTW, aEIt.Value());
          }
          aTW.Closed(true);
          aTB.Add(aTF, aTW);
          GProp_GProps aGP;
          BRepGProp::SurfaceProperties(aTF, aGP);
          anAreas[iw] = aGP.Mass();
        }
        catch (Standard_Failure const&)
        {
          anAreas[iw] = 1.;
        }
        hasBad = hasBad || anAreas[iw] <= Precision::SquareConfusion();
      }
      for (int iw = 1; hasBad && iw <= NewWires.Extent(); ++iw)
      {
        if (anAreas[iw] > Precision::SquareConfusion())
        {
          continue;
        }
        const WireInfo& anInfo = NewWires(iw);
        for (int jw = 1; jw <= NewWires.Extent() && !aLosers.Contains(iw); ++jw)
        {
          if (jw == iw || anAreas[jw] <= Precision::SquareConfusion() || NewWires(jw).HasSeam)
          {
            continue;
          }
          for (const auto& v : anInfo.Edges)
          {
            TopoDS_Vertex aV1, aV2;
            TopExp::Vertices(TopoDS::Edge(v.first), aV1, aV2);
            if (!aV1.IsNull() && !aV1.IsSame(aV2)
                && !BRep_Tool::IsClosed(TopoDS::Edge(v.first), myFace)
                && NewWires(jw).Contains(TopoDS::Edge(v.first)))
            {
              aLosers.Add(iw);
              SHOW_TOPO_SHAPE(anInfo.aWire, "WireBoundsNothing");
              break;
            }
          }
        }
      }
    }
    int aWireIndex = 0;
    for (itW = MapIteratorOfMapOfWire(NewWires); itW.More(); itW.Next())
    {
      const WireInfo& anInfo = itW.Value();
      if (aLosers.Contains(++aWireIndex))
      {
        continue;
      }
      if (!anInfo.HasSeam && !aSeamWireEdges.IsEmpty())
      {
        bool isInSeamWire = true;
        for (const auto& v : anInfo.Edges)
        {
          if (!aSeamWireEdges.Contains(v.first))
          {
            isInSeamWire = false;
            break;
          }
        }
        if (isInSeamWire)
        {
          SHOW_TOPO_SHAPE(anInfo.aWire, "SeamWireAgain");
          continue;
        }
        // A wire running once round the period does not close in UV: on its
        // own it bounds nothing. It closes a face only with the seam, and the
        // seam wires are all built -- the circle where a vanishing face's
        // offset cut this one, beyond the seam's span.
        if (!IsClosedInUV(anInfo.aWire, myFace))
        {
          SHOW_TOPO_SHAPE(anInfo.aWire, "WrapWireAlone");
          continue;
        }
      }
      myNewWires.Append(anInfo.aWire);
    }
  }

  // Re-insert the degenerated edges: each one belongs to the wire that passes
  // through the singularity vertex - without it that wire does not close in UV
  // space and the resulting face is unorientable.
  if (!DegenEdges.IsEmpty())
  {
    BRep_Builder    aBB;
    TopExp_Explorer aVExp;
    for (itl.Initialize(DegenEdges); itl.More(); itl.Next())
    {
      const TopoDS_Edge&  DE  = TopoDS::Edge(itl.Value());
      const TopoDS_Vertex aDV = TopExp::FirstVertex(TopoDS::Edge(DE.Oriented(TopAbs_FORWARD)));
      bool                added = false;
      // Already in a seam band closed at the pole.
      for (itl1.Initialize(myNewWires); itl1.More() && !added; itl1.Next())
      {
        for (TopoDS_Iterator aEIt(itl1.Value()); aEIt.More() && !added; aEIt.Next())
        {
          added = aEIt.Value().IsSame(DE);
        }
      }
      for (itl1.Initialize(myNewWires); itl1.More() && !added; itl1.Next())
      {
        TopoDS_Wire& W = TopoDS::Wire(itl1.ChangeValue());
        for (aVExp.Init(W, TopAbs_VERTEX); aVExp.More(); aVExp.Next())
        {
          if (aVExp.Current().IsSame(aDV))
          {
            // A wire the minimal-wire search tried in a face is free no longer.
            W.Free(true);
            aBB.Add(W, DE);
            SHOW_TOPO_SHAPE(DE, "DegenEdgeKept");
            added = true;
            break;
          }
        }
      }
      // The wire can reach the pole on a vertex of its own: two meridians of
      // a face cut short of a full turn, each ending on the vertex its own
      // section gave it. The pole's edge is kept out of the vertex map, so
      // it never took the vertex the map settled on; it takes it here.
      if (!added && !aDV.IsNull())
      {
        const gp_Pnt aDP = BRep_Tool::Pnt(aDV);
        for (itl1.Initialize(myNewWires); itl1.More() && !added; itl1.Next())
        {
          TopoDS_Wire& W = TopoDS::Wire(itl1.ChangeValue());
          for (aVExp.Init(W, TopAbs_VERTEX); aVExp.More(); aVExp.Next())
          {
            const TopoDS_Vertex& aWV = TopoDS::Vertex(aVExp.Current());
            if (aDP.Distance(BRep_Tool::Pnt(aWV))
                > std::max(myTolConf,
                           std::max(BRep_Tool::Tolerance(aWV), BRep_Tool::Tolerance(aDV))))
            {
              continue;
            }
            TopoDS_Edge                    aDE = TopoDS::Edge(DE.Oriented(TopAbs_FORWARD));
            NCollection_List<TopoDS_Shape> aOldV;
            for (TopoDS_Iterator aVIt(aDE); aVIt.More(); aVIt.Next())
            {
              aOldV.Append(aVIt.Value());
            }
            aDE.Free(true);
            for (NCollection_List<TopoDS_Shape>::Iterator aOIt(aOldV); aOIt.More(); aOIt.Next())
            {
              aBB.Remove(aDE, aOIt.Value());
              aBB.Add(aDE, aWV.Oriented(aOIt.Value().Orientation()));
            }
            if (!myVerticesForSubstitute.IsBound(aDV))
            {
              myVerticesForSubstitute.Bind(aDV, aWV);
            }
            W.Free(true);
            aBB.Add(W, DE);
            SHOW_TOPO_SHAPE(DE, "DegenEdgeKeptAtPoint");
            added = true;
            break;
          }
        }
      }
      if (!added)
      {
        SHOW_TOPO_SHAPE(DE, "DegenEdgeDropped");
      }
    }
  }

  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> UsedEdges;
  TopExp_Explorer                                        aExp;
  for (itl.Initialize(myNewWires); itl.More(); itl.Next())
  {
    for (aExp.Init(itl.Value(), TopAbs_EDGE); aExp.More(); aExp.Next())
    {
      UsedEdges.Add(aExp.Current());
    }
  }

  PurgeNewEdges(myCutEdges, UsedEdges, myKeptEdges);
}

//=================================================================================================

void BRepAlgo_Loop::CutEdge(const TopoDS_Edge&                    E,
                            const NCollection_List<TopoDS_Shape>& VOnE,
                            NCollection_List<TopoDS_Shape>&       NE) const
{
  CutEdge(E, VOnE, NE, true);
}

void BRepAlgo_Loop::CutEdge(const TopoDS_Edge&                    E,
                            const NCollection_List<TopoDS_Shape>& VOnE,
                            NCollection_List<TopoDS_Shape>&       NE,
                            bool                                  KeepAll) const
{
  double       Tol     = 0.001; // 5.e-05; //5.e-07;
  TopoDS_Shape aLocalE = E.Oriented(TopAbs_FORWARD);
  TopoDS_Edge  WE      = TopoDS::Edge(aLocalE);

  double                                   U1, U2;
  TopoDS_Vertex                            V1, V2;
  NCollection_Sequence<TopoDS_Shape>       SV;
  NCollection_Sequence<double>             SU;
  NCollection_List<TopoDS_Shape>::Iterator it(VOnE);
  BRep_Builder                             B;

  for (; it.More(); it.Next())
  {
    SV.Append(it.Value());
  }
  //--------------------------------
  // Parse vertices on the edge.
  //--------------------------------
  Bubble(WE, SV, SU);

  // KeepAll = true;

  if (KeepAll)
  {
    SHOW_TOPO_SHAPE(WE, "CuttingInter");
  }
  else
  {
    SHOW_TOPO_SHAPE(WE, "Cutting");
  }

  int NbVer = SV.Length();
  //----------------------------------------------------------------
  // Construction of new edges.
  // Note :  vertices at the extremities of edges are not
  //         onligatorily in the list of vertices
  //----------------------------------------------------------------
  if (SV.IsEmpty())
  {
    NE.Append(E);
    return;
  }

  for (int ii = 1; ii <= SV.Length(); ++ii)
  {
    SHOW_TOPO_SHAPE(SV(ii), "CuttingV", 1);
  }

  TopoDS_Vertex             VF, VL;
  double                    f, l;
  occ::handle<Geom2d_Curve> C = BRep_Tool::CurveOnSurface(WE, myFace, f, l);
  TopExp::Vertices(WE, VF, VL);

  TopoDS_Iterator It(WE);
  bool            Extended = false;
  for (; It.More(); It.Next())
  {
    if (It.Value().Orientation() == TopAbs_INTERNAL)
    {
      SHOW_TOPO_SHAPE(WE, "IEdge");
      SHOW_TOPO_SHAPE(It.Value(), "Internal");
      Extended = true;
      break;
    }
  }

  if (NbVer == 2)
  {
    if (SV(1).IsEqual(VF) && SV(2).IsEqual(VL))
    {
      NE.Append(E);
      return;
    }
  }

  auto InsertVertex = [&](const TopoDS_Shape& V, double U) {
    for (int ii = 1; ii <= SV.Length(); ++ii)
    {
      if (SV(ii).IsSame(V) && SU(ii) == U)
      {
        return;
      }
      if (SU(ii) > U)
      {
        SHOW_TOPO_SHAPE(V, "CuttingInsV1", 1);
        SU.InsertBefore(ii, U);
        SV.InsertBefore(ii, V);
        return;
      }
    }
    SHOW_TOPO_SHAPE(V, "CuttingInsV2", 1);
    SU.Append(U);
    SV.Append(V);
  };

  //----------------------------------------------------
  // Processing of closed edges
  // If a vertex of intersection is on the common vertex
  // it should appear at the beginning and end of SV.
  //----------------------------------------------------
  TopoDS_Vertex VCEI;
  if (!VF.IsNull() && VF.IsSame(VL))
  {
    VCEI = UpdateClosedEdge(WE, SV, SU);
    if (!VCEI.IsNull())
    {
      TopoDS_Shape aLocalV = VCEI.Oriented(TopAbs_FORWARD);
      VF                   = TopoDS::Vertex(aLocalV);
      aLocalV              = VCEI.Oriented(TopAbs_REVERSED);
      VL                   = TopoDS::Vertex(aLocalV);
      InsertVertex(VF, f);
      InsertVertex(VL, l);
    }
    else if (!Extended)
    {
      InsertVertex(VF, f);
      InsertVertex(VL, l);
    }
  }
  else if (!Extended)
  {
    //-----------------------------------------
    // Eventually all extremities of the edge.
    //-----------------------------------------
    if (!VF.IsNull())
    {
      InsertVertex(VF, f);
    }
    if (!VL.IsNull())
    {
      InsertVertex(VL, l);
    }
  }

  while (!SV.IsEmpty())
  {
    while (!KeepAll && !SV.IsEmpty() && SV.First().Orientation() != TopAbs_FORWARD)
    {
      SHOW_TOPO_SHAPE(SV.First(), SV.First().Orientation() == TopAbs_FORWARD ? "VF1" : "VR1");
      SV.Remove(1);
    }
    if (SV.IsEmpty())
    {
      break;
    }
    V1 = TopoDS::Vertex(SV.First());
    SHOW_TOPO_SHAPE(SV.First(), SV.First().Orientation() == TopAbs_FORWARD ? "VF2" : "VR2");

    SV.Remove(1);
    if (SV.IsEmpty())
    {
      break;
    }
    SHOW_TOPO_SHAPE(SV.First(), SV.First().Orientation() == TopAbs_FORWARD ? "VF4" : "VR4");
    if (KeepAll || SV.First().Orientation() == TopAbs_REVERSED)
    {
      V2 = TopoDS::Vertex(SV.First());
      //-------------------------------------------
      // Copy the edge and restriction by V1 V2.
      //-------------------------------------------
      TopoDS_Shape NewEdge    = WE.EmptyCopied();
      TopoDS_Shape aLocalEdge = V1.Oriented(TopAbs_FORWARD);
      B.Add(NewEdge, aLocalEdge);
      aLocalEdge = V2.Oriented(TopAbs_REVERSED);
      B.Add(TopoDS::Edge(NewEdge), aLocalEdge);
      if (V1.IsSame(VF))
      {
        U1 = f;
      }
      else
      {
        TopoDS_Shape aLocalV = V1.Oriented(TopAbs_INTERNAL);
        U1                   = BRep_Tool::Parameter(TopoDS::Vertex(aLocalV), WE);
      }
      if (V2.IsSame(VL))
      {
        U2 = l;
      }
      else
      {
        TopoDS_Shape aLocalV = V2.Oriented(TopAbs_INTERNAL);
        U2                   = BRep_Tool::Parameter(TopoDS::Vertex(aLocalV), WE);
      }
      B.Range(TopoDS::Edge(NewEdge), U1, U2);
      SHOW_TOPO_SHAPE(NewEdge, "CutEdge");
      NE.Append(NewEdge.Oriented(E.Orientation()));
    }
  }

  // Remove edges with size <= tolerance
  it.Initialize(NE);
  while (it.More())
  {
    // skl : I change "E" to "EE"
    TopoDS_Edge EE = TopoDS::Edge(it.Value());
    double      fpar, lpar;
    BRep_Tool::Range(EE, fpar, lpar);
    if (lpar - fpar <= Precision::Confusion())
    {
      SHOW_TOPO_SHAPE(EE, "CutEdgeRemove1");
      NE.Remove(it);
    }
    else
    {
      gp_Pnt2d pf, pl;
      BRep_Tool::UVPoints(EE, myFace, pf, pl);
      if (pf.Distance(pl) <= Tol && !BRep_Tool::IsClosed(EE))
      {
        SHOW_TOPO_SHAPE(EE, "CutEdgeRemove2");
        NE.Remove(it);
      }
      else
      {
        it.Next();
      }
    }
  }
}

//=================================================================================================

const NCollection_List<TopoDS_Shape>& BRepAlgo_Loop::NewWires() const
{
  return myNewWires;
}

//=================================================================================================

const NCollection_List<TopoDS_Shape>& BRepAlgo_Loop::NewFaces() const
{
  return myNewFaces;
}

//=================================================================================================

void BRepAlgo_Loop::WiresToFaces()
{
  if (!myNewWires.IsEmpty())
  {
    SHOW_TOPO_SHAPE(myFace, "FaceRestrict", myNewWires);
    BRepAlgo_FaceRestrictor FR;
    TopoDS_Shape            aLocalS = myFace.Oriented(TopAbs_FORWARD);
    FR.Init(TopoDS::Face(aLocalS), false, true);
    //    FR.Init (TopoDS::Face(myFace.Oriented(TopAbs_FORWARD)),
    //	     false);
    NCollection_List<TopoDS_Shape>::Iterator it(myNewWires);
    for (; it.More(); it.Next())
    {
      FR.Add(TopoDS::Wire(it.ChangeValue()));
    }

    FR.Perform();

    if (FR.IsDone())
    {
      TopAbs_Orientation OriF = myFace.Orientation();
      for (; FR.More(); FR.Next())
      {
        SHOW_TOPO_SHAPE(FR.Current(), "NewFace");
        myNewFaces.Append(FR.Current().Oriented(OriF));
      }
    }
  }
}

//=================================================================================================

const NCollection_List<TopoDS_Shape>& BRepAlgo_Loop::NewEdges(const TopoDS_Edge& E) const
{
  return myCutEdges.FindFromKey(E);
}

//=================================================================================================

void BRepAlgo_Loop::GetVerticesForSubstitute(
  NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>& VerVerMap) const
{
  VerVerMap = myVerticesForSubstitute;
}

//=================================================================================================

void BRepAlgo_Loop::VerticesForSubstitute(
  NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>& VerVerMap)
{
  myVerticesForSubstitute = VerVerMap;
}

//=================================================================================================

void BRepAlgo_Loop::UpdateVEmap(
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>&
    theVEmap)
{
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    VerLver;

  for (int ii = 1; ii <= theVEmap.Extent(); ii++)
  {
    const TopoDS_Vertex&                  aVertex = TopoDS::Vertex(theVEmap.FindKey(ii));
    const NCollection_List<TopoDS_Shape>& aElist  = theVEmap(ii);
    if (aElist.Extent() == 1 && myImageVV.IsImage(aVertex))
    {
      const TopoDS_Vertex& aProVertex = TopoDS::Vertex(myImageVV.ImageFrom(aVertex));
      SHOW_TOPO_SHAPE(aElist.First(), "VEMapE");
      SHOW_TOPO_SHAPE(aVertex, "VEMapV");
      SHOW_TOPO_SHAPE(aProVertex, "VEMapVI");
      if (VerLver.Contains(aProVertex))
      {
        NCollection_List<TopoDS_Shape>& aVlist = VerLver.ChangeFromKey(aProVertex);
        aVlist.Append(aVertex.Oriented(TopAbs_FORWARD));
      }
      else
      {
        NCollection_List<TopoDS_Shape> aVlist;
        aVlist.Append(aVertex.Oriented(TopAbs_FORWARD));
        VerLver.Add(aProVertex, aVlist);
      }
    }
  }

  if (VerLver.IsEmpty())
  {
    return;
  }

  BRep_Builder aBB;
  for (int ii = 1; ii <= VerLver.Extent(); ii++)
  {
    // In some cases (concave faces), a vertex may have images in more than one
    // distinct location. So we can't always assume all images to converge into
    // one vertex.
    const NCollection_List<TopoDS_Shape>& aVertexList = VerLver(ii);
    if (aVertexList.Extent() == 1)
    {
      continue;
    }

    SHOW_TOPO_SHAPE(TopoDS_Shape(), "VEMap", aVertexList);

    NCollection_List<TopoDS_Shape> aVlist = aVertexList;
    NCollection_List<TopoDS_Shape> OutLiers;

    gp_Pnt                                   aCentre;
    double                                   Tol2 = myTolConf * myTolConf;
    NCollection_Array1<gp_Pnt>               Points(1, aVlist.Extent());
    NCollection_List<TopoDS_Shape>::Iterator itl;

    auto GetCenter = [&](const NCollection_List<TopoDS_Shape>& Vertices) {
      double aMaxTol = 0.;

      int Count = 0;
      for (itl.Initialize(Vertices); itl.More();)
      {
        const TopoDS_Vertex& aVertex = TopoDS::Vertex(itl.Value());
        double               aTol    = BRep_Tool::Tolerance(aVertex);
        aMaxTol                      = std::max(aMaxTol, aTol);
        gp_Pnt aPnt                  = BRep_Tool::Pnt(aVertex);
        if (Count == 0 || Points(1).SquareDistance(aPnt) < Tol2)
        {
          Points(++Count) = aPnt;
          itl.Next();
        }
        else
        {
          OutLiers.Append(aVertex);
          aVlist.Remove(itl);
        }
      }

      gp_Ax2 anAxis;
      bool   IsSingular;
      // Only the first Count slots were filled on this pass; the tail holds
      // zeros or points of a previously processed cluster.
      NCollection_Array1<gp_Pnt> aFilledPoints(1, Count);
      for (int jj = 1; jj <= Count; jj++)
      {
        aFilledPoints(jj) = Points(jj);
      }
      GeomLib::AxeOfInertia(aFilledPoints, anAxis, IsSingular);
      aCentre         = anAxis.Location();
      double aMaxDist = 0.;
      for (int jj = 1; jj <= Count; jj++)
      {
        double aSqDist = aCentre.SquareDistance(Points(jj));
        aMaxDist       = std::max(aMaxDist, aSqDist);
      }
      aMaxDist = std::sqrt(aMaxDist);
      return std::max(aMaxTol, aMaxDist);
    };

    while (true)
    {
      double aMaxTol = GetCenter(aVlist);

      // Find constant vertex
      TopoDS_Vertex aConstVertex;
      for (itl.Initialize(aVlist); itl.More(); itl.Next())
      {
        const TopoDS_Vertex&                     aVertex = TopoDS::Vertex(itl.Value());
        const NCollection_List<TopoDS_Shape>&    aElist  = theVEmap.FindFromKey(aVertex);
        const TopoDS_Shape&                      anEdge  = aElist.First();
        NCollection_List<TopoDS_Shape>::Iterator itcedges(myConstEdges);
        for (; itcedges.More(); itcedges.Next())
        {
          if (anEdge.IsSame(itcedges.Value()))
          {
            aConstVertex = aVertex;
            SHOW_TOPO_SHAPE(aConstVertex, "ConstV_");
            break;
          }
        }
        if (!aConstVertex.IsNull())
        {
          break;
        }
      }
      if (aConstVertex.IsNull())
      {
        aConstVertex = TopoDS::Vertex(aVlist.First());
        SHOW_TOPO_SHAPE(aConstVertex, "ConstVN");
      }
      aBB.UpdateVertex(aConstVertex, aCentre, aMaxTol);

      for (itl.Initialize(aVlist); itl.More(); itl.Next())
      {
        const TopoDS_Vertex& aVertex = TopoDS::Vertex(itl.Value());
        if (aVertex.IsSame(aConstVertex))
        {
          continue;
        }

        const NCollection_List<TopoDS_Shape>& aElist = theVEmap.FindFromKey(aVertex);
        for (NCollection_List<TopoDS_Shape>::Iterator itl1(aElist); itl1.More(); itl1.Next())
        {
          TopoDS_Edge anEdge = TopoDS::Edge(itl1.Value());
          SHOW_TOPO_SHAPE(anEdge, "ReplaceConstV");
          anEdge.Orientation(TopAbs_FORWARD);
          TopoDS_Vertex aV1, aV2;
          TopExp::Vertices(anEdge, aV1, aV2);
          TopoDS_Vertex aVertexToRemove = (aV1.IsSame(aVertex)) ? aV1 : aV2;
          anEdge.Free(true);
          aBB.Remove(anEdge, aVertexToRemove);
          aBB.Add(anEdge, aConstVertex.Oriented(aVertexToRemove.Orientation()));
        }
      }

      if (OutLiers.Extent() <= 1)
      {
        break;
      }
      aVlist.Clear();
      aVlist.Append(OutLiers);
    }
  }

  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> Emap;
  for (int ii = 1; ii <= theVEmap.Extent(); ii++)
  {
    const NCollection_List<TopoDS_Shape>&    aElist = theVEmap(ii);
    NCollection_List<TopoDS_Shape>::Iterator itl(aElist);
    for (; itl.More(); itl.Next())
    {
      Emap.Add(itl.Value());
    }
  }

  theVEmap.Clear();
  for (int ii = 1; ii <= Emap.Extent(); ii++)
  {
    TopExp::MapShapesAndAncestors(Emap(ii), TopAbs_VERTEX, TopAbs_EDGE, theVEmap);
  }
}
