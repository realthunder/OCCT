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
#include <BRepCheck_Analyzer.hxx>
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
  const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& UsedEdges)
{
  for (int ii = 1; ii <= NewEdges.Extent(); ++ii)
  {
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
  if (_CollectingEdges)
  {
    NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> EMap;
    NCollection_List<TopoDS_Shape>                         theVerts;
    NCollection_List<TopoDS_Shape>                         LV;
    NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> MV;

    for (itl.Initialize(myConstEdges); itl.More(); itl.Next())
    {
      theEdges.Append(itl.Value());
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
        TopExp::Vertices(otherEdge, OV1, OV2);
        theVerts.Append(OV1);
        theVerts.Append(OV2);
        for (itl2.Initialize(theVerts); itl2.More(); itl2.Next())
        {
          TopoDS_Vertex aVertex = TopoDS::Vertex(itl2.Value());
          if (!MV.Add(aVertex))
          {
            continue;
          }
          double Tol = BRep_Tool::Tolerance(aVertex);
          gp_Pnt OP  = BRep_Tool::Pnt(aVertex);
          if (OP.Distance(BRep_Tool::Pnt(V1)) < Tol || OP.Distance(BRep_Tool::Pnt(V2)) < Tol)
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
                else if (P < (FP + LP) / 2)
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
        ConstEdges.Append(anEdge);
        SHOW_TOPO_SHAPE(anEdge, "ConstEdge");
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
        HasSeam(other.HasSeam)
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
                  const double&                                                         Tol)
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
      if (NW.Closed() && !NewWires.Contains(NW) && NewWires.Add(NW) > 0)
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
        FindAllLoops(NV, NE, CurrentVEMap, CurrentEdgeList, MVE, NewWires, aFace, Tol);
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
        ShapeFix_Shape aFix(aNF);
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
      if (Emap.Add(E))
      {
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
    if (DejaVu.Add(E))
    {
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
  DejaVu.Clear();

  for (int ii = 1; ii <= MVE.Extent(); ++ii)
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

      FindAllLoops(VF, CE, CurrentVEMap, CurrentEdgeList, MVE, NewWires, myFace, myTolConf);

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

        for (itl2.Initialize(MVE.FindFromKey(NNV)); itl2.More(); itl2.Next())
        {
          TopoDS_Edge NNE = TopoDS::Edge(itl2.Value());
          if (NNE.IsSame(CE))
          {
            continue;
          }
          TopExp::Vertices(NNE, V1, V2, true);
          if (V1.IsSame(V2))
          {
            // Look for an already-built seam wire first: only when none exists
            // may the plain wires be removed and replaced. Removing before the
            // scan completes would destroy wires without a replacement whenever
            // a seam wire happens to come later in map-iteration order.
            bool hasSeamWire = false;
            for (itW = MapIteratorOfMapOfWire(NewWires); itW.More(); itW.Next())
            {
              if (itW.Value().HasSeam)
              {
                hasSeamWire = true;
                break;
              }
            }
            if (hasSeamWire)
            {
              break;
            }
            for (int iw = NewWires.Extent(); iw >= 1; --iw)
            {
              const WireInfo& info = NewWires(iw);
              if (info.Contains(CE) || info.Contains(NE) || info.Contains(NNE))
              {
                SHOW_TOPO_SHAPE(info.aWire, "SeamRemove");
                NewWires.RemoveFromIndex(iw);
              }
            }

            DejaVu.Add(NE);
            DejaVu.Add(NNE);

            SHOW_TOPO_SHAPE(NE, "SeamEdge2", true);
            SHOW_TOPO_SHAPE(NNV, "SeamEdgeV2", true);
            SHOW_TOPO_SHAPE(NNE, "SeamEdge3", true);

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
            TopoDS_Wire NW = aMakeWire.Wire();

            ShapeFix_Wire aFixer;
            aFixer.Load(NW);
            aFixer.SetFace(myFace);
            aFixer.FixReorder();
            aFixer.FixConnected();
            aFixer.FixSeam(0);
            aFixer.FixEdgeCurves();
            aFixer.FixDegenerated();
            NW = aFixer.Wire();

            WireInfo aSeamInfo(NW, true);
            if (NW.Closed() && !NewWires.Contains(aSeamInfo) && NewWires.Add(aSeamInfo) > 0)
            {
              SHOW_TOPO_SHAPE(NW, "NewWire2");
            }
            else
            {
              SHOW_TOPO_SHAPE(NW, "DiscardWire2");
            }
          }
        }
      }
    }
  }

  if (!IsPeriodic)
  {
    //-----------------------------------------------
    // Split wires
    //----------------------------------------------
    SplitWires(myNewWires, NewWires, MVE, myFace);
  }
  else
  {
    for (itW = MapIteratorOfMapOfWire(NewWires); itW.More(); itW.Next())
    {
      const TopoDS_Wire& aWire = itW.Value().aWire;
      myNewWires.Append(aWire);
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

  PurgeNewEdges(myCutEdges, UsedEdges);
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
