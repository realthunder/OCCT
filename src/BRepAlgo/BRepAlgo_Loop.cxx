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
#include <BRepTools_History.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>
#include <IntTools_FClass2d.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom_Surface.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomLib.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Ax2.hxx>
#include <Precision.hxx>
#include <ShapeBuild_ReShape.hxx>
#include <ShapeFix_Shape.hxx>
#include <TColStd_SequenceOfReal.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <TopTools.hxx>
#include <TopTools_DataMapIteratorOfDataMapOfShapeShape.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_MapOfShape.hxx>
#include <TopTools_SequenceOfShape.hxx>
#include <BRepBuilderAPI_Copy.hxx>

#include <stdio.h>
//#define OCCT_DEBUG_ALGO
//#define DRAW
#ifdef DRAW
#include <DBRep.hxx>
#pragma comment(lib,"TKDraw")
#endif

#ifdef OCCT_DEBUG_ALGO
Standard_Boolean AffichLoop  = Standard_True;
Standard_Integer NbLoops     = 0;
Standard_Integer NbWires     = 1;
static char* name = new char[100];
#endif

static thread_local Standard_Integer _CollectingEdges;
static thread_local TopTools_DataMapOfShapeListOfShape _EdgeMap;

BRepAlgo_LoopIntersectingEdgeMap::BRepAlgo_LoopIntersectingEdgeMap()
{
  ++_CollectingEdges;
}
BRepAlgo_LoopIntersectingEdgeMap::~BRepAlgo_LoopIntersectingEdgeMap()
{
  if (--_CollectingEdges == 0) {
    _EdgeMap.Clear();
  }
}

TopTools_DataMapOfShapeListOfShape& BRepAlgo_LoopIntersectingEdgeMap::EdgeMap()
{
  return _EdgeMap;
}

//=======================================================================
//function : BRepAlgo_Loop
//purpose  : 
//=======================================================================

BRepAlgo_Loop::BRepAlgo_Loop():
  myTolConf (0.001)
{
}


//=======================================================================
//function : Init
//purpose  : 
//=======================================================================

void BRepAlgo_Loop::Init(const TopoDS_Face& F)
{
  myConstEdges.Clear(); 
  myEdges     .Clear();
  myVerOnEdges.Clear();
  myNewWires  .Clear();
  myNewFaces  .Clear();
  myCutEdges  .Clear();
  myFace = F;
}


//=======================================================================
//function : Bubble
//purpose  : Orders the sequence of vertices by increasing parameter. 
//=======================================================================

static void Bubble(const TopoDS_Edge&        E,
		   TopTools_SequenceOfShape& Seq,
		   TColStd_SequenceOfReal& SeqU) 
{
  Standard_Real    U;
  TopoDS_Vertex    V;

  if (Seq.IsEmpty())
    return;

  //Remove duplicates
  for (Standard_Integer i = 1; i < Seq.Length(); i++) {
    for (Standard_Integer j = i+1; j <= Seq.Length(); j++) {
      if (Seq(i) == Seq(j))
      {
        Seq.Remove(j);
        j--;
      }
    }
    TopoDS_Shape aLocalV = Seq(i)  .Oriented(TopAbs_INTERNAL);
    V = TopoDS::Vertex(aLocalV);
    U = BRep_Tool::Parameter(V,E);
    SeqU.Append(U);
  }
  V = TopoDS::Vertex(Seq.Last().Oriented(TopAbs_INTERNAL));
  U = BRep_Tool::Parameter(V,E);
  SeqU.Append(U);

  Standard_Boolean Invert   = Standard_True;
  Standard_Integer NbPoints = Seq.Length();

  while (Invert) {
    Invert = Standard_False;
    for ( Standard_Integer i = 1; i < NbPoints; i++) {
      if (SeqU(i+1) < SeqU(i)) {
	Seq.Exchange(i,i+1);
        SeqU.Exchange(i,i+1);
	Invert = Standard_True;
      }
    }
  }
}



//=======================================================================
//function : AddEdges
//purpose  : 
//=======================================================================

void BRepAlgo_Loop::AddEdge (TopoDS_Edge&                E, 
			     const TopTools_ListOfShape& LV)
{
  myEdges.Append(E);
  myVerOnEdges.Bind(E,LV);
  SHOW_TOPO_SHAPE(E, "AddEdge", LV);
}


//=======================================================================
//function : AddConstEdges
//purpose  : 
//=======================================================================

void BRepAlgo_Loop::AddConstEdge (const TopoDS_Edge& E)
{
  myConstEdges.Append(E);
}

//=======================================================================
//function : AddConstEdges
//purpose  : 
//=======================================================================

void BRepAlgo_Loop::AddConstEdges(const TopTools_ListOfShape& LE)
{
  TopTools_ListIteratorOfListOfShape itl(LE);
  for (; itl.More(); itl.Next()) {
    myConstEdges.Append(itl.Value());
  }
}

//=======================================================================
//function : SetImageVV
//purpose  : 
//=======================================================================

void BRepAlgo_Loop::SetImageVV (const BRepAlgo_Image& theImageVV)
{
  myImageVV = theImageVV;
}

//=======================================================================
//function : UpdateClosedEdge
//purpose  : If the first or the last vertex of intersection
//           coincides with the closing vertex, it is removed from SV.
//           it will be added at the beginning and the end of SV by the caller.
//=======================================================================

static TopoDS_Vertex  UpdateClosedEdge(const TopoDS_Edge&         E,
				       TopTools_SequenceOfShape&  SV)
{
  TopoDS_Vertex    VB [2], V1, V2, VRes;
  gp_Pnt           P,PC;
  Standard_Boolean OnStart = 0, OnEnd = 0;
  //// modified by jgv, 13.04.04 for OCC5634 ////
  TopExp::Vertices (E,V1,V2);
  Standard_Real    Tol = BRep_Tool::Tolerance( V1 );
  ///////////////////////////////////////////////
  
  if (SV.IsEmpty()) return VRes;

  VB[0] = TopoDS::Vertex(SV.First());
  VB[1] = TopoDS::Vertex(SV.Last ());
  PC = BRep_Tool::Pnt(V1);

  for ( Standard_Integer i = 0 ; i < 2 ; i++) {
    P = BRep_Tool::Pnt(VB [i]);
    if (P.IsEqual(PC,Tol)) {
      VRes = VB [i];
      if (i == 0) OnStart = Standard_True;
      else        OnEnd   = Standard_True;
    }
  }
  if (OnStart && OnEnd) {
    if (!VB[0].IsSame(VB[1])) {
#ifdef OCCT_DEBUG_ALGO
      if (AffichLoop)
	std::cout <<"Two different vertices on the closing vertex"<<std::endl;
#endif
    }
    else {
      SV.Remove(1);
      if (!SV.IsEmpty()) SV.Remove(SV.Length());
    }
  }
  else if (OnStart) SV.Remove(1);
  else if (OnEnd  ) SV.Remove(SV.Length());

  return VRes;
}

//=======================================================================
//function : PurgeNewEdges
//purpose  : 
//=======================================================================

static void  PurgeNewEdges(TopTools_DataMapOfShapeListOfShape& NewEdges,
			   const TopTools_MapOfShape&          UsedEdges)
{
  TopTools_DataMapIteratorOfDataMapOfShapeListOfShape it(NewEdges);
  for (; it.More(); it.Next()) {
    TopTools_ListOfShape& LNE = NewEdges.ChangeFind(it.Key());
    TopTools_ListIteratorOfListOfShape itL(LNE);
    while (itL.More()) {
      const TopoDS_Shape& NE = itL.Value();
      if (!UsedEdges.Contains(NE)) {
	LNE.Remove(itL);
      }
      else {
	itL.Next();
      }
    }
  }
  
}

//=======================================================================
//function : StoreInMVE
//purpose  : 
//=======================================================================

static void StoreInMVE (const TopoDS_Face&                  F,
			TopoDS_Edge&                  E,
			TopTools_IndexedDataMapOfShapeListOfShape& MVE,
			Standard_Boolean&                   YaCouture,
			TopTools_DataMapOfShapeShape& VerticesForSubstitute,
            const Standard_Real theTolConf)
{      
  TopoDS_Vertex V1, V2, V;
  TopTools_ListOfShape Empty;

  gp_Pnt P1, P;
  BRep_Builder BB;
  for (Standard_Integer iV = 1; iV <= MVE.Extent(); iV++)
    {
      V = TopoDS::Vertex(MVE.FindKey(iV));
      P = BRep_Tool::Pnt( V );
      TopTools_ListOfShape VList;
      TopoDS_Iterator VerExp( E );
      for (; VerExp.More(); VerExp.Next())
	VList.Append( VerExp.Value() );
      TopTools_ListIteratorOfListOfShape itl( VList );
      for (; itl.More(); itl.Next())
	{
	  V1 = TopoDS::Vertex( itl.Value() );
	  P1 = BRep_Tool::Pnt( V1 );
	  if (P.IsEqual( P1, theTolConf ) && !V.IsSame(V1))
	    {
	      V.Orientation( V1.Orientation() );
	      if (VerticesForSubstitute.IsBound( V1 ))
		{
		  TopoDS_Shape OldNewV = VerticesForSubstitute( V1 );
		  if (! OldNewV.IsSame( V ))
		    {
		      VerticesForSubstitute.Bind( OldNewV, V );
		      VerticesForSubstitute( V1 ) = V;
		    }
		}
	      else
		{
		  if (VerticesForSubstitute.IsBound( V ))
		    {
		      TopoDS_Shape NewNewV = VerticesForSubstitute( V );
		      if (! NewNewV.IsSame( V1 ))
			VerticesForSubstitute.Bind( V1, NewNewV );
		    }
		  else
		    {
		      VerticesForSubstitute.Bind( V1, V );
		      TopTools_DataMapIteratorOfDataMapOfShapeShape mapit( VerticesForSubstitute );
		      for (; mapit.More(); mapit.Next())
			if (mapit.Value().IsSame( V1 ))
			  VerticesForSubstitute( mapit.Key() ) = V;
		    }
		}
	      E.Free( Standard_True );
	      BB.Remove( E, V1 );
	      BB.Add( E, V );
              SHOW_TOPO_SHAPE(E, "ReplaceMVEE");
              SHOW_TOPO_SHAPE(V, "ReplaceMVEV");
	    }
	}
    }

  TopExp::Vertices(E,V1,V2);
  if( V1.IsNull() && V2.IsNull() ){ YaCouture = Standard_False; return; }
  if (!MVE.Contains(V1)) {
    MVE.Add(V1,Empty);
  }
  MVE.ChangeFromKey(V1).Append(E);
  // SHOW_TOPO_SHAPE(V1, "MVE_V1_");
  if (!V1.IsSame(V2)) {
     if (!MVE.Contains(V2)) {
       MVE.Add(V2,Empty);
     }
     MVE.ChangeFromKey(V2).Append(E);
     // SHOW_TOPO_SHAPE(V2, "MVE_V2_");
  }
  TopLoc_Location L ;
  Handle(Geom_Surface) S = BRep_Tool::Surface(F,L);
  if (BRep_Tool::IsClosed(E,S,L)) {
    MVE.ChangeFromKey(V2).Append(E.Reversed());
    // SHOW_TOPO_SHAPE(V2, "MVE_ClosedV2_");
    if (!V1.IsSame(V2)) {
      MVE.ChangeFromKey(V1).Append(E.Reversed());
      // SHOW_TOPO_SHAPE(V1, "MVE_ClosedV1_");
    }
    YaCouture = Standard_True;
  }
  // SHOW_TOPO_SHAPE(E, "MVE");
}

//=======================================================================
//function : Perform
//purpose  : 
//=======================================================================

void BRepAlgo_Loop::Perform()
{
    Perform(nullptr);
}

void BRepAlgo_Loop::Perform(const TopTools_ListOfShape* ContextFaces,
                            const Handle(BRepAlgo_AsDes)& AsDes)
{
  TopTools_ListIteratorOfListOfShape                  itl, itl1, itl2, itl3;
  TopoDS_Vertex                                       V1,V2,OV1,OV2;
  BRep_Builder                                        B;
  TopExp_Explorer                                     aExp;

  //------------------------------------------------
  // Check intersection in myConstEdges which is possible when make thick solid
  // with concave removed face
  //------------------------------------------------
  TopTools_ListOfShape theEdges = myEdges;
  TopTools_ListOfShape ConstEdges;
  TopTools_ListOfShape IntersectingEdges;
  if (_CollectingEdges) {
    TopTools_MapOfShape EMap;
    TopTools_ListOfShape theVerts;
    TopTools_ListOfShape LV;
    TopTools_MapOfShape MV;

    for (itl.Initialize(myConstEdges); itl.More(); itl.Next())
      theEdges.Append(itl.Value());

    for (itl.Initialize(theEdges); itl.More(); itl.Next())
    {
      TopoDS_Edge anEdge = TopoDS::Edge(itl.Value().Oriented(TopAbs_FORWARD));
      // Sewn edges can be doubled or not in myConstEdges
      if (!EMap.Add(anEdge))
          continue;

      SHOW_TOPO_SHAPE(anEdge, "InterEdge");

      LV.Clear();
      MV.Clear();

      Standard_Boolean Bounded = Standard_False;
      Standard_Real FP, LP;
      TopoDS_Shape VF, VL;
      if (myVerOnEdges.IsBound(anEdge)) {
        const TopTools_ListOfShape& LE = myVerOnEdges(anEdge);
        Bounded = Standard_True;
        for (itl1.Initialize(LE); itl1.More(); itl1.Next()) {
          if (!MV.Add(itl1.Value()))
            continue;
          if (itl1.Value().Orientation() == TopAbs_FORWARD ||
              itl1.Value().Orientation() == TopAbs_REVERSED) {
            const TopoDS_Vertex& aVertex = TopoDS::Vertex(itl1.Value());
	    Standard_Real P = BRep_Tool::Parameter(aVertex, anEdge);
            if (VF.IsNull()) {
              FP = LP = P;
              VF = aVertex.Oriented(TopAbs_FORWARD);
              VL = aVertex.Oriented(TopAbs_REVERSED);
            }
            else if (FP > P) {
              VF = aVertex.Oriented(TopAbs_FORWARD);
              FP = P;
            }
            else if (LP < P) {
              VL = aVertex.Oriented(TopAbs_REVERSED);
              LP = P;
            }
          }
        }
      }

      MV.Clear();
      if (!VF.IsNull()){
        MV.Add(VF);
        LV.Append(VF);
      }
      if (!VL.IsNull() && !VL.IsSame(VF)) {
        MV.Add(VL);
        LV.Append(VL);
      }

      Standard_Boolean Extended = Standard_False;
      Standard_Real aF, aL;
      const Handle(Geom_Curve) C = BRep_Tool::Curve(anEdge, aF, aL);
      TopExp::Vertices(anEdge, V1, V2);

      for (TopoDS_Iterator It(anEdge); It.More(); It.Next()) {
        if (It.Value().Orientation() == TopAbs_INTERNAL) {
          Extended = Standard_True;
          break;
        }
      }

      for (itl1.Initialize(LV); itl1.More(); itl1.Next()) {
        SHOW_TOPO_SHAPE(itl1.Value(), "InterV1", 1);
      }

      for (itl1.Initialize(theEdges); itl1.More(); itl1.Next()) {
        const TopoDS_Edge& otherEdge = TopoDS::Edge(itl1.Value());
        // YES, we will check vertices from the same edge just like any other
        // edges in order to determine the orientation of the internal
        // vertices. So no need to the check IsSame() here.
        //
        // if (otherEdge.IsSame(anEdge))
        //   continue;

        theVerts.Clear();
        const TopTools_ListOfShape *pLV = myVerOnEdges.Seek(otherEdge);
        if (pLV)
          theVerts = *pLV;
        TopExp::Vertices(otherEdge, OV1, OV2);
        theVerts.Append(OV1);
        theVerts.Append(OV2);
        for (itl2.Initialize(theVerts); itl2.More(); itl2.Next()) {
          TopoDS_Vertex aVertex = TopoDS::Vertex(itl2.Value());
          if (!MV.Add(aVertex))
            continue;
          Standard_Real Tol = BRep_Tool::Tolerance(aVertex);
          gp_Pnt OP = BRep_Tool::Pnt(aVertex);
          if (OP.Distance(BRep_Tool::Pnt(V1)) < Tol
              || OP.Distance(BRep_Tool::Pnt(V2)) < Tol)
            continue;
          if (Extended) {
            if (OP.Distance(BRep_Tool::Pnt(TopoDS::Vertex(VF))) < Tol
                || OP.Distance(BRep_Tool::Pnt(TopoDS::Vertex(VL))) < Tol)
                continue;
          }
          if (C.IsNull())
              continue;
          GeomAPI_ProjectPointOnCurve Proj(BRep_Tool::Pnt(aVertex), C);
          if (Proj.NbPoints() > 0) {
            Standard_Real D = Proj.LowerDistance();
            Standard_Real P = Proj.LowerDistanceParameter();
            auto ReorientVertex = [&](TopoDS_Shape &V, TopAbs_Orientation Ori) {
              if (V.Orientation() == Ori)
                return;
              V.Orientation(Ori);
              for (itl3.Initialize(LV); itl3.More(); itl3.Next()) {
                if (itl3.Value().IsSame(V)) {
                  itl3.Value().Orientation(Ori);
                  SHOW_TOPO_SHAPE(V, "VertexFlip");
                  break;
                }
              }
            };
            if (C->IsPeriodic()) {
              while (P < aF) {
                P += C->Period();
              }
            }
            if (D < Tol  && P > aF && P < aL) {
              TopoDS_Shape aLocalShape;
              if (Extended) {
                if (P < FP) {
                  aLocalShape = aVertex.Oriented(TopAbs_FORWARD);
                  SHOW_TOPO_SHAPE(aLocalShape, "VertexOverF");
                  ReorientVertex(VF, TopAbs_REVERSED);
                }
                else if (P > LP) {
                  aLocalShape = aVertex.Oriented(TopAbs_REVERSED);
                  SHOW_TOPO_SHAPE(aLocalShape, "VertexOverR");
                  ReorientVertex(VL, TopAbs_FORWARD);
                }
                else if (P < (FP + LP)/2) {
                  aLocalShape = aVertex.Oriented(TopAbs_REVERSED);
                  ReorientVertex(VF, TopAbs_FORWARD);
                }
                else {
                  aLocalShape = aVertex.Oriented(TopAbs_FORWARD);
                  ReorientVertex(VL, TopAbs_REVERSED);
                }
              }
              else if (P < (aF + aL)/2) {
                aLocalShape = aVertex.Oriented(TopAbs_REVERSED);
              }
              else {
                aLocalShape = aVertex.Oriented(TopAbs_FORWARD);
              }
              B.UpdateVertex(TopoDS::Vertex(aLocalShape),P,anEdge,Tol);
              LV.Append(aLocalShape);
              SHOW_TOPO_SHAPE(aLocalShape, "InterV", 1);
            }
            else if (D < Tol) {
              SHOW_TOPO_SHAPE(aVertex, "InterVSkip", 1);
            }
          }
        }
      }
      if (LV.Extent()) {
        if (!Bounded) {
          IntersectingEdges.Append(anEdge);
          myVerOnEdges.Bind(anEdge, LV);
        }
        else {
          *myVerOnEdges.ChangeSeek(anEdge) = LV;
        }
      }
      else {
        ConstEdges.Append(anEdge);
        SHOW_TOPO_SHAPE(anEdge, "ConstEdge");
      }
    }
    if (ConstEdges.Extent() != myConstEdges.Extent()) {
      myConstEdges = ConstEdges;
    }
  }

#ifdef OCCT_DEBUG_ALGO
  if (AffichLoop) {
    std::cout <<"NewLoop"<<std::endl;
    NbLoops++;
#ifdef DRAW
    sprintf(name,"FLoop_%d",NbLoops);
    DBRep::Set(name,myFace);
    Standard_Integer NbEdges = 1;
    for (itl.Initialize(myEdges); itl.More(); itl.Next()) { 
      const TopoDS_Edge& E = TopoDS::Edge(itl.Value());
      sprintf(name,"EEE_%d_%d",NbLoops,NbEdges++);
      DBRep::Set(name,E);
    }
    for (itl.Initialize(myConstEdges); itl.More(); itl.Next()) {
      const TopoDS_Edge& E = TopoDS::Edge(itl.Value());    
      sprintf(name,"EEC_%d_%d_",NbLoops,NbEdges++);
      DBRep::Set(name,E);
    }
#endif
  }
#endif

  //------------------------------------------------
  // Cut edges
  //------------------------------------------------
  for (itl.Initialize(theEdges); itl.More(); itl.Next())
  {
    const TopoDS_Edge& anEdge = TopoDS::Edge(itl.Value());
    if (myCutEdges.Seek(anEdge))
      continue;
    TopTools_ListOfShape LCE;
    const TopTools_ListOfShape* pVertices = myVerOnEdges.Seek (anEdge);
    if (pVertices)
    {
      Standard_Boolean KeepAll = Standard_True;
      if (ContextFaces && AsDes && AsDes->HasAscendant(anEdge)) {
        const TopTools_ListOfShape LF = AsDes->Ascendant(anEdge);
        Standard_Integer Count = 0;
        for (itl1.Initialize(LF); itl1.More(); itl1.Next()) {
          for (itl2.Initialize(*ContextFaces); itl2.More(); itl2.Next()) {
            if (itl2.Value().IsSame(itl1.Value()) && ++Count > 1) {
              KeepAll = Standard_False;
              SHOW_TOPO_SHAPE(anEdge, "NoKeepAll");
              break;
            }
          }
          if (!KeepAll)
            break;
        }
      }
      CutEdge (anEdge, *pVertices, LCE, KeepAll);
      myCutEdges.Bind(anEdge, LCE);
    }
  }

  FindLoop();

  if (_CollectingEdges) {
#if 1
    TopTools_DataMapIteratorOfDataMapOfShapeListOfShape itM(myCutEdges);
    for (; itM.More(); itM.Next()) {
      SHOW_TOPO_SHAPE(itM.Key(), "InterNewEdge", itM.Value());
      TopTools_ListOfShape *pLE = _EdgeMap.ChangeSeek(itM.Key());
      if (!pLE)
        pLE = _EdgeMap.Bound(itM.Key(), TopTools_ListOfShape());
      for (itl.Initialize(itM.Value()); itl.More(); itl.Next())
        pLE->Append(itl.Value());
    }
#else
    for (itl.Initialize(IntersectingEdges); itl.More(); itl.Next()) {
      const TopTools_ListOfShape &aList = NewEdges(TopoDS::Edge(itl.Value()));
      _EdgeMap.Bind(itl.Value(), aList);
      SHOW_TOPO_SHAPE(itl.Value(), "InterNewEdge", aList);
    }
#endif
  }
}

namespace {

typedef NCollection_DataMap<TopoDS_Shape,TopTools_MapOfShape,TopTools_ShapeMapHasher> DataMapOfMapOfShape;
typedef NCollection_DataMap<TopoDS_Shape,TopTools_MapOfShape,TopTools_ShapeMapHasher>::Iterator DataMapIteratorOfMapOfShape;

// boost::hash_combine
inline Standard_Integer combine(Standard_Integer seed, Standard_Integer h) noexcept
{
    seed ^= h + 0x9e3779b9 + (seed << 6U) + (seed >> 2U);
    return seed;
}

struct WireInfo
{
  mutable TopoDS_Wire aWire;
  std::vector<std::pair<TopoDS_Shape,Standard_Integer>> Edges;
  Standard_Integer aHashCode;
  mutable TopoDS_Face aFace;

  WireInfo(const TopoDS_Wire &W)
    : aWire(W)
  {
    TopoDS_Iterator aIt(W);
    for (; aIt.More(); aIt.Next()) {
      Edges.emplace_back(aIt.Value(), aIt.Value().HashCode(INT_MAX));
    }
    std::sort(Edges.begin(), Edges.end(),
      [](const std::pair<TopoDS_Shape, Standard_Integer>& a,
         const std::pair<TopoDS_Shape, Standard_Integer>& b)
      {
        if (a.first.TShape() < b.first.TShape())
          return true;
        if (b.first.TShape() < a.first.TShape())
          return false;
        return a.first.Location().HashCode(INT_MAX) < b.first.Location().HashCode(INT_MAX);
      });
    aHashCode = 0;
    for (const auto& s : Edges)
      aHashCode = combine(aHashCode, s.second);
  }

  WireInfo(const WireInfo &other)
    :aWire(other.aWire)
    ,Edges(other.Edges)
    ,aHashCode(other.aHashCode)
  {}
};

struct WireInfoHasher
{
  static Standard_Integer HashCode(const WireInfo &info, Standard_Integer theUpperBound)
  {
    return ::HashCode(info.aHashCode, theUpperBound);
  }

  static Standard_Boolean IsEqual (const WireInfo& a, const WireInfo& b)
  {
    if (a.Edges.size() != b.Edges.size())
      return Standard_False;
    std::size_t i;
    for (i=0; i<a.Edges.size(); ++i)
    {
      if (!a.Edges[i].first.IsSame(b.Edges[i].first))
        return Standard_False;
    }
    return Standard_True;
  }
};

typedef NCollection_Map<WireInfo,WireInfoHasher> MapOfWire;
typedef NCollection_Map<WireInfo,WireInfoHasher>::Iterator MapIteratorOfMapOfWire;

void FindAllLoops(const TopoDS_Vertex& CV,
                  const TopoDS_Edge& CE,
                  TopTools_DataMapOfShapeShape& CurrentVEMap, 
                  TopTools_ListOfShape& CurrentEdgeList, 
                  const TopTools_IndexedDataMapOfShapeListOfShape& MVE,
                  MapOfWire& NewWires,
                  const TopoDS_Face& aFace,
                  const Standard_Real &Tol)
{
  TopTools_ListIteratorOfListOfShape itl;
  TopoDS_Vertex                      V1, V2, NV;

  CurrentVEMap.Bind(CV, CE);
  CurrentEdgeList.Prepend(CE);
  // SHOW_TOPO_SHAPE(CV, "CV");
  // SHOW_TOPO_SHAPE(CE, "CE");

  TopExp::Vertices(CE, V1, V2);
  if (CV.IsSame(V1))
    NV = V2;
  else
    NV = V1;

  TopoDS_Edge EF;
  const TopoDS_Shape *pE = CurrentVEMap.Seek(NV);
  if (pE)
    EF = TopoDS::Edge(*pE);
  if (pE && !EF.IsSame(CE)) {
    TopoDS_Edge E = EF;
    TopoDS_Vertex VF, VL;
    BRepLib_MakeWire aMakeWire;
    for (;;)
    {
      TopExp::Vertices(E, V1, V2, Standard_False);
      if (VF.IsNull()) {
        VF = NV;
        VL = NV.IsSame(V1) ? V2 : V1;
        // SHOW_TOPO_SHAPE(VF, "NWVF");
      }
      else if (VL.IsSame(V1))
        VL = V2;
      else
        VL = V1;
      // SHOW_TOPO_SHAPE(VL, "NWV");
      // SHOW_TOPO_SHAPE(E, "NWE");
      aMakeWire.Add(E);
      if (VL.IsSame(VF))
        break;
      E = TopoDS::Edge(CurrentVEMap.Find(VL));
    }
    TopoDS_Wire NW = aMakeWire.Wire();
    if (NW.Closed() && NewWires.Add(NW)) {
        SHOW_TOPO_SHAPE(NW, "NewWire");
    }
    else {
        // SHOW_TOPO_SHAPE(NW, "DiscardWire");
    }
  }
  else
  {
    for (itl.Initialize(MVE.FindFromKey(NV)); itl.More(); itl.Next())
    {
      const TopoDS_Edge& NE = TopoDS::Edge(itl.Value());
      if (!NE.IsSame(CE) && !EF.IsSame(NE))
        FindAllLoops(NV, NE, CurrentVEMap, CurrentEdgeList, MVE, NewWires, aFace, Tol);
    }
  }

  CurrentVEMap.UnBind(CV);
  CurrentEdgeList.RemoveFirst();
  // SHOW_TOPO_SHAPE(CE, "Pop");
}

void SplitWires(TopTools_ListOfShape& OutputWires,
                MapOfWire& InputWires,
                TopTools_IndexedDataMapOfShapeListOfShape& MVE,
                const TopoDS_Face& aFace)
{
  TopTools_MapOfShape UsedEdges;
  TopTools_DataMapOfShapeListOfShape aEdgeWireMap;
  TopTools_IndexedMapOfShape aWireEdgeMap;
  TopTools_ListIteratorOfListOfShape itl;
  MapIteratorOfMapOfWire itM;
  TopTools_MapOfShape aPrunedMap;
  TopTools_MapOfShape aCheckMap;
  TopExp_Explorer aExp;
  BRep_Builder B;
  TopLoc_Location L;
  Standard_Real Tol = BRep_Tool::Tolerance(aFace);

  for (itM.Initialize(InputWires); itM.More(); itM.Next())
  {
    for (const auto &v : itM.Value().Edges)
      UsedEdges.Add(v.first);
  }

  for (itM.Initialize(InputWires); itM.More(); itM.Next())
  {
    IntTools_FClass2d aClassifier;
    const WireInfo& Info = itM.Value();
    TopoDS_Wire& aWire = TopoDS::Wire(Info.aWire);
    if (aPrunedMap.Contains(aWire))
      continue;

    TopoDS_Face& aNF = Info.aFace;
    if (aNF.IsNull())
    {
      Handle(Geom_Surface) S = BRep_Tool::Surface(aFace, L);
      B.MakeFace(aNF, S, L, Tol);
      B.Add(aNF, aWire);
      B.NaturalRestriction(aNF, Standard_False);

      BRepTopAdaptor_FClass2d FClass2d(aNF,Precision::PConfusion());
      if(FClass2d.PerformInfinitePoint() != TopAbs_OUT) { 
        aNF.EmptyCopy();
        aWire.Reverse();
        B.Add(aNF, aWire);
        B.NaturalRestriction(aNF, Standard_False);
      }

      BRepCheck_Analyzer anAnalyzer(aNF);
      if (!anAnalyzer.IsValid())
      {
        ShapeFix_Shape aFix(aNF);
        aFix.Perform();
        aNF = TopoDS::Face(aFix.Shape());
      }
      aNF = TopoDS::Face(aNF.Oriented(TopAbs_FORWARD));
      SHOW_TOPO_SHAPE(aNF, "MakeFace");

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
      const TopoDS_Vertex& aVertex(TopoDS::Vertex(aExp.Current()));
      const TopTools_ListOfShape& EL = MVE.FindFromKey(aVertex);
      if (EL.Extent() <= 2)
        continue;

      aCheckMap.Clear();

      aWireEdgeMap.Clear();
      TopExp::MapShapes(aWire, TopAbs_EDGE, aWireEdgeMap);

      for (itl.Initialize(EL); itl.More(); itl.Next())
      {
        const TopoDS_Edge &aEdge = TopoDS::Edge(itl.Value());

        if (aWireEdgeMap.Contains(aEdge))
          continue;

        if (!UsedEdges.Contains(aEdge))
          continue;

        if (!aCheckMap.Add(aEdge))
          continue;

        SHOW_TOPO_SHAPE(aEdge, "CheckSplit");

        // Get 2d curve of the edge on the face
        Standard_Real aT1, aT2;
        const Handle(Geom2d_Curve)& aC2D = BRep_Tool::CurveOnSurface(aEdge, aNF, aT1, aT2);
        if (aC2D.IsNull()) {
          SHOW_TOPO_SHAPE(aEdge, "Prune_NoCurve_");
          continue;
        }

        // Get middle point on the curve
        gp_Pnt2d aP2D = aC2D->Value((aT1 + aT2) / 2.);

        // Classify the point
        TopAbs_State aState = aClassifier.Perform(aP2D);

        static int index;
        if (index++ == 10) {
          Standard_Real aT11, aT12;
          const Handle(Geom2d_Curve)& aC2D1 = BRep_Tool::CurveOnSurface(aEdge, aFace, aT11, aT12);
          gp_Pnt2d aP2D1 = aC2D1->Value((aT11 + aT12) / 2.);
          aClassifier.Perform(aP2D1);
        }

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
        break;
    }
  }

  for (itM.Initialize(InputWires); itM.More(); itM.Next())
  {
    if (!aPrunedMap.Contains(itM.Value().aWire))
    {
      for (const auto &v : itM.Value().Edges)
      {
        TopTools_ListOfShape *pLE = aEdgeWireMap.ChangeSeek(v.first);
        if (!pLE)
          pLE = aEdgeWireMap.Bound(v.first, TopTools_ListOfShape());
        pLE->Append(itM.Value().aWire);
      }
    }
  }

  for (itM.Initialize(InputWires); itM.More(); itM.Next())
  {
    const TopoDS_Wire& aWire = itM.Value().aWire;
    if (aPrunedMap.Contains(aWire))
      continue;

    //-----------------------------------------------
    // Prune the wire if all of its edges are shared by some other wire, in
    // which case the wire can be interpreted as inner hole)
    //----------------------------------------------

    aCheckMap.Clear();
    Standard_Boolean Pruned = Standard_True;
    for (const auto &v : itM.Value().Edges)
    {
      const TopTools_ListOfShape& LE = aEdgeWireMap.Find(v.first);
      if (LE.Extent() == 1)
      {
        OutputWires.Append(aWire);
        Pruned = Standard_False;
        break;
      }
      for (itl.Initialize(LE); itl.More(); itl.Next())
      {
        if (!itl.Value().IsSame(aWire))
          aCheckMap.Add(itl.Value());
      }
    }
    if (Pruned)
    {
      SHOW_TOPO_SHAPE(aWire, "Prune1_");
      TopTools_MapIteratorOfMapOfShape itM1(aCheckMap);
      for (; itM1.More(); itM1.Next())
        SHOW_TOPO_SHAPE(itM1.Value(), "Prune1W_");
      // aPrunedMap.Add(aWire);
    }
  }
}

} // Anonymous namespace

void BRepAlgo_Loop::FindLoop()
{
  TopTools_ListIteratorOfListOfShape itl,  itl1;
  Standard_Boolean   YaCouture = Standard_False;

  myNewWires.Clear();
  myNewFaces.Clear();

  //-----------------------------------
  // Construction map vertex => edges
  //-----------------------------------
  TopTools_IndexedDataMapOfShapeListOfShape MVE;

  // add cut edges.
  TopTools_MapOfShape Emap;
  TopTools_DataMapIteratorOfDataMapOfShapeListOfShape itM(myCutEdges);
  for (; itM.More(); itM.Next()) {
    for (itl1.Initialize(itM.Value()); itl1.More(); itl1.Next()) {
      TopoDS_Edge& E = TopoDS::Edge(itl1.Value());
      if (Emap.Add(E))
        StoreInMVE(myFace,E,MVE,YaCouture,myVerticesForSubstitute, myTolConf);
    }
  }
  
  // add const edges
  // Sewn edges can be doubled or not in myConstEdges
  // => call only once StoreInMVE which should double them
  TopTools_MapOfShape DejaVu;
  for (itl.Initialize(myConstEdges); itl.More(); itl.Next()) {
    TopoDS_Edge& E = TopoDS::Edge(itl.Value());
    if (DejaVu.Add(E))
      StoreInMVE(myFace,E,MVE,YaCouture,myVerticesForSubstitute, myTolConf);
  }

#ifdef DRAW
  if (AffichLoop) {
    std::cout <<"NewLoop"<<std::endl;
    Standard_Integer NbEdges = 1;
    TopTools_MapOfShape Done;
    for (Standard_Integer iV = 1; iV <= MVE.Extent(); iV++) {
      for (itl.Initialize(MVE(iV)); itl.More(); itl.Next()) {
        TopoDS_Edge& E = TopoDS::Edge(itl.Value());
        if (Done.Add(E)) {
          sprintf(name,"EEC2_%d_%d_",NbLoops,NbEdges++);
          DBRep::Set(name,E);
        }
      }
    }
  }
#endif

  UpdateVEmap (MVE);
  if (MVE.IsEmpty())
    return;

  TopTools_DataMapOfShapeShape CurrentVEMap;
  TopTools_ListOfShape CurrentEdgeList;
  MapOfWire NewWires;

  //-----------------------------------------------
  // Find all possible closed wires
  //----------------------------------------------

  const TopoDS_Vertex& VF = TopoDS::Vertex(MVE.FindKey(1));
#if 0
  for (Standard_Integer ii = 1; ii <= MVE.Extent(); ++ii)
  {
    SHOW_TOPO_SHAPE(V, "MVEV");
    for (itl.Initialize(MVE(ii)); itl.More(); itl.Next())
      SHOW_TOPO_SHAPE(E, "MVE");
  }
#endif

  for (Standard_Integer ii = 1; ii <= MVE.Extent(); ++ii)
  {
    for (itl.Initialize(MVE(ii)); itl.More(); itl.Next())
    {
      FindAllLoops(VF, TopoDS::Edge(itl.Value()), CurrentVEMap,
          CurrentEdgeList, MVE, NewWires, myFace, myTolConf);
    }
  }

  //-----------------------------------------------
  // Split wires
  //----------------------------------------------
  SplitWires(myNewWires, NewWires, MVE, myFace);

  TopTools_MapOfShape UsedEdges;
  TopExp_Explorer aExp;
  for (itl.Initialize(myNewWires); itl.More(); itl.Next())
  {
    for (aExp.Init(itl.Value(), TopAbs_EDGE); aExp.More(); aExp.Next())
      UsedEdges.Add(aExp.Current());
  }

  PurgeNewEdges(myCutEdges,UsedEdges);
}

//=======================================================================
//function : CutEdges
//purpose  : 
//=======================================================================

void BRepAlgo_Loop::CutEdge (const TopoDS_Edge&          E,
			     const TopTools_ListOfShape& VOnE,
			             TopTools_ListOfShape& NE   ) const 
{
  CutEdge(E, VOnE, NE, Standard_True);
}


void BRepAlgo_Loop::CutEdge (const TopoDS_Edge&          E,
                             const TopTools_ListOfShape& VOnE,
                             TopTools_ListOfShape& NE,
                             Standard_Boolean KeepAll) const 
{
  Standard_Real Tol = 0.001; //5.e-05; //5.e-07;
  TopoDS_Shape aLocalE  = E.Oriented(TopAbs_FORWARD);
  TopoDS_Edge WE = TopoDS::Edge(aLocalE);

  Standard_Real                      U1,U2;
  TopoDS_Vertex                      V1,V2;
  TopTools_SequenceOfShape           SV;
  TColStd_SequenceOfReal             SU;
  TopTools_ListIteratorOfListOfShape it(VOnE);
  BRep_Builder                       B;

  for ( ; it.More(); it.Next()) {
    SV.Append(it.Value());
  }
  //--------------------------------
  // Parse vertices on the edge.
  //--------------------------------
  Bubble (WE,SV,SU);

  // KeepAll = Standard_True;

  if (KeepAll) {
    SHOW_TOPO_SHAPE(WE, "CuttingInter");
  }
  else {
    SHOW_TOPO_SHAPE(WE, "Cutting");
  }

  Standard_Integer NbVer = SV.Length();
  //----------------------------------------------------------------
  // Construction of new edges.
  // Note :  vertices at the extremities of edges are not 
  //         onligatorily in the list of vertices
  //----------------------------------------------------------------
  if (SV.IsEmpty()) {
    NE.Append(E);
    return;
  }

  for (Standard_Integer ii = 1; ii <= SV.Length(); ++ii) {
    SHOW_TOPO_SHAPE(SV(ii), "CuttingV", 1);
  }

  TopoDS_Vertex    VF,VL;
  Standard_Real    f,l;
  Handle(Geom2d_Curve) C = BRep_Tool::CurveOnSurface(WE,myFace,f,l);  
  TopExp::Vertices(WE,VF,VL);

  TopoDS_Iterator It(WE);
  Standard_Boolean Extended = Standard_False;
  for (; It.More(); It.Next())
  {
    if (It.Value().Orientation() == TopAbs_INTERNAL) {
      Extended = Standard_True;
      break;
    }
  }

  if (NbVer == 2) {
    if (SV(1).IsEqual(VF) && SV(2).IsEqual(VL)) {
      NE.Append(E);
#ifdef DRAW
      if (AffichLoop) {  
      DBRep::Set("ECOpied",E);
    }      
#endif
      return;
    }
  }

  auto InsertVertex = [&](const TopoDS_Shape &V, Standard_Real U) {
    for (Standard_Integer ii = 1; ii <= SV.Length(); ++ii) {
      if (SV(ii).IsSame(V))
        return;
      if (SU(ii) > U) {
        SHOW_TOPO_SHAPE(V, "CuttingInsV", 1);
        SU.InsertBefore(ii, U);
        SV.InsertBefore(ii, V);
        return;
      }
    }
    SHOW_TOPO_SHAPE(V, "CuttingInsV", 1);
    SU.Append(U);
    SV.Append(V);
  };

  //----------------------------------------------------
  // Processing of closed edges 
  // If a vertex of intersection is on the common vertex
  // it should appear at the beginning and end of SV.
  //----------------------------------------------------
  TopoDS_Vertex VCEI;
  if (!VF.IsNull() && VF.IsSame(VL)) {
    VCEI = UpdateClosedEdge(WE,SV);    
    if (!VCEI.IsNull()) {
      TopoDS_Shape aLocalV = VCEI.Oriented(TopAbs_FORWARD);
      VF = TopoDS::Vertex(aLocalV);
      aLocalV = VCEI.Oriented(TopAbs_REVERSED); 
      VL = TopoDS::Vertex(aLocalV);
      InsertVertex(VF, f);
      InsertVertex(VL, l);
    }
    else if (!Extended) {
      InsertVertex(VF, f);
      InsertVertex(VL, l);
    }
  }
  else if (!Extended) {
    //-----------------------------------------
    // Eventually all extremities of the edge.
    //-----------------------------------------
    if (!VF.IsNull())
      InsertVertex(VF, f);
    if (!VL.IsNull())
      InsertVertex(VL, l);
  }

  while (!SV.IsEmpty()) {
    while (!KeepAll && !SV.IsEmpty() && 
	   SV.First().Orientation() != TopAbs_FORWARD) {
      SV.Remove(1);
    }
    if (SV.IsEmpty())
      break;
    V1  = TopoDS::Vertex(SV.First());
    SV.Remove(1);
    if (SV.IsEmpty())
      break;
    if (SV.First().IsSame(V1))
      continue;
    if (KeepAll || SV.First().Orientation() == TopAbs_REVERSED) {
      V2  = TopoDS::Vertex(SV.First());
      //-------------------------------------------
      // Copy the edge and restriction by V1 V2.
      //-------------------------------------------
      TopoDS_Shape NewEdge = WE.EmptyCopied();
      TopoDS_Shape aLocalEdge = V1.Oriented(TopAbs_FORWARD);
      B.Add  (NewEdge,aLocalEdge);
      aLocalEdge = V2.Oriented(TopAbs_REVERSED);
      B.Add  (TopoDS::Edge(NewEdge),aLocalEdge);
      if (V1.IsSame(VF)) 
	U1 = f;
      else 
	{
	  TopoDS_Shape aLocalV = V1.Oriented(TopAbs_INTERNAL);
	  U1=BRep_Tool::Parameter(TopoDS::Vertex(aLocalV),WE);
	}
      if (V2.IsSame(VL))
	U2 = l;
      else
	{
	  TopoDS_Shape aLocalV = V2.Oriented(TopAbs_INTERNAL);
	  U2=BRep_Tool::Parameter(TopoDS::Vertex(aLocalV),WE);
	}
      B.Range (TopoDS::Edge(NewEdge),U1,U2);
#ifdef DRAW
    if (AffichLoop) {  
      DBRep::Set("Cut",NewEdge);
    }
#endif
      SHOW_TOPO_SHAPE(NewEdge, "CutEdge");
      NE.Append(NewEdge.Oriented(E.Orientation()));
    }
  }

  //Remove edges with size <= tolerance
  it.Initialize(NE);
  while (it.More())
    {
      // skl : I change "E" to "EE"
      TopoDS_Edge EE = TopoDS::Edge( it.Value() );
      Standard_Real fpar, lpar;
      BRep_Tool::Range( EE, fpar, lpar );
      if (lpar - fpar <= Precision::Confusion()) {
        SHOW_TOPO_SHAPE(EE, "CutEdgeRemove");
	NE.Remove(it);
      }
      else
	{
	  gp_Pnt2d pf, pl;
	  BRep_Tool::UVPoints( EE, myFace, pf, pl );
	  if (pf.Distance(pl) <= Tol && !BRep_Tool::IsClosed(EE)) {
            SHOW_TOPO_SHAPE(EE, "CutEdgeRemove");
	    NE.Remove(it);
          }
	  else
	    it.Next();
	}
    }
}

//=======================================================================
//function : NewWires
//purpose  : 
//=======================================================================

const TopTools_ListOfShape&  BRepAlgo_Loop::NewWires() const 
{  
  return myNewWires;
}

//=======================================================================
//function : NewFaces
//purpose  : 
//=======================================================================

const TopTools_ListOfShape&  BRepAlgo_Loop::NewFaces() const 
{  
  return myNewFaces;
}
 
//=======================================================================
//function : WiresToFaces
//purpose  : 
//=======================================================================

void  BRepAlgo_Loop::WiresToFaces() 
{  
  if (!myNewWires.IsEmpty()) {
    BRepAlgo_FaceRestrictor FR;
    TopoDS_Shape aLocalS = myFace.Oriented(TopAbs_FORWARD);
    FR.Init (TopoDS::Face(aLocalS),Standard_False, Standard_True);
//    FR.Init (TopoDS::Face(myFace.Oriented(TopAbs_FORWARD)),
//	     Standard_False);
    TopTools_ListIteratorOfListOfShape it(myNewWires);
    for (; it.More(); it.Next()) {
      FR.Add(TopoDS::Wire(it.Value()));
    }

    FR.Perform();
    
    if (FR.IsDone()) {
      TopAbs_Orientation OriF = myFace.Orientation();
      for (; FR.More(); FR.Next()) {
	myNewFaces.Append(FR.Current().Oriented(OriF));
      }
    }
    SHOW_TOPO_SHAPE(TopoDS_Shape(), "NewFace", myNewFaces);
  }
}


//=======================================================================
//function : NewEdges
//purpose  : 
//=======================================================================

const TopTools_ListOfShape&  BRepAlgo_Loop::NewEdges(const TopoDS_Edge& E) const 
{
  return myCutEdges(E);
}

//=======================================================================
//function : GetVerticesForSubstitute
//purpose  : 
//=======================================================================

void  BRepAlgo_Loop::GetVerticesForSubstitute( TopTools_DataMapOfShapeShape& VerVerMap ) const
{
  VerVerMap = myVerticesForSubstitute;
}

//=======================================================================
//function : VerticesForSubstitute
//purpose  : 
//=======================================================================

void  BRepAlgo_Loop::VerticesForSubstitute( TopTools_DataMapOfShapeShape& VerVerMap )
{
  myVerticesForSubstitute = VerVerMap;
}

//=======================================================================
//function : UpdateVEmap
//purpose  : 
//=======================================================================

void  BRepAlgo_Loop::UpdateVEmap (TopTools_IndexedDataMapOfShapeListOfShape& theVEmap)
{
  TopTools_IndexedDataMapOfShapeListOfShape VerLver;

  for (Standard_Integer ii = 1; ii <= theVEmap.Extent(); ii++)
  {
    const TopoDS_Vertex& aVertex = TopoDS::Vertex (theVEmap.FindKey(ii));
    const TopTools_ListOfShape& aElist = theVEmap(ii);
    if (aElist.Extent() == 1 && myImageVV.IsImage(aVertex))
    {
      const TopoDS_Vertex& aProVertex = TopoDS::Vertex (myImageVV.ImageFrom(aVertex));
      SHOW_TOPO_SHAPE(aElist.First(), "VEMapE");
      SHOW_TOPO_SHAPE(aVertex, "VEMapV");
      SHOW_TOPO_SHAPE(aProVertex, "VEMapVI");
      if (VerLver.Contains(aProVertex))
      {
        TopTools_ListOfShape& aVlist = VerLver.ChangeFromKey(aProVertex);
        aVlist.Append (aVertex.Oriented(TopAbs_FORWARD));
      }
      else
      {
        TopTools_ListOfShape aVlist;
        aVlist.Append (aVertex.Oriented(TopAbs_FORWARD));
        VerLver.Add (aProVertex,  aVlist);
      }
    }
  }

  if (VerLver.IsEmpty())
    return;

  BRep_Builder aBB;
  for (Standard_Integer ii = 1; ii <= VerLver.Extent(); ii++)
  {
    // In some cases (concave faces), a vertex may have images in more than one
    // distinct location. So we can't always assume all images to converge into
    // one vertex.
    const TopTools_ListOfShape& aVertexList = VerLver(ii);
    if (aVertexList.Extent() == 1)
      continue;
    
    SHOW_TOPO_SHAPE(TopoDS_Shape(), "VEMap", aVertexList);

    TopTools_ListOfShape aVlist = aVertexList;
    TopTools_ListOfShape OutLiers;

    gp_Pnt aCentre;
    Standard_Real Tol2 = myTolConf * myTolConf;
    TColgp_Array1OfPnt Points (1, aVlist.Extent());
    TopTools_ListIteratorOfListOfShape itl;

    auto GetCenter = [&](const TopTools_ListOfShape &Vertices) {

      Standard_Real aMaxTol = 0.;

      Standard_Integer Count = 0;
      for (itl.Initialize(Vertices); itl.More(); ) {
        const TopoDS_Vertex& aVertex = TopoDS::Vertex (itl.Value());
        Standard_Real aTol = BRep_Tool::Tolerance(aVertex);
        aMaxTol = Max (aMaxTol, aTol);
        gp_Pnt aPnt = BRep_Tool::Pnt(aVertex);
        if (Count == 0 || Points(1).SquareDistance(aPnt) < Tol2) {
          Points(++Count) = aPnt;
          itl.Next();
        }
        else {
          OutLiers.Append(aVertex);
          aVlist.Remove(itl);
        }
      }

      gp_Ax2 anAxis;
      Standard_Boolean IsSingular;
      GeomLib::AxeOfInertia (Points, anAxis, IsSingular);
      aCentre = anAxis.Location();
      Standard_Real aMaxDist = 0.;
      for (Standard_Integer jj = 1; jj <= Count; jj++)
      {
        Standard_Real aSqDist = aCentre.SquareDistance (Points(jj));
        aMaxDist = Max (aMaxDist, aSqDist);
      }
      aMaxDist = Sqrt(aMaxDist);
      return Max (aMaxTol, aMaxDist);
    };


    while (true) {
      Standard_Real aMaxTol = GetCenter(aVlist);

      //Find constant vertex
      TopoDS_Vertex aConstVertex;
      for (itl.Initialize(aVlist); itl.More(); itl.Next())
      {
        const TopoDS_Vertex& aVertex = TopoDS::Vertex (itl.Value());
        const TopTools_ListOfShape& aElist = theVEmap.FindFromKey(aVertex);
        const TopoDS_Shape& anEdge = aElist.First();
        TopTools_ListIteratorOfListOfShape itcedges (myConstEdges);
        for (; itcedges.More(); itcedges.Next())
          if (anEdge.IsSame (itcedges.Value()))
          {
            aConstVertex = aVertex;
            SHOW_TOPO_SHAPE(aConstVertex, "ConstV_");
            break;
          }
        if (!aConstVertex.IsNull())
          break;
      }
      if (aConstVertex.IsNull())
      {
        aConstVertex = TopoDS::Vertex(aVlist.First());
        SHOW_TOPO_SHAPE(aConstVertex, "ConstVN");
      }
      aBB.UpdateVertex (aConstVertex, aCentre, aMaxTol);

      for (itl.Initialize(aVlist); itl.More(); itl.Next())
      {
        const TopoDS_Vertex& aVertex = TopoDS::Vertex (itl.Value());
        if (aVertex.IsSame(aConstVertex))
          continue;
        
        const TopTools_ListOfShape& aElist = theVEmap.FindFromKey (aVertex);
        for (TopTools_ListIteratorOfListOfShape itl1(aElist); itl1.More(); itl1.Next())
        {
          TopoDS_Edge anEdge = TopoDS::Edge (itl1.Value());
          SHOW_TOPO_SHAPE(anEdge, "ReplaceConstV");
          anEdge.Orientation(TopAbs_FORWARD);
          TopoDS_Vertex aV1, aV2;
          TopExp::Vertices (anEdge, aV1, aV2);
          TopoDS_Vertex aVertexToRemove = (aV1.IsSame(aVertex))? aV1 : aV2;
          anEdge.Free(Standard_True);
          aBB.Remove (anEdge, aVertexToRemove);
          aBB.Add (anEdge, aConstVertex.Oriented (aVertexToRemove.Orientation()));
        }
      }

      if (OutLiers.Extent() <= 1)
        break;
      aVlist.Clear();
      aVlist.Append(OutLiers);
    }
  }

  TopTools_IndexedMapOfShape Emap;
  for (Standard_Integer ii = 1; ii <= theVEmap.Extent(); ii++)
  {
    const TopTools_ListOfShape& aElist = theVEmap(ii);
    TopTools_ListIteratorOfListOfShape itl (aElist);
    for (; itl.More(); itl.Next())
      Emap.Add (itl.Value());
  }

  theVEmap.Clear();
  for (Standard_Integer ii = 1; ii <= Emap.Extent(); ii++)
    TopExp::MapShapesAndAncestors (Emap(ii), TopAbs_VERTEX, TopAbs_EDGE, theVEmap);
}
