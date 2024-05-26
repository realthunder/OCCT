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
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
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

#define OCCT_DEBUG_ALGO

#ifdef OCCT_DEBUG_ALGO
Standard_Boolean AffichLoop  = Standard_True;
Standard_Integer NbLoops     = 0;
Standard_Integer NbWires     = 1;
static char* name = new char[100];
#endif

extern "C" {
#if 1
void showTopoShape(const TopoDS_Shape &s, const char *name);
void showTopoShapes(const TopoDS_Shape &s, const char *name, const TopTools_ListOfShape &shapes);
#else
static void showTopoShape(const TopoDS_Shape &s, const char *name)
{
}
static void showTopoShapes(const TopoDS_Shape &s, const char *name, const TopTools_ListOfShape &shapes)
{
}
#endif
}
// #define DRAW

namespace DBRep {
static void Set(const char *name, const TopoDS_Shape &shape) {
    showTopoShape(shape, name);
}
}

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
		   TopTools_SequenceOfShape& Seq) 
{
  //Remove duplicates
  for (Standard_Integer i = 1; i < Seq.Length(); i++)
    for (Standard_Integer j = i+1; j <= Seq.Length(); j++)
      if (Seq(i) == Seq(j))
      {
        Seq.Remove(j);
        j--;
      }
  
  Standard_Boolean Invert   = Standard_True;
  Standard_Integer NbPoints = Seq.Length();
  Standard_Real    U1,U2;
  TopoDS_Vertex    V1,V2;

  while (Invert) {
    Invert = Standard_False;
    for ( Standard_Integer i = 1; i < NbPoints; i++) {
      TopoDS_Shape aLocalV = Seq.Value(i)  .Oriented(TopAbs_INTERNAL);
      V1 = TopoDS::Vertex(aLocalV);
      aLocalV = Seq.Value(i+1).Oriented(TopAbs_INTERNAL);
      V2 = TopoDS::Vertex(aLocalV);
//      V1 = TopoDS::Vertex(Seq.Value(i)  .Oriented(TopAbs_INTERNAL));
//      V2 = TopoDS::Vertex(Seq.Value(i+1).Oriented(TopAbs_INTERNAL));

      U1 = BRep_Tool::Parameter(V1,E);
      U2 = BRep_Tool::Parameter(V2,E);
      if (U2 < U1) {
	Seq.Exchange(i,i+1);
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
  showTopoShapes(E, "AddEdge", LV);
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
//function : SamePnt2d
//purpose  : 
//=======================================================================

static Standard_Boolean  SamePnt2d(const TopoDS_Vertex&  V,
				   TopoDS_Edge&   E1,
				   TopoDS_Edge&   E2,
				   const TopoDS_Face&   F)
{
  Standard_Real   f1,f2,l1,l2;
  gp_Pnt2d        P1,P2;
  TopoDS_Shape aLocalF = F.Oriented(TopAbs_FORWARD);
  TopoDS_Face FF = TopoDS::Face(aLocalF);
//  TopoDS_Face FF = TopoDS::Face(F.Oriented(TopAbs_FORWARD));
  Handle(Geom2d_Curve) C1 = BRep_Tool::CurveOnSurface(E1,FF,f1,l1);  
  Handle(Geom2d_Curve) C2 = BRep_Tool::CurveOnSurface(E2,FF,f2,l2);  
  if (E1.Orientation () == TopAbs_FORWARD) P1 = C1->Value(f1);
  else                                     P1 = C1->Value(l1);
  
  if (E2.Orientation () == TopAbs_FORWARD) P2 = C2->Value(l2);
  else                                     P2 = C2->Value(f2);
  Standard_Real Tol  = 100*BRep_Tool::Tolerance(V);
  Standard_Real Dist = P1.Distance(P2);
  if (Dist < Tol)
    return Standard_True;

//  TopoDS_Face FF = TopoDS::Face(F.Oriented(TopAbs_FORWARD));
  gp_Pnt2d P11 = C1->Value(f1);
  gp_Pnt2d P12 = C1->Value(l1);
  gp_Pnt2d P21 = C2->Value(f2);
  gp_Pnt2d P22 = C2->Value(l2);
  Standard_Boolean Result = Standard_False;
  Standard_Real D;
  if ((D = P11.Distance(P21)) < Tol) {
    Result = Standard_True;
    char name[128];
    snprintf(name, sizeof(name), "Same1121_%g_", D);
    // showTopoShape(V, name);
  }
  else if ((D = P12.Distance(P21)) < Tol) {
    char name[128];
    snprintf(name, sizeof(name), "Same1221_%g_", D);
    // showTopoShape(V, name);
    Result = Standard_True;
  }
  else if ((D = P11.Distance(P22)) < Tol) {
    char name[128];
    snprintf(name, sizeof(name), "Same1122_%g_", D);
    // showTopoShape(V, name);
    Result = Standard_True;
  }
  else if ((D = P12.Distance(P22)) < Tol) {
    char name[128];
    snprintf(name, sizeof(name), "Same1221_%g_", D);
    // showTopoShape(V, name);
    Result = Standard_True;
  }

  if (Result) {
    // showTopoShape(E1, "SameE1_");
    // showTopoShape(E2, "SameE2_");
  }

  return Result;
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
              showTopoShape(E, "ReplaceMVEE");
              showTopoShape(V, "ReplaceMVEV");
	    }
	}
    }

  TopExp::Vertices(E,V1,V2);
  if( V1.IsNull() && V2.IsNull() ){ YaCouture = Standard_False; return; }
  if (!MVE.Contains(V1)) {
    MVE.Add(V1,Empty);
  }
  MVE.ChangeFromKey(V1).Append(E);
  // showTopoShape(V1, "MVE_V1_");
  if (!V1.IsSame(V2)) {
     if (!MVE.Contains(V2)) {
       MVE.Add(V2,Empty);
     }
     MVE.ChangeFromKey(V2).Append(E);
     // showTopoShape(V2, "MVE_V2_");
  }
  TopLoc_Location L ;
  Handle(Geom_Surface) S = BRep_Tool::Surface(F,L);
  if (BRep_Tool::IsClosed(E,S,L)) {
    MVE.ChangeFromKey(V2).Append(E.Reversed());
    // showTopoShape(V2, "MVE_ClosedV2_");
    if (!V1.IsSame(V2)) {
      MVE.ChangeFromKey(V1).Append(E.Reversed());
      // showTopoShape(V1, "MVE_ClosedV1_");
    }
    YaCouture = Standard_True;
  }
  showTopoShape(E, "MVE");
}

//=======================================================================
//function : Perform
//purpose  : 
//=======================================================================

void BRepAlgo_Loop::Perform()
{
  TopTools_ListIteratorOfListOfShape                  itl, itl1;
  TopoDS_Vertex                                       V1,V2;
  BRep_Builder                                        B;

  showTopoShape(myFace, "LoopFace");

  //------------------------------------------------
  // Check intersection in myConstEdges which is possible when make thick solid
  // with concave removed face
  //------------------------------------------------
  TopTools_ListOfShape ConstEdges;
  TopTools_ListOfShape IntersectingEdges;
  if (_CollectingEdges) {
    TopTools_MapOfShape EMap;
    for (itl.Initialize(myConstEdges); itl.More(); itl.Next())
    {
      const TopoDS_Edge& anEdge = TopoDS::Edge(itl.Value());
      // Sewn edges can be doubled or not in myConstEdges
      if (!EMap.Add(anEdge))
        continue;

      Standard_Real aF, aL;
      const Handle(Geom_Curve) C = BRep_Tool::Curve(anEdge, aF, aL);
      TopExp::Vertices(anEdge, V1, V2);
      TopTools_ListOfShape LV;
      for (itl1.Initialize(myConstEdges); itl1.More(); itl1.Next()) {
        const TopoDS_Edge& otherEdge = TopoDS::Edge(itl1.Value());
        if (otherEdge.IsSame(anEdge))
          continue;
        TopoDS_Vertex OV[2];
        TopExp::Vertices (otherEdge,OV[0],OV[1]);
        for (Standard_Integer i = 0; i < 2; i++) {
          Standard_Real Tol = BRep_Tool::Tolerance(OV[i]);
          gp_Pnt OP = BRep_Tool::Pnt(OV[i]);
          if (OP.Distance(BRep_Tool::Pnt(V1)) < Tol
              || OP.Distance(BRep_Tool::Pnt(V2)) < Tol)
            continue;
          GeomAPI_ProjectPointOnCurve Proj(BRep_Tool::Pnt(OV[i]), C);
          if (Proj.NbPoints() > 0) {
            Standard_Real D = Proj.LowerDistance();
            Standard_Real P = Proj.LowerDistanceParameter();
            if (D < Tol  && P > aF && P < aL) {
              TopoDS_Shape aLocalShape = OV[i].Oriented(TopAbs_INTERNAL);
              B.UpdateVertex(TopoDS::Vertex(aLocalShape),P,anEdge,Tol);
              showTopoShape(aLocalShape, "InterVertex");
              showTopoShape(anEdge, "InterEdge");
              LV.Append(aLocalShape);
            }
          }
        }
      }
      if (!LV.IsEmpty()) {
        IntersectingEdges.Append(anEdge);
        myVerOnEdges.Bind(anEdge, LV);
      }
      else {
        ConstEdges.Append(anEdge);
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
  for (itl.Initialize(myEdges); itl.More(); itl.Next())
  {
    const TopoDS_Edge& anEdge = TopoDS::Edge(itl.Value());
    if (myCutEdges.Seek(anEdge))
      continue;
    TopTools_ListOfShape LCE;
    const TopTools_ListOfShape* pVertices = myVerOnEdges.Seek (anEdge);
    if (pVertices)
    {
      CutEdge (anEdge, *pVertices, LCE);
      myCutEdges.Bind(anEdge, LCE);
    }
  }

  for (itl.Initialize(IntersectingEdges); itl.More(); itl.Next())
  {
    const TopoDS_Edge& anEdge = TopoDS::Edge(itl.Value());
    if (myCutEdges.Seek(anEdge))
      continue;
    TopTools_ListOfShape LCE;
    const TopTools_ListOfShape* pVertices = myVerOnEdges.Seek (anEdge);
    if (pVertices)
    {
      CutEdge (anEdge, *pVertices, LCE, Standard_True/*Keep all cut edge*/);
      myCutEdges.Bind(anEdge, LCE);
    }

  }

  FindLoop();

  for (itl.Initialize(IntersectingEdges); itl.More(); itl.Next()) {
    const TopTools_ListOfShape &aList = NewEdges(TopoDS::Edge(itl.Value()));
    _EdgeMap.Bind(itl.Value(), aList);
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
                  TopTools_IndexedDataMapOfShapeListOfShape& MVE,
                  MapOfWire& NewWires,
                  const TopoDS_Face& aFace,
                  const Standard_Real &Tol)
{
  TopTools_ListIteratorOfListOfShape itl;
  TopoDS_Vertex                      V1, V2, NV;

  CurrentVEMap.Bind(CV, CE);
  CurrentEdgeList.Prepend(CE);
  // showTopoShape(CV, "CV");
  showTopoShape(CE, "CE");

  TopExp::Vertices(CE, V1, V2);
  if (CV.IsSame(V1))
    NV = V2;
  else
    NV = V1;

  const TopoDS_Shape *pE = CurrentVEMap.Seek(NV);
  if (pE && !pE->IsSame(CE))
  {
    TopoDS_Edge EF = TopoDS::Edge(*pE);
    TopoDS_Edge E = EF;
#if 1
    TopoDS_Vertex VF, VL;
    BRepLib_MakeWire aMakeWire;
    for (;;)
    {
      TopExp::Vertices(E, V1, V2, Standard_False);
      // showTopoShape(V, "NWV");
      // showTopoShape(E, "NWE");
      if (VF.IsNull())
      {
        VF = NV;
        VL = NV.IsSame(V1) ? V2 : V1;
      }
      else if (VL.IsSame(V1))
        VL = V2;
      else
        VL = V1;
      aMakeWire.Add(E);
      if (VL.IsSame(VF))
        break;
      E = TopoDS::Edge(CurrentVEMap.Find(VL));
    }
    TopoDS_Wire NW = aMakeWire.Wire();
    if (NW.Closed() && NewWires.Add(NW))
    {
        showTopoShape(NW, "NewWire");
    }
#else
    TopoDS_Wire NW;
    BRep_Builder::MakeWire(NW);
    TopoDS_Vertex V = NV;
    for (;;)
    {
      TopExp::Vertices(E, V1, V2, Standard_False);
      // showTopoShape(V, "NWV");
      // showTopoShape(E, "NWE");
      if (V.IsSame(V1))
      {
        V = V2;
      }
      else
      {
        E = TopoDS::Edge(E.Reversed());
        V = V1;
      }
      B.Add(NW, E);
      if (!V.IsSame(NV))
      {
        E = TopoDS::Edge(CurrentVEMap.Find(V));
      }
      else if (SamePnt2d(V, EF, E, aFace))
      {
        NW.Closed(Standard_True);
        if (NewWires.Add(NW))
        {
          showTopoShape(NW, "NewWire");
        }
        break;
      }
      else if(BRep_Tool::Tolerance(V) < Tol)
      {
        BRep_Builder aBB;
        aBB.UpdateVertex(V, Tol);
        if (SamePnt2d(V, EF, E, aFace))
        {
          NW.Closed(Standard_True);
          if (NewWires.Add(NW))
            showTopoShape(NW, "NewWire1_");
        }
        break;
      }
      else
        break;
    }
#endif
  }
  else
  {
    for (itl.Initialize(MVE.FindFromKey(NV)); itl.More(); itl.Next())
    {
      const TopoDS_Edge& NE = TopoDS::Edge(itl.Value());
      if (!NE.IsSame(CE) && (!pE || !pE->IsSame(NE)))
        FindAllLoops(NV, NE, CurrentVEMap, CurrentEdgeList, MVE, NewWires, aFace, Tol);
    }
  }

  CurrentVEMap.UnBind(CV);
  CurrentEdgeList.RemoveFirst();
  // showTopoShapes(CE, "Pop", CurrentEdgeList);
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
      showTopoShape(aNF, "MakeFace");

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
      // the wire, we if there is any edge sharing this vertex has its middle
      // point inside the wire.
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

        showTopoShape(aEdge, "CheckSplit");

        // Get 2d curve of the edge on the face
        Standard_Real aT1, aT2;
        const Handle(Geom2d_Curve)& aC2D = BRep_Tool::CurveOnSurface(aEdge, aNF, aT1, aT2);
        if (aC2D.IsNull()) {
          showTopoShape(aEdge, "Prune_NoCurve_");
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
          showTopoShape(aNF, "Prune2_");
          aPrunedMap.Add(aWire);
          break;
        }
        else if (!aClassifier.IsHole() && aState == TopAbs_IN)
        {
          showTopoShape(aNF, "Prune3_");
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
      showTopoShape(aWire, "Prune1_");
      TopTools_MapIteratorOfMapOfShape itM1(aCheckMap);
      for (; itM1.More(); itM1.Next())
        showTopoShape(itM1.Value(), "Prune1W_");
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
  CutEdge(E, VOnE, NE, Standard_False);
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
  TopTools_ListIteratorOfListOfShape it(VOnE);
  BRep_Builder                       B;

  for ( ; it.More(); it.Next()) {
    SV.Append(it.Value());
  }
  //--------------------------------
  // Parse vertices on the edge.
  //--------------------------------
  Bubble (WE,SV);

  if (KeepAll)
    showTopoShapes(WE, "CuttingInter", VOnE);
  else {
    showTopoShapes(WE, "Cutting", VOnE);
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

  TopoDS_Vertex    VF,VL;
  Standard_Real    f,l;
  BRep_Tool::Range(WE,f,l);
  TopExp::Vertices(WE,VF,VL);

  TopoDS_Iterator It(WE);
  Standard_Boolean Extended = Standard_False;
  for (; It.More(); It.Next())
  {
    if (It.Value().Orientation() == TopAbs_INTERNAL)
    {
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
    }
    if (!Extended) {
      SV.Prepend(VF);
      SV.Append(VL);
    }
  }
  else if (!Extended) {
    //-----------------------------------------
    // Eventually all extremities of the edge.
    //-----------------------------------------
    if (!VF.IsNull() && !VF.IsSame(SV.First())) SV.Prepend(VF);
    if (!VL.IsNull() && !VL.IsSame(SV.Last ())) SV.Append (VL);
  }

  if (!KeepAll && Extended)
  {
    SV.ChangeFirst().Orientation(TopAbs_FORWARD);
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
    if (!SV.First().IsSame(V1)) {
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
      showTopoShape(NewEdge, "CutEdge");
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
        showTopoShape(EE, "CutEdgeRemove");
	NE.Remove(it);
      }
      else
	{
	  gp_Pnt2d pf, pl;
	  BRep_Tool::UVPoints( EE, myFace, pf, pl );
	  if (pf.Distance(pl) <= Tol && !BRep_Tool::IsClosed(EE)) {
            showTopoShape(EE, "CutEdgeRemove");
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
    showTopoShapes(myFace, "NewFace", myNewFaces);
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
    const TopTools_ListOfShape& aVlist = VerLver(ii);
    if (aVlist.Extent() == 1)
      continue;
    
    Standard_Real aMaxTol = 0.;
    TColgp_Array1OfPnt Points (1, aVlist.Extent());

    TopTools_ListIteratorOfListOfShape itl (aVlist);
    Standard_Integer jj = 0;
    for (; itl.More(); itl.Next())
    {
      const TopoDS_Vertex& aVertex = TopoDS::Vertex (itl.Value());
      Standard_Real aTol = BRep_Tool::Tolerance(aVertex);
      aMaxTol = Max (aMaxTol, aTol);
      gp_Pnt aPnt = BRep_Tool::Pnt(aVertex);
      Points(++jj) = aPnt;
    }

    gp_Ax2 anAxis;
    Standard_Boolean IsSingular;
    GeomLib::AxeOfInertia (Points, anAxis, IsSingular);
    gp_Pnt aCentre = anAxis.Location();
    Standard_Real aMaxDist = 0.;
    for (jj = 1; jj <= Points.Upper(); jj++)
    {
      Standard_Real aSqDist = aCentre.SquareDistance (Points(jj));
      aMaxDist = Max (aMaxDist, aSqDist);
    }
    aMaxDist = Sqrt(aMaxDist);
    aMaxTol = Max (aMaxTol, aMaxDist);

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
          showTopoShape(aConstVertex, "ConstV_");
          break;
        }
      if (!aConstVertex.IsNull())
        break;
    }
    if (aConstVertex.IsNull())
    {
      aConstVertex = TopoDS::Vertex(aVlist.First());
      showTopoShape(aConstVertex, "ConstVN");
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
        showTopoShape(anEdge, "ReplaceConstV");
        anEdge.Orientation(TopAbs_FORWARD);
        TopoDS_Vertex aV1, aV2;
        TopExp::Vertices (anEdge, aV1, aV2);
        TopoDS_Vertex aVertexToRemove = (aV1.IsSame(aVertex))? aV1 : aV2;
        anEdge.Free(Standard_True);
        aBB.Remove (anEdge, aVertexToRemove);
        aBB.Add (anEdge, aConstVertex.Oriented (aVertexToRemove.Orientation()));
      }
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
