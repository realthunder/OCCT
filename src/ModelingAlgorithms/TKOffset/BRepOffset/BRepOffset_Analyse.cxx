// Created on: 1995-10-20
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

#include <Adaptor3d_Surface.hxx>
#include <BOPTools_AlgoTools.hxx>
#include <BOPTools_AlgoTools3D.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepLib_MakeEdge.hxx>
#include <BRepLib_MakeFace.hxx>
#include <BRepLib_MakeWire.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Precision.hxx>
#include <BRepOffset_Analyse.hxx>
#include <BRepOffset_Interval.hxx>
#include <BRepOffset_Tool.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepTools.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom_Curve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomLProp_SLProps.hxx>
#include <gp.hxx>
#include <gp_Ax3.hxx>
#include <gp_Cone.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Sphere.hxx>
#include <gp_Vec.hxx>
#include <IntTools_Context.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <NCollection_Map.hxx>
#include <ChFi3d.hxx>
#include <LocalAnalysis_SurfaceContinuity.hxx>

static void CorrectOrientationOfTangent(gp_Vec&              TangVec,
                                        const TopoDS_Vertex& aVertex,
                                        const TopoDS_Edge&   anEdge)
{
  TopoDS_Vertex Vlast = TopExp::LastVertex(anEdge);
  if (aVertex.IsSame(Vlast))
  {
    TangVec.Reverse();
  }
}

static bool CheckMixedContinuity(const TopoDS_Edge& theEdge,
                                 const TopoDS_Face& theFace1,
                                 const TopoDS_Face& theFace2,
                                 const double       theAngTol);

//=================================================================================================

BRepOffset_Analyse::BRepOffset_Analyse()
    : myOffset(0.0),
      myDone(false)
{
}

//=================================================================================================

BRepOffset_Analyse::BRepOffset_Analyse(const TopoDS_Shape& S, const double Angle)
    : myOffset(0.0),
      myDone(false)
{
  Perform(S, Angle);
}

//=================================================================================================

static void EdgeAnalyse(const TopoDS_Edge&                     E,
                        const TopoDS_Face&                     F1,
                        const TopoDS_Face&                     F2,
                        const double                           SinTol,
                        NCollection_List<BRepOffset_Interval>& LI)
{
  double f, l;
  BRep_Tool::Range(E, F1, f, l);
  BRepOffset_Interval I;
  I.First(f);
  I.Last(l);
  //
  BRepAdaptor_Surface aBAsurf1(F1, false);
  GeomAbs_SurfaceType aSurfType1 = aBAsurf1.GetType();

  BRepAdaptor_Surface aBAsurf2(F2, false);
  GeomAbs_SurfaceType aSurfType2 = aBAsurf2.GetType();

  bool isTwoPlanes = (aSurfType1 == GeomAbs_Plane && aSurfType2 == GeomAbs_Plane);

  ChFiDS_TypeOfConcavity ConnectType = ChFiDS_Other;

  if (isTwoPlanes) // then use only strong condition
  {
    if (BRep_Tool::Continuity(E, F1, F2) > GeomAbs_C0)
    {
      ConnectType = ChFiDS_Tangential;
    }
    else
    {
      ConnectType = ChFi3d::DefineConnectType(E, F1, F2, SinTol, false);
    }
  }
  else
  {
    bool isTwoSplines =
      (aSurfType1 == GeomAbs_BSplineSurface || aSurfType1 == GeomAbs_BezierSurface)
      && (aSurfType2 == GeomAbs_BSplineSurface || aSurfType2 == GeomAbs_BezierSurface);
    bool isMixedConcavity = false;
    if (isTwoSplines)
    {
      double anAngTol  = 0.1;
      isMixedConcavity = CheckMixedContinuity(E, F1, F2, anAngTol);
    }

    if (!isMixedConcavity)
    {
      if (ChFi3d::IsTangentFaces(E, F1, F2)) // weak condition
      {
        ConnectType = ChFiDS_Tangential;
      }
      else
      {
        ConnectType = ChFi3d::DefineConnectType(E, F1, F2, SinTol, false);
      }
    }
    else
    {
      ConnectType = ChFiDS_Mixed;
    }
  }

  I.Type(ConnectType);
  LI.Append(I);
}

//=================================================================================================

bool CheckMixedContinuity(const TopoDS_Edge& theEdge,
                          const TopoDS_Face& theFace1,
                          const TopoDS_Face& theFace2,
                          const double       theAngTol)
{
  bool          aMixedCont = false;
  GeomAbs_Shape aCurrOrder = BRep_Tool::Continuity(theEdge, theFace1, theFace2);
  if (aCurrOrder > GeomAbs_C0)
  {
    // Method BRep_Tool::Continuity(...) always returns minimal continuity between faces
    // so, if aCurrOrder > C0 it means that faces are tangent along whole edge.
    return aMixedCont;
  }
  // But we cannot trust result, if it is C0, because this value is set by default.
  double TolC0 = std::max(0.001, 1.5 * BRep_Tool::Tolerance(theEdge));

  double aFirst;
  double aLast;

  occ::handle<Geom2d_Curve> aC2d1, aC2d2;

  if (!theFace1.IsSame(theFace2) && BRep_Tool::IsClosed(theEdge, theFace1)
      && BRep_Tool::IsClosed(theEdge, theFace2))
  {
    // Find the edge in the face 1: this edge will have correct orientation
    TopoDS_Edge anEdgeInFace1;
    TopoDS_Face aFace1 = theFace1;
    aFace1.Orientation(TopAbs_FORWARD);
    TopExp_Explorer anExplo(aFace1, TopAbs_EDGE);
    for (; anExplo.More(); anExplo.Next())
    {
      const TopoDS_Edge& anEdge = TopoDS::Edge(anExplo.Current());
      if (anEdge.IsSame(theEdge))
      {
        anEdgeInFace1 = anEdge;
        break;
      }
    }
    if (anEdgeInFace1.IsNull())
    {
      return aMixedCont;
    }

    aC2d1              = BRep_Tool::CurveOnSurface(anEdgeInFace1, aFace1, aFirst, aLast);
    TopoDS_Face aFace2 = theFace2;
    aFace2.Orientation(TopAbs_FORWARD);
    anEdgeInFace1.Reverse();
    aC2d2 = BRep_Tool::CurveOnSurface(anEdgeInFace1, aFace2, aFirst, aLast);
  }
  else
  {
    // Obtaining of pcurves of edge on two faces.
    aC2d1 = BRep_Tool::CurveOnSurface(theEdge, theFace1, aFirst, aLast);
    // For the case of seam edge
    TopoDS_Edge EE = theEdge;
    if (theFace1.IsSame(theFace2))
    {
      EE.Reverse();
    }
    aC2d2 = BRep_Tool::CurveOnSurface(EE, theFace2, aFirst, aLast);
  }

  if (aC2d1.IsNull() || aC2d2.IsNull())
  {
    return aMixedCont;
  }

  // Obtaining of two surfaces from adjacent faces.
  occ::handle<Geom_Surface> aSurf1 = BRep_Tool::Surface(theFace1);
  occ::handle<Geom_Surface> aSurf2 = BRep_Tool::Surface(theFace2);

  if (aSurf1.IsNull() || aSurf2.IsNull())
  {
    return aMixedCont;
  }

  int aNbSamples = 23;

  // Check for mixed concavity: convex in some regions, concave in others.
  const double aDelta      = (aLast - aFirst) / (aNbSamples - 1);
  bool         aHasConvex  = false;
  bool         aHasConcave = false;
  int          aNbValid    = 0;

  for (int i = 1; i <= aNbSamples; i++)
  {
    const double aPar = (i == aNbSamples) ? aLast : aFirst + (i - 1) * aDelta;

    LocalAnalysis_SurfaceContinuity aCont(aC2d1,
                                          aC2d2,
                                          aPar,
                                          aSurf1,
                                          aSurf2,
                                          GeomAbs_G1,
                                          0.001,
                                          TolC0,
                                          theAngTol,
                                          theAngTol,
                                          theAngTol);
    if (!aCont.IsDone())
    {
      continue;
    }

    aNbValid++;

    if (!aCont.IsG1() && (!aHasConvex || !aHasConcave))
    {
      const double anAngle = aCont.C0Value();
      aHasConvex           = aHasConvex || (anAngle > M_PI_2 + theAngTol);
      aHasConcave          = aHasConcave || (anAngle < M_PI_2 - theAngTol);
    }
  }

  if (aNbValid < aNbSamples / 2)
  {
    return aMixedCont;
  }

  // Mixed connectivity: both convex and concave regions exist.
  aMixedCont = aHasConvex && aHasConcave;

  return aMixedCont;
}

//=================================================================================================

static void BuildAncestors(
  const TopoDS_Shape& S,
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>&
    MA)
{
  MA.Clear();
  TopExp::MapShapesAndUniqueAncestors(S, TopAbs_VERTEX, TopAbs_EDGE, MA);
  TopExp::MapShapesAndUniqueAncestors(S, TopAbs_EDGE, TopAbs_FACE, MA);
}

//=================================================================================================

void BRepOffset_Analyse::Perform(const TopoDS_Shape&          S,
                                 const double                 Angle,
                                 const Message_ProgressRange& theRange)
{
  myShape = S;
  myNewFaces.Clear();
  myGenerated.Clear();
  myReplacement.Clear();
  myDescendants.Clear();

  myAngle       = Angle;
  double SinTol = std::abs(std::sin(Angle));

  // Build ancestors.
  BuildAncestors(S, myAncestors);

  NCollection_List<TopoDS_Shape> aLETang;
  TopExp_Explorer                Exp(S.Oriented(TopAbs_FORWARD), TopAbs_EDGE);
  Message_ProgressScope          aPSOuter(theRange, nullptr, 2);
  Message_ProgressScope          aPS(aPSOuter.Next(), "Performing edges analysis", 1, true);
  for (; Exp.More(); Exp.Next(), aPS.Next())
  {
    if (!aPS.More())
    {
      return;
    }
    const TopoDS_Edge& E = TopoDS::Edge(Exp.Current());
    if (!myMapEdgeType.IsBound(E))
    {
      NCollection_List<BRepOffset_Interval> LI;
      myMapEdgeType.Bind(E, LI);

      const NCollection_List<TopoDS_Shape>& L = Ancestors(E);
      if (L.IsEmpty())
      {
        continue;
      }

      if (L.Extent() == 2)
      {
        const TopoDS_Face& F1 = TopoDS::Face(L.First());
        const TopoDS_Face& F2 = TopoDS::Face(L.Last());
        EdgeAnalyse(E, F1, F2, SinTol, myMapEdgeType(E));

        // For tangent faces add artificial perpendicular face
        // to close the gap between them (if they have different offset values)
        if (myMapEdgeType(E).Last().Type() == ChFiDS_Tangential)
        {
          aLETang.Append(E);
        }
      }
      else if (L.Extent() == 1)
      {
        double             U1, U2;
        const TopoDS_Face& F = TopoDS::Face(L.First());
        BRep_Tool::Range(E, F, U1, U2);
        BRepOffset_Interval Inter(U1, U2, ChFiDS_Other);

        if (!BRepTools::IsReallyClosed(E, F))
        {
          Inter.Type(ChFiDS_FreeBound);
        }
        myMapEdgeType(E).Append(Inter);
      }
      else
      {
#ifdef OCCT_DEBUG
        std::cout << "edge shared by more than two faces" << std::endl;
#endif
      }
    }
  }

  TreatTangentFaces(aLETang, aPSOuter.Next());
  if (!aPSOuter.More())
  {
    return;
  }
  myDone = true;
}

//=================================================================================================

void BRepOffset_Analyse::TreatTangentFaces(const NCollection_List<TopoDS_Shape>& theLE,
                                           const Message_ProgressRange&          theRange)
{
  if (theLE.IsEmpty() || myFaceOffsetMap.IsEmpty())
  {
    // Noting to do: either there are no tangent faces in the shape or
    //               the face offset map has not been provided
    return;
  }

  // Select the edges which connect faces with different offset values
  TopoDS_Compound aCETangent;
  BRep_Builder().MakeCompound(aCETangent);
  // Bind to each tangent edge a max offset value of its faces
  NCollection_DataMap<TopoDS_Shape, double, TopTools_ShapeMapHasher> anEdgeOffsetMap;
  // Bind vertices of the tangent edges with connected edges
  // of the face with smaller offset value
  NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher> aDMVEMin;
  Message_ProgressScope aPSOuter(theRange, nullptr, 3);
  Message_ProgressScope aPS1(aPSOuter.Next(),
                             "Binding vertices with connected edges",
                             theLE.Extent());
  for (NCollection_List<TopoDS_Shape>::Iterator it(theLE); it.More(); it.Next(), aPS1.Next())
  {
    if (!aPS1.More())
    {
      return;
    }
    const TopoDS_Shape&                   aE  = it.Value();
    const NCollection_List<TopoDS_Shape>& aLA = Ancestors(aE);

    const TopoDS_Shape &aF1 = aLA.First(), aF2 = aLA.Last();

    const double* pOffsetVal1  = myFaceOffsetMap.Seek(aF1);
    const double* pOffsetVal2  = myFaceOffsetMap.Seek(aF2);
    const double  anOffsetVal1 = pOffsetVal1 ? std::abs(*pOffsetVal1) : myOffset;
    const double  anOffsetVal2 = pOffsetVal2 ? std::abs(*pOffsetVal2) : myOffset;
    if (anOffsetVal1 != anOffsetVal2)
    {
      BRep_Builder().Add(aCETangent, aE);
      anEdgeOffsetMap.Bind(aE, std::max(anOffsetVal1, anOffsetVal2));

      const TopoDS_Shape& aFMin = anOffsetVal1 < anOffsetVal2 ? aF1 : aF2;
      for (TopoDS_Iterator itV(aE); itV.More(); itV.Next())
      {
        const TopoDS_Shape& aV = itV.Value();
        if (Ancestors(aV).Extent() == 3)
        {
          for (TopExp_Explorer expE(aFMin, TopAbs_EDGE); expE.More(); expE.Next())
          {
            const TopoDS_Shape& aEMin = expE.Current();
            if (aEMin.IsSame(aE))
            {
              continue;
            }
            for (TopoDS_Iterator itV1(aEMin); itV1.More(); itV1.Next())
            {
              const TopoDS_Shape& aVx = itV1.Value();
              if (aV.IsSame(aVx))
              {
                aDMVEMin.Bind(aV, aEMin);
              }
            }
          }
        }
      }
    }
  }

  if (anEdgeOffsetMap.IsEmpty())
  {
    return;
  }

  // Create map of Face ancestors for the vertices on tangent edges
  NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    aDMVFAnc;

  Message_ProgressScope aPS2(aPSOuter.Next(), "Creating map of Face ancestors", theLE.Extent());
  for (NCollection_List<TopoDS_Shape>::Iterator itE(theLE); itE.More(); itE.Next(), aPS2.Next())
  {
    if (!aPS2.More())
    {
      return;
    }
    const TopoDS_Shape& aE = itE.Value();
    if (!anEdgeOffsetMap.IsBound(aE))
    {
      continue;
    }

    NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aMFence;
    {
      const NCollection_List<TopoDS_Shape>& aLEA = Ancestors(aE);
      for (NCollection_List<TopoDS_Shape>::Iterator itLEA(aLEA); itLEA.More(); itLEA.Next())
      {
        aMFence.Add(itLEA.Value());
      }
    }

    for (TopoDS_Iterator itV(aE); itV.More(); itV.Next())
    {
      const TopoDS_Shape&             aV   = itV.Value();
      NCollection_List<TopoDS_Shape>* pLFA = aDMVFAnc.Bound(aV, NCollection_List<TopoDS_Shape>());
      const NCollection_List<TopoDS_Shape>& aLVA = Ancestors(aV);
      for (NCollection_List<TopoDS_Shape>::Iterator itLVA(aLVA); itLVA.More(); itLVA.Next())
      {
        const TopoDS_Edge&                           aEA        = TopoDS::Edge(itLVA.Value());
        const NCollection_List<BRepOffset_Interval>* pIntervals = myMapEdgeType.Seek(aEA);
        if (!pIntervals || pIntervals->IsEmpty())
        {
          continue;
        }
        if (pIntervals->First().Type() == ChFiDS_Tangential)
        {
          continue;
        }

        const NCollection_List<TopoDS_Shape>& aLEA = Ancestors(aEA);
        for (NCollection_List<TopoDS_Shape>::Iterator itLEA(aLEA); itLEA.More(); itLEA.Next())
        {
          const TopoDS_Shape& aFA = itLEA.Value();
          if (aMFence.Add(aFA))
          {
            pLFA->Append(aFA);
          }
        }
      }
    }
  }

  occ::handle<IntTools_Context> aCtx = new IntTools_Context();
  // Tangency criteria
  double aSinTol = std::abs(std::sin(myAngle));

  // Make blocks of connected edges
  NCollection_List<NCollection_List<TopoDS_Shape>> aLCB;
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    aMVEMap;

  BOPTools_AlgoTools::MakeConnexityBlocks(aCETangent, TopAbs_VERTEX, TopAbs_EDGE, aLCB, aMVEMap);

  // Analyze each block to find co-planar edges
  Message_ProgressScope aPS3(aPSOuter.Next(),
                             "Analyzing blocks to find co-planar edges",
                             aLCB.Extent());
  for (NCollection_List<NCollection_List<TopoDS_Shape>>::Iterator itLCB(aLCB); itLCB.More();
       itLCB.Next(), aPS3.Next())
  {
    if (!aPS3.More())
    {
      return;
    }
    const NCollection_List<TopoDS_Shape>& aCB = itLCB.Value();

    NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aMFence;
    for (NCollection_List<TopoDS_Shape>::Iterator itCB1(aCB); itCB1.More(); itCB1.Next())
    {
      const TopoDS_Edge& aE1 = TopoDS::Edge(itCB1.Value());
      if (!aMFence.Add(aE1))
      {
        continue;
      }

      TopoDS_Compound aBlock;
      BRep_Builder().MakeCompound(aBlock);
      BRep_Builder().Add(aBlock, aE1.Oriented(TopAbs_FORWARD));

      double                                anOffset = anEdgeOffsetMap.Find(aE1);
      const NCollection_List<TopoDS_Shape>& aLF1     = Ancestors(aE1);

      gp_Dir aDN1;
      BOPTools_AlgoTools3D::GetNormalToFaceOnEdge(aE1, TopoDS::Face(aLF1.First()), aDN1);

      NCollection_List<TopoDS_Shape>::Iterator itCB2 = itCB1;
      for (itCB2.Next(); itCB2.More(); itCB2.Next())
      {
        const TopoDS_Edge& aE2 = TopoDS::Edge(itCB2.Value());
        if (aMFence.Contains(aE2))
        {
          continue;
        }

        const NCollection_List<TopoDS_Shape>& aLF2 = Ancestors(aE2);

        gp_Dir aDN2;
        BOPTools_AlgoTools3D::GetNormalToFaceOnEdge(aE2, TopoDS::Face(aLF2.First()), aDN2);

        if (aDN1.XYZ().Crossed(aDN2.XYZ()).Modulus() < aSinTol)
        {
          BRep_Builder().Add(aBlock, aE2.Oriented(TopAbs_FORWARD));
          aMFence.Add(aE2);
          anOffset = std::max(anOffset, anEdgeOffsetMap.Find(aE2));
        }
      }

      // Make the prism
      BRepPrimAPI_MakePrism aMP(aBlock, gp_Vec(aDN1.XYZ()) * anOffset);
      if (!aMP.IsDone())
      {
        continue;
      }

      NCollection_IndexedDataMap<TopoDS_Shape,
                                 NCollection_List<TopoDS_Shape>,
                                 TopTools_ShapeMapHasher>
        aPrismAncestors;
      TopExp::MapShapesAndAncestors(aMP.Shape(), TopAbs_EDGE, TopAbs_FACE, aPrismAncestors);
      TopExp::MapShapesAndAncestors(aMP.Shape(), TopAbs_VERTEX, TopAbs_EDGE, aPrismAncestors);

      for (TopoDS_Iterator itE(aBlock); itE.More(); itE.Next())
      {
        const TopoDS_Edge&                    aE    = TopoDS::Edge(itE.Value());
        const NCollection_List<TopoDS_Shape>& aLG   = aMP.Generated(aE);
        TopoDS_Face                           aFNew = TopoDS::Face(aLG.First());

        NCollection_List<TopoDS_Shape>& aLA = myAncestors.ChangeFromKey(aE);

        TopoDS_Shape aF1 = aLA.First();
        TopoDS_Shape aF2 = aLA.Last();

        const double* pOffsetVal1  = myFaceOffsetMap.Seek(aF1);
        const double* pOffsetVal2  = myFaceOffsetMap.Seek(aF2);
        const double  anOffsetVal1 = pOffsetVal1 ? std::abs(*pOffsetVal1) : myOffset;
        const double  anOffsetVal2 = pOffsetVal2 ? std::abs(*pOffsetVal2) : myOffset;

        const TopoDS_Shape& aFToRemove = anOffsetVal1 > anOffsetVal2 ? aF1 : aF2;
        const TopoDS_Shape& aFOpposite = anOffsetVal1 > anOffsetVal2 ? aF2 : aF1;

        // Orient the face so its normal is directed to smaller offset face
        {
          // get normal of the new face
          gp_Dir aDN;
          BOPTools_AlgoTools3D::GetNormalToFaceOnEdge(aE, aFNew, aDN);

          // get bi-normal for the aFOpposite
          TopoDS_Edge aEInF;
          for (TopExp_Explorer aExpE(aFOpposite, TopAbs_EDGE); aExpE.More(); aExpE.Next())
          {
            if (aE.IsSame(aExpE.Current()))
            {
              aEInF = TopoDS::Edge(aExpE.Current());
              break;
            }
          }

          gp_Pnt2d                       aP2d;
          gp_Pnt                         aPInF;
          double                         f, l;
          const occ::handle<Geom_Curve>& aC3D  = BRep_Tool::Curve(aEInF, f, l);
          gp_Pnt                         aPOnE = aC3D->Value((f + l) / 2.);
          BOPTools_AlgoTools3D::PointNearEdge(aEInF,
                                              TopoDS::Face(aFOpposite),
                                              (f + l) / 2.,
                                              1.e-5,
                                              aP2d,
                                              aPInF);

          gp_Vec aBN(aPOnE, aPInF);

          if (aBN.Dot(aDN) < 0)
          {
            aFNew.Reverse();
          }
        }

        // Remove the face with bigger offset value from edge ancestors
        for (NCollection_List<TopoDS_Shape>::Iterator itA(aLA); itA.More(); itA.Next())
        {
          if (itA.Value().IsSame(aFToRemove))
          {
            aLA.Remove(itA);
            break;
          }
        }
        aLA.Append(aFNew);

        myMapEdgeType(aE).Clear();
        // Analyze edge again
        EdgeAnalyse(aE, TopoDS::Face(aFOpposite), aFNew, aSinTol, myMapEdgeType(aE));

        // Analyze vertices
        NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aFNewEdgeMap;
        aFNewEdgeMap.Add(aE);
        for (TopoDS_Iterator itV(aE); itV.More(); itV.Next())
        {
          const TopoDS_Shape& aV = itV.Value();
          // Add Side edge to map of Ancestors with the correct orientation
          TopoDS_Edge aEG = TopoDS::Edge(aMP.Generated(aV).First());
          myGenerated.Bind(aV, aEG);
          {
            for (TopExp_Explorer anExpEg(aFNew, TopAbs_EDGE); anExpEg.More(); anExpEg.Next())
            {
              if (anExpEg.Current().IsSame(aEG))
              {
                aEG = TopoDS::Edge(anExpEg.Current());
                break;
              }
            }
          }

          if (aDMVEMin.IsBound(aV))
          {
            const NCollection_List<TopoDS_Shape>* pSA = aDMVFAnc.Seek(aV);
            if (pSA && pSA->Extent() == 1)
            {
              // Adjust orientation of generated edge to its new ancestor
              TopoDS_Edge aEMin = TopoDS::Edge(aDMVEMin.Find(aV));
              for (TopExp_Explorer expEx(pSA->First(), TopAbs_EDGE); expEx.More(); expEx.Next())
              {
                if (expEx.Current().IsSame(aEMin))
                {
                  aEMin = TopoDS::Edge(expEx.Current());
                  break;
                }
              }

              TopAbs_Orientation anOriInEMin(TopAbs_FORWARD), anOriInEG(TopAbs_FORWARD);

              for (TopoDS_Iterator itx(aEMin); itx.More(); itx.Next())
              {
                if (itx.Value().IsSame(aV))
                {
                  anOriInEMin = itx.Value().Orientation();
                  break;
                }
              }

              for (TopoDS_Iterator itx(aEG); itx.More(); itx.Next())
              {
                if (itx.Value().IsSame(aV))
                {
                  anOriInEG = itx.Value().Orientation();
                  break;
                }
              }

              if (anOriInEG == anOriInEMin)
              {
                aEG.Reverse();
              }
            }
          }

          NCollection_List<TopoDS_Shape>& aLVA = myAncestors.ChangeFromKey(aV);
          if (!aLVA.Contains(aEG))
          {
            aLVA.Append(aEG);
          }
          aFNewEdgeMap.Add(aEG);

          NCollection_List<TopoDS_Shape>& aLEGA =
            myAncestors(myAncestors.Add(aEG, aPrismAncestors.FindFromKey(aEG)));
          {
            // Add ancestors from the shape
            const NCollection_List<TopoDS_Shape>* pSA = aDMVFAnc.Seek(aV);
            if (pSA && !pSA->IsEmpty())
            {
              NCollection_List<TopoDS_Shape> aLSA = *pSA;
              aLEGA.Append(aLSA);
            }
          }

          myMapEdgeType.Bind(aEG, NCollection_List<BRepOffset_Interval>());
          if (aLEGA.Extent() == 2)
          {
            EdgeAnalyse(aEG,
                        TopoDS::Face(aLEGA.First()),
                        TopoDS::Face(aLEGA.Last()),
                        aSinTol,
                        myMapEdgeType(aEG));
          }
        }

        // Find an edge opposite to tangential one and add ancestors for it
        TopoDS_Edge aEOpposite;
        for (TopExp_Explorer anExpE(aFNew, TopAbs_EDGE); anExpE.More(); anExpE.Next())
        {
          if (!aFNewEdgeMap.Contains(anExpE.Current()))
          {
            aEOpposite = TopoDS::Edge(anExpE.Current());
            break;
          }
        }

        {
          // Find it in aFOpposite
          for (TopExp_Explorer anExpE(aFToRemove, TopAbs_EDGE); anExpE.More(); anExpE.Next())
          {
            const TopoDS_Shape& aEInFToRem = anExpE.Current();
            if (aE.IsSame(aEInFToRem))
            {
              if (BOPTools_AlgoTools::IsSplitToReverse(aEOpposite, aEInFToRem, aCtx))
              {
                aEOpposite.Reverse();
              }
              break;
            }
          }
        }

        NCollection_List<TopoDS_Shape> aLFOpposite;
        aLFOpposite.Append(aFNew);
        aLFOpposite.Append(aFToRemove);
        myAncestors.Add(aEOpposite, aLFOpposite);
        myMapEdgeType.Bind(aEOpposite, NCollection_List<BRepOffset_Interval>());
        EdgeAnalyse(aEOpposite,
                    aFNew,
                    TopoDS::Face(aFToRemove),
                    aSinTol,
                    myMapEdgeType(aEOpposite));

        NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>* pEEMap =
          myReplacement.ChangeSeek(aFToRemove);
        if (!pEEMap)
        {
          pEEMap = myReplacement.Bound(
            aFToRemove,
            NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>());
        }
        pEEMap->Bind(aE, aEOpposite);

        // Add ancestors for the vertices
        for (TopoDS_Iterator itV(aEOpposite); itV.More(); itV.Next())
        {
          const TopoDS_Shape&                   aV   = itV.Value();
          const NCollection_List<TopoDS_Shape>& aLVA = aPrismAncestors.FindFromKey(aV);
          myAncestors.Add(aV, aLVA);
        }

        myNewFaces.Append(aFNew);
        myGenerated.Bind(aE, aFNew);
      }
    }
  }
}

//=================================================================================================

// The normal of theF at the parameter theT of its edge theE, the way the face
// faces.
static gp_Dir FaceNormalOnEdge(const TopoDS_Edge& theE, const TopoDS_Face& theF, const double theT)
{
  gp_Dir aN;
  BOPTools_AlgoTools3D::GetNormalToFaceOnEdge(theE, theF, theT, aN);
  if (theF.Orientation() == TopAbs_REVERSED)
  {
    aN.Reverse();
  }
  return aN;
}

// A planar face on the closed wire through theEdges, its normal along theDir.
static TopoDS_Face PlanarFaceOn(const NCollection_List<TopoDS_Shape>& theEdges,
                                const gp_Dir&                         theDir)
{
  TopoDS_Face      aF;
  BRepLib_MakeWire aMW;
  for (NCollection_List<TopoDS_Shape>::Iterator it(theEdges); it.More(); it.Next())
  {
    aMW.Add(TopoDS::Edge(it.Value()));
  }
  if (!aMW.IsDone() || !aMW.Wire().Closed())
  {
    return aF;
  }
  BRepLib_MakeFace aMF(aMW.Wire(), true);
  if (!aMF.IsDone())
  {
    return aF;
  }
  aF = aMF.Face();
  BRepAdaptor_Surface aBAS(aF, false);
  if (aBAS.GetType() != GeomAbs_Plane)
  {
    return TopoDS_Face();
  }
  gp_Dir aN = aBAS.Plane().Axis().Direction();
  if (aF.Orientation() == TopAbs_REVERSED)
  {
    aN.Reverse();
  }
  if (aN.Dot(theDir) < 0.)
  {
    aF.Reverse();
  }
  return aF;
}

// Both faces on one sphere: the sphere.
static bool AreOnOneSphere(const TopoDS_Face& theF1, const TopoDS_Face& theF2, gp_Sphere& theSph)
{
  BRepAdaptor_Surface aS1(theF1, false), aS2(theF2, false);
  if (aS1.GetType() != GeomAbs_Sphere || aS2.GetType() != GeomAbs_Sphere)
  {
    return false;
  }
  const gp_Sphere aSph1 = aS1.Sphere(), aSph2 = aS2.Sphere();
  if (aSph1.Location().Distance(aSph2.Location()) > Precision::Confusion()
      || std::abs(aSph1.Radius() - aSph2.Radius()) > Precision::Confusion())
  {
    return false;
  }
  theSph = aSph1;
  return true;
}

// The wall closing a tangent edge between two faces of one sphere. The edge
// is a circle; carried a thickness round the sphere into the cap (theD, at
// the edge's middle) it is the circle of the same axis there, and the wall
// runs from that circle across the thickness on the cone from the sphere's
// centre through it -- square to the sphere, as a plane's wall is square to
// the plane. Its edges: on the cap (theEp), across the thickness at each end
// (theTe, from theVp to theVpn), and where the kept face's offset meets it
// (theEpn). Its normal runs into the cap.
static TopoDS_Face SphereCapWall(const gp_Sphere&         theSph,
                                 const BRepAdaptor_Curve& theC,
                                 const TopoDS_Vertex      theV[2],
                                 const gp_Dir&            theD,
                                 const gp_Dir&            theNrm,
                                 const double             theOffset,
                                 TopoDS_Vertex            theVp[2],
                                 TopoDS_Vertex            theVpn[2],
                                 TopoDS_Edge              theTe[2],
                                 TopoDS_Edge&             theEp,
                                 TopoDS_Edge&             theEpn)
{
  const gp_Pnt aO  = theSph.Location();
  const double aR  = theSph.Radius();
  const gp_Dir aA  = theC.Circle().Axis().Direction();
  const double f   = theC.FirstParameter(), l = theC.LastParameter();
  const gp_Pnt aPm = theC.Value((f + l) / 2.);
  // The polar angle of a point from the axis, and its way off the axis.
  auto aRadial = [&](const gp_Pnt& theP) {
    const gp_Vec aW(aO, theP);
    const gp_Vec aU = aW - gp_Vec(aA) * aW.Dot(gp_Vec(aA));
    return aU.Magnitude() < Precision::Confusion() ? gp_Vec() : aU.Normalized();
  };
  const gp_Vec aWm(aO, aPm);
  const gp_Vec aUm = aRadial(aPm);
  if (aUm.Magnitude() < 0.5)
  {
    return TopoDS_Face();
  }
  const double aTh0 = std::atan2(aUm.Dot(aWm), aWm.Dot(gp_Vec(aA)));
  const gp_Vec aETh = aUm * std::cos(aTh0) - gp_Vec(aA) * std::sin(aTh0);
  const double aS   = gp_Vec(theD).Dot(aETh) > 0. ? 1. : -1.;
  const double aTh1 = aTh0 + aS * std::abs(theOffset) / aR;
  if (aTh1 < 1.e-3 || aTh1 > M_PI - 1.e-3 || std::abs(aTh1 - M_PI / 2.) < 1.e-6)
  {
    return TopoDS_Face();
  }
  // The kept face's offset: off the sphere or into it.
  const double aRn = aR + (gp_Vec(theNrm).Dot(aWm) > 0. ? theOffset : -theOffset);
  if (aRn < Precision::Confusion())
  {
    return TopoDS_Face();
  }
  const gp_Vec aU0 = aRadial(BRep_Tool::Pnt(theV[0])), aU1 = aRadial(BRep_Tool::Pnt(theV[1]));
  if (aU0.Magnitude() < 0.5 || aU1.Magnitude() < 0.5)
  {
    return TopoDS_Face();
  }
  // The cone: its apex at the centre, its generators aTh1 off the axis; v
  // runs along them from the sphere, u round the axis from V's side.
  const bool   isUp  = aTh1 < M_PI / 2.;
  const gp_Dir aZ    = isUp ? aA : aA.Reversed();
  const gp_Ax3 anAx(aO.Translated(gp_Vec(aA) * (aR * std::cos(aTh1))), aZ, gp_Dir(aU0));
  const gp_Cone aCone(anAx, isUp ? aTh1 : M_PI - aTh1, aR * std::sin(aTh1));
  // V's partner round the axis, the way the edge runs.
  double aU1Par = std::atan2(aU0.Crossed(aU1).Dot(gp_Vec(aZ)), aU0.Dot(aU1));
  if (theC.DN(f, 1).Dot(gp_Vec(aZ).Crossed(aU0)) > 0.)
  {
    if (aU1Par <= Precision::Angular())
    {
      aU1Par += 2. * M_PI;
    }
  }
  else if (aU1Par >= -Precision::Angular())
  {
    aU1Par -= 2. * M_PI;
  }
  const double     aVn = aRn - aR;
  BRepLib_MakeFace aMF(aCone,
                       std::min(0., aU1Par),
                       std::max(0., aU1Par),
                       std::min(0., aVn),
                       std::max(0., aVn));
  if (!aMF.IsDone())
  {
    return TopoDS_Face();
  }
  TopoDS_Face aWall = aMF.Face();
  // Its edges, told apart by where their vertices lie.
  const double aTol = 1.e-6 * aR;
  gp_Pnt       aPp[2];
  for (int i = 0; i < 2; ++i)
  {
    const gp_Vec aU = i == 0 ? aU0 : aU1;
    aPp[i] = aO.Translated((aU * std::sin(aTh1) + gp_Vec(aA) * std::cos(aTh1)) * aR);
  }
  for (TopExp_Explorer anExp(aWall, TopAbs_EDGE); anExp.More(); anExp.Next())
  {
    const TopoDS_Edge& anE = TopoDS::Edge(anExp.Current());
    TopoDS_Vertex      aV1, aV2;
    TopExp::Vertices(anE, aV1, aV2);
    const double aD1 = aO.Distance(BRep_Tool::Pnt(aV1)), aD2 = aO.Distance(BRep_Tool::Pnt(aV2));
    const bool   isR1 = std::abs(aD1 - aR) < aTol, isR2 = std::abs(aD2 - aR) < aTol;
    if (isR1 && isR2)
    {
      theEp = anE;
    }
    else if (!isR1 && !isR2)
    {
      theEpn = anE;
    }
    else
    {
      const TopoDS_Vertex& aVR = isR1 ? aV1 : aV2;
      const int i = BRep_Tool::Pnt(aVR).Distance(aPp[0]) < BRep_Tool::Pnt(aVR).Distance(aPp[1]) ? 0 : 1;
      theTe[i]  = anE;
      theVp[i]  = aVR;
      theVpn[i] = isR1 ? aV2 : aV1;
    }
  }
  if (theEp.IsNull() || theEpn.IsNull() || theTe[0].IsNull() || theTe[1].IsNull()
      || theTe[0].IsSame(theTe[1]))
  {
    return TopoDS_Face();
  }
  // Into the cap.
  BRepAdaptor_Surface aBAS(aWall, false);
  gp_Pnt              aP;
  gp_Vec              aDU, aDV;
  aBAS.D1(0.5 * aU1Par, 0.5 * aVn, aP, aDU, aDV);
  gp_Vec aNW = aDU ^ aDV;
  if (aWall.Orientation() == TopAbs_REVERSED)
  {
    aNW.Reverse();
  }
  const gp_Vec aUmid = aRadial(aP);
  const gp_Vec anInto =
    (aUmid * std::cos(aTh1) - gp_Vec(aA) * std::sin(aTh1)) * aS;
  if (aNW.Dot(anInto) < 0.)
  {
    aWall.Reverse();
  }
  return aWall;
}

// The type of the edge where the kept face's offset meets the wall of
// SphereCapWall, as ChFi3d::DefineConnectType gives it for a plane's wall:
// the edge lies off the kept face's sphere and has no pcurve on it, and the
// face's normal there is the sphere's, where the edge points.
static ChFiDS_TypeOfConcavity SphereCapWallEdgeType(const TopoDS_Edge& theE,
                                                    const TopoDS_Face& theWall,
                                                    const gp_Sphere&   theSph,
                                                    const gp_Dir&      theNrm,
                                                    const gp_Pnt&      theOnKept)
{
  double                          f, l;
  const occ::handle<Geom2d_Curve> aC2 = BRep_Tool::CurveOnSurface(theE, theWall, f, l);
  if (aC2.IsNull())
  {
    return ChFiDS_Other;
  }
  BRepAdaptor_Curve aBAC(theE);
  const double      aMid = 0.5 * (f + l);
  gp_Pnt            aP;
  gp_Vec            aT;
  aBAC.D1(aMid, aP, aT);
  if (aT.Magnitude() < gp::Resolution())
  {
    return ChFiDS_Other;
  }
  aT.Normalize();
  if (BRepTools::OriEdgeInFace(theE, theWall) == TopAbs_REVERSED)
  {
    aT.Reverse();
  }
  if (theWall.Orientation() == TopAbs_REVERSED)
  {
    aT.Reverse();
  }
  BRepAdaptor_Surface aBAS(theWall, false);
  const gp_Pnt2d      aUV = aC2->Value(aMid);
  gp_Pnt              aPS;
  gp_Vec              aDU, aDV;
  aBAS.D1(aUV.X(), aUV.Y(), aPS, aDU, aDV);
  gp_Vec aN1 = aDU ^ aDV;
  if (theWall.Orientation() == TopAbs_REVERSED)
  {
    aN1.Reverse();
  }
  gp_Vec aN2(theSph.Location(), aP);
  if (gp_Vec(theNrm).Dot(gp_Vec(theSph.Location(), theOnKept)) < 0.)
  {
    aN2.Reverse();
  }
  if (aN1.Magnitude() < gp::Resolution() || aN2.Magnitude() < gp::Resolution())
  {
    return ChFiDS_Other;
  }
  return aT.Dot(aN1.Normalized() ^ aN2.Normalized()) > 0. ? ChFiDS_Convex : ChFiDS_Concave;
}

void BRepOffset_Analyse::TreatTangentCaps(
  const NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>& theCaps,
  const double                                                         theOffset)
{
  const double aT = std::abs(theOffset);
  if (aT < Precision::Confusion())
  {
    return;
  }
  const double aSign   = theOffset > 0. ? 1. : -1.;
  const double aSinTol = std::abs(std::sin(myAngle));
  // Tangency as the Arc join's tube judges it (BRepOffset_MakeOffset::ToContext).
  const double aSinTang = std::abs(std::sin(Precision::Angular()));
  BRep_Builder aBB;
  for (int iC = 1; iC <= theCaps.Extent(); ++iC)
  {
    const TopoDS_Face& aCF = TopoDS::Face(theCaps(iC));
    for (TopExp_Explorer anExp(aCF, TopAbs_EDGE); anExp.More(); anExp.Next())
    {
      const TopoDS_Edge& aEC = TopoDS::Edge(anExp.Current());
      if (!myAncestors.Contains(aEC))
      {
        continue;
      }
      const NCollection_List<TopoDS_Shape>& aLA = myAncestors.FindFromKey(aEC);
      if (aLA.Extent() != 1 || theCaps.Contains(aLA.First()))
      {
        continue;
      }
      const TopoDS_Face aN = TopoDS::Face(aLA.First());
      if (ChFi3d::DefineConnectType(aEC, aCF, aN, aSinTang, false) != ChFiDS_Tangential)
      {
        continue;
      }
      BRepAdaptor_Curve aBAC(aEC);
      // Or a circle between two faces of one sphere -- the half ball with its
      // sphere in two domes or two lunes, one removed. The kept face is its
      // own strip, as a plane is, and its offset runs on round the sphere a
      // thickness into the cap, to a wall square to the sphere: a cone from
      // its centre.
      gp_Sphere  aSph;
      const bool isOnSphere = aBAC.GetType() == GeomAbs_Circle && AreOnOneSphere(aCF, aN, aSph);
      if (aBAC.GetType() != GeomAbs_Line && !isOnSphere)
      {
        continue;
      }
      TopoDS_Vertex aV[2];
      TopExp::Vertices(TopoDS::Edge(aEC.Oriented(TopAbs_FORWARD)), aV[0], aV[1]);
      if (aV[0].IsNull() || aV[1].IsNull() || aV[0].IsSame(aV[1]))
      {
        continue;
      }
      // The kept face's normal, the same all along the edge: its tangent
      // plane is the strip's.
      const double f = aBAC.FirstParameter(), l = aBAC.LastParameter();
      const gp_Dir aNrm = FaceNormalOnEdge(aEC, aN, (f + l) / 2.);
      if (!isOnSphere
          && (!aNrm.IsEqual(FaceNormalOnEdge(aEC, aN, f), 10. * Precision::Angular())
              || !aNrm.IsEqual(FaceNormalOnEdge(aEC, aN, l), 10. * Precision::Angular())))
      {
        continue;
      }
      // Into the cap: its interior lies to the left of its edge, seen along
      // its normal.
      gp_Vec aTan = aBAC.DN((f + l) / 2., 1);
      if (aEC.Orientation() == TopAbs_REVERSED)
      {
        aTan.Reverse();
      }
      gp_Vec anIn = gp_Vec(FaceNormalOnEdge(aEC, aCF, (f + l) / 2.)) ^ aTan;
      if (anIn.Magnitude() < gp::Resolution())
      {
        continue;
      }
      const gp_Dir aD(anIn);
      // The faces at each end of the edge besides the kept one: the strip
      // and the wall meet the same face there.
      TopoDS_Face aFV[2];
      bool        isOk = true;
      for (int i = 0; i < 2 && isOk; ++i)
      {
        NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aMF;
        for (NCollection_List<TopoDS_Shape>::Iterator itE(Ancestors(aV[i])); itE.More(); itE.Next())
        {
          if (itE.Value().IsSame(aEC) || !myAncestors.Contains(itE.Value()))
          {
            continue;
          }
          for (NCollection_List<TopoDS_Shape>::Iterator itF(Ancestors(itE.Value())); itF.More();
               itF.Next())
          {
            if (!itF.Value().IsSame(aN) && !theCaps.Contains(itF.Value()))
            {
              aMF.Add(itF.Value());
            }
          }
        }
        isOk = aMF.Extent() == 1;
        if (isOk)
        {
          aFV[i] = TopoDS::Face(NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>::Iterator(aMF).Value());
        }
        else if (aMF.Extent() > 1)
        {
          // More than one: the face at the end is split too -- a box fused of
          // two, not refined, every face across the joint in two pieces. The
          // closure runs into the cap, so it meets the piece beside the cap:
          // the one sharing the cap's other edge at this vertex.
          for (TopExp_Explorer anExpC(aCF, TopAbs_EDGE); anExpC.More() && !isOk; anExpC.Next())
          {
            const TopoDS_Edge& aEO = TopoDS::Edge(anExpC.Current());
            if (aEO.IsSame(aEC) || !myAncestors.Contains(aEO))
            {
              continue;
            }
            TopoDS_Vertex aO1, aO2;
            TopExp::Vertices(aEO, aO1, aO2);
            if (!aV[i].IsSame(aO1) && !aV[i].IsSame(aO2))
            {
              continue;
            }
            for (NCollection_List<TopoDS_Shape>::Iterator itF(Ancestors(aEO)); itF.More();
                 itF.Next())
            {
              if (aMF.Contains(itF.Value()))
              {
                aFV[i] = TopoDS::Face(itF.Value());
                isOk   = true;
                break;
              }
            }
          }
        }
      }
      if (!isOk)
      {
        continue;
      }
      // V -> V' (a thickness into the cap) -> V'n (across the thickness).
      TopoDS_Vertex aVp[2], aVpn[2];
      TopoDS_Edge   aS[2], aTe[2], aEp, aEpn;
      // A planar kept face is its own tangent plane: its offset runs on to
      // the wall, and a strip would only lie on it (its sections with the
      // faces at the ends doubling the kept face's). The wall then meets the
      // kept face itself. So does a face of the cap's own sphere.
      const bool isPlanarN =
        isOnSphere || BRepAdaptor_Surface(aN, false).GetType() == GeomAbs_Plane;
      TopoDS_Face aStrip, aWall;
      if (isOnSphere)
      {
        aWall = SphereCapWall(aSph, aBAC, aV, aD, aNrm, aSign * aT, aVp, aVpn, aTe, aEp, aEpn);
        if (aWall.IsNull())
        {
          continue;
        }
      }
      else
      {
        for (int i = 0; i < 2; ++i)
        {
          const gp_Pnt aP  = BRep_Tool::Pnt(aV[i]);
          const gp_Pnt aPp = aP.Translated(gp_Vec(aD) * aT);
          aBB.MakeVertex(aVp[i], aPp, Precision::Confusion());
          aBB.MakeVertex(aVpn[i],
                         aPp.Translated(gp_Vec(aNrm) * (aSign * aT)),
                         Precision::Confusion());
        }
        BRepLib_MakeEdge aMS1(aV[0], aVp[0]), aMS2(aV[1], aVp[1]), aMEp(aVp[0], aVp[1]),
          aMT1(aVp[0], aVpn[0]), aMT2(aVp[1], aVpn[1]), aMEpn(aVpn[0], aVpn[1]);
        if (!aMS1.IsDone() || !aMS2.IsDone() || !aMEp.IsDone() || !aMT1.IsDone()
            || !aMT2.IsDone() || !aMEpn.IsDone())
        {
          continue;
        }
        aS[0]  = aMS1.Edge();
        aS[1]  = aMS2.Edge();
        aTe[0] = aMT1.Edge();
        aTe[1] = aMT2.Edge();
        aEp    = aMEp.Edge();
        aEpn   = aMEpn.Edge();
        NCollection_List<TopoDS_Shape> aLStrip, aLWall;
        aLStrip.Append(aEC.Oriented(TopAbs_FORWARD));
        aLStrip.Append(aS[1]);
        aLStrip.Append(aEp);
        aLStrip.Append(aS[0]);
        aLWall.Append(aEp);
        aLWall.Append(aTe[1]);
        aLWall.Append(aEpn);
        aLWall.Append(aTe[0]);
        aStrip = isPlanarN ? TopoDS_Face() : PlanarFaceOn(aLStrip, aNrm);
        aWall  = PlanarFaceOn(aLWall, aD);
        if ((!isPlanarN && aStrip.IsNull()) || aWall.IsNull())
        {
          continue;
        }
      }
      // The face the wall's far edge meets: the strip, or the kept face.
      const TopoDS_Face& aNear = isPlanarN ? aN : aStrip;
      // The edges as the new faces hold them.
      auto anInFace = [](const TopoDS_Shape& theE, const TopoDS_Face& theF) {
        for (TopExp_Explorer anExpE(theF, TopAbs_EDGE); anExpE.More(); anExpE.Next())
        {
          if (anExpE.Current().IsSame(theE))
          {
            return anExpE.Current();
          }
        }
        return theE;
      };
      auto aReplace = [&](const TopoDS_Shape& theF, const TopoDS_Shape& theE,
                          const TopoDS_Shape& theBy) {
        NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>* pEEMap =
          myReplacement.ChangeSeek(theF);
        if (!pEEMap)
        {
          pEEMap = myReplacement.Bound(
            theF,
            NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>());
        }
        pEEMap->Bind(theE, theBy);
      };

      if (!isPlanarN)
      {
        // The tangent edge: the kept face and the strip.
        myAncestors.ChangeFromKey(aEC).Append(aStrip);
        myMapEdgeType(aEC).Clear();
        EdgeAnalyse(TopoDS::Edge(anInFace(aEC, aStrip)), aN, aStrip, aSinTol, myMapEdgeType(aEC));
      }
      // The strip's ends and the wall's: each with the face at its end.
      for (int i = 0; i < 2; ++i)
      {
        NCollection_List<TopoDS_Shape> aLVp, aLVpn;
        if (!isPlanarN)
        {
          const TopoDS_Edge aSi = TopoDS::Edge(anInFace(aS[i], aStrip));
          NCollection_List<TopoDS_Shape> aLS;
          aLS.Append(aStrip);
          aLS.Append(aFV[i]);
          myAncestors.Add(aSi, aLS);
          myMapEdgeType.Bind(aSi, NCollection_List<BRepOffset_Interval>());
          EdgeAnalyse(aSi, aStrip, aFV[i], aSinTol, myMapEdgeType(aSi));
          myAncestors.ChangeFromKey(aV[i]).Append(aSi);
          aLVp.Append(aSi);
        }

        const TopoDS_Edge aTi = TopoDS::Edge(anInFace(aTe[i], aWall));
        NCollection_List<TopoDS_Shape> aLT;
        aLT.Append(aWall);
        aLT.Append(aFV[i]);
        myAncestors.Add(aTi, aLT);
        myMapEdgeType.Bind(aTi, NCollection_List<BRepOffset_Interval>());
        EdgeAnalyse(aTi, aWall, aFV[i], aSinTol, myMapEdgeType(aTi));
        // The wall's end edge is a line square to the kept face, and lies on
        // the face at the end only where that face is a plane. On a curved
        // one -- the sphere of half a dome, the tangent edge ending at its
        // pole -- it has no pcurve to be analysed by: the edge is convex or
        // concave as the wall runs behind the face or in front of it, seen
        // from where the edge's middle falls on the face.
        if (BRepAdaptor_Surface(aFV[i], false).GetType() != GeomAbs_Plane)
        {
          const gp_Pnt aPi  = BRep_Tool::Pnt(aVp[i]), aPn = BRep_Tool::Pnt(aVpn[i]);
          const gp_Pnt aMid = aPi.XYZ() * 0.5 + aPn.XYZ() * 0.5;
          const occ::handle<Geom_Surface> aSurf = BRep_Tool::Surface(aFV[i]);
          GeomAPI_ProjectPointOnSurf      aProj(aMid, aSurf);
          if (aProj.NbPoints() > 0)
          {
            double aU, aVv;
            aProj.LowerDistanceParameters(aU, aVv);
            GeomLProp_SLProps aProps(aSurf, aU, aVv, 1, Precision::Confusion());
            if (aProps.IsNormalDefined())
            {
              gp_Dir aNF = aProps.Normal();
              if (aFV[i].Orientation() == TopAbs_REVERSED)
              {
                aNF.Reverse();
              }
              const gp_Vec aAlong(BRep_Tool::Pnt(aV[i]), BRep_Tool::Pnt(aV[1 - i]));
              double aF1, aL1;
              BRep_Tool::Range(aTi, aF1, aL1);
              NCollection_List<BRepOffset_Interval>& aLI = myMapEdgeType(aTi);
              aLI.Clear();
              aLI.Append(BRepOffset_Interval(aF1,
                                             aL1,
                                             aAlong.Dot(gp_Vec(aNF)) < 0. ? ChFiDS_Convex
                                                                          : ChFiDS_Concave));
            }
          }
        }

        aLVp.Append(anInFace(aEp, aWall));
        aLVp.Append(aTi);
        myAncestors.Add(aVp[i], aLVp);
        // The face at the edge's end meets the closure there: its edges on
        // that face are the vertex's too (BRepOffset_Inter2d::ConnexIntByInt).
        myAncestors.ChangeFromKey(aV[i]).Append(aTi);
        aLVpn.Append(aTi);
        aLVpn.Append(anInFace(aEpn, aWall));
        myAncestors.Add(aVpn[i], aLVpn);
      }
      // The wall's edge on the cap: a free border, as the cap's own edges.
      {
        NCollection_List<TopoDS_Shape> aLE;
        aLE.Append(aWall);
        const TopoDS_Edge aEpW = TopoDS::Edge(anInFace(aEp, aWall));
        myAncestors.Add(aEpW, aLE);
        double aF1, aL1;
        BRep_Tool::Range(aEpW, aF1, aL1);
        NCollection_List<BRepOffset_Interval> aLI;
        aLI.Append(BRepOffset_Interval(aF1, aL1, ChFiDS_FreeBound));
        myMapEdgeType.Bind(aEpW, aLI);
      }
      // The wall's far edge, where the offset of the strip (or of the kept
      // face) meets it: it replaces the strip's far edge (the kept face's
      // tangent edge).
      {
        const TopoDS_Edge aEpnW = TopoDS::Edge(anInFace(aEpn, aWall));
        NCollection_List<TopoDS_Shape> aLE;
        aLE.Append(aWall);
        aLE.Append(aNear);
        myAncestors.Add(aEpnW, aLE);
        myMapEdgeType.Bind(aEpnW, NCollection_List<BRepOffset_Interval>());
        if (isOnSphere)
        {
          // Off the kept face's sphere, the edge has no pcurve on it: the
          // kept face's normal is the sphere's, where the edge points.
          double aF1, aL1;
          BRep_Tool::Range(aEpnW, aF1, aL1);
          myMapEdgeType(aEpnW).Append(
            BRepOffset_Interval(aF1,
                                aL1,
                                SphereCapWallEdgeType(aEpnW, aWall, aSph, aNrm, aBAC.Value((f + l) / 2.))));
        }
        else
        {
          EdgeAnalyse(aEpnW, aWall, aNear, aSinTol, myMapEdgeType(aEpnW));
        }
        aReplace(aNear, isPlanarN ? aEC : aEp, aEpnW);
      }
      // In the cap, the tangent edge gives way to the wall's edge.
      aReplace(aCF, aEC, anInFace(aEp, aWall));
      if (!isPlanarN)
      {
        myNewFaces.Append(aStrip);
        myFaceOffsetMap.Bind(aStrip, theOffset);
      }
      myNewFaces.Append(aWall);
    }
  }
  myDescendants.Clear();
}

//=================================================================================================

const TopoDS_Edge& BRepOffset_Analyse::EdgeReplacement(const TopoDS_Face& theF,
                                                       const TopoDS_Edge& theE) const
{
  const NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>* pEE =
    myReplacement.Seek(theF);
  if (!pEE)
  {
    return theE;
  }

  const TopoDS_Shape* pE = pEE->Seek(theE);
  if (!pE)
  {
    return theE;
  }

  return TopoDS::Edge(*pE);
}

//=================================================================================================

TopoDS_Shape BRepOffset_Analyse::Generated(const TopoDS_Shape& theS) const
{
  static TopoDS_Shape aNullShape;
  const TopoDS_Shape* pGenS = myGenerated.Seek(theS);
  return pGenS ? *pGenS : aNullShape;
}

//=================================================================================================

const NCollection_List<TopoDS_Shape>* BRepOffset_Analyse::Descendants(const TopoDS_Shape& theS,
                                                                      const bool theUpdate) const
{
  if (myDescendants.IsEmpty() || theUpdate)
  {
    myDescendants.Clear();
    const int aNbA = myAncestors.Extent();
    for (int i = 1; i <= aNbA; ++i)
    {
      const TopoDS_Shape&                   aSS = myAncestors.FindKey(i);
      const NCollection_List<TopoDS_Shape>& aLA = myAncestors(i);

      for (NCollection_List<TopoDS_Shape>::Iterator it(aLA); it.More(); it.Next())
      {
        const TopoDS_Shape& aSA = it.Value();

        NCollection_List<TopoDS_Shape>* pLD = myDescendants.ChangeSeek(aSA);
        if (!pLD)
        {
          pLD = myDescendants.Bound(aSA, NCollection_List<TopoDS_Shape>());
        }
        if (!pLD->Contains(aSS))
        {
          pLD->Append(aSS);
        }
      }
    }
  }

  return myDescendants.Seek(theS);
}

//=================================================================================================

void BRepOffset_Analyse::Clear()
{
  myDone = false;
  myShape.Nullify();
  myMapEdgeType.Clear();
  myAncestors.Clear();
  myFaceOffsetMap.Clear();
  myReplacement.Clear();
  myDescendants.Clear();
  myNewFaces.Clear();
  myGenerated.Clear();
}

//=======================================================================
// function : NCollection_List<BRepOffset_Interval>&
// purpose  :
//=======================================================================
const NCollection_List<BRepOffset_Interval>& BRepOffset_Analyse::Type(const TopoDS_Edge& E) const
{
  return myMapEdgeType(E);
}

//=================================================================================================

void BRepOffset_Analyse::Edges(const TopoDS_Vertex&            V,
                               const ChFiDS_TypeOfConcavity    T,
                               NCollection_List<TopoDS_Shape>& LE) const
{
  LE.Clear();
  const NCollection_List<TopoDS_Shape>&    L = Ancestors(V);
  NCollection_List<TopoDS_Shape>::Iterator it(L);

  for (; it.More(); it.Next())
  {
    const TopoDS_Edge&                           E          = TopoDS::Edge(it.Value());
    const NCollection_List<BRepOffset_Interval>* pIntervals = myMapEdgeType.Seek(E);
    if (pIntervals && pIntervals->Extent() > 0)
    {
      TopoDS_Vertex V1, V2;
      BRepOffset_Tool::EdgeVertices(E, V1, V2);
      if (V1.IsSame(V))
      {
        if (pIntervals->Last().Type() == T)
        {
          LE.Append(E);
        }
      }
      if (V2.IsSame(V))
      {
        if (pIntervals->First().Type() == T)
        {
          LE.Append(E);
        }
      }
    }
  }
}

//=================================================================================================

void BRepOffset_Analyse::Edges(const TopoDS_Face&              F,
                               const ChFiDS_TypeOfConcavity    T,
                               NCollection_List<TopoDS_Shape>& LE) const
{
  LE.Clear();
  TopExp_Explorer exp(F, TopAbs_EDGE);

  for (; exp.More(); exp.Next())
  {
    const TopoDS_Edge& E = TopoDS::Edge(exp.Current());

    const NCollection_List<BRepOffset_Interval>&    Lint = Type(E);
    NCollection_List<BRepOffset_Interval>::Iterator it(Lint);
    for (; it.More(); it.Next())
    {
      if (it.Value().Type() == T)
      {
        LE.Append(E);
      }
    }
  }
}

//=================================================================================================

void BRepOffset_Analyse::TangentEdges(const TopoDS_Edge&              Edge,
                                      const TopoDS_Vertex&            Vertex,
                                      NCollection_List<TopoDS_Shape>& Edges) const
{
  gp_Vec V, VRef;

  double            U, URef;
  BRepAdaptor_Curve C3d, C3dRef;

  URef   = BRep_Tool::Parameter(Vertex, Edge);
  C3dRef = BRepAdaptor_Curve(Edge);
  VRef   = C3dRef.DN(URef, 1);
  CorrectOrientationOfTangent(VRef, Vertex, Edge);
  if (VRef.SquareMagnitude() < gp::Resolution())
  {
    return;
  }

  Edges.Clear();

  const NCollection_List<TopoDS_Shape>&    Anc = Ancestors(Vertex);
  NCollection_List<TopoDS_Shape>::Iterator it(Anc);
  for (; it.More(); it.Next())
  {
    const TopoDS_Edge& CurE = TopoDS::Edge(it.Value());
    if (CurE.IsSame(Edge))
    {
      continue;
    }
    U   = BRep_Tool::Parameter(Vertex, CurE);
    C3d = BRepAdaptor_Curve(CurE);
    V   = C3d.DN(U, 1);
    CorrectOrientationOfTangent(V, Vertex, CurE);
    if (V.SquareMagnitude() < gp::Resolution())
    {
      continue;
    }
    if (V.IsOpposite(VRef, myAngle))
    {
      Edges.Append(CurE);
    }
  }
}

//=================================================================================================

void BRepOffset_Analyse::Explode(NCollection_List<TopoDS_Shape>& List,
                                 const ChFiDS_TypeOfConcavity    T) const
{
  List.Clear();
  BRep_Builder                                           B;
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> Map;

  TopExp_Explorer Fexp;
  for (Fexp.Init(myShape, TopAbs_FACE); Fexp.More(); Fexp.Next())
  {
    if (Map.Add(Fexp.Current()))
    {
      TopoDS_Face     Face = TopoDS::Face(Fexp.Current());
      TopoDS_Compound Co;
      B.MakeCompound(Co);
      B.Add(Co, Face);
      // add to Co all faces from the cloud of faces
      // G1 created from <Face>
      AddFaces(Face, Co, Map, T);
      List.Append(Co);
    }
  }
}

//=================================================================================================

void BRepOffset_Analyse::Explode(NCollection_List<TopoDS_Shape>& List,
                                 const ChFiDS_TypeOfConcavity    T1,
                                 const ChFiDS_TypeOfConcavity    T2) const
{
  List.Clear();
  BRep_Builder                                           B;
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> Map;

  TopExp_Explorer Fexp;
  for (Fexp.Init(myShape, TopAbs_FACE); Fexp.More(); Fexp.Next())
  {
    if (Map.Add(Fexp.Current()))
    {
      TopoDS_Face     Face = TopoDS::Face(Fexp.Current());
      TopoDS_Compound Co;
      B.MakeCompound(Co);
      B.Add(Co, Face);
      // add to Co all faces from the cloud of faces
      // G1 created from  <Face>
      AddFaces(Face, Co, Map, T1, T2);
      List.Append(Co);
    }
  }
}

//=================================================================================================

void BRepOffset_Analyse::AddFaces(const TopoDS_Face&                                      Face,
                                  TopoDS_Compound&                                        Co,
                                  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& Map,
                                  const ChFiDS_TypeOfConcavity                            T) const
{
  BRep_Builder                          B;
  const NCollection_List<TopoDS_Shape>* pLE = Descendants(Face);
  if (!pLE)
  {
    return;
  }
  for (NCollection_List<TopoDS_Shape>::Iterator it(*pLE); it.More(); it.Next())
  {
    const TopoDS_Edge&                           E  = TopoDS::Edge(it.Value());
    const NCollection_List<BRepOffset_Interval>& LI = Type(E);
    if (!LI.IsEmpty() && LI.First().Type() == T)
    {
      // so <NewFace> is attached to G1 by <Face>
      const NCollection_List<TopoDS_Shape>& L = Ancestors(E);
      if (L.Extent() == 2)
      {
        TopoDS_Face F1 = TopoDS::Face(L.First());
        if (F1.IsSame(Face))
        {
          F1 = TopoDS::Face(L.Last());
        }
        if (Map.Add(F1))
        {
          B.Add(Co, F1);
          AddFaces(F1, Co, Map, T);
        }
      }
    }
  }
}

//=================================================================================================

void BRepOffset_Analyse::AddFaces(const TopoDS_Face&                                      Face,
                                  TopoDS_Compound&                                        Co,
                                  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& Map,
                                  const ChFiDS_TypeOfConcavity                            T1,
                                  const ChFiDS_TypeOfConcavity                            T2) const
{
  BRep_Builder                          B;
  const NCollection_List<TopoDS_Shape>* pLE = Descendants(Face);
  if (!pLE)
  {
    return;
  }
  for (NCollection_List<TopoDS_Shape>::Iterator it(*pLE); it.More(); it.Next())
  {
    const TopoDS_Edge&                           E  = TopoDS::Edge(it.Value());
    const NCollection_List<BRepOffset_Interval>& LI = Type(E);
    if (!LI.IsEmpty() && (LI.First().Type() == T1 || LI.First().Type() == T2))
    {
      // so <NewFace> is attached to G1 by <Face>
      const NCollection_List<TopoDS_Shape>& L = Ancestors(E);
      if (L.Extent() == 2)
      {
        TopoDS_Face F1 = TopoDS::Face(L.First());
        if (F1.IsSame(Face))
        {
          F1 = TopoDS::Face(L.Last());
        }
        if (Map.Add(F1))
        {
          B.Add(Co, F1);
          AddFaces(F1, Co, Map, T1, T2);
        }
      }
    }
  }
}
