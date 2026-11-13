// Created by: Peter KURNEV
// Copyright (c) 2010-2014 OPEN CASCADE SAS
// Copyright (c) 2007-2010 CEA/DEN, EDF R&D, OPEN CASCADE
// Copyright (c) 2003-2007 OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN, CEDRAT,
//                         EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
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

#include <BOPAlgo_PaveFiller.hxx>
#include <Bnd_Box.hxx>
#include <BOPAlgo_Alerts.hxx>
#include <BOPAlgo_SectionAttribute.hxx>
#include <BOPAlgo_Tools.hxx>
#include <BOPDS_CoupleOfPaveBlocks.hxx>
#include <BOPDS_Curve.hxx>
#include <NCollection_DataMap.hxx>
#include <BOPDS_Pair.hxx>
#include <BOPDS_PaveBlock.hxx>
#include <Standard_Handle.hxx>
#include <BOPDS_DS.hxx>
#include <BOPDS_FaceInfo.hxx>
#include <BOPDS_Interf.hxx>
#include <BOPDS_Iterator.hxx>
#include <NCollection_List.hxx>
#include <BOPDS_Pave.hxx>
#include <NCollection_Map.hxx>
#include <BOPDS_Point.hxx>
#include <BOPDS_ShapeInfo.hxx>
#include <NCollection_DynamicArray.hxx>
#include <BOPTools_AlgoTools.hxx>
#include <BOPTools_AlgoTools3D.hxx>
#include <BOPTools_BoxSelector.hxx>
#include <BOPTools_Parallel.hxx>
#include <Bnd_Tools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_TEdge.hxx>
#include <BRep_Tool.hxx>
#include <BRepBndLib.hxx>
#include <BRepLib.hxx>
#include <BRepTools.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Plane.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_OffsetSurface.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <gp_Pnt.hxx>
#include <IntSurf_PntOn2S.hxx>
#include <Standard_ErrorHandler.hxx>
#include <IntTools.hxx>
#include <IntTools_Context.hxx>
#include <IntTools_Curve.hxx>
#include <IntTools_EdgeFace.hxx>
#include <IntTools_FaceFace.hxx>
#include <IntTools_PntOn2Faces.hxx>
#include <NCollection_Sequence.hxx>
#include <IntTools_Tools.hxx>
#include <NCollection_IncAllocator.hxx>
#include <Precision.hxx>
#include <Standard_Integer.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Shape.hxx>
#include <TopTools_ShapeMapHasher.hxx>
class BOPDS_PaveBlock;

static double ToleranceFF(const BRepAdaptor_Surface& aBAS1, const BRepAdaptor_Surface& aBAS2);
static void   getBoundPaves(const BOPDS_DS* theDS, const BOPDS_Curve& theNC, int theNV[2]);
static int    tangentEndVertex(const BOPDS_PDS&                     theDS,
                               const occ::handle<IntTools_Context>& theCtx,
                               const TopoDS_Face&                   theF1,
                               const TopoDS_Face&                   theF2,
                               const BOPDS_Curve&                   theNC,
                               const int                            theEnd,
                               const NCollection_Map<int>&          theMV,
                               double&                              theDist);
static int    continuedEndVertex(const BOPDS_PDS&                     theDS,
                                 const occ::handle<IntTools_Context>& theCtx,
                                 const TopoDS_Face&                   theF1,
                                 const TopoDS_Face&                   theF2,
                                 const BOPDS_Curve&                   theNC,
                                 const int                            theEnd,
                                 const NCollection_Map<int>&          theMV,
                                 double&                              theDist);
static bool   evenParameters(const occ::handle<Geom_Curve>&    theC,
                             const double                      theT1,
                             const double                      theT2,
                             const int                         theNb,
                             NCollection_DynamicArray<double>& thePars);
static void   dropStretchedPaves(const BOPDS_PDS&                             theDS,
                                 const occ::handle<IntTools_Context>&         theCtx,
                                 const NCollection_DynamicArray<BOPDS_Curve>& theVC,
                                 const NCollection_Map<int>&                  theMI,
                                 const double                                 theFuzzy,
                                 NCollection_DataMap<int, double>&            theMVTol);
static bool   curveAlongEdges(
  const BOPDS_PDS&                                            theDS,
  const occ::handle<IntTools_Context>&                        theCtx,
  const BOPDS_Curve&                                          theNC,
  const NCollection_IndexedMap<occ::handle<BOPDS_PaveBlock>>& theMPBOnIn,
  const NCollection_Map<occ::handle<BOPDS_PaveBlock>>&        theMPBCommon);
static bool   facesTogether(const occ::handle<IntTools_Context>& theCtx,
                            const TopoDS_Face&                   theF1,
                            const TopoDS_Face&                   theF2,
                            const gp_Pnt&                        thePA,
                            const gp_Pnt&                        thePB,
                            const double                         theTol,
                            NCollection_Array1<gp_Pnt>*          thePnts = nullptr,
                            double*                              theDist = nullptr,
                            NCollection_Array1<gp_Pnt2d>*        theUV1  = nullptr,
                            NCollection_Array1<gp_Pnt2d>*        theUV2  = nullptr);
static void   sharedEdges(const BOPDS_PDS&                                            theDS,
                          const NCollection_IndexedMap<occ::handle<BOPDS_PaveBlock>>& theMPBOnIn,
                          const NCollection_Map<occ::handle<BOPDS_PaveBlock>>&        theMPBCommon,
                          NCollection_List<int>&                                      theLE);
static bool   onSharedEdge(const BOPDS_PDS&                     theDS,
                           const occ::handle<IntTools_Context>& theCtx,
                           const TopoDS_Face&                   theF1,
                           const TopoDS_Face&                   theF2,
                           const gp_Pnt&                        theP,
                           const NCollection_List<int>&         theLE,
                           const double                         theTol,
                           const double                         theDMax);
static bool   lostAlongEdges(const BOPDS_PDS&                     theDS,
                             const occ::handle<IntTools_Context>& theCtx,
                             const TopoDS_Face&                   theF1,
                             const TopoDS_Face&                   theF2,
                             const BOPDS_Curve&                   theNC,
                             const NCollection_List<int>&         theLE);
static bool   trimEndsOnEdges(const BOPDS_PDS&                     theDS,
                              const occ::handle<IntTools_Context>& theCtx,
                              const TopoDS_Face&                   theF1,
                              const TopoDS_Face&                   theF2,
                              BOPDS_Curve&                         theNC,
                              const NCollection_List<int>&         theLE);
static bool   bridgeToVertex(const BOPDS_PDS&                       theDS,
                             const occ::handle<IntTools_Context>&   theCtx,
                             const TopoDS_Face&                     theF1,
                             const TopoDS_Face&                     theF2,
                             NCollection_DynamicArray<BOPDS_Curve>& theVC,
                             const int                              theJ,
                             const int                              theNV,
                             const double                           theTolV,
                             const double                           theTolFF,
                             const NCollection_List<int>&           theLE,
                             occ::handle<Geom_Curve>&               theBridge,
                             occ::handle<Geom2d_Curve>              theBridge2d[2],
                             double&                                theTolB,
                             NCollection_DataMap<int, int>&         theOwners);
static void   tangentVerticesNearCurve(const BOPDS_PDS&                        theDS,
                                       const occ::handle<IntTools_Context>&    theCtx,
                                       const TopoDS_Face&                      theF1,
                                       const TopoDS_Face&                      theF2,
                                       const BOPDS_Curve&                      theNC,
                                       const NCollection_Map<int>&             theMV,
                                       const NCollection_DataMap<int, double>& theMVTol,
                                       NCollection_List<int>&                  theLV,
                                       NCollection_List<double>&               theLD);

//=================================================================================================

static bool IsPlaneFF(const occ::handle<Geom_Surface>& theSurface)
{
  if (theSurface->IsKind(STANDARD_TYPE(Geom_Plane)))
  {
    return true;
  }
  else if (const occ::handle<Geom_OffsetSurface> anOffsetSurface =
             occ::down_cast<Geom_OffsetSurface>(theSurface))
  {
    return anOffsetSurface->BasisSurface()->IsKind(STANDARD_TYPE(Geom_Plane));
  }
  else if (const occ::handle<Geom_RectangularTrimmedSurface> aTrimmedSurface =
             occ::down_cast<Geom_RectangularTrimmedSurface>(theSurface))
  {
    return aTrimmedSurface->BasisSurface()->IsKind(STANDARD_TYPE(Geom_Plane));
  }

  return false;
}

//=================================================================================================

static bool IsClosedFF(const TopoDS_Edge&                     theEdge,
                       const occ::handle<Geom_Surface>&       theSurface,
                       const occ::handle<Poly_Triangulation>& theTriangulation,
                       const TopLoc_Location&                 theLocation,
                       bool                                   theIsPlane)
{
  if (!theIsPlane)
  {
    // Check surface
    const TopLoc_Location aLocation = theLocation.Predivided(theEdge.Location());

    // find the representation
    const BRep_TEdge* aTEdge = static_cast<const BRep_TEdge*>(theEdge.TShape().get());
    for (NCollection_List<occ::handle<BRep_CurveRepresentation>>::Iterator anIter(aTEdge->Curves());
         anIter.More();
         anIter.Next())
    {
      const occ::handle<BRep_CurveRepresentation>& aCurveRepresentation = anIter.Value();
      if (aCurveRepresentation->IsCurveOnSurface(theSurface, aLocation)
          && aCurveRepresentation->IsCurveOnClosedSurface())
      {
        return true;
      }
    }
  }

  // Check triangulation
  return BRep_Tool::IsClosed(theEdge, theTriangulation, theLocation);
}

//=================================================================================================

class BOPAlgo_FaceFace : public IntTools_FaceFace, public BOPAlgo_ParallelAlgo
{

public:
  DEFINE_STANDARD_ALLOC

  BOPAlgo_FaceFace()
      : BOPAlgo_ParallelAlgo(),
        myIF1(-1),
        myIF2(-1),
        myTolFF(1.e-7)
  {
  }

  //
  ~BOPAlgo_FaceFace() override = default;

  //
  void SetIndices(const int nF1, const int nF2)
  {
    myIF1 = nF1;
    myIF2 = nF2;
  }

  //
  void Indices(int& nF1, int& nF2) const
  {
    nF1 = myIF1;
    nF2 = myIF2;
  }

  //
  void SetFaces(const TopoDS_Face& aF1, const TopoDS_Face& aF2)
  {
    myF1 = aF1;
    myF2 = aF2;
  }

  //
  void SetBoxes(const Bnd_Box& theBox1, const Bnd_Box& theBox2)
  {
    myBox1 = theBox1;
    myBox2 = theBox2;
  }

  //
  const TopoDS_Face& Face1() const { return myF1; }

  //
  const TopoDS_Face& Face2() const { return myF2; }

  //
  void SetTolFF(const double aTolFF) { myTolFF = aTolFF; }

  //
  double TolFF() const { return myTolFF; }

  //
  void SetFuzzyValue(const double theFuzz) { IntTools_FaceFace::SetFuzzyValue(theFuzz); }

  //
  const gp_Trsf& Trsf() const { return myTrsf; }

  //
  void Perform() override
  {
    Message_ProgressScope aPS(myProgressRange, nullptr, 1);
    if (UserBreak(aPS))
    {
      return;
    }
    try
    {
      OCC_CATCH_SIGNALS

      gp_Trsf     aTrsf;
      TopoDS_Face aF1 = myF1, aF2 = myF2;
      if (BOPAlgo_Tools::TrsfToPoint(myBox1, myBox2, aTrsf))
      {
        // Shapes are located far from origin, move the shapes to the origin,
        // to increase the accuracy of intersection.
        TopLoc_Location aLoc(aTrsf);
        aF1.Move(aLoc);
        aF2.Move(aLoc);

        // The starting point is initialized only with the UV parameters
        // on the faces - 3D point is not set (see GetEFPnts method),
        // so no need to transform anything.
        // for (NCollection_List<IntSurf_PntOn2S>::Iterator it (myListOfPnts); it.More(); it.Next())
        //{
        //  IntSurf_PntOn2S& aP2S = it.ChangeValue();
        //  aP2S.SetValue (aP2S.Value().Transformed (aTrsf));
        //}

        myTrsf = aTrsf.Inverted();
      }

      IntTools_FaceFace::Perform(aF1, aF2, myRunParallel);
    }
    catch (Standard_Failure const&)
    {
      AddError(new BOPAlgo_AlertIntersectionFailed);
    }
  }

  //
  void ApplyTrsf()
  {
    if (IsDone())
    {
      // Update curves
      for (int i = 1; i <= mySeqOfCurve.Length(); ++i)
      {
        IntTools_Curve& aIC = mySeqOfCurve(i);
        aIC.Curve()->Transform(myTrsf);
      }
      // Update points
      for (int i = 1; i <= myPnts.Length(); ++i)
      {
        IntTools_PntOn2Faces& aP2F   = myPnts(i);
        IntTools_PntOnFace    aPOnF1 = aP2F.P1(), aPOnF2 = aP2F.P2();
        aPOnF1.SetPnt(aPOnF1.Pnt().Transformed(myTrsf));
        aPOnF2.SetPnt(aPOnF2.Pnt().Transformed(myTrsf));
        aP2F.SetP1(aPOnF1);
        aP2F.SetP2(aPOnF2);
      }
    }
  }

  //
protected:
  int         myIF1;
  int         myIF2;
  double      myTolFF;
  TopoDS_Face myF1;
  TopoDS_Face myF2;
  Bnd_Box     myBox1;
  Bnd_Box     myBox2;
  gp_Trsf     myTrsf;
};

//=================================================================================================

typedef NCollection_DynamicArray<BOPAlgo_FaceFace> BOPAlgo_VectorOfFaceFace;

//=================================================================================================

void BOPAlgo_PaveFiller::PerformFF(const Message_ProgressRange& theRange)
{
  // Update face info for all Face/Face intersection pairs
  // and also for the rest of the faces with FaceInfo already initialized,
  // i.e. anyhow touched faces.
  myIterator->Initialize(TopAbs_FACE, TopAbs_FACE);
  int iSize = myIterator->ExpectedLength();

  // Collect faces from intersection pairs
  NCollection_Map<int> aMIFence;
  int                  nF1, nF2;
  for (; myIterator->More(); myIterator->Next())
  {
    myIterator->Value(nF1, nF2);
    aMIFence.Add(nF1);
    aMIFence.Add(nF2);
  }
  // Collect the rest of the touched faces
  for (int i = 0; i < myDS->NbSourceShapes(); ++i)
  {
    const BOPDS_ShapeInfo& aSI = myDS->ShapeInfo(i);
    if (aSI.ShapeType() == TopAbs_FACE && aSI.HasReference())
    {
      aMIFence.Add(i);
    }
  }
  // Update face info
  myDS->UpdateFaceInfoOn(aMIFence);
  myDS->UpdateFaceInfoIn(aMIFence);

  if (!iSize)
  {
    // no intersection pairs found
    return;
  }

  Message_ProgressScope aPSOuter(theRange, nullptr, 1);

  NCollection_DynamicArray<BOPDS_InterfFF>& aFFs = myDS->InterfFF();
  aFFs.SetIncrement(iSize);
  //
  // Options for the intersection algorithm
  bool bApprox = mySectionAttribute.Approximation(), bCompC2D1 = mySectionAttribute.PCurveOnS1(),
       bCompC2D2     = mySectionAttribute.PCurveOnS2();
  double anApproxTol = 1.e-7;
  // Post-processing options
  bool bSplitCurve = false;
  //
  // Collect all pairs of Edge/Edge interferences to check if
  // some faces have to be moved to obtain more precise intersection
  NCollection_DataMap<BOPDS_Pair, NCollection_List<int>> aEEMap;
  const NCollection_DynamicArray<BOPDS_InterfEE>&        aVEEs = myDS->InterfEE();
  for (int iEE = 0; iEE < aVEEs.Length(); ++iEE)
  {
    const BOPDS_Interf& aEE = aVEEs(iEE);
    if (!aEE.HasIndexNew())
    {
      continue;
    }
    int nE1, nE2;
    aEE.Indices(nE1, nE2);

    const int nVN = aEE.IndexNew();

    BOPDS_Pair             aPair(nE1, nE2);
    NCollection_List<int>* pPoints = aEEMap.ChangeSeek(aPair);
    if (pPoints)
    {
      pPoints->Append(nVN);
    }
    else
    {
      pPoints = aEEMap.Bound(BOPDS_Pair(nE1, nE2), NCollection_List<int>());
      pPoints->Append(nVN);
    }
  }

  // Prepare the pairs of faces for intersection
  BOPAlgo_VectorOfFaceFace aVFaceFace;
  myIterator->Initialize(TopAbs_FACE, TopAbs_FACE);
  for (; myIterator->More(); myIterator->Next())
  {
    if (UserBreak(aPSOuter))
    {
      return;
    }
    myIterator->Value(nF1, nF2);

    if (myGlue == BOPAlgo_GlueOff)
    {
      const TopoDS_Face& aF1 = (*(TopoDS_Face*)(&myDS->Shape(nF1)));
      const TopoDS_Face& aF2 = (*(TopoDS_Face*)(&myDS->Shape(nF2)));
      //
      const BRepAdaptor_Surface& aBAS1 = myContext->SurfaceAdaptor(aF1);
      const BRepAdaptor_Surface& aBAS2 = myContext->SurfaceAdaptor(aF2);
      if (aBAS1.GetType() == GeomAbs_Plane && aBAS2.GetType() == GeomAbs_Plane)
      {
        // Check if the planes are really interfering
        bool bToIntersect = CheckPlanes(nF1, nF2);
        if (!bToIntersect)
        {
          BOPDS_InterfFF& aFF = aFFs.Appended();
          aFF.SetIndices(nF1, nF2);
          aFF.Init(0, 0);
          continue;
        }
      }

      // Check if there is an intersection between edges of the faces.
      // If there is an intersection, check if there is a shift between the edges
      // (intersection point is on some distance from the edges), and move one of
      // the faces to the point of exact edges intersection. This should allow
      // obtaining more precise intersection curves between the faces
      // (at least the curves should reach the boundary).
      // Note, that currently this check considers only closed edges (seam edges).
      TopoDS_Face aFShifted1 = aF1, aFShifted2 = aF2;
      // Keep shift value to use it as the tolerance for intersection curves
      double aShiftValue = 0.;

      if (aBAS1.GetType() != GeomAbs_Plane || aBAS2.GetType() != GeomAbs_Plane)
      {
        TopLoc_Location                        aLocation1;
        const occ::handle<Geom_Surface>&       aSurface1 = BRep_Tool::Surface(aF1, aLocation1);
        const occ::handle<Poly_Triangulation>& aTriangulation1 =
          BRep_Tool::Triangulation(aF1, aLocation1);
        const bool anIsPlane1 = IsPlaneFF(aSurface1);

        TopLoc_Location                        aLocation2;
        const occ::handle<Geom_Surface>&       aSurface2 = BRep_Tool::Surface(aF2, aLocation2);
        const occ::handle<Poly_Triangulation>& aTriangulation2 =
          BRep_Tool::Triangulation(aF2, aLocation2);
        const bool anIsPlane2 = IsPlaneFF(aSurface2);

        bool anIsFound = false;
        for (TopoDS_Iterator aItW1(aF1); !anIsFound && aItW1.More(); aItW1.Next())
        {
          for (TopoDS_Iterator aItE1(aItW1.Value()); !anIsFound && aItE1.More(); aItE1.Next())
          {
            const TopoDS_Edge& anEdge1      = TopoDS::Edge(aItE1.Value());
            const int          anEdgeIndex1 = myDS->Index(anEdge1);

            for (TopoDS_Iterator aItW2(aF2); !anIsFound && aItW2.More(); aItW2.Next())
            {
              for (TopoDS_Iterator aItE2(aItW2.Value()); !anIsFound && aItE2.More(); aItE2.Next())
              {
                const TopoDS_Edge& anEdge2      = TopoDS::Edge(aItE2.Value());
                const int          anEdgeIndex2 = myDS->Index(anEdge2);

                const bool anIsClosed1 =
                  IsClosedFF(anEdge1, aSurface1, aTriangulation1, aLocation1, anIsPlane1);
                const bool anIsClosed2 =
                  IsClosedFF(anEdge2, aSurface2, aTriangulation2, aLocation2, anIsPlane2);
                if (!anIsClosed1 && !anIsClosed2)
                {
                  continue;
                }

                const NCollection_List<int>* aVertexIndices =
                  aEEMap.Seek(BOPDS_Pair(anEdgeIndex1, anEdgeIndex2));
                if (!aVertexIndices)
                {
                  continue;
                }

                for (NCollection_List<int>::Iterator itEEP(*aVertexIndices); itEEP.More();
                     itEEP.Next())
                {
                  const int            aVertexIndex = itEEP.Value();
                  const TopoDS_Vertex& aVertex      = TopoDS::Vertex(myDS->Shape(aVertexIndex));
                  const gp_Pnt&        aVertexPoint = BRep_Tool::Pnt(aVertex);

                  // Compute points exactly on the edges
                  GeomAPI_ProjectPointOnCurve& aProjPC1 = myContext->ProjPC(anEdge1);
                  GeomAPI_ProjectPointOnCurve& aProjPC2 = myContext->ProjPC(anEdge2);
                  aProjPC1.Perform(aVertexPoint);
                  aProjPC2.Perform(aVertexPoint);
                  if (!aProjPC1.NbPoints() && !aProjPC2.NbPoints())
                  {
                    continue;
                  }
                  const gp_Pnt aP1 =
                    aProjPC1.NbPoints() > 0 ? aProjPC1.NearestPoint() : aVertexPoint;
                  const gp_Pnt aP2 =
                    aProjPC2.NbPoints() > 0 ? aProjPC2.NearestPoint() : aVertexPoint;

                  const double aShiftDist = aP1.Distance(aP2);
                  if (aShiftDist > BRep_Tool::Tolerance(aVertex))
                  {
                    // Move one of the faces to the point of exact intersection of edges
                    gp_Trsf aTrsf;
                    aTrsf.SetTranslation(anIsClosed1 ? gp_Vec(aP1, aP2) : gp_Vec(aP2, aP1));
                    TopLoc_Location aLoc(aTrsf);
                    (anIsClosed1 ? &aFShifted1 : &aFShifted2)->Move(aLoc);
                    aShiftValue = aShiftDist;
                    anIsFound   = true;
                  }
                }
              }
            }
          }
        }
      }
      //
      BOPAlgo_FaceFace& aFaceFace = aVFaceFace.Appended();
      //
      aFaceFace.SetRunParallel(myRunParallel);
      aFaceFace.SetIndices(nF1, nF2);
      aFaceFace.SetFaces(aFShifted1, aFShifted2);
      aFaceFace.SetBoxes(myDS->ShapeInfo(nF1).Box(), myDS->ShapeInfo(nF2).Box());
      // Note: in case of faces with closed edges it should not be less than value of the shift
      double aTolFF = std::max(aShiftValue, ToleranceFF(aBAS1, aBAS2));
      aFaceFace.SetTolFF(aTolFF);
      //
      NCollection_List<IntSurf_PntOn2S> aListOfPnts;
      GetEFPnts(nF1, nF2, aListOfPnts);
      int aNbLP = aListOfPnts.Extent();
      if (aNbLP)
      {
        aFaceFace.SetList(aListOfPnts);
      }
      //
      aFaceFace.SetParameters(bApprox, bCompC2D1, bCompC2D2, anApproxTol);
      aFaceFace.SetFuzzyValue(myFuzzyValue);
    }
    else
    {
      // for the Glue mode just add all interferences of that type
      BOPDS_InterfFF& aFF = aFFs.Appended();
      aFF.SetIndices(nF1, nF2);
      aFF.SetTangentFaces(false);
      aFF.Init(0, 0);
    }
  } // for (; myIterator->More(); myIterator->Next()) {
  //
  int k, aNbFaceFace = aVFaceFace.Length();
  ;
  Message_ProgressScope aPS(aPSOuter.Next(), "Performing Face-Face intersection", aNbFaceFace);
  for (k = 0; k < aNbFaceFace; k++)
  {
    BOPAlgo_FaceFace& aFaceFace = aVFaceFace.ChangeValue(k);
    aFaceFace.SetProgressRange(aPS.Next());
  }
  //======================================================
  // Perform intersection
  BOPTools_Parallel::Perform(myRunParallel, aVFaceFace);
  if (UserBreak(aPSOuter))
  {
    return;
  }
  //======================================================
  // Treatment of the results

  for (k = 0; k < aNbFaceFace; ++k)
  {
    if (UserBreak(aPSOuter))
    {
      return;
    }
    BOPAlgo_FaceFace& aFaceFace = aVFaceFace(k);
    aFaceFace.Indices(nF1, nF2);
    if (!aFaceFace.IsDone() || aFaceFace.HasErrors())
    {
      BOPDS_InterfFF& aFF = aFFs.Appended();
      aFF.SetIndices(nF1, nF2);
      aFF.Init(0, 0);
      // Warn about failed intersection of faces
      AddIntersectionFailedWarning(aFaceFace.Face1(), aFaceFace.Face2());
      continue;
    }
    //
    bool   bTangentFaces = aFaceFace.TangentFaces();
    double aTolFF        = aFaceFace.TolFF();
    //
    aFaceFace.PrepareLines3D(bSplitCurve);
    //
    aFaceFace.ApplyTrsf();
    //
    const NCollection_Sequence<IntTools_Curve>&       aCvsX  = aFaceFace.Lines();
    const NCollection_Sequence<IntTools_PntOn2Faces>& aPntsX = aFaceFace.Points();
    //
    int aNbCurves = aCvsX.Length();
    int aNbPoints = aPntsX.Length();
    //
    if (aNbCurves || aNbPoints)
    {
      myDS->AddInterf(nF1, nF2);
    }
    //
    BOPDS_InterfFF& aFF = aFFs.Appended();
    aFF.SetIndices(nF1, nF2);
    aFF.SetTangentFaces(bTangentFaces);
    aFF.Init(aNbCurves, aNbPoints);
    //
    // Curves
    // Fix bounding box expanding coefficient.
    double aBoxExpandValue = aTolFF;
    if (aNbCurves > 0)
    {
      // Modify geometric expanding coefficient by topology value,
      // since this bounding box used in sharing (vertex or edge).
      double aMaxVertexTol = std::max(BRep_Tool::MaxTolerance(aFaceFace.Face1(), TopAbs_VERTEX),
                                      BRep_Tool::MaxTolerance(aFaceFace.Face2(), TopAbs_VERTEX));
      aBoxExpandValue += aMaxVertexTol;
    }
    //
    NCollection_DynamicArray<BOPDS_Curve>& aVNC = aFF.ChangeCurves();
    for (int i = 1; i <= aNbCurves; ++i)
    {
      if (UserBreak(aPSOuter))
      {
        return;
      }
      Bnd_Box               aBox;
      const IntTools_Curve& aIC      = aCvsX(i);
      bool                  bIsValid = IntTools_Tools::CheckCurve(aIC, aBox);
      if (bIsValid)
      {
        BOPDS_Curve& aNC = aVNC.Appended();
        aNC.SetCurve(aIC);
        // make sure that the bounding box has the maximal gap
        aBox.Enlarge(aBoxExpandValue);
        aNC.SetBox(aBox);
        aNC.SetTolerance(std::max(aIC.Tolerance(), aTolFF));
      }
    }
    //
    // Points
    NCollection_DynamicArray<BOPDS_Point>& aVNP = aFF.ChangePoints();
    for (int i = 1; i <= aNbPoints; ++i)
    {
      const IntTools_PntOn2Faces& aPi = aPntsX(i);
      const gp_Pnt&               aP  = aPi.P1().Pnt();
      //
      BOPDS_Point& aNP = aVNP.Appended();
      aNP.SetPnt(aP);
    }
  }
}

//=======================================================================
// function : UpdateSavedTolerance
// purpose  : Updates the saved tolerance of the vertices of the edge
//           with new tolerance of edge
//=======================================================================
static void UpdateSavedTolerance(const BOPDS_PDS&                  theDS,
                                 const int                         theNE,
                                 const double                      theTolNew,
                                 NCollection_DataMap<int, double>& theMVTol)
{
  const NCollection_List<int>&    aSubShapes = theDS->ShapeInfo(theNE).SubShapes();
  NCollection_List<int>::Iterator itSS(aSubShapes);
  for (; itSS.More(); itSS.Next())
  {
    const int nV        = itSS.Value();
    double*   pTolSaved = theMVTol.ChangeSeek(nV);
    if (pTolSaved && *pTolSaved < theTolNew)
    {
      *pTolSaved = theTolNew;
    }
  }
}

//=================================================================================================

void BOPAlgo_PaveFiller::MakeBlocks(const Message_ProgressRange& theRange)
{
  Message_ProgressScope aPSOuter(theRange, nullptr, 4);
  if (myGlue != BOPAlgo_GlueOff)
  {
    return;
  }
  //
  NCollection_DynamicArray<BOPDS_InterfFF>& aFFs  = myDS->InterfFF();
  int                                       aNbFF = aFFs.Length();
  Message_ProgressScope                     aPS(aPSOuter.Next(), "Building section edges", aNbFF);
  if (!aNbFF)
  {
    return;
  }
  //
  bool   bExist, bValid2D;
  int    i, nF1, nF2, aNbC, aNbP, j;
  int    nV1, nV2;
  double aT1, aT2;
  // Main allocator for cross-iteration data (IncAllocator for performance)
  occ::handle<NCollection_BaseAllocator> aAllocator = new NCollection_IncAllocator;
  // Temporary allocator for per-iteration collections that are cleared each iteration.
  // Using separate allocator allows memory reuse via Reset(false) at each iteration.
  // Collections are kept outside the loop to reuse their bucket arrays after Clear().
  occ::handle<NCollection_IncAllocator> aTmpAllocator = new NCollection_IncAllocator;
  // Default allocator for collections that require proper Free() behavior (Remove/UnBind).
  // IncAllocator::Free() is a no-op, so collections with Remove/UnBind must use default.
  occ::handle<NCollection_BaseAllocator> aDefaultAllocator =
    NCollection_BaseAllocator::CommonBaseAllocator();
  NCollection_List<occ::handle<BOPDS_PaveBlock>>::Iterator aItLPB;
  TopoDS_Edge                                              aES;
  occ::handle<BOPDS_PaveBlock>                             aPBOut;
  //
  //-----------------------------------------------------scope f
  //
  // Per-iteration collections (use temporary allocator, cleared and reset each iteration).
  // Kept outside loop to reuse bucket arrays - Clear() preserves allocated buckets.
  NCollection_List<int>                                aLSE(aTmpAllocator);
  NCollection_List<int>                                aLBV(aTmpAllocator);
  NCollection_Map<int>                                 aMVOnIn(100, aTmpAllocator);
  NCollection_Map<int>                                 aMVCommon(100, aTmpAllocator);
  NCollection_Map<int>                                 aMVStick(100, aTmpAllocator);
  NCollection_Map<int>                                 aMVEF(100, aTmpAllocator);
  NCollection_Map<int>                                 aMVBounds(100, aTmpAllocator);
  NCollection_Map<int>                                 aMI(100, aTmpAllocator);
  NCollection_IndexedMap<occ::handle<BOPDS_PaveBlock>> aMPBOnIn(100, aTmpAllocator);
  NCollection_Map<occ::handle<BOPDS_PaveBlock>>        aMPBCommon(100, aTmpAllocator);
  NCollection_DataMap<int, NCollection_List<int>>      aDMBV(100, aTmpAllocator);
  // aMVTol has UnBind() operations - must use default allocator for proper Free()
  NCollection_DataMap<int, double> aMVTol(100, aDefaultAllocator);
  //
  // Cross-iteration collections (use main allocator, persist through entire loop)
  NCollection_Map<occ::handle<BOPDS_PaveBlock>> aMPBAdd(100, aAllocator);
  // aLPB has Remove() operations - must use default allocator for proper Free()
  NCollection_List<occ::handle<BOPDS_PaveBlock>> aLPB(aDefaultAllocator);
  NCollection_IndexedDataMap<TopoDS_Shape, BOPDS_CoupleOfPaveBlocks, TopTools_ShapeMapHasher>
                                                                  aMSCPB(100, aAllocator);
  NCollection_DataMap<TopoDS_Shape, int, TopTools_ShapeMapHasher> aMVI(100, aAllocator);
  NCollection_DataMap<occ::handle<BOPDS_PaveBlock>, NCollection_List<occ::handle<BOPDS_PaveBlock>>>
                                aDMExEdges(100, aAllocator);
  NCollection_DataMap<int, int> aDMNewSD(100, aAllocator);
  // aDMVLV has UnBind() operations - must use default allocator for proper Free()
  NCollection_DataMap<int, NCollection_List<int>>               aDMVLV(100, aDefaultAllocator);
  NCollection_DataMap<int, double>::Iterator                    aItMV;
  NCollection_IndexedMap<occ::handle<BOPDS_PaveBlock>>          aMicroPB(100, aAllocator);
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aVertsOnRejectedPB(100, aAllocator);
  // Map of PaveBlocks with the faces to which it has to be added
  BOPAlgo_DataMapOfPaveBlockListOfInteger aPBFacesMap(100, aAllocator);
  //
  // The vector aFFToRecheck contains indices of potentially problematic Face-Face intersections
  NCollection_DynamicArray<int> aFFToRecheck;
  // aNbFF may be increased while processing this loop, because it is necessary to recheck
  // some of Face-Face intersections to avoid missing section edges
  // aNbFF will be increased to the number of potentially problematic Face-Face intersections
  const int aNbFFPrev = aNbFF;
  // the vertices the edge/face pass put where an edge lying on a face leaves
  // it (a vertex interference of an edge with a face it also lies on)
  NCollection_Map<int> aMVTangent(100, aAllocator);
  {
    NCollection_DynamicArray<BOPDS_InterfEF>& aEFs = myDS->InterfEF();
    NCollection_Map<BOPDS_Pair>               aMEOnF(100, aAllocator);
    for (int k = 0; k < aEFs.Length(); ++k)
    {
      if (aEFs(k).CommonPart().Type() == TopAbs_EDGE)
      {
        int nE, nF;
        aEFs(k).Indices(nE, nF);
        aMEOnF.Add(BOPDS_Pair(nE, nF));
      }
    }
    for (int k = 0; k < aEFs.Length() && !aMEOnF.IsEmpty(); ++k)
    {
      int nE, nF, nV;
      aEFs(k).Indices(nE, nF);
      if (aEFs(k).CommonPart().Type() == TopAbs_VERTEX && aEFs(k).HasIndexNew(nV)
          && aMEOnF.Contains(BOPDS_Pair(nE, nF)))
      {
        myDS->HasShapeSD(nV, nV);
        aMVTangent.Add(nV);
      }
    }
  }
  // the number of curves each pair's intersection gave
  NCollection_DataMap<int, int> aDMNbCInt;
  for (i = 0; i < aNbFF; ++i, aPS.Next())
  {
    if (UserBreak(aPS))
    {
      return;
    }
    // after passing through all of Face-Face intersections it is necessary to return
    // to potentially problematic Face-Face intersections and process them one more time
    const int aCurInd = i < aNbFFPrev ? i : aFFToRecheck[i - aNbFFPrev];
    //
    BOPDS_InterfFF& aFF = aFFs(aCurInd);
    aFF.Indices(nF1, nF2);
    //
    NCollection_DynamicArray<BOPDS_Point>& aVP = aFF.ChangePoints();
    aNbP                                       = aVP.Length();
    NCollection_DynamicArray<BOPDS_Curve>& aVC = aFF.ChangeCurves();
    aNbC                                       = aVC.Length();
    if (!aNbP && !aNbC)
    {
      continue;
    }
    //
    const TopoDS_Face& aF1 = (*(TopoDS_Face*)(&myDS->Shape(nF1)));
    const TopoDS_Face& aF2 = (*(TopoDS_Face*)(&myDS->Shape(nF2)));
    //
    double aTolFF = std::max(BRep_Tool::Tolerance(aF1), BRep_Tool::Tolerance(aF2));
    //
    BOPDS_FaceInfo& aFI1 = myDS->ChangeFaceInfo(nF1);
    BOPDS_FaceInfo& aFI2 = myDS->ChangeFaceInfo(nF2);
    //
    // Clear per-iteration collections and reset allocator for memory reuse
    aMVOnIn.Clear();
    aMVCommon.Clear();
    aMPBOnIn.Clear();
    aMPBCommon.Clear();
    aDMBV.Clear();
    aMVTol.Clear();
    aLSE.Clear();
    aLBV.Clear();
    aMVStick.Clear();
    aMVEF.Clear();
    aMVBounds.Clear();
    aMI.Clear();
    aTmpAllocator->Reset(false);
    //
    myDS->SubShapesOnIn(nF1, nF2, aMVOnIn, aMVCommon, aMPBOnIn, aMPBCommon);
    myDS->SharedEdges(nF1, nF2, aLSE, aTmpAllocator);
    //
    // 1. Treat Points
    for (j = 0; j < aNbP; ++j)
    {
      TopoDS_Vertex            aV;
      BOPDS_CoupleOfPaveBlocks aCPB;
      //
      BOPDS_Point&  aNP = aVP.ChangeValue(j);
      const gp_Pnt& aP  = aNP.Pnt();
      //
      bExist = IsExistingVertex(aP, aTolFF, aMVOnIn);
      if (!bExist)
      {
        BOPTools_AlgoTools::MakeNewVertex(aP, aTolFF, aV);
        //
        aCPB.SetIndexInterf(i);
        aCPB.SetIndex(j);
        aMSCPB.Add(aV, aCPB);
      }
    }
    //
    // 2. Treat Curves
    GetStickVertices(nF1, nF2, aMVStick, aMVEF, aMI);
    //
    for (j = 0; j < aNbC; ++j)
    {
      BOPDS_Curve& aNC = aVC.ChangeValue(j);
      // DEBt
      aNC.InitPaveBlock1();
      //
      // In order to avoid problems connected with
      // extending tolerance of vertex while putting
      // (e.g. see "bugs modalg_6 bug26789_1" test case),
      // all not-common vertices will be checked by
      // BndBoxes before putting. For common-vertices,
      // filtering by BndBoxes is not necessary.
      PutPavesOnCurve(aMVOnIn, aMVCommon, aNC, aMI, aMVEF, aMVTol, aDMVLV);
    }

    // if some E-F vertex was put on a curve due to large E-F intersection range,
    // and it also was put on another curve correctly then remove this vertex from
    // the first curve. Detect such case if the distance to curve exceeds aTolR3D.
    FilterPavesOnCurves(aVC, aMVTol);
    // nor a vertex its stretched tolerance put on a curve that only passes by
    dropStretchedPaves(myDS, myContext, aVC, aMI, myFuzzyValue, aMVTol);

    // the curves of the intersection; those made here to reach a vertex
    // follow, and are there already when the pair is gone through again
    const int* pNbCInt = aDMNbCInt.Seek(aCurInd);
    const int  aNbCInt = pNbCInt ? *pNbCInt : aNbC;
    if (!pNbCInt)
    {
      aDMNbCInt.Bind(aCurInd, aNbC);
    }
    NCollection_Map<int> aMCLost;

    // the vertices where an edge on a face leaves it, which a curve passes near
    if (!aMVTangent.IsEmpty())
    {
      NCollection_Map<int> aMVT;
      for (const NCollection_Map<int>* pM : {&aMVCommon, &aMVStick, &aMVEF})
      {
        for (NCollection_Map<int>::Iterator aItS(*pM); aItS.More(); aItS.Next())
        {
          int nVSD = aItS.Value();
          myDS->HasShapeSD(nVSD, nVSD);
          if (aMVTangent.Contains(nVSD))
          {
            aMVT.Add(aItS.Value());
          }
        }
      }
      NCollection_List<int> aLEShared;
      if (!aMVT.IsEmpty())
      {
        sharedEdges(myDS, aMPBOnIn, aMPBCommon, aLEShared);
      }
      NCollection_DataMap<int, int> aMVOwner; // the curve that reaches a vertex
      // the curve cut short of where it loses its way near the vertex, and a
      // curve of its own from there to the vertex
      const auto aJoin = [&](const int theJ, const int theNV, const double theTolV) -> bool {
        occ::handle<Geom_Curve>   aBridge;
        occ::handle<Geom2d_Curve> aBridge2d[2];
        double                    aTolB = 0.;
        if (!bridgeToVertex(myDS,
                            myContext,
                            aF1,
                            aF2,
                            aVC,
                            theJ,
                            theNV,
                            theTolV,
                            ToleranceFF(myContext->SurfaceAdaptor(aF1),
                                        myContext->SurfaceAdaptor(aF2)),
                            aLEShared,
                            aBridge,
                            aBridge2d,
                            aTolB,
                            aMVOwner))
        {
          return false;
        }
        if (!aBridge.IsNull())
        {
          const IntTools_Curve aICB(aBridge, aBridge2d[0], aBridge2d[1], aTolB, 0.);
          Bnd_Box              aBoxB;
          if (IntTools_Tools::CheckCurve(aICB, aBoxB))
          {
            BOPDS_Curve& aNB = aVC.Appended();
            aNB.SetCurve(aICB);
            aBoxB.Enlarge(aTolB + theTolV + myFuzzyValue);
            aNB.SetBox(aBoxB);
            aNB.SetTolerance(aTolB);
            aNB.InitPaveBlock1();
            PutPaveOnCurve(theNV, aTolB, aNB, aMI, aMVTol, aDMVLV);
          }
        }
        return true;
      };
      NCollection_DynamicArray<int>    aLaterJ, aLaterV;
      NCollection_DynamicArray<double> aLaterD;
      for (j = 0; j < aNbCInt && !aMVT.IsEmpty(); ++j)
      {
        BOPDS_Curve& aNC = aVC.ChangeValue(j);
        // not a curve along an edge the faces share already: its pieces go to
        // that edge, and a vertex grown to reach its end would stay grown
        if (curveAlongEdges(myDS, myContext, aNC, aMPBOnIn, aMPBCommon))
        {
          // and where it ends on one of them it takes no vertex of its own
          if (trimEndsOnEdges(myDS, myContext, aF1, aF2, aNC, aLEShared))
          {
            aMCLost.Add(j);
          }
          continue;
        }
        // nor one lost beside such an edge: it takes no vertex at all
        if (lostAlongEdges(myDS, myContext, aF1, aF2, aNC, aLEShared))
        {
          aMCLost.Add(j);
          continue;
        }
        NCollection_List<int>    aLV;
        NCollection_List<double> aLD;
        tangentVerticesNearCurve(myDS, myContext, aF1, aF2, aNC, aMVT, aMVTol, aLV, aLD);
        NCollection_List<double>::Iterator aItD(aLD);
        for (NCollection_List<int>::Iterator aItV(aLV); aItV.More(); aItV.Next(), aItD.Next())
        {
          const int    nV    = aItV.Value();
          const double aTolV = aMVTol.IsBound(nV)
                                 ? aMVTol(nV)
                                 : BRep_Tool::Tolerance(TopoDS::Vertex(myDS->Shape(nV)));
          if (!aJoin(j, nV, aTolV))
          {
            aLaterJ.Append(j);
            aLaterV.Append(nV);
            aLaterD.Append(aItD.Value());
          }
        }
      }
      // a line the vertex was not joined to is tried once more if another
      // line owns the vertex by now; else the vertex is grown to the line
      for (int k = 0; k < aLaterJ.Length(); ++k)
      {
        const int    nV    = aLaterV(k);
        const double aTolV = aMVTol.IsBound(nV)
                               ? aMVTol(nV)
                               : BRep_Tool::Tolerance(TopoDS::Vertex(myDS->Shape(nV)));
        if (!aMVOwner.IsBound(nV) || !aJoin(aLaterJ(k), nV, aTolV))
        {
          PutPaveOnCurve(nV,
                         aLaterD(k) * (1. + 1.e-6),
                         aVC.ChangeValue(aLaterJ(k)),
                         aMI,
                         aMVTol,
                         aDMVLV);
        }
      }
      aNbC = aVC.Length();
    }

    for (j = 0; j < aNbC; ++j)
    {
      BOPDS_Curve&          aNC = aVC.ChangeValue(j);
      const IntTools_Curve& aIC = aNC.Curve();
      //
      PutStickPavesOnCurve(aF1, aF2, aMI, aVC, j, aMVStick, aMVTol, aDMVLV);
      // 904/F7
      if (aNbCInt == 1)
      {
        PutEFPavesOnCurve(aVC, j, aMI, aMVEF, aMVTol, aDMVLV);
      }
      //
      if (aIC.HasBounds() && !aMCLost.Contains(j))
      {
        aLBV.Clear();
        //
        // an end short of a vertex of both faces that the section runs on to
        {
          int aEndNV[2];
          getBoundPaves(myDS, aNC, aEndNV);
          for (int k = 0; k < 2; ++k)
          {
            if (aEndNV[k] >= 0)
            {
              continue;
            }
            NCollection_Map<int> aMV;
            aMV.Assign(aMVCommon);
            for (NCollection_Map<int>::Iterator aItS(aMVStick); aItS.More(); aItS.Next())
            {
              aMV.Add(aItS.Value());
            }
            RemoveUsedVertices(aVC, aMV);
            double    aD;
            const int nV = continuedEndVertex(myDS, myContext, aF1, aF2, aNC, k, aMV, aD);
            if (nV >= 0)
            {
              PutPaveOnCurve(nV, aD * (1. + 1.e-6), aNC, aMI, aMVTol, aDMVLV);
            }
          }
        }
        //
        PutBoundPaveOnCurve(aF1, aF2, aNC, aLBV);
        //
        // an end outside a face, at a corner where the faces are tangent
        int aBndNV[2];
        getBoundPaves(myDS, aNC, aBndNV);
        for (int k = 0; k < 2; ++k)
        {
          if (aBndNV[k] >= 0)
          {
            continue;
          }
          NCollection_Map<int> aMV;
          aMV.Assign(aMVCommon);
          for (NCollection_Map<int>::Iterator aItS(aMVStick); aItS.More(); aItS.Next())
          {
            aMV.Add(aItS.Value());
          }
          RemoveUsedVertices(aVC, aMV);
          double    aD;
          const int nV = tangentEndVertex(myDS, myContext, aF1, aF2, aNC, k, aMV, aD);
          if (nV >= 0)
          {
            PutPaveOnCurve(nV, aD, aNC, aMI, aMVTol, aDMVLV);
          }
        }
        //
        if (!aLBV.IsEmpty())
        {
          aDMBV.Bind(j, aLBV);
          NCollection_List<int>::Iterator aItI(aLBV);
          for (; aItI.More(); aItI.Next())
          {
            aMVBounds.Add(aItI.Value());
          }
        }
      }
    } // for (j=0; j<aNbC; ++j) {

    // Put closing pave if needed
    for (j = 0; j < aNbC; ++j)
    {
      BOPDS_Curve& aNC = aVC.ChangeValue(j);
      PutClosingPaveOnCurve(aNC);
    }
    //

    BOPTools_BoxTree aPBTree;
    {
      // Fill the tree with boxes of pave blocks ON/IN
      // Tree will be build on first selection from the tree.
      const int aNbPB = aMPBOnIn.Extent();
      aPBTree.SetSize(aNbPB);
      for (int iPB = 1; iPB <= aNbPB; ++iPB)
      {
        const occ::handle<BOPDS_PaveBlock>& aPB = aMPBOnIn(iPB);
        if (!aPB->HasEdge())
        {
          continue;
        }

        if (myDS->ShapeInfo(aPB->OriginalEdge()).HasFlag())
        {
          continue;
        }

        aPBTree.Add(iPB, Bnd_Tools::Bnd2BVH(myDS->ShapeInfo(aPB->Edge()).Box()));
      }
    }

    // Added additional check of Face-Face intersection to avoid missing section edges
    // because of sequence of Face-Face interference processing
    bool isToRecheck = aNbC > 0 && i < aNbFFPrev;
    //
    // 3. Make section edges
    for (j = 0; j < aNbC; ++j)
    {
      BOPDS_Curve&          aNC     = aVC.ChangeValue(j);
      const IntTools_Curve& aIC     = aNC.Curve();
      double                aTolR3D = std::max(aNC.Tolerance(), aNC.TangentialTolerance());
      //
      NCollection_List<occ::handle<BOPDS_PaveBlock>>& aLPBC = aNC.ChangePaveBlocks();
      occ::handle<BOPDS_PaveBlock>&                   aPB1  = aNC.ChangePaveBlock1();
      //
      aLPB.Clear();
      aPB1->Update(aLPB, false);
      //
      if (aLPB.Extent() != 0)
      {
        isToRecheck = false;
      }

      aItLPB.Initialize(aLPB);
      for (; aItLPB.More(); aItLPB.Next())
      {
        occ::handle<BOPDS_PaveBlock>& aPB = aItLPB.ChangeValue();
        aPB->Indices(nV1, nV2);
        aPB->Range(aT1, aT2);
        //
        if (fabs(aT1 - aT2) < Precision::PConfusion())
        {
          continue;
        }

        // Check validity of the block for the faces:
        // classify bounding and middle points on the curve
        // relatively both faces
        bValid2D = myContext->IsValidBlockForFaces(aT1, aT2, aIC, aF1, aF2, aTolR3D);
        if (!bValid2D)
        {
          continue;
        }
        //
        int    nEOut;
        double aTolNew = -1.;
        bExist         = IsExistingPaveBlock(aPB, aNC, aLSE, nEOut, aTolNew);
        if (bExist)
        {
          // Update edge with new tolerance
          UpdateEdgeTolerance(nEOut, aTolNew);
          // Update aMVTol map with new tolerances of vertices
          UpdateSavedTolerance(myDS, nEOut, aTolNew, aMVTol);
          continue;
        }
        //
        const TopoDS_Vertex& aV1 = (*(TopoDS_Vertex*)(&myDS->Shape(nV1)));
        const TopoDS_Vertex& aV2 = (*(TopoDS_Vertex*)(&myDS->Shape(nV2)));
        //
        // check if the pave block has a valid range
        double aFirst, aLast;
        if (!BRepLib::FindValidRange(GeomAdaptor_Curve(aIC.Curve()),
                                     aTolR3D,
                                     aT1,
                                     BRep_Tool::Pnt(aV1),
                                     std::max(aTolR3D, BRep_Tool::Tolerance(aV1)),
                                     aT2,
                                     BRep_Tool::Pnt(aV2),
                                     std::max(aTolR3D, BRep_Tool::Tolerance(aV2)),
                                     aFirst,
                                     aLast))
        {
          // If the pave block does not have valid range, i.e. it is completely
          // covered by the tolerance spheres of its vertices, it will be
          // passed into post treatment process to fuse its vertices.
          // The pave block itself will not be kept.
          if (!aMVBounds.Contains(nV1) && !aMVBounds.Contains(nV2))
          {
            aMicroPB.Add(aPB);
            // keep vertices for post treatment
            aMVI.Bind(aV1, nV1);
            aMVI.Bind(aV2, nV2);
          }
          continue;
        }
        //
        bExist =
          IsExistingPaveBlock(aPB, aNC, aTolR3D, aMPBOnIn, aPBTree, aMPBCommon, aPBOut, aTolNew);
        if (bExist)
        {
          bool bInF1 =
            (aFI1.PaveBlocksOn().Contains(aPBOut) || aFI1.PaveBlocksIn().Contains(aPBOut));
          bool bInF2 =
            (aFI2.PaveBlocksOn().Contains(aPBOut) || aFI2.PaveBlocksIn().Contains(aPBOut));
          if (!bInF1 || !bInF2)
          {
            // Update edge to touch both faces
            int                nE    = aPBOut->Edge();
            const TopoDS_Edge& aE    = *(TopoDS_Edge*)&myDS->Shape(nE);
            double             aTolE = BRep_Tool::Tolerance(aE);
            if (aTolNew < aNC.Tolerance())
            {
              aTolNew = aNC.Tolerance(); // use real tolerance of intersection
            }
            if (aTolNew > aTolE)
            {
              UpdateEdgeTolerance(nE, aTolNew);
              // Update aMVTol map with new tolerances of vertices
              UpdateSavedTolerance(myDS, nE, aTolNew, aMVTol);
            }

            // Face without pave block
            const int              nF     = bInF1 ? nF2 : nF1;
            NCollection_List<int>* pFaces = aPBFacesMap.ChangeSeek(aPBOut);
            if (!pFaces)
            {
              pFaces = aPBFacesMap.Bound(aPBOut, NCollection_List<int>());
            }
            // List is expected to be short, so we allow the check here
            if (pFaces->IsEmpty() || !pFaces->Contains(nF))
            {
              pFaces->Append(nF);
            }

            // Try fusing the vertices of the existing pave block
            // with the vertices put on the real section curve (except
            // for technological vertices, which will be removed)
            int nVOut1, nVOut2;
            aPBOut->Indices(nVOut1, nVOut2);
            if (nV1 != nVOut1 && nV1 != nVOut2 && !aMVBounds.Contains(nV1))
            {
              aVertsOnRejectedPB.Add(aV1);
            }
            if (nV2 != nVOut1 && nV2 != nVOut2 && !aMVBounds.Contains(nV2))
            {
              aVertsOnRejectedPB.Add(aV2);
            }

            if (aMPBAdd.Add(aPBOut))
            {
              // Add edge for processing as the section edge
              PreparePostTreatFF(aCurInd, j, aPBOut, aMSCPB, aMVI, aLPBC);
            }
          }
          continue;
        }
        //
        // Make Edge
        BOPTools_AlgoTools::MakeEdge(aIC, aV1, aT1, aV2, aT2, aTolR3D, aES);
        // Make p-curves
        BOPTools_AlgoTools::MakePCurve(aES,
                                       aF1,
                                       aF2,
                                       aIC,
                                       mySectionAttribute.PCurveOnS1(),
                                       mySectionAttribute.PCurveOnS2(),
                                       myContext);
        //
        // Append the Pave Block to the Curve j
        aLPBC.Append(aPB);
        //
        // Keep info for post treatment
        BOPDS_CoupleOfPaveBlocks aCPB;
        aCPB.SetIndexInterf(aCurInd);
        aCPB.SetIndex(j);
        aCPB.SetPaveBlock1(aPB);
        //
        aMSCPB.Add(aES, aCPB);
        aMVI.Bind(aV1, nV1);
        aMVI.Bind(aV2, nV2);
        //
        aMVTol.UnBind(nV1);
        aMVTol.UnBind(nV2);

        // Add existing pave blocks for post treatment
        ProcessExistingPaveBlocks(aCurInd,
                                  j,
                                  nF1,
                                  nF2,
                                  aES,
                                  aMPBOnIn,
                                  aPBTree,
                                  aMSCPB,
                                  aMVI,
                                  aLPBC,
                                  aPBFacesMap,
                                  aMPBAdd);
      }
      //
      aLPBC.RemoveFirst();
    } // for (j=0; j<aNbC; ++j) {
    if (isToRecheck)
    {
      aFFToRecheck.Append(aCurInd);
      ++aNbFF;
    }
    //
    // back to previous tolerance values for unused vertices
    // and forget about SD groups of such vertices
    aItMV.Initialize(aMVTol);
    for (; aItMV.More(); aItMV.Next())
    {
      nV1         = aItMV.Key();
      double aTol = aItMV.Value();
      //
      const TopoDS_Vertex&             aV = *(TopoDS_Vertex*)&myDS->Shape(nV1);
      const occ::handle<BRep_TVertex>& TV = *((occ::handle<BRep_TVertex>*)&aV.TShape());
      TV->Tolerance(aTol);
      // reset bnd box
      BOPDS_ShapeInfo& aSIDS  = myDS->ChangeShapeInfo(nV1);
      Bnd_Box&         aBoxDS = aSIDS.ChangeBox();
      aBoxDS                  = Bnd_Box();
      BRepBndLib::Add(aV, aBoxDS);
      aBoxDS.SetGap(aBoxDS.GetGap() + Precision::Confusion());
      //
      if (aDMVLV.IsBound(nV1))
      {
        aDMVLV.UnBind(nV1);
      }
    }
    //
    ProcessExistingPaveBlocks(aCurInd,
                              nF1,
                              nF2,
                              aMPBOnIn,
                              aPBTree,
                              aDMBV,
                              aMSCPB,
                              aMVI,
                              aPBFacesMap,
                              aMPBAdd);
  } // for (i=0; i<aNbFF; ++i) {

  // Remove "micro" section edges
  RemoveMicroSectionEdges(aMSCPB, aMicroPB);

  // post treatment
  MakeSDVerticesFF(aDMVLV, aDMNewSD);
  PostTreatFF(aMSCPB,
              aDMExEdges,
              aDMNewSD,
              aMicroPB,
              aVertsOnRejectedPB,
              aAllocator,
              aPSOuter.Next(2));
  if (HasErrors())
  {
    return;
  }
  // reduce tolerances of section edges where it is appropriate
  CorrectToleranceOfSE();
  //
  // update face info
  UpdateFaceInfo(aDMExEdges, aDMNewSD, aPBFacesMap);
  // Update all pave blocks
  UpdatePaveBlocks(aDMNewSD);
  //
  // Treat possible common zones by trying to put each section edge
  // into all faces, not participated in creation of that edge, as IN edge

  PutSEInOtherFaces(aPSOuter.Next());
}

//=================================================================================================

void BOPAlgo_PaveFiller::MakeSDVerticesFF(
  const NCollection_DataMap<int, NCollection_List<int>>& theDMVLV,
  NCollection_DataMap<int, int>&                         theDMNewSD)
{
  // Create a new SD vertex for each group of coinciding vertices
  // and put new substitutions to theDMNewSD.
  NCollection_DataMap<int, NCollection_List<int>>::Iterator aItG(theDMVLV);
  for (; aItG.More(); aItG.Next())
  {
    const NCollection_List<int>& aList = aItG.Value();
    // make SD vertices w/o creation of interfs
    int nSD = MakeSDVertices(aList, false);
    // update theDMNewSD
    NCollection_List<int>::Iterator aItL(aList);
    for (; aItL.More(); aItL.Next())
    {
      int nV = aItL.Value();
      theDMNewSD.Bind(nV, nSD);
    }
  }
}

//=================================================================================================

void BOPAlgo_PaveFiller::PostTreatFF(
  NCollection_IndexedDataMap<TopoDS_Shape, BOPDS_CoupleOfPaveBlocks, TopTools_ShapeMapHasher>&
    theMSCPB,
  NCollection_DataMap<occ::handle<BOPDS_PaveBlock>, NCollection_List<occ::handle<BOPDS_PaveBlock>>>&
                                                                       aDMExEdges,
  NCollection_DataMap<int, int>&                                       aDMNewSD,
  const NCollection_IndexedMap<occ::handle<BOPDS_PaveBlock>>&          theMicroPB,
  const NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>& theVertsOnRejectedPB,
  const occ::handle<NCollection_BaseAllocator>&                        theAllocator,
  const Message_ProgressRange&                                         theRange)
{
  int aNbS = theMSCPB.Extent();
  if (!aNbS)
  {
    return;
  }
  //
  bool                                     bHasPaveBlocks, bOld;
  int                                      nSx, nVSD, iX, iP, iC, j, nV, iV = 0, iE, k;
  int                                      aNbLPBx;
  TopAbs_ShapeEnum                         aType;
  TopoDS_Shape                             aV, aE;
  NCollection_List<TopoDS_Shape>::Iterator aItLS;
  NCollection_List<occ::handle<BOPDS_PaveBlock>>::Iterator aItLPB;
  BOPDS_PDS                                                aPDS;
  occ::handle<BOPDS_PaveBlock>                             aPB1;
  BOPDS_Pave                                               aPave[2];
  BOPDS_ShapeInfo                                          aSI;
  //
  NCollection_List<TopoDS_Shape> aLS(theAllocator);
  BOPAlgo_PaveFiller             aPF(theAllocator);
  aPF.SetIsPrimary(false);
  aPF.SetNonDestructive(myNonDestructive);
  //
  NCollection_DynamicArray<BOPDS_InterfFF>& aFFs  = myDS->InterfFF();
  int                                       aNbFF = aFFs.Length();
  //

  // Find unused vertices
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> VertsUnused;
  NCollection_Map<int>                                          IndMap;
  for (int i = 0; i < aNbFF; i++)
  {
    BOPDS_InterfFF& aFF = aFFs(i);
    int             nF1, nF2;
    aFF.Indices(nF1, nF2);

    NCollection_Map<int> aMV, aMVEF, aMI;
    GetStickVertices(nF1, nF2, aMV, aMVEF, aMI);
    NCollection_DynamicArray<BOPDS_Curve>& aVC = aFF.ChangeCurves();
    RemoveUsedVertices(aVC, aMV);

    NCollection_Map<int>::Iterator itmap(aMV);
    for (; itmap.More(); itmap.Next())
    {
      int                 indV    = itmap.Value();
      const TopoDS_Shape& aVertex = myDS->Shape(indV);
      if (IndMap.Add(indV))
      {
        VertsUnused.Add(aVertex);
      }
      else
      {
        VertsUnused.RemoveKey(aVertex);
      }
    }
  }
  /////////////////////

  int aNbME     = theMicroPB.Extent();
  int aNbVOnRPB = theVertsOnRejectedPB.Extent();
  // 0
  if (aNbS == 1 && (aNbME == 0) && (aNbVOnRPB == 0) && VertsUnused.IsEmpty())
  {
    const TopoDS_Shape&             aS   = theMSCPB.FindKey(1);
    const BOPDS_CoupleOfPaveBlocks& aCPB = theMSCPB.FindFromIndex(1);
    //
    aType = aS.ShapeType();
    if (aType == TopAbs_VERTEX)
    {
      aSI.SetShapeType(aType);
      aSI.SetShape(aS);
      iV = myDS->Append(aSI);
      //
      iX                                          = aCPB.IndexInterf();
      iP                                          = aCPB.Index();
      BOPDS_InterfFF&                        aFF  = aFFs(iX);
      NCollection_DynamicArray<BOPDS_Point>& aVNP = aFF.ChangePoints();
      BOPDS_Point&                           aNP  = aVNP(iP);
      aNP.SetIndex(iV);
    }
    else if (aType == TopAbs_EDGE)
    {
      aPB1 = aCPB.PaveBlock1();
      //
      if (aPB1->HasEdge())
      {
        NCollection_List<occ::handle<BOPDS_PaveBlock>> aLPBx;
        aLPBx.Append(aPB1);
        aDMExEdges.Bind(aPB1, aLPBx);
      }
      else
      {
        aSI.SetShapeType(aType);
        aSI.SetShape(aS);
        iE = myDS->Append(aSI);
        //
        aPB1->SetEdge(iE);
      }
    }
    return;
  }
  //
  // 1 prepare arguments
  // Avoid intersection of existing edges among themselves
  TopoDS_Compound anExistingEdges;
  BRep_Builder().MakeCompound(anExistingEdges);
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> anAddedSD;
  for (k = aNbS; k > 0; --k)
  {
    const TopoDS_Shape&                 aS  = theMSCPB.FindKey(k);
    const occ::handle<BOPDS_PaveBlock>& aPB = theMSCPB(k).PaveBlock1();
    if (!aPB.IsNull() && aPB->HasEdge())
    {
      BRep_Builder().Add(anExistingEdges, aS);
    }
    else
    {
      aLS.Append(aS);
    }
    // add vertices-candidates for SD from the map aDMNewSD,
    // so that they took part in fuse operation.
    TopoDS_Iterator itV(aS);
    for (; itV.More(); itV.Next())
    {
      const TopoDS_Shape& aVer = itV.Value();
      int                 iVer = myDS->Index(aVer);
      const int*          pSD  = aDMNewSD.Seek(iVer);
      if (pSD)
      {
        const TopoDS_Shape& aVSD = myDS->Shape(*pSD);
        if (anAddedSD.Add(aVSD))
        {
          aLS.Append(aVSD);
        }
      }
    }
  }
  if (anExistingEdges.NbChildren() > 0)
  {
    aLS.Append(anExistingEdges);
  }
  //
  // The section edges considered as a micro should be
  // specially treated - their vertices should be united and
  // the edge itself should be removed. Thus, we add only
  // its vertices into operation.
  //
  BRep_Builder aBB;
  for (k = 1; k <= aNbME; ++k)
  {
    int nVerts[2];
    theMicroPB(k)->Indices(nVerts[0], nVerts[1]);
    TopoDS_Vertex aVerts[2];
    for (int i = 0; i < 2; ++i)
    {
      const int* pSD = aDMNewSD.Seek(nVerts[i]);
      aVerts[i]      = TopoDS::Vertex(myDS->Shape(pSD ? *pSD : nVerts[i]));
      if (anAddedSD.Add(aVerts[i]))
      {
        aLS.Append(aVerts[i]);
      }
    }
    //
    if (aVerts[0].IsSame(aVerts[1]))
    {
      continue;
    }
    //
    // make sure these vertices will be united
    const gp_Pnt& aP1 = BRep_Tool::Pnt(aVerts[0]);
    const gp_Pnt& aP2 = BRep_Tool::Pnt(aVerts[1]);
    //
    double aDist  = aP1.Distance(aP2);
    double aTolV1 = BRep_Tool::Tolerance(aVerts[0]);
    double aTolV2 = BRep_Tool::Tolerance(aVerts[1]);
    //
    aDist -= (aTolV1 + aTolV2);
    if (aDist > 0.)
    {
      aDist /= 2.;
      aBB.UpdateVertex(aVerts[0], aTolV1 + aDist);
      aBB.UpdateVertex(aVerts[1], aTolV2 + aDist);
    }
  }

  // Add vertices put on the real section curves to unify them with the
  // vertices of the edges, by which these sections curves have been rejected
  // and with unused vertices
  const NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>* VerMap[2] = {
    &theVertsOnRejectedPB,
    &VertsUnused};
  for (int imap = 0; imap < 2; imap++)
  {
    int NbVer = VerMap[imap]->Extent();
    for (int i = 1; i <= NbVer; ++i)
    {
      TopoDS_Shape aVer = VerMap[imap]->FindKey(i);
      int          iVer = myDS->Index(aVer);
      const int*   pSD  = aDMNewSD.Seek(iVer);
      if (pSD)
      {
        aVer = myDS->Shape(*pSD);
      }

      if (anAddedSD.Add(aVer))
      {
        aLS.Append(aVer);
      }
    }
  }
  //
  Message_ProgressScope aPS(theRange, "Intersection of section edges", 1);

  // 2 Fuse shapes
  aPF.SetRunParallel(myRunParallel);
  aPF.SetArguments(aLS);
  aPF.Perform(aPS.Next());
  if (aPF.HasErrors())
  {
    AddError(new BOPAlgo_AlertPostTreatFF);
    return;
  }
  aPDS = aPF.PDS();
  //
  // Map to store the real tolerance of the common block
  // and avoid repeated computation of it
  NCollection_DataMap<occ::handle<BOPDS_CommonBlock>, double> aMCBTol;
  // Map to avoid creation of different pave blocks for
  // the same intersection edge
  NCollection_DataMap<int, occ::handle<BOPDS_PaveBlock>> aMEPB;
  //
  aItLS.Initialize(aLS);
  for (; aItLS.More(); aItLS.Next())
  {
    const TopoDS_Shape& aSx = aItLS.Value();
    if (aSx.ShapeType() == TopAbs_COMPOUND)
    {
      for (TopoDS_Iterator itC(aSx); itC.More(); itC.Next())
      {
        aLS.Append(itC.Value());
      }
      continue;
    }
    nSx                         = aPDS->Index(aSx);
    const BOPDS_ShapeInfo& aSIx = aPDS->ShapeInfo(nSx);
    //
    aType = aSIx.ShapeType();
    //
    if (aType == TopAbs_VERTEX)
    {
      bool bIntersectionPoint = theMSCPB.Contains(aSx);
      //
      if (aPDS->HasShapeSD(nSx, nVSD))
      {
        aV = aPDS->Shape(nVSD);
      }
      else
      {
        aV = aSx;
      }
      // index of new vertex in theDS -> iV
      iV = myDS->Index(aV);
      if (iV < 0)
      {
        aSI.SetShapeType(aType);
        aSI.SetShape(aV);
        iV = myDS->Append(aSI);
      }
      //
      if (!bIntersectionPoint)
      {
        // save SD connection
        nSx = myDS->Index(aSx);
        if (nSx != iV)
        {
          aDMNewSD.Bind(nSx, iV);
          myDS->AddShapeSD(nSx, iV);
        }
      }
      else
      {
        // update FF interference
        const BOPDS_CoupleOfPaveBlocks& aCPB        = theMSCPB.FindFromKey(aSx);
        iX                                          = aCPB.IndexInterf();
        iP                                          = aCPB.Index();
        BOPDS_InterfFF&                        aFF  = aFFs(iX);
        NCollection_DynamicArray<BOPDS_Point>& aVNP = aFF.ChangePoints();
        BOPDS_Point&                           aNP  = aVNP(iP);
        aNP.SetIndex(iV);
      }
    } // if (aType==TopAbs_VERTEX) {
    //
    else if (aType == TopAbs_EDGE)
    {
      bHasPaveBlocks                       = aPDS->HasPaveBlocks(nSx);
      const BOPDS_CoupleOfPaveBlocks& aCPB = theMSCPB.FindFromKey(aSx);
      iX                                   = aCPB.IndexInterf();
      iC                                   = aCPB.Index();
      aPB1                                 = aCPB.PaveBlock1();
      //
      bOld = aPB1->HasEdge();
      if (bOld)
      {
        NCollection_List<occ::handle<BOPDS_PaveBlock>> aLPBx;
        aDMExEdges.Bind(aPB1, aLPBx);
      }
      //
      if (!bHasPaveBlocks)
      {
        if (bOld)
        {
          aDMExEdges.ChangeFind(aPB1).Append(aPB1);
        }
        else
        {
          aSI.SetShapeType(aType);
          aSI.SetShape(aSx);
          iE = myDS->Append(aSI);
          //
          aPB1->SetEdge(iE);
        }
      }
      else
      {
        BOPDS_InterfFF&                                 aFF   = aFFs(iX);
        NCollection_DynamicArray<BOPDS_Curve>&          aVNC  = aFF.ChangeCurves();
        BOPDS_Curve&                                    aNC   = aVNC(iC);
        NCollection_List<occ::handle<BOPDS_PaveBlock>>& aLPBC = aNC.ChangePaveBlocks();
        //
        // check if edge occurred to be micro edge;
        // note we check not the edge aSx itself, but its image in aPDS
        const NCollection_List<occ::handle<BOPDS_PaveBlock>>& aLPBx = aPDS->PaveBlocks(nSx);
        aNbLPBx                                                     = aLPBx.Extent();
        if (aNbLPBx == 0 || (aNbLPBx == 1 && !aLPBx.First()->HasShrunkData()))
        {
          NCollection_List<occ::handle<BOPDS_PaveBlock>>::Iterator it(aLPBC);
          for (; it.More(); it.Next())
          {
            if (it.Value() == aPB1)
            {
              aLPBC.Remove(it);
              break;
            }
          }

          // The edge became micro edge, check vertices for SD
          TopoDS_Iterator itV(aSx);
          for (; itV.More(); itV.Next())
          {
            aLS.Append(itV.Value());
          }

          continue;
        }
        //
        if (bOld && !aNbLPBx)
        {
          aDMExEdges.ChangeFind(aPB1).Append(aPB1);
          continue;
        }
        //
        if (!bOld)
        {
          aItLPB.Initialize(aLPBC);
          for (; aItLPB.More(); aItLPB.Next())
          {
            const occ::handle<BOPDS_PaveBlock>& aPBC = aItLPB.Value();
            if (aPBC == aPB1)
            {
              aLPBC.Remove(aItLPB);
              break;
            }
          }
        }
        //
        if (aNbLPBx)
        {
          aItLPB.Initialize(aLPBx);
          for (; aItLPB.More(); aItLPB.Next())
          {
            const occ::handle<BOPDS_PaveBlock>& aPBx  = aItLPB.Value();
            const occ::handle<BOPDS_PaveBlock>  aPBRx = aPDS->RealPaveBlock(aPBx);
            //
            // update vertices of paves
            aPave[0] = aPBx->Pave1();
            aPave[1] = aPBx->Pave2();
            for (j = 0; j < 2; ++j)
            {
              nV = aPave[j].Index();
              aV = aPDS->Shape(nV);
              // index of new vertex in myDS -> iV
              iV = myDS->Index(aV);
              if (iV < 0)
              {
                aSI.SetShapeType(TopAbs_VERTEX);
                aSI.SetShape(aV);
                iV = myDS->Append(aSI);
              }
              const BOPDS_Pave& aP1 = !j ? aPB1->Pave1() : aPB1->Pave2();
              if (aP1.Index() != iV)
              {
                if (aP1.Parameter() == aPave[j].Parameter())
                {
                  aDMNewSD.Bind(aP1.Index(), iV);
                  myDS->AddShapeSD(aP1.Index(), iV);
                }
                else
                {
                  // check aPDS to have the SD connection between these vertices
                  const TopoDS_Shape& aVPave = myDS->Shape(aP1.Index());
                  int                 nVnewSD, nVnew = aPDS->Index(aVPave);
                  if (aPDS->HasShapeSD(nVnew, nVnewSD))
                  {
                    if (nVnewSD == nV)
                    {
                      aDMNewSD.Bind(aP1.Index(), iV);
                      myDS->AddShapeSD(aP1.Index(), iV);
                    }
                  }
                }
              }
              //
              aPave[j].SetIndex(iV);
            }
            //
            // add edge
            aE = aPDS->Shape(aPBRx->Edge());
            iE = myDS->Index(aE);
            if (iE < 0)
            {
              aSI.SetShapeType(aType);
              aSI.SetShape(aE);
              iE = myDS->Append(aSI);
            }
            //
            // update real edge tolerance according to distances in common block if any
            if (aPDS->IsCommonBlock(aPBRx))
            {
              const occ::handle<BOPDS_CommonBlock>& aCB  = aPDS->CommonBlock(aPBRx);
              double*                               pTol = aMCBTol.ChangeSeek(aCB);
              if (!pTol)
              {
                double aTol = BOPAlgo_Tools::ComputeToleranceOfCB(aCB, aPDS, aPF.Context());
                pTol        = aMCBTol.Bound(aCB, aTol);
              }
              //
              if (aNC.Tolerance() < *pTol)
              {
                aNC.SetTolerance(*pTol);
              }
            }
            // append new PaveBlock to aLPBC
            occ::handle<BOPDS_PaveBlock>* pPBC = aMEPB.ChangeSeek(iE);
            if (!pPBC)
            {
              pPBC = aMEPB.Bound(iE, new BOPDS_PaveBlock());
              BOPDS_Pave aPaveR1, aPaveR2;
              aPaveR1 = aPBRx->Pave1();
              aPaveR2 = aPBRx->Pave2();
              aPaveR1.SetIndex(myDS->Index(aPDS->Shape(aPaveR1.Index())));
              aPaveR2.SetIndex(myDS->Index(aPDS->Shape(aPaveR2.Index())));
              //
              (*pPBC)->SetPave1(aPaveR1);
              (*pPBC)->SetPave2(aPaveR2);
              (*pPBC)->SetEdge(iE);
            }
            //
            if (bOld)
            {
              (*pPBC)->SetOriginalEdge(aPB1->OriginalEdge());
              aDMExEdges.ChangeFind(aPB1).Append(*pPBC);
            }
            else
            {
              aLPBC.Append(*pPBC);
            }
          }
        }
      }
    } // else if (aType==TopAbs_EDGE)
  } // for (; aItLS.More(); aItLS.Next()) {

  // Update SD for vertices that did not participate in operation
  NCollection_DataMap<int, int>::Iterator itDM(aDMNewSD);
  for (; itDM.More(); itDM.Next())
  {
    const int* pSD = aDMNewSD.Seek(itDM.Value());
    if (pSD)
    {
      itDM.ChangeValue() = *pSD;
      myDS->AddShapeSD(itDM.Key(), *pSD);
    }
  }
}

//=================================================================================================

void BOPAlgo_PaveFiller::UpdateFaceInfo(
  NCollection_DataMap<occ::handle<BOPDS_PaveBlock>, NCollection_List<occ::handle<BOPDS_PaveBlock>>>&
                                                 theDME,
  const NCollection_DataMap<int, int>&           theDMV,
  const BOPAlgo_DataMapOfPaveBlockListOfInteger& thePBFacesMap)
{
  int                                                      i, j, nV1, nF1, nF2, aNbFF, aNbC, aNbP;
  NCollection_List<occ::handle<BOPDS_PaveBlock>>::Iterator aItLPB;
  NCollection_Map<int>                                     aMF;

  // Unify pave blocks of the existing edges united on the post-treat stage
  NCollection_DataMap<int, NCollection_List<occ::handle<BOPDS_PaveBlock>>> anEdgeLPB;

  NCollection_DynamicArray<BOPDS_InterfFF>& aFFs = myDS->InterfFF();
  aNbFF                                          = aFFs.Length();
  // 1. Sections (curves, points);
  for (i = 0; i < aNbFF; ++i)
  {
    BOPDS_InterfFF& aFF = aFFs(i);
    aFF.Indices(nF1, nF2);
    //
    BOPDS_FaceInfo& aFI1 = myDS->ChangeFaceInfo(nF1);
    BOPDS_FaceInfo& aFI2 = myDS->ChangeFaceInfo(nF2);
    //
    // 1.1. Section edges
    NCollection_DynamicArray<BOPDS_Curve>& aVNC = aFF.ChangeCurves();
    aNbC                                        = aVNC.Length();
    for (j = 0; j < aNbC; ++j)
    {
      BOPDS_Curve&                                    aNC   = aVNC(j);
      NCollection_List<occ::handle<BOPDS_PaveBlock>>& aLPBC = aNC.ChangePaveBlocks();
      //
      // Add section edges to face info
      aItLPB.Initialize(aLPBC);
      for (; aItLPB.More();)
      {
        const occ::handle<BOPDS_PaveBlock>& aPB = aItLPB.Value();
        //
        // Treat existing pave blocks
        if (theDME.IsBound(aPB))
        {
          NCollection_List<occ::handle<BOPDS_PaveBlock>>& aLPB = theDME.ChangeFind(aPB);
          UpdateExistingPaveBlocks(aPB, aLPB, thePBFacesMap);

          NCollection_List<occ::handle<BOPDS_PaveBlock>>::Iterator itLPB(aLPB);
          for (; itLPB.More(); itLPB.Next())
          {
            const int                                       nE      = itLPB.Value()->Edge();
            NCollection_List<occ::handle<BOPDS_PaveBlock>>* pLPBOnE = anEdgeLPB.ChangeSeek(nE);
            if (!pLPBOnE)
            {
              pLPBOnE = anEdgeLPB.Bound(nE, NCollection_List<occ::handle<BOPDS_PaveBlock>>());
            }
            pLPBOnE->Append(itLPB.Value());
          }

          aLPBC.Remove(aItLPB);
          continue;
        }
        //
        aFI1.ChangePaveBlocksSc().Add(aPB);
        aFI2.ChangePaveBlocksSc().Add(aPB);
        // Add edge-PB connection
        const int                                       nE      = aPB->Edge();
        NCollection_List<occ::handle<BOPDS_PaveBlock>>* pLPBOnE = anEdgeLPB.ChangeSeek(nE);
        if (!pLPBOnE)
        {
          pLPBOnE = anEdgeLPB.Bound(nE, NCollection_List<occ::handle<BOPDS_PaveBlock>>());
        }
        pLPBOnE->Append(aPB);

        aItLPB.Next();
      }
    }
    //
    // 1.2. Section vertices
    const NCollection_DynamicArray<BOPDS_Point>& aVNP = aFF.Points();
    aNbP                                              = aVNP.Length();
    for (j = 0; j < aNbP; ++j)
    {
      const BOPDS_Point& aNP = aVNP(j);
      nV1                    = aNP.Index();
      if (nV1 < 0)
      {
        continue;
      }
      aFI1.ChangeVerticesSc().Add(nV1);
      aFI2.ChangeVerticesSc().Add(nV1);
    }
    //
    aMF.Add(nF1);
    aMF.Add(nF2);
  }

  bool bNewCB = false;
  {
    // Unify pave blocks of the existing edges united on the post-treat stage
    // by making new Common blocks from them
    NCollection_DataMap<int, NCollection_List<occ::handle<BOPDS_PaveBlock>>>::Iterator itDM(
      anEdgeLPB);
    for (; itDM.More(); itDM.Next())
    {
      const NCollection_List<occ::handle<BOPDS_PaveBlock>>& aLPB = itDM.Value();
      if (aLPB.Extent() == 1)
      {
        continue;
      }

      bNewCB = true;

      // Find or create common block
      occ::handle<BOPDS_CommonBlock> aCB;
      // Collect all faces attached to common blocks
      NCollection_Map<int> aMFaces;
      // Collect all pave blocks attached to common blocks
      NCollection_IndexedMap<occ::handle<BOPDS_PaveBlock>> aMPaveBlocks;

      NCollection_List<occ::handle<BOPDS_PaveBlock>>::Iterator itLPB(aLPB);
      for (; itLPB.More(); itLPB.Next())
      {
        const occ::handle<BOPDS_PaveBlock>& aPB = itLPB.Value();
        aMPaveBlocks.Add(aPB);

        if (myDS->IsCommonBlock(aPB))
        {
          occ::handle<BOPDS_CommonBlock> aPBCB = myDS->CommonBlock(aPB);
          // Get pave blocks
          const NCollection_List<occ::handle<BOPDS_PaveBlock>>& aLPBOnCB = aPBCB->PaveBlocks();
          for (NCollection_List<occ::handle<BOPDS_PaveBlock>>::Iterator it(aLPBOnCB); it.More();
               it.Next())
          {
            aMPaveBlocks.Add(it.Value());
          }

          // Get faces
          const NCollection_List<int>& aLFacesOnCB = aPBCB->Faces();
          for (NCollection_List<int>::Iterator it(aLFacesOnCB); it.More(); it.Next())
          {
            aMFaces.Add(it.Value());
          }

          if (aCB.IsNull())
          {
            aCB = aPBCB;
          }
        }
      }

      if (aCB.IsNull())
      {
        // None of the pave blocks in the list is a common block,
        // so create the new one.
        aCB = new BOPDS_CommonBlock;
        aCB->SetPaveBlocks(aLPB);
        itLPB.Initialize(aLPB);
        for (; itLPB.More(); itLPB.Next())
        {
          const occ::handle<BOPDS_PaveBlock>& aPB = itLPB.Value();
          myDS->SetCommonBlock(aPB, aCB);
        }
      }
      else
      {
        // Update common block with new pave blocks
        NCollection_List<occ::handle<BOPDS_PaveBlock>> aLPBNew;
        {
          const int aNbPB = aMPaveBlocks.Extent();
          for (int iPB = 1; iPB <= aNbPB; ++iPB)
          {
            const occ::handle<BOPDS_PaveBlock>& aPB = aMPaveBlocks(iPB);
            myDS->SetCommonBlock(aPB, aCB);
            aLPBNew.Append(aPB);
          }
        }
        aCB->SetPaveBlocks(aLPBNew);

        // Update faces of the common block
        NCollection_List<int> aLFaces;
        for (NCollection_Map<int>::Iterator it(aMFaces); it.More(); it.Next())
        {
          aLFaces.Append(it.Value());
        }
        aCB->SetFaces(aLFaces);
      }
    }
  }

  bool bVerts = theDMV.Extent() > 0;
  bool bEdges = theDME.Extent() > 0 || bNewCB;
  //
  if (!bVerts && !bEdges)
  {
    return;
  }
  //
  // 2. Update Face Info information with new vertices and new
  //    pave blocks created in PostTreatFF from existing ones
  int                                     nV2;
  NCollection_Map<int>::Iterator          aItMF;
  NCollection_DataMap<int, int>::Iterator aItMV;
  //
  aItMF.Initialize(aMF);
  for (; aItMF.More(); aItMF.Next())
  {
    nF1 = aItMF.Value();
    //
    BOPDS_FaceInfo& aFI = myDS->ChangeFaceInfo(nF1);
    //
    // 2.1. Update information about vertices
    if (bVerts)
    {
      NCollection_Map<int>& aMVOn = aFI.ChangeVerticesOn();
      NCollection_Map<int>& aMVIn = aFI.ChangeVerticesIn();
      //
      aItMV.Initialize(theDMV);
      for (; aItMV.More(); aItMV.Next())
      {
        nV1 = aItMV.Key();
        nV2 = aItMV.Value();
        //
        if (aMVOn.Remove(nV1))
        {
          aMVOn.Add(nV2);
        }
        //
        if (aMVIn.Remove(nV1))
        {
          aMVIn.Add(nV2);
        }
      } // for (; aItMV.More(); aItMV.Next()) {
    } // if (bVerts) {
    //
    // 2.2. Update information about pave blocks
    if (bEdges)
    {
      NCollection_Map<occ::handle<BOPDS_PaveBlock>>         aMPBFence;
      NCollection_IndexedMap<occ::handle<BOPDS_PaveBlock>>* pMPB[] = {&aFI.ChangePaveBlocksOn(),
                                                                      &aFI.ChangePaveBlocksIn(),
                                                                      &aFI.ChangePaveBlocksSc()};
      for (i = 0; i < 3; ++i)
      {
        NCollection_IndexedMap<occ::handle<BOPDS_PaveBlock>> aMPBCopy = *pMPB[i];
        pMPB[i]->Clear();
        const int aNbPB = aMPBCopy.Extent();
        for (j = 1; j <= aNbPB; ++j)
        {
          const occ::handle<BOPDS_PaveBlock>&                   aPB  = aMPBCopy(j);
          const NCollection_List<occ::handle<BOPDS_PaveBlock>>* pLPB = theDME.Seek(aPB);
          if (pLPB && !pLPB->IsEmpty())
          {
            aItLPB.Initialize(*pLPB);
            for (; aItLPB.More(); aItLPB.Next())
            {
              const occ::handle<BOPDS_PaveBlock>& aPB1 = aItLPB.Value();
              const occ::handle<BOPDS_PaveBlock>& aPBR = myDS->RealPaveBlock(aPB1);
              if (aMPBFence.Add(aPBR))
              {
                pMPB[i]->Add(aPBR);
              }
            }
          }
          else
          {
            const occ::handle<BOPDS_PaveBlock>& aPBR = myDS->RealPaveBlock(aPB);
            if (aMPBFence.Add(aPBR))
            {
              pMPB[i]->Add(aPBR);
            }
          }
        } // for (j = 1; j <= aNbPB; ++j) {
      } // for (i = 0; i < 2; ++i) {
    } // if (bEdges) {
  }
}

//=================================================================================================

bool BOPAlgo_PaveFiller::IsExistingVertex(const gp_Pnt&               aP,
                                          const double                theTolR3D,
                                          const NCollection_Map<int>& aMVOnIn) const
{
  bool                           bRet;
  int                            nV, iFlag;
  double                         aTolCheck;
  Bnd_Box                        aBoxP;
  NCollection_Map<int>::Iterator aIt;
  //
  aTolCheck = theTolR3D + myFuzzyValue;
  bRet      = true;
  //
  aBoxP.Add(aP);
  aBoxP.Enlarge(theTolR3D);
  //
  aIt.Initialize(aMVOnIn);
  for (; aIt.More(); aIt.Next())
  {
    nV                           = aIt.Value();
    const BOPDS_ShapeInfo& aSIV  = myDS->ShapeInfo(nV);
    const TopoDS_Vertex&   aV    = (*(TopoDS_Vertex*)(&aSIV.Shape()));
    const Bnd_Box&         aBoxV = aSIV.Box();
    //
    if (!aBoxP.IsOut(aBoxV))
    {
      iFlag = BOPTools_AlgoTools::ComputeVV(aV, aP, aTolCheck);
      if (!iFlag)
      {
        return bRet;
      }
    }
  }
  return !bRet;
}

//=================================================================================================

bool BOPAlgo_PaveFiller::IsExistingPaveBlock(const occ::handle<BOPDS_PaveBlock>& thePB,
                                             const BOPDS_Curve&                  theNC,
                                             const NCollection_List<int>&        theLSE,
                                             int&                                theNEOut,
                                             double&                             theTolNew)
{
  if (theLSE.IsEmpty())
  {
    return false;
  }

  double                          aT1, aT2, aTm, aTx, aTolE, aTolCheck, aTol, aDist;
  int                             nE, iFlag, nV1, nV2;
  gp_Pnt                          aPm;
  Bnd_Box                         aBoxPm;
  NCollection_List<int>::Iterator aItLI;
  //
  thePB->Range(aT1, aT2);
  thePB->Indices(nV1, nV2);
  const TopoDS_Vertex &aV1 = TopoDS::Vertex(myDS->Shape(nV1)),
                      &aV2 = TopoDS::Vertex(myDS->Shape(nV2));
  const double aTolV1 = BRep_Tool::Tolerance(aV1), aTolV2 = BRep_Tool::Tolerance(aV2);

  aTol = std::max(aTolV1, aTolV2);

  aTm = IntTools_Tools::IntermediatePoint(aT1, aT2);
  theNC.Curve().D0(aTm, aPm);
  aBoxPm.Add(aPm);
  aBoxPm.Enlarge(aTol);
  //
  aItLI.Initialize(theLSE);
  for (; aItLI.More(); aItLI.Next())
  {
    nE = aItLI.Value();
    if (nE < 0)
    {
      continue;
    }
    const BOPDS_ShapeInfo& aSIE  = myDS->ChangeShapeInfo(nE);
    const Bnd_Box&         aBoxE = aSIE.Box();
    if (!aBoxE.IsOut(aBoxPm))
    {
      const TopoDS_Edge& aE = (*(TopoDS_Edge*)(&aSIE.Shape()));
      aTolE                 = BRep_Tool::Tolerance(aE);
      aTolCheck             = std::max(aTolE, aTol) + myFuzzyValue;
      iFlag                 = myContext->ComputePE(aPm, aTolCheck, aE, aTx, aDist);
      if (!iFlag)
      {
        theNEOut  = nE;
        theTolNew = aDist;
        return true;
      }
    }
  }
  return false;
}

//=================================================================================================

bool BOPAlgo_PaveFiller::IsExistingPaveBlock(
  const occ::handle<BOPDS_PaveBlock>&                         thePB,
  const BOPDS_Curve&                                          theNC,
  const double                                                theTolR3D,
  const NCollection_IndexedMap<occ::handle<BOPDS_PaveBlock>>& theMPBOnIn,
  BOPTools_BoxTree&                                           thePBTree,
  const NCollection_Map<occ::handle<BOPDS_PaveBlock>>&        theMPBCommon,
  occ::handle<BOPDS_PaveBlock>&                               aPBOut,
  double&                                                     theTolNew)
{
  const IntTools_Curve& aIC = theNC.Curve();

  double aT1, aT2;
  thePB->Range(aT1, aT2);

  int nV11, nV12;
  thePB->Indices(nV11, nV12);

  // first point
  Bnd_Box aBoxP1;
  gp_Pnt  aP1;
  aIC.D0(aT1, aP1);
  aBoxP1.Add(aP1);
  const double aTolV11 = BRep_Tool::Tolerance(TopoDS::Vertex(myDS->Shape(nV11)));
  aBoxP1.Enlarge(aTolV11);

  // Find edges intersecting by AABB with the first point
  BOPTools_BoxTreeSelector aSelector;
  aSelector.SetBox(Bnd_Tools::Bnd2BVH(aBoxP1));
  aSelector.SetBVHSet(&thePBTree);
  if (!aSelector.Select())
  {
    return false;
  }

  // intermediate point
  Bnd_Box                        aBoxPm;
  double                         aTm = IntTools_Tools::IntermediatePoint(aT1, aT2);
  gp_Pnt                         aPm;
  gp_Vec                         aVTgt1;
  const occ::handle<Geom_Curve>& aC3d = aIC.Curve();
  aC3d->D1(aTm, aPm, aVTgt1);
  aBoxPm.Add(aPm);
  bool isVtgt1Valid = aVTgt1.SquareMagnitude() > gp::Resolution();
  if (isVtgt1Valid)
  {
    aVTgt1.Normalize();
  }

  // last point
  Bnd_Box aBoxP2;
  gp_Pnt  aP2;
  aIC.D0(aT2, aP2);
  aBoxP2.Add(aP2);
  const double aTolV12 = BRep_Tool::Tolerance(TopoDS::Vertex(myDS->Shape(nV12)));
  aBoxP2.Enlarge(aTolV12);

  const double aTolV1 = std::max(aTolV11, aTolV12) + myFuzzyValue;

  double aTolCheck = theTolR3D + myFuzzyValue;

  // Some limit values to define "thin" face when iflag1=iflag2=2 and
  // edge has no common block with any face
  double       aMaxTolAdd   = 0.001; // Maximal tolerance of edge allowed
  const double aCoeffTolAdd = 10.;   // Coeff to define max. tolerance with help of aTolCheck
  aMaxTolAdd                = std::min(aMaxTolAdd, aCoeffTolAdd * aTolCheck);

  // Look for the existing pave block closest to the section curve
  bool bFound = false;
  theTolNew   = ::RealLast();

  for (NCollection_List<int>::Iterator it(aSelector.Indices()); it.More(); it.Next())
  {
    const occ::handle<BOPDS_PaveBlock>& aPB = theMPBOnIn(it.Value());

    int nV21, nV22;
    aPB->Indices(nV21, nV22);

    const double aTolV21 = BRep_Tool::Tolerance(TopoDS::Vertex(myDS->Shape(nV21)));
    const double aTolV22 = BRep_Tool::Tolerance(TopoDS::Vertex(myDS->Shape(nV22)));
    const double aTolV2  = std::max(aTolV21, aTolV22) + myFuzzyValue;

    const BOPDS_ShapeInfo& aSISp  = myDS->ChangeShapeInfo(aPB->Edge());
    const TopoDS_Edge&     aSp    = (*(TopoDS_Edge*)(&aSISp.Shape()));
    const Bnd_Box&         aBoxSp = aSISp.Box();

    int iFlag1 = (nV11 == nV21 || nV11 == nV22) ? 2 : 1;
    int iFlag2 = (nV12 == nV21 || nV12 == nV22) ? 2 : (!aBoxSp.IsOut(aBoxP2) ? 1 : 0);
    if (!iFlag2)
    {
      continue;
    }

    double aDist  = 0.;
    double aCoeff = 1.; // Coeff for taking in account deflections between edge and theNC
    // when aPB is not common block
    double aDistm1m2 = 0.;
    int    aPEStatus = 1;

    double aRealTol = aTolCheck;
    if (myDS->IsCommonBlock(aPB))
    {
      aRealTol = std::max(aRealTol, std::max(aTolV1, aTolV2));
      if (theMPBCommon.Contains(aPB))
      {
        // for an edge, which is a common block with a face,
        // increase the chance to coincide with section curve
        aRealTol *= 2.;
      }
    }
    else if (iFlag1 == 2 && iFlag2 == 2)
    {
      // Check, if edge could be common block with section curve
      //  and increase the chance to coincide with section curve
      // skip processing if one edge is closed, but other is not closed
      // such configurations can give iFlag1 == 2 && iFlag2 == 2
      bool bSkipProcessing =
        ((nV11 == nV12) && (nV21 != nV22)) || ((nV11 != nV12) && (nV21 == nV22));

      if (!bSkipProcessing)
      {

        if (isVtgt1Valid)
        {
          BRepAdaptor_Curve aBAC2(aSp);
          if (aIC.Type() != GeomAbs_Line || aBAC2.GetType() != GeomAbs_Line)
          {
            double aTldp;
            double aTolAdd =
              2. * std::min(aMaxTolAdd, std::max(aRealTol, std::max(aTolV1, aTolV2)));
            aPEStatus = myContext->ComputePE(aPm, aTolAdd, aSp, aTldp, aDistm1m2);

            if (aPEStatus == 0)
            {
              gp_Pnt aPm2;
              gp_Vec aVTgt2;
              aBAC2.D1(aTldp, aPm2, aVTgt2);
              if (aVTgt2.SquareMagnitude() > gp::Resolution())
              {
                // The angle should be close to zero
                double aCos = aVTgt1.Dot(aVTgt2.Normalized());
                if (std::abs(aCos) >= 0.9063)
                {
                  aRealTol = aTolAdd;
                  aCoeff   = 2.;
                }
              }
            }
          }
        }
      }
    }

    Bnd_Box aBoxTmp = aBoxPm;
    aBoxTmp.Enlarge(aRealTol);

    double aDistToSp = 0.;
    double aTx;
    if (aBoxSp.IsOut(aBoxTmp) || aPEStatus < 0)
    {
      continue;
    }
    else if (aPEStatus == 0) // aPEStatus == 0 for case iflag1 == iflag2 == 2
    {
      aDistToSp = aDistm1m2;
    }
    else if (aPEStatus == 1) // Projection has not been done yet
    {
      aPEStatus = myContext->ComputePE(aPm, aRealTol, aSp, aTx, aDistToSp);
      if (aPEStatus < 0)
      {
        continue;
      }
    }
    //
    if (iFlag1 == 1)
    {
      iFlag1 = !myContext->ComputePE(aP1, aRealTol, aSp, aTx, aDist);
      if (iFlag1 && aDistToSp < aDist)
      {
        aDistToSp = aDist;
      }
    }
    //
    if (iFlag2 == 1)
    {
      iFlag2 = !myContext->ComputePE(aP2, aRealTol, aSp, aTx, aDist);
      if (iFlag2 && aDistToSp < aDist)
      {
        aDistToSp = aDist;
      }
    }
    //
    if (iFlag1 && iFlag2)
    {
      if (aDistToSp < theTolNew)
      {
        aPBOut    = aPB;
        theTolNew = aCoeff * aDistToSp;
        bFound    = true;
      }
    }
  }
  return bFound;
}

//=================================================================================================

static void getBoundPaves(const BOPDS_DS* theDS, const BOPDS_Curve& theNC, int theNV[2])
{
  theNV[0] = theNV[1] = -1;

  // get extreme paves
  const occ::handle<BOPDS_PaveBlock>& aPB   = theNC.PaveBlocks().First();
  const NCollection_List<BOPDS_Pave>& aLP   = aPB->ExtPaves();
  int                                 aNbEP = aLP.Extent();
  if (aNbEP == 0)
  {
    return;
  }
  double aTmin = RealLast();
  double aTmax = -aTmin;
  for (NCollection_List<BOPDS_Pave>::Iterator aItLP(aLP); aItLP.More(); aItLP.Next())
  {
    const BOPDS_Pave& aPv = aItLP.Value();
    int               nV;
    double            aTV;
    aPv.Contents(nV, aTV);
    if (aTV < aTmin)
    {
      theNV[0] = aPv.Index();
      aTmin    = aTV;
    }
    if (aTV > aTmax)
    {
      theNV[1] = aPv.Index();
      aTmax    = aTV;
    }
  }

  // compare extreme vertices with ends of the curve
  const IntTools_Curve& aIC = theNC.Curve();
  double                aT[2];
  gp_Pnt                aP[2];
  aIC.Bounds(aT[0], aT[1], aP[0], aP[1]);
  double aTol = std::max(theNC.Tolerance(), theNC.TangentialTolerance());
  aTol += Precision::Confusion();
  for (int j = 0; j < 2; ++j)
  {
    const BOPDS_ShapeInfo& aSIV  = theDS->ShapeInfo(theNV[j]);
    const TopoDS_Vertex&   aV    = (*(TopoDS_Vertex*)(&aSIV.Shape()));
    int                    iFlag = BOPTools_AlgoTools::ComputeVV(aV, aP[j], aTol);
    if (iFlag != 0)
    {
      theNV[j] = -1;
    }
  }
}

//=================================================================================================

// The normal of the face's surface at the projection of the point; false if
// it does not project.
static bool normalAt(const occ::handle<IntTools_Context>& theCtx,
                     const TopoDS_Face&                   theF,
                     const gp_Pnt&                        theP,
                     gp_Dir&                              theN)
{
  GeomAPI_ProjectPointOnSurf& aPS = theCtx->ProjPS(theF);
  aPS.Perform(theP);
  if (!aPS.IsDone() || !aPS.NbPoints())
  {
    return false;
  }
  double u, v;
  aPS.LowerDistanceParameters(u, v);
  BOPTools_AlgoTools3D::GetNormalToSurface(BRep_Tool::Surface(theF), u, v, theN);
  return true;
}

//=================================================================================================

// A vertex of both faces for the end <theEnd> of the section curve, which
// lies outside one of them, at a tangent end. Two faces tangent at a corner
// -- the offsets of two fillets each meeting a third face tangentially
// there -- cut each other in a curve that passes the corner at a hair's
// distance, not through it, and leaves the faces there across the boundary
// of one: the end is not on both faces and takes no vertex, and with one
// end only the section has no edge. The corner's vertex is the end if the
// faces are tangent at it and stay within the curve's tolerance of each other
// all the way from the end to it (the distance times the sine of the angle
// between them at the end), and it is near the end, not as far as the curve
// is long. Returns the vertex and its distance, or -1.
static int tangentEndVertex(const BOPDS_PDS&                     theDS,
                            const occ::handle<IntTools_Context>& theCtx,
                            const TopoDS_Face&                   theF1,
                            const TopoDS_Face&                   theF2,
                            const BOPDS_Curve&                   theNC,
                            const int                            theEnd,
                            const NCollection_Map<int>&          theMV,
                            double&                              theDist)
{
  const IntTools_Curve& aIC = theNC.Curve();
  double                aT[2];
  gp_Pnt                aP[2];
  aIC.Bounds(aT[0], aT[1], aP[0], aP[1]);
  const double aTolR3D = std::max(theNC.Tolerance(), theNC.TangentialTolerance());
  const double aDMax   = 0.1 * aP[0].Distance(aP[1]);
  gp_Dir       aDN[2];
  if (!normalAt(theCtx, theF1, aP[theEnd], aDN[0]) || !normalAt(theCtx, theF2, aP[theEnd], aDN[1]))
  {
    return -1;
  }
  // the sine of the angle between the faces at the end
  const double aSin  = gp_Vec(aDN[0]).Crossed(gp_Vec(aDN[1])).Magnitude();
  int          nVEnd = -1;
  theDist            = RealLast();
  for (NCollection_Map<int>::Iterator aIt(theMV); aIt.More(); aIt.Next())
  {
    const int    nV  = aIt.Value();
    const gp_Pnt aPV = BRep_Tool::Pnt(TopoDS::Vertex(theDS->Shape(nV)));
    const double aD  = aPV.Distance(aP[theEnd]);
    if (aD >= theDist || aD > aDMax || aD * aSin > aTolR3D)
    {
      continue;
    }
    gp_Dir aDNV[2];
    if (!normalAt(theCtx, theF1, aPV, aDNV[0]) || !normalAt(theCtx, theF2, aPV, aDNV[1])
        || 1. - std::abs(aDNV[0] * aDNV[1]) > 5.e-9) // the creasing criteria
    {
      continue;
    }
    nVEnd   = nV;
    theDist = aD;
  }
  return nVEnd;
}

//=================================================================================================

// A vertex of both faces for the end <theEnd> of the section curve, which
// stops short of it. Where a section runs into a corner of both faces along
// their boundaries -- the rim of a cylinder lying in the plane the other
// face ends in, that face's bounding arc tangent to the rim at the corner
// -- the two faces and their two boundaries are all within tolerance of
// each other well before the corner: in realthunder/FreeCAD#280, a pocket's
// cylinder against a revolved end of the part, the arc is within 1e-7 of
// the cylinder for 0.25 before the corner. The walking line stops anywhere
// in there (0.0051 short). Its end is on both faces and would take a vertex
// of its own, and the pieces of the two boundaries from that vertex to the
// corner, a hair apart (2e-7), are not found the same: a face of no area
// between them, edges no two faces share, no solid split, and the cut came
// back the part unchanged, the common and the fuse empty, with no error.
// The corner's vertex is the end if the section goes on to it: it lies ahead
// of the end, the way the curve runs, the points from the end to it are on
// both surfaces within the curve's tolerance, and it is near, not as far as
// a tenth of the curve. Returns the vertex and its distance, or -1.
static int continuedEndVertex(const BOPDS_PDS&                     theDS,
                              const occ::handle<IntTools_Context>& theCtx,
                              const TopoDS_Face&                   theF1,
                              const TopoDS_Face&                   theF2,
                              const BOPDS_Curve&                   theNC,
                              const int                            theEnd,
                              const NCollection_Map<int>&          theMV,
                              double&                              theDist)
{
  const IntTools_Curve& aIC = theNC.Curve();
  double                aT[2];
  gp_Pnt                aP[2];
  aIC.Bounds(aT[0], aT[1], aP[0], aP[1]);
  const double aTolR3D = std::max(theNC.Tolerance(), theNC.TangentialTolerance());
  const double aDMax   = 0.1 * aP[0].Distance(aP[1]);
  // the way the curve runs out of its end
  gp_Pnt aPE;
  gp_Vec aVT;
  aIC.Curve()->D1(aT[theEnd], aPE, aVT);
  if (aVT.Magnitude() < gp::Resolution())
  {
    return -1;
  }
  if (theEnd == 0)
  {
    aVT.Reverse();
  }
  const gp_Dir aDT(aVT);
  int          nVEnd = -1;
  theDist            = RealLast();
  for (NCollection_Map<int>::Iterator aIt(theMV); aIt.More(); aIt.Next())
  {
    const int            nV  = aIt.Value();
    const TopoDS_Vertex& aV  = TopoDS::Vertex(theDS->Shape(nV));
    const gp_Pnt         aPV = BRep_Tool::Pnt(aV);
    const double         aD  = aPV.Distance(aP[theEnd]);
    // one within tolerance of the end is the end's already
    if (aD >= theDist || aD > aDMax || aD <= aTolR3D + BRep_Tool::Tolerance(aV))
    {
      continue;
    }
    const gp_Vec aVPV(aP[theEnd], aPV);
    if (aVPV * gp_Vec(aDT) < 0.9 * aD)
    {
      continue;
    }
    bool bOn = true;
    for (int i = 1; i <= 4 && bOn; ++i)
    {
      const gp_Pnt aPi = aP[theEnd].Translated(aVPV * (0.25 * i));
      for (const TopoDS_Face* pF : {&theF1, &theF2})
      {
        GeomAPI_ProjectPointOnSurf& aPS = theCtx->ProjPS(*pF);
        aPS.Perform(aPi);
        if (!aPS.IsDone() || !aPS.NbPoints() || aPS.LowerDistance() > aTolR3D)
        {
          bOn = false;
          break;
        }
      }
    }
    if (!bOn)
    {
      continue;
    }
    nVEnd   = nV;
    theDist = aD;
  }
  return nVEnd;
}

//=================================================================================================

// The part of an edge the tolerance of the vertex theNV is stretched along
// (BOPAlgo_PaveFiller::ExtendedTolerance(), which gives the distance only):
// the edge of the edge/edge or edge/face interference that made the vertex,
// and its range there.
static bool stretchedPart(const BOPDS_PDS&            theDS,
                          const int                   theNV,
                          const NCollection_Map<int>& theMI,
                          const int                   theType,
                          TopoDS_Edge&                theE,
                          double&                     theT1,
                          double&                     theT2)
{
  if (!theDS->IsNewShape(theNV))
  {
    return false;
  }
  NCollection_DynamicArray<BOPDS_InterfEE>& aEEs = theDS->InterfEE();
  NCollection_DynamicArray<BOPDS_InterfEF>& aEFs = theDS->InterfEF();
  for (int k = theType == 2 ? 1 : 0; k < (theType == 1 ? 1 : 2); ++k)
  {
    const int aNb = !k ? aEEs.Length() : aEFs.Length();
    for (int i = 0; i < aNb; ++i)
    {
      const BOPDS_Interf* aInt =
        !k ? (const BOPDS_Interf*)(&aEEs(i)) : (const BOPDS_Interf*)(&aEFs(i));
      if (aInt->IndexNew() == theNV && theMI.Contains(aInt->Index1())
          && theMI.Contains(aInt->Index2()))
      {
        const IntTools_CommonPrt& aCP = !k ? aEEs(i).CommonPart() : aEFs(i).CommonPart();
        theE                          = aCP.Edge1();
        aCP.Range1(theT1, theT2);
        return !theE.IsNull();
      }
    }
  }
  return false;
}

//=================================================================================================

// Whether the section curve comes within theTol of the part [theT1, theT2]
// of the edge anywhere within theR of the point thePV, which is nearest the
// curve at theT: the curve's points each way from there, a sixteenth of
// theR apart, and between the nearest of them and its neighbours.
static bool curveMeetsEdgePart(const occ::handle<IntTools_Context>& theCtx,
                               const IntTools_Curve&                theIC,
                               const double                         theT,
                               const gp_Pnt&                        thePV,
                               const double                         theR,
                               const TopoDS_Edge&                   theE,
                               const double                         theT1,
                               const double                         theT2,
                               const double                         theTol)
{
  const occ::handle<Geom_Curve>& aC = theIC.Curve();
  if (aC.IsNull() || theR <= gp::Resolution())
  {
    return true;
  }
  double aTC[2] = {aC->FirstParameter(), aC->LastParameter()};
  if (theIC.HasBounds())
  {
    gp_Pnt aPC[2];
    theIC.Bounds(aTC[0], aTC[1], aPC[0], aPC[1]);
  }
  GeomAPI_ProjectPointOnCurve& aProj  = theCtx->ProjPC(theE);
  const double                 aTolT  = 1.e-3 * (theT2 - theT1) + Precision::PConfusion();
  const auto                   toPart = [&](const double theTc) -> double {
    aProj.Perform(aC->Value(theTc));
    if (!aProj.NbPoints())
    {
      return RealLast();
    }
    const double aTE = aProj.LowerDistanceParameter();
    return aTE >= theT1 - aTolT && aTE <= theT2 + aTolT ? aProj.LowerDistance() : RealLast();
  };
  if (toPart(theT) <= theTol)
  {
    return true;
  }
  for (int k = 0; k < 2; ++k)
  {
    // along the curve this way, while within theR of the point
    double aDt    = (aTC[k] - theT) / 4096.;
    double aTc    = theT, aTPrev = theT, aTBest = theT, aTBestPrev = theT, aTBestNext = theT;
    double aDBest = RealLast();
    gp_Pnt aPc    = aC->Value(theT);
    for (int n = 0; n < 512 && std::abs(aDt) > 1.e-12 * std::abs(aTC[k] - theT); ++n)
    {
      double aTn = aTc + aDt;
      if ((aTC[k] - aTn) * aDt <= 0.)
      {
        aTn = aTC[k];
      }
      const gp_Pnt aPn = aC->Value(aTn);
      if (aPn.Distance(aPc) > theR / 16. && std::abs(aDt) > 1.e-9 * std::abs(aTC[k] - theT))
      {
        aDt *= 0.5;
        continue;
      }
      if (aPn.Distance(thePV) > theR)
      {
        break;
      }
      const double aD = toPart(aTn);
      if (aD <= theTol)
      {
        return true;
      }
      if (aD < aDBest)
      {
        aDBest     = aD;
        aTBest     = aTn;
        aTBestPrev = aTc;
        aTBestNext = aTn; // set when the next point is made
      }
      else if (aTBest == aTc)
      {
        aTBestNext = aTn;
      }
      if (aPn.Distance(aPc) < theR / 64.)
      {
        aDt *= 2.;
      }
      aTPrev = aTc;
      aTc    = aTn;
      aPc    = aPn;
      if (aTn == aTC[k])
      {
        break;
      }
    }
    (void)aTPrev;
    if (aDBest >= RealLast())
    {
      continue;
    }
    // nearer between the nearest point and its neighbours
    double aA = std::min(aTBestPrev, aTBestNext), aB = std::max(aTBestPrev, aTBestNext);
    for (int n = 0; n < 48 && aB - aA > 1.e-14 * std::max(1., std::abs(aB)); ++n)
    {
      const double aM1 = aA + (aB - aA) / 3., aM2 = aB - (aB - aA) / 3.;
      const double aD1 = toPart(aM1), aD2 = toPart(aM2);
      if (std::min(aD1, aD2) <= theTol)
      {
        return true;
      }
      if (aD1 < aD2)
      {
        aB = aM2;
      }
      else
      {
        aA = aM1;
      }
    }
  }
  return false;
}

//=================================================================================================

// Takes a vertex off a section curve that it was put on by its stretched
// tolerance and that does not come to it. A vertex made where an edge meets
// a face (or another edge) has its tolerance stretched, for the curves of
// that face, to the ends of the part of the edge that lies on the face round
// it (BOPAlgo_PaveFiller::ExtendedTolerance()), and a curve within that
// distance takes the vertex and grows it to reach
// (BOPAlgo_PaveFiller::PutPaveOnCurve()). FilterPavesOnCurves() takes such a
// vertex off a curve again where another curve has it nearer, unless the
// vertex is square to the curve -- the curve passing by it. Passing by is not
// coming to it: a prism standing 0.01 inside a torus on a profile that is the
// torus' own cut has its bottom edge cross the torus at a grazing angle,
// lying on it for 0.3 each way. One section of the wall with the torus ends
// on that vertex, as it should; another, a line 0.01 from the edge all along
// which turns away 0.046 from the vertex, took the vertex too and grew it to
// 0.046. Five edges then end in one vertex, a face of the torus is not made
// (BOPAlgo_AlertFaceBuilderUnusedEdges), and the common of the two solids is
// nothing, with no error. So: where another curve of the pair has the vertex
// within its own tolerance, a curve that has it further off than its
// tolerance, passes on both ways beyond that distance, and nowhere on the
// way comes to the part of the edge the tolerance was stretched along,
// gives the vertex up. The vertex's tolerance goes back to what the curves
// that keep it need.
static void dropStretchedPaves(const BOPDS_PDS&                             theDS,
                               const occ::handle<IntTools_Context>&         theCtx,
                               const NCollection_DynamicArray<BOPDS_Curve>& theVC,
                               const NCollection_Map<int>&                  theMI,
                               const double                                 theFuzzy,
                               NCollection_DataMap<int, double>&            theMVTol)
{
  for (int i = 0; i < theVC.Length(); ++i)
  {
    const BOPDS_Curve&    aNC = theVC(i);
    const IntTools_Curve& aIC = aNC.Curve();
    if (!aIC.HasBounds() || aNC.PaveBlocks().IsEmpty())
    {
      continue;
    }
    double aTB[2];
    gp_Pnt aPB[2];
    aIC.Bounds(aTB[0], aTB[1], aPB[0], aPB[1]);
    const double aTolR3D = std::max(aNC.Tolerance(), aNC.TangentialTolerance());
    const occ::handle<BOPDS_PaveBlock>& aPB1 = aNC.PaveBlocks().First();
    NCollection_List<int>               aLOff;
    for (NCollection_List<BOPDS_Pave>::Iterator aItP(aPB1->ExtPaves()); aItP.More(); aItP.Next())
    {
      const int     nV   = aItP.Value().Index();
      const double* pTol = theMVTol.Seek(nV); // what it was before a curve grew it
      if (!pTol || !theDS->IsNewShape(nV))
      {
        continue;
      }
      const TopoDS_Vertex& aV  = TopoDS::Vertex(theDS->Shape(nV));
      const gp_Pnt         aPV = BRep_Tool::Pnt(aV);
      const double         aD  = aPV.Distance(aIC.Curve()->Value(aItP.Value().Parameter()));
      if (aD <= aTolR3D + *pTol + theFuzzy || aPB[0].Distance(aPV) <= 2. * aD
          || aPB[1].Distance(aPV) <= 2. * aD)
      {
        continue;
      }
      // another curve of the pair has it
      bool bOther = false;
      for (int j = 0; j < theVC.Length() && !bOther; ++j)
      {
        const BOPDS_Curve& aNCj = theVC(j);
        if (j == i || aNCj.PaveBlocks().IsEmpty())
        {
          continue;
        }
        const double aTolJ = std::max(aNCj.Tolerance(), aNCj.TangentialTolerance());
        for (NCollection_List<BOPDS_Pave>::Iterator aItJ(aNCj.PaveBlocks().First()->ExtPaves());
             aItJ.More() && !bOther;
             aItJ.Next())
        {
          bOther = aItJ.Value().Index() == nV
                   && aPV.Distance(aNCj.Curve().Curve()->Value(aItJ.Value().Parameter()))
                        <= aTolJ + *pTol + theFuzzy;
        }
      }
      if (!bOther)
      {
        continue;
      }
      TopoDS_Edge aEPart;
      double      aTPart[2] = {0., 0.};
      if (stretchedPart(theDS, nV, theMI, 0, aEPart, aTPart[0], aTPart[1])
          && !curveMeetsEdgePart(theCtx,
                                 aIC,
                                 aItP.Value().Parameter(),
                                 aPV,
                                 2. * aD,
                                 aEPart,
                                 aTPart[0],
                                 aTPart[1],
                                 aTolR3D + std::max(*pTol, BRep_Tool::Tolerance(aEPart))
                                   + theFuzzy))
      {
        aLOff.Append(nV);
      }
    }
    for (NCollection_List<int>::Iterator aItV(aLOff); aItV.More(); aItV.Next())
    {
      const int nV = aItV.Value();
      aPB1->RemoveExtPave(nV);
      // as wide as the curves that keep it need
      const TopoDS_Vertex& aV   = TopoDS::Vertex(theDS->Shape(nV));
      const gp_Pnt         aPV  = BRep_Tool::Pnt(aV);
      double               aTol = theMVTol(nV);
      for (int j = 0; j < theVC.Length(); ++j)
      {
        if (theVC(j).PaveBlocks().IsEmpty())
        {
          continue;
        }
        for (NCollection_List<BOPDS_Pave>::Iterator aItJ(theVC(j).PaveBlocks().First()->ExtPaves());
             aItJ.More();
             aItJ.Next())
        {
          if (aItJ.Value().Index() == nV)
          {
            aTol = std::max(aTol,
                            aPV.Distance(theVC(j).Curve().Curve()->Value(aItJ.Value().Parameter()))
                              + Precision::Confusion());
          }
        }
      }
      if (aTol < BRep_Tool::Tolerance(aV))
      {
        (*(occ::handle<BRep_TVertex>*)&aV.TShape())->Tolerance(aTol);
      }
    }
  }
}

//=================================================================================================

// The parameters of theNb points evenly along the length of the curve from
// theT1 to theT2, the two among them. The length is the polygon's through the
// curve's knots and the points between them (512 pieces at least, four a
// span). A walking line's parameter is far from its length, and
// GCPnts_UniformAbscissa integrates each span's length down to the curve's
// resolution, which for one such line is 5e-24 (the pad on a torus cut
// through its plane of symmetry, away from it: more than 100 s a boolean in
// curveAlongEdges() alone).
static bool evenParameters(const occ::handle<Geom_Curve>&    theC,
                           const double                      theT1,
                           const double                      theT2,
                           const int                         theNb,
                           NCollection_DynamicArray<double>& thePars)
{
  thePars.Clear();
  if (theC.IsNull() || theNb < 2 || theT2 - theT1 <= Precision::PConfusion())
  {
    return false;
  }
  NCollection_DynamicArray<double> aTK;
  aTK.Append(theT1);
  const occ::handle<Geom_BSplineCurve> aBS = occ::down_cast<Geom_BSplineCurve>(theC);
  for (int i = 1; !aBS.IsNull() && i <= aBS->NbKnots(); ++i)
  {
    const double aK = aBS->Knot(i);
    if (aK > aTK.Last() && aK < theT2)
    {
      aTK.Append(aK);
    }
  }
  aTK.Append(theT2);
  const int                        aNbSub = std::max(4, 512 / (aTK.Length() - 1) + 1);
  NCollection_DynamicArray<double> aT, aL;
  gp_Pnt                           aPPrev = theC->Value(theT1);
  aT.Append(theT1);
  aL.Append(0.);
  for (int i = 1; i < aTK.Length(); ++i)
  {
    for (int k = 1; k <= aNbSub; ++k)
    {
      const double aTi = k == aNbSub ? aTK(i) : aTK(i - 1) + (aTK(i) - aTK(i - 1)) * k / aNbSub;
      const gp_Pnt aPi = theC->Value(aTi);
      aT.Append(aTi);
      aL.Append(aL.Last() + aPi.Distance(aPPrev));
      aPPrev = aPi;
    }
  }
  const double aLen = aL.Last();
  if (aLen <= gp::Resolution())
  {
    return false;
  }
  int i = 1;
  thePars.Append(theT1);
  for (int k = 1; k < theNb - 1; ++k)
  {
    const double aLk = aLen * k / (theNb - 1);
    while (i < aL.Length() - 1 && aL(i) < aLk)
    {
      ++i;
    }
    const double aDL = aL(i) - aL(i - 1);
    const double aDT = aT(i) - aT(i - 1);
    thePars.Append(aDL > gp::Resolution() ? aT(i - 1) + aDT * (aLk - aL(i - 1)) / aDL : aT(i));
  }
  thePars.Append(theT2);
  return true;
}

//=================================================================================================

// Whether the section curve runs along edges already on both faces: points
// inside it, evenly along its length (a walking line's parameter is not --
// #617's ran along an edge for 99% of its range, then climbed the faces'
// other branch for half its length), each on one of them.
static bool curveAlongEdges(
  const BOPDS_PDS&                                            theDS,
  const occ::handle<IntTools_Context>&                        theCtx,
  const BOPDS_Curve&                                          theNC,
  const NCollection_IndexedMap<occ::handle<BOPDS_PaveBlock>>& theMPBOnIn,
  const NCollection_Map<occ::handle<BOPDS_PaveBlock>>&        theMPBCommon)
{
  const IntTools_Curve& aIC = theNC.Curve();
  if (!aIC.HasBounds())
  {
    return false;
  }
  double aT[2];
  gp_Pnt aP[2];
  aIC.Bounds(aT[0], aT[1], aP[0], aP[1]);
  const double aTolR3D = std::max(theNC.Tolerance(), theNC.TangentialTolerance());
  NCollection_List<int> aLE;
  for (int i = 1; i <= theMPBOnIn.Extent(); ++i)
  {
    if (theMPBOnIn(i)->HasEdge())
    {
      aLE.Append(theMPBOnIn(i)->Edge());
    }
  }
  for (NCollection_Map<occ::handle<BOPDS_PaveBlock>>::Iterator aIt(theMPBCommon); aIt.More();
       aIt.Next())
  {
    if (aIt.Value()->HasEdge())
    {
      aLE.Append(aIt.Value()->Edge());
    }
  }
  NCollection_DynamicArray<double> aTE;
  if (!evenParameters(aIC.Curve(), aT[0], aT[1], 21, aTE))
  {
    return false;
  }
  for (int k = 1; k < aTE.Length() - 1; ++k)
  {
    const gp_Pnt aP  = aIC.Curve()->Value(aTE(k));
    bool         bOn = false;
    for (NCollection_List<int>::Iterator aItE(aLE); aItE.More() && !bOn; aItE.Next())
    {
      const TopoDS_Edge& aE = TopoDS::Edge(theDS->Shape(aItE.Value()));
      if (BRep_Tool::Degenerated(aE))
      {
        continue;
      }
      GeomAPI_ProjectPointOnCurve& aProj = theCtx->ProjPC(aE);
      aProj.Perform(aP);
      bOn = aProj.NbPoints() && aProj.LowerDistance() <= 10. * (aTolR3D + BRep_Tool::Tolerance(aE));
    }
    if (!bOn)
    {
      return false;
    }
  }
  return true;
}

// Whether the two faces stay within theTol of each other from thePA to thePB:
// points along the line between the two in the first face's parameters, each
// against the second face. The angle between the faces at one end times the
// distance bounds the same thing for faces that are flat between the two, and
// overstates it where they close on each other, as towards a point where
// they are tangent. The points start at thePA and end at thePB (what the two
// are off the first face is shared out along the line); thePnts takes them,
// theUV1 and theUV2 where they are on the two faces, theDist the farthest
// any is from either face.
static bool facesTogether(const occ::handle<IntTools_Context>& theCtx,
                          const TopoDS_Face&                   theF1,
                          const TopoDS_Face&                   theF2,
                          const gp_Pnt&                        thePA,
                          const gp_Pnt&                        thePB,
                          const double                         theTol,
                          NCollection_Array1<gp_Pnt>*          thePnts,
                          double*                              theDist,
                          NCollection_Array1<gp_Pnt2d>*        theUV1,
                          NCollection_Array1<gp_Pnt2d>*        theUV2)
{
  const occ::handle<Geom_Surface> aS1 = BRep_Tool::Surface(theF1);
  double                          aU[2], aV[2];
  gp_Vec                          aOff[2];
  GeomAPI_ProjectPointOnSurf&     aPS1 = theCtx->ProjPS(theF1);
  for (int i = 0; i < 2; ++i)
  {
    const gp_Pnt& aP = i == 0 ? thePA : thePB;
    aPS1.Perform(aP);
    if (!aPS1.IsDone() || !aPS1.NbPoints() || aPS1.LowerDistance() > theTol)
    {
      return false;
    }
    aPS1.LowerDistanceParameters(aU[i], aV[i]);
    aOff[i] = gp_Vec(aS1->Value(aU[i], aV[i]), aP);
  }
  if (aS1->IsUPeriodic())
  {
    aU[1] += aS1->UPeriod() * std::round((aU[0] - aU[1]) / aS1->UPeriod());
  }
  if (aS1->IsVPeriodic())
  {
    aV[1] += aS1->VPeriod() * std::round((aV[0] - aV[1]) / aS1->VPeriod());
  }
  GeomAPI_ProjectPointOnSurf& aPS2 = theCtx->ProjPS(theF2);
  const int                   aNb  = thePnts ? thePnts->Length() - 1 : 8;
  double                      aDM  = std::max(aOff[0].Magnitude(), aOff[1].Magnitude());
  for (int i = 0; i <= aNb; ++i)
  {
    const double aW = double(i) / aNb;
    const gp_Pnt aP =
      aS1->Value(aU[0] + (aU[1] - aU[0]) * aW, aV[0] + (aV[1] - aV[0]) * aW)
        .Translated(aOff[0] * (1. - aW) + aOff[1] * aW);
    aPS2.Perform(aP);
    if (!aPS2.IsDone() || !aPS2.NbPoints() || aPS2.LowerDistance() > theTol)
    {
      return false;
    }
    aDM = std::max(aDM, aPS2.LowerDistance());
    if (thePnts)
    {
      thePnts->ChangeValue(thePnts->Lower() + i) = aP;
    }
    if (theUV1)
    {
      theUV1->ChangeValue(theUV1->Lower() + i) =
        gp_Pnt2d(aU[0] + (aU[1] - aU[0]) * aW, aV[0] + (aV[1] - aV[0]) * aW);
    }
    if (theUV2)
    {
      double aU2, aV2;
      aPS2.LowerDistanceParameters(aU2, aV2);
      theUV2->ChangeValue(theUV2->Lower() + i) = gp_Pnt2d(aU2, aV2);
    }
  }
  if (theDist)
  {
    *theDist = aDM;
  }
  return true;
}

//=================================================================================================

// The edges already on both faces (curveAlongEdges()).
static void sharedEdges(const BOPDS_PDS&                                            theDS,
                        const NCollection_IndexedMap<occ::handle<BOPDS_PaveBlock>>& theMPBOnIn,
                        const NCollection_Map<occ::handle<BOPDS_PaveBlock>>&        theMPBCommon,
                        NCollection_List<int>&                                      theLE)
{
  for (int i = 1; i <= theMPBOnIn.Extent(); ++i)
  {
    if (theMPBOnIn(i)->HasEdge()
        && !BRep_Tool::Degenerated(TopoDS::Edge(theDS->Shape(theMPBOnIn(i)->Edge()))))
    {
      theLE.Append(theMPBOnIn(i)->Edge());
    }
  }
  for (NCollection_Map<occ::handle<BOPDS_PaveBlock>>::Iterator aIt(theMPBCommon); aIt.More();
       aIt.Next())
  {
    if (aIt.Value()->HasEdge()
        && !BRep_Tool::Degenerated(TopoDS::Edge(theDS->Shape(aIt.Value()->Edge()))))
    {
      theLE.Append(aIt.Value()->Edge());
    }
  }
}

//=================================================================================================

// Whether the point is on an edge both faces share as far as the faces tell:
// it is as near the edge as curveAlongEdges() asks, or the faces stay
// together (within the edge's tolerance, or theTol) from the point to the
// edge, which is no further than theDMax.
static bool onSharedEdge(const BOPDS_PDS&                     theDS,
                         const occ::handle<IntTools_Context>& theCtx,
                         const TopoDS_Face&                   theF1,
                         const TopoDS_Face&                   theF2,
                         const gp_Pnt&                        theP,
                         const NCollection_List<int>&         theLE,
                         const double                         theTol,
                         const double                         theDMax)
{
  for (NCollection_List<int>::Iterator aItE(theLE); aItE.More(); aItE.Next())
  {
    const TopoDS_Edge&           aE    = TopoDS::Edge(theDS->Shape(aItE.Value()));
    GeomAPI_ProjectPointOnCurve& aProj = theCtx->ProjPC(aE);
    aProj.Perform(theP);
    if (aProj.NbPoints() && aProj.LowerDistance() <= theDMax
        && (aProj.LowerDistance() <= 10. * (theTol + BRep_Tool::Tolerance(aE))
            || facesTogether(theCtx,
                             theF1,
                             theF2,
                             theP,
                             aProj.NearestPoint(),
                             std::max(theTol, BRep_Tool::Tolerance(aE)))))
    {
      return true;
    }
  }
  return false;
}

//=================================================================================================

// Whether the section curve, nothing put on it, is lost beside the edges both
// faces share: every point of it on one of them as far as the faces tell
// (onSharedEdge()), none further from it than a tenth of the curve. Round a
// point where the faces are tangent a walking line may be a stray piece
// beside the edge that runs through the point (#617's torus: 0.31 long,
// 0.002 off the pad's bottom edge, half a millimetre from the point); it took
// a vertex at each end and stayed in its face as an edge going nowhere.
static bool lostAlongEdges(const BOPDS_PDS&                     theDS,
                           const occ::handle<IntTools_Context>& theCtx,
                           const TopoDS_Face&                   theF1,
                           const TopoDS_Face&                   theF2,
                           const BOPDS_Curve&                   theNC,
                           const NCollection_List<int>&         theLE)
{
  const IntTools_Curve& aIC = theNC.Curve();
  if (!aIC.HasBounds() || theLE.IsEmpty() || !theNC.PaveBlocks().First()->ExtPaves().IsEmpty())
  {
    return false;
  }
  double aT[2];
  gp_Pnt aP[2];
  aIC.Bounds(aT[0], aT[1], aP[0], aP[1]);
  const double      aTolR3D = std::max(theNC.Tolerance(), theNC.TangentialTolerance());
  NCollection_DynamicArray<double> aTE;
  if (!evenParameters(aIC.Curve(), aT[0], aT[1], 9, aTE))
  {
    return false;
  }
  double aLen = 0.;
  for (int k = 1; k < aTE.Length(); ++k)
  {
    aLen += aIC.Curve()->Value(aTE(k)).Distance(aIC.Curve()->Value(aTE(k - 1)));
  }
  for (int k = 0; k < aTE.Length(); ++k)
  {
    if (!onSharedEdge(theDS,
                      theCtx,
                      theF1,
                      theF2,
                      aIC.Curve()->Value(aTE(k)),
                      theLE,
                      aTolR3D,
                      0.1 * aLen))
    {
      return false;
    }
  }
  return true;
}

//=================================================================================================

// Ends a section curve that runs along the edges both faces share
// (curveAlongEdges()) at the last vertex put on it, at each end that has
// none and is itself on one of those edges as far as the faces tell
// (onSharedEdge()). Such an end took a vertex of its own, and the piece from
// it to the next vertex became an edge beside the one already there (#617's
// pad: 6.79 long from 5e-4 off the tangent point, which no longer grows to
// take it in, see bridgeToVertex()). True if nothing of the curve is left.
static bool trimEndsOnEdges(const BOPDS_PDS&                     theDS,
                            const occ::handle<IntTools_Context>& theCtx,
                            const TopoDS_Face&                   theF1,
                            const TopoDS_Face&                   theF2,
                            BOPDS_Curve&                         theNC,
                            const NCollection_List<int>&         theLE)
{
  const IntTools_Curve& aIC = theNC.Curve();
  if (!aIC.HasBounds())
  {
    return false;
  }
  double aT[2];
  gp_Pnt aP[2];
  aIC.Bounds(aT[0], aT[1], aP[0], aP[1]);
  const NCollection_List<BOPDS_Pave>& aLP     = theNC.PaveBlocks().First()->ExtPaves();
  const double                        aTolR3D = std::max(theNC.Tolerance(), theNC.TangentialTolerance());
  const double                        aDMax   = 0.1 * aP[0].Distance(aP[1]);
  int                                 aEndNV[2] = {-1, -1};
  if (!aLP.IsEmpty())
  {
    getBoundPaves(theDS, theNC, aEndNV);
  }
  double aTN[2] = {aT[0], aT[1]};
  bool   bTrim  = false;
  for (int k = 0; k < 2; ++k)
  {
    if (aEndNV[k] >= 0
        || !onSharedEdge(theDS, theCtx, theF1, theF2, aP[k], theLE, aTolR3D, aDMax))
    {
      continue;
    }
    if (aLP.IsEmpty())
    {
      return true;
    }
    aTN[k] = k == 0 ? RealLast() : -RealLast();
    for (NCollection_List<BOPDS_Pave>::Iterator aIt(aLP); aIt.More(); aIt.Next())
    {
      aTN[k] = k == 0 ? std::min(aTN[k], aIt.Value().Parameter())
                      : std::max(aTN[k], aIt.Value().Parameter());
    }
    bTrim = true;
  }
  if (!bTrim)
  {
    return false;
  }
  if (aTN[1] - aTN[0] < Precision::PConfusion())
  {
    return true;
  }
  const occ::handle<Geom_BSplineCurve> aC3D = occ::down_cast<Geom_BSplineCurve>(aIC.Curve());
  occ::handle<Geom2d_BSplineCurve>     aC2D[2];
  aC2D[0] = occ::down_cast<Geom2d_BSplineCurve>(aIC.FirstCurve2d());
  aC2D[1] = occ::down_cast<Geom2d_BSplineCurve>(aIC.SecondCurve2d());
  if (aC3D.IsNull() || aC3D->IsPeriodic() || (aC2D[0].IsNull() && !aIC.FirstCurve2d().IsNull())
      || (aC2D[1].IsNull() && !aIC.SecondCurve2d().IsNull()))
  {
    return false;
  }
  occ::handle<Geom_BSplineCurve> aC3 = occ::down_cast<Geom_BSplineCurve>(aC3D->Copy());
  aC3->Segment(aTN[0], aTN[1]);
  occ::handle<Geom2d_BSplineCurve> aC2[2];
  for (int i = 0; i < 2; ++i)
  {
    if (!aC2D[i].IsNull())
    {
      aC2[i] = occ::down_cast<Geom2d_BSplineCurve>(aC2D[i]->Copy());
      aC2[i]->Segment(aTN[0], aTN[1]);
    }
  }
  theNC.SetCurve(IntTools_Curve(aC3, aC2[0], aC2[1], aIC.Tolerance(), aIC.TangentialTolerance()));
  return false;
}

//=================================================================================================

// A curve of its own from the section curve theVC(theJ) to the vertex theNV,
// which tangentVerticesNearCurve() found it passing near, the section curve
// ending where the new one starts. A walking line that comes to a point where
// the faces are tangent has them within tolerance of each other all round it
// and loses its way there. #617's pad on a torus, the tangent points on the
// pad's bottom edge, which lies on the torus: one line turns 0.07 short of the
// point, runs on beside the edge (0.005 to 0.03 off it, the nearest 0.04 from
// the point) and stops 0.17 past it on nothing; the other comes along the edge
// 1e-3 off it, crosses it 0.007 the wrong side of the point, hooks back over
// it 0.12 further on and leaves on the other branch. Putting the vertex on
// such a line where it is nearest grows the vertex to reach it (0.04 and
// 0.0075 there), keeps the hook as an edge and the run past the point as an
// edge going nowhere, and the faces' wires close only within those
// tolerances: areas and volumes read off them come out by that much times
// the size of the face. So the line is cut where it is still sound and the
// rest is made: the points on the first face between the cut and the vertex,
// every one within tolerance of the second (theTolV: the vertex's tolerance
// before any curve grew it; the edge it is on is known no better), a curve
// through them. The cut is the first point of the line, past its last sharp
// turn within 64 times the nearest, from which that holds and which is not on
// a shared edge itself (the curve from there to the vertex would run along
// that edge). One side of the line from the point taken must be lost, or be
// nothing, the line ending there. Lost is: nothing put on it, its end near,
// and its points on the vertex, on an edge both faces share as far as the
// faces tell (together from the point to there) or off a face; or, however
// far it goes and whatever is on it, every point of it on such an edge or
// off a face -- the line running on beside the edge it came to, which is
// that edge's already (the vertices on it go with it). The point taken is
// the line's nearest to the vertex, or another foot of the vertex on the line
// if that leaves less of it: a line may come by the vertex, run on, turn
// round and come by again (0.44 on, the pad on a torus cut through its plane
// of symmetry), and all of that is lost. And several lines may go by one
// vertex, the same line more than once, each lost its own way (three there,
// up to 2.4e-3 apart within 0.12 of the vertex, so that cut and joined each
// by itself they cross): the first joined to the vertex owns it (theOwners),
// and a later one is cut at the first of its points from which on it lies on
// the owner, and ends there with no curve of its own. So is one no side of
// which is lost, if one side lies on the owner: left as it is, it would be
// there twice between the vertex and where the owner was cut (the caller
// asks again for a line once its vertex has an owner). Returns false to have
// the vertex put on the line as it is; true with no curve if the line was
// cut for this vertex already (a pair of faces gone through again) or ends
// on the owner.
static bool bridgeToVertex(const BOPDS_PDS&                       theDS,
                           const occ::handle<IntTools_Context>&   theCtx,
                           const TopoDS_Face&                     theF1,
                           const TopoDS_Face&                     theF2,
                           NCollection_DynamicArray<BOPDS_Curve>& theVC,
                           const int                              theJ,
                           const int                              theNV,
                           const double                           theTolV,
                           const double                           theTolFF,
                           const NCollection_List<int>&           theLE,
                           occ::handle<Geom_Curve>&               theBridge,
                           occ::handle<Geom2d_Curve>              theBridge2d[2],
                           double&                                theTolB,
                           NCollection_DataMap<int, int>&         theOwners)
{
  BOPDS_Curve&          aNC = theVC.ChangeValue(theJ);
  const IntTools_Curve& aIC = aNC.Curve();
  const occ::handle<Geom_BSplineCurve> aC3D = occ::down_cast<Geom_BSplineCurve>(aIC.Curve());
  occ::handle<Geom2d_BSplineCurve>     aC2D[2];
  aC2D[0] = occ::down_cast<Geom2d_BSplineCurve>(aIC.FirstCurve2d());
  aC2D[1] = occ::down_cast<Geom2d_BSplineCurve>(aIC.SecondCurve2d());
  if (!aIC.HasBounds() || aC3D.IsNull() || aC3D->IsPeriodic()
      || (aC2D[0].IsNull() && !aIC.FirstCurve2d().IsNull())
      || (aC2D[1].IsNull() && !aIC.SecondCurve2d().IsNull()))
  {
    return false;
  }
  double aT[2];
  gp_Pnt aP[2];
  aIC.Bounds(aT[0], aT[1], aP[0], aP[1]);
  const gp_Pnt aPV     = BRep_Tool::Pnt(TopoDS::Vertex(theDS->Shape(theNV)));
  const double aTolR3D = std::max(aNC.Tolerance(), aNC.TangentialTolerance());
  const double aTol    = std::max(aTolR3D, theTolV);
  const double aDMax   = 0.1 * aP[0].Distance(aP[1]);
  // the curve's point nearest the vertex
  double aD = aPV.Distance(aP[0]), aTV = aT[0];
  if (aPV.Distance(aP[1]) < aD)
  {
    aD  = aPV.Distance(aP[1]);
    aTV = aT[1];
  }
  {
    GeomAPI_ProjectPointOnCurve aProj(aPV, aC3D, aT[0], aT[1]);
    if (aProj.NbPoints() && aProj.LowerDistance() < aD)
    {
      aD  = aProj.LowerDistance();
      aTV = aProj.LowerDistanceParameter();
    }
  }
  gp_Pnt aPC = aC3D->Value(aTV);
  // done already?
  for (int i = 0; i < theVC.Length(); ++i)
  {
    const IntTools_Curve& aICi = theVC(i).Curve();
    double                aTi[2];
    gp_Pnt                aPi[2];
    if (i != theJ && aICi.HasBounds() && aICi.Bounds(aTi[0], aTi[1], aPi[0], aPi[1])
        && ((aPi[0].Distance(aPC) <= Precision::Confusion()
             && aPi[1].Distance(aPV) <= Precision::Confusion())
            || (aPi[1].Distance(aPC) <= Precision::Confusion()
                && aPi[0].Distance(aPV) <= Precision::Confusion())))
    {
      if (!theOwners.IsBound(theNV))
      {
        theOwners.Bind(theNV, theJ);
      }
      return true;
    }
  }
  // the side that is lost, or nothing: its end near the vertex and nothing put
  // on it, or all of it beside the shared edges or off a face
  const NCollection_List<BOPDS_Pave>& aLP    = aNC.PaveBlocks().First()->ExtPaves();
  const auto                          isLost = [&](const double theTV, const int k) -> bool {
    const double aDT = aT[k] - theTV;
    if (std::abs(aDT) < Precision::PConfusion())
    {
      return true;
    }
    const bool bNear = aP[k].Distance(aPV) <= aDMax;
    bool       bLost = true;
    for (NCollection_List<BOPDS_Pave>::Iterator aIt(aLP); aIt.More() && bLost && bNear; aIt.Next())
    {
      bLost = (aIt.Value().Parameter() - theTV) * aDT <= 0.;
    }
    const int aNbS = bNear ? 8 : 64;
    for (int i = 1; i <= aNbS && bLost; ++i)
    {
      const gp_Pnt aPi = aC3D->Value(theTV + aDT * i / aNbS);
      if (bNear && facesTogether(theCtx, theF1, theF2, aPi, aPV, aTol))
      {
        continue;
      }
      bLost = onSharedEdge(theDS, theCtx, theF1, theF2, aPi, theLE, aTolR3D, aDMax)
              || !theCtx->IsValidPointForFaces(aPi, theF1, theF2, aTolR3D);
    }
    return bLost;
  };
  int kTail = -1;
  {
    const bool bL[2] = {isLost(aTV, 0), isLost(aTV, 1)};
    if (bL[0] && bL[1])
    {
      return false; // both sides: nothing of the curve to keep
    }
    if (bL[0] != bL[1])
    {
      kTail = bL[0] ? 0 : 1;
    }
  }
  // the line may come by the vertex once more, after a hook: what it does
  // before the last time is lost with the side
  {
    const double                aTV0  = aTV, aD0 = aD;
    double                      aKept = kTail < 0 ? RealLast() : std::abs(aT[1 - kTail] - aTV);
    GeomAPI_ProjectPointOnCurve aProj(aPV, aC3D, aT[0], aT[1]);
    for (int i = 1; i <= aProj.NbPoints(); ++i)
    {
      const double aTi = aProj.Parameter(i), aDi = aProj.Distance(i);
      if (std::abs(aTi - aTV0) < Precision::PConfusion() || aDi > aDMax)
      {
        continue;
      }
      const gp_Pnt aPi = aC3D->Value(aTi);
      if (!facesTogether(theCtx, theF1, theF2, aPi, aPV, aTol))
      {
        continue;
      }
      const bool bL[2] = {isLost(aTi, 0), isLost(aTi, 1)};
      if (bL[0] == bL[1])
      {
        continue;
      }
      const int k = bL[0] ? 0 : 1;
      if (std::abs(aT[1 - k] - aTi) < aKept)
      {
        aKept = std::abs(aT[1 - k] - aTi);
        kTail = k;
        aTV   = aTi;
        aD    = std::max(aDi, aD0);
        aPC   = aPi;
      }
    }
  }
  // the points of a side, out to 64 times the nearest, none past a vertex put
  // on the curve
  const double aDFar  = std::min(64. * aD, aDMax);
  const auto   sample = [&](const int                         theKB,
                            NCollection_DynamicArray<double>& theTS,
                            NCollection_DynamicArray<gp_Pnt>& thePS) {
    double       aTEnd = aT[theKB];
    const double aSign = aTEnd > aTV ? 1. : -1.;
    for (NCollection_List<BOPDS_Pave>::Iterator aIt(aLP); aIt.More(); aIt.Next())
    {
      const double aTP = aIt.Value().Parameter();
      if ((aTP - aTV) * aSign > 0. && (aTEnd - aTP) * aSign > 0.)
      {
        aTEnd = aTP;
      }
    }
    theTS.Clear();
    thePS.Clear();
    theTS.Append(aTV);
    thePS.Append(aPC);
    double aTc = aTV, aDt = (aTEnd - aTV) / 4096.;
    while (theTS.Length() < 1024 && std::abs(aDt) > 1.e-12 * std::abs(aTEnd - aTV))
    {
      const double aTn = aTc + aDt;
      if ((aTEnd - aTn) * aSign <= 0.)
      {
        break;
      }
      const gp_Pnt aPn   = aC3D->Value(aTn);
      const double aStep = aPn.Distance(thePS.Last());
      if (aStep > 0.5 * aD)
      {
        aDt *= 0.5;
        continue;
      }
      if (aPn.Distance(aPV) > aDFar)
      {
        break;
      }
      aTc = aTn;
      theTS.Append(aTc);
      thePS.Append(aPn);
      if (aStep < 0.125 * aD)
      {
        aDt *= 2.;
      }
    }
  };
  // another curve reaches this vertex already: this one is the same line once
  // more if a side of it lies on that one, and where the two have lost their
  // way they are not the same. The first of the side's points from which on
  // it is on the other, or -1
  const int* pOwner  = theOwners.Seek(theNV);
  const auto onOwner = [&](const NCollection_DynamicArray<gp_Pnt>& thePS) -> int {
    if (!pOwner || *pOwner == theJ)
    {
      return -1;
    }
    const BOPDS_Curve&    aNO = theVC(*pOwner);
    const IntTools_Curve& aIO = aNO.Curve();
    double                aTO[2];
    gp_Pnt                aPO[2];
    if (!aIO.HasBounds() || !aIO.Bounds(aTO[0], aTO[1], aPO[0], aPO[1]))
    {
      return -1;
    }
    const double aTolO = aTolR3D + std::max(aNO.Tolerance(), aNO.TangentialTolerance());
    GeomAPI_ProjectPointOnCurve aProjO;
    aProjO.Init(aIO.Curve(), aTO[0], aTO[1]);
    int iOn = thePS.Length();
    for (int i = thePS.Length() - 1; i >= 0; --i)
    {
      aProjO.Perform(thePS(i));
      if (!aProjO.NbPoints() || aProjO.LowerDistance() > aTolO)
      {
        break;
      }
      iOn = i;
    }
    return iOn < thePS.Length() - 1 ? iOn : -1;
  };
  NCollection_DynamicArray<double> aTS;
  NCollection_DynamicArray<gp_Pnt> aPS;
  if (kTail < 0)
  {
    // no side of it is lost. If one lies on the owner the other goes all the
    // same: left as it is, the line would be there twice between the vertex
    // and where the owner was cut
    for (int k = 0; k < 2; ++k)
    {
      sample(k, aTS, aPS);
      if (onOwner(aPS) >= 0)
      {
        if (kTail >= 0)
        {
          return false;
        }
        kTail = 1 - k;
      }
    }
    if (kTail < 0)
    {
      return false;
    }
  }
  const int kB = 1 - kTail;
  sample(kB, aTS, aPS);
  // sound from its last sharp turn on (20 degrees in a step, the steps under
  // half the nearest distance): the first point after it the faces are
  // together from
  int iFirst = 0;
  for (int i = 1; i < aTS.Length() - 1; ++i)
  {
    const gp_Vec aV1(aPS(i - 1), aPS(i)), aV2(aPS(i), aPS(i + 1));
    if (aV1.Magnitude() > gp::Resolution() && aV2.Magnitude() > gp::Resolution()
        && aV1 * aV2 < 0.94 * aV1.Magnitude() * aV2.Magnitude())
    {
      iFirst = i + 1;
    }
  }
  NCollection_Array1<gp_Pnt>   aPnts(1, 17);
  NCollection_Array1<gp_Pnt2d> aUV1(1, 17), aUV2(1, 17);
  double                       aTB = aTV;
  bool                         bOK = false;
  // cut where it is on the owner from there on, to end on it
  const int  iOn      = onOwner(aPS);
  const bool bOnOwner = iOn >= 0;
  if (bOnOwner)
  {
    aTB = aTS(iOn);
  }
  for (int i = iFirst; i < aTS.Length() && !bOK && !bOnOwner; ++i)
  {
    // not from a point on a shared edge: the curve from there to the vertex
    // would run along that edge
    if (onSharedEdge(theDS,
                     theCtx,
                     theF1,
                     theF2,
                     aPS(i),
                     theLE,
                     aTolR3D,
                     10. * (aTolR3D + Precision::Confusion())))
    {
      continue;
    }
    aTB = aTS(i);
    bOK =
      facesTogether(theCtx, theF1, theF2, aPS(i), aPV, aTol, &aPnts, &theTolB, &aUV1, &aUV2);
  }
  if (!bOnOwner && (!bOK || aPnts.First().Distance(aPnts.Last()) <= aTol))
  {
    return false;
  }
  occ::handle<Geom_Curve>   aBr3;
  occ::handle<Geom2d_Curve> aBr2[2];
  if (!bOnOwner)
  {
    // the curve through them
    occ::handle<NCollection_HArray1<gp_Pnt>> aHP = new NCollection_HArray1<gp_Pnt>(1, 17);
    occ::handle<NCollection_HArray1<double>> aHT = new NCollection_HArray1<double>(1, 17);
    const double                             aL  = aPnts.First().Distance(aPnts.Last());
    for (int i = 1; i <= 17; ++i)
    {
      aHP->SetValue(i, aPnts(i));
      aHT->SetValue(i, aL * (i - 1) / 16.);
    }
    GeomAPI_Interpolate anInt(aHP, aHT, false, Precision::Confusion());
    anInt.Perform();
    if (!anInt.IsDone())
    {
      return false;
    }
    // and where it is on the faces: the line on the first, through the points'
    // places on the second, kept on one side of a seam
    NCollection_Array1<gp_Pnt2d> aPoles1(1, 2);
    NCollection_Array1<double>   aKnots1(1, 2);
    NCollection_Array1<int>      aMults1(1, 2);
    aPoles1(1) = aUV1.First();
    aPoles1(2) = aUV1.Last();
    aKnots1(1) = 0.;
    aKnots1(2) = aL;
    aMults1.Init(2);
    const occ::handle<Geom_Surface>            aS2  = BRep_Tool::Surface(theF2);
    occ::handle<NCollection_HArray1<gp_Pnt2d>> aHUV = new NCollection_HArray1<gp_Pnt2d>(1, 17);
    for (int i = 1; i <= 17; ++i)
    {
      gp_Pnt2d aUV = aUV2(i);
      if (i > 1 && aS2->IsUPeriodic())
      {
        aUV.SetX(aUV.X()
                 + aS2->UPeriod()
                     * std::round((aHUV->Value(i - 1).X() - aUV.X()) / aS2->UPeriod()));
      }
      if (i > 1 && aS2->IsVPeriodic())
      {
        aUV.SetY(aUV.Y()
                 + aS2->VPeriod()
                     * std::round((aHUV->Value(i - 1).Y() - aUV.Y()) / aS2->VPeriod()));
      }
      aHUV->SetValue(i, aUV);
    }
    Geom2dAPI_Interpolate anInt2(aHUV, aHT, false, 1.e-13);
    anInt2.Perform();
    if (!anInt2.IsDone())
    {
      return false;
    }
    aBr3    = anInt.Curve();
    aBr2[0] = new Geom2d_BSplineCurve(aPoles1, aKnots1, aMults1, 1);
    aBr2[1] = anInt2.Curve();
  }
  // the section curve from the cut on
  aT[kTail] = aTB;
  if (aT[1] - aT[0] < Precision::PConfusion())
  {
    return false;
  }
  {
    NCollection_List<int> aLOff;
    for (NCollection_List<BOPDS_Pave>::Iterator aIt(aLP); aIt.More(); aIt.Next())
    {
      if (aIt.Value().Parameter() < aT[0] - Precision::PConfusion()
          || aIt.Value().Parameter() > aT[1] + Precision::PConfusion())
      {
        aLOff.Append(aIt.Value().Index());
      }
    }
    for (NCollection_List<int>::Iterator aIt(aLOff); aIt.More(); aIt.Next())
    {
      aNC.ChangePaveBlock1()->RemoveExtPave(aIt.Value());
    }
  }
  occ::handle<Geom_BSplineCurve> aC3 = occ::down_cast<Geom_BSplineCurve>(aC3D->Copy());
  aC3->Segment(aT[0], aT[1]);
  occ::handle<Geom2d_BSplineCurve> aC2[2];
  for (int i = 0; i < 2; ++i)
  {
    if (!aC2D[i].IsNull())
    {
      aC2[i] = occ::down_cast<Geom2d_BSplineCurve>(aC2D[i]->Copy());
      aC2[i]->Segment(aT[0], aT[1]);
    }
  }
  // its tolerance may have been the lost part's: what is left, measured
  double aTolC = aIC.Tolerance();
  {
    NCollection_DynamicArray<double> aTE;
    double aDev = evenParameters(aC3, aT[0], aT[1], 65, aTE) ? 0. : aTolC;
    for (int k = 0; aDev < aTolC && k < aTE.Length(); ++k)
    {
      const gp_Pnt aPk = aC3->Value(aTE(k));
      for (const TopoDS_Face* pF : {&theF1, &theF2})
      {
        GeomAPI_ProjectPointOnSurf& aPS = theCtx->ProjPS(*pF);
        aPS.Perform(aPk);
        aDev = aPS.IsDone() && aPS.NbPoints() ? std::max(aDev, aPS.LowerDistance()) : aTolC;
      }
    }
    aTolC = std::min(aTolC, std::max(2. * aDev, theTolFF));
  }
  if (aTolC < aIC.Tolerance() && aNC.Tolerance() <= aIC.Tolerance())
  {
    aNC.SetTolerance(aTolC);
  }
  aNC.SetCurve(IntTools_Curve(aC3, aC2[0], aC2[1], aTolC, aIC.TangentialTolerance()));
  if (bOnOwner)
  {
    return true;
  }
  theBridge      = aBr3;
  theBridge2d[0] = aBr2[0];
  theBridge2d[1] = aBr2[1];
  theTolB        = std::max(theTolB, theTolFF);
  if (!theOwners.IsBound(theNV))
  {
    theOwners.Bind(theNV, theJ);
  }
  return true;
}

//=================================================================================================

// The vertices of theMV the section curve passes near where the faces are
// tangent at them: each with the curve's point nearest to it, if the faces
// stay within the curve's tolerance of each other from there to the vertex
// (the distance times the sine of the angle between them at the point), as
// tangentEndVertex() takes a tangent end. A walking line through a point where
// two branches of the section cross -- a pad's wall on a revolution, tangent
// to it where the revolution turns vertical on the pad's bottom edge, which
// lies on it (#617) -- loses its way a few 1e-3 off the point and runs on the
// other branch or stops; the edge/face pass has put a vertex at the point.
// Where the angle says too far, the faces are measured (facesTogether()),
// against the vertex's own tolerance too (theMVTol: what it was before a
// curve grew it): the edge it is on is known no better than that. #617's pad
// on a torus, its profile a copy of the torus' section 6e-5 off it, has the
// line 0.04 from the vertex with the faces 2.9e-5 apart on the way, where the
// angle times the distance says 5.9e-5 and the line's own tolerance is 1.7e-5.
static void tangentVerticesNearCurve(const BOPDS_PDS&                        theDS,
                                     const occ::handle<IntTools_Context>&    theCtx,
                                     const TopoDS_Face&                      theF1,
                                     const TopoDS_Face&                      theF2,
                                     const BOPDS_Curve&                      theNC,
                                     const NCollection_Map<int>&             theMV,
                                     const NCollection_DataMap<int, double>& theMVTol,
                                     NCollection_List<int>&                  theLV,
                                     NCollection_List<double>&               theLD)
{
  const IntTools_Curve& aIC = theNC.Curve();
  if (!aIC.HasBounds())
  {
    return;
  }
  double aT[2];
  gp_Pnt aP[2];
  aIC.Bounds(aT[0], aT[1], aP[0], aP[1]);
  const double                   aTolR3D = std::max(theNC.Tolerance(), theNC.TangentialTolerance());
  const double                   aDMax   = 0.1 * aP[0].Distance(aP[1]);
  const occ::handle<Geom_Curve>& aC      = aIC.Curve();
  for (NCollection_Map<int>::Iterator aIt(theMV); aIt.More(); aIt.Next())
  {
    const int            nV  = aIt.Value();
    const TopoDS_Vertex& aV  = TopoDS::Vertex(theDS->Shape(nV));
    const gp_Pnt         aPV = BRep_Tool::Pnt(aV);
    GeomAPI_ProjectPointOnCurve aProj(aPV, aC, aT[0], aT[1]);
    double                      aD = aPV.Distance(aP[0]), aTV = aT[0];
    if (aPV.Distance(aP[1]) < aD)
    {
      aD  = aPV.Distance(aP[1]);
      aTV = aT[1];
    }
    if (aProj.NbPoints() && aProj.LowerDistance() < aD)
    {
      aD  = aProj.LowerDistance();
      aTV = aProj.LowerDistanceParameter();
    }
    if (aD <= BRep_Tool::Tolerance(aV) + aTolR3D || aD > aDMax)
    {
      continue; // on the curve already, or too far
    }
    const gp_Pnt aPC = aC->Value(aTV);
    gp_Dir       aDN[2], aDNV[2];
    if (!normalAt(theCtx, theF1, aPC, aDN[0]) || !normalAt(theCtx, theF2, aPC, aDN[1]))
    {
      continue;
    }
    if (gp_Vec(aDN[0]).Crossed(gp_Vec(aDN[1])).Magnitude() * aD > aTolR3D)
    {
      const double aTolV = theMVTol.IsBound(nV) ? theMVTol(nV) : BRep_Tool::Tolerance(aV);
      if (!facesTogether(theCtx, theF1, theF2, aPC, aPV, std::max(aTolR3D, aTolV)))
      {
        continue;
      }
    }
    if (!normalAt(theCtx, theF1, aPV, aDNV[0]) || !normalAt(theCtx, theF2, aPV, aDNV[1])
        || 1. - std::abs(aDNV[0] * aDNV[1]) > 5.e-9) // the creasing criteria
    {
      continue;
    }
    theLV.Append(nV);
    theLD.Append(aD);
  }
}

//=================================================================================================

void BOPAlgo_PaveFiller::PutBoundPaveOnCurve(const TopoDS_Face&     aF1,
                                             const TopoDS_Face&     aF2,
                                             BOPDS_Curve&           aNC,
                                             NCollection_List<int>& aLVB)
{
  const IntTools_Curve& aIC = aNC.Curve();
  double                aT[2];
  gp_Pnt                aP[2];
  aIC.Bounds(aT[0], aT[1], aP[0], aP[1]);
  double                        aTolR3D = std::max(aNC.Tolerance(), aNC.TangentialTolerance());
  occ::handle<BOPDS_PaveBlock>& aPB     = aNC.ChangePaveBlock1();
  // Get numbers of vertices assigned to the ends of the curve
  int aBndNV[2];
  getBoundPaves(myDS, aNC, aBndNV);
  //
  double aTolVnew = Precision::Confusion();
  bool   isClosed = aP[1].IsEqual(aP[0], aTolVnew);
  if (isClosed && (aBndNV[0] > 0 || aBndNV[1] > 0))
  {
    return;
  }

  for (int j = 0; j < 2; ++j)
  {
    if (aBndNV[j] < 0)
    {
      // no vertex on this end
      if (j && isClosed)
      {
        // if curve is closed, process only one bound
        continue;
      }
      bool bVF = myContext->IsValidPointForFaces(aP[j], aF1, aF2, aTolR3D);
      if (!bVF)
      {
        continue;
      }
      TopoDS_Vertex aVn;
      BOPTools_AlgoTools::MakeNewVertex(aP[j], aTolR3D, aVn);
      BOPTools_AlgoTools::UpdateVertex(aIC, aT[j], aVn);
      aTolVnew = BRep_Tool::Tolerance(aVn);

      BOPDS_ShapeInfo aSIVn;
      aSIVn.SetShapeType(TopAbs_VERTEX);
      aSIVn.SetShape(aVn);

      Bnd_Box& aBox = aSIVn.ChangeBox();
      BRepBndLib::Add(aVn, aBox);
      aBox.SetGap(aBox.GetGap() + Precision::Confusion());

      int nVn = myDS->Append(aSIVn);

      BOPDS_Pave aPn;
      aPn.SetIndex(nVn);
      aPn.SetParameter(aT[j]);
      aPB->AppendExtPave(aPn);

      aLVB.Append(nVn);
    }
  }
}

//=================================================================================================

void BOPAlgo_PaveFiller::PutPavesOnCurve(const NCollection_Map<int>&       theMVOnIn,
                                         const NCollection_Map<int>&       theMVCommon,
                                         BOPDS_Curve&                      theNC,
                                         const NCollection_Map<int>&       theMI,
                                         const NCollection_Map<int>&       theMVEF,
                                         NCollection_DataMap<int, double>& theMVTol,
                                         NCollection_DataMap<int, NCollection_List<int>>& theDMVLV)
{
  int                            nV;
  NCollection_Map<int>::Iterator aIt;
  //
  const Bnd_Box& aBoxC   = theNC.Box();
  const double   aTolR3D = std::max(theNC.Tolerance(), theNC.TangentialTolerance());
  //
  // Put EF vertices first
  aIt.Initialize(theMVEF);
  for (; aIt.More(); aIt.Next())
  {
    nV = aIt.Value();
    PutPaveOnCurve(nV, aTolR3D, theNC, theMI, theMVTol, theDMVLV, 2);
  }

  // Put all other vertices
  aIt.Initialize(theMVOnIn);
  for (; aIt.More(); aIt.Next())
  {
    nV = aIt.Value();
    if (theMVEF.Contains(nV))
    {
      continue;
    }

    if (!theMVCommon.Contains(nV))
    {
      const BOPDS_ShapeInfo& aSIV  = myDS->ShapeInfo(nV);
      const Bnd_Box&         aBoxV = aSIV.Box();
      //
      if (aBoxC.IsOut(aBoxV))
      {
        continue;
      }
      if (!myDS->IsNewShape(nV))
      {
        continue;
      }
    }
    //
    PutPaveOnCurve(nV, aTolR3D, theNC, theMI, theMVTol, theDMVLV, 1);
  }
}

//=================================================================================================

namespace
{
struct PaveBlockDist
{
  occ::handle<BOPDS_PaveBlock> PB;
  double                       SquareDist; // square distance from vertex to the paveblock
  double                       SinAngle;   // sinus of angle between projection vector
  // and tangent at projection point
  double Tolerance; // tolerance of the section curve
};
} // namespace

void BOPAlgo_PaveFiller::FilterPavesOnCurves(const NCollection_DynamicArray<BOPDS_Curve>& theVNC,
                                             NCollection_DataMap<int, double>&            theMVTol)
{
  // For each vertex found in ExtPaves of pave blocks of section curves
  // collect list of pave blocks with distance to the curve
  NCollection_IndexedDataMap<int, NCollection_List<PaveBlockDist>> aIDMVertPBs;
  int                                                              i;
  const double                                                     anEps = gp::Resolution();
  for (i = 0; i < theVNC.Length(); ++i)
  {
    const BOPDS_Curve&    aNC     = theVNC(i);
    const IntTools_Curve& aIC     = aNC.Curve();
    const double          aTolR3D = std::max(aNC.Tolerance(), aNC.TangentialTolerance());
    GeomAdaptor_Curve     aGAC(aIC.Curve());
    const occ::handle<BOPDS_PaveBlock>&    aPB    = aNC.PaveBlocks().First();
    const NCollection_List<BOPDS_Pave>&    aPaves = aPB->ExtPaves();
    NCollection_List<BOPDS_Pave>::Iterator itPaves(aPaves);
    for (; itPaves.More(); itPaves.Next())
    {
      const BOPDS_Pave&    aPave = itPaves.Value();
      int                  nV    = aPave.Index();
      const TopoDS_Vertex& aV    = TopoDS::Vertex(myDS->Shape(nV));
      // compute distance from vertex to the point on curve with vertex parameter
      gp_Pnt aPV  = BRep_Tool::Pnt(aV);
      double aPar = aPave.Parameter();
      gp_Pnt aPonC;
      gp_Vec aD1;
      aGAC.D1(aPar, aPonC, aD1);
      gp_Vec aProjVec(aPV, aPonC);
      double aSqDist  = aProjVec.SquareMagnitude();
      double aSqD1Mod = aD1.SquareMagnitude();
      double aSin     = aProjVec.CrossSquareMagnitude(aD1);
      if (aSqDist > anEps && aSqD1Mod > anEps)
      {
        aSin = sqrt(aSin / aSqDist / aSqD1Mod);
      }
      NCollection_List<PaveBlockDist>* pList = aIDMVertPBs.ChangeSeek(nV);
      if (!pList)
      {
        pList =
          &aIDMVertPBs.ChangeFromIndex(aIDMVertPBs.Add(nV, NCollection_List<PaveBlockDist>()));
      }
      PaveBlockDist aPBD = {aPB, aSqDist, aSin, aTolR3D};
      pList->Append(aPBD);
    }
  }

  // Process each vertex
  const double aSinAngleMin = 0.5; // angle below which projection is suspicious
  for (i = 1; i <= aIDMVertPBs.Extent(); i++)
  {
    int                                    nV    = aIDMVertPBs.FindKey(i);
    const NCollection_List<PaveBlockDist>& aList = aIDMVertPBs(i);
    // Find a pave with minimal distance
    double                                    aMinDist = RealLast();
    occ::handle<BOPDS_PaveBlock>              aPBMinDist;
    NCollection_List<PaveBlockDist>::Iterator itL(aList);
    for (; itL.More(); itL.Next())
    {
      const PaveBlockDist& aPBD = itL.Value();
      if (aPBD.SquareDist < aMinDist)
      {
        aMinDist   = aPBD.SquareDist;
        aPBMinDist = aPBD.PB;
      }
    }
    // Remove a vertex from a pave block if the distance is greater than the tolerance
    // and there are other pave blocks for which the distance is less than the current.
    // Do not remove a vertex if it is projected on the curve with quite large angle
    // (see test bugs modalg_6 bug27761).

    // Reduce tolerance for the vertex to the value of maximal distance to
    // to section curve on which it will be kept.
    double aMaxDistKept = -1;
    bool   isRemoved    = false;
    for (itL.Init(aList); itL.More(); itL.Next())
    {
      const PaveBlockDist& aPBD       = itL.Value();
      double               aCheckDist = 100. * std::max(aPBD.Tolerance * aPBD.Tolerance, aMinDist);
      if (aPBD.SquareDist > aCheckDist && aPBD.SinAngle < aSinAngleMin)
      {
        aPBD.PB->RemoveExtPave(nV);
        isRemoved = true;
      }
      else if (aPBD.SquareDist > aMaxDistKept)
      {
        aMaxDistKept = aPBD.SquareDist;
      }
    }

    if (isRemoved && aMaxDistKept > 0)
    {
      const double* pTol = theMVTol.Seek(nV);
      if (pTol)
      {
        const TopoDS_Vertex& aV = *(TopoDS_Vertex*)&myDS->Shape(nV);
        const double aRealTol   = std::max(*pTol, sqrt(aMaxDistKept) + Precision::Confusion());
        (*(occ::handle<BRep_TVertex>*)&aV.TShape())->Tolerance(aRealTol);
      }
    }
  }
}

//=================================================================================================

bool BOPAlgo_PaveFiller::ExtendedTolerance(const int                   nV,
                                           const NCollection_Map<int>& aMI,
                                           double&                     aTolVExt,
                                           const int                   aType)
{
  bool bFound = false;
  if (!(myDS->IsNewShape(nV)))
  {
    return bFound;
  }
  //
  int           i, k, aNbLines, aNbInt;
  double        aT11, aT12, aD1, aD2, aD;
  TopoDS_Vertex aV;
  gp_Pnt        aPV, aP11, aP12;
  //
  k      = 0;
  aNbInt = 2;
  if (aType == 1)
  {
    aNbInt = 1;
  }
  else if (aType == 2)
  {
    k = 1;
  }
  //
  aV  = (*(TopoDS_Vertex*)(&myDS->Shape(nV)));
  aPV = BRep_Tool::Pnt(aV);
  //
  NCollection_DynamicArray<BOPDS_InterfEE>& aEEs = myDS->InterfEE();
  NCollection_DynamicArray<BOPDS_InterfEF>& aEFs = myDS->InterfEF();
  //
  for (; k < aNbInt; ++k)
  {
    aNbLines = !k ? aEEs.Length() : aEFs.Length();
    for (i = 0; i < aNbLines; ++i)
    {
      BOPDS_Interf* aInt = !k ? (BOPDS_Interf*)(&aEEs(i)) : (BOPDS_Interf*)(&aEFs(i));
      if (aInt->IndexNew() == nV)
      {
        if (aMI.Contains(aInt->Index1()) && aMI.Contains(aInt->Index2()))
        {
          const IntTools_CommonPrt& aComPrt = !k ? aEEs(i).CommonPart() : aEFs(i).CommonPart();
          //
          const TopoDS_Edge& aE1 = aComPrt.Edge1();
          aComPrt.Range1(aT11, aT12);
          BOPTools_AlgoTools::PointOnEdge(aE1, aT11, aP11);
          BOPTools_AlgoTools::PointOnEdge(aE1, aT12, aP12);
          aD1 = aPV.Distance(aP11);
          aD2 = aPV.Distance(aP12);
          aD  = (aD1 > aD2) ? aD1 : aD2;
          if (aD > aTolVExt)
          {
            aTolVExt = aD;
          }
          return !bFound;
        } // if (aMI.Contains(aEF.Index1()) && aMI.Contains(aEF.Index2())) {
      } // if (aInt->IndexNew() == nV) {
    } // for (i = 0; i < aNbLines; ++i) {
  } // for (k=0; k<2; ++k) {
  return bFound;
}

//=================================================================================================

void BOPAlgo_PaveFiller::GetEFPnts(const int                          nF1,
                                   const int                          nF2,
                                   NCollection_List<IntSurf_PntOn2S>& aListOfPnts)
{
  int                  nE, nF, nFOpposite, aNbEFs, i;
  double               U1, U2, V1, V2, f, l;
  NCollection_Map<int> aMI;
  //
  // collect indexes of all shapes from nF1 and nF2.
  GetFullShapeMap(nF1, aMI);
  GetFullShapeMap(nF2, aMI);
  //
  NCollection_DynamicArray<BOPDS_InterfEF>& aEFs = myDS->InterfEF();
  aNbEFs                                         = aEFs.Length();
  //
  for (i = 0; i < aNbEFs; ++i)
  {
    const BOPDS_InterfEF& aEF = aEFs(i);
    if (aEF.HasIndexNew())
    {
      aEF.Indices(nE, nFOpposite);
      if (aMI.Contains(nE) && aMI.Contains(nFOpposite))
      {
        const IntTools_CommonPrt& aCP        = aEF.CommonPart();
        double                    aPar       = aCP.VertexParameter1();
        const TopoDS_Edge&        aE         = (*(TopoDS_Edge*)(&myDS->Shape(nE)));
        const TopoDS_Face&        aFOpposite = (*(TopoDS_Face*)(&myDS->Shape(nFOpposite)));
        //
        const occ::handle<Geom_Curve>& aCurve = BRep_Tool::Curve(aE, f, l);
        //
        nF                                = (nFOpposite == nF1) ? nF2 : nF1;
        const TopoDS_Face&        aF      = (*(TopoDS_Face*)(&myDS->Shape(nF)));
        occ::handle<Geom2d_Curve> aPCurve = BRep_Tool::CurveOnSurface(aE, aF, f, l);
        //
        GeomAPI_ProjectPointOnSurf& aProj = myContext->ProjPS(aFOpposite);
        //
        gp_Pnt aPoint;
        aCurve->D0(aPar, aPoint);
        IntSurf_PntOn2S aPnt;
        if (!aPCurve.IsNull())
        {
          gp_Pnt2d aP2d = aPCurve->Value(aPar);
          aProj.Perform(aPoint);
          if (aProj.IsDone())
          {
            aProj.LowerDistanceParameters(U1, V1);
            if (nF == nF1)
            {
              aPnt.SetValue(aP2d.X(), aP2d.Y(), U1, V1);
            }
            else
            {
              aPnt.SetValue(U1, V1, aP2d.X(), aP2d.Y());
            }
            aListOfPnts.Append(aPnt);
          }
        }
        else
        {
          GeomAPI_ProjectPointOnSurf& aProj1 = myContext->ProjPS(aF);
          aProj1.Perform(aPoint);
          aProj.Perform(aPoint);
          if (aProj1.IsDone() && aProj.IsDone())
          {
            aProj1.LowerDistanceParameters(U1, V1);
            aProj.LowerDistanceParameters(U2, V2);
            if (nF == nF1)
            {
              aPnt.SetValue(U1, V1, U2, V2);
            }
            else
            {
              aPnt.SetValue(U2, V2, U1, V1);
            }
            aListOfPnts.Append(aPnt);
          }
        }
      }
    }
  }
}

//=================================================================================================

void BOPAlgo_PaveFiller::PutEFPavesOnCurve(const NCollection_DynamicArray<BOPDS_Curve>& theVC,
                                           const int                                    theIndex,
                                           const NCollection_Map<int>&                  aMI,
                                           const NCollection_Map<int>&                  aMVEF,
                                           NCollection_DataMap<int, double>&            aMVTol,
                                           NCollection_DataMap<int, NCollection_List<int>>& aDMVLV)
{
  if (!aMVEF.Extent())
  {
    return;
  }
  //
  const BOPDS_Curve&    aNC = theVC.Value(theIndex);
  const IntTools_Curve& aIC = aNC.Curve();
  GeomAbs_CurveType     aTypeC;
  aTypeC = aIC.Type();
  if (aTypeC != GeomAbs_BezierCurve && aTypeC != GeomAbs_BSplineCurve)
  {
    return;
  }
  //
  int                  nV;
  NCollection_Map<int> aMV;
  //
  aMV.Assign(aMVEF);
  RemoveUsedVertices(theVC, aMV);
  if (!aMV.Extent())
  {
    return;
  }
  //
  double     aDist;
  BOPDS_Pave aPave;
  //
  const occ::handle<Geom_Curve>& aC3D    = aIC.Curve();
  GeomAPI_ProjectPointOnCurve&   aProjPT = myContext->ProjPT(aC3D);
  //
  NCollection_Map<int>::Iterator aItMI;
  aItMI.Initialize(aMV);
  for (; aItMI.More(); aItMI.Next())
  {
    nV                       = aItMI.Value();
    const TopoDS_Vertex& aV  = (*(TopoDS_Vertex*)(&myDS->Shape(nV)));
    gp_Pnt               aPV = BRep_Tool::Pnt(aV);
    aProjPT.Perform(aPV);
    int aNbPoints = aProjPT.NbPoints();
    if (aNbPoints)
    {
      aDist = aProjPT.LowerDistance();
      PutPaveOnCurve(nV, aDist, aNC, aMI, aMVTol, aDMVLV);
    }
  }
}

//=================================================================================================

void BOPAlgo_PaveFiller::PutStickPavesOnCurve(
  const TopoDS_Face&                               aF1,
  const TopoDS_Face&                               aF2,
  const NCollection_Map<int>&                      aMI,
  const NCollection_DynamicArray<BOPDS_Curve>&     theVC,
  const int                                        theIndex,
  const NCollection_Map<int>&                      aMVStick,
  NCollection_DataMap<int, double>&                aMVTol,
  NCollection_DataMap<int, NCollection_List<int>>& aDMVLV)
{
  const BOPDS_Curve& aNC = theVC.Value(theIndex);
  // Get numbers of vertices assigned to the ends of the curve
  int aBndNV[2];
  getBoundPaves(myDS, aNC, aBndNV);
  if (aBndNV[0] >= 0 && aBndNV[1] >= 0)
  {
    // both curve ends already have assigned vertices
    return;
  }
  NCollection_Map<int> aMV;
  aMV.Assign(aMVStick);
  RemoveUsedVertices(theVC, aMV);
  //
  if (!aMV.Extent())
  {
    return;
  }
  //
  occ::handle<Geom_Surface> aS1 = BRep_Tool::Surface(aF1);
  occ::handle<Geom_Surface> aS2 = BRep_Tool::Surface(aF2);
  //
  const IntTools_Curve&     aIC = aNC.Curve();
  occ::handle<Geom2d_Curve> aC2D[2];
  //
  aC2D[0] = aIC.FirstCurve2d();
  aC2D[1] = aIC.SecondCurve2d();
  if (!aC2D[0].IsNull() && !aC2D[1].IsNull())
  {
    int                            nV, m, n;
    double                         aTC[2], aD, aD2, u, v, aDT2, aScPr, aDScPr;
    gp_Pnt                         aPC[2], aPV;
    gp_Dir                         aDN[2];
    gp_Pnt2d                       aP2D;
    NCollection_Map<int>::Iterator aItMI, aItMI1;
    //
    aDT2   = 2e-7;  // the rich criteria
    aDScPr = 5.e-9; // the creasing criteria
    aIC.Bounds(aTC[0], aTC[1], aPC[0], aPC[1]);
    //
    aItMI.Initialize(aMV);
    for (; aItMI.More(); aItMI.Next())
    {
      nV                      = aItMI.Value();
      const TopoDS_Vertex& aV = *((TopoDS_Vertex*)&myDS->Shape(nV));
      aPV                     = BRep_Tool::Pnt(aV);
      //
      for (m = 0; m < 2; ++m)
      {
        if (aBndNV[m] >= 0)
        {
          continue;
        }
        aD2 = aPC[m].SquareDistance(aPV);
        if (aD2 > aDT2)
        { // no rich
          continue;
        }
        //
        for (n = 0; n < 2; ++n)
        {
          occ::handle<Geom_Surface>& aS = (!n) ? aS1 : aS2;
          aC2D[n]->D0(aTC[m], aP2D);
          aP2D.Coord(u, v);
          BOPTools_AlgoTools3D::GetNormalToSurface(aS, u, v, aDN[n]);
        }
        //
        aScPr = aDN[0] * aDN[1];
        if (aScPr < 0.)
        {
          aScPr = -aScPr;
        }
        aScPr = 1. - aScPr;
        //
        if (aScPr > aDScPr)
        {
          continue;
        }
        //
        // The intersection curve aIC is vanishing curve (the crease)
        aD = sqrt(aD2);
        //
        PutPaveOnCurve(nV, aD, aNC, aMI, aMVTol, aDMVLV);
      }
    } // for (jVU=1; jVU=aNbVU; ++jVU) {
  }
}

//=================================================================================================

void BOPAlgo_PaveFiller::GetStickVertices(const int             nF1,
                                          const int             nF2,
                                          NCollection_Map<int>& aMVStick,
                                          NCollection_Map<int>& aMVEF,
                                          NCollection_Map<int>& aMI)
{
  int nS1, nS2, nVNew, aTypeInt, i;
  //
  NCollection_DynamicArray<BOPDS_InterfVV>& aVVs = myDS->InterfVV();
  NCollection_DynamicArray<BOPDS_InterfVE>& aVEs = myDS->InterfVE();
  NCollection_DynamicArray<BOPDS_InterfEE>& aEEs = myDS->InterfEE();
  NCollection_DynamicArray<BOPDS_InterfVF>& aVFs = myDS->InterfVF();
  NCollection_DynamicArray<BOPDS_InterfEF>& aEFs = myDS->InterfEF();
  //
  int aNbLines[5] = {aVVs.Length(), aVEs.Length(), aEEs.Length(), aVFs.Length(), aEFs.Length()};
  // collect indices of all shapes from nF1 and nF2.
  aMI.Clear();
  GetFullShapeMap(nF1, aMI);
  GetFullShapeMap(nF2, aMI);
  //
  // collect VV, VE, EE, VF interferences
  for (aTypeInt = 0; aTypeInt < 4; ++aTypeInt)
  {
    for (i = 0; i < aNbLines[aTypeInt]; ++i)
    {
      BOPDS_Interf* aInt = (aTypeInt == 0)
                             ? (BOPDS_Interf*)(&aVVs(i))
                             : ((aTypeInt == 1) ? (BOPDS_Interf*)(&aVEs(i))
                                                : ((aTypeInt == 2) ? (BOPDS_Interf*)(&aEEs(i))
                                                                   : (BOPDS_Interf*)(&aVFs(i))));
      if (aInt->HasIndexNew())
      {
        aInt->Indices(nS1, nS2);
        if (aMI.Contains(nS1) && aMI.Contains(nS2))
        {
          nVNew = aInt->IndexNew();
          myDS->HasShapeSD(nVNew, nVNew);
          aMVStick.Add(nVNew);
        }
      }
    }
  }
  // collect EF interferences
  for (i = 0; i < aNbLines[4]; ++i)
  {
    const BOPDS_InterfEF& aInt = aEFs(i);
    if (aInt.HasIndexNew())
    {
      aInt.Indices(nS1, nS2);
      if (aMI.Contains(nS1) && aMI.Contains(nS2))
      {
        nVNew = aInt.IndexNew();
        myDS->HasShapeSD(nVNew, nVNew);
        aMVStick.Add(nVNew);
        aMVEF.Add(nVNew);
      }
    }
  }
}

//=================================================================================================

void BOPAlgo_PaveFiller::GetFullShapeMap(const int nF, NCollection_Map<int>& aMI)
{
  NCollection_List<int>::Iterator aIt;
  int                             nS;
  //
  const BOPDS_ShapeInfo&       aSI = myDS->ShapeInfo(nF);
  const NCollection_List<int>& aLI = aSI.SubShapes();
  //
  aMI.Add(nF);
  aIt.Initialize(aLI);
  for (; aIt.More(); aIt.Next())
  {
    nS = aIt.Value();
    aMI.Add(nS);
  }
}

//=================================================================================================

void BOPAlgo_PaveFiller::RemoveUsedVertices(const NCollection_DynamicArray<BOPDS_Curve>& aVC,
                                            NCollection_Map<int>&                        aMV)
{
  if (aMV.IsEmpty())
  {
    return;
  }

  for (int i = 0; i < aVC.Length(); ++i)
  {
    const BOPDS_Curve&                                       aNC   = aVC.Value(i);
    const NCollection_List<occ::handle<BOPDS_PaveBlock>>&    aLPBC = aNC.PaveBlocks();
    NCollection_List<occ::handle<BOPDS_PaveBlock>>::Iterator itPB(aLPBC);
    for (; itPB.More(); itPB.Next())
    {
      const occ::handle<BOPDS_PaveBlock>&    aPB = itPB.Value();
      const NCollection_List<BOPDS_Pave>&    aLP = aPB->ExtPaves();
      NCollection_List<BOPDS_Pave>::Iterator itLP(aLP);
      for (; itLP.More(); itLP.Next())
      {
        aMV.Remove(itLP.Value().Index());
      }

      aMV.Remove(aPB->Pave1().Index());
      aMV.Remove(aPB->Pave2().Index());
    }
  }
}

//=================================================================================================

void BOPAlgo_PaveFiller::PutPaveOnCurve(const int                                        nV,
                                        const double                                     aTolR3D,
                                        const BOPDS_Curve&                               aNC,
                                        const NCollection_Map<int>&                      aMI,
                                        NCollection_DataMap<int, double>&                aMVTol,
                                        NCollection_DataMap<int, NCollection_List<int>>& aDMVLV,
                                        const int iCheckExtend)
{
  bool   bIsVertexOnLine;
  double aT;
  //
  const TopoDS_Vertex&                aV  = (*(TopoDS_Vertex*)(&myDS->Shape(nV)));
  const occ::handle<BOPDS_PaveBlock>& aPB = aNC.PaveBlocks().First();
  const IntTools_Curve&               aIC = aNC.Curve();
  //
  double aTolV = (aMVTol.IsBound(nV) ? aMVTol(nV) : BRep_Tool::Tolerance(aV));

  bIsVertexOnLine = myContext->IsVertexOnLine(aV, aTolV, aIC, aTolR3D + myFuzzyValue, aT);
  if (!bIsVertexOnLine && iCheckExtend && !myVertsToAvoidExtension.Contains(nV))
  {
    double anExtraTol = aTolV;
    if (ExtendedTolerance(nV, aMI, anExtraTol, iCheckExtend))
    {
      bIsVertexOnLine = myContext->IsVertexOnLine(aV, anExtraTol, aIC, aTolR3D + myFuzzyValue, aT);
      if (bIsVertexOnLine)
      {
        gp_Pnt aPOnC;
        aIC.D0(aT, aPOnC);
        aTolV = aPOnC.Distance(BRep_Tool::Pnt(aV));
      }
    }
  }
  //
  if (bIsVertexOnLine)
  {
    // check if aPB contains the parameter aT
    bool   bExist;
    int    nVUsed;
    double aPTol, aDTol;
    //
    aDTol = BOPTools_AlgoTools::DTolerance();
    //
    GeomAdaptor_Curve aGAC(aIC.Curve());
    aPTol = aGAC.Resolution(std::max(aTolR3D, aTolV));
    //
    bExist = aPB->ContainsParameter(aT, aPTol, nVUsed);
    if (bExist)
    {
      // use existing pave
      NCollection_List<int>* pList = aDMVLV.ChangeSeek(nVUsed);
      if (!pList)
      {
        pList = aDMVLV.Bound(nVUsed, NCollection_List<int>());
        pList->Append(nVUsed);
        if (!aMVTol.IsBound(nVUsed))
        {
          const TopoDS_Vertex& aVUsed = (*(TopoDS_Vertex*)(&myDS->Shape(nVUsed)));
          aTolV                       = BRep_Tool::Tolerance(aVUsed);
          aMVTol.Bind(nVUsed, aTolV);
        }
      }
      // avoid repeated elements in the list
      NCollection_List<int>::Iterator aItLI(*pList);
      for (; aItLI.More(); aItLI.Next())
      {
        if (aItLI.Value() == nV)
        {
          break;
        }
      }
      if (!aItLI.More())
      {
        pList->Append(nV);
      }
      // save initial tolerance for the vertex
      if (!aMVTol.IsBound(nV))
      {
        aTolV = BRep_Tool::Tolerance(aV);
        aMVTol.Bind(nV, aTolV);
      }
    }
    else
    {
      // add new pave
      BOPDS_Pave aPave;
      aPave.SetIndex(nV);
      aPave.SetParameter(aT);
      aPB->AppendExtPave(aPave);
      //
      gp_Pnt aP1   = aGAC.Value(aT);
      aTolV        = BRep_Tool::Tolerance(aV);
      gp_Pnt aP2   = BRep_Tool::Pnt(aV);
      double aDist = aP1.Distance(aP2);
      if (aTolV < aDist + aDTol)
      {
        BRep_Builder().UpdateVertex(aV, aDist + aDTol);
        //
        if (!aMVTol.IsBound(nV))
        {
          aMVTol.Bind(nV, aTolV);
        }
        //
        BOPDS_ShapeInfo& aSIDS  = myDS->ChangeShapeInfo(nV);
        Bnd_Box&         aBoxDS = aSIDS.ChangeBox();
        BRepBndLib::Add(aV, aBoxDS);
        aBoxDS.SetGap(aBoxDS.GetGap() + Precision::Confusion());
      }
    }
  }
}

//=================================================================================================

void BOPAlgo_PaveFiller::ProcessExistingPaveBlocks(
  const int                                                   theInt,
  const int                                                   theCur,
  const int                                                   nF1,
  const int                                                   nF2,
  const TopoDS_Edge&                                          theES,
  const NCollection_IndexedMap<occ::handle<BOPDS_PaveBlock>>& theMPBOnIn,
  BOPTools_BoxTree&                                           thePBTree,
  NCollection_IndexedDataMap<TopoDS_Shape, BOPDS_CoupleOfPaveBlocks, TopTools_ShapeMapHasher>&
                                                                   theMSCPB,
  NCollection_DataMap<TopoDS_Shape, int, TopTools_ShapeMapHasher>& theMVI,
  NCollection_List<occ::handle<BOPDS_PaveBlock>>&                  theLPBC,
  BOPAlgo_DataMapOfPaveBlockListOfInteger&                         thePBFacesMap,
  NCollection_Map<occ::handle<BOPDS_PaveBlock>>&                   theMPB)
{
  Bnd_Box aBoxES;
  BRepBndLib::Add(theES, aBoxES, false);

  BOPTools_BoxTreeSelector aSelector;
  aSelector.SetBox(Bnd_Tools::Bnd2BVH(aBoxES));
  aSelector.SetBVHSet(&thePBTree);
  if (!aSelector.Select())
  {
    return;
  }

  const double aTolES = BRep_Tool::Tolerance(theES);

  const BOPDS_FaceInfo& aFI1 = myDS->FaceInfo(nF1);
  const BOPDS_FaceInfo& aFI2 = myDS->FaceInfo(nF2);

  for (NCollection_List<int>::Iterator itPB(aSelector.Indices()); itPB.More(); itPB.Next())
  {
    const occ::handle<BOPDS_PaveBlock>& aPBF = theMPBOnIn(itPB.Value());
    if (theMPB.Contains(aPBF))
    {
      continue;
    }

    bool bInF1 = (aFI1.PaveBlocksOn().Contains(aPBF) || aFI1.PaveBlocksIn().Contains(aPBF));
    bool bInF2 = (aFI2.PaveBlocksOn().Contains(aPBF) || aFI2.PaveBlocksIn().Contains(aPBF));
    if (bInF1 && bInF2)
    {
      // Add all common edges for post treatment
      theMPB.Add(aPBF);
      PreparePostTreatFF(theInt, theCur, aPBF, theMSCPB, theMVI, theLPBC);
      continue;
    }

    const int                                  nF = bInF1 ? nF2 : nF1;
    const NCollection_List<EdgeRangeDistance>* pList =
      myDistances.Seek(BOPDS_Pair(aPBF->OriginalEdge(), nF));
    if (!pList)
    {
      continue;
    }

    double aT1, aT2;
    aPBF->Range(aT1, aT2);

    double aDist = RealLast();
    for (NCollection_List<EdgeRangeDistance>::Iterator itR(*pList); itR.More(); itR.Next())
    {
      const EdgeRangeDistance& aRangeDist = itR.Value();
      if ((aT1 <= aRangeDist.First && aRangeDist.First <= aT2)
          || (aT1 <= aRangeDist.Last && aRangeDist.Last <= aT2)
          || (aRangeDist.First <= aT1 && aT1 <= aRangeDist.Last)
          || (aRangeDist.First <= aT2 && aT2 <= aRangeDist.Last))
      {
        aDist = aRangeDist.Distance;
        break;
      }
    }
    if (aDist < RealLast())
    {
      const TopoDS_Edge& aEF     = TopoDS::Edge(myDS->Shape(aPBF->Edge()));
      const double       aTolSum = aTolES + BRep_Tool::Tolerance(aEF);

      if (aDist <= aTolSum)
      {
        theMPB.Add(aPBF);
        PreparePostTreatFF(theInt, theCur, aPBF, theMSCPB, theMVI, theLPBC);

        NCollection_List<int>* pFaces = thePBFacesMap.ChangeSeek(aPBF);
        if (!pFaces)
        {
          pFaces = thePBFacesMap.Bound(aPBF, NCollection_List<int>());
        }
        if (pFaces->IsEmpty() || !pFaces->Contains(nF))
        {
          pFaces->Append(nF);
        }
      }
    }
  }
}

//=================================================================================================

void BOPAlgo_PaveFiller::ProcessExistingPaveBlocks(
  const int                                                   theInt,
  const int                                                   nF1,
  const int                                                   nF2,
  const NCollection_IndexedMap<occ::handle<BOPDS_PaveBlock>>& aMPBOnIn,
  BOPTools_BoxTree&                                           thePBTree,
  const NCollection_DataMap<int, NCollection_List<int>>&      aDMBV,
  NCollection_IndexedDataMap<TopoDS_Shape, BOPDS_CoupleOfPaveBlocks, TopTools_ShapeMapHasher>&
                                                                   aMSCPB,
  NCollection_DataMap<TopoDS_Shape, int, TopTools_ShapeMapHasher>& aMVI,
  BOPAlgo_DataMapOfPaveBlockListOfInteger&                         thePBFacesMap,
  NCollection_Map<occ::handle<BOPDS_PaveBlock>>&                   aMPB)
{
  if (aDMBV.IsEmpty())
  {
    return;
  }
  //
  double                                                    aT, dummy;
  int                                                       nV, nE, iC, iFlag;
  NCollection_List<int>::Iterator                           aItLI;
  NCollection_DataMap<int, NCollection_List<int>>::Iterator aItBV;
  //
  NCollection_DynamicArray<BOPDS_InterfFF>& aFFs = myDS->InterfFF();
  BOPDS_InterfFF&                           aFF  = aFFs(theInt);
  NCollection_DynamicArray<BOPDS_Curve>&    aVC  = aFF.ChangeCurves();
  //
  const BOPDS_FaceInfo& aFI1 = myDS->FaceInfo(nF1);
  const BOPDS_FaceInfo& aFI2 = myDS->FaceInfo(nF2);
  //
  aItBV.Initialize(aDMBV);
  for (; aItBV.More(); aItBV.Next())
  {
    iC                                = aItBV.Key();
    const NCollection_List<int>& aLBV = aItBV.Value();
    //
    BOPDS_Curve&                                    aNC   = aVC.ChangeValue(iC);
    NCollection_List<occ::handle<BOPDS_PaveBlock>>& aLPBC = aNC.ChangePaveBlocks();
    //
    aItLI.Initialize(aLBV);
    for (; aItLI.More(); aItLI.Next())
    {
      nV                           = aItLI.Value();
      const BOPDS_ShapeInfo& aSIV  = myDS->ShapeInfo(nV);
      const Bnd_Box&         aBoxV = aSIV.Box();
      const TopoDS_Vertex&   aV    = *(TopoDS_Vertex*)&aSIV.Shape();
      if (!aMVI.IsBound(aV))
      {
        continue;
      }
      //
      BOPTools_BoxTreeSelector aSelector;
      aSelector.SetBox(Bnd_Tools::Bnd2BVH(aBoxV));
      aSelector.SetBVHSet(&thePBTree);
      if (!aSelector.Select())
      {
        continue;
      }

      for (NCollection_List<int>::Iterator it(aSelector.Indices()); it.More(); it.Next())
      {
        const occ::handle<BOPDS_PaveBlock>& aPB = aMPBOnIn(it.Value());
        if (aPB->Pave1().Index() == nV || aPB->Pave2().Index() == nV)
        {
          continue;
        }
        //
        if (aMPB.Contains(aPB))
        {
          continue;
        }
        nE                          = aPB->Edge();
        const BOPDS_ShapeInfo& aSIE = myDS->ShapeInfo(nE);
        const TopoDS_Edge&     aE   = *(TopoDS_Edge*)&aSIE.Shape();
        //
        iFlag = myContext->ComputeVE(aV, aE, aT, dummy, myFuzzyValue);
        if (!iFlag)
        {
          aMPB.Add(aPB);
          PreparePostTreatFF(theInt, iC, aPB, aMSCPB, aMVI, aLPBC);

          // Add faces to PB
          bool bInF1 = (aFI1.PaveBlocksOn().Contains(aPB) || aFI1.PaveBlocksIn().Contains(aPB));
          bool bInF2 = (aFI2.PaveBlocksOn().Contains(aPB) || aFI2.PaveBlocksIn().Contains(aPB));
          if (!bInF1 || !bInF2)
          {
            // Face without pave block
            const int              nF     = bInF1 ? nF2 : nF1;
            NCollection_List<int>* pFaces = thePBFacesMap.ChangeSeek(aPB);
            if (!pFaces)
            {
              pFaces = thePBFacesMap.Bound(aPB, NCollection_List<int>());
            }
            // List is expected to be short, so we allow the check here
            if (pFaces->IsEmpty() || !pFaces->Contains(nF))
            {
              pFaces->Append(nF);
            }
          }
        }
      }
    }
  }
}

//=================================================================================================

void BOPAlgo_PaveFiller::UpdateExistingPaveBlocks(
  const occ::handle<BOPDS_PaveBlock>&             aPBf,
  NCollection_List<occ::handle<BOPDS_PaveBlock>>& aLPB,
  const BOPAlgo_DataMapOfPaveBlockListOfInteger&  thePBFacesMap)
{
  if (!aLPB.Extent())
  {
    return;
  }
  //
  int                                                      nE;
  bool                                                     bCB;
  occ::handle<BOPDS_PaveBlock>                             aPB, aPB1, aPB2, aPB2n;
  occ::handle<BOPDS_CommonBlock>                           aCB;
  NCollection_List<occ::handle<BOPDS_PaveBlock>>::Iterator aIt, aIt1, aIt2;
  //
  // 1. Remove old pave blocks
  const occ::handle<BOPDS_CommonBlock>& aCB1 = myDS->CommonBlock(aPBf);
  bCB                                        = !aCB1.IsNull();
  NCollection_List<occ::handle<BOPDS_PaveBlock>> aLPB1;
  //
  if (bCB)
  {
    aLPB1.Assign(aCB1->PaveBlocks());
  }
  else
  {
    aLPB1.Append(aPBf);
  }
  aIt1.Initialize(aLPB1);
  for (; aIt1.More(); aIt1.Next())
  {
    aPB1 = aIt1.Value();
    nE   = aPB1->OriginalEdge();
    //
    NCollection_List<occ::handle<BOPDS_PaveBlock>>& aLPB2 = myDS->ChangePaveBlocks(nE);
    aIt2.Initialize(aLPB2);
    for (; aIt2.More(); aIt2.Next())
    {
      aPB2 = aIt2.Value();
      if (aPB1 == aPB2)
      {
        aLPB2.Remove(aIt2);
        break;
      }
    }
  }
  //
  // 2. Update pave blocks
  if (bCB)
  {
    // Create new common blocks
    NCollection_List<occ::handle<BOPDS_PaveBlock>> aLPBNew;
    const NCollection_List<int>&                   aFaces = aCB1->Faces();
    aIt.Initialize(aLPB);
    for (; aIt.More(); aIt.Next())
    {
      const occ::handle<BOPDS_PaveBlock>& aPBValue         = aIt.Value();
      BOPDS_Pave                          aPBValuePaves[2] = {aPBValue->Pave1(), aPBValue->Pave2()};
      //
      aCB = new BOPDS_CommonBlock;
      aIt1.Initialize(aLPB1);
      for (; aIt1.More(); aIt1.Next())
      {
        aPB2 = aIt1.Value();
        nE   = aPB2->OriginalEdge();
        //
        // Create new pave block
        aPB2n = new BOPDS_PaveBlock;
        if (aPBValue->OriginalEdge() == nE)
        {
          aPB2n->SetPave1(aPBValuePaves[0]);
          aPB2n->SetPave2(aPBValuePaves[1]);
        }
        else
        {
          // For the different original edge compute the parameters of paves
          BOPDS_Pave aPave[2];

          if (aPBValuePaves[0].Index() == aPBValuePaves[1].Index()
              && aPB2->Pave1().Index() == aPB2->Pave2().Index())
          {
            // still closed
            aPave[0].SetIndex(aPBValuePaves[0].Index());
            aPave[0].SetParameter(aPB2->Pave1().Parameter());
            aPave[1].SetIndex(aPBValuePaves[1].Index());
            aPave[1].SetParameter(aPB2->Pave2().Parameter());
          }
          else
          {
            for (int i = 0; i < 2; ++i)
            {
              int nV = aPBValuePaves[i].Index();
              aPave[i].SetIndex(nV);
              if (nV == aPB2->Pave1().Index())
              {
                aPave[i].SetParameter(aPB2->Pave1().Parameter());
              }
              else if (nV == aPB2->Pave2().Index())
              {
                aPave[i].SetParameter(aPB2->Pave2().Parameter());
              }
              else
              {
                // Compute the parameter by projecting the point
                const TopoDS_Vertex& aV   = TopoDS::Vertex(myDS->Shape(nV));
                const TopoDS_Edge&   aEOr = TopoDS::Edge(myDS->Shape(nE));
                double               aTOut, aDist;
                int iErr = myContext->ComputeVE(aV, aEOr, aTOut, aDist, myFuzzyValue);
                if (!iErr)
                {
                  aPave[i].SetParameter(aTOut);
                }
                else
                {
                  // Unable to project - set the parameter of the closest boundary
                  const TopoDS_Vertex& aV1 = TopoDS::Vertex(myDS->Shape(aPB2->Pave1().Index()));
                  const TopoDS_Vertex& aV2 = TopoDS::Vertex(myDS->Shape(aPB2->Pave2().Index()));
                  //
                  gp_Pnt aP  = BRep_Tool::Pnt(aV);
                  gp_Pnt aP1 = BRep_Tool::Pnt(aV1);
                  gp_Pnt aP2 = BRep_Tool::Pnt(aV2);
                  //
                  double aDist1 = aP.SquareDistance(aP1);
                  double aDist2 = aP.SquareDistance(aP2);
                  //
                  aPave[i].SetParameter(aDist1 < aDist2 ? aPB2->Pave1().Parameter()
                                                        : aPB2->Pave2().Parameter());
                }
              }
            }
          }
          //
          if (aPave[1].Parameter() < aPave[0].Parameter())
          {
            BOPDS_Pave aPaveTmp = aPave[0];
            aPave[0]            = aPave[1];
            aPave[1]            = aPaveTmp;
          }
          //
          aPB2n->SetPave1(aPave[0]);
          aPB2n->SetPave2(aPave[1]);
        }
        //
        aPB2n->SetEdge(aPBValue->Edge());
        aPB2n->SetOriginalEdge(nE);
        aCB->AddPaveBlock(aPB2n);
        myDS->SetCommonBlock(aPB2n, aCB);
        myDS->ChangePaveBlocks(nE).Append(aPB2n);
      }
      aCB->SetFaces(aFaces);
      //
      const occ::handle<BOPDS_PaveBlock>& aPBNew = aCB->PaveBlocks().First();
      aLPBNew.Append(aPBNew);
    }
    //
    aLPB = aLPBNew;
  }
  else
  {
    nE                                                    = aPBf->OriginalEdge();
    NCollection_List<occ::handle<BOPDS_PaveBlock>>& aLPBE = myDS->ChangePaveBlocks(nE);
    aIt.Initialize(aLPB);
    for (; aIt.More(); aIt.Next())
    {
      aPB = aIt.Value();
      aLPBE.Append(aPB);
    }
  }

  // Try to project the edge on the faces
  const NCollection_List<int>* pLFaces = thePBFacesMap.Seek(aPBf);
  if (!pLFaces)
  {
    return;
  }
  NCollection_List<int>::Iterator itLF(*pLFaces);
  for (; itLF.More(); itLF.Next())
  {
    const int          nF  = itLF.Value();
    BOPDS_FaceInfo&    aFI = myDS->ChangeFaceInfo(nF);
    const TopoDS_Face& aF  = TopoDS::Face(myDS->Shape(nF));

    aIt.Initialize(aLPB);
    for (; aIt.More(); aIt.Next())
    {
      aPB = aIt.ChangeValue();
      if (aFI.PaveBlocksOn().Contains(aPB) || aFI.PaveBlocksIn().Contains(aPB))
      {
        continue;
      }

      const TopoDS_Edge& aE = *(TopoDS_Edge*)&myDS->Shape(aPB->Edge());
      //
      IntTools_EdgeFace anEF;
      anEF.SetEdge(aE);
      anEF.SetFace(aF);
      anEF.SetFuzzyValue(myFuzzyValue);
      anEF.SetRange(aPB->Pave1().Parameter(), aPB->Pave2().Parameter());
      anEF.SetContext(myContext);
      anEF.Perform();
      //
      const NCollection_Sequence<IntTools_CommonPrt>& aCPrts = anEF.CommonParts();
      bool bCoincide = (aCPrts.Length() == 1 && aCPrts(1).Type() == TopAbs_EDGE);
      if (bCoincide)
      {
        aCB = myDS->CommonBlock(aPB);
        if (aCB.IsNull())
        {
          aCB = new BOPDS_CommonBlock;
          aCB->AddPaveBlock(aPB);
          myDS->SetCommonBlock(aPB, aCB);
        }
        aCB->AddFace(nF);
        aFI.ChangePaveBlocksIn().Add(aPB);
      }
    }
  }
}

//=================================================================================================

void BOPAlgo_PaveFiller::PutClosingPaveOnCurve(BOPDS_Curve& aNC)
{
  const IntTools_Curve&          aIC  = aNC.Curve();
  const occ::handle<Geom_Curve>& aC3D = aIC.Curve();
  // check 3d curve
  if (aC3D.IsNull())
  {
    return;
  }

  // check bounds
  if (!aIC.HasBounds())
  {
    return;
  }

  // check closeness
  double aT[2];
  gp_Pnt aP[2];
  aIC.Bounds(aT[0], aT[1], aP[0], aP[1]);

  // Find the pave which has been put at one of the ends
  BOPDS_Pave aPave;
  // Index of the vertex put at one of the ends
  int nV = -1;
  // Keep the opposite parameter
  double aTOp = 0.;
  // Keep the opposite bounding point
  gp_Pnt aPOp;

  occ::handle<BOPDS_PaveBlock>&          aPB = aNC.ChangePaveBlock1();
  NCollection_List<BOPDS_Pave>&          aLP = aPB->ChangeExtPaves();
  NCollection_List<BOPDS_Pave>::Iterator aItLP(aLP);
  for (; aItLP.More() && (nV < 0); aItLP.Next())
  {
    aPave      = aItLP.Value();
    double aTC = aPave.Parameter();
    for (int j = 0; j < 2; ++j)
    {
      if (std::abs(aTC - aT[j]) < Precision::PConfusion())
      {
        nV   = aPave.Index();
        aTOp = (!j) ? aT[1] : aT[0];
        aPOp = (!j) ? aP[1] : aP[0];
        break;
      }
    }
  }

  if (nV < 0)
  {
    // A vertex at an end in space, put on the curve a little way in: a
    // cylinder's seam 1e-8 from where its circle section starts is 2e-7 from
    // the circle's end, within the vertex's tolerance, so no vertex was made
    // at the end -- the seam's stands for it, at 1e-8, and no pave was found
    // at a bound. Not closed, the circle had one pave and no block: the
    // section was lost. The closing pave goes at the opposite bound.
    const double aTolC = std::max(aNC.Tolerance(), aNC.TangentialTolerance()) + Precision::Confusion();
    for (aItLP.Initialize(aLP); aItLP.More() && (nV < 0); aItLP.Next())
    {
      const BOPDS_Pave&    aPv  = aItLP.Value();
      const TopoDS_Vertex& aVx  = TopoDS::Vertex(myDS->Shape(aPv.Index()));
      const gp_Pnt         aPVx = BRep_Tool::Pnt(aVx);
      const double         aTol = BRep_Tool::Tolerance(aVx) + aTolC;
      for (int j = 0; j < 2; ++j)
      {
        if (aPVx.Distance(aP[j]) <= aTol)
        {
          aPave = aPv;
          nV    = aPv.Index();
          aTOp  = (!j) ? aT[1] : aT[0];
          aPOp  = (!j) ? aP[1] : aP[0];
          break;
        }
      }
    }
    // a pave at the opposite bound already: nothing to close
    for (aItLP.Initialize(aLP); aItLP.More() && nV >= 0; aItLP.Next())
    {
      if (std::abs(aItLP.Value().Parameter() - aTOp) < Precision::PConfusion())
      {
        nV = -1;
      }
    }
  }

  if (nV < 0)
  {
    // No paves on the bounds of the curve
    return;
  }

  // Check if the curve is closed using the tolerance
  // of found vertex
  const TopoDS_Vertex& aV    = TopoDS::Vertex(myDS->Shape(nV));
  double               aTolV = BRep_Tool::Tolerance(aV);
  gp_Pnt               aPV   = BRep_Tool::Pnt(aV);
  // Tolerance for the point on the curve
  double aTolP = std::max(aNC.Tolerance(), aNC.TangentialTolerance());
  aTolP += Precision::Confusion();

  const double aDistVP = aPV.Distance(aPOp);
  if (aDistVP > aTolV + aTolP)
  {
    // Curve is not closed
    return;
  }

  // Check if there will be valid range on the curve
  double aFirst, aLast;
  double aNewTolV = std::max(aTolV, aDistVP + BOPTools_AlgoTools::DTolerance());
  if (!BRepLib::FindValidRange(GeomAdaptor_Curve(aIC.Curve()),
                               aIC.Tolerance(),
                               aT[0],
                               aP[0],
                               aNewTolV,
                               aT[1],
                               aP[1],
                               aNewTolV,
                               aFirst,
                               aLast))
  {
    // No valid range
    return;
  }

  if (aNewTolV > aTolV)
  {
    int nVn = UpdateVertex(nV, aNewTolV);
    if (nVn != nV)
    {
      aPave.SetIndex(nVn);
      nV = nVn;
    }
    aTolV = BRep_Tool::Tolerance(TopoDS::Vertex(myDS->Shape(nV)));
  }

  // Add closing pave to the curve
  BOPDS_Pave aNewPave;
  aNewPave.SetIndex(nV);
  aNewPave.SetParameter(aTOp);
  aLP.Append(aNewPave);
}

//=================================================================================================

void BOPAlgo_PaveFiller::PreparePostTreatFF(
  const int                           aInt,
  const int                           aCur,
  const occ::handle<BOPDS_PaveBlock>& aPB,
  NCollection_IndexedDataMap<TopoDS_Shape, BOPDS_CoupleOfPaveBlocks, TopTools_ShapeMapHasher>&
                                                                   aMSCPB,
  NCollection_DataMap<TopoDS_Shape, int, TopTools_ShapeMapHasher>& aMVI,
  NCollection_List<occ::handle<BOPDS_PaveBlock>>&                  aLPBC)
{
  int nV1, nV2;
  //
  aLPBC.Append(aPB);
  //
  aPB->Indices(nV1, nV2);
  const TopoDS_Vertex& aV1 = (*(TopoDS_Vertex*)(&myDS->Shape(nV1)));
  const TopoDS_Vertex& aV2 = (*(TopoDS_Vertex*)(&myDS->Shape(nV2)));
  const TopoDS_Edge&   aE  = *(TopoDS_Edge*)&myDS->Shape(aPB->Edge());
  // Keep info for post treatment
  BOPDS_CoupleOfPaveBlocks aCPB;
  aCPB.SetIndexInterf(aInt);
  aCPB.SetIndex(aCur);
  aCPB.SetPaveBlock1(aPB);
  //
  aMSCPB.Add(aE, aCPB);
  aMVI.Bind(aV1, nV1);
  aMVI.Bind(aV2, nV2);
}

//=================================================================================================

bool BOPAlgo_PaveFiller::CheckPlanes(const int nF1, const int nF2) const
{
  bool                           bToIntersect;
  int                            i, nV2, iCnt;
  NCollection_Map<int>::Iterator aIt;
  //
  bToIntersect = false;
  //
  const BOPDS_FaceInfo& aFI1 = myDS->ChangeFaceInfo(nF1);
  const BOPDS_FaceInfo& aFI2 = myDS->ChangeFaceInfo(nF2);
  //
  const NCollection_Map<int>& aMVIn1 = aFI1.VerticesIn();
  const NCollection_Map<int>& aMVOn1 = aFI1.VerticesOn();
  //
  iCnt = 0;
  for (i = 0; (i < 2 && !bToIntersect); ++i)
  {
    const NCollection_Map<int>& aMV2 = (!i) ? aFI2.VerticesIn() : aFI2.VerticesOn();
    //
    aIt.Initialize(aMV2);
    for (; aIt.More(); aIt.Next())
    {
      nV2 = aIt.Value();
      if (aMVIn1.Contains(nV2) || aMVOn1.Contains(nV2))
      {
        ++iCnt;
        if (iCnt > 1)
        {
          bToIntersect = !bToIntersect;
          break;
        }
      }
    }
  }
  //
  return bToIntersect;
}

//=================================================================================================

void BOPAlgo_PaveFiller::UpdatePaveBlocks(const NCollection_DataMap<int, int>& aDMNewSD)
{
  if (aDMNewSD.IsEmpty())
  {
    return;
  }
  //
  int                                                      nSp, aNbPBP, nV[2], i, j;
  double                                                   aT[2];
  bool                                                     bCB, bRebuild;
  NCollection_List<occ::handle<BOPDS_PaveBlock>>::Iterator aItPB;
  NCollection_Map<occ::handle<BOPDS_PaveBlock>>            aMPB;
  NCollection_Map<int>                                     aMicroEdges;
  //
  NCollection_List<occ::handle<BOPDS_PaveBlock>> anAllPBs;

  // Get pave blocks of section edges
  NCollection_DynamicArray<BOPDS_InterfFF>& aFFs  = myDS->InterfFF();
  int                                       aNbFF = aFFs.Length();
  for (i = 0; i < aNbFF; ++i)
  {
    const BOPDS_InterfFF&                        aFF  = aFFs(i);
    const NCollection_DynamicArray<BOPDS_Curve>& aVNC = aFF.Curves();
    int                                          aNbC = aVNC.Length();
    for (j = 0; j < aNbC; ++j)
    {
      const BOPDS_Curve&                                    aNC   = aVNC(j);
      const NCollection_List<occ::handle<BOPDS_PaveBlock>>& aLPBC = aNC.PaveBlocks();
      aItPB.Initialize(aLPBC);
      for (; aItPB.More(); aItPB.Next())
      {
        anAllPBs.Append(aItPB.Value());
      }
    }
  }

  // Get pave blocks from the pool
  NCollection_DynamicArray<NCollection_List<occ::handle<BOPDS_PaveBlock>>>& aPBP =
    myDS->ChangePaveBlocksPool();
  aNbPBP = aPBP.Length();
  for (i = 0; i < aNbPBP; ++i)
  {
    NCollection_List<occ::handle<BOPDS_PaveBlock>>& aLPB = aPBP(i);
    aItPB.Initialize(aLPB);
    for (; aItPB.More(); aItPB.Next())
    {
      anAllPBs.Append(aItPB.Value());
    }
  }

  // Process all pave blocks
  aItPB.Initialize(anAllPBs);
  for (; aItPB.More(); aItPB.Next())
  {
    occ::handle<BOPDS_PaveBlock>          aPB = aItPB.Value();
    const occ::handle<BOPDS_CommonBlock>& aCB = myDS->CommonBlock(aPB);
    bCB                                       = !aCB.IsNull();
    if (bCB)
    {
      aPB = aCB->PaveBlock1();
    }
    //
    if (aMPB.Add(aPB))
    {
      bRebuild = false;
      aPB->Indices(nV[0], nV[1]);
      aPB->Range(aT[0], aT[1]);
      // remember the fact if the edge had different vertices before substitution
      bool wasRegularEdge = (nV[0] != nV[1]);
      //
      for (j = 0; j < 2; ++j)
      {
        if (aDMNewSD.IsBound(nV[j]))
        {
          BOPDS_Pave aPave;
          //
          nV[j] = aDMNewSD.Find(nV[j]);
          aPave.SetIndex(nV[j]);
          aPave.SetParameter(aT[j]);
          //
          bRebuild = true;
          if (!j)
          {
            aPB->SetPave1(aPave);
          }
          else
          {
            aPB->SetPave2(aPave);
          }
        }
      }
      //
      if (bRebuild)
      {
        int nE = aPB->Edge();
        // Check if the Pave Block has the edge set
        if (nE < 0)
        {
          // untouched edge
          nE = aPB->OriginalEdge();
        }
        bool isDegEdge = myDS->ShapeInfo(nE).HasFlag();
        if (wasRegularEdge && !isDegEdge && nV[0] == nV[1])
        {
          // now edge has the same vertex on both ends;
          // check if it is not a regular closed curve.
          FillShrunkData(aPB);
          if (!aPB->HasShrunkData())
          {
            // micro edge, so mark it for removal
            aMicroEdges.Add(nE);
            continue;
          }
        }
        nSp = SplitEdge(nE, nV[0], aT[0], nV[1], aT[1]);
        if (bCB)
        {
          aCB->SetEdge(nSp);
        }
        else
        {
          aPB->SetEdge(nSp);
        }
      } // if (bRebuild) {
    } // if (aMPB.Add(aPB)) {
  } // for (; aItPB.More(); aItPB.Next()) {
  aMPB.Clear();

  if (aMicroEdges.Extent())
  {
    RemovePaveBlocks(aMicroEdges);
  }
}

//=================================================================================================

void BOPAlgo_PaveFiller::RemovePaveBlocks(const NCollection_Map<int>& theEdges)
{
  // Remove all pave blocks referring to input edges:
  //
  // 1. from the Pave Blocks Pool
  NCollection_DynamicArray<NCollection_List<occ::handle<BOPDS_PaveBlock>>>& aPBP =
    myDS->ChangePaveBlocksPool();
  int aNbPBP = aPBP.Length(), i;
  for (i = 0; i < aNbPBP; ++i)
  {
    NCollection_List<occ::handle<BOPDS_PaveBlock>>& aLPB = aPBP(i);
    //
    NCollection_List<occ::handle<BOPDS_PaveBlock>>::Iterator aItPB(aLPB);
    while (aItPB.More())
    {
      const occ::handle<BOPDS_PaveBlock>& aPB = aItPB.Value();
      if (theEdges.Contains(aPB->Edge()))
      {
        aLPB.Remove(aItPB);
      }
      else
      {
        aItPB.Next();
      }
    }
  }

  // 2. from section curves
  NCollection_Map<int>                      aMPassed;
  NCollection_DynamicArray<BOPDS_InterfFF>& aFFs  = myDS->InterfFF();
  int                                       aNbFF = aFFs.Length(), j;
  for (i = 0; i < aNbFF; ++i)
  {
    BOPDS_InterfFF& aFF = aFFs(i);
    // remove from Section pave blocks
    NCollection_DynamicArray<BOPDS_Curve>& aVNC = aFF.ChangeCurves();
    int                                    aNbC = aVNC.Length();
    for (j = 0; j < aNbC; ++j)
    {
      BOPDS_Curve&                                             aNC  = aVNC(j);
      NCollection_List<occ::handle<BOPDS_PaveBlock>>&          aLPB = aNC.ChangePaveBlocks();
      NCollection_List<occ::handle<BOPDS_PaveBlock>>::Iterator aItPB(aLPB);
      while (aItPB.More())
      {
        const occ::handle<BOPDS_PaveBlock>& aPB = aItPB.Value();
        if (theEdges.Contains(aPB->Edge()))
        {
          aLPB.Remove(aItPB);
        }
        else
        {
          aItPB.Next();
        }
      }
    }
  }

  // 3. From Face Info
  for (i = 0; i < myDS->NbSourceShapes(); ++i)
  {
    const BOPDS_ShapeInfo& aSI = myDS->ShapeInfo(i);
    if (aSI.ShapeType() != TopAbs_FACE)
    {
      continue;
    }
    if (!aSI.HasReference())
    {
      continue;
    }

    BOPDS_FaceInfo&                                       aFI     = myDS->ChangeFaceInfo(i);
    NCollection_IndexedMap<occ::handle<BOPDS_PaveBlock>>* aIMPB[] = {&aFI.ChangePaveBlocksIn(),
                                                                     &aFI.ChangePaveBlocksOn(),
                                                                     &aFI.ChangePaveBlocksSc()};
    for (int k = 0; k < 3; k++)
    {
      int aNbPB = aIMPB[k]->Extent(), m;
      for (m = 1; m <= aNbPB; ++m)
      {
        const occ::handle<BOPDS_PaveBlock>& aPB = aIMPB[k]->FindKey(m);
        if (theEdges.Contains(aPB->Edge()))
        {
          break;
        }
      }
      if (m <= aNbPB)
      {
        NCollection_IndexedMap<occ::handle<BOPDS_PaveBlock>> aMPBCopy = *aIMPB[k];
        aIMPB[k]->Clear();
        for (m = 1; m <= aNbPB; ++m)
        {
          const occ::handle<BOPDS_PaveBlock>& aPB = aMPBCopy(m);
          if (!theEdges.Contains(aPB->Edge()))
          {
            aIMPB[k]->Add(aPB);
          }
        }
      }
    }
  }
}

//=======================================================================
// function : ToleranceFF
// purpose  : Computes the TolFF according to the tolerance value and
//           types of the faces.
//=======================================================================
double ToleranceFF(const BRepAdaptor_Surface& aBAS1, const BRepAdaptor_Surface& aBAS2)
{
  double aTol1  = aBAS1.Tolerance();
  double aTol2  = aBAS2.Tolerance();
  double aTolFF = std::max(aTol1, aTol2);
  //
  bool isAna1, isAna2;
  isAna1 = (aBAS1.GetType() == GeomAbs_Plane || aBAS1.GetType() == GeomAbs_Cylinder
            || aBAS1.GetType() == GeomAbs_Cone || aBAS1.GetType() == GeomAbs_Sphere
            || aBAS1.GetType() == GeomAbs_Torus);
  //
  isAna2 = (aBAS2.GetType() == GeomAbs_Plane || aBAS2.GetType() == GeomAbs_Cylinder
            || aBAS2.GetType() == GeomAbs_Cone || aBAS2.GetType() == GeomAbs_Sphere
            || aBAS2.GetType() == GeomAbs_Torus);
  //
  if (!isAna1 || !isAna2)
  {
    aTolFF = std::max(aTolFF, 5.e-6);
  }
  return aTolFF;
}

//=================================================================================================

void BOPAlgo_PaveFiller::UpdateBlocksWithSharedVertices()
{
  if (!myNonDestructive)
  {
    return;
  }
  //
  int aNbFF;
  //
  NCollection_DynamicArray<BOPDS_InterfFF>& aFFs = myDS->InterfFF();
  aNbFF                                          = aFFs.Length();
  if (!aNbFF)
  {
    return;
  }
  //
  bool                 bOnCurve, bHasShapeSD;
  int                  i, nF1, nF2, aNbC, j, nV, nVSD;
  double               aTolV;
  NCollection_Map<int> aMF;
  // Old vertices UpdateVertex() has replaced by new ones
  NCollection_Map<int> aMVReplaced;
  //
  for (i = 0; i < aNbFF; ++i)
  {
    BOPDS_InterfFF& aFF = aFFs(i);
    //
    NCollection_DynamicArray<BOPDS_Curve>& aVC = aFF.ChangeCurves();
    aNbC                                       = aVC.Length();
    if (!aNbC)
    {
      continue;
    }
    //
    aFF.Indices(nF1, nF2);
    //
    if (aMF.Add(nF1))
    {
      myDS->UpdateFaceInfoOn(nF1);
    }
    if (aMF.Add(nF2))
    {
      myDS->UpdateFaceInfoOn(nF2);
    }
    //
    // Collect old vertices that are shared for nF1, nF2 ->aMI;
    NCollection_Map<int>           aMI;
    NCollection_Map<int>::Iterator aItMI;
    //
    BOPDS_FaceInfo& aFI1 = myDS->ChangeFaceInfo(nF1);
    BOPDS_FaceInfo& aFI2 = myDS->ChangeFaceInfo(nF2);
    //
    const NCollection_Map<int>& aMVOn1 = aFI1.VerticesOn();
    const NCollection_Map<int>& aMVIn1 = aFI1.VerticesIn();
    const NCollection_Map<int>& aMVOn2 = aFI2.VerticesOn();
    const NCollection_Map<int>& aMVIn2 = aFI2.VerticesIn();
    //
    for (j = 0; j < 2; ++j)
    {
      const NCollection_Map<int>& aMV1 = (!j) ? aMVOn1 : aMVIn1;
      aItMI.Initialize(aMV1);
      for (; aItMI.More(); aItMI.Next())
      {
        nV = aItMI.Value();
        if (myDS->IsNewShape(nV))
        {
          continue;
        }
        if (aMVOn2.Contains(nV) || aMVIn2.Contains(nV))
        {
          aMI.Add(nV);
        }
      }
    }
    //
    // Try to put vertices aMI on curves
    for (j = 0; j < aNbC; ++j)
    {
      BOPDS_Curve& aNC     = aVC.ChangeValue(j);
      double       aTolR3D = std::max(aNC.Tolerance(), aNC.TangentialTolerance());
      //
      aItMI.Initialize(aMI);
      for (; aItMI.More(); aItMI.Next())
      {
        nV = aItMI.Value();
        //
        bHasShapeSD = myDS->HasShapeSD(nV, nVSD);
        if (bHasShapeSD)
        {
          continue;
        }
        //
        bOnCurve = EstimatePaveOnCurve(nV, aNC, aTolR3D);
        if (!bOnCurve)
        {
          continue;
        }
        //
        const TopoDS_Vertex& aV = *((TopoDS_Vertex*)&myDS->Shape(nV));
        aTolV                   = BRep_Tool::Tolerance(aV);
        //
        if (UpdateVertex(nV, aTolV) != nV)
        {
          aMVReplaced.Add(nV);
        }
        myDS->InitPaveBlocksForVertex(nV);
      }
    } // for (j=0; j<aNbC; ++j) {
  } // for (i=0; i<aNbFF; ++i) {
  //
  UpdateCommonBlocksWithSDVertices();
  //
  // The faces' ON vertices were taken before these vertices were replaced,
  // and nothing takes them again before MakeBlocks() puts them on the
  // section curves: a section edge would end on the old vertex while the
  // split edges of the faces end on the new one, the section edge would
  // not connect, and the faces it should split would stay whole: a fuse of
  // two solids sharing a face came out empty. Whether a shared vertex is
  // replaced here can turn on 1e-14 in the input (tests/fork/occ-issues local02).
  if (!aMVReplaced.IsEmpty())
  {
    const int aNbS = myDS->NbSourceShapes();
    for (int nF = 0; nF < aNbS; ++nF)
    {
      if (myDS->ShapeInfo(nF).ShapeType() != TopAbs_FACE || !myDS->HasFaceInfo(nF))
      {
        continue;
      }
      NCollection_Map<int>& aMVOn = myDS->ChangeFaceInfo(nF).ChangeVerticesOn();
      for (NCollection_Map<int>::Iterator aItV(aMVReplaced); aItV.More(); aItV.Next())
      {
        nV = aItV.Value();
        if (aMVOn.Remove(nV))
        {
          aMVOn.Add(myDS->GetSameDomainIndex(nV));
        }
      }
    }
  }
}

//=================================================================================================

bool BOPAlgo_PaveFiller::EstimatePaveOnCurve(const int          nV,
                                             const BOPDS_Curve& aNC,
                                             const double       aTolR3D)
{
  bool   bIsVertexOnLine;
  double aT;
  //
  const TopoDS_Vertex&  aV  = *((TopoDS_Vertex*)&myDS->Shape(nV));
  const IntTools_Curve& aIC = aNC.Curve();
  //
  // The tolerance PutPaveOnCurve() puts the vertex on the curve with. With
  // less, a vertex at a distance between the two is not replaced here, yet
  // PutPaveOnCurve() takes it onto the curve and raises its tolerance in
  // place -- an input's vertex, in non-destructive mode (tests/fork/occ-issues
  // local04).
  bIsVertexOnLine = myContext->IsVertexOnLine(aV, aIC, aTolR3D + myFuzzyValue, aT);
  return bIsVertexOnLine;
}

//=================================================================================================

void BOPAlgo_PaveFiller::CorrectToleranceOfSE()
{
  NCollection_DynamicArray<BOPDS_InterfFF>& aFFs = myDS->InterfFF();
  NCollection_IndexedDataMap<int, NCollection_List<occ::handle<BOPDS_PaveBlock>>> aMVIPBs;
  NCollection_Map<int>                                                            aMVIToReduce;
  // Fence map to avoid repeated checking of the same edge
  NCollection_Map<occ::handle<BOPDS_PaveBlock>> aMPB;
  //
  // 1. iterate on all sections F-F
  int aNb = aFFs.Length(), i;
  for (i = 0; i < aNb; ++i)
  {
    BOPDS_InterfFF& aFF = aFFs(i);
    //
    NCollection_DynamicArray<BOPDS_Curve>& aVNC = aFF.ChangeCurves();
    int                                    aNbC = aVNC.Length(), k;
    for (k = 0; k < aNbC; ++k)
    {
      BOPDS_Curve&                                             aNC  = aVNC(k);
      NCollection_List<occ::handle<BOPDS_PaveBlock>>&          aLPB = aNC.ChangePaveBlocks();
      NCollection_List<occ::handle<BOPDS_PaveBlock>>::Iterator aItLPB(aLPB);
      for (; aItLPB.More();)
      {
        const occ::handle<BOPDS_PaveBlock>& aPB = aItLPB.Value();
        int                                 nE;
        if (!aPB->HasEdge(nE))
        {
          aLPB.Remove(aItLPB);
          continue;
        }
        //
        if (!aMPB.Add(aPB))
        {
          aItLPB.Next();
          continue;
        }
        //
        bool bIsReduced = false;
        if (aPB->OriginalEdge() < 0)
        {
          // It is possible that due to small angle between faces the
          // common zone between faces can be large and the tangential
          // tolerance of the curve will be large as well.
          // Here we're trying to reduce the tolerance of the section
          // edge using the valid tolerance of the edge.
          // Note, that if the pave block has created common block with
          // other edges its valid tolerance could have been changed to
          // cover all edges in common block (see PostTreatFF() method).
          double aTolC    = aNC.Tolerance();
          double aTolTang = aNC.TangentialTolerance();
          if (aTolC < aTolTang)
          {
            const TopoDS_Edge& aE    = TopoDS::Edge(myDS->Shape(nE));
            double             aTolE = BRep_Tool::Tolerance(aE);
            if (aTolC < aTolE)
            {
              // reduce edge tolerance
              static_cast<BRep_TEdge*>(aE.TShape().get())->Tolerance(aTolC);
              bIsReduced = true;
            }
          }
        }
        //
        // fill in the map vertex index - pave blocks
        for (int j = 0; j < 2; j++)
        {
          int nV = (j == 0 ? aPB->Pave1().Index() : aPB->Pave2().Index());
          myDS->HasShapeSD(nV, nV);
          NCollection_List<occ::handle<BOPDS_PaveBlock>>* pPBList = aMVIPBs.ChangeSeek(nV);
          if (!pPBList)
          {
            pPBList = &aMVIPBs.ChangeFromIndex(
              aMVIPBs.Add(nV, NCollection_List<occ::handle<BOPDS_PaveBlock>>()));
          }
          pPBList->Append(aPB);
          if (bIsReduced)
          {
            aMVIToReduce.Add(nV);
          }
        }
        aItLPB.Next();
      }
    }
  }
  //
  if (aMVIToReduce.IsEmpty())
  {
    return;
  }
  //
  // 2. try to reduce tolerances of connected vertices
  // 2.1 find all other edges containing these connected vertices to avoid
  //     reducing the tolerance to the value less than the tolerances of edges,
  //     i.e. minimal tolerance for the vertex is the max tolerance of the
  //     edges containing this vertex
  NCollection_DataMap<int, double>                                          aMVITol;
  NCollection_DynamicArray<NCollection_List<occ::handle<BOPDS_PaveBlock>>>& aPBP =
    myDS->ChangePaveBlocksPool();
  aNb = aPBP.Length();
  for (i = 0; i < aNb; ++i)
  {
    const NCollection_List<occ::handle<BOPDS_PaveBlock>>&    aLPB = aPBP(i);
    NCollection_List<occ::handle<BOPDS_PaveBlock>>::Iterator aItLPB(aLPB);
    for (; aItLPB.More(); aItLPB.Next())
    {
      const occ::handle<BOPDS_PaveBlock>& aPB = aItLPB.Value();
      int                                 nE;
      if (!aPB->HasEdge(nE))
      {
        continue;
      }
      const TopoDS_Edge& aE    = TopoDS::Edge(myDS->Shape(nE));
      double             aTolE = BRep_Tool::Tolerance(aE);
      //
      int nV[2];
      aPB->Indices(nV[0], nV[1]);
      //
      for (int j = 0; j < 2; j++)
      {
        if (aMVIToReduce.Contains(nV[j]))
        {
          double* aMaxTol = aMVITol.ChangeSeek(nV[j]);
          if (!aMaxTol)
          {
            aMVITol.Bind(nV[j], aTolE);
          }
          else if (aTolE > *aMaxTol)
          {
            *aMaxTol = aTolE;
          }
          NCollection_List<occ::handle<BOPDS_PaveBlock>>& aPBList = aMVIPBs.ChangeFromKey(nV[j]);
          aPBList.Append(aPB);
        }
      }
    }
  }
  //
  // 2.2 reduce tolerances if possible
  aNb = aMVIPBs.Extent();
  for (i = 1; i <= aNb; ++i)
  {
    int nV = aMVIPBs.FindKey(i);
    if (!aMVIToReduce.Contains(nV))
    {
      continue;
    }
    //
    const TopoDS_Vertex& aV      = TopoDS::Vertex(myDS->Shape(nV));
    double               aTolV   = BRep_Tool::Tolerance(aV);
    const double*        pMaxTol = aMVITol.Seek(nV);
    double               aMaxTol = pMaxTol ? *pMaxTol : 0.;
    // it makes no sense to compute the real tolerance if it is
    // impossible to reduce the tolerance at least 0.1% of the current value
    if (aTolV - aMaxTol < 0.001 * aTolV)
    {
      continue;
    }
    //
    // compute the maximal distance from the vertex to the adjacent edges
    gp_Pnt aP = BRep_Tool::Pnt(aV);
    //
    // Avoid repeated checks
    NCollection_Map<occ::handle<BOPDS_PaveBlock>> aMPBFence;
    //
    const NCollection_List<occ::handle<BOPDS_PaveBlock>>&    aLPB = aMVIPBs.FindFromIndex(i);
    NCollection_List<occ::handle<BOPDS_PaveBlock>>::Iterator aItLPB(aLPB);
    for (; aItLPB.More(); aItLPB.Next())
    {
      const occ::handle<BOPDS_PaveBlock>& aPB = aItLPB.Value();
      if (!aMPBFence.Add(aPB))
      {
        continue;
      }
      int                nE = aPB->Edge();
      const TopoDS_Edge& aE = TopoDS::Edge(myDS->Shape(nE));
      BRepAdaptor_Curve  aC(aE);
      for (int iPave = 0; iPave < 2; ++iPave)
      {
        const BOPDS_Pave& aPave = !iPave ? aPB->Pave1() : aPB->Pave2();
        int               nVSD  = aPave.Index();
        myDS->HasShapeSD(nVSD, nVSD);
        if (nVSD != nV)
        {
          continue;
        }
        //
        gp_Pnt aPonE = aC.Value(aPave.Parameter());
        double aDist = aP.Distance(aPonE);
        aDist += BRep_Tool::Tolerance(aE);
        if (aDist > aMaxTol)
        {
          aMaxTol = aDist;
        }
      }
    }
    //
    if (aMaxTol < aTolV)
    {
      static_cast<BRep_TVertex*>(aV.TShape().get())->Tolerance(aMaxTol);
    }
  }
}

//=================================================================================================

void BOPAlgo_PaveFiller::PutSEInOtherFaces(const Message_ProgressRange& theRange)
{
  // Try to intersect each section edge with the faces
  // not participated in its creation

  // Get all section edges
  NCollection_IndexedMap<occ::handle<BOPDS_PaveBlock>> aMPBScAll;

  NCollection_DynamicArray<BOPDS_InterfFF>& aFFs  = myDS->InterfFF();
  const int                                 aNbFF = aFFs.Length();
  Message_ProgressScope                     aPS(theRange, nullptr, 1);
  for (int i = 0; i < aNbFF; ++i)
  {
    const NCollection_DynamicArray<BOPDS_Curve>& aVNC = aFFs(i).Curves();
    const int                                    aNbC = aVNC.Length();
    for (int j = 0; j < aNbC; ++j)
    {
      const NCollection_List<occ::handle<BOPDS_PaveBlock>>&    aLPBC = aVNC(j).PaveBlocks();
      NCollection_List<occ::handle<BOPDS_PaveBlock>>::Iterator aItPB(aLPBC);
      for (; aItPB.More(); aItPB.Next())
      {
        aMPBScAll.Add(aItPB.Value());
      }
    }
  }
  // Perform intersection of collected pave blocks
  ForceInterfEF(aMPBScAll, aPS.Next(), false);
}

//=================================================================================================

void BOPAlgo_PaveFiller::RemoveMicroSectionEdges(
  NCollection_IndexedDataMap<TopoDS_Shape, BOPDS_CoupleOfPaveBlocks, TopTools_ShapeMapHasher>&
                                                        theMSCPB,
  NCollection_IndexedMap<occ::handle<BOPDS_PaveBlock>>& theMicroPB)
{
  if (theMSCPB.IsEmpty())
  {
    // no section edges
    return;
  }

  // Get all F/F interferences
  NCollection_DynamicArray<BOPDS_InterfFF>& aFFs = myDS->InterfFF();

  // Build the new map of section edges avoiding the micro edges
  NCollection_IndexedDataMap<TopoDS_Shape, BOPDS_CoupleOfPaveBlocks, TopTools_ShapeMapHasher>
    aSEPBMap;
  // Analyze all section edges
  int aNbCPB = theMSCPB.Extent();
  for (int i = 1; i <= aNbCPB; ++i)
  {
    const TopoDS_Shape&             aSI  = theMSCPB.FindKey(i);
    const BOPDS_CoupleOfPaveBlocks& aCPB = theMSCPB(i);

    if (aSI.ShapeType() != TopAbs_EDGE)
    {
      // Not an edge
      aSEPBMap.Add(aSI, aCPB);
      continue;
    }

    // Get pave block for analysis
    const occ::handle<BOPDS_PaveBlock>& aPB = aCPB.PaveBlock1();
    if (aPB->HasEdge())
    {
      // Not a real section edge
      aSEPBMap.Add(aSI, aCPB);
      continue;
    }

    if (!BOPTools_AlgoTools::IsMicroEdge(TopoDS::Edge(aSI), myContext, false))
    {
      // Normal edge
      aSEPBMap.Add(aSI, aCPB);
      continue;
    }

    // Micro edge is found, avoid it in the <theMSCPB> map
    // and remove from the F/F Intersection info structure

    // Get F/F interference which created this micro edge
    BOPDS_InterfFF& aFF = aFFs(aCPB.IndexInterf());
    // Get curve from which this edge has been created
    BOPDS_Curve& aCurve = aFF.ChangeCurves().ChangeValue(aCPB.Index());
    // Get all section pave blocks created from this curve
    NCollection_List<occ::handle<BOPDS_PaveBlock>>& aLPBC = aCurve.ChangePaveBlocks();
    // Remove pave block from the list
    for (NCollection_List<occ::handle<BOPDS_PaveBlock>>::Iterator it(aLPBC); it.More(); it.Next())
    {
      if (it.Value() == aPB)
      {
        aLPBC.Remove(it);
        break;
      }
    }

    // Add the pave block of "micro" edge into outgoing map for
    // unification of its vertices in the PostTreatFF method
    theMicroPB.Add(aPB);
  }

  // Overwrite the old map if necessary
  if (aSEPBMap.Extent() != theMSCPB.Extent())
  {
    theMSCPB = aSEPBMap;
  }
}

//=================================================================================================

void BOPAlgo_PaveFiller::RemoveMicroEdges()
{
  // Fence map
  NCollection_Map<occ::handle<BOPDS_PaveBlock>> aMPBFence;
  // Resulting map of micro edges
  NCollection_Map<int> aMicroEdges;
  // Check all pave blocks from the pool to find the micro edges
  NCollection_DynamicArray<NCollection_List<occ::handle<BOPDS_PaveBlock>>>& aPBP =
    myDS->ChangePaveBlocksPool();
  int aNbPBP = aPBP.Length();
  for (int i = 0; i < aNbPBP; ++i)
  {
    NCollection_List<occ::handle<BOPDS_PaveBlock>>& aLPB = aPBP(i);
    if (aLPB.Extent() < 2)
    {
      // No splits
      continue;
    }

    if (myDS->ShapeInfo(aLPB.First()->OriginalEdge()).HasFlag())
    {
      continue;
    }

    NCollection_List<occ::handle<BOPDS_PaveBlock>>::Iterator it(aLPB);
    for (; it.More(); it.Next())
    {
      const occ::handle<BOPDS_PaveBlock>& aPB  = it.Value();
      occ::handle<BOPDS_PaveBlock>        aPBR = myDS->RealPaveBlock(aPB);

      if (aMPBFence.Add(aPBR))
      {
        int nV1, nV2;
        aPBR->Indices(nV1, nV2);
        if (nV1 == nV2)
        {
          // Check if it has the valid range
          FillShrunkData(aPBR);
          if (!aPBR->HasShrunkData())
          {
            aMicroEdges.Add(aPBR->Edge());
          }
        }
      }
    }
  }
  RemovePaveBlocks(aMicroEdges);
}
