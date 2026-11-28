// Created on: 1995-06-22
// Created by: Flore Lantheaume
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

#include <BRepFilletAPI_MakeChamfer.hxx>
#include <ChFiDS_Spine.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopOpeBRepBuild_HBuilder.hxx>
#include <TopOpeBRepDS_HDataStructure.hxx>

#include <BOPAlgo_ArgumentAnalyzer.hxx>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <Bnd_Box.hxx>
#include <Extrema_ExtPC.hxx>
#include <NCollection_DataMap.hxx>
#include <NCollection_IndexedMap.hxx>
#include <Standard_Failure.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>

#include <vector>

namespace
{
using ShapeMap     = NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>;
using ShapeListMap =
  NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>;

//! The history of a result made of several computations, one contour after
//! another (BuildOneAfterAnother), kept where the builder keeps its builder
//! of topologic operations: BRepFilletAPI_MakeChamfer can take no member.
class ChamferChain : public TopOpeBRepBuild_HBuilder
{
public:
  ChamferChain(const TopOpeBRepDS_BuildTool& theBT)
      : TopOpeBRepBuild_HBuilder(theBT)
  {
  }

  //! the faces of the input that are not in the result as they were, with
  //! what they became (an empty list: nothing)
  ShapeListMap myModified;
  //! the edges and vertices of the input with the faces made from them
  ShapeListMap myGenerated;

  DEFINE_STANDARD_RTTI_INLINE(ChamferChain, TopOpeBRepBuild_HBuilder)
};

//! true inside a computation of BuildOneAfterAnother: no chain within a chain
thread_local bool THE_IN_CHAIN = false;

//! One contour as it was given: its edges and its size.
struct ChainContour
{
  NCollection_List<TopoDS_Shape> myEdges;
  bool                           myIsSym = true;
  double                         myDist  = 0.;
};

void AppendNew(NCollection_List<TopoDS_Shape>& theList, const TopoDS_Shape& theS)
{
  for (NCollection_List<TopoDS_Shape>::Iterator it(theList); it.More(); it.Next())
  {
    if (it.Value().IsSame(theS))
    {
      return;
    }
  }
  theList.Append(theS);
}

//! What is left in <theShape> of the edge <theE> of the input: the edge
//! itself, or the edges that lie along it -- their ends and middle on its
//! curve, within its range.
void EdgesLeftOf(const TopoDS_Edge&              theE,
                 const TopoDS_Shape&             theShape,
                 const ShapeMap&                 theSubShapes,
                 NCollection_List<TopoDS_Shape>& theLeft)
{
  if (theSubShapes.Contains(theE))
  {
    AppendNew(theLeft, theE);
    return;
  }
  BRepAdaptor_Curve aC(theE);
  Bnd_Box           aBox;
  BRepBndLib::Add(theE, aBox);
  const double aTol = std::max(10. * BRep_Tool::Tolerance(theE), 1.e-6);
  aBox.Enlarge(aTol);
  const auto IsOn = [&](const gp_Pnt& theP) {
    if (aBox.IsOut(theP))
    {
      return false;
    }
    Extrema_ExtPC anExt(theP, aC);
    if (anExt.IsDone())
    {
      for (int k = 1; k <= anExt.NbExt(); k++)
      {
        if (anExt.SquareDistance(k) < aTol * aTol)
        {
          return true;
        }
      }
    }
    return theP.Distance(aC.Value(aC.FirstParameter())) < aTol
           || theP.Distance(aC.Value(aC.LastParameter())) < aTol;
  };
  ShapeMap aSeen;
  for (TopExp_Explorer ex(theShape, TopAbs_EDGE); ex.More(); ex.Next())
  {
    const TopoDS_Edge& anX = TopoDS::Edge(ex.Current());
    if (!aSeen.Add(anX) || BRep_Tool::Degenerated(anX))
    {
      continue;
    }
    BRepAdaptor_Curve aCX(anX);
    const double      aF = aCX.FirstParameter(), aL = aCX.LastParameter();
    if (IsOn(aCX.Value(0.5 * (aF + aL))) && IsOn(aCX.Value(aF)) && IsOn(aCX.Value(aL)))
    {
      AppendNew(theLeft, anX);
    }
  }
}

//! Where a chain is: the shape so far, the faces of the input with what
//! they are in it (none bound: itself), and the edges and vertices of the
//! input with the faces made from them.
struct ChainState
{
  TopoDS_Shape myShape;
  ShapeListMap myModified;
  ShapeListMap myGenerated;
};

//! What a chain is made on.
struct ChainInput
{
  TopoDS_Shape                                                  myShape;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> myFaces;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> mySources;
};

double VolumeOf(const TopoDS_Shape& theS)
{
  GProp_GProps aProps;
  BRepGProp::VolumeProperties(theS, aProps);
  return aProps.Mass();
}

//! The loosest edge or vertex of <theS>.
double MaxToleranceOf(const TopoDS_Shape& theS)
{
  double aMax = 0.;
  for (TopExp_Explorer ex(theS, TopAbs_EDGE); ex.More(); ex.Next())
  {
    aMax = std::max(aMax, BRep_Tool::Tolerance(TopoDS::Edge(ex.Current())));
  }
  for (TopExp_Explorer ex(theS, TopAbs_VERTEX); ex.More(); ex.Next())
  {
    aMax = std::max(aMax, BRep_Tool::Tolerance(TopoDS::Vertex(ex.Current())));
  }
  return aMax;
}

//! One contour made on what <theState> left of its edges. False where
//! nothing is left of them, the computation fails, or an edge left is not
//! chamfered by it.
bool ChainStep(const ChainInput& theInput, const ChainContour& theContour, ChainState& theState)
{
  const TopoDS_Shape& aCur = theState.myShape;
  ShapeMap            aSub;
  for (TopExp_Explorer ex(aCur, TopAbs_FACE); ex.More(); ex.Next())
  {
    aSub.Add(ex.Current());
  }
  for (TopExp_Explorer ex(aCur, TopAbs_EDGE); ex.More(); ex.Next())
  {
    aSub.Add(ex.Current());
  }
  for (TopExp_Explorer ex(aCur, TopAbs_VERTEX); ex.More(); ex.Next())
  {
    aSub.Add(ex.Current());
  }
  // what is left of the contour's edges
  ShapeListMap                   aLeftOf;
  NCollection_List<TopoDS_Shape> aLeft;
  for (NCollection_List<TopoDS_Shape>::Iterator itE(theContour.myEdges); itE.More(); itE.Next())
  {
    NCollection_List<TopoDS_Shape> aL;
    EdgesLeftOf(TopoDS::Edge(itE.Value()), aCur, aSub, aL);
    for (NCollection_List<TopoDS_Shape>::Iterator it(aL); it.More(); it.Next())
    {
      AppendNew(aLeft, it.Value());
    }
    aLeftOf.Bind(itE.Value(), aL);
  }
  if (aLeft.IsEmpty())
  {
    return false;
  }
  BRepFilletAPI_MakeChamfer aMk(aCur);
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    anEF;
  TopExp::MapShapesAndAncestors(aCur, TopAbs_EDGE, TopAbs_FACE, anEF);
  for (NCollection_List<TopoDS_Shape>::Iterator it(aLeft); it.More(); it.Next())
  {
    const TopoDS_Edge& anE = TopoDS::Edge(it.Value());
    if (theContour.myIsSym)
    {
      aMk.Add(theContour.myDist, anE);
    }
    else
    {
      if (!anEF.Contains(anE) || anEF.FindFromKey(anE).IsEmpty())
      {
        return false;
      }
      aMk.Add(theContour.myDist,
              theContour.myDist,
              anE,
              TopoDS::Face(anEF.FindFromKey(anE).First()));
    }
  }
  if (aMk.NbContours() == 0)
  {
    return false;
  }
  try
  {
    aMk.Build();
  }
  catch (Standard_Failure const&)
  {
    return false;
  }
  if (!aMk.IsDone())
  {
    return false;
  }
  for (NCollection_List<TopoDS_Shape>::Iterator it(aLeft); it.More(); it.Next())
  {
    if (aMk.Generated(it.Value()).IsEmpty())
    {
      return false;
    }
  }
  const TopoDS_Shape aNext = aMk.Shape();
  ShapeMap           aNextFaces;
  for (TopExp_Explorer ex(aNext, TopAbs_FACE); ex.More(); ex.Next())
  {
    aNextFaces.Add(ex.Current());
  }
  // a face of aCur in aNext
  const auto Image = [&](const TopoDS_Shape& theF, NCollection_List<TopoDS_Shape>& theInto) {
    const NCollection_List<TopoDS_Shape>& aM = aMk.Modified(theF);
    if (!aM.IsEmpty())
    {
      for (NCollection_List<TopoDS_Shape>::Iterator it(aM); it.More(); it.Next())
      {
        AppendNew(theInto, it.Value());
      }
    }
    else if (aNextFaces.Contains(theF))
    {
      AppendNew(theInto, theF);
    }
  };
  ShapeListMap& aMod = theState.myModified;
  ShapeListMap& aGen = theState.myGenerated;
  for (int i = 1; i <= theInput.myFaces.Extent(); i++)
  {
    const TopoDS_Shape&            aF = theInput.myFaces(i);
    NCollection_List<TopoDS_Shape> aNew;
    if (aMod.IsBound(aF))
    {
      for (NCollection_List<TopoDS_Shape>::Iterator it(aMod.Find(aF)); it.More(); it.Next())
      {
        Image(it.Value(), aNew);
      }
    }
    else if (aSub.Contains(aF))
    {
      Image(aF, aNew);
    }
    if (aNew.Extent() != 1 || !aNew.First().IsSame(aF))
    {
      aMod.Bind(aF, aNew);
    }
    else
    {
      aMod.UnBind(aF);
    }
  }
  for (int i = 1; i <= theInput.mySources.Extent(); i++)
  {
    const TopoDS_Shape&            aS = theInput.mySources(i);
    NCollection_List<TopoDS_Shape> aNew;
    if (aGen.IsBound(aS))
    {
      for (NCollection_List<TopoDS_Shape>::Iterator it(aGen.Find(aS)); it.More(); it.Next())
      {
        Image(it.Value(), aNew);
      }
    }
    NCollection_List<TopoDS_Shape> aNow;
    if (aLeftOf.IsBound(aS))
    {
      aNow = aLeftOf.Find(aS);
    }
    else if (aSub.Contains(aS))
    {
      aNow.Append(aS);
    }
    for (NCollection_List<TopoDS_Shape>::Iterator itN(aNow); itN.More(); itN.Next())
    {
      const NCollection_List<TopoDS_Shape>& aG = aMk.Generated(itN.Value());
      for (NCollection_List<TopoDS_Shape>::Iterator it(aG); it.More(); it.Next())
      {
        AppendNew(aNew, it.Value());
      }
    }
    if (!aNew.IsEmpty() || aGen.IsBound(aS))
    {
      aGen.Bind(aS, aNew);
    }
  }
  theState.myShape = aNext;
  return true;
}

//! Where the contours of <theBuilder> were not made together: one after
//! another, each on what the ones before it left of its edges -- in the
//! order given and, where that makes nothing, in the reverse. What two
//! chamfer features one after the other make; at a vertex of four sharp
//! edges and more the builder has no corner for two of them (#962's tooth,
//! #631's wedge). At the latest fix level, whatever the caller's: what a
//! level refuses goes on to the levels above it.
//!
//! Only where each contour is made by itself, and the result's volume is
//! the input's with what each takes or adds by itself to within a quarter
//! of those amounts: each ends on the others where they meet, a little less
//! than by itself (#962's, 1.07 of 13.0; #273's, 0.66 of 22.6) -- one that
//! took 139 where the two by themselves take 6.5 (#474's coil at its block)
//! is no chamfer of them; nor, where all take or all add, one that takes or
//! adds more than they do by themselves (#962's arm, an edge of 0.05
//! chamfered at 1 with its neighbour: 1.88 for 1.55). And a result BOP's check finds crossing itself is
//! none, nor one with an edge looser than the input's or a twentieth of the
//! smallest chamfer (the terms the builder's own fallbacks keep a result
//! on; #273's handle, a knot of facets with a vertex of 0.094 at size 1).
//!
//! Two contours only. Three and more at a vertex (the sweeps' every edge of
//! a vertex of four) come out valid and within these terms with slits and
//! undercuts where the third meets what the first two left: no chamfer
//! anyone would take.
bool BuildOneAfterAnother(ChFi3d_ChBuilder& theBuilder,
                          TopoDS_Shape&     theResult,
                          ShapeListMap&     theModified,
                          ShapeListMap&     theGenerated)
{
  const int aNb = theBuilder.NbElements();
  if (THE_IN_CHAIN || aNb != 2)
  {
    return false;
  }
  std::vector<ChainContour> aContours(aNb);
  for (int anIC = 1; anIC <= aNb; anIC++)
  {
    ChainContour& aContour = aContours[anIC - 1];
    if (theBuilder.IsChamfer(anIC) == ChFiDS_Sym)
    {
      theBuilder.GetDist(anIC, aContour.myDist);
    }
    else if (theBuilder.IsChamfer(anIC) == ChFiDS_TwoDist)
    {
      double aD1 = 0., aD2 = 0.;
      theBuilder.Dists(anIC, aD1, aD2);
      if (std::abs(aD1 - aD2) > 1.e-9 * std::max(aD1, aD2))
      {
        return false;
      }
      aContour.myDist  = aD1;
      aContour.myIsSym = false;
    }
    else
    {
      return false;
    }
    const occ::handle<ChFiDS_Spine> aSpine = theBuilder.Value(anIC);
    if (aSpine.IsNull() || aContour.myDist <= 0.)
    {
      return false;
    }
    for (int j = 1; j <= aSpine->NbEdges(); j++)
    {
      aContour.myEdges.Append(aSpine->Edges(j));
    }
  }
  struct ChainGuard
  {
    ChainGuard()
        : myLevel(ChFi3d_SetFixLevel(ChFi3d_FixLevelLatest()))
    {
      THE_IN_CHAIN = true;
    }

    ~ChainGuard()
    {
      THE_IN_CHAIN = false;
      ChFi3d_SetFixLevel(myLevel);
    }

    int myLevel;
  } aGuard;
  ChainInput anInput;
  anInput.myShape = theBuilder.InputShape();
  TopExp::MapShapes(anInput.myShape, TopAbs_FACE, anInput.myFaces);
  TopExp::MapShapes(anInput.myShape, TopAbs_EDGE, anInput.mySources);
  TopExp::MapShapes(anInput.myShape, TopAbs_VERTEX, anInput.mySources);
  // each by itself
  std::vector<ChainState> anAlone(aNb);
  for (int i = 0; i < aNb; i++)
  {
    anAlone[i].myShape = anInput.myShape;
    if (!ChainStep(anInput, aContours[i], anAlone[i]))
    {
      return false;
    }
  }
  const double aVolume    = VolumeOf(anInput.myShape);
  const double aTolerance = MaxToleranceOf(anInput.myShape);
  double       aSmallest  = aContours[0].myDist;
  for (int i = 1; i < aNb; i++)
  {
    aSmallest = std::min(aSmallest, aContours[i].myDist);
  }
  double       aSum = 0., anAbs = 0.;
  for (int i = 0; i < aNb; i++)
  {
    const double aD = VolumeOf(anAlone[i].myShape) - aVolume;
    aSum  += aD;
    anAbs += std::abs(aD);
  }
  for (int aTry = 0; aTry < 2; aTry++)
  {
    ChainState aState = anAlone[aTry == 0 ? 0 : aNb - 1];
    bool       isMade = true;
    for (int k = 1; k < aNb && isMade; k++)
    {
      isMade = ChainStep(anInput, aContours[aTry == 0 ? k : aNb - 1 - k], aState);
    }
    if (!isMade)
    {
      continue;
    }
    const double aChange = VolumeOf(aState.myShape) - aVolume;
    if (std::abs(aChange - aSum) > 0.25 * anAbs
        || (std::abs(aSum) > anAbs * (1. - 1.e-9) && std::abs(aChange) > anAbs * (1. + 1.e-6)))
    {
      // (the second: all taking or all adding, and more of it than they
      // do by themselves)
      continue;
    }
    if (MaxToleranceOf(aState.myShape) > std::max(aTolerance, 0.05 * aSmallest))
    {
      continue;
    }
    BOPAlgo_ArgumentAnalyzer aCheck;
    aCheck.SetShape1(aState.myShape);
    aCheck.ArgumentTypeMode()   = false;
    aCheck.SelfInterMode()      = true;
    aCheck.SmallEdgeMode()      = false;
    aCheck.RebuildFaceMode()    = false;
    aCheck.TangentMode()        = false;
    aCheck.MergeVertexMode()    = false;
    aCheck.MergeEdgeMode()      = false;
    aCheck.ContinuityMode()     = false;
    aCheck.CurveOnSurfaceMode() = false;
    aCheck.StopOnFirstFaulty()  = true;
    aCheck.Perform();
    if (aCheck.HasFaulty())
    {
      continue;
    }
    theResult    = aState.myShape;
    theModified  = aState.myModified;
    theGenerated = aState.myGenerated;
    return true;
  }
  return false;
}
} // namespace

//=================================================================================================

BRepFilletAPI_MakeChamfer::BRepFilletAPI_MakeChamfer(const TopoDS_Shape& S)
    : myBuilder(S)
{
}

//=================================================================================================

void BRepFilletAPI_MakeChamfer::Add(const TopoDS_Edge& E)
{
  myBuilder.Add(E);
}

//=================================================================================================

void BRepFilletAPI_MakeChamfer::Add(const double Dis, const TopoDS_Edge& E)
{
  myBuilder.Add(Dis, E);
}

//=================================================================================================

void BRepFilletAPI_MakeChamfer::SetDist(const double Dis, const int IC, const TopoDS_Face& F)
{
  myBuilder.SetDist(Dis, IC, F);
}

//=================================================================================================

void BRepFilletAPI_MakeChamfer::GetDist(const int IC, double& Dis) const
{
  myBuilder.GetDist(IC, Dis);
}

//=================================================================================================

void BRepFilletAPI_MakeChamfer::Add(const double       Dis1,
                                    const double       Dis2,
                                    const TopoDS_Edge& E,
                                    const TopoDS_Face& F)
{
  myBuilder.Add(Dis1, Dis2, E, F);
}

//=================================================================================================

void BRepFilletAPI_MakeChamfer::SetDists(const double       Dis1,
                                         const double       Dis2,
                                         const int          IC,
                                         const TopoDS_Face& F)
{
  myBuilder.SetDists(Dis1, Dis2, IC, F);
}

//=================================================================================================

void BRepFilletAPI_MakeChamfer::Dists(const int IC, double& Dis1, double& Dis2) const
{
  double temp1, temp2;
  myBuilder.Dists(IC, temp1, temp2);
  Dis1 = temp1;
  Dis2 = temp2;
}

//=================================================================================================

void BRepFilletAPI_MakeChamfer::AddDA(const double       Dis,
                                      const double       Angle,
                                      const TopoDS_Edge& E,
                                      const TopoDS_Face& F)
{
  myBuilder.AddDA(Dis, Angle, E, F);
}

//=================================================================================================

void BRepFilletAPI_MakeChamfer::SetDistAngle(const double       Dis,
                                             const double       Angle,
                                             const int          IC,
                                             const TopoDS_Face& F)
{
  myBuilder.SetDistAngle(Dis, Angle, IC, F);
}

//=================================================================================================

void BRepFilletAPI_MakeChamfer::GetDistAngle(const int IC, double& Dis, double& Angle) const
{
  myBuilder.GetDistAngle(IC, Dis, Angle);
}

//=================================================================================================

void BRepFilletAPI_MakeChamfer::SetMode(const ChFiDS_ChamfMode theMode)
{
  myBuilder.SetMode(theMode);
}

//=================================================================================================

bool BRepFilletAPI_MakeChamfer::IsSymetric(const int IC) const
{
  ChFiDS_ChamfMethod ChamfMeth = myBuilder.IsChamfer(IC);
  bool               ret       = false;

  if (ChamfMeth == ChFiDS_Sym)
  {
    ret = true;
  }

  return ret;
}

//=================================================================================================

bool BRepFilletAPI_MakeChamfer::IsTwoDistances(const int IC) const
{
  ChFiDS_ChamfMethod ChamfMeth = myBuilder.IsChamfer(IC);
  bool               ret       = false;

  if (ChamfMeth == ChFiDS_TwoDist)
  {
    ret = true;
  }

  return ret;
}

//=================================================================================================

bool BRepFilletAPI_MakeChamfer::IsDistanceAngle(const int IC) const
{
  ChFiDS_ChamfMethod ChamfMeth = myBuilder.IsChamfer(IC);
  bool               ret       = false;

  if (ChamfMeth == ChFiDS_DistAngle)
  {
    ret = true;
  }

  return ret;
}

//=================================================================================================

void BRepFilletAPI_MakeChamfer::ResetContour(const int IC)
{
  myBuilder.ResetContour(IC);
}

//=================================================================================================

int BRepFilletAPI_MakeChamfer::NbContours() const
{
  return myBuilder.NbElements();
}

//=================================================================================================

int BRepFilletAPI_MakeChamfer::Contour(const TopoDS_Edge& E) const
{
  return myBuilder.Contains(E);
}

//=================================================================================================

int BRepFilletAPI_MakeChamfer::NbEdges(const int I) const

{
  const occ::handle<ChFiDS_Spine>& Spine = myBuilder.Value(I);
  return (Spine->NbEdges());
}

//=================================================================================================

const TopoDS_Edge& BRepFilletAPI_MakeChamfer::Edge(const int I, const int J) const
{
  const occ::handle<ChFiDS_Spine>& Spine = myBuilder.Value(I);
  const TopoDS_Edge&               E     = Spine->Edges(J);
  return E;
}

//=================================================================================================

void BRepFilletAPI_MakeChamfer::Remove(const TopoDS_Edge& E)
{
  myBuilder.Remove(E);
}

//=================================================================================================

double BRepFilletAPI_MakeChamfer::Length(const int IC) const
{
  return myBuilder.Length(IC);
}

//=================================================================================================

TopoDS_Vertex BRepFilletAPI_MakeChamfer::FirstVertex(const int IC) const
{
  return myBuilder.FirstVertex(IC);
}

//=================================================================================================

TopoDS_Vertex BRepFilletAPI_MakeChamfer::LastVertex(const int IC) const
{
  return myBuilder.LastVertex(IC);
}

//=================================================================================================

double BRepFilletAPI_MakeChamfer::Abscissa(const int IC, const TopoDS_Vertex& V) const
{
  return myBuilder.Abscissa(IC, V);
}

//=================================================================================================

double BRepFilletAPI_MakeChamfer::RelativeAbscissa(const int IC, const TopoDS_Vertex& V) const
{
  return myBuilder.RelativeAbscissa(IC, V);
}

//=================================================================================================

bool BRepFilletAPI_MakeChamfer::ClosedAndTangent(const int IC) const
{
  return myBuilder.ClosedAndTangent(IC);
}

//=================================================================================================

bool BRepFilletAPI_MakeChamfer::Closed(const int IC) const
{
  return myBuilder.Closed(IC);
}

//=================================================================================================

occ::handle<TopOpeBRepBuild_HBuilder> BRepFilletAPI_MakeChamfer::Builder() const
{
  return myBuilder.Builder();
}

//=================================================================================================

void BRepFilletAPI_MakeChamfer::Build(const Message_ProgressRange& /*theRange*/)
{
  if (!occ::down_cast<ChamferChain>(myBuilder.Builder()).IsNull())
  {
    // a chain's history from a computation before
    myBuilder.SetBuilder(new TopOpeBRepBuild_HBuilder(myBuilder.Builder()->BuildTool()));
  }
  std::exception_ptr aFailure;
  try
  {
    myBuilder.Compute();
  }
  catch (Standard_Failure const&)
  {
    if (THE_IN_CHAIN || myBuilder.NbElements() != 2)
    {
      throw;
    }
    aFailure = std::current_exception();
  }
  bool isMade = !aFailure && myBuilder.IsDone();
  if (isMade)
  {
    myShape = myBuilder.Shape();
  }
  else
  {
    // the contours one after another (BuildOneAfterAnother)
    occ::handle<ChamferChain> aChain = new ChamferChain(myBuilder.Builder()->BuildTool());
    TopoDS_Shape              aResult;
    bool                      isChained = false;
    try
    {
      isChained = BuildOneAfterAnother(myBuilder, aResult, aChain->myModified, aChain->myGenerated);
    }
    catch (Standard_Failure const&)
    {
      isChained = false;
    }
    if (isChained)
    {
      myBuilder.SetBuilder(aChain);
      myShape = aResult;
      isMade  = true;
    }
    else if (aFailure)
    {
      std::rethrow_exception(aFailure);
    }
  }
  if (isMade)
  {
    Done();

    // creation of the Map; a thawed copy stands for the face it was copied
    // from (BRepFilletAPI_MakeFillet::Build)
    TopExp_Explorer ex;
    for (ex.Init(myShape, TopAbs_FACE); ex.More(); ex.Next())
    {
      myMap.Add(ex.Current());
      TopoDS_Shape anOriginal = ex.Current();
      for (occ::handle<TopoDS_TShape> aFrom =
             TopoDS_TShape::ThawedFrom(ex.Current().TShape().get());
           !aFrom.IsNull();
           aFrom = TopoDS_TShape::ThawedFrom(aFrom.get()))
      {
        anOriginal.TShape(aFrom);
        myMap.Add(anOriginal);
      }
    }
  }
}

//=================================================================================================

void BRepFilletAPI_MakeChamfer::Reset()
{
  NotDone();
  myBuilder.Reset();
  myMap.Clear();
}

//=================================================================================================

const NCollection_List<TopoDS_Shape>& BRepFilletAPI_MakeChamfer::Generated(const TopoDS_Shape& EorV)
{
  const occ::handle<ChamferChain> aChain = occ::down_cast<ChamferChain>(myBuilder.Builder());
  if (!aChain.IsNull())
  {
    myGenerated.Clear();
    if (aChain->myGenerated.IsBound(EorV))
    {
      myGenerated = aChain->myGenerated.Find(EorV);
    }
    return myGenerated;
  }
  return myBuilder.Generated(EorV);
}

//=================================================================================================

const NCollection_List<TopoDS_Shape>& BRepFilletAPI_MakeChamfer::Modified(const TopoDS_Shape& F)
{
  myGenerated.Clear();

  const occ::handle<ChamferChain> aChain = occ::down_cast<ChamferChain>(myBuilder.Builder());
  if (!aChain.IsNull())
  {
    if (aChain->myModified.IsBound(F))
    {
      myGenerated = aChain->myModified.Find(F);
    }
    return myGenerated;
  }

  if (myBuilder.Builder()->IsSplit(F, TopAbs_OUT))
  {
    NCollection_List<TopoDS_Shape>::Iterator It(myBuilder.Builder()->Splits(F, TopAbs_OUT));
    for (; It.More(); It.Next())
    {
      myGenerated.Append(It.Value());
    }
  }
  if (myBuilder.Builder()->IsSplit(F, TopAbs_IN))
  {
    NCollection_List<TopoDS_Shape>::Iterator It(myBuilder.Builder()->Splits(F, TopAbs_IN));
    for (; It.More(); It.Next())
    {
      myGenerated.Append(It.Value());
    }
  }
  if (myBuilder.Builder()->IsSplit(F, TopAbs_ON))
  {
    NCollection_List<TopoDS_Shape>::Iterator It(myBuilder.Builder()->Splits(F, TopAbs_ON));
    for (; It.More(); It.Next())
    {
      myGenerated.Append(It.Value());
    }
  }
  return myGenerated;
}

//=================================================================================================

bool BRepFilletAPI_MakeChamfer::IsDeleted(const TopoDS_Shape& F)
{
  const occ::handle<ChamferChain> aChain = occ::down_cast<ChamferChain>(myBuilder.Builder());
  if (!aChain.IsNull())
  {
    return !(myMap.Contains(F)
             || (aChain->myModified.IsBound(F) && !aChain->myModified.Find(F).IsEmpty()));
  }
  return !(myMap.Contains(F) || myBuilder.Builder()->IsSplit(F, TopAbs_OUT)
           || myBuilder.Builder()->IsSplit(F, TopAbs_IN)
           || myBuilder.Builder()->IsSplit(F, TopAbs_ON));
}

//=================================================================================================

void BRepFilletAPI_MakeChamfer::Simulate(const int IC)
{
  myBuilder.Simulate(IC);
}

//=================================================================================================

int BRepFilletAPI_MakeChamfer::NbSurf(const int IC) const
{
  return myBuilder.NbSurf(IC);
}

//=================================================================================================

occ::handle<NCollection_HArray1<ChFiDS_CircSection>> BRepFilletAPI_MakeChamfer::Sect(
  const int IC,
  const int IS) const
{
  return myBuilder.Sect(IC, IS);
}
