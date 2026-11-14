// Created on: 1996-07-26
// Created by: Stagiaire Xuan Trang PHAMPHU
// Copyright (c) 1996-1999 Matra Datavision
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

#include <ChFi3d_ChBuilder.hxx>

#include <ChFi3d_Builder_0.hxx>
#include <ChFiDS_ChamfSpine.hxx>
#include <ChFiDS_Spine.hxx>
#include <ChFiDS_Stripe.hxx>
#include <ChFiDS_SurfData.hxx>
#include <NCollection_Map.hxx>
#include <TopTools_ShapeMapHasher.hxx>

#include <algorithm>
#include <cmath>

//=================================================================================================

// The larger of a chamfer's two distances from its edge (0 where the spine
// is no chamfer's).
static double ChamferSize(const occ::handle<ChFiDS_Spine>& theSpine)
{
  const occ::handle<ChFiDS_ChamfSpine> aSpine = occ::down_cast<ChFiDS_ChamfSpine>(theSpine);
  if (aSpine.IsNull())
  {
    return 0.;
  }
  double aD1 = 0., aD2 = 0.;
  switch (aSpine->IsChamfer())
  {
    case ChFiDS_Sym:
      aSpine->GetDist(aD1);
      aD2 = aD1;
      break;
    case ChFiDS_TwoDist:
      aSpine->Dists(aD1, aD2);
      break;
    case ChFiDS_DistAngle: {
      double anAngle = 0.;
      aSpine->GetDistAngle(aD1, anAngle);
      aD2 = aD1 * std::tan(std::min(std::abs(anAngle), 1.4));
      break;
    }
  }
  return std::max(aD1, aD2);
}

//=================================================================================================

void ChFi3d_ChBuilder::PerformTwoCorner(const int Index)
{
  if (PerformTwoCornerbyInter(Index))
  {
    return;
  }
  // The two chamfers' corner cannot be made by their intersection -- one
  // edge convex and the other concave, their lines on the common face not
  // crossing, an end not on an edge. The intersection gives up before it
  // has put anything in the data structure; its failure was not looked at,
  // and the corner was left with nothing. The corner of several edges
  // makes it, as it does for two fillets
  // (ChFi3d_FilBuilder::PerformTwoCorner) -- where the vertex has the two
  // edges and one more, and neither chamfer's surface runs on past the
  // vertex by more than one and a half times its size: that corner takes each surface's
  // end as it is, and one carried on along the next edge stayed in the
  // result, a chamfer of an edge nobody asked for (#962's arm, 13 along
  // the block's foot; valid, and clean by BOP's check). And the result is
  // refused where that corner left an edge looser than a twentieth of the
  // chamfer (ChFi3d_TwoChamferCorner): a chamfer of 1 on a piece of an
  // edge 0.7 long came out with vertices of 0.23.
  const TopoDS_Vertex& aVtx = myVDataMap.FindKey(Index);
  // (the map lists an edge once for each face it is in)
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> anEdges;
  if (myVEMap.Contains(aVtx))
  {
    for (NCollection_List<TopoDS_Shape>::Iterator itE(myVEMap(aVtx)); itE.More(); itE.Next())
    {
      anEdges.Add(itE.Value());
    }
  }
  bool isAtVertex = anEdges.Extent() == 3;
  for (NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator itS(myVDataMap(Index));
       itS.More() && isAtVertex;
       itS.Next())
  {
    const occ::handle<ChFiDS_Stripe>& aStripe = itS.Value();
    const occ::handle<ChFiDS_Spine>&  aSpine  = aStripe->Spine();
    int                               aSens   = 0;
    const int                         anIndex = ChFi3d_IndexOfSurfData(aVtx, aStripe, aSens);
    if (anIndex < 1 || aStripe->SetOfSurfData().IsNull()
        || anIndex > aStripe->SetOfSurfData()->Length())
    {
      isAtVertex = false;
      break;
    }
    const occ::handle<ChFiDS_SurfData>& aSD     = aStripe->SetOfSurfData()->Value(anIndex);
    const bool                          isFirst = aSens == 1;
    const double aPast  = isFirst ? aSpine->FirstParameter(1) - aSD->FirstSpineParam()
                                  : aSD->LastSpineParam() - aSpine->LastParameter(aSpine->NbEdges());
    isAtVertex = aPast <= 1.5 * ChamferSize(aSpine);
  }
  if (isAtVertex)
  {
    // the result is held to this corner's terms (ChFi3d_Builder::Compute)
    for (NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator itS(myVDataMap(Index)); itS.More();
         itS.Next())
    {
      const double aSize = ChamferSize(itS.Value()->Spine());
      double&      aMin  = ChFi3d_TwoChamferCorner();
      if (aSize > 0. && (aMin <= 0. || aSize < aMin))
      {
        aMin = aSize;
      }
    }
    PerformMoreThreeCorner(Index, 2);
    done = true;
  }
}
