// Created on: 1994-03-29
// Created by: Isabelle GRIGNON
// Copyright (c) 1994-1999 Matra Datavision
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

#include <Adaptor3d_CurveOnSurface.hxx>
#include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBlend_Line.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepLProp_CLProps.hxx>
#include <BRepTopAdaptor_TopolTool.hxx>
#include <BSplCLib.hxx>
#include <ChFi3d.hxx>
#include <ChFi3d_Builder_0.hxx>
#include <ChFi3d_FilBuilder.hxx>
#include <ChFiDS_CommonPoint.hxx>
#include <ChFiDS_FaceInterference.hxx>
#include <ChFiDS_FilSpine.hxx>
#include <ChFiDS_SurfData.hxx>
#include <NCollection_Sequence.hxx>
#include <NCollection_HSequence.hxx>
#include <ChFiDS_Regul.hxx>
#include <ChFiDS_Spine.hxx>
#include <ChFiDS_State.hxx>
#include <ChFiDS_Stripe.hxx>
#include <ChFiKPart_ComputeData.hxx>
#include <ElCLib.hxx>
#include <Extrema_ExtCC.hxx>
#include <Extrema_ExtPC.hxx>
#include <Extrema_POnCurv.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_ElementarySurface.hxx>
#include <Geom_Line.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <Geom2dConvert.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <GeomFill_ConstrainedFilling.hxx>
#include <GeomFill_SimpleBound.hxx>
#include <GeomProjLib.hxx>
#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Lin.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <Standard_ErrorHandler.hxx>
#include <Standard_NotImplemented.hxx>
#include <Precision.hxx>
#include <StdFail_NotDone.hxx>
#include <NCollection_Array1.hxx>
#include <Standard_Integer.hxx>
#include <NCollection_List.hxx>
#include <TopAbs.hxx>
#include <TopAbs_Orientation.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopOpeBRepDS_DataStructure.hxx>
#include <TopOpeBRepDS_HDataStructure.hxx>

#include <exception>

#ifdef OCCT_DEBUG
  #include <Geom_TrimmedCurve.hxx>
extern bool ChFi3d_GettraceDRAWSPINE();
extern bool ChFi3d_GetcontextFORCEFILLING();
  #include <OSD_Chronometer.hxx>

extern double t_t2cornerinit, t_perf2cornerbyinter, t_chfikpartcompdata, t_cheminement,
  t_remplissage, t_t2cornerDS;
extern void ChFi3d_InitChron(OSD_Chronometer& ch);
extern void ChFi3d_ResultChron(OSD_Chronometer& ch, double& time);
#endif

//=======================================================================
// function : ToricRotule
// purpose  : Test if it is a particular case of torus routine.
//           Three planes with two constant incident fillets
//           of the same radius and the third face perpendicular to
//           two others are required.
//=======================================================================

static bool ToricRotule(const BRepAdaptor_Surface&        fac,
                        const BRepAdaptor_Surface&        s1,
                        const BRepAdaptor_Surface&        s2,
                        const occ::handle<ChFiDS_Stripe>& c1,
                        const occ::handle<ChFiDS_Stripe>& c2)

{
  double tolesp = 1.e-7;

  occ::handle<ChFiDS_FilSpine> sp1 = occ::down_cast<ChFiDS_FilSpine>(c1->Spine());
  occ::handle<ChFiDS_FilSpine> sp2 = occ::down_cast<ChFiDS_FilSpine>(c2->Spine());
  if (sp1.IsNull() || sp2.IsNull())
  {
    return false;
  }
  if (!sp1->IsConstant() || !sp2->IsConstant())
  {
    return false;
  }
  if ((fac.GetType() != GeomAbs_Plane) || (s1.GetType() != GeomAbs_Plane)
      || (s2.GetType() != GeomAbs_Plane))
  {
    return false;
  }
  gp_Dir df  = fac.Plane().Position().Direction();
  gp_Dir ds1 = s1.Plane().Position().Direction();
  gp_Dir ds2 = s2.Plane().Position().Direction();
  if (std::abs(df.Dot(ds1)) >= tolesp || std::abs(df.Dot(ds2)) >= tolesp)
  {
    return false;
  }
  double r1 = sp1->Radius();
  double r2 = sp2->Radius();
  return std::abs(r1 - r2) < tolesp;
}

//=======================================================================
// function : PieceOfArc
// purpose  : <Arc>, an edge of <Vtx>, is not an edge of <F> but of a piece
//           of the same plane beside it, across an edge of <Vtx> the two
//           share -- a wall kept in coplanar pieces (a body without
//           Refine), the split between two of them, <Split>, ending at
//           <Vtx>: that piece, <Fn>, and <F> as the shell has it, <Fsh>.
//           Or <F> is itself a piece beyond the one at <Vtx>, the split
//           an edge of <F> away from <Vtx>: <AtVtx> is then false.
//=======================================================================

static bool PieceOfArc(const TopoDS_Vertex& Vtx,
                       const TopoDS_Face&   F,
                       const TopoDS_Edge&   Arc,
                       const ChFiDS_Map&    VEMap,
                       const ChFiDS_Map&    EFMap,
                       const double         Tol,
                       TopoDS_Face&         Fsh,
                       TopoDS_Face&         Fn,
                       TopoDS_Edge&         Split,
                       bool&                AtVtx)
{
  AtVtx = true;
  if (!VEMap.Contains(Vtx) || !EFMap.Contains(Arc) || EFMap(Arc).Extent() != 2)
  {
    return false;
  }
  bool isOfVtx = false;
  for (NCollection_List<TopoDS_Shape>::Iterator itE(VEMap(Vtx)); itE.More(); itE.Next())
  {
    if (Arc.IsSame(itE.Value()))
    {
      isOfVtx = true;
    }
  }
  for (NCollection_List<TopoDS_Shape>::Iterator itF(EFMap(Arc)); itF.More(); itF.Next())
  {
    if (F.IsSame(itF.Value()))
    {
      return false;
    }
  }
  if (!isOfVtx)
  {
    return false;
  }
  for (NCollection_List<TopoDS_Shape>::Iterator itE(VEMap(Vtx)); itE.More(); itE.Next())
  {
    const TopoDS_Edge& aSplit = TopoDS::Edge(itE.Value());
    if (aSplit.IsSame(Arc) || BRep_Tool::Degenerated(aSplit) || EFMap(aSplit).Extent() != 2)
    {
      continue;
    }
    TopoDS_Face aF1, aF2;
    for (NCollection_List<TopoDS_Shape>::Iterator itF(EFMap(aSplit)); itF.More(); itF.Next())
    {
      if (F.IsSame(itF.Value()))
      {
        aF1 = TopoDS::Face(itF.Value());
      }
      else
      {
        aF2 = TopoDS::Face(itF.Value());
      }
    }
    if (aF1.IsNull() || aF2.IsNull() || !ChFi3d_CoplanarPieces(aF1, aF2, Tol))
    {
      continue;
    }
    for (NCollection_List<TopoDS_Shape>::Iterator itF(EFMap(Arc)); itF.More(); itF.Next())
    {
      if (aF2.IsSame(itF.Value()))
      {
        Fsh   = aF1;
        Fn    = aF2;
        Split = aSplit;
        return true;
      }
    }
  }
  // <F> does not come to <Vtx>
  for (TopExp_Explorer exV(F, TopAbs_VERTEX); exV.More(); exV.Next())
  {
    if (exV.Current().IsSame(Vtx))
    {
      return false;
    }
  }
  for (TopExp_Explorer exE(F, TopAbs_EDGE); exE.More(); exE.Next())
  {
    const TopoDS_Edge& aSplit = TopoDS::Edge(exE.Current());
    if (BRep_Tool::Degenerated(aSplit) || !EFMap.Contains(aSplit) || EFMap(aSplit).Extent() != 2)
    {
      continue;
    }
    TopoDS_Face aF1, aF2;
    for (NCollection_List<TopoDS_Shape>::Iterator itF(EFMap(aSplit)); itF.More(); itF.Next())
    {
      if (F.IsSame(itF.Value()))
      {
        aF1 = TopoDS::Face(itF.Value());
      }
      else
      {
        aF2 = TopoDS::Face(itF.Value());
      }
    }
    if (aF1.IsNull() || aF2.IsNull() || !ChFi3d_CoplanarPieces(aF1, aF2, Tol))
    {
      continue;
    }
    for (NCollection_List<TopoDS_Shape>::Iterator itF(EFMap(Arc)); itF.More(); itF.Next())
    {
      if (aF2.IsSame(itF.Value()))
      {
        Fsh   = aF1;
        Fn    = aF2;
        Split = aSplit;
        AtVtx = false;
        return true;
      }
    }
  }
  return false;
}

static void RemoveSD(occ::handle<ChFiDS_Stripe>& Stripe, const int num1, const int num2)
{
  NCollection_Sequence<occ::handle<ChFiDS_SurfData>>& Seq =
    Stripe->ChangeSetOfSurfData()->ChangeSequence();
  if (Seq.IsEmpty())
  {
    return;
  }
  if (num1 == num2)
  {
    Seq.Remove(num1);
  }
  else
  {
    Seq.Remove(num1, num2);
  }
}

//=================================================================================================

void ChFi3d_FilBuilder::PerformTwoCorner(const int Index)
{
#ifdef OCCT_DEBUG
  OSD_Chronometer ch;
  ChFi3d_InitChron(ch); // init perf initialisation
#endif

  done                                                        = false;
  const TopoDS_Vertex&                                   Vtx  = myVDataMap.FindKey(Index);
  TopOpeBRepDS_DataStructure&                            DStr = myDS->ChangeDS();
  NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator It;
  It.Initialize(myVDataMap(Index));
  occ::handle<ChFiDS_Stripe>                         st1, st2;
  int                                                Sens1, Sens2;
  int                                                Isd1, Isd2, i1, i2;
  occ::handle<ChFiDS_SurfData>                       sd1, sd2;
  NCollection_Sequence<occ::handle<ChFiDS_SurfData>> SeqFil1, SeqFil2;
  occ::handle<Geom_Surface>                          surf1, surf2;
  bool                                               OkinterCC, Okvisavis, SameSide;
  int                                                IFaCo1, IFaCo2;
  double                                             UIntPC1, UIntPC2;
  TopoDS_Face                                        FaCo;
  TopoDS_Edge                                        E1, E2, E;
  TopoDS_Vertex                                      V1, V2;
  //  gp_Pnt P1,P2;
  int    nbsurf1, nbsurf2, deb1, fin1, deb2, fin2;
  double parE1, parE2;
  // Necessary information on fillets is extracted
  //------------------------------------------------------

  // the first
  //----------

  st1  = It.Value();
  Isd1 = ChFi3d_IndexOfSurfData(Vtx, st1, Sens1);

  // the second
  //----------
  It.Next();
  st2 = It.Value();
  if (st2 == st1)
  {
    Sens2 = -1;
    Isd2  = st2->SetOfSurfData()->Length();
  }
  else
  {
    Isd2 = ChFi3d_IndexOfSurfData(Vtx, st2, Sens2);
  }

  // If two edges to rounded are tangent GeomPlate is called

  if (Sens1 == 1)
  {
    E1 = st1->Spine()->Edges(1);
  }
  else
  {
    E1 = st1->Spine()->Edges(st1->Spine()->NbEdges());
  }

  if (Sens2 == 1)
  {
    E2 = st2->Spine()->Edges(1);
  }
  else
  {
    E2 = st2->Spine()->Edges(st2->Spine()->NbEdges());
  }

  BRepAdaptor_Curve BCurv1(E1);
  BRepAdaptor_Curve BCurv2(E2);
  parE1 = BRep_Tool::Parameter(Vtx, E1);
  parE2 = BRep_Tool::Parameter(Vtx, E2);
  BRepLProp_CLProps CL1(BCurv1, parE1, 1, 1.e-4);
  BRepLProp_CLProps CL2(BCurv2, parE2, 1, 1.e-4);
  gp_Dir            dir1, dir2;
  CL1.Tangent(dir1);
  CL2.Tangent(dir2);
  if (Sens1 == -1)
  {
    dir1.Reverse();
  }
  if (Sens2 == -1)
  {
    dir2.Reverse();
  }
  double ang1;
  ang1 = std::abs(dir1.Angle(dir2));
  if (ang1 < M_PI / 180.)
  {
    PerformMoreThreeCorner(Index, 2);
    done = true;
    return;
  }

  OkinterCC = ChFi3d_IsInFront(DStr,
                               st1,
                               st2,
                               Isd1,
                               Isd2,
                               Sens1,
                               Sens2,
                               UIntPC1,
                               UIntPC2,
                               FaCo,
                               SameSide,
                               IFaCo1,
                               IFaCo2,
                               Okvisavis,
                               Vtx,
                               true);

  // The two fillets on the two pieces of a wall kept in coplanar pieces (a
  // body without Refine), the split between them ending at Vtx -- the
  // reentrant corner of an L-shaped wall made of two faces: they have no
  // face in common, and the corner was a plate. On the wall in one face the
  // two lines on it meet past the split, on one of the pieces, and the
  // fillets are cut by one another there. The same is done here: while they
  // are intersected, the line that stops on the split has the other's piece
  // for its face; it then ends on the split, and its piece from there to
  // the meeting point is put on the other piece.
  // The same where the split passes Vtx by, one fillet as wide as or wider
  // than the piece its spine is on, its line on the piece beyond (a pocket's
  // mouth beside a split): the other's line stops short of the split or on
  // it, and crosses it to the meeting point.
  if (!Okvisavis && st1 != st2)
  {
    const bool                          isf[2] = {Sens1 == 1, Sens2 == 1};
    const occ::handle<ChFiDS_Stripe>    aSt[2] = {st1, st2};
    const occ::handle<ChFiDS_SurfData> aSd[2] = {st1->SetOfSurfData()->Value(Isd1),
                                                  st2->SetOfSurfData()->Value(Isd2)};
    int         aJ[2] = {0, 0};
    TopoDS_Face aFsh[2];
    TopoDS_Edge aSplit;
    for (int j1 = 1; j1 <= 2 && aSplit.IsNull(); j1++)
    {
      for (int j2 = 1; j2 <= 2 && aSplit.IsNull(); j2++)
      {
        if (aSd[0]->Index(j1) <= 0 || aSd[1]->Index(j2) <= 0)
        {
          continue;
        }
        const TopoDS_Shape& aS1 = DStr.Shape(aSd[0]->Index(j1));
        const TopoDS_Shape& aS2 = DStr.Shape(aSd[1]->Index(j2));
        if (aS1.IsNull() || aS2.IsNull() || aS1.ShapeType() != TopAbs_FACE
            || aS2.ShapeType() != TopAbs_FACE || aS1.IsSame(aS2))
        {
          continue;
        }
        // an edge of Vtx, or else one of the face
        NCollection_List<TopoDS_Shape> aCand;
        for (NCollection_List<TopoDS_Shape>::Iterator itE(myVEMap(Vtx)); itE.More(); itE.Next())
        {
          aCand.Append(itE.Value());
        }
        for (TopExp_Explorer exE(aS1, TopAbs_EDGE); exE.More(); exE.Next())
        {
          aCand.Append(exE.Current());
        }
        for (NCollection_List<TopoDS_Shape>::Iterator itE(aCand); itE.More(); itE.Next())
        {
          const TopoDS_Edge& anE = TopoDS::Edge(itE.Value());
          if (BRep_Tool::Degenerated(anE) || !myEFMap.Contains(anE)
              || myEFMap(anE).Extent() != 2)
          {
            continue;
          }
          TopoDS_Face aFa, aFb;
          for (NCollection_List<TopoDS_Shape>::Iterator itF(myEFMap(anE)); itF.More(); itF.Next())
          {
            if (aS1.IsSame(itF.Value()))
            {
              aFa = TopoDS::Face(itF.Value());
            }
            else if (aS2.IsSame(itF.Value()))
            {
              aFb = TopoDS::Face(itF.Value());
            }
          }
          if (aFa.IsNull() || aFb.IsNull() || !ChFi3d_CoplanarPieces(aFa, aFb, 10 * tolapp3d))
          {
            continue;
          }
          aSplit  = anE;
          aJ[0]   = j1;
          aJ[1]   = j2;
          aFsh[0] = aFa;
          aFsh[1] = aFb;
          break;
        }
      }
    }
    // the one line stops on the split, the other does not
    int ic = -1;
    if (!aSplit.IsNull())
    {
      for (int k = 0; k < 2; k++)
      {
        const ChFiDS_CommonPoint& aCP = aSd[k]->Vertex(isf[k], aJ[k]);
        if (aCP.IsOnArc() && aCP.Arc().IsSame(aSplit))
        {
          ic = (ic == -1 ? k : -2);
        }
      }
    }
    if (!aSplit.IsNull() && ic != -2 && st1->Spine()->Status(isf[0]) == ChFiDS_OnSame
        && st2->Spine()->Status(isf[1]) == ChFiDS_OnSame)
    {
      // both fillets end on one edge on their other sides: at one point of
      // it, or at two (radii or faces that differ), as the fillets on a wall
      // in one face that PerformTwoCornerbyInter takes
      const ChFiDS_CommonPoint& aCPa1 = aSd[0]->Vertex(isf[0], 3 - aJ[0]);
      const ChFiDS_CommonPoint& aCPa2 = aSd[1]->Vertex(isf[1], 3 - aJ[1]);
      bool isOk = aCPa1.IsOnArc() && aCPa2.IsOnArc() && aCPa1.Arc().IsSame(aCPa2.Arc());
      // two straight lines that meet; a walked fillet's is a B-spline, and
      // straight all the same on a plane along a straight spine
      occ::handle<Geom_Line> aLin[2];
      double                 aT[2] = {0., 0.};
      gp_Pnt                 aX;
      if (isOk)
      {
        gp_Lin aL[2];
        for (int k = 0; k < 2 && isOk; k++)
        {
          const ChFiDS_FaceInterference& aFik = aSd[k]->Interference(aJ[k]);
          isOk = aFik.LineIndex() > 0 && !DStr.Curve(aFik.LineIndex()).Curve().IsNull();
          if (!isOk)
          {
            break;
          }
          const occ::handle<Geom_Curve>& aCk = DStr.Curve(aFik.LineIndex()).Curve();
          const GeomAdaptor_Curve        aGk(aCk);
          if (aGk.GetType() == GeomAbs_Line)
          {
            aL[k]   = aGk.Line();
            aLin[k] = new Geom_Line(aL[k]);
            continue;
          }
          const double aWa = aFik.FirstParameter(), aWb = aFik.LastParameter();
          const gp_Pnt aPa = aCk->Value(aWa), aPb = aCk->Value(aWb);
          isOk             = aPa.Distance(aPb) > 10. * tolapp3d;
          if (isOk)
          {
            aL[k] = gp_Lin(aPa, gp_Dir(gp_Vec(aPa, aPb)));
            for (int i = 1; i <= 7 && isOk; i++)
            {
              isOk = aL[k].Distance(aCk->Value(aWa + i * (aWb - aWa) / 8.)) <= tolapp3d;
            }
          }
        }
        if (isOk)
        {
          const gp_Vec aD0(aL[0].Direction()), aD1(aL[1].Direction());
          const gp_Vec aN  = aD0.Crossed(aD1);
          const double aN2 = aN.SquareMagnitude();
          isOk             = aN2 > 1.e-12;
          if (isOk)
          {
            const gp_Vec aPP(aL[0].Location(), aL[1].Location());
            aT[0] = aPP.Crossed(aD1).Dot(aN) / aN2;
            aT[1] = aPP.Crossed(aD0).Dot(aN) / aN2;
            aX    = aL[0].Location().Translated(aT[0] * aD0);
            isOk  = aX.Distance(aL[1].Location().Translated(aT[1] * aD1)) <= 10. * tolapp3d;
          }
        }
      }
      // The line that stops on the split; or, neither on it, the one that
      // is to cross it: the two meet on the other's piece.
      double aWQ = 0., aParQ = 0.;
      if (isOk && ic >= 0)
      {
        aWQ   = aSd[ic]->Interference(aJ[ic]).Parameter(isf[ic]);
        aParQ = aSd[ic]->Vertex(isf[ic], aJ[ic]).ParameterOnArc();
      }
      else if (isOk)
      {
        for (int k = 0; k < 2 && ic < 0; k++)
        {
          BRepClass_FaceClassifier aCls(aFsh[k], aX, tolapp3d);
          if (aCls.State() == TopAbs_IN)
          {
            ic = 1 - k;
          }
        }
        isOk = false;
        if (ic >= 0 && !aLin[ic].IsNull())
        {
          const GeomAdaptor_Curve aGL(aLin[ic]);
          const BRepAdaptor_Curve aSC(aSplit);
          const double            aRes = aSC.Resolution(10. * tolapp3d);
          Extrema_ExtCC           anExt(aGL, aSC);
          for (int i = 1; anExt.IsDone() && !anExt.IsParallel() && i <= anExt.NbExt() && !isOk;
               i++)
          {
            Extrema_POnCurv aP1, aP2;
            anExt.Points(i, aP1, aP2);
            aWQ   = aP1.Parameter();
            aParQ = aP2.Parameter();
            isOk  = anExt.SquareDistance(i) <= 100. * tolapp3d * tolapp3d
                    && aParQ > aSC.FirstParameter() + aRes && aParQ < aSC.LastParameter() - aRes;
          }
        }
      }
      // the line that is carried on is a line as the DS has it
      isOk                       = isOk && !aLin[ic].IsNull();
      const int              io  = isOk ? 1 - ic : 0;
      occ::handle<Geom_Line> aLc = isOk ? aLin[ic] : occ::handle<Geom_Line>();
      const double           aWX = isOk ? aT[ic] : 0.;
      if (isOk)
      {
        // past the split, in the other's piece; the fillet is carried on
        // past its end to there, on a surface that goes on as it is
        isOk = isf[ic] ? aWX < aWQ - tolapp3d : aWX > aWQ + tolapp3d;
        if (isOk)
        {
          BRepClass_FaceClassifier aCls(aFsh[io], aX, tolapp3d);
          isOk = aCls.State() == TopAbs_IN
                 && DStr.Surface(aSd[ic]->Surf()).Surface()->IsKind(
                   STANDARD_TYPE(Geom_ElementarySurface));
        }
      }
      if (isOk)
      {
        ChFiDS_FaceInterference&        aFic     = aSd[ic]->ChangeInterference(aJ[ic]);
        const int                       aIndSave = aSd[ic]->Index(aJ[ic]);
        const occ::handle<Geom2d_Curve> aPcSave  = aFic.PCurveOnFace();
        const double                    aWSave   = aFic.Parameter(isf[ic]);
        // the two pieces the other way up on their surfaces: the line's
        // transition, told by its face's surface, is the other one there
        const bool                      isTurned = aFsh[0].Orientation() != aFsh[1].Orientation();
        const occ::handle<Geom2d_Curve> aPcOn    = GeomProjLib::Curve2d(
          aLc,
          BRep_Tool::Surface(TopoDS::Face(DStr.Shape(aSd[io]->Index(aJ[io])))));
        isOk = !aPcOn.IsNull();
        if (isOk)
        {
          auto aSetFace = [&](const int                        theInd,
                              const occ::handle<Geom2d_Curve>& thePc,
                              const double                     theW) {
            if (aJ[ic] == 1)
            {
              aSd[ic]->ChangeIndexOfS1(theInd);
            }
            else
            {
              aSd[ic]->ChangeIndexOfS2(theInd);
            }
            aFic.ChangePCurveOnFace() = thePc;
            aFic.SetParameter(theW, isf[ic]);
            if (isTurned)
            {
              aFic.SetTransition(TopAbs::Reverse(aFic.Transition()));
            }
          };
          // the line as long as to pass the meeting point, as it is on the
          // wall in one face
          aSetFace(aSd[io]->Index(aJ[io]), aPcOn, 2. * aWX - aWQ);
          bool isDone = false;
          try
          {
            OCC_CATCH_SIGNALS
            isDone = PerformTwoCornerbyInter(Index);
          }
          catch (Standard_Failure const&)
          {
            aSetFace(aIndSave, aPcSave, aWSave);
            throw;
          }
          // the line is cut where they meet, when they are
          const double aWP = isDone ? aFic.Parameter(isf[ic]) : aWSave;
          aSetFace(aIndSave, aPcSave, aWP);
          if (isDone)
          {
            // The end of that fillet -- the curve the two are cut along --
            // goes on its surface here, as FILDS would put it; its line
            // ends on the split, and FILDS ends it there.
            const int          Isurf = aSd[ic]->Surf();
            const int          indP  = aSt[ic]->IndexPoint(isf[ic], aJ[ic]);
            TopAbs_Orientation aTra  = TopAbs_FORWARD;
            if (aSd[ic]->IndexOfS1() > 0)
            {
              aTra = DStr.Shape(aSd[ic]->IndexOfS1()).Orientation();
            }
            aTra = TopAbs::Compose(aTra, aSd[ic]->Orientation());
            aTra =
              TopAbs::Compose(TopAbs::Reverse(aSd[ic]->InterferenceOnS1().Transition()), aTra);
            DStr.ChangeSurfaceInterferences(Isurf).Append(
              ChFi3d_FilCurveInDS(aSt[ic]->Curve(isf[ic]),
                                  Isurf,
                                  aSt[ic]->PCurve(isf[ic]),
                                  TopAbs::Compose(isf[ic] ? TopAbs::Reverse(aTra) : aTra,
                                                  aSt[ic]->Orientation(isf[ic]))));
            ChFi3d_StoreLineOverSplit(DStr,
                                      aSt[ic],
                                      aSd[ic],
                                      isf[ic],
                                      aJ[ic],
                                      Isurf,
                                      new GeomAdaptor_Surface(DStr.Surface(Isurf).Surface()),
                                      aLc,
                                      aSplit,
                                      aFsh[io],
                                      aWQ,
                                      aWP,
                                      aParQ,
                                      indP,
                                      10 * tolapp3d,
                                      myRegul);
            aSt[ic]->InDS(isf[ic]);
            done = true;
            return;
          }
        }
      }
    }
  }

  // Two fillets of opposite concavity facing one another on their last
  // surfaces, their lines on the common face not meeting there: one of
  // them ends in two surfaces -- its line crosses a split of its wall (a
  // body without Refine) a little before the spine's end -- and the lines
  // meet on the surface before, where the wall in one face has them meet.
  // The corner is made of that surface, the one past it gone as what the
  // corner cuts off the fillet anyway.
  if (Okvisavis && !OkinterCC && st1 != st2
      && ((st1->Spine()->Status(Sens1 == 1) == ChFiDS_OnSame
           && st2->Spine()->Status(Sens2 == 1) == ChFiDS_OnDiff)
          || (st1->Spine()->Status(Sens1 == 1) == ChFiDS_OnDiff
              && st2->Spine()->Status(Sens2 == 1) == ChFiDS_OnSame)))
  {
    for (int k = 1; k <= 2; k++)
    {
      occ::handle<ChFiDS_Stripe>& aSt   = k == 1 ? st1 : st2;
      const int                   aSens = k == 1 ? Sens1 : Sens2;
      const int                   aNb   = aSt->SetOfSurfData()->Length();
      if (aNb < 2 || aSt->Spine()->Status(aSens == 1) != ChFiDS_OnSame)
      {
        continue;
      }
      const int   aPrev = aSens == 1 ? 2 : aNb - 1;
      double      aU1 = 0., aU2 = 0.;
      TopoDS_Face aFa;
      bool        aSS = false, aVis = false;
      int         aJ1 = 0, aJ2 = 0;
      if (ChFi3d_IsInFront(DStr,
                           st1,
                           st2,
                           k == 1 ? aPrev : Isd1,
                           k == 2 ? aPrev : Isd2,
                           Sens1,
                           Sens2,
                           aU1,
                           aU2,
                           aFa,
                           aSS,
                           aJ1,
                           aJ2,
                           aVis,
                           Vtx,
                           true)
          && aVis)
      {
        RemoveSD(aSt, aSens == 1 ? 1 : aNb, aSens == 1 ? 1 : aNb);
        Isd1      = ChFi3d_IndexOfSurfData(Vtx, st1, Sens1);
        Isd2      = ChFi3d_IndexOfSurfData(Vtx, st2, Sens2);
        OkinterCC = true;
        Okvisavis = true;
        UIntPC1   = aU1;
        UIntPC2   = aU2;
        FaCo      = aFa;
        SameSide  = aSS;
        IFaCo1    = aJ1;
        IFaCo2    = aJ2;
        break;
      }
    }
  }

  // The fillet OnDiff in two surfaces at its end, for the two pieces of its
  // edge: the wall beside it is kept in coplanar pieces (a body without
  // Refine), and the split's foot is on the edge, nearer the vertex than the
  // fillet is wide. The lines meet on the common face past that foot, on
  // the surface before the last, and the corner runs from there to the last
  // one's end on an edge of the vertex: over both. On the wall in one face
  // they are one surface. They are made one here: the surface before the
  // last carries on to the vertex -- its line on the common face to the
  // last one's end, its line on the wall over the split to the last one's
  // end on the piece at the vertex -- and the last goes. While the corner
  // is made the line is shown on that piece; then it ends on the split
  // again and its piece beyond is put on the piece at the vertex and on the
  // fillet (StoreLineOverSplit), as a line over a split wall.
  struct MergedEnd
  {
    occ::handle<ChFiDS_Stripe>   Stripe;
    occ::handle<ChFiDS_SurfData> Sd, SdSave, SdGone;
    occ::handle<Geom_Line>       Line;
    occ::handle<Geom2d_Curve>    PcSave;
    TopoDS_Edge                  Split;
    TopoDS_Face                  Fn;
    double                       WQ = 0., ParQ = 0.;
    int                          Side = 0, IndSave = 0;
    bool                         isFirst = false, isTurned = false, isShown = false;
    int                          Uncaught = std::uncaught_exceptions();
    //! the two surfaces again, as they were
    void Undo()
    {
      if (Sd.IsNull())
      {
        return;
      }
      Sd->Copy(SdSave);
      NCollection_Sequence<occ::handle<ChFiDS_SurfData>>& aSeq =
        Stripe->ChangeSetOfSurfData()->ChangeSequence();
      if (isFirst)
      {
        aSeq.Prepend(SdGone);
      }
      else
      {
        aSeq.Append(SdGone);
      }
      Sd.Nullify();
    }
    //! a corner that fails leaves the fillet as it found it
    ~MergedEnd()
    {
      if (std::uncaught_exceptions() > Uncaught)
      {
        Undo();
      }
    }
  } aMerged;
  if ((!Okvisavis || !OkinterCC) && st1 != st2
      && ((st1->Spine()->Status(Sens1 == 1) == ChFiDS_OnSame
           && st2->Spine()->Status(Sens2 == 1) == ChFiDS_OnDiff)
          || (st1->Spine()->Status(Sens1 == 1) == ChFiDS_OnDiff
              && st2->Spine()->Status(Sens2 == 1) == ChFiDS_OnSame)))
  {
    for (int k = 1; k <= 2 && aMerged.Sd.IsNull(); k++)
    {
      occ::handle<ChFiDS_Stripe>& aSt = k == 1 ? st1 : st2;
      const bool                  isf = (k == 1 ? Sens1 : Sens2) == 1;
      NCollection_Sequence<occ::handle<ChFiDS_SurfData>>& aSeq =
        aSt->ChangeSetOfSurfData()->ChangeSequence();
      const int aNb = aSeq.Length();
      if (aNb < 2 || aSt->Spine()->Status(isf) != ChFiDS_OnDiff)
      {
        continue;
      }
      const occ::handle<ChFiDS_SurfData> aSdE = aSeq(isf ? 1 : aNb);
      const occ::handle<ChFiDS_SurfData> aSdM = aSeq(isf ? 2 : aNb - 1);
      // one face on the one side, two coplanar pieces on the other
      int jc = 0;
      for (int j = 1; j <= 2; j++)
      {
        if (aSdE->Index(j) > 0 && aSdE->Index(j) == aSdM->Index(j))
        {
          jc = jc == 0 ? j : -1;
        }
      }
      if (jc <= 0)
      {
        continue;
      }
      const int jo = 3 - jc;
      if (aSdE->Index(jo) <= 0 || aSdM->Index(jo) <= 0
          || DStr.Shape(aSdE->Index(jo)).ShapeType() != TopAbs_FACE
          || DStr.Shape(aSdM->Index(jo)).ShapeType() != TopAbs_FACE
          || DStr.Shape(aSdM->Index(jc)).ShapeType() != TopAbs_FACE)
      {
        continue;
      }
      const TopoDS_Face aFm = TopoDS::Face(DStr.Shape(aSdM->Index(jo)));
      const TopoDS_Face aFe = TopoDS::Face(DStr.Shape(aSdE->Index(jo)));
      const TopoDS_Face aFc = TopoDS::Face(DStr.Shape(aSdM->Index(jc)));
      if (!ChFi3d_CoplanarPieces(aFm, aFe, 10 * tolapp3d))
      {
        continue;
      }
      // two straight lines, on a surface that goes on past its end as it is
      ChFiDS_FaceInterference&        aFimO = aSdM->ChangeInterference(jo);
      ChFiDS_FaceInterference&        aFimC = aSdM->ChangeInterference(jc);
      const ChFiDS_CommonPoint&       aCPp  = aSdE->Vertex(isf, jo);
      const ChFiDS_CommonPoint&       aCPc  = aSdE->Vertex(isf, jc);
      const occ::handle<Geom_Surface> aSm   = DStr.Surface(aSdM->Surf()).Surface();
      if (aFimO.LineIndex() == 0 || aFimC.LineIndex() == 0
          || DStr.Curve(aFimO.LineIndex()).Curve().IsNull()
          || DStr.Curve(aFimC.LineIndex()).Curve().IsNull() || aFimO.PCurveOnSurf().IsNull()
          || aFimC.PCurveOnSurf().IsNull() || aFimC.PCurveOnFace().IsNull() || aSm.IsNull()
          || !aSm->IsKind(STANDARD_TYPE(Geom_ElementarySurface)) || !aCPp.IsOnArc()
          || aCPp.IsVertex())
      {
        continue;
      }
      const GeomAdaptor_Curve aGO(DStr.Curve(aFimO.LineIndex()).Curve());
      const GeomAdaptor_Curve aGC(DStr.Curve(aFimC.LineIndex()).Curve());
      if (aGO.GetType() != GeomAbs_Line || aGC.GetType() != GeomAbs_Line)
      {
        continue;
      }
      // the last surface's two ends, on the lines of the one before
      const double aWQ = aFimO.Parameter(isf);
      const double aWP = ElCLib::Parameter(aGO.Line(), aCPp.Point());
      const double aWC = ElCLib::Parameter(aGC.Line(), aCPc.Point());
      const double aWM = aFimC.Parameter(isf);
      if (aGO.Value(aWP).Distance(aCPp.Point()) > tolapp3d
          || aGC.Value(aWC).Distance(aCPc.Point()) > tolapp3d
          || (isf ? aWP > aWQ - tolapp3d : aWP < aWQ + tolapp3d)
          || (isf ? aWC > aWM - tolapp3d : aWC < aWM + tolapp3d))
      {
        continue;
      }
      // and on its surface, and on the common face, carried on
      const gp_Pnt2d            aUVsO = aFimO.PCurveOnSurf()->Value(aWP);
      const gp_Pnt2d            aUVsC = aFimC.PCurveOnSurf()->Value(aWC);
      const gp_Pnt2d            aUVfC = aFimC.PCurveOnFace()->Value(aWC);
      const BRepAdaptor_Surface aBSc(aFc, false);
      if (aSm->Value(aUVsO.X(), aUVsO.Y()).Distance(aCPp.Point()) > tolapp3d
          || aSm->Value(aUVsC.X(), aUVsC.Y()).Distance(aCPc.Point()) > tolapp3d
          || aBSc.Value(aUVfC.X(), aUVfC.Y()).Distance(aCPc.Point()) > tolapp3d)
      {
        continue;
      }
      // the split the line stops on, between the two pieces
      const occ::handle<Geom_Line> aLc = new Geom_Line(aGO.Line());
      const gp_Pnt                 aPQ = aLc->Value(aWQ);
      TopoDS_Edge                  aEs;
      double                       aParQ = 0.;
      for (TopExp_Explorer exE(aFm, TopAbs_EDGE); exE.More() && aEs.IsNull(); exE.Next())
      {
        const TopoDS_Edge& aE = TopoDS::Edge(exE.Current());
        if (BRep_Tool::Degenerated(aE))
        {
          continue;
        }
        bool isBoth = false;
        for (TopExp_Explorer exN(aFe, TopAbs_EDGE); exN.More() && !isBoth; exN.Next())
        {
          isBoth = exN.Current().IsSame(aE);
        }
        if (!isBoth)
        {
          continue;
        }
        const BRepAdaptor_Curve aCE(aE);
        Extrema_ExtPC           anExt(aPQ, aCE);
        for (int i = 1; anExt.IsDone() && i <= anExt.NbExt() && aEs.IsNull(); i++)
        {
          if (anExt.SquareDistance(i) <= 100. * tolapp3d * tolapp3d)
          {
            aEs   = aE;
            aParQ = anExt.Point(i).Parameter();
          }
        }
      }
      const occ::handle<Geom2d_Curve> aPcOn =
        aEs.IsNull() ? occ::handle<Geom2d_Curve>()
                     : GeomProjLib::Curve2d(aLc, BRep_Tool::Surface(aFe));
      if (aPcOn.IsNull())
      {
        continue;
      }

      // made one
      aMerged.Stripe = aSt;
      aMerged.SdSave = new ChFiDS_SurfData();
      aMerged.SdSave->Copy(aSdM);
      aMerged.SdGone   = aSdE;
      aMerged.Sd       = aSdM;
      aMerged.Line     = aLc;
      aMerged.PcSave   = aFimO.PCurveOnFace();
      aMerged.Split    = aEs;
      aMerged.Fn       = aFe;
      aMerged.WQ       = aWQ;
      aMerged.ParQ     = aParQ;
      aMerged.Side     = jo;
      aMerged.IndSave  = aSdM->Index(jo);
      aMerged.isFirst  = isf;
      aMerged.isTurned = aFe.Orientation() != aFm.Orientation();
      aSeq.Remove(isf ? 1 : aNb);
      // Twice. The line shown on the piece at the vertex, to the last
      // surface's end there, for fillets that meet on the face they both
      // lie on whole. Or left on its own piece, to the split: the other
      // fillet, wider than the piece at the vertex, lies on that piece too,
      // the lines meet on it before the split, and the corner cuts the
      // line there.
      for (int iShow = 1; iShow >= 0 && !aMerged.Sd.IsNull(); iShow--)
      {
        aSdM->Copy(aMerged.SdSave);
        gp_Pnt2d a2d[4];
        aSdM->Get2dPoints(a2d[0], a2d[1], a2d[2], a2d[3]);
        aFimC.SetParameter(aWC, isf);
        aSdM->ChangeVertex(isf, jc)       = aCPc;
        a2d[2 * (jc - 1) + (isf ? 0 : 1)] = aUVsC;
        if (iShow == 1)
        {
          if (jo == 1)
          {
            aSdM->ChangeIndexOfS1(aSdE->Index(jo));
          }
          else
          {
            aSdM->ChangeIndexOfS2(aSdE->Index(jo));
          }
          aFimO.ChangePCurveOnFace() = aPcOn;
          aFimO.SetParameter(aWP, isf);
          if (aMerged.isTurned)
          {
            aFimO.SetTransition(TopAbs::Reverse(aFimO.Transition()));
          }
          aSdM->ChangeVertex(isf, jo)       = aCPp;
          a2d[2 * (jo - 1) + (isf ? 0 : 1)] = aUVsO;
        }
        aSdM->Set2dPoints(a2d[0], a2d[1], a2d[2], a2d[3]);
        if (isf)
        {
          aSdM->FirstSpineParam(aSdE->FirstSpineParam());
        }
        else
        {
          aSdM->LastSpineParam(aSdE->LastSpineParam());
        }
        // do the lines meet now, on the face they should?
        double      aU1 = 0., aU2 = 0.;
        TopoDS_Face aFa;
        bool        aSS = false, aVis = false;
        int         aJ1 = 0, aJ2 = 0;
        const int   aI1 = ChFi3d_IndexOfSurfData(Vtx, st1, Sens1);
        const int   aI2 = ChFi3d_IndexOfSurfData(Vtx, st2, Sens2);
        if (ChFi3d_IsInFront(DStr,
                             st1,
                             st2,
                             aI1,
                             aI2,
                             Sens1,
                             Sens2,
                             aU1,
                             aU2,
                             aFa,
                             aSS,
                             aJ1,
                             aJ2,
                             aVis,
                             Vtx,
                             true)
            && aVis && (k == 1 ? aJ1 : aJ2) == (iShow == 1 ? jc : jo)
            // (on its own piece: short of the split, not on it -- a fillet as
            // wide as the piece at the vertex has its line along the split)
            && (iShow == 1
                || (isf ? (k == 1 ? aU1 : aU2) > aWQ + 10 * tolapp3d
                        : (k == 1 ? aU1 : aU2) < aWQ - 10 * tolapp3d)))
        {
          aMerged.isShown = iShow == 1;
          Isd1            = aI1;
          Isd2            = aI2;
          OkinterCC       = true;
          Okvisavis       = true;
          UIntPC1         = aU1;
          UIntPC2         = aU2;
          FaCo            = aFa;
          SameSide        = aSS;
          IFaCo1          = aJ1;
          IFaCo2          = aJ2;
          break;
        }
        if (iShow == 0)
        {
          aMerged.Undo();
        }
      }
    }
  }

  bool trouve = false;
  if (!Okvisavis)
  {

    // one is not limited to the first or the last surfdata
    // to find the opposing data
    nbsurf1 = st1->SetOfSurfData()->Length();
    nbsurf2 = st2->SetOfSurfData()->Length();
    deb1    = 1;
    deb2    = 1;
    fin1    = 1;
    fin2    = 1;
    if (nbsurf1 != 1)
    {
      if (Sens1 == 1)
      {
        deb1 = 1;
        fin1 = 2;
      }
      else
      {
        deb1 = nbsurf1 - 1;
        fin1 = nbsurf1;
      }
    }
    if (nbsurf2 != 1)
    {
      if (Sens2 == 1)
      {
        deb2 = 1;
        fin2 = 2;
      }
      else
      {
        deb2 = nbsurf2 - 1;
        fin2 = nbsurf2;
      }
    }

    for (i1 = deb1; i1 <= fin1 && !trouve; i1++)
    {
      Isd1 = i1;
      for (i2 = deb2; i2 <= fin2 && !trouve; i2++)
      {
        Isd2 = i2;

        OkinterCC = ChFi3d_IsInFront(DStr,
                                     st1,
                                     st2,
                                     Isd1,
                                     Isd2,
                                     Sens1,
                                     Sens2,
                                     UIntPC1,
                                     UIntPC2,
                                     FaCo,
                                     SameSide,
                                     IFaCo1,
                                     IFaCo2,
                                     Okvisavis,
                                     Vtx,
                                     true);
        trouve    = Okvisavis;
        // Face to face on a face the vertex is not on: the surface before
        // the last of a fillet whose line crosses a split of its wall (a
        // body without Refine) before its end, and the other fillet on the
        // piece beyond. This corner is not made of a surface that does not
        // come to the vertex -- it would be a plate from where that surface
        // ends; it is refused, as it was while ChFi3d_IsInFront asked for
        // the vertex's parameters on the face.
        if (trouve)
        {
          bool isAtVtx = false;
          for (TopExp_Explorer exV(FaCo, TopAbs_VERTEX); exV.More() && !isAtVtx; exV.Next())
          {
            isAtVtx = exV.Current().IsSame(Vtx);
          }
          if (!isAtVtx)
          {
            throw Standard_Failure("TwoCorner : fillets face to face on a face off the vertex");
          }
        }
      }
    }
    if (!trouve)
    {
      PerformMoreThreeCorner(Index, 2);
      done = true;
      return;
    }
    else
    {
      if (Sens1 == 1 && Isd1 != 1)
      {
        RemoveSD(st1, 1, 1);
      }
      if (Sens1 != 1 && Isd1 != nbsurf1)
      {
        RemoveSD(st1, fin1, fin1);
      }
      if (Sens2 == 1 && Isd2 != 1)
      {
        RemoveSD(st2, 1, 1);
      }
      if (Sens2 != 1 && Isd2 != nbsurf2)
      {
        RemoveSD(st2, fin2, fin2);
      }
    }
    Isd1 = ChFi3d_IndexOfSurfData(Vtx, st1, Sens1);
    Isd2 = ChFi3d_IndexOfSurfData(Vtx, st2, Sens2);
  }
  // throw StdFail_NotDone("TwoCorner : no common face");
  int IFaArc1 = 3 - IFaCo1, IFaArc2 = 3 - IFaCo2;
  SeqFil1                  = st1->ChangeSetOfSurfData()->ChangeSequence();
  SeqFil2                  = st2->ChangeSetOfSurfData()->ChangeSequence();
  sd1                      = SeqFil1.ChangeValue(Isd1);
  surf1                    = DStr.Surface(sd1->Surf()).Surface();
  sd2                      = SeqFil2.ChangeValue(Isd2);
  surf2                    = DStr.Surface(sd2->Surf()).Surface();
  TopAbs_Orientation OFaCo = FaCo.Orientation();
  // The concavities are analyzed and the opposite face and the
  // eventual intersection of 2 pcurves on this face are found.

  ChFiDS_State Stat1, Stat2;
  bool         isfirst1 = (Sens1 == 1);
  bool         isfirst2 = (Sens2 == 1);
  Stat1                 = st1->Spine()->Status(isfirst1);
  Stat2                 = st2->Spine()->Status(isfirst2);
  bool c1biseau         = (Stat1 == ChFiDS_AllSame);
  bool c1rotule         = (Stat1 == ChFiDS_OnSame && Stat2 == ChFiDS_OnSame);

  // It is checked if the fillets have a commonpoint on a common arc.
  // This edge is the pivot of the bevel or the knee.

  ChFiDS_CommonPoint& CP1 = sd1->ChangeVertex(isfirst1, IFaArc1);
  ChFiDS_CommonPoint& CP2 = sd2->ChangeVertex(isfirst2, IFaArc2);

  // The line of one fillet on the arc's side stopped on a split of its
  // face, kept in coplanar pieces (a body without Refine), the other's on
  // an edge of the piece beyond: on the face in one piece both end on that
  // edge, the pivot, and the fillets are cut by one another
  // (PerformTwoCornerbyInter). The line is carried over the split to the
  // pivot for that, on the piece beyond; then it ends on the split again,
  // and its piece from there to the pivot is put on the piece beyond and
  // on its own fillet (StoreLineOverSplit), as a line over a split wall.
  // The same for a line that ends in its face, short of the split, where
  // the fillet's spine ends: straight on, it crosses the split.
  if (st1 != st2 && (CP1.IsOnArc() || CP2.IsOnArc())
      && !(CP1.IsOnArc() && CP2.IsOnArc() && CP1.Arc().IsSame(CP2.Arc())))
  {
    auto hasEdge = [](const TopoDS_Face& theF, const TopoDS_Edge& theE) {
      for (TopExp_Explorer exE(theF, TopAbs_EDGE); exE.More(); exE.Next())
      {
        if (exE.Current().IsSame(theE))
        {
          return true;
        }
      }
      return false;
    };
    for (int k = 0; k < 2; k++)
    {
      const occ::handle<ChFiDS_SurfData>& aSdk = k == 0 ? sd1 : sd2;
      const occ::handle<ChFiDS_SurfData>& aSdo = k == 0 ? sd2 : sd1;
      const occ::handle<ChFiDS_Stripe>&   aStk = k == 0 ? st1 : st2;
      const bool                          isfk = k == 0 ? isfirst1 : isfirst2;
      const int                           jk   = k == 0 ? IFaArc1 : IFaArc2;
      const int                           jo   = k == 0 ? IFaArc2 : IFaArc1;
      ChFiDS_CommonPoint&                 aCPk = k == 0 ? CP1 : CP2;
      const ChFiDS_CommonPoint&           aCPo = k == 0 ? CP2 : CP1;
      if (!aCPo.IsOnArc() || aSdk->Index(jk) <= 0 || aSdo->Index(jo) <= 0
          || DStr.Shape(aSdk->Index(jk)).ShapeType() != TopAbs_FACE
          || DStr.Shape(aSdo->Index(jo)).ShapeType() != TopAbs_FACE)
      {
        continue;
      }
      const TopoDS_Face aFk = TopoDS::Face(DStr.Shape(aSdk->Index(jk)));
      const TopoDS_Face aFo = TopoDS::Face(DStr.Shape(aSdo->Index(jo)));
      // a straight line, on a surface that goes on past its end as it is
      ChFiDS_FaceInterference& aFik = aSdk->ChangeInterference(jk);
      if (aFik.LineIndex() == 0 || DStr.Curve(aFik.LineIndex()).Curve().IsNull()
          || !DStr.Surface(aSdk->Surf()).Surface()->IsKind(STANDARD_TYPE(Geom_ElementarySurface)))
      {
        continue;
      }
      const GeomAdaptor_Curve aGk(DStr.Curve(aFik.LineIndex()).Curve());
      if (aGk.GetType() != GeomAbs_Line)
      {
        continue;
      }
      const occ::handle<Geom_Line> aLc = new Geom_Line(aGk.Line());
      const GeomAdaptor_Curve      aGL(aLc);
      // where the line ends now, and where it is on the split
      const double aWEnd = aFik.Parameter(isfk);
      double       aWQ   = aWEnd;
      double       aParQ = 0.;
      TopoDS_Edge  aEs;
      if (aCPk.IsOnArc())
      {
        aEs   = aCPk.Arc();
        aParQ = aCPk.ParameterOnArc();
      }
      else
      {
        // the split of its face the line crosses first, past its end
        for (TopExp_Explorer exE(aFk, TopAbs_EDGE); exE.More(); exE.Next())
        {
          const TopoDS_Edge& aE = TopoDS::Edge(exE.Current());
          if (!myEFMap.Contains(aE) || BRep_Tool::Degenerated(aE))
          {
            continue;
          }
          bool isSplit = false;
          for (NCollection_List<TopoDS_Shape>::Iterator itF(myEFMap(aE)); itF.More() && !isSplit;
               itF.Next())
          {
            const TopoDS_Face& aF = TopoDS::Face(itF.Value());
            isSplit = !aF.IsSame(aFk) && ChFi3d_CoplanarPieces(aFk, aF, 10 * tolapp3d);
          }
          if (!isSplit)
          {
            continue;
          }
          const BRepAdaptor_Curve aCE(aE);
          Extrema_ExtCC           anExtE(aGL, aCE);
          for (int i = 1; anExtE.IsDone() && !anExtE.IsParallel() && i <= anExtE.NbExt(); i++)
          {
            Extrema_POnCurv aP1, aP2;
            anExtE.Points(i, aP1, aP2);
            const double aW = aP1.Parameter(), aPar = aP2.Parameter();
            if (anExtE.SquareDistance(i) <= 100. * tolapp3d * tolapp3d
                && aPar > aCE.FirstParameter() + aCE.Resolution(tolapp3d)
                && aPar < aCE.LastParameter() - aCE.Resolution(tolapp3d)
                && (isfk ? aW < aWEnd - tolapp3d : aW > aWEnd + tolapp3d)
                && (aEs.IsNull() || std::abs(aW - aWEnd) < std::abs(aWQ - aWEnd)))
            {
              aEs   = aE;
              aWQ   = aW;
              aParQ = aPar;
            }
          }
        }
        if (aEs.IsNull())
        {
          continue;
        }
      }
      TopoDS_Edge aEp = aCPo.Arc();
      // the piece beyond the split, with the pivot for an edge
      TopoDS_Face aFn;
      if (myEFMap.Contains(aEs) && hasEdge(aFk, aEs))
      {
        for (NCollection_List<TopoDS_Shape>::Iterator itF(myEFMap(aEs)); itF.More(); itF.Next())
        {
          const TopoDS_Face& aF = TopoDS::Face(itF.Value());
          if (!aF.IsSame(aFk) && ChFi3d_CoplanarPieces(aFk, aF, 10 * tolapp3d))
          {
            aFn = aF;
          }
        }
      }
      // The pivot in two pieces itself, the split's foot between them: the
      // other fillet ends on the piece at the vertex, an edge of this
      // fillet's face, and the line carried over the split meets the piece
      // beyond, an edge of the piece of the face there. Both are edges of
      // the other fillet's face, one tangent to the other at the foot
      // (ChFi3d_PivotInTwoPieces); PerformTwoCornerbyInter takes the two for
      // one pivot.
      if (!aFn.IsNull() && !hasEdge(aFn, aEp) && hasEdge(aFk, aEp) && hasEdge(aFo, aEp))
      {
        for (TopExp_Explorer exE(aFn, TopAbs_EDGE); exE.More(); exE.Next())
        {
          const TopoDS_Edge& aE = TopoDS::Edge(exE.Current());
          if (hasEdge(aFo, aE) && ChFi3d_PivotInTwoPieces(Vtx, aEp, aE, aEs))
          {
            aEp = aE;
            break;
          }
        }
      }
      if (aFn.IsNull() || !hasEdge(aFn, aEp) || !hasEdge(aFo, aEp))
      {
        continue;
      }
      // it meets the pivot past the split
      double                  aWP = 0., aParP = 0.;
      bool                    isMet = false;
      const BRepAdaptor_Curve aPC(aEp);
      Extrema_ExtCC           anExt(aGL, aPC);
      for (int i = 1; anExt.IsDone() && !anExt.IsParallel() && i <= anExt.NbExt() && !isMet; i++)
      {
        Extrema_POnCurv aP1, aP2;
        anExt.Points(i, aP1, aP2);
        aWP   = aP1.Parameter();
        aParP = aP2.Parameter();
        isMet = anExt.SquareDistance(i) <= 100. * tolapp3d * tolapp3d
                && aParP >= aPC.FirstParameter() && aParP <= aPC.LastParameter()
                && (isfk ? aWP < aWQ - tolapp3d : aWP > aWQ + tolapp3d);
      }
      const occ::handle<Geom2d_Curve> aPcOn =
        isMet ? GeomProjLib::Curve2d(aLc, BRep_Tool::Surface(aFn)) : occ::handle<Geom2d_Curve>();
      if (aPcOn.IsNull())
      {
        continue;
      }
      const int                       aIndSave = aSdk->Index(jk);
      const occ::handle<Geom2d_Curve> aPcSave  = aFik.PCurveOnFace();
      const ChFiDS_CommonPoint        aCPSave  = aCPk;
      const bool                      isTurned = aFn.Orientation() != aFk.Orientation();
      auto aSetFace = [&](const int theInd, const occ::handle<Geom2d_Curve>& thePc, const double theW) {
        if (jk == 1)
        {
          aSdk->ChangeIndexOfS1(theInd);
        }
        else
        {
          aSdk->ChangeIndexOfS2(theInd);
        }
        aFik.ChangePCurveOnFace() = thePc;
        aFik.SetParameter(theW, isfk);
        if (isTurned)
        {
          aFik.SetTransition(TopAbs::Reverse(aFik.Transition()));
        }
      };
      aSetFace(DStr.AddShape(aFn), aPcOn, aWP);
      aCPk.Reset();
      aCPk.SetPoint(aLc->Value(aWP));
      aCPk.SetArc(std::max(aCPSave.Tolerance(), tolapp3d),
                  aEp,
                  aParP,
                  ChFi3d_TransitionOnEdge(aFik, aFn, aEp, isfk));
      bool isDone = false;
      try
      {
        OCC_CATCH_SIGNALS
        isDone = PerformTwoCornerbyInter(Index);
      }
      catch (Standard_Failure const&)
      {
        aSetFace(aIndSave, aPcSave, aWEnd);
        aCPk = aCPSave;
        throw;
      }
      if (!isDone)
      {
        aSetFace(aIndSave, aPcSave, aWEnd);
        aCPk = aCPSave;
        break;
      }
      const double aWPn = aFik.Parameter(isfk);
      // (the pivot is cut where the line ends on it: the line's end is the
      // split's from here on, and FILDS cuts that)
      const TopAbs_Orientation aTraP = aCPk.TransitionOnArc();
      aSetFace(aIndSave, aPcSave, aWPn);
      const int          Isurf = aSdk->Surf();
      const int          indP  = aStk->IndexPoint(isfk, jk);
      const int          IarcP = DStr.AddShape(aEp);
      DStr.ChangeShapeInterferences(IarcP).Append(
        ChFi3d_FilPointInDS(aTraP, IarcP, indP, aParP));
      TopAbs_Orientation aTra  = TopAbs_FORWARD;
      if (aSdk->IndexOfS1() > 0)
      {
        aTra = DStr.Shape(aSdk->IndexOfS1()).Orientation();
      }
      aTra = TopAbs::Compose(aTra, aSdk->Orientation());
      aTra = TopAbs::Compose(TopAbs::Reverse(aSdk->InterferenceOnS1().Transition()), aTra);
      DStr.ChangeSurfaceInterferences(Isurf).Append(
        ChFi3d_FilCurveInDS(aStk->Curve(isfk),
                            Isurf,
                            aStk->PCurve(isfk),
                            TopAbs::Compose(isfk ? TopAbs::Reverse(aTra) : aTra,
                                            aStk->Orientation(isfk))));
      ChFi3d_StoreLineOverSplit(DStr,
                                aStk,
                                aSdk,
                                isfk,
                                jk,
                                Isurf,
                                new GeomAdaptor_Surface(DStr.Surface(Isurf).Surface()),
                                aLc,
                                aEs,
                                aFn,
                                aWQ,
                                aWPn,
                                aParQ,
                                indP,
                                10 * tolapp3d,
                                myRegul);
      aStk->InDS(isfk);
      done = true;
      return;
    }
  }

  bool resetcp1 = false;
  bool resetcp2 = false;

  TopoDS_Edge pivot;
  bool        yapiv = false;
  if (CP1.IsOnArc())
  {
    pivot = CP1.Arc();
  }
  else
  {
    // No pivot is needed by two fillets of opposite concavity: the corner
    // is made of the surface OnSame cut where the lines meet on the common
    // face, and of the arc the surface OnDiff ends on. Where the first
    // surface is a line's that stops in its face (a walked fillet, or one
    // stopped by a split of its wall), that corner is made as where the
    // line reaches the face's edge.
    bool isInv = false;
    if (OkinterCC && !SameSide
        && ((Stat1 == ChFiDS_OnSame && Stat2 == ChFiDS_OnDiff)
            || (Stat1 == ChFiDS_OnDiff && Stat2 == ChFiDS_OnSame)))
    {
      const bool                is1Sam = Stat1 == ChFiDS_OnSame;
      const ChFiDS_CommonPoint& aCPdif = is1Sam ? CP2 : CP1;
      const int                 aIop   = is1Sam ? sd1->Index(IFaArc1) : sd2->Index(IFaArc2);
      if (aCPdif.IsOnArc() && aIop > 0 && DStr.Shape(aIop).ShapeType() == TopAbs_FACE)
      {
        const TopoDS_Face aFop = TopoDS::Face(DStr.Shape(aIop));
        for (TopExp_Explorer exF(aFop, TopAbs_EDGE); exF.More() && !isInv; exF.Next())
        {
          isInv = exF.Current().IsSame(aCPdif.Arc());
        }
        if (!isInv && !aCPdif.IsVertex())
        {
          TopoDS_Face aFa, aFb;
          TopoDS_Edge aSp;
          bool        isAt = true;
          isInv = PieceOfArc(Vtx,
                             aFop,
                             aCPdif.Arc(),
                             myVEMap,
                             myEFMap,
                             10 * tolapp3d,
                             aFa,
                             aFb,
                             aSp,
                             isAt);
        }
      }
    }
    if (!isInv)
    {
      PerformMoreThreeCorner(Index, 2);
      done = true;
      return;
    }
  }
  if (CP1.IsOnArc() && CP2.IsOnArc())
  {
    yapiv = (pivot.IsSame(CP2.Arc()));
  }
  occ::handle<BRepAdaptor_Curve> Hpivot;
  bool                           sameparam = false;
  double                         parCP1 = 0., parCP2 = 0.;
  if (yapiv)
  {
    Hpivot      = new BRepAdaptor_Curve(pivot);
    parCP1      = CP1.ParameterOnArc();
    parCP2      = CP2.ParameterOnArc();
    gp_Pnt tst1 = Hpivot->Value(parCP1);
    gp_Pnt tst2 = Hpivot->Value(parCP2);
    sameparam   = tst1.Distance(tst2) <= tolapp3d;
  }
  occ::handle<BRepAdaptor_Surface> HFaCo = new BRepAdaptor_Surface();
  occ::handle<BRepAdaptor_Surface> HFaPiv;
  occ::handle<BRepAdaptor_Surface> HBRS1 = new BRepAdaptor_Surface();
  occ::handle<BRepAdaptor_Surface> HBRS2 = new BRepAdaptor_Surface();

  BRepAdaptor_Surface& BRS1   = *HBRS1;
  BRepAdaptor_Surface& BRS2   = *HBRS2;
  BRepAdaptor_Surface& BRFaCo = *HFaCo;
  BRFaCo.Initialize(FaCo);

  TopoDS_Face                      FF1, FF2, F, FaPiv;
  TopAbs_Orientation               pctrans = TopAbs_FORWARD;
  occ::handle<Geom2d_BSplineCurve> PCurveOnPiv;
  FF1 = TopoDS::Face(DStr.Shape(sd1->Index(IFaArc1)));
  FF2 = TopoDS::Face(DStr.Shape(sd2->Index(IFaArc2)));
  if (FF1.IsNull() || FF2.IsNull())
  {
    PerformMoreThreeCorner(Index, 2);
    done = true;
    return;
  }
  BRS1.Initialize(FF1);
  BRS2.Initialize(FF2);

  if (yapiv)
  {
    NCollection_List<TopoDS_Shape>::Iterator Kt;
    bool                                     ok1 = false, ok2 = false;
    for (Kt.Initialize(myEFMap(pivot)); Kt.More(); Kt.Next())
    {
      F = TopoDS::Face(Kt.Value());
      if (!ok1 && FF1.IsSame(F))
      {
        ok1 = true;
      }
      if (!ok2 && FF2.IsSame(F))
      {
        ok2 = true;
      }
    }
    if (!ok1 || !ok2)
    {
      PerformMoreThreeCorner(Index, 2);
      done = true;
      return;
    }
  }

#ifdef OCCT_DEBUG
  ChFi3d_ResultChron(ch, t_t2cornerinit); // result perf initialisation
#endif

  // bevel
  //------
  ChFiDS_CommonPoint      cp11, cp12, cp21, cp22;
  ChFiDS_FaceInterference intf11, intf12, intf21, intf22;

  if (c1biseau)
  {
#ifdef OCCT_DEBUG
    ChFi3d_InitChron(ch); // init perf PerformTwoCornerbyInter
#endif

    done = PerformTwoCornerbyInter(Index);

#ifdef OCCT_DEBUG
    ChFi3d_ResultChron(ch, t_perf2cornerbyinter); // result perf  PerformTwoCornerbyInter
#endif

    if (!done)
    {
      PerformMoreThreeCorner(Index, 2);
      done = true;
      return;
    }
  }
  else if (c1rotule)
  { // save.
    cp11   = sd1->Vertex(isfirst1, 1);
    cp12   = sd1->Vertex(isfirst1, 2);
    cp21   = sd2->Vertex(isfirst2, 1);
    cp22   = sd2->Vertex(isfirst2, 2);
    intf11 = sd1->InterferenceOnS1();
    intf12 = sd1->InterferenceOnS2();
    intf21 = sd2->InterferenceOnS1();
    intf22 = sd2->InterferenceOnS2();
#ifdef OCCT_DEBUG
    ChFi3d_InitChron(ch); // init perf PerformTwoCornerbyInter
#endif

    done = PerformTwoCornerbyInter(Index);

#ifdef OCCT_DEBUG
    ChFi3d_ResultChron(ch, t_perf2cornerbyinter); // result perf  PerformTwoCornerbyInter
#endif
    if (!done)
    {
      // restore
      sd1->ChangeVertex(isfirst1, 1) = cp11;
      sd1->ChangeVertex(isfirst1, 2) = cp12;
      sd2->ChangeVertex(isfirst2, 1) = cp21;
      sd2->ChangeVertex(isfirst2, 2) = cp22;
      sd1->ChangeInterferenceOnS1()  = intf11;
      sd1->ChangeInterferenceOnS2()  = intf12;
      sd2->ChangeInterferenceOnS1()  = intf21;
      sd2->ChangeInterferenceOnS2()  = intf22;
      done                           = false;
    }
  }

  if (!c1biseau && !done)
  {
    // new cornerdata is created
    //-------------------------------
    occ::handle<ChFiDS_Stripe>                                        corner = new ChFiDS_Stripe();
    occ::handle<NCollection_HSequence<occ::handle<ChFiDS_SurfData>>>& cornerset =
      corner->ChangeSetOfSurfData();
    cornerset                         = new NCollection_HSequence<occ::handle<ChFiDS_SurfData>>();
    occ::handle<ChFiDS_SurfData> coin = new ChFiDS_SurfData();
    cornerset->Append(coin);

    if (SameSide)
    {
      if (ToricRotule(BRFaCo, BRS1, BRS2, st1, st2))
      {
        // Direct construction.
        // ---------------------

        int                bid;
        TopAbs_Orientation ori     = OFaCo;
        TopAbs_Orientation oriS    = st1->Orientation(IFaCo1);
        TopAbs_Orientation OFF1    = FF1.Orientation();
        TopAbs_Orientation oriSFF1 = st1->Orientation(IFaArc1);
        bid                        = 1;
        bid                        = ChFi3d::NextSide(ori, OFF1, oriS, oriSFF1, bid);
        TopAbs_Orientation op1 = TopAbs_FORWARD, op2 = TopAbs_FORWARD;
        if (yapiv)
        {
          bid = ChFi3d::ConcaveSide(BRS1, BRS2, pivot, op1, op2);
        }
        op1 = TopAbs::Reverse(op1);
        op2 = TopAbs::Reverse(op2);
#ifdef OCCT_DEBUG
        ChFi3d_InitChron(ch); // init perf ChFiKPart_ComputeData
#endif
        double radius = occ::down_cast<ChFiDS_FilSpine>(st1->Spine())->Radius();
        done          = ChFiKPart_ComputeData::ComputeCorner(DStr,
                                                    coin,
                                                    HFaCo,
                                                    HBRS1,
                                                    HBRS2,
                                                    OFaCo,
                                                    ori,
                                                    op1,
                                                    op2,
                                                    radius);
#ifdef OCCT_DEBUG
        ChFi3d_ResultChron(ch, t_chfikpartcompdata); // result perf ChFiKPart_ComputeData
#endif
      }
      else
      {
        // Construction by filling remplissage
        // ----------------------------
        if (!sameparam && Hpivot.IsNull())
        {
          // the fill between two distinct points runs along the pivot, and
          // the stripes do not end on one edge here: there is none
          PerformMoreThreeCorner(Index, 2);
          done = true;
          return;
        }
        double                         uPCArc1, uPCArc2;
        gp_Pnt2d                       p2da1, p2df1, p2da2, p2df2, p2dfac1, p2dfac2;
        gp_Vec2d                       v2dfac1, v2dfac2;
        occ::handle<GeomFill_Boundary> B1, B2, Bpiv, Bfac;
        uPCArc1 = sd1->Interference(IFaArc1).Parameter(isfirst1);
        p2da1   = sd1->Interference(IFaArc1).PCurveOnSurf()->Value(uPCArc1);
        p2df1   = sd1->Interference(IFaCo1).PCurveOnSurf()->Value(uPCArc1);
        sd1->Interference(IFaCo1).PCurveOnFace()->D1(uPCArc1, p2dfac1, v2dfac1);
        uPCArc2 = sd2->Interference(IFaArc2).Parameter(isfirst2);
        p2da2   = sd2->Interference(IFaArc2).PCurveOnSurf()->Value(uPCArc2);
        p2df2   = sd2->Interference(IFaCo2).PCurveOnSurf()->Value(uPCArc2);
        sd2->Interference(IFaCo2).PCurveOnFace()->D1(uPCArc2, p2dfac2, v2dfac2);
#ifdef OCCT_DEBUG
        ChFi3d_InitChron(ch); // init perf filling
#endif
        B1 = ChFi3d_mkbound(surf1, p2df1, p2da1, tolapp3d, 2.e-4);
        B2 = ChFi3d_mkbound(surf2, p2df2, p2da2, tolapp3d, 2.e-4);
        occ::handle<Geom2d_Curve> PCurveOnFace;
        Bfac = ChFi3d_mkbound(HFaCo,
                              PCurveOnFace,
                              Sens1,
                              p2dfac1,
                              v2dfac1,
                              Sens2,
                              p2dfac2,
                              v2dfac2,
                              tolapp3d,
                              2.e-4);
        GeomFill_ConstrainedFilling fil(8, 20);
        if (sameparam)
        {
          fil.Init(Bfac, B2, B1, true);
        }
        else
        {
          occ::handle<Adaptor3d_Curve> HPivTrim =
            Hpivot->Trim(std::min(parCP1, parCP2), std::max(parCP1, parCP2), tolesp);
          Bpiv = new GeomFill_SimpleBound(HPivTrim, tolapp3d, 2.e-4);
          fil.Init(Bfac, B2, Bpiv, B1, true);
          BRepAdaptor_Curve2d pcpivot;
          gp_Vec              dArc, dcf;
          gp_Pnt              bidon;
          Hpivot->D1(parCP1, bidon, dArc);
          double fb1, lb1;
          B1->Bounds(fb1, lb1);
          B1->D1(lb1, bidon, dcf);
          bool pivotverslebas = dArc.Dot(dcf) <= 0.;
          bool pcfalenvers    = (parCP1 > parCP2);
          if ((pivotverslebas && !pcfalenvers) || (!pivotverslebas && pcfalenvers))
          {
            FaPiv    = FF2;
            HFaPiv   = HBRS2;
            resetcp2 = true;
          }
          else
          {
            FaPiv    = FF1;
            HFaPiv   = HBRS1;
            resetcp1 = true;
          }
          FaPiv.Orientation(TopAbs_FORWARD);
          pcpivot.Initialize(pivot, FaPiv);
          TopExp_Explorer Expl;
          for (Expl.Init(FaPiv, TopAbs_EDGE); Expl.More(); Expl.Next())
          {
            if (Expl.Current().IsSame(pivot))
            {
              pctrans = Expl.Current().Orientation();
              break;
            }
          }
          if (pcpivot.GetType() != GeomAbs_BSplineCurve)
          {
            occ::handle<Geom2d_TrimmedCurve> trc =
              new Geom2d_TrimmedCurve(pcpivot.Curve(),
                                      std::min(parCP1, parCP2),
                                      std::max(parCP1, parCP2));
            PCurveOnPiv = Geom2dConvert::CurveToBSplineCurve(trc);
          }
          else
          {
            PCurveOnPiv =
              Geom2dConvert::SplitBSplineCurve(occ::down_cast<Geom2d_BSplineCurve>(pcpivot.Curve()),
                                               std::min(parCP1, parCP2),
                                               std::max(parCP1, parCP2),
                                               tol2d);
          }
          NCollection_Array1<double> kk(PCurveOnPiv->Knots());
          BSplCLib::Reparametrize(0., 1., kk);
          PCurveOnPiv->SetKnots(kk);
          if (pcfalenvers)
          {
            PCurveOnPiv->Reverse();
            pctrans = TopAbs::Reverse(pctrans);
          }
        }
        occ::handle<Geom_Surface> Surfcoin = fil.Surface();
        done                               = CompleteData(coin,
                            Surfcoin,
                            HFaCo,
                            PCurveOnFace,
                            HFaPiv,
                            PCurveOnPiv,
                            OFaCo,
                            true,
                            false,
                            false,
                            false,
                            false);
#ifdef OCCT_DEBUG
        ChFi3d_ResultChron(ch, t_remplissage); // result perf filling
#endif
      }
#ifdef OCCT_DEBUG
      ChFi3d_InitChron(ch); // init perf update DS
#endif
      if (done)
      {
        // Update 3 CornerData and the DS
        // ----------------------------------------
        if (resetcp1)
        {
          gp_Pnt pjyl   = CP1.Point();
          double tolsav = CP1.Tolerance();
          CP1.Reset();
          CP1.SetPoint(pjyl);
          CP1.SetTolerance(tolsav);
        }
        else if (resetcp2)
        {
          gp_Pnt pjyl   = CP2.Point();
          double tolsav = CP2.Tolerance();
          CP2.Reset();
          CP2.SetPoint(pjyl);
          CP2.SetTolerance(tolsav);
        }
        double                    P1deb, P2deb, P1fin, P2fin;
        int                       If1, If2, Il1, Il2, Icf, Icl;
        const ChFiDS_CommonPoint& Pf1 = coin->VertexFirstOnS1();
        ChFiDS_CommonPoint&       Pf2 = coin->ChangeVertexFirstOnS2();
        const ChFiDS_CommonPoint& Pl1 = coin->VertexLastOnS1();
        ChFiDS_CommonPoint&       Pl2 = coin->ChangeVertexLastOnS2();
        Pf2                           = CP1;
        Pl2                           = CP2;

        // the corner to start,
        // -----------------------
        ChFiDS_Regul regdeb, regfin;
        If1 = ChFi3d_IndexPointInDS(Pf1, DStr);
        If2 = ChFi3d_IndexPointInDS(Pf2, DStr);
        Il1 = ChFi3d_IndexPointInDS(Pl1, DStr);
        if (sameparam)
        {
          Il2 = If2;
        }
        else
        {
          Il2 = ChFi3d_IndexPointInDS(Pl2, DStr);
        }

        gp_Pnt2d pp1, pp2;
        pp1 =
          coin->InterferenceOnS1().PCurveOnSurf()->Value(coin->InterferenceOnS1().FirstParameter());
        pp2 =
          coin->InterferenceOnS2().PCurveOnSurf()->Value(coin->InterferenceOnS2().FirstParameter());
        occ::handle<Geom_Curve> C3d;
        double                  tolreached;
        ChFi3d_ComputeArete(Pf1,
                            pp1,
                            Pf2,
                            pp2,
                            DStr.Surface(coin->Surf()).Surface(),
                            C3d,
                            corner->ChangeFirstPCurve(),
                            P1deb,
                            P2deb,
                            tolapp3d,
                            tol2d,
                            tolreached,
                            0);
        double par1 = sd1->Interference(IFaArc1).Parameter(isfirst1);
        pp1         = sd1->Interference(IFaCo1).PCurveOnSurf()->Value(par1);
        pp2         = sd1->Interference(IFaArc1).PCurveOnSurf()->Value(par1);
        double tolr1;
        ChFi3d_ComputePCurv(C3d,
                            pp1,
                            pp2,
                            st1->ChangePCurve(isfirst1),
                            DStr.Surface(sd1->Surf()).Surface(),
                            P1deb,
                            P2deb,
                            tolapp3d,
                            tolr1);
        tolreached = std::max(tolreached, tolr1);
        TopOpeBRepDS_Curve Tcurv1(C3d, tolreached);
        Icf = DStr.AddCurve(Tcurv1);
        regdeb.SetCurve(Icf);
        regdeb.SetS1(coin->Surf(), false);
        regdeb.SetS2(sd1->Surf(), false);
        myRegul.Append(regdeb);
        corner->ChangeFirstCurve(Icf);
        corner->ChangeFirstParameters(P1deb, P2deb);
        corner->ChangeIndexFirstPointOnS1(If1);
        corner->ChangeIndexFirstPointOnS2(If2);

        pp1 =
          coin->InterferenceOnS1().PCurveOnSurf()->Value(coin->InterferenceOnS1().LastParameter());
        pp2 =
          coin->InterferenceOnS2().PCurveOnSurf()->Value(coin->InterferenceOnS2().LastParameter());
        ChFi3d_ComputeArete(Pl1,
                            pp1,
                            Pl2,
                            pp2,
                            DStr.Surface(coin->Surf()).Surface(),
                            C3d,
                            corner->ChangeLastPCurve(),
                            P1fin,
                            P2fin,
                            tolapp3d,
                            tol2d,
                            tolreached,
                            0);
        double par2 = sd2->Interference(IFaArc2).Parameter(isfirst2);
        pp1         = sd2->Interference(IFaCo2).PCurveOnSurf()->Value(par2);
        pp2         = sd2->Interference(IFaArc2).PCurveOnSurf()->Value(par2);
        double tolr2;
        ChFi3d_ComputePCurv(C3d,
                            pp1,
                            pp2,
                            st2->ChangePCurve(isfirst2),
                            DStr.Surface(sd2->Surf()).Surface(),
                            P1deb,
                            P2deb,
                            tolapp3d,
                            tolr2);
        tolreached = std::max(tolreached, tolr2);
        TopOpeBRepDS_Curve Tcurv2(C3d, tolreached);
        Icl = DStr.AddCurve(Tcurv2);
        regfin.SetCurve(Icl);
        regfin.SetS1(coin->Surf(), false);
        regfin.SetS2(sd2->Surf(), false);
        myRegul.Append(regfin);
        corner->ChangeLastCurve(Icl);
        corner->ChangeLastParameters(P1fin, P2fin);
        corner->ChangeIndexLastPointOnS1(Il1);
        corner->ChangeIndexLastPointOnS2(Il2);

        coin->ChangeIndexOfS1(DStr.AddShape(FaCo));
        if (sameparam)
        {
          coin->ChangeIndexOfS2(0);
        }
        else
        {
          coin->ChangeIndexOfS2(DStr.AddShape(FaPiv));
          coin->ChangeInterferenceOnS2().SetTransition(pctrans);
        }
        corner->SetSolidIndex(st1->SolidIndex());

        // then the starting Stripe,
        // ------------------------
        st1->SetCurve(Icf, isfirst1);
        st1->SetIndexPoint(If1, isfirst1, IFaCo1);
        st1->SetIndexPoint(If2, isfirst1, IFaArc1);
        st1->SetParameters(isfirst1, P1deb, P2deb);
        sd1->ChangeVertex(isfirst1, IFaCo1)  = Pf1;
        sd1->ChangeVertex(isfirst1, IFaArc1) = Pf2;
        sd1->ChangeInterference(IFaCo1).SetParameter(par1, isfirst1);
        if (IFaCo1 == 2)
        {
          st1->SetOrientation(TopAbs_REVERSED, isfirst1);
        }

        // then the end Stripe,
        // -------------------------
        st2->SetCurve(Icl, isfirst2);
        st2->SetIndexPoint(Il1, isfirst2, IFaCo2);
        st2->SetIndexPoint(Il2, isfirst2, IFaArc2);
        st2->SetParameters(isfirst2, P1fin, P2fin);
        sd2->ChangeVertex(isfirst2, IFaCo2)  = Pl1;
        sd2->ChangeVertex(isfirst2, IFaArc2) = Pl2;
        sd2->ChangeInterference(IFaCo2).SetParameter(par2, isfirst2);
        if (IFaCo2 == 2)
        {
          st2->SetOrientation(TopAbs_REVERSED, isfirst2);
        }
      }
#ifdef OCCT_DEBUG
      ChFi3d_ResultChron(ch, t_t2cornerDS); // result perf update DS
#endif
    }
    else
    {
      // it is necessary to make difference with
      if (!OkinterCC)
      {
        throw Standard_Failure("TwoCorner : No intersection pc pc");
      }
      occ::handle<ChFiDS_Stripe>   stsam, stdif;
      occ::handle<ChFiDS_SurfData> sdsam, sddif;
      double                       uintpcsam = 0., uintpcdif = 0.;
      int                          ifacosam = 0, ifacodif = 0, ifaopsam = 0, ifaopdif = 0;
      bool                         isfirstsam = false, isfirstdif = false;
      if (Stat1 == ChFiDS_OnSame && Stat2 == ChFiDS_OnDiff)
      {
        stsam      = st1;
        sdsam      = sd1;
        uintpcsam  = UIntPC1;
        ifacosam   = IFaCo1;
        ifaopsam   = IFaArc1;
        isfirstsam = isfirst1;
        stdif      = st2;
        sddif      = sd2;
        uintpcdif  = UIntPC2;
        ifacodif   = IFaCo2;
        ifaopdif   = IFaArc2;
        isfirstdif = isfirst2;
      }
      else if (Stat1 == ChFiDS_OnDiff && Stat2 == ChFiDS_OnSame)
      {
        stsam      = st2;
        sdsam      = sd2;
        uintpcsam  = UIntPC2;
        ifacosam   = IFaCo2;
        ifaopsam   = IFaArc2;
        isfirstsam = isfirst2;
        stdif      = st1;
        sddif      = sd1;
        uintpcdif  = UIntPC1;
        ifacodif   = IFaCo1;
        ifaopdif   = IFaArc1;
        isfirstdif = isfirst1;
      }
      else
      {
        throw Standard_Failure("TwoCorner : Config unknown");
      }
      // It is checked if surface ondiff has a point on arc from the side opposed
      // to the common face and if this arc is connected to the base face
      // opposed to common face of the surface onsame.
      ChFiDS_CommonPoint& cpopdif = sddif->ChangeVertex(isfirstdif, ifaopdif);
      if (!cpopdif.IsOnArc())
      {
        throw Standard_Failure("TwoCorner : No point on restriction on surface OnDiff");
      }
      const TopoDS_Edge& Arcopdif = cpopdif.Arc();
      const TopoDS_Face& Fopsam   = TopoDS::Face(DStr.Shape(sdsam->Index(ifaopsam)));
      TopExp_Explorer    ex;
      for (ex.Init(Fopsam, TopAbs_EDGE); ex.More(); ex.Next())
      {
        if (ex.Current().IsSame(Arcopdif))
        {
          break;
        }
        else if (!ex.More())
        {
          throw Standard_Failure("TwoCorner : No common face to loop the contour");
        }
      }
#ifdef OCCT_DEBUG
      ChFi3d_InitChron(ch); // init perf filling
#endif
      occ::handle<GeomFill_Boundary> Bsam, Bdif, Bfac;
      gp_Pnt2d ppopsam = sdsam->Interference(ifaopsam).PCurveOnSurf()->Value(uintpcsam);
      gp_Pnt2d ppcosam = sdsam->Interference(ifacosam).PCurveOnSurf()->Value(uintpcsam);
      occ::handle<Geom_Surface>        surfsam  = DStr.Surface(sdsam->Surf()).Surface();
      occ::handle<GeomAdaptor_Surface> Hsurfsam = new GeomAdaptor_Surface(surfsam);
      occ::handle<Geom2d_Curve>        pcsurfsam;
      Bsam              = ChFi3d_mkbound(Hsurfsam, pcsurfsam, ppopsam, ppcosam, tolapp3d, 2.e-4);
      double   upcopdif = sddif->Interference(ifaopdif).Parameter(isfirstdif);
      gp_Pnt2d ppopdif  = sddif->Interference(ifaopdif).PCurveOnSurf()->Value(upcopdif);
      gp_Pnt2d ppcodif  = sddif->Interference(ifacodif).PCurveOnSurf()->Value(uintpcdif);
      occ::handle<Geom_Surface>        surfdif  = DStr.Surface(sddif->Surf()).Surface();
      occ::handle<GeomAdaptor_Surface> Hsurfdif = new GeomAdaptor_Surface(surfdif);
      occ::handle<Geom2d_Curve>        pcsurfdif;
      Bdif = ChFi3d_mkbound(Hsurfdif, pcsurfdif, ppcodif, ppopdif, tolapp3d, 2.e-4);
      gp_Pnt2d ppfacsam, ppfacdif;
      gp_Pnt   PPfacsam, PPfacdif;
      gp_Vec   VVfacsam, VVfacdif;
      sdsam->Interference(ifaopsam).PCurveOnFace()->D0(uintpcsam, ppfacsam);
      const occ::handle<Geom_Curve>& curvopsam =
        DStr.Curve(sdsam->Interference(ifaopsam).LineIndex()).Curve();
      curvopsam->D1(uintpcsam, PPfacsam, VVfacsam);
      BRepAdaptor_Curve2d PCArcFac(Arcopdif, Fopsam);
      PCArcFac.D0(cpopdif.ParameterOnArc(), ppfacdif);
      // jgv for OCC26173
      BRepAdaptor_Surface SurFopsam(Fopsam);
      if (SurFopsam.IsUClosed())
      {
        double Uperiod = SurFopsam.LastUParameter() - SurFopsam.FirstUParameter();
        if (std::abs(ppfacsam.X() - ppfacdif.X()) > Uperiod / 2)
        {
          if (ppfacdif.X() < ppfacsam.X())
          {
            ppfacdif.SetX(ppfacdif.X() + Uperiod);
          }
          else
          {
            ppfacdif.SetX(ppfacdif.X() - Uperiod);
          }
        }
      }
      //////////////////
      BRepAdaptor_Curve CArcFac(Arcopdif);
      CArcFac.D1(cpopdif.ParameterOnArc(), PPfacdif, VVfacdif);
      occ::handle<BRepAdaptor_Surface> HBRFopsam = new BRepAdaptor_Surface();
      BRepAdaptor_Surface&             BRFopsam  = *HBRFopsam;
      BRFopsam.Initialize(Fopsam, false);
      occ::handle<Geom2d_Curve> pcFopsam =
        ChFi3d_BuildPCurve(HBRFopsam, ppfacsam, VVfacsam, ppfacdif, VVfacdif, true);
      Bfac = ChFi3d_mkbound(HBRFopsam, pcFopsam, tolapp3d, 2.e-4);
      GeomFill_ConstrainedFilling fil(8, 20);
      fil.Init(Bsam, Bdif, Bfac, true);
      occ::handle<Geom_Surface> Surfcoin = fil.Surface();
      TopAbs_Orientation        Osurfsam = sdsam->Orientation();
      occ::handle<Geom2d_Curve> pcnul;
      done = CompleteData(coin,
                          Surfcoin,
                          Hsurfsam,
                          pcsurfsam,
                          HBRFopsam,
                          pcnul,
                          Osurfsam,
                          true,
                          false,
                          false,
                          false,
                          false);
#ifdef OCCT_DEBUG
      ChFi3d_ResultChron(ch, t_remplissage); // result perf filling
#endif
      if (!done)
      {
        throw Standard_Failure("concavites inverted : fail");
      }
#ifdef OCCT_DEBUG
      ChFi3d_InitChron(ch); // init perf update DS
#endif
      // Update 3 CornerData and the DS
      // ----------------------------------------
      // the corner to start,
      // -----------------------
      double                    P1deb, P2deb, P1fin, P2fin;
      int                       If1, If2, Il1, Il2, Icf, Icl;
      const ChFiDS_CommonPoint& Pf1 = coin->VertexFirstOnS1();
      ChFiDS_CommonPoint&       Pf2 = coin->ChangeVertexFirstOnS2();
      const ChFiDS_CommonPoint& Pl1 = coin->VertexLastOnS1();
      ChFiDS_CommonPoint&       Pl2 = coin->ChangeVertexLastOnS2();
      Pf2 = Pl2 = cpopdif;

      ChFiDS_Regul regdeb, regfin;
      If1 = ChFi3d_IndexPointInDS(Pf1, DStr);
      If2 = ChFi3d_IndexPointInDS(Pf2, DStr);
      Il1 = ChFi3d_IndexPointInDS(Pl1, DStr);
      Il2 = If2;

      gp_Pnt2d pp1, pp2;
      pp1 =
        coin->InterferenceOnS1().PCurveOnSurf()->Value(coin->InterferenceOnS1().FirstParameter());
      pp2 =
        coin->InterferenceOnS2().PCurveOnSurf()->Value(coin->InterferenceOnS2().FirstParameter());
      occ::handle<Geom_Curve> C3d;
      double                  tolreached;
      ChFi3d_ComputeArete(Pf1,
                          pp1,
                          Pf2,
                          pp2,
                          DStr.Surface(coin->Surf()).Surface(),
                          C3d,
                          corner->ChangeFirstPCurve(),
                          P1deb,
                          P2deb,
                          tolapp3d,
                          tol2d,
                          tolreached,
                          0);
      double                         tolr1;
      occ::handle<GeomAdaptor_Curve> HC3d = new GeomAdaptor_Curve(C3d);
      ChFi3d_SameParameter(HC3d, pcFopsam, HBRFopsam, tolapp3d, tolr1);
      tolreached = std::max(tolreached, tolr1);
      // The wall kept in coplanar pieces (a body without Refine), the split
      // between two of them ending at Vtx: the arc the surface OnDiff ends
      // on bounds the piece beside the one the surface OnSame runs on, and
      // the corner's edge on the wall, drawn in the plane of both, ends on an
      // edge its face does not have.
      // Or the surface OnSame runs on the piece beyond the one at Vtx, its
      // spine's piece under it whole, the split away from Vtx: the edge
      // crosses it on its way to the arc.
      TopoDS_Face  aFopSh, aFarc;
      TopoDS_Edge  aSplit;
      bool         isAtVtx = true;
      const bool   isPiece = !cpopdif.IsVertex()
                           && PieceOfArc(Vtx,
                                         Fopsam,
                                         Arcopdif,
                                         myVEMap,
                                         myEFMap,
                                         10 * tolapp3d,
                                         aFopSh,
                                         aFarc,
                                         aSplit,
                                         isAtVtx);
      // Where that edge crosses the split, inside both: it runs from the one
      // piece into the other there.
      bool   isCross = false;
      double aWQ = 0., aParQ = 0.;
      if (isPiece)
      {
        const double            aTolX = std::max(tolreached, BRep_Tool::Tolerance(aSplit));
        const GeomAdaptor_Curve aGC(C3d, P1deb, P2deb);
        const BRepAdaptor_Curve aSC(aSplit);
        const double            aTolW = aGC.Resolution(aTolX), aTolS = aSC.Resolution(aTolX);
        Extrema_ExtCC           anExt(aGC, aSC);
        if (!anExt.IsDone() || anExt.IsParallel())
        {
          throw Standard_Failure("TwoCorner : the corner's edge and the split of its wall");
        }
        int aNbX = 0;
        for (int i = 1; i <= anExt.NbExt(); i++)
        {
          if (anExt.SquareDistance(i) > aTolX * aTolX)
          {
            continue;
          }
          Extrema_POnCurv aP1, aP2;
          anExt.Points(i, aP1, aP2);
          if (aP1.Parameter() <= P1deb + aTolW || aP1.Parameter() >= P2deb - aTolW
              || aP2.Parameter() <= aSC.FirstParameter() + aTolS
              || aP2.Parameter() >= aSC.LastParameter() - aTolS)
          {
            throw Standard_Failure("TwoCorner : the corner's edge through an end of the split");
          }
          aNbX++;
          aWQ   = aP1.Parameter();
          aParQ = aP2.Parameter();
        }
        if (aNbX > 1)
        {
          throw Standard_Failure("TwoCorner : the corner's edge crosses the split twice");
        }
        isCross = aNbX == 1;
        if (!isCross && !isAtVtx)
        {
          throw Standard_Failure("TwoCorner : the corner's edge does not cross the split");
        }
      }
      TopOpeBRepDS_Curve Tcurv1(C3d, tolreached);
      if (isCross)
      {
        Tcurv1.ChangeCurve() = new Geom_TrimmedCurve(C3d, P1deb, aWQ);
      }
      Icf = DStr.AddCurve(Tcurv1);
      // place the pcurve on face in the DS
      TopAbs_Orientation OpcFopsam = sdsam->Interference(ifaopsam).Transition();
      int                IFopsam   = sdsam->Index(ifaopsam);
      if (isfirstsam)
      {
        OpcFopsam = TopAbs::Reverse(OpcFopsam);
      }
      occ::handle<TopOpeBRepDS_SurfaceCurveInterference> interf =
        ChFi3d_FilCurveInDS(Icf, IFopsam, pcFopsam, OpcFopsam);
      DStr.ChangeShapeInterferences(IFopsam).Append(interf);

      if (isCross)
      {
        // The edge is cut where it crosses the split, as the split would
        // cut it: to there on the piece of the surface OnSame, from there on
        // the piece beside, which the arc bounds; two edges of the corner.
        // The split is cut there too, its piece from Vtx gone with Vtx.
        const gp_Pnt aPQ   = C3d->Value(aWQ);
        const double aTolQ = std::max(tolreached,
                                      aPQ.Distance(BRepAdaptor_Curve(aSplit).Value(aParQ)));
        const int    IQ    = DStr.AddPoint(TopOpeBRepDS_Point(aPQ, aTolQ));
        DStr.ChangeCurveInterferences(Icf).Append(
          ChFi3d_FilPointInDS(TopAbs_FORWARD, Icf, If1, P1deb));
        DStr.ChangeCurveInterferences(Icf).Append(
          ChFi3d_FilPointInDS(TopAbs_REVERSED, Icf, IQ, aWQ));
        const int Icf2 =
          DStr.AddCurve(TopOpeBRepDS_Curve(new Geom_TrimmedCurve(C3d, aWQ, P2deb), tolreached));
        DStr.ChangeCurveInterferences(Icf2).Append(
          ChFi3d_FilPointInDS(TopAbs_FORWARD, Icf2, IQ, aWQ));
        DStr.ChangeCurveInterferences(Icf2).Append(
          ChFi3d_FilPointInDS(TopAbs_REVERSED, Icf2, If2, P2deb));
        // on the piece beside: the same way round or, that piece the other
        // way up on its surface, the other
        const occ::handle<Geom2d_Curve> aPc2 =
          GeomProjLib::Curve2d(C3d, aWQ, P2deb, BRep_Tool::Surface(aFarc));
        if (aPc2.IsNull())
        {
          throw Standard_Failure("TwoCorner : the corner's edge over a split has no pcurve");
        }
        const int IFarc = DStr.AddShape(aFarc);
        DStr.ChangeShapeInterferences(IFarc).Append(
          ChFi3d_FilCurveInDS(Icf2,
                              IFarc,
                              aPc2,
                              aFarc.Orientation() == aFopSh.Orientation()
                                ? OpcFopsam
                                : TopAbs::Reverse(OpcFopsam)));
        // on the corner, as FILDS puts its first edge
        TopAbs_Orientation aTraCoin = TopAbs::Compose(sdsam->Orientation(), coin->Orientation());
        aTraCoin =
          TopAbs::Compose(TopAbs::Reverse(coin->InterferenceOnS1().Transition()), aTraCoin);
        DStr.ChangeSurfaceInterferences(coin->Surf())
          .Append(ChFi3d_FilCurveInDS(Icf2,
                                      coin->Surf(),
                                      corner->FirstPCurve(),
                                      TopAbs::Compose(TopAbs::Reverse(aTraCoin),
                                                      corner->FirstPCurveOrientation())));
        ChFiDS_Regul aRegul;
        aRegul.SetCurve(Icf2);
        aRegul.SetS1(coin->Surf(), false);
        aRegul.SetS2(IFarc, true);
        myRegul.Append(aRegul);
        // away from Vtx, its piece along the fillet OnSame gone under it
        bool isFromQ = true;
        if (isAtVtx)
        {
          isFromQ = BRep_Tool::Parameter(Vtx, aSplit) < aParQ;
        }
        else
        {
          const double  aWIn = uintpcsam
                               + 0.25 * (sdsam->Interference(ifaopsam).Parameter(!isfirstsam)
                                         - uintpcsam);
          const BRepAdaptor_Curve aSCin(aSplit);
          Extrema_ExtPC           anExtIn(curvopsam->Value(aWIn), aSCin);
          double        aD2 = Precision::Infinite(), aParIn = aParQ;
          for (int i = 1; anExtIn.IsDone() && i <= anExtIn.NbExt(); i++)
          {
            if (anExtIn.SquareDistance(i) < aD2)
            {
              aD2    = anExtIn.SquareDistance(i);
              aParIn = anExtIn.Point(i).Parameter();
            }
          }
          if (std::abs(aParIn - aParQ) < Precision::PConfusion())
          {
            throw Standard_Failure("TwoCorner : the side of the split under the fillet");
          }
          isFromQ = aParIn < aParQ;
        }
        const int IS = DStr.AddShape(aSplit);
        DStr.ChangeShapeInterferences(IS).Append(
          ChFi3d_FilPointInDS(isFromQ ? TopAbs_FORWARD : TopAbs_REVERSED, IS, IQ, aParQ));
      }
      else if (isPiece)
      {
        // It passes the split by, outside the two. On the wall in one face
        // the arc is cut at its end and the piece of the arc from Vtx is
        // gone, the face closing along the corner; in pieces that piece of
        // the arc stays, an edge between the two: the face of the corner's
        // edge closes along it to Vtx and the split, and the piece beside
        // keeps its outline.
        double                        aFirst, aLast;
        const occ::handle<Geom_Curve> aCarc = BRep_Tool::Curve(Arcopdif, aFirst, aLast);
        const double                  aParV = BRep_Tool::Parameter(Vtx, Arcopdif);
        const double                  aParP = cpopdif.ParameterOnArc();
        const double                  aWa = std::min(aParV, aParP), aWb = std::max(aParV, aParP);
        if (!aCarc.IsNull() && aWb - aWa > Precision::PConfusion())
        {
          TopoDS_Face aFarcF = aFarc;
          aFarcF.Orientation(TopAbs_FORWARD);
          TopAbs_Orientation anOarc = TopAbs_FORWARD;
          for (TopExp_Explorer exa(aFarcF, TopAbs_EDGE); exa.More(); exa.Next())
          {
            if (Arcopdif.IsSame(exa.Current()))
            {
              anOarc = exa.Current().Orientation();
              break;
            }
          }
          // on the piece beside as the arc lies in it; on the corner's, the
          // other side of it, the other way -- the same, that one the other
          // way up on its surface
          const TopAbs_Orientation anOop =
            aFarc.Orientation() == aFopSh.Orientation() ? TopAbs::Reverse(anOarc) : anOarc;
          const occ::handle<Geom2d_Curve> aPcArc =
            GeomProjLib::Curve2d(aCarc, aWa, aWb, BRep_Tool::Surface(aFarc));
          const occ::handle<Geom2d_Curve> aPcOp =
            GeomProjLib::Curve2d(aCarc, aWa, aWb, BRep_Tool::Surface(Fopsam));
          if (aPcArc.IsNull() || aPcOp.IsNull())
          {
            throw Standard_Failure("TwoCorner : the arc between two pieces has no pcurve");
          }
          const int Iarcp = DStr.AddCurve(
            TopOpeBRepDS_Curve(new Geom_TrimmedCurve(aCarc, aWa, aWb),
                               std::max(BRep_Tool::Tolerance(Arcopdif), Precision::Confusion())));
          const int IVtx  = DStr.AddShape(Vtx);
          const int IFarc = DStr.AddShape(aFarc);
          DStr.ChangeShapeInterferences(IFarc).Append(
            ChFi3d_FilCurveInDS(Iarcp, IFarc, aPcArc, anOarc));
          DStr.ChangeShapeInterferences(IFopsam).Append(
            ChFi3d_FilCurveInDS(Iarcp, IFopsam, aPcOp, anOop));
          if (aParV < aParP)
          {
            DStr.ChangeCurveInterferences(Iarcp).Append(
              ChFi3d_FilVertexInDS(TopAbs_FORWARD, Iarcp, IVtx, aWa));
            DStr.ChangeCurveInterferences(Iarcp).Append(
              ChFi3d_FilPointInDS(TopAbs_REVERSED, Iarcp, If2, aWb));
          }
          else
          {
            DStr.ChangeCurveInterferences(Iarcp).Append(
              ChFi3d_FilPointInDS(TopAbs_FORWARD, Iarcp, If2, aWa));
            DStr.ChangeCurveInterferences(Iarcp).Append(
              ChFi3d_FilVertexInDS(TopAbs_REVERSED, Iarcp, IVtx, aWb));
          }
          // Vtx stays, a vertex of that face, and the spine OnSame, an edge
          // of it under its fillet, is removed from there
          const TopoDS_Edge& aSpineE =
            stsam->Spine()->Edges(isfirstsam ? 1 : stsam->Spine()->NbEdges());
          TopAbs_Orientation anOV = TopAbs_FORWARD;
          for (TopExp_Explorer exv(aSpineE.Oriented(TopAbs_FORWARD), TopAbs_VERTEX);
               exv.More();
               exv.Next())
          {
            if (Vtx.IsSame(exv.Current()))
            {
              anOV = exv.Current().Orientation();
              break;
            }
          }
          const int ISpineE = DStr.AddShape(aSpineE);
          DStr.ChangeShapeInterferences(ISpineE).Append(
            ChFi3d_FilVertexInDS(TopAbs::Reverse(anOV),
                                 ISpineE,
                                 IVtx,
                                 BRep_Tool::Parameter(Vtx, aSpineE)));
        }
      }

      regdeb.SetCurve(Icf);
      regdeb.SetS1(coin->Surf(), false);
      regdeb.SetS2(IFopsam, true);
      myRegul.Append(regdeb);
      corner->ChangeFirstCurve(Icf);
      corner->ChangeFirstParameters(P1deb, P2deb);
      corner->ChangeIndexFirstPointOnS1(If1);
      corner->ChangeIndexFirstPointOnS2(If2);

      pp1 =
        coin->InterferenceOnS1().PCurveOnSurf()->Value(coin->InterferenceOnS1().LastParameter());
      pp2 =
        coin->InterferenceOnS2().PCurveOnSurf()->Value(coin->InterferenceOnS2().LastParameter());
      ChFi3d_ComputeArete(Pl1,
                          pp1,
                          Pl2,
                          pp2,
                          DStr.Surface(coin->Surf()).Surface(),
                          C3d,
                          corner->ChangeLastPCurve(),
                          P1fin,
                          P2fin,
                          tolapp3d,
                          tol2d,
                          tolreached,
                          0);
      double tolr2;
      HC3d->Load(C3d);
      ChFi3d_SameParameter(HC3d, pcsurfdif, Hsurfdif, tolapp3d, tolr2);
      tolreached = std::max(tolreached, tolr2);
      TopOpeBRepDS_Curve Tcurv2(C3d, tolreached);
      Icl = DStr.AddCurve(Tcurv2);
      regfin.SetCurve(Icl);
      regfin.SetS1(coin->Surf(), false);
      regfin.SetS2(sddif->Surf(), false);
      myRegul.Append(regfin);
      corner->ChangeLastCurve(Icl);
      corner->ChangeLastParameters(P1fin, P2fin);
      corner->ChangeIndexLastPointOnS1(Il1);
      corner->ChangeIndexLastPointOnS2(Il2);

      coin->ChangeIndexOfS1(-sdsam->Surf());
      coin->ChangeIndexOfS2(0);

      corner->SetSolidIndex(stsam->SolidIndex());

      // then Stripe OnSame
      // ---------------------
      const ChFiDS_FaceInterference& intcoin1 = coin->InterferenceOnS1();
      stsam->SetCurve(intcoin1.LineIndex(), isfirstsam);
      stsam->InDS(isfirstsam); // filDS already works from the corner.
      stsam->ChangePCurve(isfirstsam) = coin->InterferenceOnS1().PCurveOnFace();
      stsam->SetIndexPoint(If1, isfirstsam, ifaopsam);
      stsam->SetIndexPoint(Il1, isfirstsam, ifacosam);
      stsam->SetParameters(isfirstsam, intcoin1.FirstParameter(), intcoin1.LastParameter());
      sdsam->ChangeVertex(isfirstsam, ifaopsam) = Pf1;
      sdsam->ChangeVertex(isfirstsam, ifacosam) = Pl1;
      sdsam->ChangeInterferenceOnS1().SetParameter(uintpcsam, isfirstsam);
      sdsam->ChangeInterferenceOnS2().SetParameter(uintpcsam, isfirstsam);
      if (ifaopsam == 2)
      {
        stsam->SetOrientation(TopAbs_REVERSED, isfirstsam);
      }

      // then Stripe OnDiff
      // ---------------------
      stdif->SetCurve(Icl, isfirstdif);
      stdif->ChangePCurve(isfirstdif) = pcsurfdif;
      stdif->SetIndexPoint(Il2, isfirstdif, ifaopdif);
      stdif->SetIndexPoint(Il1, isfirstdif, ifacodif);
      stdif->SetParameters(isfirstdif, P1fin, P2fin);
      sddif->ChangeVertex(isfirstdif, ifaopdif) = Pl2;
      sddif->ChangeVertex(isfirstdif, ifacodif) = Pl1;
      sddif->ChangeInterference(ifacodif).SetParameter(uintpcdif, isfirstdif);
      if (ifaopdif == 1)
      {
        stdif->SetOrientation(TopAbs_REVERSED, isfirstdif);
      }
      // The fillet OnDiff made of its last two surfaces, above: its line on
      // the wall ends on the split again, and its piece beyond, to the
      // corner's point, is put on the piece at the vertex and on the fillet.
      // The arc it ended on is cut here, the line's end being the split's
      // from now on, and the fillet's end is put in the DS with it.
      if (!aMerged.Sd.IsNull() && aMerged.isShown && aMerged.Sd == sddif
          && aMerged.Side == ifaopdif)
      {
        ChFiDS_FaceInterference& aFik = sddif->ChangeInterference(ifaopdif);
        const ChFiDS_CommonPoint aCPp = sddif->Vertex(isfirstdif, ifaopdif);
        const double             aWPn = aFik.Parameter(isfirstdif);
        if (ifaopdif == 1)
        {
          sddif->ChangeIndexOfS1(aMerged.IndSave);
        }
        else
        {
          sddif->ChangeIndexOfS2(aMerged.IndSave);
        }
        aFik.ChangePCurveOnFace() = aMerged.PcSave;
        if (aMerged.isTurned)
        {
          aFik.SetTransition(TopAbs::Reverse(aFik.Transition()));
        }
        const int Isurf = sddif->Surf();
        const int indP  = stdif->IndexPoint(isfirstdif, ifaopdif);
        const int IarcP = DStr.AddShape(aCPp.Arc());
        DStr.ChangeShapeInterferences(IarcP).Append(
          ChFi3d_FilPointInDS(aCPp.TransitionOnArc(), IarcP, indP, aCPp.ParameterOnArc()));
        TopAbs_Orientation aTra = TopAbs_FORWARD;
        if (sddif->IndexOfS1() > 0)
        {
          aTra = DStr.Shape(sddif->IndexOfS1()).Orientation();
        }
        aTra = TopAbs::Compose(aTra, sddif->Orientation());
        aTra = TopAbs::Compose(TopAbs::Reverse(sddif->InterferenceOnS1().Transition()), aTra);
        DStr.ChangeSurfaceInterferences(Isurf).Append(
          ChFi3d_FilCurveInDS(stdif->Curve(isfirstdif),
                              Isurf,
                              stdif->PCurve(isfirstdif),
                              TopAbs::Compose(isfirstdif ? TopAbs::Reverse(aTra) : aTra,
                                              stdif->Orientation(isfirstdif))));
        ChFi3d_StoreLineOverSplit(DStr,
                                  stdif,
                                  sddif,
                                  isfirstdif,
                                  ifaopdif,
                                  Isurf,
                                  new GeomAdaptor_Surface(DStr.Surface(Isurf).Surface()),
                                  aMerged.Line,
                                  aMerged.Split,
                                  aMerged.Fn,
                                  aMerged.WQ,
                                  aWPn,
                                  aMerged.ParQ,
                                  indP,
                                  10 * tolapp3d,
                                  myRegul);
        stdif->InDS(isfirstdif);
        aMerged.Sd.Nullify();
      }
#ifdef OCCT_DEBUG
      ChFi3d_ResultChron(ch, t_t2cornerDS); // result perf update DS
#endif
    }
    if (!myEVIMap.IsBound(Vtx))
    {
      NCollection_List<int> li;
      myEVIMap.Bind(Vtx, li);
    }
    myEVIMap.ChangeFind(Vtx).Append(coin->Surf());
    myListStripe.Append(corner);
  }
}
