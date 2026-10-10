// Created on: 1994-03-09
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

//  Modified by skv - Mon Jun  7 18:38:57 2004 OCC5898
//  Modified by skv - Thu Aug 21 11:55:58 2008 OCC20222

#include <Adaptor2d_Curve2d.hxx>
#include <Blend_FuncInv.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAlgo_NormalProjection.hxx>
#include <BRepBlend_Line.hxx>
#include <BRepExtrema_ExtCC.hxx>
#include <BRepLib_MakeEdge.hxx>
#include <BRepTools.hxx>
#include <BRepTopAdaptor_TopolTool.hxx>
#include <ChFi3d.hxx>
#include <ChFi3d_Builder.hxx>
#include <ChFi3d_Builder_0.hxx>
#include <ChFiDS_CommonPoint.hxx>
#include <ChFiDS_FaceInterference.hxx>
#include <ChFiDS_SurfData.hxx>
#include <NCollection_Sequence.hxx>
#include <NCollection_HSequence.hxx>
#include <ChFiDS_Stripe.hxx>
#include <NCollection_List.hxx>
#include <ChFiDS_Map.hxx>
#include <ChFiDS_Spine.hxx>
#include <ElCLib.hxx>
#include <Extrema_ExtCC.hxx>
#include <Extrema_ExtPC.hxx>
#include <Extrema_ExtPS.hxx>
#include <Extrema_LocateExtCC.hxx>
#include <Extrema_POnCurv.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2d_Line.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <Geom2dInt_GInter.hxx>
#include <Geom_BezierSurface.hxx>
#include <Geom_Line.hxx>
#include <GeomAPI_IntCS.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <Geom_BoundedCurve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Curve.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_Surface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <GeomAbs_Shape.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <GeomInt_IntSS.hxx>
#include <GeomLib.hxx>
#include <GeomProjLib.hxx>
#include <gp_Lin.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Vec2d.hxx>
#include <IntCurveSurface_HInter.hxx>
#include <IntCurveSurface_IntersectionPoint.hxx>
#include <IntRes2d_IntersectionPoint.hxx>
#include <Precision.hxx>
#include <Standard_ConstructionError.hxx>
#include <Standard_Failure.hxx>
#include <Standard_NotImplemented.hxx>
#include <StdFail_NotDone.hxx>
#include <NCollection_Array1.hxx>
#include <TopAbs.hxx>
#include <TopAbs_Orientation.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopOpeBRepBuild_HBuilder.hxx>
#include <TopOpeBRepDS_Curve.hxx>
#include <TopOpeBRepDS_CurvePointInterference.hxx>
#include <TopOpeBRepDS_DataStructure.hxx>
#include <TopOpeBRepDS_HDataStructure.hxx>
#include <TopOpeBRepDS_Kind.hxx>
#include <TopOpeBRepDS_Interference.hxx>
#include <TopOpeBRepDS_Point.hxx>
#include <TopOpeBRepDS_SolidSurfaceInterference.hxx>
#include <TopOpeBRepDS_Surface.hxx>
#include <TopOpeBRepDS_SurfaceCurveInterference.hxx>
#include <TopOpeBRepDS_Transition.hxx>

#ifdef OCCT_DEBUG
//  Modified by Sergey KHROMOV - Thu Apr 11 12:23:40 2002 Begin
// The method
// ChFi3d_Builder::PerformMoreSurfdata(const int Index)
// is totally rewroted.
//  Modified by Sergey KHROMOV - Thu Apr 11 12:23:40 2002 End

extern double t_same, t_inter, t_sameinter;
extern void   ChFi3d_InitChron(OSD_Chronometer& ch);
extern void   ChFi3d_ResultChron(OSD_Chronometer& ch, double& time);
#endif
#include <Geom2dAPI_ProjectPointOnCurve.hxx>
#include <math_FunctionSample.hxx>
#include <IntRes2d_IntersectionSegment.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_BoundedSurface.hxx>

static double recadre(const double p,
                      const double ref,
                      const bool   isfirst,
                      const double first,
                      const double last)
{
  double pp = p;
  if (isfirst)
  {
    pp -= (last - first);
  }
  else
  {
    pp += (last - first);
  }
  if (std::abs(pp - ref) < std::abs(p - ref))
  {
    return pp;
  }
  return p;
}

//=======================================================================
// function : Update
// purpose  : Calculate the intersection of the face at the end of
//           the tangency line to update CommonPoint and its
//           parameter in FaceInterference.
//=======================================================================

static bool Update(const occ::handle<Adaptor3d_Surface>& fb,
                   const occ::handle<Adaptor2d_Curve2d>& pcfb,
                   const occ::handle<Adaptor3d_Surface>& surf,
                   ChFiDS_FaceInterference&              fi,
                   ChFiDS_CommonPoint&                   cp,
                   gp_Pnt2d&                             p2dbout,
                   const bool                            isfirst,
                   double&                               pared,
                   double&                               wop,
                   const double                          tol)
{
  Adaptor3d_CurveOnSurface         c1(pcfb, fb);
  occ::handle<Geom2d_Curve>        pc  = fi.PCurveOnSurf();
  occ::handle<Geom2dAdaptor_Curve> hpc = new Geom2dAdaptor_Curve(pc);
  Adaptor3d_CurveOnSurface         c2(hpc, surf);
  Extrema_LocateExtCC              ext(c1, c2, pared, wop);
  if (ext.IsDone())
  {
    double dist2 = ext.SquareDistance();
    if (dist2 < tol * tol)
    {
      Extrema_POnCurv ponc1, ponc2;
      ext.Point(ponc1, ponc2);
      double parfb = ponc1.Parameter();
      p2dbout      = pcfb->Value(parfb);
      pared        = ponc1.Parameter();
      wop          = ponc2.Parameter();
      fi.SetParameter(wop, isfirst);
      cp.Reset();
      cp.SetPoint(ponc1.Value());
      return true;
    }
  }
  return false;
}

//=======================================================================
// function : Update
// purpose  : Intersect surface <fb> and 3d curve <ct>
//           Update <isfirst> parameter of FaceInterference <fi> and point of
//           CommonPoint <cp>. Return new intersection parameters in <wop>
//           and <p2dbout>
//=======================================================================

static bool Update(const occ::handle<Adaptor3d_Surface>& fb,
                   const occ::handle<Adaptor3d_Curve>&   ct,
                   ChFiDS_FaceInterference&              fi,
                   ChFiDS_CommonPoint&                   cp,
                   gp_Pnt2d&                             p2dbout,
                   const bool                            isfirst,
                   double&                               wop)
{
  IntCurveSurface_HInter Intersection;
  // check if in KPart the limits of the tangency line
  // are already in place at this stage.
  // Modif lvt : the periodic cases are reframed, espercially if nothing was found.
  double w, uf = ct->FirstParameter(), ul = ct->LastParameter();

  double wbis = 0.;

  bool isperiodic = ct->IsPeriodic(), recadrebis = false;
  Intersection.Perform(ct, fb);
  if (Intersection.IsDone())
  {
    int    nbp = Intersection.NbPoints(), i, isol = 0, isolbis = 0;
    double dist    = Precision::Infinite();
    double distbis = Precision::Infinite();
    for (i = 1; i <= nbp; i++)
    {
      w = Intersection.Point(i).W();
      if (isperiodic)
      {
        w = recadre(w, wop, isfirst, uf, ul);
      }
      if (uf <= w && ul >= w && std::abs(w - wop) < dist)
      {
        isol = i;
        dist = std::abs(w - wop);
      }
    }
    // None on the curve: a point at an end of it, found a rounding error
    // past the end, is the end when the curve is tangent to the surface
    // there (the point is ill-conditioned along it). A crossing found past
    // the end is past the end.
    if (isol == 0 && !isperiodic)
    {
      const double tolw = Precision::PConfusion();
      for (i = 1; i <= nbp; i++)
      {
        const IntCurveSurface_IntersectionPoint& ip = Intersection.Point(i);
        w                                           = ip.W();
        if (uf - tolw > w || ul + tolw < w || std::abs(w - wop) >= dist)
        {
          continue;
        }
        gp_Pnt P;
        gp_Vec T, DU, DV;
        ct->D1(w, P, T);
        fb->D1(ip.U(), ip.V(), P, DU, DV);
        const gp_Vec N = DU.Crossed(DV);
        if (T.Magnitude() > gp::Resolution() && N.Magnitude() > gp::Resolution()
            && std::abs(T.Dot(N)) <= 1.e-6 * T.Magnitude() * N.Magnitude())
        {
          isol = i;
          dist = std::abs(w - wop);
        }
      }
    }
    if (isperiodic)
    {
      for (i = 1; i <= nbp; i++)
      {
        w = Intersection.Point(i).W();
        if (uf <= w && ul >= w && std::abs(w - wop) < distbis
            && (std::abs(w - ul) <= 0.01 || std::abs(w - uf) <= 0.01))
        {
          isolbis    = i;
          wbis       = recadre(w, wop, isfirst, uf, ul);
          distbis    = std::abs(wbis - wop);
          recadrebis = true;
        }
      }
    }
    if (isol == 0 && isolbis == 0)
    {
      return false;
    }
    if (!recadrebis)
    {
      IntCurveSurface_IntersectionPoint pint = Intersection.Point(isol);
      p2dbout.SetCoord(pint.U(), pint.V());
      w = pint.W();
      if (isperiodic)
      {
        w = ElCLib::InPeriod(w, uf, ul);
      }
      else
      {
        w = std::min(std::max(w, uf), ul);
      }
    }
    else
    {
      if (dist > distbis)
      {
        IntCurveSurface_IntersectionPoint pint = Intersection.Point(isolbis);
        p2dbout.SetCoord(pint.U(), pint.V());
        w = wbis;
      }
      else
      {
        IntCurveSurface_IntersectionPoint pint = Intersection.Point(isol);
        p2dbout.SetCoord(pint.U(), pint.V());
        w = pint.W();
        w = ElCLib::InPeriod(w, uf, ul);
      }
    }
    fi.SetParameter(w, isfirst);
    cp.Reset();
    cp.SetPoint(ct->Value(w));
    wop = w;
    return true;
  }
  return false;
}

//=======================================================================
// function : IntersUpdateOnSame
// purpose  : Intersect  curve <c3dFI> of ChFi-<Fop> interference with extended
//           surface <HBs> of <Fprol> . Return intersection parameters in
//           <FprolUV>, <c3dU> and updating <FIop> and <CPop>
//           <HGs> is a surface of ChFi
//           <Fop> is a face having 2 edges at corner with OnSame state
//           <Fprol> is a face non-adjacent to spine edge
//           <Vtx> is a corner vertex
//=======================================================================

static bool IntersUpdateOnSame(occ::handle<GeomAdaptor_Surface>& HGs,
                               occ::handle<BRepAdaptor_Surface>& HBs,
                               const occ::handle<Geom_Curve>&    c3dFI,
                               const TopoDS_Face&                Fop,
                               const TopoDS_Face&                Fprol,
                               const TopoDS_Edge&                Eprol,
                               const TopoDS_Vertex&              Vtx,
                               const bool                        isFirst,
                               const double                      Tol,
                               ChFiDS_FaceInterference&          FIop,
                               ChFiDS_CommonPoint&               CPop,
                               gp_Pnt2d&                         FprolUV,
                               double&                           c3dU)
{
  // add more or less restrictive criterions to
  // decide if the intersection is done with the face at
  // extended end or if the end is sharp.
  // A fillet whose line on <Fop> collapsed to a point (its radius that of
  // the face's curvature) has no curve there to intersect.
  if (c3dFI.IsNull())
  {
    return false;
  }
  double                         uf = FIop.FirstParameter();
  double                         ul = FIop.LastParameter();
  occ::handle<GeomAdaptor_Curve> Hc3df;
  if (c3dFI->IsPeriodic())
  {
    Hc3df = new GeomAdaptor_Curve(c3dFI);
  }
  else
  {
    Hc3df = new GeomAdaptor_Curve(c3dFI, uf, ul);
  }

  if (Update(HBs, Hc3df, FIop, CPop, FprolUV, isFirst, c3dU))
  {
    return true;
  }

  if (!ChFi3d::IsTangentFaces(Eprol, Fprol, Fop))
  {
    return false;
  }

  occ::handle<Geom2d_Curve> gpcprol = BRep_Tool::CurveOnSurface(Eprol, Fprol, uf, ul);
  if (gpcprol.IsNull())
  {
    throw Standard_ConstructionError("Failed to get p-curve of edge");
  }
  occ::handle<Geom2dAdaptor_Curve> pcprol  = new Geom2dAdaptor_Curve(gpcprol);
  double                           partemp = BRep_Tool::Parameter(Vtx, Eprol);

  return Update(HBs, pcprol, HGs, FIop, CPop, FprolUV, isFirst, partemp, c3dU, Tol);
}

//=======================================================================
// function : Update
// purpose  : Calculate the extrema curveonsurf/curveonsurf to prefer
//           the values concerning the trace on surf and the pcurve on the
//           face at end.
//=======================================================================

static bool Update(const occ::handle<Adaptor3d_Surface>& face,
                   const occ::handle<Adaptor2d_Curve2d>& edonface,
                   const occ::handle<Adaptor3d_Surface>& surf,
                   ChFiDS_FaceInterference&              fi,
                   ChFiDS_CommonPoint&                   cp,
                   const bool                            isfirst)
{
  if (!cp.IsOnArc())
  {
    return false;
  }
  Adaptor3d_CurveOnSurface  c1(edonface, face);
  double                    pared      = cp.ParameterOnArc();
  double                    parltg     = fi.Parameter(isfirst);
  occ::handle<Geom2d_Curve> pc         = fi.PCurveOnSurf();
  double                    f          = fi.FirstParameter();
  double                    l          = fi.LastParameter();
  double                    delta      = 0.1 * (l - f);
  f                                    = std::max(f - delta, pc->FirstParameter());
  l                                    = std::min(l + delta, pc->LastParameter());
  occ::handle<Geom2dAdaptor_Curve> hpc = new Geom2dAdaptor_Curve(pc, f, l);
  Adaptor3d_CurveOnSurface         c2(hpc, surf);

  Extrema_LocateExtCC ext(c1, c2, pared, parltg);
  if (ext.IsDone())
  {
    Extrema_POnCurv ponc1, ponc2;
    ext.Point(ponc1, ponc2);
    pared  = ponc1.Parameter();
    parltg = ponc2.Parameter();
    if ((parltg > f) && (parltg < l))
    {
      ////modified by jgv, 10.05.2012 for the bug 23139, 25657////
      occ::handle<Geom2d_Curve> PConF = fi.PCurveOnFace();
      if (!PConF.IsNull())
      {
        occ::handle<Geom2d_TrimmedCurve> aTrCurve = occ::down_cast<Geom2d_TrimmedCurve>(PConF);
        if (!aTrCurve.IsNull())
        {
          PConF = aTrCurve->BasisCurve();
        }
        if (!PConF->IsPeriodic())
        {
          if (isfirst)
          {
            double fpar = PConF->FirstParameter();
            if (parltg < fpar)
            {
              parltg = fpar;
            }
          }
          else
          {
            double lpar = PConF->LastParameter();
            if (parltg > lpar)
            {
              parltg = lpar;
            }
          }
        }
      }
      /////////////////////////////////////////////////////
      fi.SetParameter(parltg, isfirst);
      cp.SetArc(cp.Tolerance(), cp.Arc(), pared, cp.TransitionOnArc());
      return true;
    }
  }
  return false;
}

//=================================================================================================

static void ChFi3d_ExtendSurface(occ::handle<Geom_Surface>& S, int& prol)
{
  if (prol)
  {
    return;
  }

  prol = (S->IsKind(STANDARD_TYPE(Geom_BSplineSurface))  ? 1
          : S->IsKind(STANDARD_TYPE(Geom_BezierSurface)) ? 2
                                                         : 0);
  if (!prol)
  {
    return;
  }

  double length, umin, umax, vmin, vmax;
  gp_Pnt P1, P2;
  S->Bounds(umin, umax, vmin, vmax);
  S->D0(umin, vmin, P1);
  S->D0(umax, vmax, P2);
  length = P1.Distance(P2);

  occ::handle<Geom_BoundedSurface> aBS = occ::down_cast<Geom_BoundedSurface>(S);
  GeomLib::ExtendSurfByLength(aBS, length, 1, false, true);
  GeomLib::ExtendSurfByLength(aBS, length, 1, true, true);
  GeomLib::ExtendSurfByLength(aBS, length, 1, false, false);
  GeomLib::ExtendSurfByLength(aBS, length, 1, true, false);
  S = aBS;
}

//=======================================================================
// function : ComputeCurve2d
// purpose  : calculate the 2d of the curve Ct on face Face
//=======================================================================

static void ComputeCurve2d(const occ::handle<Geom_Curve>& Ct,
                           TopoDS_Face&                   Face,
                           occ::handle<Geom2d_Curve>&     C2d)
{
  TopoDS_Edge                                                   E1;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> MapE1;
  BRepLib_MakeEdge                                              Bedge(Ct);
  TopoDS_Edge                                                   edg = Bedge.Edge();
  BRepAlgo_NormalProjection                                     OrtProj;
  OrtProj.Init(Face);
  OrtProj.Add(edg);
  OrtProj.SetParams(1.e-6, 1.e-6, GeomAbs_C1, 14, 16);
  OrtProj.SetLimit(false);
  OrtProj.Compute3d(false);
  OrtProj.Build();
  double up1, up2;
  if (OrtProj.IsDone())
  {
    TopExp::MapShapes(OrtProj.Projection(), TopAbs_EDGE, MapE1);
    if (MapE1.Extent() != 0)
    {
      TopoDS_Shape aLocalShape = TopoDS_Shape(MapE1(1));
      E1                       = TopoDS::Edge(aLocalShape);
      //      E1=TopoDS::Edge( TopoDS_Shape (MapE1(1)));
      C2d = BRep_Tool::CurveOnSurface(E1, Face, up1, up2);
    }
  }
}

//=======================================================================
// function : PCurveInFace
// purpose  : the pcurve of <E> on <F>, for <E> as it lies in <F>. An edge
//           with two pcurves on the surface of <F> that occurs in <F> only
//           once -- a seam left on a piece of a split closed face -- would
//           otherwise give the pcurve of its own orientation, which may lie
//           a period away from the face.
//=======================================================================

static occ::handle<Geom2d_Curve> PCurveInFace(const TopoDS_Edge& E,
                                              const TopoDS_Face& F,
                                              double&            f,
                                              double&            l)
{
  if (BRep_Tool::IsClosed(E, F) && !BRepTools::IsReallyClosed(E, F))
  {
    for (TopExp_Explorer ex(F, TopAbs_EDGE); ex.More(); ex.Next())
    {
      if (ex.Current().IsSame(E))
      {
        return BRep_Tool::CurveOnSurface(TopoDS::Edge(ex.Current()), F, f, l);
      }
    }
  }
  return BRep_Tool::CurveOnSurface(E, F, f, l);
}

//=================================================================================================

static void ChFi3d_Recale(const BRepAdaptor_Surface& Bs,
                          gp_Pnt2d&                  p1,
                          gp_Pnt2d&                  p2,
                          const bool                 refon1)
{
  occ::handle<Geom_Surface>                   surf = Bs.GeomSurfaceOriginal();
  occ::handle<Geom_RectangularTrimmedSurface> ts =
    occ::down_cast<Geom_RectangularTrimmedSurface>(surf);
  if (!ts.IsNull())
  {
    surf = ts->BasisSurface();
  }
  if (surf->IsUPeriodic())
  {
    double u1 = p1.X(), u2 = p2.X();
    double uper = surf->UPeriod();
    if (fabs(u2 - u1) > 0.5 * uper)
    {
      if (u2 < u1 && refon1)
      {
        u2 += uper;
      }
      else if (u2 < u1 && !refon1)
      {
        u1 -= uper;
      }
      else if (u1 < u2 && refon1)
      {
        u2 -= uper;
      }
      else if (u1 < u2 && !refon1)
      {
        u1 += uper;
      }
    }
    p1.SetX(u1);
    p2.SetX(u2);
  }
  if (surf->IsVPeriodic())
  {
    double v1 = p1.Y(), v2 = p2.Y();
    double vper = surf->VPeriod();
    if (fabs(v2 - v1) > 0.5 * vper)
    {
      if (v2 < v1 && refon1)
      {
        v2 += vper;
      }
      else if (v2 < v1 && !refon1)
      {
        v1 -= vper;
      }
      else if (v1 < v2 && refon1)
      {
        v2 -= vper;
      }
      else if (v1 < v2 && !refon1)
      {
        v1 += vper;
      }
    }
    p1.SetY(v1);
    p2.SetY(v2);
  }
}

//=======================================================================
// function : ChFi3d_SelectStripe
// purpose  : find stripe with ChFiDS_OnSame state if <thePrepareOnSame> is True
//=======================================================================

bool ChFi3d_SelectStripe(NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator& It,
                         const TopoDS_Vertex&                                    Vtx,
                         const bool                                              thePrepareOnSame)
{
  if (!thePrepareOnSame)
  {
    return true;
  }

  for (; It.More(); It.Next())
  {
    int                        sens   = 0;
    occ::handle<ChFiDS_Stripe> stripe = It.Value();
    ChFi3d_IndexOfSurfData(Vtx, stripe, sens);
    ChFiDS_State stat;
    if (sens == 1)
    {
      stat = stripe->Spine()->FirstStatus();
    }
    else
    {
      stat = stripe->Spine()->LastStatus();
    }
    if (stat == ChFiDS_OnSame)
    {
      return true;
    }
  }

  return false;
}

static bool containV(const TopoDS_Face& F1, const TopoDS_Vertex& V);

//=======================================================================
// function : TangentNeighbour
// purpose  : The edge of <V> between <F> and another face <Fn>, across
//           which the two are tangent (one wall in two faces), or a null
//           edge. <OF> is <F>'s orientation in the shell.
//=======================================================================

static TopoDS_Edge TangentNeighbour(const TopoDS_Vertex& V,
                                    const TopoDS_Face&   F,
                                    const TopoDS_Face&   Fn,
                                    const ChFiDS_Map&    VEMap,
                                    const ChFiDS_Map&    EFMap,
                                    TopAbs_Orientation&  OF)
{
  for (NCollection_List<TopoDS_Shape>::Iterator itE(VEMap(V)); itE.More(); itE.Next())
  {
    const TopoDS_Edge& anE = TopoDS::Edge(itE.Value());
    TopoDS_Face        aF1, aF2;
    for (NCollection_List<TopoDS_Shape>::Iterator itF(EFMap(anE)); itF.More(); itF.Next())
    {
      if (F.IsSame(itF.Value()))
      {
        aF1 = TopoDS::Face(itF.Value());
      }
      else if (Fn.IsSame(itF.Value()))
      {
        aF2 = TopoDS::Face(itF.Value());
      }
    }
    if (!aF1.IsNull() && !aF2.IsNull() && ChFi3d::IsTangentFaces(anE, aF1, aF2))
    {
      OF = aF1.Orientation();
      return anE;
    }
  }
  return TopoDS_Edge();
}

static bool hasVertex(const TopoDS_Edge& E, const TopoDS_Vertex& V)
{
  TopoDS_Vertex aV1, aV2;
  TopExp::Vertices(E, aV1, aV2);
  return V.IsSame(aV1) || V.IsSame(aV2);
}

//=======================================================================
// function : EdgesToArc
// purpose  : The edges of <Fv> on the way from <V> to <Arc>, setting off
//           away from <Away> (an edge of <V>) and never along it, each
//           with its vertex toward <V>; empty when the way does not reach
//           <Arc>. When the fillet is wider than the face beside <V>
//           there, its common point lies on <Arc> beyond it, and these
//           edges lie under the fillet.
//=======================================================================

static NCollection_List<TopoDS_Shape> EdgesToArc(const TopoDS_Face&   Fv,
                                                 const TopoDS_Vertex& V,
                                                 const TopoDS_Edge&   Away,
                                                 const TopoDS_Edge&   Arc)
{
  NCollection_List<TopoDS_Shape> aWay;
  TopoDS_Vertex                  aV    = V;
  TopoDS_Edge                    aPrev = Away;
  for (int aStep = 0; aStep < 10; aStep++)
  {
    TopoDS_Edge aNext;
    for (TopExp_Explorer ex(Fv, TopAbs_EDGE); ex.More() && aNext.IsNull(); ex.Next())
    {
      const TopoDS_Edge& anE = TopoDS::Edge(ex.Current());
      if (anE.IsSame(aPrev) || BRep_Tool::Degenerated(anE))
      {
        continue;
      }
      TopoDS_Vertex aV1, aV2;
      TopExp::Vertices(anE, aV1, aV2);
      if (aV1.IsSame(aV) || aV2.IsSame(aV))
      {
        aNext = anE;
      }
    }
    if (aNext.IsNull() || aNext.IsSame(Away))
    {
      break;
    }
    if (aNext.IsSame(Arc))
    {
      return aWay;
    }
    aWay.Append(aNext);
    aWay.Append(aV);
    TopoDS_Vertex aV1, aV2;
    TopExp::Vertices(aNext, aV1, aV2);
    aV    = aV1.IsSame(aV) ? aV2 : aV1;
    aPrev = aNext;
  }
  return NCollection_List<TopoDS_Shape>();
}

//=======================================================================
// function : SplitUnderLine
// purpose  : The edge of the face <F> the fillet's line on it runs along,
//           or a null edge: the line, from the vertex <V1> to the vertex
//           <V2>, lies on an edge joining the two, the split toward the
//           piece of a wall kept in tangent pieces the spine is on (see
//           ChFi3d_SplitPieceOfSpine). The radius is then that piece's
//           width: the piece goes under the fillet, and the split with it.
//=======================================================================

static TopoDS_Edge SplitUnderLine(const TopoDS_Face&               F,
                                  const TopoDS_Vertex&             V1,
                                  const TopoDS_Vertex&             V2,
                                  const ChFiDS_FaceInterference&   Fi,
                                  const occ::handle<ChFiDS_Spine>& Spine,
                                  const ChFiDS_Map&                EFMap,
                                  const double                     Tol)
{
  if (V1.IsSame(V2) || Fi.PCurveOnFace().IsNull())
  {
    return TopoDS_Edge();
  }
  for (TopExp_Explorer ex(F, TopAbs_EDGE); ex.More(); ex.Next())
  {
    const TopoDS_Edge& anE = TopoDS::Edge(ex.Current());
    TopoDS_Vertex      aV1, aV2;
    TopExp::Vertices(anE, aV1, aV2);
    if (!((aV1.IsSame(V1) && aV2.IsSame(V2)) || (aV1.IsSame(V2) && aV2.IsSame(V1)))
        || ChFi3d_SplitPieceOfSpine(anE, F, Spine, EFMap).IsNull())
    {
      continue;
    }
    // the line's points between its ends on the edge
    BRepAdaptor_Surface aS(F, false);
    BRepAdaptor_Curve   aC(anE);
    bool                isOn = true;
    for (int i = 1; i <= 3 && isOn; i++)
    {
      const double   t   = Fi.FirstParameter() + i * (Fi.LastParameter() - Fi.FirstParameter()) / 4;
      const gp_Pnt2d aUV = Fi.PCurveOnFace()->Value(t);
      Extrema_ExtPC  anExt(aS.Value(aUV.X(), aUV.Y()), aC);
      double         aD2 = Precision::Infinite();
      for (int j = 1; anExt.IsDone() && j <= anExt.NbExt(); j++)
      {
        aD2 = std::min(aD2, anExt.SquareDistance(j));
      }
      isOn = aD2 <= Tol * Tol;
    }
    if (isOn)
    {
      return anE;
    }
  }
  return TopoDS_Edge();
}

static bool containE(const TopoDS_Face& F1, const TopoDS_Edge& E);

//=======================================================================
// function : StraightFromVertex
// purpose  : True when the edge <E> of the vertex <V> is a straight
//           segment; <L> is then its line, from <V> on away from <E>.
//=======================================================================

static bool StraightFromVertex(const TopoDS_Edge& E, const TopoDS_Vertex& V, gp_Lin& L)
{
  TopoDS_Vertex aV1, aV2;
  TopExp::Vertices(E, aV1, aV2);
  if (aV1.IsSame(aV2) || BRep_Tool::Degenerated(E) || !hasVertex(E, V))
  {
    return false;
  }
  const gp_Pnt aP = BRep_Tool::Pnt(V);
  const gp_Pnt aQ = BRep_Tool::Pnt(V.IsSame(aV1) ? aV2 : aV1);
  if (aP.Distance(aQ) <= 10 * Precision::Confusion())
  {
    return false;
  }
  L = gp_Lin(aP, gp_Dir(gp_Vec(aQ, aP)));
  BRepAdaptor_Curve aC(E);
  const double      aTol = std::max(BRep_Tool::Tolerance(E), Precision::Confusion());
  for (int i = 1; i < 8; i++)
  {
    const double t = aC.FirstParameter() + i * (aC.LastParameter() - aC.FirstParameter()) / 8;
    if (L.Distance(aC.Value(t)) > aTol)
    {
      return false;
    }
  }
  return true;
}

//=======================================================================
// function : CurveOnFacePeriod
// purpose  : The pcurve of <C> on <S>, the surface of <F> extended, on the
//           period of <F>'s domain at its vertex <V> (where <C> starts),
//           or a null curve when <C> does not lie on <S> within <Tol>.
//=======================================================================

static occ::handle<Geom2d_Curve> CurveOnFacePeriod(const occ::handle<Geom_Curve>&   C,
                                                   const occ::handle<Geom_Surface>& S,
                                                   const TopoDS_Face&               F,
                                                   const TopoDS_Vertex&             V,
                                                   const double                     Tol)
{
  occ::handle<Geom_Surface> aS = S;
  const occ::handle<Geom_RectangularTrimmedSurface> aTrimmed =
    occ::down_cast<Geom_RectangularTrimmedSurface>(aS);
  if (!aTrimmed.IsNull())
  {
    aS = aTrimmed->BasisSurface();
  }
  const double              f = C->FirstParameter(), l = C->LastParameter();
  occ::handle<Geom2d_Curve> aC2d = GeomProjLib::Curve2d(C, f, l, aS);
  if (aC2d.IsNull())
  {
    return aC2d;
  }
  for (int i = 0; i <= 4; i++)
  {
    const double   t   = f + i * (l - f) / 4;
    const gp_Pnt2d aUV = aC2d->Value(t);
    if (aS->Value(aUV.X(), aUV.Y()).Distance(C->Value(t)) > Tol)
    {
      return occ::handle<Geom2d_Curve>();
    }
  }
  const gp_Pnt2d aUVv = BRep_Tool::Parameters(V, F);
  const gp_Pnt2d aUVc = aC2d->Value(f);
  double         aDu = 0., aDv = 0.;
  if (aS->IsUPeriodic())
  {
    aDu = aS->UPeriod() * std::floor((aUVv.X() - aUVc.X()) / aS->UPeriod() + 0.5);
  }
  if (aS->IsVPeriodic())
  {
    aDv = aS->VPeriod() * std::floor((aUVv.Y() - aUVc.Y()) / aS->VPeriod() + 0.5);
  }
  if (aDu != 0. || aDv != 0.)
  {
    aC2d->Translate(gp_Vec2d(aDu, aDv));
  }
  return aC2d;
}

//=======================================================================
// function : CrossingsOnFace
// purpose  : Where the segment of <L> from <W1> to <W2> (and <Tol> past
//           its ends) meets the edges of <F> within <Tol> (or an edge's
//           tolerance): each as the edge, the parameter on <L> and the one
//           on the edge. False when the segment runs along an edge.
//=======================================================================

namespace
{
struct LineCrossing
{
  TopoDS_Edge Edge;
  double      W;
  double      Par;
};
} // namespace

static bool CrossingsOnFace(const occ::handle<Geom_Line>&   L,
                            const double                    W1,
                            const double                    W2,
                            const TopoDS_Face&              F,
                            const double                    Tol,
                            NCollection_List<LineCrossing>& Crossings)
{
  // a crossing at an end of the segment is one
  const GeomAdaptor_Curve aSeg(L, std::min(W1, W2) - Tol, std::max(W1, W2) + Tol);
  for (TopExp_Explorer ex(F, TopAbs_EDGE); ex.More(); ex.Next())
  {
    const TopoDS_Edge& anE = TopoDS::Edge(ex.Current());
    if (BRep_Tool::Degenerated(anE))
    {
      continue;
    }
    const double      aTol = std::max(Tol, BRep_Tool::Tolerance(anE));
    BRepAdaptor_Curve aC(anE);
    Extrema_ExtCC     anExt(aSeg, aC);
    if (!anExt.IsDone())
    {
      return false;
    }
    if (anExt.IsParallel())
    {
      if (anExt.SquareDistance(1) <= aTol * aTol)
      {
        return false;
      }
      continue;
    }
    for (int i = 1; i <= anExt.NbExt(); i++)
    {
      if (anExt.SquareDistance(i) > aTol * aTol)
      {
        continue;
      }
      Extrema_POnCurv aP1, aP2;
      anExt.Points(i, aP1, aP2);
      bool isKnown = false;
      for (NCollection_List<LineCrossing>::Iterator it(Crossings); it.More() && !isKnown; it.Next())
      {
        isKnown = it.Value().Edge.IsSame(anE) && std::abs(it.Value().W - aP1.Parameter()) <= aTol;
      }
      if (!isKnown)
      {
        Crossings.Append({anE, aP1.Parameter(), aP2.Parameter()});
      }
    }
  }
  return true;
}

//=======================================================================
// function : CoplanarPieces
// purpose  : True when the faces <F1> and <F2>, oriented as in the shell,
//           are pieces of one plane the same way up: a wall kept in
//           coplanar pieces (a body without Refine).
//=======================================================================

static bool CoplanarPieces(const TopoDS_Face& F1, const TopoDS_Face& F2, const double Tol)
{
  const BRepAdaptor_Surface aS1(F1, false), aS2(F2, false);
  if (aS1.GetType() != GeomAbs_Plane || aS2.GetType() != GeomAbs_Plane)
  {
    return false;
  }
  const gp_Pln aPl1 = aS1.Plane(), aPl2 = aS2.Plane();
  // the surface's normal, whichever way its axes turn
  gp_Dir aN1 = aPl1.Position().XDirection().Crossed(aPl1.Position().YDirection());
  gp_Dir aN2 = aPl2.Position().XDirection().Crossed(aPl2.Position().YDirection());
  if (F1.Orientation() == TopAbs_REVERSED)
  {
    aN1.Reverse();
  }
  if (F2.Orientation() == TopAbs_REVERSED)
  {
    aN2.Reverse();
  }
  return aN1.Dot(aN2) >= 1. - 1.e-9 && aPl1.Distance(aPl2.Location()) <= Tol;
}

//=======================================================================
// function : LineOverSplit
// purpose  : The fillet's line on the plane face <Fop>, the straight line
//           <L3d> from its end at <W0> (the walk stopped at the spine's
//           end, in the face), carried on away from the fillet (back when
//           <isFirst>) to the face at end <Fv>: it crosses one edge of
//           <Fop>, <Esplit>, into <Fn> -- a piece of the same plane, the
//           same way up: a wall kept in coplanar pieces -- and meets <Fv>
//           on an edge of <Fn> and <Fv>, <Earc>, crossing no other. <L> is
//           the line, <WQ>, <WP> its parameters at <Esplit> and <Earc>,
//           <ParQ>, <ParP> theirs. The line then ends on <Fv> as it would
//           on the wall in one face, and runs over both pieces.
//=======================================================================

static bool LineOverSplit(const occ::handle<Geom_Curve>& L3d,
                          const double                   W0,
                          const bool                     isFirst,
                          const TopoDS_Face&             Fop,
                          const TopoDS_Face&             Fv,
                          const ChFiDS_Map&              EFMap,
                          const double                   Tol,
                          occ::handle<Geom_Line>&        L,
                          TopoDS_Edge&                   Esplit,
                          TopoDS_Face&                   Fn,
                          TopoDS_Edge&                   Earc,
                          double&                        WQ,
                          double&                        WP,
                          double&                        ParQ,
                          double&                        ParP)
{
  if (L3d.IsNull())
  {
    return false;
  }
  const GeomAdaptor_Curve aC(L3d);
  if (aC.GetType() != GeomAbs_Line)
  {
    return false;
  }
  L = new Geom_Line(aC.Line());
  if (L->Value(W0).Distance(L3d->Value(W0)) > Tol)
  {
    return false;
  }
  const BRepAdaptor_Surface aSop(Fop, false);
  if (aSop.GetType() != GeomAbs_Plane)
  {
    return false;
  }
  // where the line meets Fv's surface, the nearest past its end
  const double              aSense = isFirst ? -1. : 1.;
  occ::handle<Geom_Surface> aSv    = BRep_Tool::Surface(Fv);
  GeomAPI_IntCS             anInt(L, aSv);
  bool                      isFound = false;
  for (int i = 1; anInt.IsDone() && i <= anInt.NbPoints(); i++)
  {
    double u, v, w;
    anInt.Parameters(i, u, v, w);
    if ((w - W0) * aSense > Tol && (!isFound || (w - WP) * aSense < 0.))
    {
      WP      = w;
      isFound = true;
    }
  }
  if (!isFound)
  {
    return false;
  }
  // on the way there, on Fop: one crossing, an edge into a coplanar piece
  NCollection_List<LineCrossing> aOnOp;
  if (!CrossingsOnFace(L, W0 + 2. * aSense * Tol, WP, Fop, Tol, aOnOp) || aOnOp.Extent() != 1)
  {
    return false;
  }
  Esplit = aOnOp.First().Edge;
  WQ     = aOnOp.First().W;
  ParQ   = aOnOp.First().Par;
  if (std::abs(WP - WQ) <= Tol || EFMap(Esplit).Extent() != 2)
  {
    return false;
  }
  for (NCollection_List<TopoDS_Shape>::Iterator it(EFMap(Esplit)); it.More(); it.Next())
  {
    if (!Fop.IsSame(it.Value()))
    {
      Fn = TopoDS::Face(it.Value());
    }
  }
  if (Fn.IsNull() || Fn.IsSame(Fv))
  {
    return false;
  }
  if (!CoplanarPieces(Fop, Fn, Tol))
  {
    return false;
  }
  // on Fn: from Esplit to an edge Fn shares with Fv, nothing between
  NCollection_List<LineCrossing> aOnN;
  if (!CrossingsOnFace(L, WQ, WP, Fn, Tol, aOnN))
  {
    return false;
  }
  const double aTolW = 10. * Tol;
  for (NCollection_List<LineCrossing>::Iterator it(aOnN); it.More(); it.Next())
  {
    const LineCrossing& aX = it.Value();
    if (aX.Edge.IsSame(Esplit) && std::abs(aX.W - WQ) <= aTolW)
    {
      continue;
    }
    if (Earc.IsNull() && std::abs(aX.W - WP) <= aTolW && containE(Fv, aX.Edge))
    {
      Earc = aX.Edge;
      ParP = aX.Par;
      continue;
    }
    return false;
  }
  if (Earc.IsNull())
  {
    return false;
  }
  // both points inside their edges
  for (int k = 0; k < 2; k++)
  {
    const TopoDS_Edge& anE = k == 0 ? Esplit : Earc;
    const double       aP  = k == 0 ? ParQ : ParP;
    double             f, l;
    BRep_Tool::Range(anE, f, l);
    const double aTolP =
      BRepAdaptor_Curve(anE).Resolution(std::max(Tol, BRep_Tool::Tolerance(anE)));
    if (aP <= f + aTolP || aP >= l - aTolP)
    {
      return false;
    }
  }
  return true;
}

//=======================================================================
// function : TransitionOnEdge
// purpose  : The transition of the fillet's line with the interference
//           <Fi> on the face <Fi>'s face leaving (arriving through, when
//           <isFirst>) the face <F> across its edge <E>, as the walk gives a
//           line ending on an arc. <F> is <Fi>'s face or a piece of the same
//           plane beside it, the same way up in the shell: one the other way
//           round on its surface (<isReversed>) has the line the other way.
//=======================================================================

static TopAbs_Orientation TransitionOnEdge(const ChFiDS_FaceInterference& Fi,
                                           const TopoDS_Face&             F,
                                           const TopoDS_Edge&             E,
                                           const bool                     isFirst,
                                           const bool                     isReversed = false)
{
  TopAbs_Orientation anOr = E.Orientation();
  for (TopExp_Explorer ex(F.Oriented(TopAbs_FORWARD), TopAbs_EDGE); ex.More(); ex.Next())
  {
    if (E.IsSame(ex.Current()))
    {
      anOr = ex.Current().Orientation();
      break;
    }
  }
  const TopAbs_Orientation aTr =
    TopAbs::Compose(isReversed ? TopAbs::Reverse(Fi.Transition()) : Fi.Transition(), anOr);
  return isFirst ? TopAbs::Reverse(aTr) : aTr;
}

//=======================================================================
// function : StoreLineOverSplit
// purpose  : The fillet's line on side <Ons> of <Fd>, carried over a wall
//           in coplanar pieces (LineOverSplit) to P (<IndP>, at <WP> on
//           <L>) on the far piece <Fn>: the line is made to end at Q (at
//           <WQ>), on the split <Esplit>, where FILDS ends it, and its piece
//           from Q to P is put in the DS on <Fn> and on the fillet's
//           surface <Isurf> (<HGs>), as FILDS puts a line.
//=======================================================================

static void StoreLineOverSplit(TopOpeBRepDS_DataStructure&             DStr,
                               const occ::handle<ChFiDS_Stripe>&       Stripe,
                               const occ::handle<ChFiDS_SurfData>&     Fd,
                               const bool                              isFirst,
                               const int                               Ons,
                               const int                               Isurf,
                               const occ::handle<GeomAdaptor_Surface>& HGs,
                               const occ::handle<Geom_Line>&           L,
                               const TopoDS_Edge&                      Esplit,
                               const TopoDS_Face&                      Fn,
                               const double                            WQ,
                               const double                            WP,
                               const double                            ParQ,
                               const int                               IndP,
                               const double                            Tol,
                               NCollection_List<ChFiDS_Regul>&         Regul)
{
  // on the fillet's surface, as FILDS puts the lines
  TopAbs_Orientation aTrafil1 = TopAbs_FORWARD;
  if (Fd->IndexOfS1() > 0)
  {
    aTrafil1 = DStr.Shape(Fd->IndexOfS1()).Orientation();
  }
  aTrafil1 = TopAbs::Compose(aTrafil1, Fd->Orientation());
  aTrafil1 = TopAbs::Compose(TopAbs::Reverse(Fd->InterferenceOnS1().Transition()), aTrafil1);

  ChFiDS_FaceInterference& aFi  = Fd->ChangeInterference(Ons);
  const TopoDS_Face        aFOp = TopoDS::Face(DStr.Shape(Fd->Index(Ons)));
  // the far piece the other way round on its surface: the line too
  const bool isRevN = Fn.Orientation() != aFOp.Orientation();

  // the line ends at Q on the split, where FILDS takes it; Q is where
  // the line crosses the edge, within the gap between them there
  ChFiDS_CommonPoint& aCPOp = Fd->ChangeVertex(isFirst, Ons);
  const double        aTolQ =
    std::max({aCPOp.Tolerance(),
              BRepAdaptor_Curve(Esplit).Value(ParQ).Distance(L->Value(WQ)),
              Precision::Confusion()});
  aCPOp.Reset();
  aCPOp.SetPoint(L->Value(WQ));
  aCPOp.SetArc(aTolQ, Esplit, ParQ, TransitionOnEdge(aFi, aFOp, Esplit, isFirst));
  aFi.SetParameter(WQ, isFirst);
  const int indQ = ChFi3d_IndexPointInDS(aCPOp, DStr);
  Stripe->SetIndexPoint(indQ, isFirst, Ons);

  // the piece from Q to P, the way the line runs, on the far piece and on
  // the fillet's surface
  const double                  aWa  = std::min(WQ, WP), aWb = std::max(WQ, WP);
  const occ::handle<Geom_Curve> aCQP = new Geom_TrimmedCurve(L, aWa, aWb);
  occ::handle<Geom2d_Curve>     aPcN =
    CurveOnFacePeriod(aCQP, BRep_Tool::Surface(Fn), Fn, TopExp::FirstVertex(Esplit), Tol);
  occ::handle<Geom2d_Curve>                aPsS  = aFi.PCurveOnSurf();
  const occ::handle<Geom2d_TrimmedCurve> aPsTr = occ::down_cast<Geom2d_TrimmedCurve>(aPsS);
  if (!aPsTr.IsNull())
  {
    aPsS = aPsTr->BasisCurve();
  }
  if (aPcN.IsNull() || aPsS.IsNull())
  {
    throw Standard_Failure("ChFi3d : the line over a split wall has no pcurve");
  }
  aPsS = new Geom2d_TrimmedCurve(aPsS, aWa, aWb);
  const occ::handle<GeomAdaptor_Surface> aHSn = new GeomAdaptor_Surface(BRep_Tool::Surface(Fn));
  const double                           aTolQP =
    std::max(ChFi3d_EvalTolReached(HGs, aPsS, aHSn, aPcN, aCQP), Precision::Confusion());
  if (aTolQP > Tol)
  {
    throw Standard_Failure("ChFi3d : the line over a split wall is off the fillet");
  }
  const int IQP = DStr.AddCurve(TopOpeBRepDS_Curve(aCQP, aTolQP));
  const int IFn = DStr.AddShape(Fn);
  DStr.ChangeShapeInterferences(IFn).Append(
    ChFi3d_FilCurveInDS(IQP,
                        IFn,
                        aPcN,
                        isRevN ? TopAbs::Reverse(aFi.Transition()) : aFi.Transition()));
  DStr.ChangeSurfaceInterferences(Isurf).Append(
    ChFi3d_FilCurveInDS(IQP, Isurf, aPsS, Ons == 1 ? aTrafil1 : TopAbs::Reverse(aTrafil1)));
  DStr.ChangeCurveInterferences(IQP).Append(
    ChFi3d_FilPointInDS(TopAbs_FORWARD, IQP, isFirst ? IndP : indQ, aWa));
  DStr.ChangeCurveInterferences(IQP).Append(
    ChFi3d_FilPointInDS(TopAbs_REVERSED, IQP, isFirst ? indQ : IndP, aWb));
  ChFiDS_Regul aRegul;
  aRegul.SetCurve(IQP);
  aRegul.SetS1(Isurf, false);
  aRegul.SetS2(IFn, true);
  Regul.Append(aRegul);
}

//=======================================================================
// function : KeepCurveOnFace
// purpose  : The corner's curve on <F> (<C>, from <First> to <Last>) can
//           end tangent to an iso-line bounding <F> -- a fillet of the
//           radius of the cylinder it ends on leaves the cylinder's circle
//           at fourth order -- and its approximation then runs a rounding
//           error past the line near the end: the face's wire crosses
//           itself there. The curve is refined towards that end (simple
//           knots, the curve as smooth as it was) until its poles close in
//           on it, the poles past the line are put on it and the one next
//           to the end a margin inside, so the curve leaves the line at an
//           angle: within the hull of its poles, it is inside but at its
//           end. Left alone when that moves the curve more than
//           Precision::Confusion(). Returns how far it moved the curve in
//           3D, at most (0 when left alone).
//=======================================================================

static double KeepCurveOnFace(occ::handle<Geom2d_Curve>& C,
                              const TopoDS_Face&         F,
                              const double               First,
                              const double               Last)
{
  occ::handle<Geom2d_BSplineCurve> aBS = occ::down_cast<Geom2d_BSplineCurve>(C);
  if (aBS.IsNull() || F.IsNull())
  {
    return 0.;
  }
  BRepAdaptor_Surface aS(F, false);
  double              aU1, aU2, aV1, aV2;
  BRepTools::UVBounds(F, aU1, aU2, aV1, aV2);
  const double aMaxMove = Precision::Confusion();
  const double aMargin  = 0.1 * Precision::Confusion();
  double       aMoved   = 0.;
  for (int anEnd = 0; anEnd < 2; ++anEnd)
  {
    const bool   isLast = anEnd == 1;
    const double aT     = isLast ? Last : First;
    if (std::abs(aT - (isLast ? aBS->LastParameter() : aBS->FirstParameter()))
        > Precision::PConfusion())
    {
      continue;
    }
    const gp_Pnt2d aEnd = aBS->Value(aT);
    for (int aBound = 0; aBound < 4; ++aBound)
    {
      const bool onU = aBound < 2;
      if (onU ? aS.IsUPeriodic() : aS.IsVPeriodic())
      {
        continue;
      }
      // parameter per unit of length across the line
      const double aRes  = onU ? aS.UResolution(1.) : aS.VResolution(1.);
      const double aLine = aBound == 0 ? aU1 : aBound == 1 ? aU2 : aBound == 2 ? aV1 : aV2;
      const double aSide = aBound % 2 == 0 ? -1. : 1.;
      auto         past  = [&](const gp_Pnt2d& P) {
        return aSide * ((onU ? P.X() : P.Y()) - aLine);
      };
      if (std::abs(past(aEnd)) > aRes * Precision::Confusion())
      {
        continue;
      }
      bool isPast = false;
      for (int i = 1; i <= aBS->NbPoles() && !isPast; ++i)
      {
        isPast = past(aBS->Pole(i)) > 0.;
      }
      if (!isPast)
      {
        continue;
      }
      // simple knots nearer the end each time: the poles close in on the
      // curve there, and the curve keeps its continuity. The refinement
      // that moves the curve least, once within twice the margin.
      const double                     aRange = aBS->LastParameter() - aBS->FirstParameter();
      occ::handle<Geom2d_BSplineCurve> aRef   = occ::down_cast<Geom2d_BSplineCurve>(aBS->Copy());
      occ::handle<Geom2d_BSplineCurve> aBest;
      double                           aBestMove = aMaxMove;
      for (int k = 1; k <= 12 && aBestMove > 2. * aMargin; ++k)
      {
        // a knot there already stays as it is
        NCollection_Array1<double> aKnot(1, 1);
        NCollection_Array1<int>    aMult(1, 1);
        aKnot(1) = isLast ? aT - std::ldexp(aRange, -k) : aT + std::ldexp(aRange, -k);
        aMult(1) = 1;
        aRef->InsertKnots(aKnot, aMult, 0., false);
        const int aNb   = aRef->NbPoles();
        const int aNext = isLast ? aNb - 1 : 2;
        // how far past the line, or short of the margin next to the end
        auto   over  = [&](const int i) {
          return past(aRef->Pole(i)) + (i == aNext ? aMargin * aRes : 0.);
        };
        double aMove = 0.;
        for (int i = 1; i <= aNb; ++i)
        {
          aMove = std::max(aMove, over(i) / aRes);
        }
        if (aMove > aBestMove)
        {
          continue;
        }
        aBest     = occ::down_cast<Geom2d_BSplineCurve>(aRef->Copy());
        aBestMove = aMove;
        for (int i = 1; i <= aNb; ++i)
        {
          const double anOver = over(i);
          if (anOver > 0.)
          {
            gp_Pnt2d aP = aBest->Pole(i);
            if (onU)
            {
              aP.SetX(aP.X() - aSide * anOver);
            }
            else
            {
              aP.SetY(aP.Y() - aSide * anOver);
            }
            aBest->SetPole(i, aP);
          }
        }
      }
      if (!aBest.IsNull())
      {
        aBS    = aBest;
        aMoved = std::max(aMoved, aBestMove);
      }
    }
  }
  C = aBS;
  return aMoved;
}

//=======================================================================
// function : PerformOneCorner
// purpose  : Calculate a corner with three edges and a fillet.
//           3 separate case: (22/07/94 only 1st is implemented)
//
//           - same concavity on three edges, intersection with the
//             face at end,
//           - concavity of 2 outgoing edges is opposite to the one of the fillet,
//             if the face at end is ready for that, the same in  case 1 on extended face,
//             otherwise a small cap is done with GeomFill,
//           - only one outgoing edge has concavity opposed to the edge of the
//             fillet and the third edge, the top of the corner is reread
//             in the empty of the fillet and closed, either by extending the face
//             at end if it is plane and orthogonal to the
//             guiding edge, or by a cap of type GeomFill.
//
//           <thePrepareOnSame> means that only needed thing is redefinition
//           of intersection pameter of OnSame-Stripe with <Arcprol>
//           (eap, Arp 9 2002, occ266)
//=======================================================================

void ChFi3d_Builder::PerformOneCorner(const int Index, const bool thePrepareOnSame)
{
  TopOpeBRepDS_DataStructure& DStr = myDS->ChangeDS();

#ifdef OCCT_DEBUG
  OSD_Chronometer ch; // init perf for PerformSetOfKPart
#endif
  // the top,
  const TopoDS_Vertex& Vtx = myVDataMap.FindKey(Index);
  // The fillet is returned,
  NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator StrIt;
  StrIt.Initialize(myVDataMap(Index));
  if (!ChFi3d_SelectStripe(StrIt, Vtx, thePrepareOnSame))
  {
    return;
  }
  occ::handle<ChFiDS_Stripe>                          stripe = StrIt.Value();
  const occ::handle<ChFiDS_Spine>                     spine  = stripe->Spine();
  NCollection_Sequence<occ::handle<ChFiDS_SurfData>>& SeqFil =
    stripe->ChangeSetOfSurfData()->ChangeSequence();
  // SurfData and its CommonPoints,
  int sens = 0;

  // Choose proper SurfData
  int  num     = ChFi3d_IndexOfSurfData(Vtx, stripe, sens);
  bool isfirst = (sens == 1);
  if (isfirst)
  {
    for (; num < SeqFil.Length()
           && ((SeqFil.Value(num)->IndexOfS1() == 0) || (SeqFil.Value(num)->IndexOfS2() == 0));)
    {
      SeqFil.Remove(num); // The surplus is removed
    }
  }
  else
  {
    for (; num > 1
           && ((SeqFil.Value(num)->IndexOfS1() == 0) || (SeqFil.Value(num)->IndexOfS2() == 0));)
    {
      SeqFil.Remove(num); // The surplus is removed
      num--;
    }
  }

  occ::handle<ChFiDS_SurfData>& Fd  = SeqFil.ChangeValue(num);
  ChFiDS_CommonPoint&           CV1 = Fd->ChangeVertex(isfirst, 1);
  ChFiDS_CommonPoint&           CV2 = Fd->ChangeVertex(isfirst, 2);
  // To evaluate the new points.
  Bnd_Box box1, box2;

  // The cases of cap and intersection are processed separately.
  // ----------------------------------------------------------
  ChFiDS_State stat;
  if (isfirst)
  {
    stat = spine->FirstStatus();
  }
  else
  {
    stat = spine->LastStatus();
  }
  bool        onsame = (stat == ChFiDS_OnSame);
  TopoDS_Face Fv, Fad, Fop;
  TopoDS_Edge Arcpiv, Arcprol, Arcspine;
  // edges of Fv from Vtx to Arcpiv when Arcpiv is not an edge of Vtx
  NCollection_List<TopoDS_Shape> Swallowed;
  // Arcprol in a face tangent to Fop across an edge of Vtx (one wall in two
  // faces): that face, the edge, and whether the extension runs along it
  TopoDS_Face Fprol;
  TopoDS_Edge Etan;
  bool        zobOnEtan = false;
  // The cut over a wall in two tangent faces (see below): the face beyond
  // the one at Vtx, the straight edge between the two from Arcpiv's vertex
  // Vp, its line carried on past Vp, the edge of Fad from Vp to Vtx, and
  // the two faces' orientations in the shell
  TopoDS_Face        FvT;
  TopoDS_Edge        Etg, Eunder;
  TopoDS_Vertex      Vp;
  gp_Lin             LinTg;
  TopAbs_Orientation OFvT = TopAbs_FORWARD, OFv = TopAbs_FORWARD;
  // The fillet's line over a wall in coplanar pieces (see LineOverSplit):
  // its side, the line, the edge it crosses and the one it ends on, the
  // piece beyond, and the parameters of the two points
  int                    onsOS = 0;
  occ::handle<Geom_Line> LinOS;
  TopoDS_Edge            EsplitOS, EarcOS;
  TopoDS_Face            FnOS;
  double                 wQOS = 0., wPOS = 0., parQOS = 0., parPOS = 0.;
  if (isfirst)
  {
    Arcspine = spine->Edges(1);
  }
  else
  {
    Arcspine = spine->Edges(spine->NbEdges());
  }
  TopAbs_Orientation               OArcprolv = TopAbs_FORWARD, OArcprolop = TopAbs_FORWARD;
  int                              ICurve;
  occ::handle<BRepAdaptor_Surface> HBs  = new BRepAdaptor_Surface();
  occ::handle<BRepAdaptor_Surface> HBad = new BRepAdaptor_Surface();
  occ::handle<BRepAdaptor_Surface> HBop = new BRepAdaptor_Surface();
  occ::handle<BRepAdaptor_Surface> HBsT = new BRepAdaptor_Surface();
  BRepAdaptor_Surface&             Bs   = *HBs;
  BRepAdaptor_Surface&             Bad  = *HBad;
  BRepAdaptor_Surface&             Bop  = *HBop;
  occ::handle<Geom_Curve>          Cc;
  occ::handle<Geom2d_Curve>        Pc, Ps;
  double                           Ubid, Vbid; //,mu,Mu,mv,Mv;
  double                           Udeb = 0., Ufin = 0.;
  //  gp_Pnt2d UVf1,UVl1,UVf2,UVl2;
  //  double Du,Dv,Step;
  bool inters  = true;
  int  IFadArc = 1, IFopArc = 2;
  Fop = TopoDS::Face(DStr.Shape(Fd->Index(IFopArc)));
  TopExp_Explorer ex;

#ifdef OCCT_DEBUG
  ChFi3d_InitChron(ch); // init perf condition  if (onsame)
#endif

  if (onsame)
  {
    if (!CV1.IsOnArc() && !CV2.IsOnArc())
    {
      throw Standard_ConstructionError("Corner OnSame : no point on arc");
    }
    else if (CV1.IsOnArc() && CV2.IsOnArc())
    {
      bool sur1 = false, sur2 = false;
      for (ex.Init(CV1.Arc(), TopAbs_VERTEX); ex.More(); ex.Next())
      {
        if (Vtx.IsSame(ex.Current()))
        {
          sur1 = true;
          break;
        }
      }
      for (ex.Init(CV2.Arc(), TopAbs_VERTEX); ex.More(); ex.Next())
      {
        if (Vtx.IsSame(ex.Current()))
        {
          sur2 = true;
          break;
        }
      }
      // An arc of a wall kept in tangent pieces reaches Vtx along its other
      // face through the pieces: it is on Vtx as the whole wall's would be,
      // and the edge of Vtx it continues as stands for it.
      TopoDS_Edge EV1 = CV1.Arc(), EV2 = CV2.Arc();
      for (int ons = 1; ons <= 2 && sur1 != sur2; ons++)
      {
        if (ons == 1 ? sur1 : sur2)
        {
          continue;
        }
        const TopoDS_Edge& anArc = ons == 1 ? CV1.Arc() : CV2.Arc();
        const TopoDS_Face  aF    = TopoDS::Face(DStr.Shape(Fd->Index(ons)));
        TopoDS_Face        anOther;
        for (NCollection_List<TopoDS_Shape>::Iterator itF(myEFMap(anArc)); itF.More(); itF.Next())
        {
          if (!aF.IsSame(itF.Value()))
          {
            anOther = TopoDS::Face(itF.Value());
          }
        }
        if (aF.IsNull() || anOther.IsNull())
        {
          continue;
        }
        const TopoDS_Edge anEV =
          ChFi3d_EdgeOnSplitToVertex(anArc, aF, anOther, Vtx, myEFMap, myVEMap);
        if (!anEV.IsNull())
        {
          (ons == 1 ? EV1 : EV2)   = anEV;
          (ons == 1 ? sur1 : sur2) = true;
        }
      }
      if (sur1 && sur2)
      {
        TopoDS_Edge E[3];
        E[0] = EV1;
        E[1] = EV2;
        E[2] = Arcspine;
        if (ChFi3d_EdgeState(E, myEFMap) != ChFiDS_OnDiff)
        {
          IFadArc = 2;
        }
      }
      else if (sur2)
      {
        IFadArc = 2;
      }
    }
    else if (CV2.IsOnArc())
    {
      IFadArc = 2;
    }
    IFopArc = 3 - IFadArc;

    Arcpiv = Fd->Vertex(isfirst, IFadArc).Arc();
    Fad    = TopoDS::Face(DStr.Shape(Fd->Index(IFadArc)));
    Fop    = TopoDS::Face(DStr.Shape(Fd->Index(IFopArc)));
    NCollection_List<TopoDS_Shape>::Iterator It;
    // The face at end is returned without check of its unicity.
    for (It.Initialize(myEFMap(Arcpiv)); It.More(); It.Next())
    {
      if (!Fad.IsSame(It.Value()))
      {
        Fv = TopoDS::Face(It.Value());
        break;
      }
    }

    // Does the face at bout contain the Vertex ?
    bool isinface = false;
    for (ex.Init(Fv, TopAbs_VERTEX); ex.More(); ex.Next())
    {
      if (ex.Current().IsSame(Vtx))
      {
        isinface = true;
        break;
      }
    }
    if (!isinface && Fd->Vertex(isfirst, 3 - IFadArc).IsOnArc())
    {
      IFadArc = 3 - IFadArc;
      IFopArc = 3 - IFopArc;
      Arcpiv  = Fd->Vertex(isfirst, IFadArc).Arc();
      Fad     = TopoDS::Face(DStr.Shape(Fd->Index(IFadArc)));
      Fop     = TopoDS::Face(DStr.Shape(Fd->Index(IFopArc)));
      // NCollection_List<TopoDS_Shape>::Iterator It;
      // The face at end is returned without check of its unicity.
      for (It.Initialize(myEFMap(Arcpiv)); It.More(); It.Next())
      {
        if (!Fad.IsSame(It.Value()))
        {
          Fv = TopoDS::Face(It.Value());
          break;
        }
      }
    }

    if (Fv.IsNull())
    {
      throw StdFail_NotDone("OneCorner : face at end not found");
    }

    // The fillet's line on Fad passes the end of the edge between Fad and
    // the face at Vtx (a drafted wall's rounded corner, a cone) onto an
    // edge of the face beyond it (the wall's plane), into which the face at
    // Vtx runs on, tangent, across a straight edge from Arcpiv's vertex:
    // one wall in two faces, and Fv does not hold Vtx. The cut runs over
    // both, from Arcpiv on the face beyond (FvT) across that edge carried
    // on past the vertex, to Fop on the face at Vtx, which is Fv from here.
    if (!containV(Fv, Vtx))
    {
      for (ex.Init(Arcpiv, TopAbs_VERTEX); ex.More() && FvT.IsNull(); ex.Next())
      {
        const TopoDS_Vertex& aVp = TopoDS::Vertex(ex.Current());
        for (It.Initialize(myVEMap(aVp)); It.More() && FvT.IsNull(); It.Next())
        {
          const TopoDS_Edge& anE = TopoDS::Edge(It.Value());
          if (anE.IsSame(Arcpiv) || !hasVertex(anE, Vtx) || !containE(Fad, anE))
          {
            continue;
          }
          TopoDS_Face aFV;
          for (NCollection_List<TopoDS_Shape>::Iterator itF(myEFMap(anE)); itF.More(); itF.Next())
          {
            if (!Fad.IsSame(itF.Value()))
            {
              aFV = TopoDS::Face(itF.Value());
            }
          }
          if (aFV.IsNull() || aFV.IsSame(Fv) || aFV.IsSame(Fop))
          {
            continue;
          }
          TopAbs_Orientation anOT = TopAbs_FORWARD, anOV = TopAbs_FORWARD;
          const TopoDS_Edge  anEtg = TangentNeighbour(aVp, Fv, aFV, myVEMap, myEFMap, anOT);
          gp_Lin             aLin;
          if (anEtg.IsNull() || !StraightFromVertex(anEtg, aVp, aLin))
          {
            continue;
          }
          TangentNeighbour(aVp, aFV, Fv, myVEMap, myEFMap, anOV);
          FvT    = Fv;
          Fv     = aFV;
          Etg    = anEtg;
          Eunder = anE;
          Vp     = aVp;
          LinTg  = aLin;
          OFvT   = anOT;
          OFv    = anOV;
        }
      }
      if (!FvT.IsNull())
      {
        FvT.Orientation(TopAbs_FORWARD);
      }
    }

    Fv.Orientation(TopAbs_FORWARD);
    Fad.Orientation(TopAbs_FORWARD);
    Fop.Orientation(TopAbs_FORWARD);

    // The edge that will be extended is returned: the edge of the vertex
    // between Fv and Fop -- or a face tangent to Fop across an edge of the
    // vertex (one wall in two faces). Arcpiv need not be an edge of the
    // vertex: when the fillet is wider than Fad's side there, its common
    // point lies on a face beyond it, and Fv has another edge at the vertex
    // (Fad's old side), which is not the one to extend.
    for (int pass = 0; pass < 3 && Arcprol.IsNull(); pass++)
    {
      for (It.Initialize(myVEMap(Vtx)); It.More() && Arcprol.IsNull(); It.Next())
      {
        if (Arcpiv.IsSame(It.Value()))
        {
          continue;
        }
        if (pass < 2)
        {
          bool isOp = false;
          for (NCollection_List<TopoDS_Shape>::Iterator itF(myEFMap(It.Value()));
               itF.More() && !isOp;
               itF.Next())
          {
            const TopoDS_Face& aF   = TopoDS::Face(itF.Value());
            TopAbs_Orientation aBid = TopAbs_FORWARD;
            isOp = !aF.IsSame(Fv)
                   && (pass == 0 ? aF.IsSame(Fop)
                                 : !TangentNeighbour(Vtx, Fop, aF, myVEMap, myEFMap, aBid).IsNull());
          }
          if (!isOp)
          {
            continue;
          }
        }
        for (ex.Init(Fv, TopAbs_EDGE); ex.More(); ex.Next())
        {
          if (It.Value().IsSame(ex.Current()))
          {
            Arcprol   = TopoDS::Edge(It.Value());
            OArcprolv = ex.Current().Orientation();
            break;
          }
        }
      }
    }
    if (Arcprol.IsNull()) /*throw StdFail_NotDone("OneCorner : edge a prolonger non trouve");*/
    {
      PerformIntersectionAtEnd(Index);
      return;
    }
    // The fillet is wider than Fad's side at Vtx (coplanar faces of one wall,
    // say): the edges of Fv on the way from Vtx to Arcpiv lie under the
    // fillet, and nothing else cuts them away.
    if (!containV(Fad, Vtx))
    {
      Swallowed = EdgesToArc(Fv, Vtx, Arcprol, Arcpiv);
    }
    // The edge of Fad from Arcpiv's vertex to Vtx lies under the fillet.
    if (!FvT.IsNull())
    {
      Swallowed.Append(Eunder);
      Swallowed.Append(Vtx);
    }
    for (ex.Init(Fop, TopAbs_EDGE); ex.More(); ex.Next())
    {
      if (Arcprol.IsSame(ex.Current()))
      {
        OArcprolop = ex.Current().Orientation();
        break;
      }
    }
    if (!ex.More())
    {
      // Arcprol is an edge of a face tangent to Fop: its orientation there,
      // in the shell, stands for the one it would have in Fop.
      TopAbs_Orientation aOFop = TopAbs_FORWARD;
      for (It.Initialize(myEFMap(Arcprol)); It.More() && Etan.IsNull(); It.Next())
      {
        if (!Fv.IsSame(It.Value()))
        {
          Fprol = TopoDS::Face(It.Value());
          Etan  = TangentNeighbour(Vtx, Fop, Fprol, myVEMap, myEFMap, aOFop);
        }
      }
      if (!Etan.IsNull())
      {
        for (ex.Init(Fprol, TopAbs_EDGE); ex.More(); ex.Next())
        {
          if (Arcprol.IsSame(ex.Current()))
          {
            OArcprolop = aOFop == TopAbs_REVERSED ? TopAbs::Reverse(ex.Current().Orientation())
                                                  : ex.Current().Orientation();
            break;
          }
        }
      }
    }
    TopoDS_Face               FFv;
    double                    tol;
    int                       prol = 0;
    BRep_Builder              BRE;
    occ::handle<Geom_Surface> Sface;
    Sface = BRep_Tool::Surface(Fv);
    // A trimmed surface (a drafted wall's cone, say) stops at the face no
    // less than the face does: the cut and the extension lie on its basis
    // past the trim, in the same parameters.
    const occ::handle<Geom_RectangularTrimmedSurface> aTrimmed =
      occ::down_cast<Geom_RectangularTrimmedSurface>(Sface);
    if (!aTrimmed.IsNull())
    {
      Sface = aTrimmed->BasisSurface();
    }
    ChFi3d_ExtendSurface(Sface, prol);
    tol = BRep_Tool::Tolerance(Fv);
    BRE.MakeFace(FFv, Sface, tol);
    if (prol || !aTrimmed.IsNull())
    {
      Bs.Initialize(FFv, false);
      DStr.SetNewSurface(Fv, Sface);
    }
    else
    {
      Bs.Initialize(Fv, false);
    }
    if (!FvT.IsNull())
    {
      occ::handle<Geom_Surface> SfaceT = BRep_Tool::Surface(FvT);
      const occ::handle<Geom_RectangularTrimmedSurface> aTrimmedT =
        occ::down_cast<Geom_RectangularTrimmedSurface>(SfaceT);
      if (!aTrimmedT.IsNull())
      {
        SfaceT = aTrimmedT->BasisSurface();
      }
      int prolT = 0;
      ChFi3d_ExtendSurface(SfaceT, prolT);
      TopoDS_Face FFvT;
      BRE.MakeFace(FFvT, SfaceT, BRep_Tool::Tolerance(FvT));
      HBsT->Initialize(FFvT, false);
      if (prolT || !aTrimmedT.IsNull())
      {
        DStr.SetNewSurface(FvT, SfaceT);
      }
    }
    Bad.Initialize(Fad);
    Bop.Initialize(Fop);
  }
  // in case of OnSame it is necessary to modify the CommonPoint
  // in the empty and its parameter in the FaceInterference.
  // They are both returned in non const references. Attention the modifications are done behind
  // de CV1,CV2,Fi1,Fi2.
  ChFiDS_CommonPoint&      CPopArc = Fd->ChangeVertex(isfirst, IFopArc);
  ChFiDS_FaceInterference& FiopArc = Fd->ChangeInterference(IFopArc);
  ChFiDS_CommonPoint&      CPadArc = Fd->ChangeVertex(isfirst, IFadArc);
  ChFiDS_FaceInterference& FiadArc = Fd->ChangeInterference(IFadArc);
  // the parameter of the vertex in the air is initialiced with the value of
  // its opposite (point on arc).
  double                           wop = Fd->ChangeInterference(IFadArc).Parameter(isfirst);
  occ::handle<Geom_Curve>          c3df;
  occ::handle<GeomAdaptor_Surface> HGs =
    new GeomAdaptor_Surface(DStr.Surface(Fd->Surf()).Surface());
  gp_Pnt2d p2dbout;

  if (onsame)
  {

    ChFiDS_CommonPoint saveCPopArc = CPopArc;
    c3df                           = DStr.Curve(FiopArc.LineIndex()).Curve();

    inters = IntersUpdateOnSame(HGs,
                                HBs,
                                c3df,
                                Fop,
                                Fv,
                                Arcprol,
                                Vtx,
                                isfirst,
                                10 * tolapp3d, // in
                                FiopArc,
                                CPopArc,
                                p2dbout,
                                wop); // out

    occ::handle<BRepAdaptor_Curve2d> pced = new BRepAdaptor_Curve2d();
    pced->Initialize(CPadArc.Arc(), FvT.IsNull() ? Fv : FvT);
    // in the case of degenerated Fi, parameter difference can be even negative (eap, occ293)
    if ((FiadArc.LastParameter() - FiadArc.FirstParameter()) > 10 * tolesp)
    {
      Update(FvT.IsNull() ? HBs : HBsT, pced, HGs, FiadArc, CPadArc, isfirst);
    }

    // The fillet's line on Fop ended on Etan, the edge between Fop and the
    // face Arcprol is extended in: its end is still there, and the
    // extension runs along Etan. Or the line ended inside Fop, and its
    // update to Fv's surface put the end on Etan: Fv's surface holds Etan
    // from there to Vtx (Fv and Fprol one smooth wall, Fop tangent to it
    // along Etan), so the extension runs along Etan just the same.
    if (inters && !Etan.IsNull() && !CPopArc.IsOnArc())
    {
      const bool        wasOnEtan = saveCPopArc.IsOnArc() && saveCPopArc.Arc().IsSame(Etan);
      BRepAdaptor_Curve aCEtan(Etan);
      Extrema_ExtPC     anExt(CPopArc.Point(), aCEtan);
      const double      aTol = std::max(wasOnEtan ? saveCPopArc.Tolerance() : CPopArc.Tolerance(),
                                   10 * tolapp3d);
      int               iMin = 0;
      double            dMin = aTol * aTol;
      for (int i = 1; anExt.IsDone() && i <= anExt.NbExt(); i++)
      {
        if (anExt.SquareDistance(i) <= dMin)
        {
          dMin = anExt.SquareDistance(i);
          iMin = i;
        }
      }
      if (iMin > 0 && wasOnEtan)
      {
        CPopArc.SetArc(aTol, Etan, anExt.Point(iMin).Parameter(), saveCPopArc.TransitionOnArc());
        zobOnEtan = true;
      }
      else if (iMin > 0)
      {
        const double aPar  = anExt.Point(iMin).Parameter();
        const double aParV = BRep_Tool::Parameter(Vtx, Etan);
        // Etan between the end and Vtx lies on Fv's surface
        bool onFv = std::abs(aParV - aPar) > Precision::PConfusion();
        for (int k = 1; k <= 3 && onFv; k++)
        {
          const gp_Pnt               aP = aCEtan.Value(aPar + (aParV - aPar) * k / 4.);
          GeomAPI_ProjectPointOnSurf aProj(aP, BRep_Tool::Surface(Fv));
          onFv = aProj.NbPoints() > 0 && aProj.LowerDistance() <= aTol;
        }
        if (onFv)
        {
          // the transition the walk gives a line ending on an arc
          TopAbs_Orientation anOr = Etan.Orientation();
          for (ex.Init(Fop, TopAbs_EDGE); ex.More(); ex.Next())
          {
            if (Etan.IsSame(ex.Current()))
            {
              anOr = ex.Current().Orientation();
              break;
            }
          }
          TopAbs_Orientation aTr = TopAbs::Compose(FiopArc.Transition(), anOr);
          if (isfirst)
          {
            aTr = TopAbs::Reverse(aTr);
          }
          CPopArc.SetArc(aTol, Etan, aPar, aTr);
          zobOnEtan = true;
        }
      }
    }

    if (thePrepareOnSame)
    {
      // saveCPopArc.SetParameter(wop);
      saveCPopArc.SetPoint(CPopArc.Point());
      CPopArc = saveCPopArc;
      return;
    }
  }
  else
  {
    inters = FindFace(Vtx, CV1, CV2, Fv, Fop);
    // The walk stopped at the spine's end with one side on an arc of Vtx
    // and the other in its face, a plane: the line there, carried on,
    // crosses into a coplanar piece of the same wall and meets the face at
    // end beyond (a wall kept in pieces, the radius wider than the piece's
    // edge at Vtx). The line ends on that edge of the piece, as it would on
    // the wall in one face, and is carried over the split below.
    if (!inters && CV1.IsOnArc() != CV2.IsOnArc())
    {
      const int                 onsArc = CV1.IsOnArc() ? 1 : 2;
      const int                 onsAir = 3 - onsArc;
      const ChFiDS_CommonPoint& aCPArc = Fd->Vertex(isfirst, onsArc);
      const TopoDS_Face         aFArc  = TopoDS::Face(DStr.Shape(Fd->Index(onsArc)));
      const TopoDS_Face         aFAir  = TopoDS::Face(DStr.Shape(Fd->Index(onsAir)));
      ChFiDS_FaceInterference&  aFi    = Fd->ChangeInterference(onsAir);
      TopoDS_Face               aFv;
      for (NCollection_List<TopoDS_Shape>::Iterator itF(myEFMap(aCPArc.Arc())); itF.More();
           itF.Next())
      {
        if (!aFArc.IsSame(itF.Value()))
        {
          aFv = TopoDS::Face(itF.Value());
        }
      }
      if (hasVertex(aCPArc.Arc(), Vtx) && !aFv.IsNull() && containV(aFv, Vtx)
          && aFi.LineIndex() != 0
          && LineOverSplit(DStr.Curve(aFi.LineIndex()).Curve(),
                           aFi.Parameter(isfirst),
                           isfirst,
                           aFAir,
                           aFv,
                           myEFMap,
                           10 * tolapp3d,
                           LinOS,
                           EsplitOS,
                           FnOS,
                           EarcOS,
                           wQOS,
                           wPOS,
                           parQOS,
                           parPOS))
      {
        ChFiDS_CommonPoint&      aCPAir  = Fd->ChangeVertex(isfirst, onsAir);
        const ChFiDS_CommonPoint aSaved  = aCPAir;
        const double             aSavedW = aFi.Parameter(isfirst);
        aCPAir.Reset();
        aCPAir.SetPoint(LinOS->Value(wPOS));
        aCPAir.SetArc(std::max(aSaved.Tolerance(), 10 * tolapp3d),
                      EarcOS,
                      parPOS,
                      TransitionOnEdge(aFi,
                                       FnOS,
                                       EarcOS,
                                       isfirst,
                                       FnOS.Orientation() != aFAir.Orientation()));
        aFi.SetParameter(wPOS, isfirst);
        inters = FindFace(Vtx, CV1, CV2, Fv, Fop);
        if (inters)
        {
          onsOS = onsAir;
        }
        else
        {
          aCPAir = aSaved;
          aFi.SetParameter(aSavedW, isfirst);
        }
      }
    }
    if (!inters)
    {
      PerformIntersectionAtEnd(Index);
      return;
    }
    // A common point on an arc that does not reach Vtx: the fillet is
    // wider than the face beside Vtx, and the edges of Fv on the way
    // there lie under it, as in the case OnSame.
    for (int ons = 1; ons <= 2; ons++)
    {
      const ChFiDS_CommonPoint& aCP    = Fd->Vertex(isfirst, ons);
      const ChFiDS_CommonPoint& aCPOpp = Fd->Vertex(isfirst, 3 - ons);
      if (aCP.IsOnArc() && aCPOpp.IsOnArc() && !hasVertex(aCP.Arc(), Vtx)
          && hasVertex(aCPOpp.Arc(), Vtx))
      {
        NCollection_List<TopoDS_Shape> aWay = EdgesToArc(Fv, Vtx, aCPOpp.Arc(), aCP.Arc());
        Swallowed.Append(aWay);
      }
    }
    Bs.Initialize(Fv);
    occ::handle<BRepAdaptor_Curve2d> pced = new BRepAdaptor_Curve2d();
    pced->Initialize(CV1.Arc(), Fv);
    Update(HBs, pced, HGs, Fd->ChangeInterferenceOnS1(), CV1, isfirst);
    pced->Initialize(CV2.Arc(), Fv);
    Update(HBs, pced, HGs, Fd->ChangeInterferenceOnS2(), CV2, isfirst);
  }

#ifdef OCCT_DEBUG
  ChFi3d_ResultChron(ch, t_same); // result perf condition if (same)
  ChFi3d_InitChron(ch);           // init perf condition if (inters)
#endif

  TopoDS_Edge                    edgecouture;
  bool                           couture, intcouture = false;
  double                         tolreached = tolapp3d;
  double                         par1 = 0., par2 = 0.;
  int                            indpt = 0, Icurv1 = 0, Icurv2 = 0;
  occ::handle<Geom_TrimmedCurve> curv1, curv2;
  occ::handle<Geom2d_Curve>      c2d1, c2d2;

  int Isurf = Fd->Surf();
  // the cut's end on Fop's side, in Fv's parameters
  gp_Pnt2d aCutOnOp;
  // The cut over two faces (FvT set): its piece on FvT (Cc is the piece on
  // Fv), the point between them, where it crosses Etg's line at wTg, and
  // that line from Vp to there with its pcurves on FvT and Fv
  occ::handle<Geom_Curve>   CcT, CTg;
  occ::handle<Geom2d_Curve> PsT, PcT, PTgT, PTgV;
  gp_Pnt                    PTg;
  double                    wTg = 0.;
  // the tolerances the two pieces reached, each its own
  double tolT = tolapp3d, tolV = tolapp3d;

  if (inters)
  {
    HGs                                = ChFi3d_BoundSurf(DStr, Fd, 1, 2);
    const ChFiDS_FaceInterference& Fi1 = Fd->InterferenceOnS1();
    const ChFiDS_FaceInterference& Fi2 = Fd->InterferenceOnS2();
    NCollection_Array1<double>     Pardeb(1, 4), Parfin(1, 4);
    gp_Pnt2d                       pfil1, pfac1, pfil2, pfac2;
    occ::handle<Geom2d_Curve>      Hc1, Hc2;
    if (onsame && IFopArc == 1)
    {
      pfac1 = p2dbout;
    }
    else
    {
      Hc1 = PCurveInFace(CV1.Arc(), FvT.IsNull() ? Fv : FvT, Ubid, Ubid);
      if (Hc1.IsNull())
      {
        throw Standard_ConstructionError("Failed to get p-curve of edge");
      }
      pfac1 = Hc1->Value(CV1.ParameterOnArc());
    }
    if (onsame && IFopArc == 2)
    {
      pfac2 = p2dbout;
    }
    else
    {
      Hc2 = PCurveInFace(CV2.Arc(), FvT.IsNull() ? Fv : FvT, Ubid, Ubid);
      if (Hc2.IsNull())
      {
        throw Standard_ConstructionError("Failed to get p-curve of edge");
      }
      pfac2 = Hc2->Value(CV2.ParameterOnArc());
    }
    if (Fi1.LineIndex() != 0)
    {
      pfil1 = Fi1.PCurveOnSurf()->Value(Fi1.Parameter(isfirst));
    }
    else
    {
      pfil1 = Fi1.PCurveOnSurf()->Value(Fi1.Parameter(!isfirst));
    }
    if (Fi2.LineIndex() != 0)
    {
      pfil2 = Fi2.PCurveOnSurf()->Value(Fi2.Parameter(isfirst));
    }
    else
    {
      pfil2 = Fi2.PCurveOnSurf()->Value(Fi2.Parameter(!isfirst));
    }
    if (!FvT.IsNull())
    {
      // The point where the cut crosses from FvT to Fv: on Etg's line
      // carried on past Vp, the nearest of its intersections with the
      // fillet's surface within the surface's bounds.
      const occ::handle<Geom_Surface>& aSfil = DStr.Surface(Fd->Surf()).Surface();
      GeomAPI_IntCS                    anInt(new Geom_Line(LinTg), aSfil);
      const double aU1 = HGs->FirstUParameter(), aU2 = HGs->LastUParameter();
      const double aV1 = HGs->FirstVParameter(), aV2 = HGs->LastVParameter();
      const double aMargU = 0.1 * (aU2 - aU1), aMargV = 0.1 * (aV2 - aV1);
      gp_Pnt2d     pfilT;
      bool         isFound = false;
      for (int i = 1; anInt.IsDone() && i <= anInt.NbPoints(); i++)
      {
        double u, v, w;
        anInt.Parameters(i, u, v, w);
        if (aSfil->IsUPeriodic())
        {
          u = ElCLib::InPeriod(u, pfil1.X() - 0.5 * aSfil->UPeriod(),
                               pfil1.X() + 0.5 * aSfil->UPeriod());
        }
        if (aSfil->IsVPeriodic())
        {
          v = ElCLib::InPeriod(v, pfil1.Y() - 0.5 * aSfil->VPeriod(),
                               pfil1.Y() + 0.5 * aSfil->VPeriod());
        }
        if (w <= Precision::Confusion() || u < aU1 - aMargU || u > aU2 + aMargU || v < aV1 - aMargV
            || v > aV2 + aMargV || (isFound && w >= wTg))
        {
          continue;
        }
        isFound = true;
        wTg     = w;
        pfilT.SetCoord(u, v);
        PTg = anInt.Point(i);
      }
      if (!isFound)
      {
        throw Standard_Failure("OneCorner : the cut does not cross the wall's tangent edge");
      }
      CTg  = new Geom_TrimmedCurve(new Geom_Line(LinTg), 0., wTg);
      PTgT = CurveOnFacePeriod(CTg, HBsT->Surface().Surface(), FvT, Vp, 10 * tolapp3d);
      PTgV = CurveOnFacePeriod(CTg, HBs->Surface().Surface(), Fv, Vp, 10 * tolapp3d);
      if (PTgT.IsNull() || PTgV.IsNull())
      {
        throw Standard_Failure("OneCorner : the wall's tangent edge does not carry on");
      }
      gp_Pnt2d       pfacT  = PTgT->Value(wTg);
      gp_Pnt2d       pfacV  = PTgV->Value(wTg);
      gp_Pnt2d&      pfacAd = IFadArc == 1 ? pfac1 : pfac2;
      gp_Pnt2d&      pfacOp = IFopArc == 1 ? pfac1 : pfac2;
      ChFi3d_Recale(*HBsT, pfacT, pfacAd, true);
      ChFi3d_Recale(Bs, pfacV, pfacOp, true);
      aCutOnOp = pfacOp;
      // The cut's end on Fop's side lies on Fv beside Eunder, between Vtx
      // and Etg: Fv's surface carried on (a B-spline's, say) can meet the
      // fillet's line again far off, nearer to where the walk ended.
      {
        double                          f, l;
        const occ::handle<Geom2d_Curve> aPcU = BRep_Tool::CurveOnSurface(Eunder, Fv, f, l);
        bool                            isBeside = false;
        if (!aPcU.IsNull())
        {
          const double                  aMarg = 0.05 * (l - f);
          Geom2dAPI_ProjectPointOnCurve aPr(pfacOp, aPcU, f - aMarg, l + aMarg);
          isBeside = aPr.NbPoints() > 0 && aPr.LowerDistanceParameter() >= f - aMarg
                     && aPr.LowerDistanceParameter() <= l + aMarg;
        }
        if (!isBeside)
        {
          throw Standard_Failure("OneCorner : the cut's end on the side is off the face at Vtx");
        }
      }

      // The piece on FvT runs from Arcpiv to the point, the one on Fv from
      // there to its end on Fop's side. A curve starting at the point, at
      // the end of the fillet's line, can come out loose (on a seam of Fv
      // there, say); the one from the end on Fop's side to the fillet's
      // line on Fv's surface carried on past Etg, cut at the point, can
      // cross Fv's seam. Of the two, the one that fits better.
      gp_Pnt2d pfacVAd;
      {
        const gp_Pnt2d&            pfilAd = IFadArc == 1 ? pfil1 : pfil2;
        GeomAPI_ProjectPointOnSurf aProj(HGs->Value(pfilAd.X(), pfilAd.Y()),
                                         HBs->Surface().Surface());
        if (aProj.NbPoints() > 0)
        {
          double u, v;
          aProj.LowerDistanceParameters(u, v);
          pfacVAd.SetCoord(u, v);
          ChFi3d_Recale(Bs, pfacV, pfacVAd, true);
        }
        else
        {
          pfacVAd = pfacV;
        }
      }
      // the piece the k-th along the cut (from CV1), on FvT or Fv, starting
      // or ending at the point, or on Fv through it to the fillet's line
      auto aPiece = [&](const int                  k,
                        const bool                 theWhole,
                        occ::handle<Geom_Curve>&   theC,
                        occ::handle<Geom2d_Curve>& thePs,
                        occ::handle<Geom2d_Curve>& thePc,
                        double&                    theTol) -> bool {
        const bool           onT    = (k == IFadArc);
        BRepAdaptor_Surface& aBs    = onT ? *HBsT : Bs;
        const gp_Pnt2d&      pfacTV = onT ? pfacT : (theWhole ? pfacVAd : pfacV);
        const gp_Pnt2d&      pfilA  = k == 2 && !theWhole ? pfilT : pfil1;
        const gp_Pnt2d&      pfilB  = k == 1 && !theWhole ? pfilT : pfil2;
        const gp_Pnt2d&      pfacA  = k == 1 ? pfac1 : pfacTV;
        const gp_Pnt2d&      pfacB  = k == 1 ? pfacTV : pfac2;
        Pardeb(1)                   = pfilA.X();
        Pardeb(2)                   = pfilA.Y();
        Pardeb(3)                   = pfacA.X();
        Pardeb(4)                   = pfacA.Y();
        Parfin(1)                   = pfilB.X();
        Parfin(2)                   = pfilB.Y();
        Parfin(3)                   = pfacB.X();
        Parfin(4)                   = pfacB.Y();
        double uu1, uu2, vv1, vv2;
        ChFi3d_Boite(pfacA, pfacB, uu1, uu2, vv1, vv2);
        ChFi3d_BoundFac(aBs, uu1, uu2, vv1, vv2);
        theTol = tolapp3d;
        if (!ChFi3d_ComputeCurves(HGs,
                                  onT ? HBsT : HBs,
                                  Pardeb,
                                  Parfin,
                                  theC,
                                  thePs,
                                  thePc,
                                  tolapp3d,
                                  tol2d,
                                  theTol))
        {
          return false;
        }
        if (theWhole)
        {
          GeomAPI_ProjectPointOnCurve aProjC(PTg, theC);
          if (aProjC.NbPoints() == 0 || aProjC.LowerDistance() > 10 * tolapp3d)
          {
            return false;
          }
          const double aPar = aProjC.LowerDistanceParameter();
          theC = k == 2 ? new Geom_TrimmedCurve(theC, aPar, theC->LastParameter())
                        : new Geom_TrimmedCurve(theC, theC->FirstParameter(), aPar);
          // the piece kept's own: past the point the curve runs on off the
          // parts of the surfaces it is computed on
          theTol = std::max(ChFi3d_EvalTolReached(HGs, thePs, HBs, thePc, theC),
                            aProjC.LowerDistance());
        }
        return true;
      };
      double aTolT = 0., aTolV = 0.;
      if (!aPiece(IFadArc, false, CcT, PsT, PcT, aTolT))
      {
        throw Standard_Failure("OneCorner : echec calcul intersection");
      }
      const bool isCut = aPiece(IFopArc, false, Cc, Ps, Pc, aTolV);
      if (!isCut || aTolV > tolapp3d)
      {
        occ::handle<Geom_Curve>   aC;
        occ::handle<Geom2d_Curve> aPs, aPc;
        double                    aTolW = 0.;
        if (aPiece(IFopArc, true, aC, aPs, aPc, aTolW) && (!isCut || aTolW < aTolV))
        {
          Cc    = aC;
          Ps    = aPs;
          Pc    = aPc;
          aTolV = aTolW;
        }
        else if (!isCut)
        {
          throw Standard_Failure("OneCorner : echec calcul intersection");
        }
      }
      tolreached = std::max(aTolT, aTolV);
      tolT       = aTolT;
      tolV       = aTolV;
      Udeb = Cc->FirstParameter();
      Ufin = Cc->LastParameter();
      couture = false;
    }
    else
    {
      if (onsame)
      {
        ChFi3d_Recale(Bs, pfac1, pfac2, (IFadArc == 1));
      }
      aCutOnOp = IFopArc == 1 ? pfac1 : pfac2;

      Pardeb(1) = pfil1.X();
      Pardeb(2) = pfil1.Y();
      Pardeb(3) = pfac1.X();
      Pardeb(4) = pfac1.Y();
      Parfin(1) = pfil2.X();
      Parfin(2) = pfil2.Y();
      Parfin(3) = pfac2.X();
      Parfin(4) = pfac2.Y();

      double uu1, uu2, vv1, vv2;
      ChFi3d_Boite(pfac1, pfac2, uu1, uu2, vv1, vv2);
      ChFi3d_BoundFac(Bs, uu1, uu2, vv1, vv2);

      if (!ChFi3d_ComputeCurves(HGs, HBs, Pardeb, Parfin, Cc, Ps, Pc, tolapp3d, tol2d, tolreached))
      {
        throw Standard_Failure("OneCorner : echec calcul intersection");
      }

      Udeb = Cc->FirstParameter();
      Ufin = Cc->LastParameter();

      //  determine if the curve has an intersection with edge of sewing

      ChFi3d_Couture(Fv, couture, edgecouture);
    }

    // a curve ending tangent to its face's boundary kept on the face
    const double aKeptV = KeepCurveOnFace(Pc, Fv, Udeb, Ufin);
    tolV += aKeptV;
    tolreached += aKeptV;
    if (!FvT.IsNull())
    {
      const double aKeptT = KeepCurveOnFace(PcT, FvT, CcT->FirstParameter(), CcT->LastParameter());
      tolT += aKeptT;
      tolreached = std::max(tolreached, tolT);
    }

    if (couture && !BRep_Tool::Degenerated(edgecouture))
    {

      // double Ubid,Vbid;
      occ::handle<Geom_Curve>        C     = BRep_Tool::Curve(edgecouture, Ubid, Vbid);
      occ::handle<Geom_TrimmedCurve> Ctrim = new Geom_TrimmedCurve(C, Ubid, Vbid);
      GeomAdaptor_Curve              cur1(Ctrim->BasisCurve());
      GeomAdaptor_Curve              cur2(Cc);
      Extrema_ExtCC                  extCC(cur1, cur2);
      if (extCC.IsDone() && extCC.NbExt() != 0)
      {
        int    imin     = 0;
        double distmin2 = RealLast();
        for (int i = 1; i <= extCC.NbExt(); i++)
        {
          if (extCC.SquareDistance(i) < distmin2)
          {
            distmin2 = extCC.SquareDistance(i);
            imin     = i;
          }
        }
        if (distmin2 <= Precision::SquareConfusion())
        {
          Extrema_POnCurv ponc1, ponc2;
          extCC.Points(imin, ponc1, ponc2);
          par1       = ponc1.Parameter();
          par2       = ponc2.Parameter();
          double Tol = 1.e-4;
          if (std::abs(par2 - Udeb) > Tol && std::abs(Ufin - par2) > Tol)
          {
            gp_Pnt             P1 = ponc1.Value();
            TopOpeBRepDS_Point tpoint(P1, Tol);
            indpt      = DStr.AddPoint(tpoint);
            intcouture = true;
            curv1      = new Geom_TrimmedCurve(Cc, Udeb, par2);
            curv2      = new Geom_TrimmedCurve(Cc, par2, Ufin);
            TopOpeBRepDS_Curve tcurv1(curv1, tolreached);
            TopOpeBRepDS_Curve tcurv2(curv2, tolreached);
            Icurv1 = DStr.AddCurve(tcurv1);
            Icurv2 = DStr.AddCurve(tcurv2);
          }
        }
      }
    }
  }
  else
  { // (!inters)
    throw Standard_NotImplemented("OneCorner : bouchon non ecrit");
  }
  int                IShape = DStr.AddShape(Fv);
  TopAbs_Orientation Et     = TopAbs_FORWARD;
  // the face Arcpiv bounds
  const TopoDS_Face& FvArc = FvT.IsNull() ? Fv : FvT;
  if (IFadArc == 1)
  {
    TopExp_Explorer Exp;
    for (Exp.Init(FvArc.Oriented(TopAbs_FORWARD), TopAbs_EDGE); Exp.More(); Exp.Next())
    {
      if (Exp.Current().IsSame(CV1.Arc()))
      {
        Et = TopAbs::Reverse(TopAbs::Compose(Exp.Current().Orientation(), CV1.TransitionOnArc()));
        break;
      }
    }
  }
  else
  {
    TopExp_Explorer Exp;
    for (Exp.Init(FvArc.Oriented(TopAbs_FORWARD), TopAbs_EDGE); Exp.More(); Exp.Next())
    {
      if (Exp.Current().IsSame(CV2.Arc()))
      {
        Et = TopAbs::Compose(Exp.Current().Orientation(), CV2.TransitionOnArc());
        break;
      }
    }

    //
  }

#ifdef OCCT_DEBUG
  ChFi3d_ResultChron(ch, t_inter); // result perf condition if (inter)
  ChFi3d_InitChron(ch);            // init perf condition  if (onsame && inters)
#endif

  stripe->SetIndexPoint(ChFi3d_IndexPointInDS(CV1, DStr), isfirst, 1);
  stripe->SetIndexPoint(ChFi3d_IndexPointInDS(CV2, DStr), isfirst, 2);

  // the boxes of the cut's points: the one between its pieces
  Bnd_Box boxT;
  if (!FvT.IsNull())
  {
    // The cut over two faces, in two curves meeting where it crosses Etg's
    // line carried on past Vp; the line from Vp to there bounds both faces.
    // The curves and their points go in the DS here, as for a cut crossing
    // a seam (below), and FILDS leaves the end alone. Each curve keeps the
    // tolerance it reached -- the piece over Fv's seam can come out looser
    // than the others -- and the point between them the largest.
    indpt           = DStr.AddPoint(TopOpeBRepDS_Point(PTg, tolreached));
    const int ICT   = DStr.AddCurve(TopOpeBRepDS_Curve(CcT, tolT));
    const int ICV   = DStr.AddCurve(TopOpeBRepDS_Curve(Cc, tolV));
    const int IShT  = DStr.AddShape(FvT);
    // Et is FvT's, from Arcpiv; Fv's is the same in the shell
    const TopAbs_Orientation EtV = OFvT == OFv ? Et : TopAbs::Reverse(Et);
    DStr.ChangeShapeInterferences(IShT).Append(ChFi3d_FilCurveInDS(ICT, IShT, PcT, Et));
    DStr.ChangeShapeInterferences(IShape).Append(ChFi3d_FilCurveInDS(ICV, IShape, Pc, EtV));
    // on the fillet's surface, as FILDS puts the end's curve
    TopAbs_Orientation aTrafil1 = TopAbs_FORWARD;
    if (Fd->IndexOfS1() > 0)
    {
      aTrafil1 = DStr.Shape(Fd->IndexOfS1()).Orientation();
    }
    aTrafil1 = TopAbs::Compose(aTrafil1, Fd->Orientation());
    aTrafil1 = TopAbs::Compose(TopAbs::Reverse(Fd->InterferenceOnS1().Transition()), aTrafil1);
    const TopAbs_Orientation EtS = isfirst ? TopAbs::Reverse(aTrafil1) : aTrafil1;
    DStr.ChangeSurfaceInterferences(Isurf).Append(ChFi3d_FilCurveInDS(ICT, Isurf, PsT, EtS));
    DStr.ChangeSurfaceInterferences(Isurf).Append(ChFi3d_FilCurveInDS(ICV, Isurf, Ps, EtS));
    stripe->InDS(isfirst);
    const int                      ind1 = stripe->IndexPoint(isfirst, 1);
    const int                      ind2 = stripe->IndexPoint(isfirst, 2);
    const int                      IC1  = IFadArc == 1 ? ICT : ICV;
    const int                      IC2  = IFadArc == 1 ? ICV : ICT;
    const occ::handle<Geom_Curve>& C1   = IFadArc == 1 ? CcT : Cc;
    const occ::handle<Geom_Curve>& C2   = IFadArc == 1 ? Cc : CcT;
    DStr.ChangeCurveInterferences(IC1).Append(
      ChFi3d_FilPointInDS(TopAbs_FORWARD, IC1, ind1, C1->FirstParameter()));
    DStr.ChangeCurveInterferences(IC1).Append(
      ChFi3d_FilPointInDS(TopAbs_REVERSED, IC1, indpt, C1->LastParameter()));
    DStr.ChangeCurveInterferences(IC2).Append(
      ChFi3d_FilPointInDS(TopAbs_FORWARD, IC2, indpt, C2->FirstParameter()));
    DStr.ChangeCurveInterferences(IC2).Append(
      ChFi3d_FilPointInDS(TopAbs_REVERSED, IC2, ind2, C2->LastParameter()));

    // Etg's line from Vp, on each face as Etg is there: the same way when
    // Etg runs to Vp
    TopoDS_Vertex aV1, aV2;
    TopExp::Vertices(TopoDS::Edge(Etg.Oriented(TopAbs_FORWARD)), aV1, aV2);
    const bool         sameDir = Vp.IsSame(aV2);
    const double tolG =
      std::max(ChFi3d_EvalTolReached(HBsT, PTgT, HBs, PTgV, CTg), Precision::Confusion());
    const int    IG   = DStr.AddCurve(TopOpeBRepDS_Curve(CTg, tolG));
    for (int k = 0; k < 2; k++)
    {
      const TopoDS_Face& aF  = k == 0 ? FvT : Fv;
      TopAbs_Orientation anO = TopAbs_FORWARD;
      for (ex.Init(aF.Oriented(TopAbs_FORWARD), TopAbs_EDGE); ex.More(); ex.Next())
      {
        if (Etg.IsSame(ex.Current()))
        {
          anO = ex.Current().Orientation();
          break;
        }
      }
      DStr.ChangeShapeInterferences(k == 0 ? IShT : IShape)
        .Append(ChFi3d_FilCurveInDS(IG,
                                    k == 0 ? IShT : IShape,
                                    k == 0 ? PTgT : PTgV,
                                    sameDir ? anO : TopAbs::Reverse(anO)));
    }
    DStr.ChangeCurveInterferences(IG).Append(
      ChFi3d_FilVertexInDS(TopAbs_FORWARD, IG, DStr.AddShape(Vp), 0.));
    DStr.ChangeCurveInterferences(IG).Append(
      ChFi3d_FilPointInDS(TopAbs_REVERSED, IG, indpt, wTg));
    ChFi3d_EnlargeBox(IFadArc == 1 ? HBsT : HBs, IFadArc == 1 ? PcT : Pc,
                      C1->FirstParameter(), C1->LastParameter(), box1, boxT);
    ChFi3d_EnlargeBox(IFadArc == 1 ? HBs : HBsT, IFadArc == 1 ? Pc : PcT,
                      C2->FirstParameter(), C2->LastParameter(), boxT, box2);
  }
  else if (onsOS != 0)
  {
    // The line over a wall in coplanar pieces: its piece on the far one,
    // from where it crosses the split (Q) to the cut's end (P), and the
    // cut, go in the DS here, as for the cut over two faces above; the
    // line on the face the spine is on ends at Q, and FILDS ends it there.
    if (intcouture)
    {
      throw Standard_Failure("OneCorner : the cut over a split wall crosses a seam");
    }
    const int onsArc = 3 - onsOS;
    const int indP   = stripe->IndexPoint(isfirst, onsOS);
    const int indArc = stripe->IndexPoint(isfirst, onsArc);
    // on the fillet's surface, as FILDS puts the end's curve and the lines
    TopAbs_Orientation aTrafil1 = TopAbs_FORWARD;
    if (Fd->IndexOfS1() > 0)
    {
      aTrafil1 = DStr.Shape(Fd->IndexOfS1()).Orientation();
    }
    aTrafil1 = TopAbs::Compose(aTrafil1, Fd->Orientation());
    aTrafil1 = TopAbs::Compose(TopAbs::Reverse(Fd->InterferenceOnS1().Transition()), aTrafil1);
    const TopAbs_Orientation EtS = isfirst ? TopAbs::Reverse(aTrafil1) : aTrafil1;

    // the cut, from CV1's side to CV2's
    const int ICut = DStr.AddCurve(TopOpeBRepDS_Curve(Cc, tolreached));
    DStr.ChangeShapeInterferences(IShape).Append(ChFi3d_FilCurveInDS(ICut, IShape, Pc, Et));
    DStr.ChangeSurfaceInterferences(Isurf).Append(ChFi3d_FilCurveInDS(ICut, Isurf, Ps, EtS));
    DStr.ChangeCurveInterferences(ICut).Append(
      ChFi3d_FilPointInDS(TopAbs_FORWARD, ICut, onsOS == 1 ? indP : indArc, Udeb));
    DStr.ChangeCurveInterferences(ICut).Append(
      ChFi3d_FilPointInDS(TopAbs_REVERSED, ICut, onsOS == 2 ? indP : indArc, Ufin));

    // P on the edge of the far piece and the face at end
    ChFiDS_FaceInterference& aFi  = Fd->ChangeInterference(onsOS);
    const TopoDS_Face        aFOp = TopoDS::Face(DStr.Shape(Fd->Index(onsOS)));
    // the far piece the other way round on its surface: the line too
    const bool isRevN = FnOS.Orientation() != aFOp.Orientation();
    DStr.ChangeShapeInterferences(DStr.AddShape(EarcOS))
      .Append(ChFi3d_FilPointInDS(TransitionOnEdge(aFi, FnOS, EarcOS, isfirst, isRevN),
                                  DStr.AddShape(EarcOS),
                                  indP,
                                  parPOS));
    Bnd_Box boxP;
    ChFi3d_EnlargeBox(HBs, Pc, Udeb, Ufin, onsOS == 1 ? boxP : box1, onsOS == 2 ? boxP : box2);
    ChFi3d_EnlargeBox(EarcOS, myEFMap(EarcOS), parPOS, boxP);
    boxP.Add(LinOS->Value(wPOS));

    // the line ends at Q on the split, and its piece from Q to P
    StoreLineOverSplit(DStr,
                       stripe,
                       Fd,
                       isfirst,
                       onsOS,
                       Isurf,
                       HGs,
                       LinOS,
                       EsplitOS,
                       FnOS,
                       wQOS,
                       wPOS,
                       parQOS,
                       indP,
                       10 * tolapp3d,
                       myRegul);
    stripe->InDS(isfirst);
    (onsOS == 1 ? box1 : box2).Add(LinOS->Value(wQOS));
    ChFi3d_SetPointTolerance(DStr, boxP, indP);
  }
  else if (!intcouture)
  {
    // there is no intersection with the sewing edge
    // the curve Cc is stored in the stripe
    // the storage in the DS is not done by FILDS.

    TopOpeBRepDS_Curve Tc(Cc, tolreached);
    ICurve = DStr.AddCurve(Tc);
    occ::handle<TopOpeBRepDS_SurfaceCurveInterference> Interfc =
      ChFi3d_FilCurveInDS(ICurve, IShape, Pc, Et);

    // 31/01/02 akm vvv : (OCC119) Prevent the builder from creating
    //                    intersecting fillets - they are bad.
    Geom2dInt_GInter    anIntersector;
    Geom2dAdaptor_Curve aCorkPCurve(Pc, Udeb, Ufin);

    // Take all the interferences with faces from all the stripes
    // and look if their pcurves intersect our cork pcurve.
    // Unfortunately, by this moment they do not exist in DStr.
    NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator aStrIt(myListStripe);
    for (; aStrIt.More(); aStrIt.Next())
    {
      occ::handle<ChFiDS_Stripe> aCheckStripe = aStrIt.Value();
      occ::handle<NCollection_HSequence<occ::handle<ChFiDS_SurfData>>> aSeqData =
        aCheckStripe->SetOfSurfData();
      // Loop on parts of the stripe
      int iPart;
      for (iPart = 1; iPart <= aSeqData->Length(); iPart++)
      {
        occ::handle<ChFiDS_SurfData> aData = aSeqData->Value(iPart);
        Geom2dAdaptor_Curve          anOtherPCurve;
        if (IShape == aData->IndexOfS1())
        {
          const occ::handle<Geom2d_Curve>& aPCurve = aData->InterferenceOnS1().PCurveOnFace();
          if (aPCurve.IsNull())
          {
            continue;
          }

          anOtherPCurve.Load(aPCurve,
                             aData->InterferenceOnS1().FirstParameter(),
                             aData->InterferenceOnS1().LastParameter());
        }
        else if (IShape == aData->IndexOfS2())
        {
          const occ::handle<Geom2d_Curve>& aPCurve = aData->InterferenceOnS2().PCurveOnFace();
          if (aPCurve.IsNull())
          {
            continue;
          }

          anOtherPCurve.Load(aPCurve,
                             aData->InterferenceOnS2().FirstParameter(),
                             aData->InterferenceOnS2().LastParameter());
        }
        else
        {
          // Normal case - no common surface
          continue;
        }
        if (IsEqual(anOtherPCurve.LastParameter(), anOtherPCurve.FirstParameter()))
        {
          // Degenerates
          continue;
        }
        anIntersector.Perform(aCorkPCurve, anOtherPCurve, tol2d, Precision::PConfusion());
        if (anIntersector.NbSegments() > 0 || anIntersector.NbPoints() > 0)
        {
          throw StdFail_NotDone("OneCorner : fillets have too big radiuses");
        }
      }
    }
    NCollection_List<occ::handle<TopOpeBRepDS_Interference>>::Iterator anIter(
      DStr.ChangeShapeInterferences(IShape));
    for (; anIter.More(); anIter.Next())
    {
      occ::handle<TopOpeBRepDS_SurfaceCurveInterference> anOtherIntrf =
        occ::down_cast<TopOpeBRepDS_SurfaceCurveInterference>(anIter.Value());
      // We need only interferences between cork face and curves
      // of intersection with another fillet surfaces
      if (anOtherIntrf.IsNull())
      {
        continue;
      }
      // Look if there is an intersection between pcurves
      occ::handle<Geom_TrimmedCurve> anOtherCur =
        occ::down_cast<Geom_TrimmedCurve>(DStr.Curve(anOtherIntrf->Geometry()).Curve());
      if (anOtherCur.IsNull())
      {
        continue;
      }
      Geom2dAdaptor_Curve anOtherPCurve(anOtherIntrf->PCurve(),
                                        anOtherCur->FirstParameter(),
                                        anOtherCur->LastParameter());
      anIntersector.Perform(aCorkPCurve, anOtherPCurve, tol2d, Precision::PConfusion());
      if (anIntersector.NbSegments() > 0 || anIntersector.NbPoints() > 0)
      {
        throw StdFail_NotDone("OneCorner : fillets have too big radiuses");
      }
    }
    // 31/01/02 akm ^^^
    DStr.ChangeShapeInterferences(IShape).Append(Interfc);
    //// modified by jgv, 26.03.02 for OCC32 ////
    ChFiDS_CommonPoint CV[2];
    CV[0] = CV1;
    CV[1] = CV2;
    for (int i = 0; i < 2; i++)
    {
      if (CV[i].IsOnArc() && ChFi3d_IsPseudoSeam(CV[i].Arc(), Fv))
      {
        gp_Pnt2d                  pfac1, PcF, PcL;
        gp_Vec2d                  DerPc, DerHc;
        double                    first, last, prm1, prm2;
        bool                      onfirst, FirstToPar;
        occ::handle<Geom2d_Curve> Hc = BRep_Tool::CurveOnSurface(CV[i].Arc(), Fv, first, last);
        if (Hc.IsNull())
        {
          throw Standard_ConstructionError("Failed to get p-curve of edge");
        }
        pfac1   = Hc->Value(CV[i].ParameterOnArc());
        PcF     = Pc->Value(Udeb);
        PcL     = Pc->Value(Ufin);
        onfirst = pfac1.Distance(PcF) < pfac1.Distance(PcL);
        if (onfirst)
        {
          Pc->D1(Udeb, PcF, DerPc);
        }
        else
        {
          Pc->D1(Ufin, PcL, DerPc);
          DerPc.Reverse();
        }
        Hc->D1(CV[i].ParameterOnArc(), pfac1, DerHc);
        if (DerHc.Dot(DerPc) > 0.)
        {
          prm1       = CV[i].ParameterOnArc();
          prm2       = last;
          FirstToPar = false;
        }
        else
        {
          prm1       = first;
          prm2       = CV[i].ParameterOnArc();
          FirstToPar = true;
        }
        occ::handle<Geom_Curve> Ct = BRep_Tool::Curve(CV[i].Arc(), first, last);
        Ct                         = new Geom_TrimmedCurve(Ct, prm1, prm2);
        double                                           toled = BRep_Tool::Tolerance(CV[i].Arc());
        TopOpeBRepDS_Curve                               tcurv(Ct, toled);
        occ::handle<TopOpeBRepDS_CurvePointInterference> Interfp1, Interfp2;
        int                                              indcurv;
        indcurv       = DStr.AddCurve(tcurv);
        int indpoint  = (isfirst) ? stripe->IndexFirstPointOnS1() : stripe->IndexLastPointOnS1();
        int indvertex = DStr.AddShape(Vtx);
        if (FirstToPar)
        {
          Interfp1 = ChFi3d_FilPointInDS(TopAbs_FORWARD, indcurv, indvertex, prm1, true);
          Interfp2 = ChFi3d_FilPointInDS(TopAbs_REVERSED, indcurv, indpoint, prm2, false);
        }
        else
        {
          Interfp1 = ChFi3d_FilPointInDS(TopAbs_FORWARD, indcurv, indpoint, prm1, false);
          Interfp2 = ChFi3d_FilPointInDS(TopAbs_REVERSED, indcurv, indvertex, prm2, true);
        }
        DStr.ChangeCurveInterferences(indcurv).Append(Interfp1);
        DStr.ChangeCurveInterferences(indcurv).Append(Interfp2);
        int indface = DStr.AddShape(Fv);
        Interfc     = ChFi3d_FilCurveInDS(indcurv, indface, Hc, CV[i].Arc().Orientation());
        DStr.ChangeShapeInterferences(indface).Append(Interfc);
        TopoDS_Edge aLocalEdge = CV[i].Arc();
        aLocalEdge.Reverse();
        occ::handle<Geom2d_Curve> HcR = BRep_Tool::CurveOnSurface(aLocalEdge, Fv, first, last);
        if (HcR.IsNull())
        {
          throw Standard_ConstructionError("Failed to get p-curve of edge");
        }
        Interfc = ChFi3d_FilCurveInDS(indcurv, indface, HcR, aLocalEdge.Orientation());
        DStr.ChangeShapeInterferences(indface).Append(Interfc);
        // modify degenerated edge
        bool            DegenExist = false;
        TopoDS_Edge     Edeg;
        TopExp_Explorer Explo(Fv, TopAbs_EDGE);
        for (; Explo.More(); Explo.Next())
        {
          const TopoDS_Edge& Ecur = TopoDS::Edge(Explo.Current());
          if (BRep_Tool::Degenerated(Ecur))
          {
            TopoDS_Vertex Vf, Vl;
            TopExp::Vertices(Ecur, Vf, Vl);
            if (Vf.IsSame(Vtx) || Vl.IsSame(Vtx))
            {
              DegenExist = true;
              Edeg       = Ecur;
              break;
            }
          }
        }
        if (DegenExist)
        {
          double                    fd, ld;
          occ::handle<Geom2d_Curve> Cd = BRep_Tool::CurveOnSurface(Edeg, Fv, fd, ld);
          if (Cd.IsNull())
          {
            throw Standard_ConstructionError("Failed to get p-curve of edge");
          }
          occ::handle<Geom2d_TrimmedCurve> tCd = occ::down_cast<Geom2d_TrimmedCurve>(Cd);
          if (!tCd.IsNull())
          {
            Cd = tCd->BasisCurve();
          }
          gp_Pnt2d                      P2d = (FirstToPar) ? Hc->Value(first) : Hc->Value(last);
          Geom2dAPI_ProjectPointOnCurve Projector(P2d, Cd);
          double                        par  = Projector.LowerDistanceParameter();
          int                           Ideg = DStr.AddShape(Edeg);
          // clang-format off
          TopAbs_Orientation ori = (par < fd)? TopAbs_FORWARD : TopAbs_REVERSED; //if par<fd => par>ld
          // clang-format on
          Interfp1 = ChFi3d_FilPointInDS(ori, Ideg, indvertex, par, true);
          DStr.ChangeShapeInterferences(Ideg).Append(Interfp1);
        }
      }
    }
    /////////////////////////////////////////////
    stripe->ChangePCurve(isfirst) = Ps;
    stripe->SetCurve(ICurve, isfirst);
    stripe->SetParameters(isfirst, Udeb, Ufin);
  }
  else
  {
    // curves curv1 are curv2 stored in the DS
    // these curves will not be reconstructed by FILDS as
    // one places stripe->InDS(isfirst);

    // interferences of curv1 and curv2 on Fv
    ComputeCurve2d(curv1, Fv, c2d1);
    occ::handle<TopOpeBRepDS_SurfaceCurveInterference> InterFv;
    InterFv = ChFi3d_FilCurveInDS(Icurv1, IShape, c2d1, Et);
    DStr.ChangeShapeInterferences(IShape).Append(InterFv);
    ComputeCurve2d(curv2, Fv, c2d2);
    InterFv = ChFi3d_FilCurveInDS(Icurv2, IShape, c2d2, Et);
    DStr.ChangeShapeInterferences(IShape).Append(InterFv);
    // interferences of curv1 and curv2 on Isurf
    if (Fd->Orientation() == Fv.Orientation())
    {
      Et = TopAbs::Reverse(Et);
    }
    c2d1    = new Geom2d_TrimmedCurve(Ps, Udeb, par2);
    InterFv = ChFi3d_FilCurveInDS(Icurv1, Isurf, c2d1, Et);
    DStr.ChangeSurfaceInterferences(Isurf).Append(InterFv);
    c2d2    = new Geom2d_TrimmedCurve(Ps, par2, Ufin);
    InterFv = ChFi3d_FilCurveInDS(Icurv2, Isurf, c2d2, Et);
    DStr.ChangeSurfaceInterferences(Isurf).Append(InterFv);

    // limitation of the sewing edge
    int                                              Iarc = DStr.AddShape(edgecouture);
    occ::handle<TopOpeBRepDS_CurvePointInterference> Interfedge;
    TopAbs_Orientation                               ori;
    TopoDS_Vertex                                    Vdeb, Vfin;
    Vdeb = TopExp::FirstVertex(edgecouture);
    Vfin = TopExp::LastVertex(edgecouture);
    double pard, parf;
    pard = BRep_Tool::Parameter(Vdeb, edgecouture);
    parf = BRep_Tool::Parameter(Vfin, edgecouture);
    if (std::abs(par1 - pard) < std::abs(parf - par1))
    {
      ori = TopAbs_FORWARD;
    }
    else
    {
      ori = TopAbs_REVERSED;
    }
    Interfedge = ChFi3d_FilPointInDS(ori, Iarc, indpt, par1);
    DStr.ChangeShapeInterferences(Iarc).Append(Interfedge);

    // creation of CurveInterferences from Icurv1 and Icurv2
    stripe->InDS(isfirst);
    int                                              ind1 = stripe->IndexPoint(isfirst, 1);
    int                                              ind2 = stripe->IndexPoint(isfirst, 2);
    occ::handle<TopOpeBRepDS_CurvePointInterference> interfprol =
      ChFi3d_FilPointInDS(TopAbs_FORWARD, Icurv1, ind1, Udeb);
    DStr.ChangeCurveInterferences(Icurv1).Append(interfprol);
    interfprol = ChFi3d_FilPointInDS(TopAbs_REVERSED, Icurv1, indpt, par2);
    DStr.ChangeCurveInterferences(Icurv1).Append(interfprol);
    interfprol = ChFi3d_FilPointInDS(TopAbs_FORWARD, Icurv2, indpt, par2);
    DStr.ChangeCurveInterferences(Icurv2).Append(interfprol);
    interfprol = ChFi3d_FilPointInDS(TopAbs_REVERSED, Icurv2, ind2, Ufin);
    DStr.ChangeCurveInterferences(Icurv2).Append(interfprol);
  }

  if (FvT.IsNull() && onsOS == 0)
  {
    ChFi3d_EnlargeBox(HBs, Pc, Udeb, Ufin, box1, box2);
  }

  // The fillet's line on a face running from vertex to vertex along the
  // split toward the piece the spine is on: the split goes with the piece.
  for (int ons = 1; ons <= 2 && inters; ons++)
  {
    const ChFiDS_CommonPoint& aCP    = Fd->Vertex(isfirst, ons);
    const ChFiDS_CommonPoint& aCPEnd = Fd->Vertex(!isfirst, ons);
    if (!aCP.IsVertex() || !aCPEnd.IsVertex())
    {
      continue;
    }
    const TopoDS_Edge aSplit = SplitUnderLine(TopoDS::Face(DStr.Shape(Fd->Index(ons))),
                                              aCP.Vertex(),
                                              aCPEnd.Vertex(),
                                              Fd->Interference(ons),
                                              spine,
                                              myEFMap,
                                              10 * tolapp3d);
    if (!aSplit.IsNull())
    {
      Swallowed.Append(aSplit);
      Swallowed.Append(aCP.Vertex());
    }
  }
  // The edges under the fillet go the way of the spine's edge: cut away
  // from their end toward Vtx.
  for (NCollection_List<TopoDS_Shape>::Iterator itS(Swallowed); itS.More() && inters; itS.Next())
  {
    const TopoDS_Edge&   anE = TopoDS::Edge(itS.Value());
    itS.Next();
    const TopoDS_Vertex& aV  = TopoDS::Vertex(itS.Value());
    TopAbs_Orientation   aOV = TopAbs_FORWARD;
    for (ex.Init(anE.Oriented(TopAbs_FORWARD), TopAbs_VERTEX); ex.More(); ex.Next())
    {
      if (aV.IsSame(ex.Current()))
      {
        aOV = ex.Current().Orientation();
        break;
      }
    }
    DStr.ChangeShapeInterferences(DStr.AddShape(anE))
      .Append(ChFi3d_FilVertexInDS(TopAbs::Reverse(aOV),
                                   DStr.AddShape(anE),
                                   DStr.AddShape(aV),
                                   BRep_Tool::Parameter(aV, anE)));
  }

  if (onsame && inters)
  {
// VARIANT 1:
// A small missing end of curve is added for the extension
// of the face at end and the limitation of the opposing face.

//   VARIANT 2 : extend Arcprol, not create new small edge
//   To do: modify for intcouture
#define VARIANT1

    // First of all the points are cut with the edge of the spine.
    int                IArcspine = DStr.AddShape(Arcspine);
    int                IVtx      = DStr.AddShape(Vtx);
    TopAbs_Orientation OVtx      = TopAbs_FORWARD;
    for (ex.Init(Arcspine.Oriented(TopAbs_FORWARD), TopAbs_VERTEX); ex.More(); ex.Next())
    {
      if (Vtx.IsSame(ex.Current()))
      {
        OVtx = ex.Current().Orientation();
        break;
      }
    }
    OVtx                                                    = TopAbs::Reverse(OVtx);
    double                                           parVtx = BRep_Tool::Parameter(Vtx, Arcspine);
    occ::handle<TopOpeBRepDS_CurvePointInterference> interfv =
      ChFi3d_FilVertexInDS(OVtx, IArcspine, IVtx, parVtx);
    DStr.ChangeShapeInterferences(IArcspine).Append(interfv);
    // Now the missing curves are constructed.
    TopoDS_Vertex V2;
    for (ex.Init(Arcprol.Oriented(TopAbs_FORWARD), TopAbs_VERTEX); ex.More(); ex.Next())
    {
      if (Vtx.IsSame(ex.Current()))
      {
        OVtx = ex.Current().Orientation();
      }
      else
      {
        V2 = TopoDS::Vertex(ex.Current());
      }
    }

    occ::handle<Geom2d_Curve> Hc;
#ifdef VARIANT1
    parVtx = BRep_Tool::Parameter(Vtx, Arcprol);
#else
    parVtx = BRep_Tool::Parameter(V2, Arcprol);
#endif
    const ChFiDS_FaceInterference& Fiop = Fd->Interference(IFopArc);
    gp_Pnt2d                       pop1, pop2, pv1, pv2;
    Hc = PCurveInFace(Arcprol, Fop, Ubid, Ubid);
    if (Hc.IsNull())
    {
      throw Standard_ConstructionError("Failed to get p-curve of edge");
    }
    pop1 = Hc->Value(parVtx);
    pop2 = Fiop.PCurveOnFace()->Value(Fiop.Parameter(isfirst));
    Hc   = PCurveInFace(Arcprol, Fv, Ubid, Ubid);
    if (Hc.IsNull())
    {
      throw Standard_ConstructionError("Failed to get p-curve of edge");
    }
    pv1 = Hc->Value(parVtx);
    // The extension starts where the cut ends, on the cut's period of a
    // periodic Fv: Vtx on Fv's seam has a parameter on either side of it,
    // and Arcprol's there need not be the cut's.
    pv2 = aCutOnOp;
    ChFi3d_Recale(Bs, pv1, pv2, false);
    NCollection_Array1<double> Pardeb(1, 4), Parfin(1, 4);
    Pardeb(1) = pop1.X();
    Pardeb(2) = pop1.Y();
    Pardeb(3) = pv1.X();
    Pardeb(4) = pv1.Y();
    Parfin(1) = pop2.X();
    Parfin(2) = pop2.Y();
    Parfin(3) = pv2.X();
    Parfin(4) = pv2.Y();
    double uu1, uu2, vv1, vv2;
    ChFi3d_Boite(pv1, pv2, uu1, uu2, vv1, vv2);
    ChFi3d_BoundFac(Bs, uu1, uu2, vv1, vv2);
    ChFi3d_Boite(pop1, pop2, uu1, uu2, vv1, vv2);
    ChFi3d_BoundFac(Bop, uu1, uu2, vv1, vv2);

    occ::handle<Geom_Curve>   zob3d;
    occ::handle<Geom2d_Curve> zob2dop, zob2dv;
    // double tolreached;
    if (!ChFi3d_ComputeCurves(HBop,
                              HBs,
                              Pardeb,
                              Parfin,
                              zob3d,
                              zob2dop,
                              zob2dv,
                              tolapp3d,
                              tol2d,
                              tolreached))
    {
      throw Standard_Failure("OneCorner : echec calcul intersection");
    }

    Udeb = zob3d->FirstParameter();
    Ufin = zob3d->LastParameter();
    TopOpeBRepDS_Curve Zob(zob3d, tolreached);
    int                IZob = DStr.AddCurve(Zob);

    // it is determined if Fop has an edge of sewing
    // it is determined if the curve has an intersection with the edge of sewing

    // TopoDS_Edge edgecouture;
    // bool couture;
    ChFi3d_Couture(Fop, couture, edgecouture);

    if (couture && !BRep_Tool::Degenerated(edgecouture))
    {
      BRepLib_MakeEdge  Bedge(zob3d);
      TopoDS_Edge       edg = Bedge.Edge();
      BRepExtrema_ExtCC extCC(edgecouture, edg);
      if (extCC.IsDone() && extCC.NbExt() != 0)
      {
        for (int i = 1; i <= extCC.NbExt() && !intcouture; i++)
        {
          if (extCC.SquareDistance(i) <= 1.e-8)
          {
            par1                  = extCC.ParameterOnE1(i);
            par2                  = extCC.ParameterOnE2(i);
            gp_Pnt             P1 = extCC.PointOnE1(i);
            TopOpeBRepDS_Point tpoint(P1, 1.e-4);
            indpt      = DStr.AddPoint(tpoint);
            intcouture = true;
            curv1      = new Geom_TrimmedCurve(zob3d, Udeb, par2);
            curv2      = new Geom_TrimmedCurve(zob3d, par2, Ufin);
            TopOpeBRepDS_Curve tcurv1(curv1, tolreached);
            TopOpeBRepDS_Curve tcurv2(curv2, tolreached);
            Icurv1 = DStr.AddCurve(tcurv1);
            Icurv2 = DStr.AddCurve(tcurv2);
          }
        }
      }
    }
    if (intcouture)
    {

      // interference of curv1 and curv2 on Ishape
      Et = TopAbs::Reverse(TopAbs::Compose(OVtx, OArcprolv));
      ComputeCurve2d(curv1, Fop, c2d1);
      occ::handle<TopOpeBRepDS_SurfaceCurveInterference> InterFv =
        ChFi3d_FilCurveInDS(Icurv1, IShape, /*zob2dv*/ c2d1, Et);
      DStr.ChangeShapeInterferences(IShape).Append(InterFv);
      ComputeCurve2d(curv2, Fop, c2d2);
      InterFv = ChFi3d_FilCurveInDS(Icurv2, IShape, /*zob2dv*/ c2d2, Et);
      DStr.ChangeShapeInterferences(IShape).Append(InterFv);

      // limitation of the sewing edge
      int                                              Iarc = DStr.AddShape(edgecouture);
      occ::handle<TopOpeBRepDS_CurvePointInterference> Interfedge;
      TopAbs_Orientation                               ori;
      TopoDS_Vertex                                    Vdeb, Vfin;
      Vdeb = TopExp::FirstVertex(edgecouture);
      Vfin = TopExp::LastVertex(edgecouture);
      double pard, parf;
      pard = BRep_Tool::Parameter(Vdeb, edgecouture);
      parf = BRep_Tool::Parameter(Vfin, edgecouture);
      if (std::abs(par1 - pard) < std::abs(parf - par1))
      {
        ori = TopAbs_REVERSED;
      }
      else
      {
        ori = TopAbs_FORWARD;
      }
      Interfedge = ChFi3d_FilPointInDS(ori, Iarc, indpt, par1);
      DStr.ChangeShapeInterferences(Iarc).Append(Interfedge);

      //  interference of curv1 and curv2 on Iop
      int Iop = DStr.AddShape(Fop);
      Et      = TopAbs::Reverse(TopAbs::Compose(OVtx, OArcprolop));
      occ::handle<TopOpeBRepDS_SurfaceCurveInterference> Interfop;
      ComputeCurve2d(curv1, Fop, c2d1);
      Interfop = ChFi3d_FilCurveInDS(Icurv1, Iop, c2d1, Et);
      DStr.ChangeShapeInterferences(Iop).Append(Interfop);
      ComputeCurve2d(curv2, Fop, c2d2);
      Interfop = ChFi3d_FilCurveInDS(Icurv2, Iop, c2d2, Et);
      DStr.ChangeShapeInterferences(Iop).Append(Interfop);
      occ::handle<TopOpeBRepDS_CurvePointInterference> interfprol =
        ChFi3d_FilVertexInDS(TopAbs_FORWARD, Icurv1, IVtx, Udeb);
      DStr.ChangeCurveInterferences(Icurv1).Append(interfprol);
      interfprol = ChFi3d_FilPointInDS(TopAbs_REVERSED, Icurv1, indpt, par2);
      DStr.ChangeCurveInterferences(Icurv1).Append(interfprol);
      int icc    = stripe->IndexPoint(isfirst, IFopArc);
      interfprol = ChFi3d_FilPointInDS(TopAbs_FORWARD, Icurv2, indpt, par2);
      DStr.ChangeCurveInterferences(Icurv2).Append(interfprol);
      interfprol = ChFi3d_FilPointInDS(TopAbs_REVERSED, Icurv2, icc, Ufin);
      DStr.ChangeCurveInterferences(Icurv2).Append(interfprol);
    }
    else
    {
      Et = TopAbs::Reverse(TopAbs::Compose(OVtx, OArcprolv));
      occ::handle<TopOpeBRepDS_SurfaceCurveInterference> InterFv =
        ChFi3d_FilCurveInDS(IZob, IShape, zob2dv, Et);
      DStr.ChangeShapeInterferences(IShape).Append(InterFv);
      Et = TopAbs::Reverse(TopAbs::Compose(OVtx, OArcprolop));
      if (zobOnEtan)
      {
        // The extension runs along Etan, which the fillet cuts at its end:
        // it bounds Fprol, not Fop.
        TopoDS_Face aFprolF = Fprol;
        aFprolF.Orientation(TopAbs_FORWARD);
        for (ex.Init(aFprolF, TopAbs_EDGE); ex.More(); ex.Next())
        {
          if (Arcprol.IsSame(ex.Current()))
          {
            Et = TopAbs::Reverse(TopAbs::Compose(OVtx, ex.Current().Orientation()));
            break;
          }
        }
        zob2dop = GeomProjLib::Curve2d(zob3d, Udeb, Ufin, BRep_Tool::Surface(aFprolF));
        // On the period of Fprol's domain at Vtx: a projection on a periodic
        // surface (a cylinder) lands on the surface's first period, which
        // need not be the face's.
        if (!zob2dop.IsNull())
        {
          const occ::handle<Geom_Surface> aSprol = BRep_Tool::Surface(aFprolF);
          const gp_Pnt2d                  aUVv  = BRep_Tool::Parameters(Vtx, aFprolF);
          const gp_Pnt2d                  aUVz  = zob2dop->Value(Udeb);
          double                          aDu = 0., aDv = 0.;
          if (aSprol->IsUPeriodic())
          {
            const double aPer = aSprol->UPeriod();
            aDu               = aPer * std::floor((aUVv.X() - aUVz.X()) / aPer + 0.5);
          }
          if (aSprol->IsVPeriodic())
          {
            const double aPer = aSprol->VPeriod();
            aDv               = aPer * std::floor((aUVv.Y() - aUVz.Y()) / aPer + 0.5);
          }
          if (aDu != 0. || aDv != 0.)
          {
            zob2dop->Translate(gp_Vec2d(aDu, aDv));
          }
        }
      }
      int Iop = DStr.AddShape(zobOnEtan ? Fprol : Fop);
      occ::handle<TopOpeBRepDS_SurfaceCurveInterference> Interfop =
        ChFi3d_FilCurveInDS(IZob, Iop, zob2dop, Et);
      DStr.ChangeShapeInterferences(Iop).Append(Interfop);
      occ::handle<TopOpeBRepDS_CurvePointInterference> interfprol;
#ifdef VARIANT1
      interfprol = ChFi3d_FilVertexInDS(TopAbs_FORWARD, IZob, IVtx, Udeb);
#else
      {
        int IV2    = DStr.AddShape(V2); // VARIANT 2
        interfprol = ChFi3d_FilVertexInDS(TopAbs_FORWARD, IZob, IV2, Udeb);
      }
#endif
      DStr.ChangeCurveInterferences(IZob).Append(interfprol);
      int icc    = stripe->IndexPoint(isfirst, IFopArc);
      interfprol = ChFi3d_FilPointInDS(TopAbs_REVERSED, IZob, icc, Ufin);
      DStr.ChangeCurveInterferences(IZob).Append(interfprol);
#ifdef VARIANT1
      {
        if (IFopArc == 1)
        {
          box1.Add(zob3d->Value(Ufin));
        }
        else
        {
          box2.Add(zob3d->Value(Ufin));
        }
      }
#else
      {
        // cut off existing Arcprol
        int iArcprol = DStr.AddShape(Arcprol);
        interfprol   = ChFi3d_FilPointInDS(OVtx, iArcprol, icc, Udeb);
        DStr.ChangeShapeInterferences(Arcprol).Append(interfprol);
      }
#endif
    }
  }
  ChFi3d_EnlargeBox(DStr, stripe, Fd, box1, box2, isfirst);
  if (CV1.IsOnArc())
  {
    ChFi3d_EnlargeBox(CV1.Arc(), myEFMap(CV1.Arc()), CV1.ParameterOnArc(), box1);
  }
  if (CV2.IsOnArc())
  {
    ChFi3d_EnlargeBox(CV2.Arc(), myEFMap(CV2.Arc()), CV2.ParameterOnArc(), box2);
  }
  if (!CV1.IsVertex())
  {
    ChFi3d_SetPointTolerance(DStr, box1, stripe->IndexPoint(isfirst, 1));
  }
  if (!CV2.IsVertex())
  {
    ChFi3d_SetPointTolerance(DStr, box2, stripe->IndexPoint(isfirst, 2));
  }
  if (!FvT.IsNull())
  {
    ChFi3d_SetPointTolerance(DStr, boxT, indpt);
  }

#ifdef OCCT_DEBUG
  ChFi3d_ResultChron(ch, t_sameinter); // result perf condition if (same &&inter)
#endif
}

//=======================================================================
// function : cherche_face
// purpose  : find face F belonging to the map, different from faces
//           F1  F2 F3 and containing edge E
//=======================================================================

static void cherche_face(const NCollection_List<TopoDS_Shape>& map,
                         const TopoDS_Edge&                    E,
                         const TopoDS_Face&                    F1,
                         const TopoDS_Face&                    F2,
                         const TopoDS_Face&                    F3,
                         TopoDS_Face&                          F)
{
  TopoDS_Face                              Fcur;
  bool                                     trouve = false;
  NCollection_List<TopoDS_Shape>::Iterator It;
  int                                      ie;
  for (It.Initialize(map); It.More() && !trouve; It.Next())
  {
    Fcur = TopoDS::Face(It.Value());
    if (!Fcur.IsSame(F1) && !Fcur.IsSame(F2) && !Fcur.IsSame(F3))
    {
      NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> MapE;
      TopExp::MapShapes(Fcur, TopAbs_EDGE, MapE);
      for (ie = 1; ie <= MapE.Extent() && !trouve; ie++)
      {
        TopoDS_Shape aLocalShape = TopoDS_Shape(MapE(ie));
        if (E.IsSame(TopoDS::Edge(aLocalShape)))
        //            if (E.IsSame(TopoDS::Edge(TopoDS_Shape (MapE(ie)))))
        {
          F      = Fcur;
          trouve = true;
        }
      }
    }
  }
  if (F.IsNull())
  {
    throw Standard_ConstructionError("Failed to find face.");
  }
}

//=======================================================================
// function : cherche_edge1
// purpose  : find common edge between faces F1 and F2
//=======================================================================

static void cherche_edge1(const TopoDS_Face& F1, const TopoDS_Face& F2, TopoDS_Edge& Edge)
{
  int                                                           i, j;
  TopoDS_Edge                                                   Ecur1, Ecur2;
  bool                                                          trouve = false;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> MapE1, MapE2;
  TopExp::MapShapes(F1, TopAbs_EDGE, MapE1);
  TopExp::MapShapes(F2, TopAbs_EDGE, MapE2);
  for (i = 1; i <= MapE1.Extent() && !trouve; i++)
  {
    TopoDS_Shape aLocalShape = TopoDS_Shape(MapE1(i));
    Ecur1                    = TopoDS::Edge(aLocalShape);
    //	Ecur1=TopoDS::Edge(TopoDS_Shape (MapE1(i)));
    for (j = 1; j <= MapE2.Extent() && !trouve; j++)
    {
      aLocalShape = TopoDS_Shape(MapE2(j));
      Ecur2       = TopoDS::Edge(aLocalShape);
      //	      Ecur2=TopoDS::Edge(TopoDS_Shape (MapE2(j)));
      if (Ecur2.IsSame(Ecur1))
      {
        Edge   = Ecur1;
        trouve = true;
      }
    }
  }
  if (Edge.IsNull())
  {
    throw Standard_ConstructionError("Failed to find edge.");
  }
}

//=======================================================================
// function : containV
// purpose  : return true if vertex V belongs to F1
//=======================================================================

static bool containV(const TopoDS_Face& F1, const TopoDS_Vertex& V)
{
  int                                                           i;
  TopoDS_Vertex                                                 Vcur;
  bool                                                          trouve  = false;
  bool                                                          contain = false;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> MapV;
  TopExp::MapShapes(F1, TopAbs_VERTEX, MapV);
  for (i = 1; i <= MapV.Extent() && !trouve; i++)
  {
    TopoDS_Shape aLocalShape = TopoDS_Shape(MapV(i));
    Vcur                     = TopoDS::Vertex(aLocalShape);
    //	Vcur=TopoDS::Vertex(TopoDS_Shape (MapV(i)));
    if (Vcur.IsSame(V))
    {
      contain = true;
      trouve  = true;
    }
  }
  return contain;
}

//=======================================================================
// function : containE
// purpose  : return true if edge E belongs to F1
//=======================================================================

static bool containE(const TopoDS_Face& F1, const TopoDS_Edge& E)
{
  int                                                           i;
  TopoDS_Edge                                                   Ecur;
  bool                                                          trouve  = false;
  bool                                                          contain = false;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> MapE;
  TopExp::MapShapes(F1, TopAbs_EDGE, MapE);
  for (i = 1; i <= MapE.Extent() && !trouve; i++)
  {
    TopoDS_Shape aLocalShape = TopoDS_Shape(MapE(i));
    Ecur                     = TopoDS::Edge(aLocalShape);
    //	Ecur=TopoDS::Edge(TopoDS_Shape (MapE(i)));
    if (Ecur.IsSame(E))
    {
      contain = true;
      trouve  = true;
    }
  }
  return contain;
}

//=======================================================================
// function : IsShrink
// purpose  : check if U (if <isU>==True) or V of points of <PC> is within
//           <tol> from <Param>, check points between <Pf> and <Pl>
//=======================================================================

static bool IsShrink(const Geom2dAdaptor_Curve& PC,
                     const double               Pf,
                     const double               Pl,
                     const double               Param,
                     const bool                 isU,
                     const double               tol)
{
  switch (PC.GetType())
  {
    case GeomAbs_Line: {
      gp_Pnt2d P1 = PC.Value(Pf);
      gp_Pnt2d P2 = PC.Value(Pl);
      return std::abs(P1.Coord(isU ? 1 : 2) - Param) <= tol
             && std::abs(P2.Coord(isU ? 1 : 2) - Param) <= tol;
    }
    case GeomAbs_BezierCurve:
    case GeomAbs_BSplineCurve: {
      math_FunctionSample aSample(Pf, Pl, 10);
      int                 i;
      for (i = 1; i <= aSample.NbPoints(); i++)
      {
        gp_Pnt2d P = PC.Value(aSample.GetParameter(i));
        if (std::abs(P.Coord(isU ? 1 : 2) - Param) > tol)
        {
          return false;
        }
      }
      return true;
    }
    default:;
  }
  return false;
}

//=======================================================================
// function : SectionCrossesEndFace
// purpose  : One side of the fillet ends at <P1> on <Arcpiv>, an edge of
//           <Vtx> between <Fad> and the face at end <Fv>; the other side
//           runs on past <Vtx> on its face, and is at <P2> in the same
//           section. IntersectMoreCorner extends <Fv> to the section's
//           curve: the piece added is the corner of the extended surface
//           between <Vtx>, <P1> and <P2>, and holds the section. True
//           when <Fv>'s other edge at <Vtx> runs into that corner -- a
//           wall stands on <Fv> there (the end of an arm fused to a taller
//           block): the section crosses <Fv>'s own edge, and what lies
//           beyond that edge faces away from <Fv>.
//=======================================================================

static bool SectionCrossesEndFace(const TopoDS_Vertex& Vtx,
                                  const TopoDS_Edge&   Arcpiv,
                                  const TopoDS_Face&   Fad,
                                  const gp_Pnt&        P1,
                                  const gp_Pnt&        P2,
                                  const ChFiDS_Map&    VEMap,
                                  const ChFiDS_Map&    EFMap)
{
  TopoDS_Face Fv;
  for (NCollection_List<TopoDS_Shape>::Iterator It(EFMap(Arcpiv)); It.More(); It.Next())
  {
    if (!Fad.IsSame(It.Value()))
    {
      Fv = TopoDS::Face(It.Value());
      break;
    }
  }
  if (Fv.IsNull())
  {
    return false;
  }
  TopoDS_Edge Arcprol;
  for (NCollection_List<TopoDS_Shape>::Iterator It(VEMap(Vtx)); It.More(); It.Next())
  {
    const TopoDS_Edge& E = TopoDS::Edge(It.Value());
    if (!E.IsSame(Arcpiv) && !BRep_Tool::Degenerated(E) && containE(Fv, E))
    {
      Arcprol = E;
      break;
    }
  }
  if (Arcprol.IsNull())
  {
    return false;
  }
  // Fv's normal at Vtx, and Arcprol's tangent there, away from Vtx
  const gp_Pnt2d      uv = BRep_Tool::Parameters(Vtx, Fv);
  BRepAdaptor_Surface Sv(Fv, false);
  gp_Pnt              PV;
  gp_Vec              DU, DV;
  Sv.D1(uv.X(), uv.Y(), PV, DU, DV);
  gp_Vec N = DU.Crossed(DV);
  if (N.Magnitude() <= gp::Resolution())
  {
    return false;
  }
  N.Normalize();
  BRepAdaptor_Curve Cprol(Arcprol);
  gp_Pnt            PE;
  gp_Vec            D;
  const double      t = BRep_Tool::Parameter(Vtx, Arcprol);
  Cprol.D1(t, PE, D);
  if (std::abs(t - Cprol.LastParameter()) < std::abs(t - Cprol.FirstParameter()))
  {
    D.Reverse();
  }
  // the three directions from Vtx, in Fv's tangent plane
  const gp_Pnt P = BRep_Tool::Pnt(Vtx);
  gp_Vec       A(P, P1), B(P, P2);
  A -= N * A.Dot(N);
  B -= N * B.Dot(N);
  D -= N * D.Dot(N);
  const double sAB = A.Crossed(B).Dot(N);
  const double tol = Precision::Angular() * A.Magnitude() * B.Magnitude();
  if (std::abs(sAB) <= tol || D.Magnitude() <= gp::Resolution())
  {
    return false;
  }
  return A.Crossed(D).Dot(N) * sAB > 0. && D.Crossed(B).Dot(N) * sAB > 0.;
}

//=======================================================================
// function : ArcPastSplit
// purpose  : The side of the fillet on <Fs> ends on <Arc>, an edge of <Fs>
//           and <Fe> away from <Vtx>, because it crossed into <Fs> before
//           its end: <Fs> is a piece of the same plane as <Fn>, a face of
//           <Vtx>, beside it across an edge at <Arc>'s vertex V, and <Arc>
//           carries on beyond V the edge <Eat> of <Vtx> between <Fn> and
//           <Fe> -- a wall kept in coplanar pieces (a body without Refine)
//           standing on <Fe>. On the wall in one face the side would end on
//           the edge through <Vtx>; the faces at the end are the same.
//=======================================================================

static bool ArcPastSplit(const TopoDS_Vertex& Vtx,
                         const TopoDS_Edge&   Arc,
                         const TopoDS_Face&   Fs,
                         const ChFiDS_Map&    VEMap,
                         const ChFiDS_Map&    EFMap,
                         const double         Tol,
                         TopoDS_Edge&         Eat)
{
  if (hasVertex(Arc, Vtx) || EFMap(Arc).Extent() != 2)
  {
    return false;
  }
  TopoDS_Face Fe;
  for (NCollection_List<TopoDS_Shape>::Iterator It(EFMap(Arc)); It.More(); It.Next())
  {
    if (!Fs.IsSame(It.Value()))
    {
      Fe = TopoDS::Face(It.Value());
    }
  }
  if (Fe.IsNull() || !containV(Fe, Vtx))
  {
    return false;
  }
  for (NCollection_List<TopoDS_Shape>::Iterator It(VEMap(Vtx)); It.More(); It.Next())
  {
    const TopoDS_Edge& anE = TopoDS::Edge(It.Value());
    if (anE.IsSame(Arc) || BRep_Tool::Degenerated(anE) || EFMap(anE).Extent() != 2
        || !containE(Fe, anE))
    {
      continue;
    }
    // V, the vertex of anE's other end, on Arc
    TopoDS_Vertex aV1, aV2;
    TopExp::Vertices(anE, aV1, aV2);
    const TopoDS_Vertex aV = Vtx.IsSame(aV1) ? aV2 : aV1;
    if (aV.IsSame(Vtx) || !hasVertex(Arc, aV))
    {
      continue;
    }
    TopoDS_Face Fn;
    for (NCollection_List<TopoDS_Shape>::Iterator itF(EFMap(anE)); itF.More(); itF.Next())
    {
      if (!Fe.IsSame(itF.Value()))
      {
        Fn = TopoDS::Face(itF.Value());
      }
    }
    if (Fn.IsNull() || Fn.IsSame(Fs) || !CoplanarPieces(Fs, Fn, Tol))
    {
      continue;
    }
    // the split: an edge of V between Fn and Fs
    for (NCollection_List<TopoDS_Shape>::Iterator itE(VEMap(aV)); itE.More(); itE.Next())
    {
      const TopoDS_Edge& aSplit = TopoDS::Edge(itE.Value());
      if (!aSplit.IsSame(anE) && !aSplit.IsSame(Arc) && containE(Fn, aSplit)
          && containE(Fs, aSplit))
      {
        Eat = anE;
        return true;
      }
    }
  }
  return false;
}

//=================================================================================================

void ChFi3d_Builder::PerformIntersectionAtEnd(const int Index)
{

  // intersection at end of fillet with at least two faces
  // process the following cases:
  // - top has n (n>3) adjacent edges
  // - top has 3 edges and fillet on one of edges touches
  //   more than one face

#ifdef OCCT_DEBUG
  OSD_Chronometer ch; // init perf
#endif

  TopOpeBRepDS_DataStructure&                            DStr = myDS->ChangeDS();
  const int                                              nn   = 15;
  NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator It;
  It.Initialize(myVDataMap(Index));
  occ::handle<ChFiDS_Stripe>                          stripe = It.Value();
  const occ::handle<ChFiDS_Spine>                     spine  = stripe->Spine();
  NCollection_Sequence<occ::handle<ChFiDS_SurfData>>& SeqFil =
    stripe->ChangeSetOfSurfData()->ChangeSequence();
  const TopoDS_Vertex& Vtx     = myVDataMap.FindKey(Index);
  int                  sens    = 0, num, num1;
  bool                 couture = false, isfirst;
  // int sense;
  TopoDS_Edge edgelibre1, edgelibre2, EdgeSpine;
  bool        bordlibre;
  // determine the number of faces and edges
  NCollection_Array1<TopoDS_Shape>         tabedg(0, nn);
  TopoDS_Face                              F1, F2;
  int                                      nface = ChFi3d_nbface(myVFMap(Vtx));
  NCollection_List<TopoDS_Shape>::Iterator ItF;
  int                                      nbarete;
  nbarete = ChFi3d_NbNotDegeneratedEdges(Vtx, myVEMap);
  ChFi3d_ChercheBordsLibres(myVEMap, Vtx, bordlibre, edgelibre1, edgelibre2);
  if (bordlibre)
  {
    nbarete = (nbarete - 2) / 2 + 2;
  }
  else
  {
    nbarete = nbarete / 2;
  }
  // it is determined if there is an edge of sewing and it face

  TopoDS_Face facecouture;
  TopoDS_Edge edgecouture;

  bool trouve = false;
  for (ItF.Initialize(myVFMap(Vtx)); ItF.More() && !couture; ItF.Next())
  {
    TopoDS_Face fcur = TopoDS::Face(ItF.Value());
    ChFi3d_CoutureOnVertex(fcur, Vtx, couture, edgecouture);
    if (couture)
    {
      facecouture = fcur;
    }
  }
  // it is determined if one of edges adjacent to the fillet is regular
  bool          reg1, reg2;
  TopoDS_Edge   Ecur, Eadj1, Eadj2;
  TopoDS_Face   Fga, Fdr;
  TopoDS_Vertex Vbid1;
  int           nbsurf, nbedge;
  reg1    = false;
  reg2    = false;
  nbsurf  = SeqFil.Length();
  nbedge  = spine->NbEdges();
  num     = ChFi3d_IndexOfSurfData(Vtx, stripe, sens);
  isfirst = (sens == 1);
  ChFiDS_State state;
  if (isfirst)
  {
    EdgeSpine = spine->Edges(1);
    num1      = num + 1;
    state     = spine->FirstStatus();
  }
  else
  {
    EdgeSpine = spine->Edges(nbedge);
    num1      = num - 1;
    state     = spine->LastStatus();
  }
  if (nbsurf != nbedge && nbsurf != 1)
  {
    ChFi3d_edge_common_faces(myEFMap(EdgeSpine), F1, F2);
    if (F1.IsSame(facecouture))
    {
      Eadj1 = edgecouture;
    }
    else
    {
      ChFi3d_cherche_element(Vtx, EdgeSpine, F1, Eadj1, Vbid1);
    }
    ChFi3d_edge_common_faces(myEFMap(Eadj1), Fga, Fdr);
    //  Modified by Sergey KHROMOV - Fri Dec 21 17:57:32 2001 Begin
    //  reg1=BRep_Tool::Continuity(Eadj1,Fga,Fdr)!=GeomAbs_C0;
    reg1 = ChFi3d::IsTangentFaces(Eadj1, Fga, Fdr);
    //  Modified by Sergey KHROMOV - Fri Dec 21 17:57:33 2001 End
    if (F2.IsSame(facecouture))
    {
      Eadj2 = edgecouture;
    }
    else
    {
      ChFi3d_cherche_element(Vtx, EdgeSpine, F2, Eadj2, Vbid1);
    }
    ChFi3d_edge_common_faces(myEFMap(Eadj2), Fga, Fdr);
    //  Modified by Sergey KHROMOV - Fri Dec 21 17:58:22 2001 Begin
    //  reg2=BRep_Tool::Continuity(Eadj2,Fga,Fdr)!=GeomAbs_C0;
    reg2 = ChFi3d::IsTangentFaces(Eadj2, Fga, Fdr);
    //  Modified by Sergey KHROMOV - Fri Dec 21 17:58:24 2001 End

    // two faces common to the edge are found
    if (reg1 || reg2)
    {
      bool               compoint1 = false;
      bool               compoint2 = false;
      ChFiDS_CommonPoint cp1, cp2;
      cp1 = SeqFil(num1)->ChangeVertex(isfirst, 1);
      cp2 = SeqFil(num1)->ChangeVertex(isfirst, 2);
      if (cp1.IsOnArc())
      {
        if (cp1.Arc().IsSame(Eadj1) || cp1.Arc().IsSame(Eadj2))
        {
          compoint1 = true;
        }
      }
      if (cp2.IsOnArc())
      {
        if (cp2.Arc().IsSame(Eadj1) || cp2.Arc().IsSame(Eadj2))
        {
          compoint2 = true;
        }
      }
      if (compoint1 && compoint2)
      {
        SeqFil.Remove(num);
        num = ChFi3d_IndexOfSurfData(Vtx, stripe, sens);
        if (isfirst)
        {
          num1 = num + 1;
        }
        else
        {
          num1 = num - 1;
        }
        reg1 = false;
        reg2 = false;
      }
    }
  }
  // One side of the fillet on an edge of Vtx, the other in a plane face at
  // the end: the line there, carried on, crosses into a coplanar piece of
  // the face and meets a face of Vtx on an edge of that piece (a wall kept
  // in pieces, see LineOverSplit). The line ends on that edge, as on the
  // wall in one face, and is carried over the split at the end. A walk past
  // the end that stopped on the split, its other side held at a point
  // (ChFi3d_Purge), is left out: the line is carried from the data before.
  int                    onsOS = 0;
  occ::handle<Geom_Line> LinOS;
  TopoDS_Edge            EsplitOS, EarcOS;
  TopoDS_Face            FnOS;
  double                 wQOS = 0., wPOS = 0., parQOS = 0., parPOS = 0.;
  // the edge of Vtx an arc away from it carries on (ArcPastSplit)
  TopoDS_Edge EdgePast;
  if (!couture && !bordlibre && !reg1 && !reg2 && nbarete == 4)
  {
    int       numOS = num;
    const int iHeld = SeqFil(num)->IndexOfS1() == 0   ? 1
                      : SeqFil(num)->IndexOfS2() == 0 ? 2
                                                      : 0;
    if (iHeld != 0 && num1 >= 1 && num1 <= SeqFil.Length())
    {
      numOS = num1;
    }
    const occ::handle<ChFiDS_SurfData>& aSD = SeqFil(numOS);
    if (aSD->IndexOfS1() > 0 && aSD->IndexOfS2() > 0
        && aSD->Vertex(isfirst, 1).IsOnArc() != aSD->Vertex(isfirst, 2).IsOnArc())
    {
      const int                      onsArc = aSD->Vertex(isfirst, 1).IsOnArc() ? 1 : 2;
      const int                      onsAir = 3 - onsArc;
      const ChFiDS_FaceInterference& aFi    = aSD->Interference(onsAir);
      const TopoDS_Face aFAir = TopoDS::Face(DStr.Shape(aSD->Index(onsAir)));
      // the face at end on that side: across the edge of Vtx in it, the
      // only one
      TopoDS_Face aFv;
      int         nbFv = 0;
      for (NCollection_List<TopoDS_Shape>::Iterator itE(myVEMap(Vtx)); itE.More(); itE.Next())
      {
        const TopoDS_Edge& anE = TopoDS::Edge(itE.Value());
        if (anE.IsSame(EdgeSpine) || BRep_Tool::Degenerated(anE) || !containE(aFAir, anE)
            || myEFMap(anE).Extent() != 2)
        {
          continue;
        }
        for (NCollection_List<TopoDS_Shape>::Iterator itF(myEFMap(anE)); itF.More(); itF.Next())
        {
          if (!aFAir.IsSame(itF.Value()) && !itF.Value().IsSame(aFv))
          {
            aFv = TopoDS::Face(itF.Value());
            nbFv++;
          }
        }
      }
      // the walk past the end held the arc's side and stopped on the split
      const bool isHeldOk =
        numOS == num
        || (iHeld == onsArc && SeqFil(num)->Vertex(isfirst, onsAir).IsOnArc()
            && SeqFil(num)->Index(onsAir) == aSD->Index(onsAir));
      if (isHeldOk && nbFv == 1 && !aFv.IsNull()
          && hasVertex(aSD->Vertex(isfirst, onsArc).Arc(), Vtx) && aFi.LineIndex() != 0
          && LineOverSplit(DStr.Curve(aFi.LineIndex()).Curve(),
                           aFi.Parameter(isfirst),
                           isfirst,
                           aFAir,
                           aFv,
                           myEFMap,
                           10 * tolapp3d,
                           LinOS,
                           EsplitOS,
                           FnOS,
                           EarcOS,
                           wQOS,
                           wPOS,
                           parQOS,
                           parPOS)
          && (numOS == num || SeqFil(num)->Vertex(isfirst, onsAir).Arc().IsSame(EsplitOS))
          && ArcPastSplit(Vtx, EarcOS, FnOS, myVEMap, myEFMap, 10 * tolapp3d, EdgePast))
      {
        onsOS = onsAir;
        if (numOS != num)
        {
          SeqFil.Remove(num);
          num  = ChFi3d_IndexOfSurfData(Vtx, stripe, sens);
          num1 = isfirst ? num + 1 : num - 1;
        }
      }
    }
  }

  // there is only one face at end if FindFace is true and if the face
  // is not the face with sewing edge
  TopoDS_Face                  face;
  occ::handle<ChFiDS_SurfData> Fd        = SeqFil.ChangeValue(num);
  ChFiDS_CommonPoint&          CV1       = Fd->ChangeVertex(isfirst, 1);
  ChFiDS_CommonPoint&          CV2       = Fd->ChangeVertex(isfirst, 2);
  bool                         onecorner = false;
  if (onsOS != 0)
  {
    // the line's side ends at P on the far piece's edge
    ChFiDS_CommonPoint&      aCP   = Fd->ChangeVertex(isfirst, onsOS);
    ChFiDS_FaceInterference& aFi   = Fd->ChangeInterference(onsOS);
    const TopoDS_Face        aFAir = TopoDS::Face(DStr.Shape(Fd->Index(onsOS)));
    const double             aTol  = aCP.Tolerance();
    aCP.Reset();
    aCP.SetPoint(LinOS->Value(wPOS));
    aCP.SetArc(aTol,
               EarcOS,
               parPOS,
               TransitionOnEdge(aFi,
                                FnOS,
                                EarcOS,
                                isfirst,
                                FnOS.Orientation() != aFAir.Orientation()));
    aFi.SetParameter(wPOS, isfirst);
  }
  if (FindFace(Vtx, CV1, CV2, face))
  {
    if (!couture)
    {
      onecorner = true;
    }
    else if (!face.IsSame(facecouture))
    {
      onecorner = true;
    }
  }
  if (onecorner)
  {
    if (ChFi3d_Builder::MoreSurfdata(Index))
    {
      ChFi3d_Builder::PerformMoreSurfdata(Index);
      return;
    }
  }
  if (!onecorner && (reg1 || reg2) && !couture && state != ChFiDS_OnSame)
  {
    PerformMoreThreeCorner(Index, 1);
    return;
  }
  occ::handle<GeomAdaptor_Surface> HGs = ChFi3d_BoundSurf(DStr, Fd, 1, 2);
  ChFiDS_FaceInterference          Fi1 = Fd->InterferenceOnS1();
  ChFiDS_FaceInterference          Fi2 = Fd->InterferenceOnS2();
  GeomAdaptor_Surface&             Gs  = *HGs;
  occ::handle<BRepAdaptor_Surface> HBs = new BRepAdaptor_Surface();
  BRepAdaptor_Surface&             Bs  = *HBs;
  occ::handle<Geom_Curve>          Cc;
  occ::handle<Geom2d_Curve>        Pc, Ps;
  double                           Ubid, Vbid;
  TopAbs_Orientation               orsurfdata;
  orsurfdata                             = Fd->Orientation();
  int                          IsurfPrev = 0, Isurf = Fd->Surf();
  occ::handle<ChFiDS_SurfData> SDprev;
  if (num1 > 0 && num1 <= SeqFil.Length())
  {
    SDprev    = SeqFil(num1);
    IsurfPrev = SDprev->Surf();
  }
  // calculate the orientation of curves at end

  double             tolpt = 1.e-4;
  double             tolreached;
  TopAbs_Orientation orcourbe, orface, orien;

  stripe->SetIndexPoint(ChFi3d_IndexPointInDS(CV1, DStr), isfirst, 1);
  stripe->SetIndexPoint(ChFi3d_IndexPointInDS(CV2, DStr), isfirst, 2);

  //  gp_Pnt p3d;
  //  gp_Pnt2d p2d;
  double             dist;
  int                Ishape1 = Fd->IndexOfS1();
  TopAbs_Orientation trafil1 = TopAbs_FORWARD;
  if (Ishape1 != 0)
  {
    if (Ishape1 > 0)
    {
      trafil1 = DStr.Shape(Ishape1).Orientation();
    }
#ifdef OCCT_DEBUG
    else
    {
      std::cout << "erreur" << std::endl;
    }
#endif
    trafil1 = TopAbs::Compose(trafil1, Fd->Orientation());

    trafil1 = TopAbs::Compose(TopAbs::Reverse(Fi1.Transition()), trafil1);
  }
#ifdef OCCT_DEBUG
  else
    std::cout << "erreur" << std::endl;
#endif
  // eap, Apr 22 2002, occ 293
  //   Fi1.PCurveOnFace()->D0(Fi1.LastParameter(),p2d);
  //   const occ::handle<Geom_Surface> Stemp =
  //     BRep_Tool::Surface(TopoDS::Face(DStr.Shape(Ishape1)));
  //   Stemp ->D0(p2d.X(),p2d.Y(),p3d);
  //   dist=p3d.Distance(CV1.Point());
  //   if (dist<tolpt) orcourbe=trafil1;
  //   else            orcourbe=TopAbs::Reverse(trafil1);
  if (!isfirst)
  {
    orcourbe = trafil1;
  }
  else
  {
    orcourbe = TopAbs::Reverse(trafil1);
  }

  // eap, Apr 22 2002, occ 293
  // variables to show OnSame situation
  bool isOnSame1, isOnSame2;
  // In OnSame situation, the case of degenerated FaceInterference curve
  // is probable when a corner cuts the ChFi3d earlier built on OnSame edge.
  // In such a case, chamfer face can partially shrink to a line and we need
  // to cut off that shrinked part
  // If <isOnSame1>, FaceInterference with F2 can be degenerated
  bool checkShrink, isShrink, isUShrink;
  isShrink = isUShrink = isOnSame1 = isOnSame2 = false;
  double   checkShrParam = 0., prevSDParam = 0.;
  gp_Pnt2d midP2d;
  int      midIpoint = 0;

  // find Fi1,Fi2 lengths used to extend ChFi surface
  // and by the way define necessity to check shrink
  gp_Pnt2d P2d1 = Fi1.PCurveOnSurf()->Value(Fi1.Parameter(isfirst));
  gp_Pnt2d P2d2 = Fi1.PCurveOnSurf()->Value(Fi1.Parameter(!isfirst));
  gp_Pnt   aP1, aP2;
  HGs->D0(P2d1.X(), P2d1.Y(), aP1);
  HGs->D0(P2d2.X(), P2d2.Y(), aP2);
  double Fi1Length = aP1.Distance(aP2);
  //  double eps = Precision::Confusion();
  checkShrink = (Fi1Length <= Precision::Confusion());

  gp_Pnt2d P2d3 = Fi2.PCurveOnSurf()->Value(Fi2.Parameter(isfirst));
  gp_Pnt2d P2d4 = Fi2.PCurveOnSurf()->Value(Fi2.Parameter(!isfirst));
  HGs->D0(P2d3.X(), P2d3.Y(), aP1);
  HGs->D0(P2d4.X(), P2d4.Y(), aP2);
  double Fi2Length = aP1.Distance(aP2);
  checkShrink      = checkShrink || (Fi2Length <= Precision::Confusion());

  if (checkShrink)
  {
    if (std::abs(P2d2.Y() - P2d4.Y()) <= Precision::PConfusion())
    {
      isUShrink     = false;
      checkShrParam = P2d2.Y();
    }
    else if (std::abs(P2d2.X() - P2d4.X()) <= Precision::PConfusion())
    {
      isUShrink     = true;
      checkShrParam = P2d2.X();
    }
    else
    {
      checkShrink = false;
    }
  }

  /***********************************************************************/
  //  find faces intersecting with the fillet and edges limiting intersections
  //  nbface is the nb of faces intersected, Face[i] contains the faces
  //  to intersect (i=0.. nbface-1). Edge[i] contains edges limiting
  //  the intersections (i=0 ..nbface)
  /**********************************************************************/

  int                                      nb = 1, nbface;
  TopoDS_Edge                              E1, E2, Edge[nn], E, Ei, edgesau;
  TopoDS_Face                              facesau;
  bool                                     oneintersection1 = false;
  bool                                     oneintersection2 = false;
  TopoDS_Face                              Face[nn], F, F3;
  TopoDS_Vertex                            V1, V2, V, Vfin;
  bool                                     findonf1 = false, findonf2 = false;
  NCollection_List<TopoDS_Shape>::Iterator It3;
  F1 = TopoDS::Face(DStr.Shape(Fd->IndexOfS1()));
  F2 = TopoDS::Face(DStr.Shape(Fd->IndexOfS2()));
  F3 = F1;
  if (couture || bordlibre)
  {
    nface = nface + 1;
  }
  if (nface == 3)
  {
    nbface = 2;
  }
  else
  {
    nbface = nface - 2;
  }
  if (!CV1.IsOnArc() || !CV2.IsOnArc())
  {
    PerformMoreThreeCorner(Index, 1);
    return;
  }

  Edge[0]      = CV1.Arc();
  Edge[nbface] = CV2.Arc();
  tabedg.SetValue(0, Edge[0]);
  tabedg.SetValue(nbface, Edge[nbface]);
  // processing of a fillet arriving on a vertex
  // edge contained in CV.Arc is not inevitably good
  // the edge concerned by the intersection is found

  double dist1, dist2;
  if (CV1.IsVertex())
  {
    trouve               = false;
    /*TopoDS_Vertex */ V = CV1.Vertex();
    for (It3.Initialize(myVEMap(V)); It3.More() && !trouve; It3.Next())
    {
      E = TopoDS::Edge(It3.Value());
      if (!E.IsSame(Edge[0]) && (containE(F1, E)))
      {
        trouve = true;
      }
    }
    TopoDS_Vertex Vt, V3, V4;
    V1 = TopExp::FirstVertex(Edge[0]);
    V2 = TopExp::LastVertex(Edge[0]);
    if (V.IsSame(V1))
    {
      Vt = V2;
    }
    else
    {
      Vt = V1;
    }
    dist1 = (BRep_Tool::Pnt(Vt)).Distance(BRep_Tool::Pnt(Vtx));
    V3    = TopExp::FirstVertex(E);
    V4    = TopExp::LastVertex(E);
    if (V.IsSame(V3))
    {
      Vt = V4;
    }
    else
    {
      Vt = V3;
    }
    dist2 = (BRep_Tool::Pnt(Vt)).Distance(BRep_Tool::Pnt(Vtx));
    if (dist2 < dist1)
    {
      Edge[0] = E;
      TopAbs_Orientation ori;
      if (V2.IsSame(V3) || V1.IsSame(V4))
      {
        ori = CV1.TransitionOnArc();
      }
      else
      {
        ori = TopAbs::Reverse(CV1.TransitionOnArc());
      }
      double par = BRep_Tool::Parameter(V, Edge[0]);
      double tol = CV1.Tolerance();
      CV1.SetArc(tol, Edge[0], par, ori);
    }
  }

  if (CV2.IsVertex())
  {
    trouve              = false;
    /*TopoDS_Vertex*/ V = CV2.Vertex();
    for (It3.Initialize(myVEMap(V)); It3.More() && !trouve; It3.Next())
    {
      E = TopoDS::Edge(It3.Value());
      if (!E.IsSame(Edge[2]) && (containE(F2, E)))
      {
        trouve = true;
      }
    }
    TopoDS_Vertex Vt, V3, V4;
    V1 = TopExp::FirstVertex(Edge[2]);
    V2 = TopExp::LastVertex(Edge[2]);
    if (V.IsSame(V1))
    {
      Vt = V2;
    }
    else
    {
      Vt = V1;
    }
    dist1 = (BRep_Tool::Pnt(Vt)).Distance(BRep_Tool::Pnt(Vtx));
    V3    = TopExp::FirstVertex(E);
    V4    = TopExp::LastVertex(E);
    if (V.IsSame(V3))
    {
      Vt = V4;
    }
    else
    {
      Vt = V3;
    }
    dist2 = (BRep_Tool::Pnt(Vt)).Distance(BRep_Tool::Pnt(Vtx));
    if (dist2 < dist1)
    {
      Edge[2] = E;
      TopAbs_Orientation ori;
      if (V2.IsSame(V3) || V1.IsSame(V4))
      {
        ori = CV2.TransitionOnArc();
      }
      else
      {
        ori = TopAbs::Reverse(CV2.TransitionOnArc());
      }
      double par = BRep_Tool::Parameter(V, Edge[2]);
      double tol = CV2.Tolerance();
      CV2.SetArc(tol, Edge[2], par, ori);
    }
  }
  if (!onecorner)
  {
    // If there is a regular edge, the faces adjacent to it
    // are not in Fd->IndexOfS1 or Fd->IndexOfS2

    //     TopoDS_Face Find1 ,Find2;
    //     if (isfirst)
    //       edge=stripe->Spine()->Edges(1);
    //     else  edge=stripe->Spine()->Edges(stripe->Spine()->NbEdges());
    //     It3.Initialize(myEFMap(edge));
    //     Find1=TopoDS::Face(It3.Value());
    //     trouve=false;
    //     for (It3.Initialize(myEFMap(edge));It3.More()&&!trouve;It3.Next()) {
    //       F=TopoDS::Face (It3.Value());
    //       if (!F.IsSame(Find1)) {
    // 	Find2=F;trouve=true;
    //       }
    //     }

    // if nface =3 there is a top with 3 edges and a fillet
    // and their common points are on different faces
    // otherwise there is a case when a top has more than 3 edges

    if (nface == 3)
    {
      if (CV1.IsVertex())
      {
        findonf1 = true;
      }
      if (CV2.IsVertex())
      {
        findonf2 = true;
      }
      if (!findonf1)
      {
        NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> MapV;
        TopExp::MapShapes(Edge[0], TopAbs_VERTEX, MapV);
        if (MapV.Extent() == 2)
        {
          if (!MapV(1).IsSame(Vtx) && !MapV(2).IsSame(Vtx))
          {
            findonf1 = true;
          }
        }
      }
      if (!findonf2)
      {
        NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> MapV;
        TopExp::MapShapes(Edge[2], TopAbs_VERTEX, MapV);
        if (MapV.Extent() == 2)
        {
          if (!MapV(1).IsSame(Vtx) && !MapV(2).IsSame(Vtx))
          {
            findonf2 = true;
          }
        }
      }

      // detect and process OnSame situatuation
      if (state == ChFiDS_OnSame)
      {
        TopoDS_Edge threeE[3];
        ChFi3d_cherche_element(Vtx, EdgeSpine, F1, threeE[0], V2);
        ChFi3d_cherche_element(Vtx, EdgeSpine, F2, threeE[1], V2);
        threeE[2] = EdgeSpine;
        if (ChFi3d_EdgeState(threeE, myEFMap) == ChFiDS_OnSame)
        {
          isOnSame1 = true;
          nb        = 1;
          Edge[0]   = threeE[0];
          ChFi3d_cherche_face1(myEFMap(Edge[0]), F1, Face[0]);
          if (findonf2)
          {
            findonf1 = true; // not to look for Face[0] again
          }
          else
          {
            Edge[1] = CV2.Arc();
          }
        }
        else
        {
          isOnSame2 = true;
        }
      }

      // findonf1 findonf2 show if F1 and/or F2 are adjacent
      // to many faces at end
      // the faces at end and intersected edges are found

      if (findonf1 && !isOnSame1)
      {
        if (CV1.TransitionOnArc() == TopAbs_FORWARD)
        {
          V1 = TopExp::FirstVertex(CV1.Arc());
        }
        else
        {
          V1 = TopExp::LastVertex(CV1.Arc());
        }
        ChFi3d_cherche_face1(myEFMap(CV1.Arc()), F1, Face[0]);
        nb = 1;
        Ei = Edge[0];
        while (!V1.IsSame(Vtx))
        {
          ChFi3d_cherche_element(V1, Ei, F1, E, V2);
          V1 = V2;
          Ei = E;
          ChFi3d_cherche_face1(myEFMap(E), F1, Face[nb]);
          cherche_edge1(Face[nb - 1], Face[nb], Edge[nb]);
          nb++;
          if (nb >= nn)
          {
            throw Standard_Failure("IntersectionAtEnd : the max number of faces reached");
          }
        }
        if (!findonf2)
        {
          Edge[nb] = CV2.Arc();
        }
      }
      if (findonf2 && !isOnSame2)
      {
        if (!findonf1)
        {
          nb = 1;
        }
        V1 = Vtx;
        if (CV2.TransitionOnArc() == TopAbs_FORWARD)
        {
          Vfin = TopExp::LastVertex(CV2.Arc());
        }
        else
        {
          Vfin = TopExp::FirstVertex(CV2.Arc());
        }
        if (!findonf1)
        {
          ChFi3d_cherche_face1(myEFMap(CV1.Arc()), F1, Face[nb - 1]);
        }
        ChFi3d_cherche_element(V1, EdgeSpine, F2, E, V2);
        Ei = E;
        V1 = V2;
        while (!V1.IsSame(Vfin))
        {
          ChFi3d_cherche_element(V1, Ei, F2, E, V2);
          Ei = E;
          V1 = V2;
          ChFi3d_cherche_face1(myEFMap(E), F2, Face[nb]);
          cherche_edge1(Face[nb - 1], Face[nb], Edge[nb]);
          nb++;
          if (nb >= nn)
          {
            throw Standard_Failure("IntersectionAtEnd : the max number of faces reached");
          }
        }
        Edge[nb] = CV2.Arc();
      }
      if (isOnSame2)
      {
        cherche_edge1(Face[nb - 1], F2, Edge[nb]);
        Face[nb] = F2;
      }

      nbface = nb;
    }

    else
    {

      //  this is the case when a top has more than three edges
      //  the faces and edges concerned are found
      bool /*trouve,*/ possible1, possible2;
      trouve = possible1 = possible2 = false;
      TopExp_Explorer ex;
      nb = 0;
      for (ex.Init(CV1.Arc(), TopAbs_VERTEX); ex.More(); ex.Next())
      {
        if (Vtx.IsSame(ex.Current()))
        {
          possible1 = true;
        }
      }
      for (ex.Init(CV2.Arc(), TopAbs_VERTEX); ex.More(); ex.Next())
      {
        if (Vtx.IsSame(ex.Current()))
        {
          possible2 = true;
        }
      }
      // A side ending on an edge away from Vtx past a split of its face in
      // coplanar pieces: the corner is that of the faces around Vtx, as on
      // the wall in one face, the edge of Vtx the arc carries on left out.
      if (possible1 != possible2 && nbarete == 4 && !couture && !bordlibre)
      {
        const int iaway = possible1 ? 2 : 1;
        if (ArcPastSplit(Vtx,
                         Fd->Vertex(isfirst, iaway).Arc(),
                         iaway == onsOS ? FnOS : (iaway == 1 ? F1 : F2),
                         myVEMap,
                         myEFMap,
                         10 * tolapp3d,
                         EdgePast))
        {
          possible1 = possible2 = true;
          tabedg.SetValue(nn, EdgePast);
        }
      }
      if ((possible1 && possible2) || (!possible1 && !possible2) || (nbarete > 4))
      {
        while (!trouve)
        {
          nb++;
          if (nb >= nn)
          {
            throw Standard_Failure("IntersectionAtEnd : the max number of faces reached");
          }
          if (nb != 1)
          {
            F3 = Face[nb - 2];
          }
          Face[nb - 1] = F3;
          if (CV1.Arc().IsSame(edgelibre1))
          {
            cherche_face(myVFMap(Vtx), edgelibre2, F1, F2, F3, Face[nb - 1]);
          }
          else if (CV1.Arc().IsSame(edgelibre2))
          {
            cherche_face(myVFMap(Vtx), edgelibre1, F1, F2, F3, Face[nb - 1]);
          }
          else
          {
            cherche_face(myVFMap(Vtx), Edge[nb - 1], F1, F2, F3, Face[nb - 1]);
          }
          ChFi3d_cherche_edge(Vtx, tabedg, Face[nb - 1], Edge[nb], V);
          tabedg.SetValue(nb, Edge[nb]);
          if (Edge[nb].IsSame(CV2.Arc()))
          {
            trouve = true;
          }
        }
        nbface = nb;
      }
      else
      {
        const int                ip  = possible1 ? 1 : 2, io = 3 - ip;
        ChFiDS_FaceInterference& Fio = Fd->ChangeInterference(io);
        const double             w   = Fd->Interference(ip).Parameter(isfirst);
        if (!couture && !bordlibre && w >= Fio.FirstParameter() && w <= Fio.LastParameter())
        {
          // The point of the line past Vtx in the section through the
          // other side's end
          const gp_Pnt2d uv = Fio.PCurveOnSurf()->Value(w);
          const gp_Pnt   Po = DStr.Surface(Fd->Surf()).Surface()->Value(uv.X(), uv.Y());
          if (SectionCrossesEndFace(Vtx,
                                    Fd->Vertex(isfirst, ip).Arc(),
                                    ip == 1 ? F1 : F2,
                                    Fd->Vertex(isfirst, ip).Point(),
                                    Po,
                                    myVEMap,
                                    myEFMap))
          {
            // The line is cut there, as an edge splitting its face there
            // would cut it, and the corner is filled as that of the split
            // face.
            Fio.SetParameter(w, isfirst);
            ChFiDS_CommonPoint& CVo = Fd->ChangeVertex(isfirst, io);
            CVo.Reset();
            CVo.SetPoint(Po);
            CVo.SetTolerance(tolapp3d);
            stripe->SetIndexPoint(ChFi3d_IndexPointInDS(CVo, DStr), isfirst, io);
            PerformMoreThreeCorner(Index, 1);
            return;
          }
        }
        IntersectMoreCorner(Index);
        return;
      }
      if (nbarete == 4)
      {
        // if two consecutive edges are G1 there is only one face of intersection
        double        ang1 = 0.0;
        TopoDS_Vertex Vcom;
        trouve = false;
        ChFi3d_cherche_vertex(Edge[0], Edge[1], Vcom, trouve);
        if (Vcom.IsSame(Vtx))
        {
          ang1 = ChFi3d_AngleEdge(Vtx, Edge[0], Edge[1]);
        }
        if (std::abs(ang1 - M_PI) < 0.01)
        {
          oneintersection1 = true;
          facesau          = Face[0];
          edgesau          = Edge[1];
          Face[0]          = Face[1];
          Edge[1]          = Edge[2];
          nbface           = 1;
        }

        if (!oneintersection1)
        {
          trouve = false;
          ChFi3d_cherche_vertex(Edge[1], Edge[2], Vcom, trouve);
          if (Vcom.IsSame(Vtx))
          {
            ang1 = ChFi3d_AngleEdge(Vtx, Edge[1], Edge[2]);
          }
          if (std::abs(ang1 - M_PI) < 0.01)
          {
            oneintersection2 = true;
            facesau          = Face[1];
            edgesau          = Edge[1];
            Edge[1]          = Edge[2];
            nbface           = 1;
          }
        }
      }
      else if (nbarete == 5)
      {
        // pro15368
        //  Modified by Sergey KHROMOV - Fri Dec 21 18:07:43 2001 End
        bool isTangent0 = ChFi3d::IsTangentFaces(Edge[0], F1, Face[0]);
        bool isTangent1 = ChFi3d::IsTangentFaces(Edge[1], Face[0], Face[1]);
        bool isTangent2 = ChFi3d::IsTangentFaces(Edge[2], Face[1], Face[2]);
        if ((isTangent0 || isTangent2) && isTangent1)
        {
          //         GeomAbs_Shape cont0,cont1,cont2;
          //         cont0=BRep_Tool::Continuity(Edge[0],F1,Face[0]);
          //         cont1=BRep_Tool::Continuity(Edge[1],Face[0],Face[1]);
          //         cont2=BRep_Tool::Continuity(Edge[2],Face[1],Face[2]);
          //         if ((cont0!=GeomAbs_C0 || cont2!=GeomAbs_C0) && cont1!=GeomAbs_C0) {
          //  Modified by Sergey KHROMOV - Fri Dec 21 18:07:49 2001 Begin
          facesau          = Face[0];
          edgesau          = Edge[0];
          nbface           = 1;
          Edge[1]          = Edge[3];
          Face[0]          = Face[2];
          oneintersection1 = true;
        }
      }
    }
  }
  else
  {
    nbface  = 1;
    Face[0] = face;
    Edge[1] = Edge[2];
  }

  NCollection_Array1<double> Pardeb(1, 4), Parfin(1, 4);
  gp_Pnt2d                   pfil1, pfac1, pfil2, pfac2, pint, pfildeb;
  occ::handle<Geom2d_Curve>  Hc1, Hc2;
  IntCurveSurface_HInter     inters;
  // all of them zero: an extension is made once, ChFi3d_ExtendSurface
  // leaving a face it is given as made alone (the last prolface[nn] is for Fd)
  int                        proledge[nn] = {}, prolface[nn + 1] = {};
  int                        shrink[nn]   = {};
  TopoDS_Face                faceprol[nn];
  int                        indcurve[nn], indpoint2 = 0, indpoint1 = 0;
  occ::handle<TopOpeBRepDS_CurvePointInterference>   Interfp1, Interfp2, Interfedge[nn];
  occ::handle<TopOpeBRepDS_SurfaceCurveInterference> Interfc, InterfPC[nn], InterfPS[nn];
  double                                             u2, v2, p1, p2, paredge1;
  double                                             paredge2 = 0., tolex = 1.e-4;
  bool                                               extend = false;
  occ::handle<Geom_Surface>                          Sfacemoins1, Sface;
  /***************************************************************************/
  // calculate intersection of the fillet and each face
  // and storage in the DS
  /***************************************************************************/
  for (nb = 1; nb <= nbface; nb++)
  {
    prolface[nb - 1] = 0;
    proledge[nb - 1] = 0;
    shrink[nb - 1]   = 0;
  }
  proledge[nbface] = 0;
  prolface[nn]     = 0;
  if (oneintersection1 || oneintersection2)
  {
    faceprol[1] = facesau;
  }
  if (!isOnSame1 && !isOnSame2)
  {
    checkShrink = false;
  }
  // in OnSame situation we need intersect Fd with Edge[0] or Edge[nbface] as well
  if (isOnSame1)
  {
    nb = 0;
  }
  else
  {
    nb = 1;
  }
  bool intersOnSameFailed = false;

  for (; nb <= nbface; nb++)
  {
    extend = false;
    // the fillet meets the line of an edge between two faces beyond the
    // edge's end at Vtx: the faces grow by the piece of that line from Vtx
    bool pastVtx = false;
    E2           = Edge[nb];
    if (!nb)
    {
      F = F1;
    }
    else
    {
      F = Face[nb - 1];
      if (!prolface[nb - 1])
      {
        faceprol[nb - 1] = F;
      }
    }

    if (F.IsNull())
    {
      throw Standard_NullObject("IntersectionAtEnd : Trying to intersect with NULL face");
    }

    Sfacemoins1 = BRep_Tool::Surface(F);
    occ::handle<Geom_Curve>   cint;
    occ::handle<Geom2d_Curve> C2dint1, C2dint2, cface, cfacemoins1;

    ///////////////////////////////////////////////////////
    // determine intersections of edges and the fillet
    // to find limitations of intersections face - fillet
    ///////////////////////////////////////////////////////

    if (nb == 1)
    {
      Hc1 = BRep_Tool::CurveOnSurface(Edge[0], Face[0], Ubid, Ubid);
      if (isOnSame1)
      {
        // update interference param on Fi1 and point of CV1
        if (prolface[0])
        {
          Bs.Initialize(faceprol[0], false);
        }
        else
        {
          Bs.Initialize(Face[0], false);
        }
        const occ::handle<Geom_Curve>& c3df = DStr.Curve(Fi1.LineIndex()).Curve();
        double                         Ufi  = Fi2.Parameter(isfirst);
        ChFiDS_FaceInterference&       Fi   = Fd->ChangeInterferenceOnS1();
        if (!IntersUpdateOnSame(HGs,
                                HBs,
                                c3df,
                                F1,
                                Face[0],
                                Edge[0],
                                Vtx,
                                isfirst,
                                10 * tolapp3d, // in
                                Fi,
                                CV1,
                                pfac1,
                                Ufi))
        { // out
          throw Standard_Failure("IntersectionAtEnd: pb intersection Face - Fi");
        }
        Fi1 = Fi;
        if (intersOnSameFailed)
        { // probable at fillet building
          // look for paredge2
          Geom2dAPI_ProjectPointOnCurve proj;
          if (C2dint2.IsNull())
          {
            proj.Init(pfac1, Hc1);
          }
          else
          {
            proj.Init(pfac1, C2dint2);
          }
          paredge2 = proj.LowerDistanceParameter();
        }
        // update stripe point
        TopOpeBRepDS_Point tpoint(CV1.Point(), tolapp3d);
        indpoint1 = DStr.AddPoint(tpoint);
        stripe->SetIndexPoint(indpoint1, isfirst, 1);
        // reset arc of CV1
        TopoDS_Vertex vert1, vert2;
        TopExp::Vertices(Edge[0], vert1, vert2);
        TopAbs_Orientation arcOri = Vtx.IsSame(vert1) ? TopAbs_FORWARD : TopAbs_REVERSED;
        CV1.SetArc(tolapp3d, Edge[0], paredge2, arcOri);
      }
      else
      {
        if (Hc1.IsNull())
        {
          // curve 2d not found. Sfacemoins1 is extended and projection is done there
          // CV1.Point ()
          ChFi3d_ExtendSurface(Sfacemoins1, prolface[0]);
          if (prolface[0])
          {
            extend = true;
            BRep_Builder BRE;
            double       tol = BRep_Tool::Tolerance(F);
            BRE.MakeFace(faceprol[0], Sfacemoins1, F.Location(), tol);
            if (!isOnSame1)
            {
              GeomAdaptor_Surface Asurf;
              Asurf.Load(Sfacemoins1);
              Extrema_ExtPS ext(CV1.Point(), Asurf, tol, tol, Extrema_ExtFlag_MIN);
              double        uc1, vc1;
              if (ext.IsDone())
              {
                ext.Point(1).Parameter(uc1, vc1);
                pfac1.SetX(uc1);
                pfac1.SetY(vc1);
              }
            }
          }
        }
        else
        {
          pfac1 = Hc1->Value(CV1.ParameterOnArc());
        }
      }
      paredge1 = CV1.ParameterOnArc();
      if (Fi1.LineIndex() != 0)
      {
        pfil1 = Fi1.PCurveOnSurf()->Value(Fi1.Parameter(isfirst));
      }
      else
      {
        pfil1 = Fi1.PCurveOnSurf()->Value(Fi1.Parameter(!isfirst));
      }
      pfildeb = pfil1;
    }
    else
    {
      pfil1    = pfil2;
      paredge1 = paredge2;
      pfac1    = pint;
    }

    if (nb != nbface || isOnSame2)
    {
      int nbp;

      occ::handle<Geom_Curve> C;
      C                                    = BRep_Tool::Curve(E2, Ubid, Vbid);
      occ::handle<Geom_TrimmedCurve> Ctrim = new Geom_TrimmedCurve(C, Ubid, Vbid);
      double                         Utrim, Vtrim;
      Utrim = Ctrim->BasisCurve()->FirstParameter();
      Vtrim = Ctrim->BasisCurve()->LastParameter();
      if (Ctrim->IsPeriodic())
      {
        if (Ubid > Ctrim->Period())
        {
          Ubid = (Utrim + Vtrim) / 2;
          Vbid = Vtrim;
        }
        else
        {
          Ubid = Utrim;
          Vbid = (Utrim + Vtrim) / 2;
        }
      }
      else
      {
        Ubid = Utrim;
        Vbid = Vtrim;
      }
      occ::handle<GeomAdaptor_Curve> HC  = new GeomAdaptor_Curve(C, Ubid, Vbid);
      GeomAdaptor_Curve&             Cad = *HC;
      inters.Perform(HC, HGs);
      if (!prolface[nn] && (!inters.IsDone() || (inters.NbPoints() == 0)))
      {
        // extend surface of conge
        occ::handle<Geom_BoundedSurface> S1 =
          occ::down_cast<Geom_BoundedSurface>(DStr.Surface(Fd->Surf()).Surface());
        if (!S1.IsNull())
        {
          double length = 0.5 * std::max(Fi1Length, Fi2Length);
          GeomLib::ExtendSurfByLength(S1, length, 1, false, !isfirst);
          prolface[nn] = 1;
          if (!stripe->IsInDS(!isfirst))
          {
            Gs.Load(S1);
            inters.Perform(HC, HGs);
            if (inters.IsDone() && inters.NbPoints() != 0)
            {
              Fd->ChangeSurf(
                DStr.AddSurface(TopOpeBRepDS_Surface(S1, DStr.ChangeSurface(Isurf).Tolerance())));
              // update history
              if (myEVIMap.IsBound(EdgeSpine))
              {
                NCollection_List<int>::Iterator itl(myEVIMap.ChangeFind(EdgeSpine));
                for (; itl.More(); itl.Next())
                {
                  if (itl.Value() == Isurf)
                  {
                    myEVIMap.ChangeFind(EdgeSpine).Remove(itl);
                    break;
                  }
                }
                myEVIMap.ChangeFind(EdgeSpine).Append(Fd->Surf());
              }
              else
              {
                NCollection_List<int> IndexList;
                IndexList.Append(Fd->Surf());
                myEVIMap.Bind(EdgeSpine, IndexList);
              }
              ////////////////
              Isurf = Fd->Surf();
            }
          }
        }
      }
      if (!inters.IsDone() || (inters.NbPoints() == 0))
      {
        occ::handle<Geom_BSplineCurve> cd  = occ::down_cast<Geom_BSplineCurve>(C);
        occ::handle<Geom_BezierCurve>  cd1 = occ::down_cast<Geom_BezierCurve>(C);
        if (!cd.IsNull() || !cd1.IsNull())
        {
          BRep_Builder BRE;
          Sface = BRep_Tool::Surface(Face[nb]);
          ChFi3d_ExtendSurface(Sface, prolface[nb]);
          double tol = BRep_Tool::Tolerance(F);
          BRE.MakeFace(faceprol[nb], Sface, Face[nb].Location(), tol);
          if (nb && !prolface[nb - 1])
          {
            ChFi3d_ExtendSurface(Sfacemoins1, prolface[nb - 1]);
            if (prolface[nb - 1])
            {
              tol = BRep_Tool::Tolerance(F);
              BRE.MakeFace(faceprol[nb - 1], Sfacemoins1, F.Location(), tol);
            }
          }
          else
          {
            int prol = 0;
            ChFi3d_ExtendSurface(Sfacemoins1, prol);
          }
          GeomInt_IntSS InterSS(Sfacemoins1, Sface, 1.e-7, true, true, true);
          if (InterSS.IsDone())
          {
            trouve = false;
            for (int i = 1; i <= InterSS.NbLines() && !trouve; i++)
            {
              extend  = true;
              cint    = InterSS.Line(i);
              C2dint1 = InterSS.LineOnS1(i);
              C2dint2 = InterSS.LineOnS2(i);
              Cad.Load(cint);
              inters.Perform(HC, HGs);
              trouve = inters.IsDone() && inters.NbPoints() != 0;
              // eap occ293, eval tolex on finally trimmed curves
              //               occ::handle<GeomAdaptor_Surface> H1=new
              //               GeomAdaptor_Surface(Sfacemoins1); occ::handle<GeomAdaptor_Surface>
              //               H2=new GeomAdaptor_Surface(Sface);
              //              tolex=ChFi3d_EvalTolReached(H1,C2dint1,H2,C2dint2,cint);
              tolex = InterSS.TolReached3d();
            }
          }
        }
      }
      if (inters.IsDone())
      {
        nbp = inters.NbPoints();
        if (nbp == 0)
        {
          if (nb == 0 || nb == nbface)
          {
            intersOnSameFailed = true;
          }
          else
          {
            PerformMoreThreeCorner(Index, 1);
            return;
          }
        }
        else
        {
          gp_Pnt P       = BRep_Tool::Pnt(Vtx);
          double distmin = P.Distance(inters.Point(1).Pnt());
          nbp            = 1;
          for (int i = 2; i <= inters.NbPoints(); i++)
          {
            dist = P.Distance(inters.Point(i).Pnt());
            if (dist < distmin)
            {
              distmin = dist;
              nbp     = i;
            }
          }
          gp_Pnt2d pt2d(inters.Point(nbp).U(), inters.Point(nbp).V());
          pfil2    = pt2d;
          paredge2 = inters.Point(nbp).W();
          if (!extend)
          {
            cfacemoins1 = BRep_Tool::CurveOnSurface(E2, F, u2, v2);
            if (cfacemoins1.IsNull())
            {
              throw Standard_ConstructionError("Failed to get p-curve of edge");
            }
            cface = BRep_Tool::CurveOnSurface(E2, Face[nb], u2, v2);
            if (cface.IsNull())
            {
              throw Standard_ConstructionError("Failed to get p-curve of edge");
            }
            cfacemoins1->D0(paredge2, pfac2);
            cface->D0(paredge2, pint);
            // A point past the edge's end, given to the edge, would lengthen
            // it -- and a second fillet cutting the same edge at its other end
            // (a seam of coplanar faces under both) would then leave both
            // versions of the edge. Across a seam of tangent faces, the piece
            // beyond Vtx is a curve of its own.
            if (nb != nbface && ChFi3d::IsTangentFaces(E2, F, Face[nb]))
            {
              double aF, aL;
              BRep_Tool::Range(E2, aF, aL);
              const gp_Pnt  aP    = C->Value(paredge2);
              const double  aTolE = BRep_Tool::Tolerance(E2);
              TopoDS_Vertex aVF, aVL;
              TopExp::Vertices(E2, aVF, aVL);
              if ((paredge2 < aF && aVF.IsSame(Vtx) && aP.Distance(BRep_Tool::Pnt(aVF)) > aTolE)
                  || (paredge2 > aL && aVL.IsSame(Vtx)
                      && aP.Distance(BRep_Tool::Pnt(aVL)) > aTolE))
              {
                pastVtx = true;
                cint    = C;
                C2dint1 = cfacemoins1;
                C2dint2 = cface;
              }
            }
          }
          else if (C2dint1.IsNull() || C2dint2.IsNull())
          {
            throw Standard_ConstructionError("Failed to get p-curve of edge");
          }
          else
          {
            C2dint1->D0(paredge2, pfac2);
            C2dint2->D0(paredge2, pint);
          }
        }
      }
      else
      {
        throw Standard_Failure("IntersectionAtEnd: pb intersection Face cb");
      }
    }
    else
    {
      Hc2 = BRep_Tool::CurveOnSurface(E2, Face[nbface - 1], Ubid, Ubid);
      if (Hc2.IsNull())
      {
        // curve 2d is not found,  Sfacemoins1 is extended CV2.Point() is projected there

        ChFi3d_ExtendSurface(Sfacemoins1, prolface[0]);
        if (prolface[0])
        {
          BRep_Builder BRE;
          extend     = true;
          double tol = BRep_Tool::Tolerance(F);
          BRE.MakeFace(faceprol[nb - 1], Sfacemoins1, F.Location(), tol);
          GeomAdaptor_Surface Asurf;
          Asurf.Load(Sfacemoins1);
          Extrema_ExtPS ext(CV2.Point(), Asurf, tol, tol, Extrema_ExtFlag_MIN);
          double        uc2, vc2;
          if (ext.IsDone())
          {
            ext.Point(1).Parameter(uc2, vc2);
            pfac2.SetX(uc2);
            pfac2.SetY(vc2);
          }
        }
      }
      else
      {
        pfac2 = Hc2->Value(CV2.ParameterOnArc());
      }
      paredge2 = CV2.ParameterOnArc();
      if (Fi2.LineIndex() != 0)
      {
        pfil2 = Fi2.PCurveOnSurf()->Value(Fi2.Parameter(isfirst));
      }
      else
      {
        pfil2 = Fi2.PCurveOnSurf()->Value(Fi2.Parameter(!isfirst));
      }
    }
    if (!nb)
    {
      continue; // found paredge1 on Edge[0] in OnSame situation on F1
    }

    if (nb == nbface && isOnSame2)
    {
      // update interference param on Fi2 and point of CV2
      if (prolface[nb - 1])
      {
        Bs.Initialize(faceprol[nb - 1]);
      }
      else
      {
        Bs.Initialize(Face[nb - 1]);
      }
      const occ::handle<Geom_Curve>& c3df = DStr.Curve(Fi2.LineIndex()).Curve();
      double                         Ufi  = Fi1.Parameter(isfirst);
      ChFiDS_FaceInterference&       Fi   = Fd->ChangeInterferenceOnS2();
      if (!IntersUpdateOnSame(HGs,
                              HBs,
                              c3df,
                              F2,
                              F,
                              Edge[nb],
                              Vtx,
                              isfirst,
                              10 * tolapp3d, // in
                              Fi,
                              CV2,
                              pfac2,
                              Ufi))
      { // out
        throw Standard_Failure("IntersectionAtEnd: pb intersection Face - Fi");
      }
      Fi2 = Fi;
      if (intersOnSameFailed)
      { // probable at fillet building
        // look for paredge2
        Geom2dAPI_ProjectPointOnCurve proj;
        if (extend)
        {
          proj.Init(pfac2, C2dint2);
        }
        else
        {
          proj.Init(pfac2, BRep_Tool::CurveOnSurface(E2, Face[nbface - 1], Ubid, Ubid));
        }
        paredge2 = proj.LowerDistanceParameter();
      }
      // update stripe point
      TopOpeBRepDS_Point tpoint(CV2.Point(), tolapp3d);
      indpoint2 = DStr.AddPoint(tpoint);
      stripe->SetIndexPoint(indpoint2, isfirst, 2);
      // reset arc of CV2
      TopoDS_Vertex vert1, vert2;
      TopExp::Vertices(Edge[nbface], vert1, vert2);
      TopAbs_Orientation arcOri = Vtx.IsSame(vert1) ? TopAbs_FORWARD : TopAbs_REVERSED;
      CV2.SetArc(tolapp3d, Edge[nbface], paredge2, arcOri);
    }

    if (prolface[nb - 1])
    {
      Bs.Initialize(faceprol[nb - 1]);
    }
    else
    {
      Bs.Initialize(Face[nb - 1]);
    }

    // offset of parameters if they are not in the same period

    // commented by eap 30 May 2002 occ354
    // the following code may cause trimming a wrong part of periodic surface

    //     double  deb,xx1,xx2;
    //     bool  moins2pi,moins2pi1,moins2pi2;
    //     if (DStr.Surface(Fd->Surf()).Surface()->IsUPeriodic()) {
    //       deb=pfildeb.X();
    //       xx1=pfil1.X();
    //       xx2=pfil2.X();
    //       moins2pi=std::abs(deb)< std::abs(std::abs(deb)-2*M_PI);
    //       moins2pi1=std::abs(xx1)< std::abs(std::abs(xx1)-2*M_PI);
    //       moins2pi2=std::abs(xx2)< std::abs(std::abs(xx2)-2*M_PI);
    //       if (moins2pi1!=moins2pi2) {
    //         if  (moins2pi) {
    //           if (!moins2pi1) xx1=xx1-2*M_PI;
    //           if (!moins2pi2) xx2=xx2-2*M_PI;
    //         }
    //         else {
    //           if (moins2pi1) xx1=xx1+2*M_PI;
    //           if (moins2pi2) xx2=xx2+2*M_PI;
    //         }
    //       }
    //       pfil1.SetX(xx1);
    //       pfil2.SetX(xx2);
    //     }
    //     if (couture || Sfacemoins1->IsUPeriodic()) {

    //       double ufmin,ufmax,vfmin,vfmax;
    //       BRepTools::UVBounds(Face[nb-1],ufmin,ufmax,vfmin,vfmax);
    //       deb=ufmin;
    //       xx1=pfac1.X();
    //       xx2=pfac2.X();
    //       moins2pi=std::abs(deb)< std::abs(std::abs(deb)-2*M_PI);
    //       moins2pi1=std::abs(xx1)< std::abs(std::abs(xx1)-2*M_PI);
    //       moins2pi2=std::abs(xx2)< std::abs(std::abs(xx2)-2*M_PI);
    //       if (moins2pi1!=moins2pi2) {
    //         if  (moins2pi) {
    //           if (!moins2pi1) xx1=xx1-2*M_PI;
    //           if (!moins2pi2) xx2=xx2-2*M_PI;
    //         }
    //         else {
    //           if (moins2pi1) xx1=xx1+2*M_PI;
    //           if (moins2pi2) xx2=xx2+2*M_PI;
    //         }
    //       }
    //       pfac1.SetX(xx1);
    //       pfac2.SetX(xx2);
    //     }

    Pardeb(1) = pfil1.X();
    Pardeb(2) = pfil1.Y();
    Pardeb(3) = pfac1.X();
    Pardeb(4) = pfac1.Y();
    Parfin(1) = pfil2.X();
    Parfin(2) = pfil2.Y();
    Parfin(3) = pfac2.X();
    Parfin(4) = pfac2.Y();

    double uu1, uu2, vv1, vv2;
    ChFi3d_Boite(pfac1, pfac2, uu1, uu2, vv1, vv2);
    ChFi3d_BoundFac(Bs, uu1, uu2, vv1, vv2);

    //////////////////////////////////////////////////////////////////////
    // calculate intersections face - fillet
    //////////////////////////////////////////////////////////////////////

    if (!ChFi3d_ComputeCurves(HGs,
                              HBs,
                              Pardeb,
                              Parfin,
                              Cc,
                              Ps,
                              Pc,
                              tolapp3d,
                              tol2d,
                              tolreached,
                              nbface == 1))
    {
      PerformMoreThreeCorner(Index, 1);
      return;
    }
    // storage of information in the data structure

    // evaluate tolerances
    p1 = Cc->FirstParameter();
    p2 = Cc->LastParameter();
    double   to1, to2;
    gp_Pnt2d p2d1, p2d2;
    gp_Pnt   P1, P2, P3, P4, P5, P6, P7, P8;
    HGs->D0(Pardeb(1), Pardeb(2), P1);
    HGs->D0(Parfin(1), Parfin(2), P2);
    HBs->D0(Pardeb(3), Pardeb(4), P3);
    HBs->D0(Parfin(3), Parfin(4), P4);
    Pc->D0(p1, p2d1);
    Pc->D0(p2, p2d2);
    HBs->D0(p2d1.X(), p2d1.Y(), P7);
    HBs->D0(p2d2.X(), p2d2.Y(), P8);
    Ps->D0(p1, p2d1);
    Ps->D0(p2, p2d2);
    HGs->D0(p2d1.X(), p2d1.Y(), P5);
    HGs->D0(p2d2.X(), p2d2.Y(), P6);
    to1 = std::max(P1.Distance(P5) + P3.Distance(P7), tolreached);
    to2 = std::max(P2.Distance(P6) + P4.Distance(P8), tolreached);

    //////////////////////////////////////////////////////////////////////
    // storage in the DS of the intersection curve
    //////////////////////////////////////////////////////////////////////

    bool Isvtx1 = false;
    bool Isvtx2 = false;
    int  indice;

    if (nb == 1)
    {
      indpoint1 = stripe->IndexPoint(isfirst, 1);
      if (!CV1.IsVertex())
      {
        TopOpeBRepDS_Point& tpt = DStr.ChangePoint(indpoint1);
        tpt.Tolerance(std::max(tpt.Tolerance(), to1));
      }
      else
      {
        Isvtx1 = true;
      }
    }
    if (nb == nbface)
    {
      indpoint2 = stripe->IndexPoint(isfirst, 2);
      if (!CV2.IsVertex())
      {
        TopOpeBRepDS_Point& tpt = DStr.ChangePoint(indpoint2);
        tpt.Tolerance(std::max(tpt.Tolerance(), to2));
      }
      else
      {
        Isvtx2 = true;
      }
    }
    else
    {
      gp_Pnt             point = Cc->Value(Cc->LastParameter());
      TopOpeBRepDS_Point tpoint(point, to2);
      indpoint2 = DStr.AddPoint(tpoint);
    }

    if (nb != 1)
    {
      TopOpeBRepDS_Point& tpt = DStr.ChangePoint(indpoint1);
      tpt.Tolerance(std::max(tpt.Tolerance(), to1));
    }
    TopOpeBRepDS_Curve tcurv3d(Cc, tolreached);
    indcurve[nb - 1] = DStr.AddCurve(tcurv3d);

    Interfp1 = ChFi3d_FilPointInDS(TopAbs_FORWARD,
                                   indcurve[nb - 1],
                                   indpoint1,
                                   Cc->FirstParameter(),
                                   Isvtx1);

    Interfp2 = ChFi3d_FilPointInDS(TopAbs_REVERSED,
                                   indcurve[nb - 1],
                                   indpoint2,
                                   Cc->LastParameter(),
                                   Isvtx2);

    DStr.ChangeCurveInterferences(indcurve[nb - 1]).Append(Interfp1);
    DStr.ChangeCurveInterferences(indcurve[nb - 1]).Append(Interfp2);

    //////////////////////////////////////////////////////////////////////
    // storage for the face
    //////////////////////////////////////////////////////////////////////

    TopAbs_Orientation ori = TopAbs_FORWARD;
    orface                 = Face[nb - 1].Orientation();
    if (orface == orsurfdata)
    {
      orien = TopAbs::Reverse(orcourbe);
    }
    else
    {
      orien = orcourbe;
    }
    // limitation of edges of faces
    if (nb == 1)
    {
      int Iarc1 = DStr.AddShape(Edge[0]);
      Interfedge[0] =
        ChFi3d_FilPointInDS(CV1.TransitionOnArc(), Iarc1, indpoint1, paredge1, Isvtx1);
      // DStr.ChangeShapeInterferences(Edge[0]).Append(Interfp1);
    }
    if (nb == nbface)
    {
      int Iarc2 = DStr.AddShape(Edge[nb]);
      Interfedge[nb] =
        ChFi3d_FilPointInDS(CV2.TransitionOnArc(), Iarc2, indpoint2, paredge2, Isvtx2);
      // DStr.ChangeShapeInterferences(Edge[nb]).Append(Interfp2);
    }

    if (nb != nbface || oneintersection1 || oneintersection2)
    {
      if (nface == 3)
      {
        V1 = TopExp::FirstVertex(Edge[nb]);
        V2 = TopExp::LastVertex(Edge[nb]);
        if (containV(F1, V1) || containV(F2, V1))
        {
          ori = TopAbs_FORWARD;
        }
        else if (containV(F1, V2) || containV(F2, V2))
        {
          ori = TopAbs_REVERSED;
        }
        else
        {
          throw Standard_Failure("IntersectionAtEnd : pb orientation");
        }

        if (containV(F1, V1) && containV(F1, V2))
        {
          dist1 = (BRep_Tool::Pnt(V1)).Distance(BRep_Tool::Pnt(Vtx));
          dist2 = (BRep_Tool::Pnt(V2)).Distance(BRep_Tool::Pnt(Vtx));
          if (dist1 < dist2)
          {
            ori = TopAbs_FORWARD;
          }
          else
          {
            ori = TopAbs_REVERSED;
          }
        }
        if (containV(F2, V1) && containV(F2, V2))
        {
          dist1 = (BRep_Tool::Pnt(V1)).Distance(BRep_Tool::Pnt(Vtx));
          dist2 = (BRep_Tool::Pnt(V2)).Distance(BRep_Tool::Pnt(Vtx));
          if (dist1 < dist2)
          {
            ori = TopAbs_FORWARD;
          }
          else
          {
            ori = TopAbs_REVERSED;
          }
        }
      }
      else
      {
        if (TopExp::FirstVertex(Edge[nb]).IsSame(Vtx))
        {
          ori = TopAbs_FORWARD;
        }
        else
        {
          ori = TopAbs_REVERSED;
        }
      }
      if (!extend && !pastVtx && !(oneintersection1 || oneintersection2))
      {
        int Iarc2      = DStr.AddShape(Edge[nb]);
        Interfedge[nb] = ChFi3d_FilPointInDS(ori, Iarc2, indpoint2, paredge2);
        //  DStr.ChangeShapeInterferences(Edge[nb]).Append(Interfp2);
      }
      else
      {
        if (!(oneintersection1 || oneintersection2))
        {
          proledge[nb] = true;
        }
        int    indp1, indp2, ind;
        gp_Pnt pext;
        double ubid, vbid;
        pext = BRep_Tool::Pnt(Vtx);
        GeomAdaptor_Curve       cad;
        occ::handle<Geom_Curve> csau;
        if (!(oneintersection1 || oneintersection2))
        {
          cad.Load(cint);
          csau = cint;
        }
        else
        {
          csau                              = BRep_Tool::Curve(edgesau, ubid, vbid);
          occ::handle<Geom_BoundedCurve> C1 = occ::down_cast<Geom_BoundedCurve>(csau);
          if (oneintersection1 && extend)
          {
            if (!C1.IsNull())
            {
              gp_Pnt Pl;
              Pl = C1->Value(C1->LastParameter());
              // bool sens;
              sens = Pl.Distance(pext) < tolpt;
              GeomLib::ExtendCurveToPoint(C1, CV1.Point(), 1, sens != 0);
              csau = C1;
            }
          }
          else if (oneintersection2 && extend)
          {
            if (!C1.IsNull())
            {
              gp_Pnt Pl;
              Pl = C1->Value(C1->LastParameter());
              // bool sens;
              sens = Pl.Distance(pext) < tolpt;
              GeomLib::ExtendCurveToPoint(C1, CV2.Point(), 1, sens != 0);
              csau = C1;
            }
          }
          cad.Load(csau);
        }
        Extrema_ExtPC ext(pext, cad, tolpt);
        double        par1, par2, par, ParVtx;
        bool          vtx1 = false;
        bool          vtx2 = false;
        par1               = ext.Point(1).Parameter();
        ParVtx             = par1;
        if (oneintersection1 || oneintersection2)
        {
          if (oneintersection2)
          {
            pext = CV2.Point();
            ind  = indpoint2;
          }
          else
          {
            pext = CV1.Point();
            ind  = indpoint1;
          }
          Extrema_ExtPC ext2(pext, cad, tolpt);
          par2 = ext2.Point(1).Parameter();
        }
        else
        {
          par2 = paredge2;
          ind  = indpoint2;
        }
        if (par1 > par2)
        {
          indp1 = ind;
          indp2 = DStr.AddShape(Vtx);
          vtx2  = true;
          par   = par1;
          par1  = par2;
          par2  = par;
        }
        else
        {
          indp1 = DStr.AddShape(Vtx);
          indp2 = ind;
          vtx1  = true;
        }
        occ::handle<Geom_Curve> Ct = new Geom_TrimmedCurve(csau, par1, par2);
        TopAbs_Orientation      orient;
        Cc->D0(Cc->FirstParameter(), P1);
        Cc->D0(Cc->LastParameter(), P2);
        Ct->D0(Ct->FirstParameter(), P3);
        Ct->D0(Ct->LastParameter(), P4);
        if (P2.Distance(P3) < tolpt || P1.Distance(P4) < tolpt)
        {
          orient = orien;
        }
        else
        {
          orient = TopAbs::Reverse(orien);
        }
        if (oneintersection1 || oneintersection2)
        {
          indice = DStr.AddShape(Face[0]);
          if (extend)
          {
            DStr.SetNewSurface(Face[0], Sfacemoins1);
            ComputeCurve2d(Ct, faceprol[0], C2dint1);
          }
          else
          {
            TopoDS_Edge aLocalEdge = edgesau;
            if (edgesau.Orientation() != orient)
            {
              aLocalEdge.Reverse();
            }
            C2dint1 = BRep_Tool::CurveOnSurface(aLocalEdge, Face[0], ubid, vbid);
          }
        }
        else
        {
          indice = DStr.AddShape(Face[nb - 1]);
          if (!pastVtx)
          {
            DStr.SetNewSurface(Face[nb - 1], Sfacemoins1);
          }
        }
        //// for periodic 3d curves ////
        if (cad.IsPeriodic() && !C2dint1.IsNull())
        {
          gp_Pnt2d                      P2d = BRep_Tool::Parameters(Vtx, Face[0]);
          Geom2dAPI_ProjectPointOnCurve Projector(P2d, C2dint1);
          par          = Projector.LowerDistanceParameter();
          double shift = par - ParVtx;
          if (std::abs(shift) > Precision::Confusion())
          {
            par1 += shift;
            par2 += shift;
          }
        }
        ////////////////////////////////

        Ct = new Geom_TrimmedCurve(csau, par1, par2);
        if (oneintersection1 || oneintersection2)
        {
          tolex = 10 * BRep_Tool::Tolerance(edgesau);
        }
        if (extend)
        {
          occ::handle<GeomAdaptor_Surface> H1, H2;
          H1 = new GeomAdaptor_Surface(Sfacemoins1);
          if (Sface.IsNull())
          {
            tolex = std::max(tolex, ChFi3d_EvalTolReached(H1, C2dint1, H1, C2dint1, Ct));
          }
          else
          {
            H2    = new GeomAdaptor_Surface(Sface);
            tolex = std::max(tolex, ChFi3d_EvalTolReached(H1, C2dint1, H2, C2dint2, Ct));
          }
        }
        TopOpeBRepDS_Curve tcurv(Ct, tolex);
        int                indcurv;
        indcurv  = DStr.AddCurve(tcurv);
        Interfp1 = ChFi3d_FilPointInDS(TopAbs_FORWARD, indcurv, indp1, par1, vtx1);
        Interfp2 = ChFi3d_FilPointInDS(TopAbs_REVERSED, indcurv, indp2, par2, vtx2);
        DStr.ChangeCurveInterferences(indcurv).Append(Interfp1);
        DStr.ChangeCurveInterferences(indcurv).Append(Interfp2);

        Interfc = ChFi3d_FilCurveInDS(indcurv, indice, C2dint1, orient);
        DStr.ChangeShapeInterferences(indice).Append(Interfc);
        if (oneintersection1 || oneintersection2)
        {
          indice = DStr.AddShape(facesau);
          if (facesau.Orientation() == Face[0].Orientation())
          {
            orient = TopAbs::Reverse(orient);
          }
          if (extend)
          {
            ComputeCurve2d(Ct, faceprol[1], C2dint2);
          }
          else
          {
            TopoDS_Edge aLocalEdge = edgesau;
            if (edgesau.Orientation() != orient)
            {
              aLocalEdge.Reverse();
            }
            C2dint2 = BRep_Tool::CurveOnSurface(aLocalEdge, facesau, ubid, vbid);
            // Reverse for case of edgesau on closed surface (Face[0] is equal to facesau)
          }
        }
        else
        {
          indice = DStr.AddShape(Face[nb]);
          if (!pastVtx)
          {
            DStr.SetNewSurface(Face[nb], Sface);
          }
          if (Face[nb].Orientation() == Face[nb - 1].Orientation())
          {
            orient = TopAbs::Reverse(orient);
          }
        }
        if (!bordlibre)
        {
          Interfc = ChFi3d_FilCurveInDS(indcurv, indice, C2dint2, orient);
          DStr.ChangeShapeInterferences(indice).Append(Interfc);
        }
      }
    }

    if (checkShrink
        && IsShrink(Ps, p1, p2, checkShrParam, isUShrink, Precision::Parametric(tolreached)))
    {
      shrink[nb - 1] = 1;
      // store section face-chamf curve for previous SurfData
      // Suppose Fd and SDprev are parametrized similarly
      if (!isShrink)
      { // first time
        const ChFiDS_FaceInterference& Fi = SDprev->InterferenceOnS1();
        gp_Pnt2d                       UV = Fi.PCurveOnSurf()->Value(Fi.Parameter(isfirst));
        prevSDParam                       = isUShrink ? UV.X() : UV.Y();
      }
      gp_Pnt2d UV1 = p2d1, UV2 = p2d2;
      UV1.SetCoord(isUShrink ? 1 : 2, prevSDParam);
      UV2.SetCoord(isUShrink ? 1 : 2, prevSDParam);
      double aTolreached;
      ChFi3d_ComputePCurv(Cc,
                          UV1,
                          UV2,
                          Ps,
                          DStr.Surface(SDprev->Surf()).Surface(),
                          p1,
                          p2,
                          tolapp3d,
                          aTolreached);
      TopOpeBRepDS_Curve& TCurv = DStr.ChangeCurve(indcurve[nb - 1]);
      TCurv.Tolerance(std::max(TCurv.Tolerance(), aTolreached));

      InterfPS[nb - 1] = ChFi3d_FilCurveInDS(indcurve[nb - 1], IsurfPrev, Ps, orcourbe);
      DStr.ChangeSurfaceInterferences(IsurfPrev).Append(InterfPS[nb - 1]);

      if (isOnSame2)
      {
        midP2d    = p2d2;
        midIpoint = indpoint2;
      }
      else if (!isShrink)
      {
        midP2d    = p2d1;
        midIpoint = indpoint1;
      }
      isShrink = true;
    } // end if shrink

    indice           = DStr.AddShape(Face[nb - 1]);
    InterfPC[nb - 1] = ChFi3d_FilCurveInDS(indcurve[nb - 1], indice, Pc, orien);
    if (!shrink[nb - 1])
    {
      InterfPS[nb - 1] = ChFi3d_FilCurveInDS(indcurve[nb - 1], Isurf, Ps, orcourbe);
    }
    indpoint1 = indpoint2;

  } // end loop on faces being intersected with ChFi

  if (isOnSame1)
  {
    CV1.Reset();
  }
  if (isOnSame2)
  {
    CV2.Reset();
  }

  for (nb = 1; nb <= nbface; nb++)
  {
    int indice = DStr.AddShape(Face[nb - 1]);
    DStr.ChangeShapeInterferences(indice).Append(InterfPC[nb - 1]);
    if (!shrink[nb - 1])
    {
      DStr.ChangeSurfaceInterferences(Isurf).Append(InterfPS[nb - 1]);
    }
    if (!proledge[nb - 1])
    {
      DStr.ChangeShapeInterferences(Edge[nb - 1]).Append(Interfedge[nb - 1]);
    }
  }
  DStr.ChangeShapeInterferences(Edge[nbface]).Append(Interfedge[nbface]);
  // The edge of Vtx the arc carries on lies whole under the fillet: it is
  // cut at its far end, the arc's vertex, and nothing of it kept.
  if (!EdgePast.IsNull())
  {
    TopoDS_Vertex aV1, aV2;
    TopExp::Vertices(EdgePast, aV1, aV2);
    const bool           isFwd  = aV1.IsSame(Vtx);
    const TopoDS_Vertex& aVfar  = isFwd ? aV2 : aV1;
    const int            IPast  = DStr.AddShape(EdgePast);
    DStr.ChangeShapeInterferences(IPast).Append(
      ChFi3d_FilPointInDS(isFwd ? TopAbs_FORWARD : TopAbs_REVERSED,
                          IPast,
                          DStr.AddShape(aVfar),
                          BRep_Tool::Parameter(aVfar, EdgePast),
                          true));
  }
  if (onsOS != 0)
  {
    StoreLineOverSplit(DStr,
                       stripe,
                       Fd,
                       isfirst,
                       onsOS,
                       Isurf,
                       HGs,
                       LinOS,
                       EsplitOS,
                       FnOS,
                       wQOS,
                       wPOS,
                       parQOS,
                       stripe->IndexPoint(isfirst, onsOS),
                       10 * tolapp3d,
                       myRegul);
  }

  if (!isShrink)
  {
    stripe->InDS(isfirst);
  }
  else
  {
    // compute curves for !<isfirst> end of <Fd> and <isfirst> end of previous <SurfData>

    // for Fd
    // Bnd_Box box;
    gp_Pnt2d UV, UV1 = midP2d, UV2 = midP2d;
    if (isOnSame1)
    {
      UV = UV2 = Fi1.PCurveOnSurf()->Value(Fi1.Parameter(!isfirst));
    }
    else
    {
      UV = UV1 = Fi2.PCurveOnSurf()->Value(Fi2.Parameter(!isfirst));
    }
    double                    aTolreached;
    occ::handle<Geom_Curve>   C3d;
    occ::handle<Geom_Surface> aSurf = DStr.Surface(Fd->Surf()).Surface();
    // box.Add(aSurf->Value(UV.X(), UV.Y()));

    ChFi3d_ComputeArete(CV1,
                        UV1,
                        CV2,
                        UV2,
                        aSurf, // in
                        C3d,
                        Ps,
                        p1,
                        p2,
                        tolapp3d,
                        tol2d,
                        aTolreached,
                        0); // out except tolers

    indpoint1 = indpoint2 = midIpoint;
    gp_Pnt point;
    if (isOnSame1)
    {
      point = C3d->Value(p2);
      TopOpeBRepDS_Point tpoint(point, aTolreached);
      indpoint2 = DStr.AddPoint(tpoint);
      UV        = Ps->Value(p2);
    }
    else
    {
      point = C3d->Value(p1);
      TopOpeBRepDS_Point tpoint(point, aTolreached);
      indpoint1 = DStr.AddPoint(tpoint);
      UV        = Ps->Value(p1);
    }
    // box.Add(point);
    // box.Add(aSurf->Value(UV.X(), UV.Y()));

    TopOpeBRepDS_Curve Crv   = TopOpeBRepDS_Curve(C3d, aTolreached);
    int                Icurv = DStr.AddCurve(Crv);
    Interfp1                 = ChFi3d_FilPointInDS(TopAbs_FORWARD, Icurv, indpoint1, p1, false);
    Interfp2                 = ChFi3d_FilPointInDS(TopAbs_REVERSED, Icurv, indpoint2, p2, false);
    Interfc                  = ChFi3d_FilCurveInDS(Icurv, Isurf, Ps, orcourbe);
    DStr.ChangeCurveInterferences(Icurv).Append(Interfp1);
    DStr.ChangeCurveInterferences(Icurv).Append(Interfp2);
    DStr.ChangeSurfaceInterferences(Isurf).Append(Interfc);

    // for SDprev
    aSurf = DStr.Surface(SDprev->Surf()).Surface();
    UV1.SetCoord(isUShrink ? 1 : 2, prevSDParam);
    UV2.SetCoord(isUShrink ? 1 : 2, prevSDParam);

    ChFi3d_ComputePCurv(C3d, UV1, UV2, Pc, aSurf, p1, p2, tolapp3d, aTolreached);

    Crv.Tolerance(std::max(Crv.Tolerance(), aTolreached));
    Interfc = ChFi3d_FilCurveInDS(Icurv, IsurfPrev, Pc, TopAbs::Reverse(orcourbe));
    DStr.ChangeSurfaceInterferences(IsurfPrev).Append(Interfc);

    // UV = isOnSame1 ? UV2 : UV1;
    // box.Add(aSurf->Value(UV.X(), UV.Y()));
    // UV = Ps->Value(isOnSame1 ? p2 : p1);
    // box.Add(aSurf->Value(UV.X(), UV.Y()));
    // ChFi3d_SetPointTolerance(DStr,box, isOnSame1 ? indpoint2 : indpoint1);

    // to process properly this case in ChFi3d_FilDS()
    stripe->InDS(isfirst, 2);
    Fd->ChangeInterference(isOnSame1 ? 2 : 1).SetLineIndex(0);
    ChFiDS_CommonPoint& CPprev1 = SDprev->ChangeVertex(isfirst, isOnSame1 ? 2 : 1);
    ChFiDS_CommonPoint& CPlast1 = Fd->ChangeVertex(isfirst, isOnSame1 ? 2 : 1);
    ChFiDS_CommonPoint& CPlast2 = Fd->ChangeVertex(!isfirst, isOnSame1 ? 2 : 1);
    if (CPprev1.IsOnArc())
    {
      CPlast1 = CPprev1;
      CPprev1.Reset();
      CPprev1.SetPoint(CPlast1.Point());
      CPlast2.Reset();
      CPlast2.SetPoint(CPlast1.Point());
    }

    // in shrink case, self intersection is possible at <midIpoint>,
    // eval its tolerance intersecting Ps and Pcurve at end.
    // Find end curves closest to shrinked part
    for (nb = 0; nb < nbface; nb++)
    {
      if (isOnSame1 ? shrink[nb + 1] : !shrink[nb])
      {
        break;
      }
    }
    occ::handle<Geom_Curve>   Cend  = DStr.Curve(indcurve[nb]).Curve();
    occ::handle<Geom2d_Curve> PCend = InterfPS[nb]->PCurve();
    // point near which self intersection may occur
    TopOpeBRepDS_Point& Pds   = DStr.ChangePoint(midIpoint);
    const gp_Pnt&       Pvert = Pds.Point();
    double              tol   = Pds.Tolerance();

    Geom2dAdaptor_Curve PC1(Ps), PC2(PCend);
    Geom2dInt_GInter    Intersector(PC1, PC2, Precision::PConfusion(), Precision::PConfusion());
    if (!Intersector.IsDone())
    {
      return;
    }
    for (nb = 1; nb <= Intersector.NbPoints(); nb++)
    {
      const IntRes2d_IntersectionPoint& ip   = Intersector.Point(nb);
      gp_Pnt                            Pint = C3d->Value(ip.ParamOnFirst());
      tol                                    = std::max(tol, Pvert.Distance(Pint));
      Pint                                   = Cend->Value(ip.ParamOnSecond());
      tol                                    = std::max(tol, Pvert.Distance(Pint));
    }
    for (nb = 1; nb <= Intersector.NbSegments(); nb++)
    {
      const IntRes2d_IntersectionSegment& is = Intersector.Segment(nb);
      if (is.HasFirstPoint())
      {
        const IntRes2d_IntersectionPoint& ip   = is.FirstPoint();
        gp_Pnt                            Pint = C3d->Value(ip.ParamOnFirst());
        tol                                    = std::max(tol, Pvert.Distance(Pint));
        Pint                                   = Cend->Value(ip.ParamOnSecond());
        tol                                    = std::max(tol, Pvert.Distance(Pint));
      }
      if (is.HasLastPoint())
      {
        const IntRes2d_IntersectionPoint& ip   = is.LastPoint();
        gp_Pnt                            Pint = C3d->Value(ip.ParamOnFirst());
        tol                                    = std::max(tol, Pvert.Distance(Pint));
        Pint                                   = Cend->Value(ip.ParamOnSecond());
        tol                                    = std::max(tol, Pvert.Distance(Pint));
      }
    }
    Pds.Tolerance(tol);
  }
}

//  Modified by Sergey KHROMOV - Thu Apr 11 12:23:40 2002 Begin

//=======================================================================
// function : PerformMoreSurfdata
// purpose  :  determine intersections at end on several surfdata
//=======================================================================
void ChFi3d_Builder::PerformMoreSurfdata(const int Index)
{
  TopOpeBRepDS_DataStructure&                         DStr       = myDS->ChangeDS();
  const NCollection_List<occ::handle<ChFiDS_Stripe>>& aLOfStripe = myVDataMap(Index);
  occ::handle<ChFiDS_Stripe>                          aStripe;
  occ::handle<ChFiDS_Spine>                           aSpine;
  double                                              aTol3d = 1.e-4;

  if (aLOfStripe.IsEmpty())
  {
    return;
  }

  aStripe = aLOfStripe.First();
  aSpine  = aStripe->Spine();

  NCollection_Sequence<occ::handle<ChFiDS_SurfData>>& aSeqSurfData =
    aStripe->ChangeSetOfSurfData()->ChangeSequence();
  const TopoDS_Vertex&         aVtx    = myVDataMap.FindKey(Index);
  int                          aSens   = 0;
  int                          anInd   = ChFi3d_IndexOfSurfData(aVtx, aStripe, aSens);
  bool                         isFirst = (aSens == 1);
  int                          anIndPrev;
  occ::handle<ChFiDS_SurfData> aSurfData;
  ChFiDS_CommonPoint           aCP1;
  ChFiDS_CommonPoint           aCP2;

  aSurfData = aSeqSurfData.Value(anInd);

  aCP1 = aSurfData->Vertex(isFirst, 1);
  aCP2 = aSurfData->Vertex(isFirst, 2);

  occ::handle<Geom_Surface> aSurfPrev;
  occ::handle<Geom_Surface> aSurf;
  TopoDS_Face               aFace;
  TopoDS_Face               aNeighborFace;

  FindFace(aVtx, aCP1, aCP2, aFace);
  aSurfPrev = BRep_Tool::Surface(aFace);

  if (aSens == 1)
  {
    anIndPrev = anInd + 1;
  }
  else
  {
    anIndPrev = anInd - 1;
  }

  TopoDS_Edge                              anArc1;
  TopoDS_Edge                              anArc2;
  NCollection_List<TopoDS_Shape>::Iterator anIter(myVEMap(aVtx));
  bool                                     isFound = false;

  for (; anIter.More() && !isFound; anIter.Next())
  {
    anArc1 = TopoDS::Edge(anIter.Value());

    if (containE(aFace, anArc1))
    {
      isFound = true;
    }
  }

  isFound = false;
  anIter.Initialize(myVEMap(aVtx));

  for (; anIter.More() && !isFound; anIter.Next())
  {
    anArc2 = TopoDS::Edge(anIter.Value());

    if (containE(aFace, anArc2) && !anArc2.IsSame(anArc1))
    {
      isFound = true;
    }
  }

  // determination of common points aCP1onArc, aCP2onArc and aCP2NotonArc
  // aCP1onArc    is the point on arc of index anInd
  // aCP2onArc    is the point on arc of index anIndPrev
  // aCP2NotonArc is the point of index anIndPrev which is not on arc.

  bool               is1stCP1OnArc;
  bool               is2ndCP1OnArc;
  ChFiDS_CommonPoint aCP1onArc;
  ChFiDS_CommonPoint aCP2onArc;
  ChFiDS_CommonPoint aCP2NotonArc;

  aSurfData = aSeqSurfData.Value(anIndPrev);
  aCP1      = aSurfData->Vertex(isFirst, 1);
  aCP2      = aSurfData->Vertex(isFirst, 2);

  if (aCP1.IsOnArc() && (aCP1.Arc().IsSame(anArc1) || aCP1.Arc().IsSame(anArc2)))
  {
    aCP2onArc     = aCP1;
    aCP2NotonArc  = aCP2;
    is2ndCP1OnArc = true;
  }
  else if (aCP2.IsOnArc() && (aCP2.Arc().IsSame(anArc1) || aCP2.Arc().IsSame(anArc2)))
  {
    aCP2onArc     = aCP2;
    aCP2NotonArc  = aCP1;
    is2ndCP1OnArc = false;
  }
  else
  {
    return;
  }

  aSurfData = aSeqSurfData.Value(anInd);
  aCP1      = aSurfData->Vertex(isFirst, 1);
  aCP2      = aSurfData->Vertex(isFirst, 2);

  if (aCP1.Point().Distance(aCP2onArc.Point()) <= aTol3d)
  {
    aCP1onArc     = aCP2;
    is1stCP1OnArc = false;
  }
  else
  {
    aCP1onArc     = aCP1;
    is1stCP1OnArc = true;
  }

  if (!aCP1onArc.IsOnArc())
  {
    return;
  }

  // determination of neighbor surface
  int indSurface;
  if (is1stCP1OnArc)
  {
    indSurface = myListStripe.First()->SetOfSurfData()->Value(anInd)->IndexOfS1();
  }
  else
  {
    indSurface = myListStripe.First()->SetOfSurfData()->Value(anInd)->IndexOfS2();
  }

  aNeighborFace = TopoDS::Face(myDS->Shape(indSurface));

  // calculation of intersections
  occ::handle<Geom_Curve>   aCracc;
  occ::handle<Geom2d_Curve> aPCurv1;
  double                    aParf;
  double                    aParl;
  double                    aTolReached;

  aSurfData = aSeqSurfData.Value(anInd);

  if (isFirst)
  {
    ChFi3d_ComputeArete(aSurfData->VertexLastOnS1(),
                        aSurfData->InterferenceOnS1().PCurveOnSurf()->Value(
                          aSurfData->InterferenceOnS1().LastParameter()),
                        aSurfData->VertexLastOnS2(),
                        aSurfData->InterferenceOnS2().PCurveOnSurf()->Value(
                          aSurfData->InterferenceOnS2().LastParameter()),
                        DStr.Surface(aSurfData->Surf()).Surface(),
                        aCracc,
                        aPCurv1,
                        aParf,
                        aParl,
                        aTol3d,
                        tol2d,
                        aTolReached,
                        0);
  }
  else
  {
    ChFi3d_ComputeArete(aSurfData->VertexFirstOnS1(),
                        aSurfData->InterferenceOnS1().PCurveOnSurf()->Value(
                          aSurfData->InterferenceOnS1().FirstParameter()),
                        aSurfData->VertexFirstOnS2(),
                        aSurfData->InterferenceOnS2().PCurveOnSurf()->Value(
                          aSurfData->InterferenceOnS2().FirstParameter()),
                        DStr.Surface(aSurfData->Surf()).Surface(),
                        aCracc,
                        aPCurv1,
                        aParf,
                        aParl,
                        aTol3d,
                        tol2d,
                        aTolReached,
                        0);
  }

  // calculation of the index of the line on anInd.
  // aPClineOnSurf is the pcurve on anInd.
  // aPClineOnFace is the pcurve on face.
  ChFiDS_FaceInterference aFI;

  if (is1stCP1OnArc)
  {
    aFI = aSurfData->InterferenceOnS1();
  }
  else
  {
    aFI = aSurfData->InterferenceOnS2();
  }

  occ::handle<Geom_Curve>   aCline;
  occ::handle<Geom2d_Curve> aPClineOnSurf;
  occ::handle<Geom2d_Curve> aPClineOnFace;
  int                       indLine;

  indLine       = aFI.LineIndex();
  aCline        = DStr.Curve(aFI.LineIndex()).Curve();
  aPClineOnSurf = aFI.PCurveOnSurf();
  aPClineOnFace = aFI.PCurveOnFace();

  // intersection between the SurfData number anInd and the Face aFace.
  // Obtaining of curves aCint1, aPCint11 and aPCint12.
  aSurf = DStr.Surface(aSurfData->Surf()).Surface();

  GeomInt_IntSS                    anInterSS(aSurfPrev, aSurf, 1.e-7, true, true, true);
  occ::handle<Geom_Curve>          aCint1;
  occ::handle<Geom2d_Curve>        aPCint11;
  occ::handle<Geom2d_Curve>        aPCint12;
  occ::handle<GeomAdaptor_Surface> H1      = new GeomAdaptor_Surface(aSurfPrev);
  occ::handle<GeomAdaptor_Surface> H2      = new GeomAdaptor_Surface(aSurf);
  double                           aTolex1 = 0.;
  int                              i;
  gp_Pnt                           aPext1;
  gp_Pnt                           aPext2;
  gp_Pnt                           aPext;
  bool                             isPextFound;

  if (!anInterSS.IsDone())
  {
    return;
  }

  isFound = false;

  for (i = 1; i <= anInterSS.NbLines() && !isFound; i++)
  {
    aCint1   = anInterSS.Line(i);
    aPCint11 = anInterSS.LineOnS1(i);
    aPCint12 = anInterSS.LineOnS2(i);
    aTolex1  = ChFi3d_EvalTolReached(H1, aPCint11, H2, aPCint12, aCint1);

    aCint1->D0(aCint1->FirstParameter(), aPext1);
    aCint1->D0(aCint1->LastParameter(), aPext2);

    //  Modified by skv - Mon Jun  7 18:38:57 2004 OCC5898 Begin
    //     if (aPext1.Distance(aCP1onArc.Point()) <= aTol3d ||
    // 	aPext2.Distance(aCP1onArc.Point()))
    if (aPext1.Distance(aCP1onArc.Point()) <= aTol3d
        || aPext2.Distance(aCP1onArc.Point()) <= aTol3d)
    {
      //  Modified by skv - Mon Jun  7 18:38:58 2004 OCC5898 End
      isFound = true;
    }
  }

  if (!isFound)
  {
    return;
  }

  if (aPext1.Distance(aCP2onArc.Point()) > aTol3d && aPext1.Distance(aCP1onArc.Point()) > aTol3d)
  {
    aPext       = aPext1;
    isPextFound = true;
  }
  else if (aPext2.Distance(aCP2onArc.Point()) > aTol3d
           && aPext2.Distance(aCP1onArc.Point()) > aTol3d)
  {
    aPext       = aPext2;
    isPextFound = true;
  }
  else
  {
    isPextFound = false;
  }

  bool   isDoSecondSection = false;
  double aPar              = 0.;

  if (isPextFound)
  {
    GeomAdaptor_Curve aCad(aCracc);
    Extrema_ExtPC     anExt(aPext, aCad, aTol3d);

    if (!anExt.IsDone())
    {
      return;
    }

    isFound = false;
    for (i = 1; i <= anExt.NbExt() && !isFound; i++)
    {
      if (anExt.IsMin(i))
      {
        gp_Pnt aProjPnt = anExt.Point(i).Value();

        if (aPext.Distance(aProjPnt) <= aTol3d)
        {
          aPar              = anExt.Point(i).Parameter();
          isDoSecondSection = true;
        }
      }
    }
  }

  occ::handle<Geom_Curve> aTrCracc;
  TopAbs_Orientation      anOrSD1;
  TopAbs_Orientation      anOrSD2;
  int                     indShape;

  anOrSD1   = aSurfData->Orientation();
  aSurfData = aSeqSurfData.Value(anIndPrev);
  anOrSD2   = aSurfData->Orientation();
  aSurf     = DStr.Surface(aSurfData->Surf()).Surface();

  // The following variables will be used if isDoSecondSection is true
  occ::handle<Geom_Curve>   aCint2;
  occ::handle<Geom2d_Curve> aPCint21;
  occ::handle<Geom2d_Curve> aPCint22;
  double                    aTolex2 = 0.;

  if (isDoSecondSection)
  {
    double aPar1;

    aCracc->D0(aCracc->FirstParameter(), aPext1);

    if (aPext1.Distance(aCP2NotonArc.Point()) <= aTol3d)
    {
      aPar1 = aCracc->FirstParameter();
    }
    else
    {
      aPar1 = aCracc->LastParameter();
    }

    if (aPar1 < aPar)
    {
      aTrCracc = new Geom_TrimmedCurve(aCracc, aPar1, aPar);
    }
    else
    {
      aTrCracc = new Geom_TrimmedCurve(aCracc, aPar, aPar1);
    }

    // Second section
    GeomInt_IntSS anInterSS2(aSurfPrev, aSurf, 1.e-7, true, true, true);

    if (!anInterSS2.IsDone())
    {
      return;
    }

    H1 = new GeomAdaptor_Surface(aSurfPrev);
    H2 = new GeomAdaptor_Surface(aSurf);

    isFound = false;

    for (i = 1; i <= anInterSS2.NbLines() && !isFound; i++)
    {
      aCint2   = anInterSS2.Line(i);
      aPCint21 = anInterSS2.LineOnS1(i);
      aPCint22 = anInterSS2.LineOnS2(i);
      aTolex2  = ChFi3d_EvalTolReached(H1, aPCint21, H2, aPCint22, aCint2);

      aCint2->D0(aCint2->FirstParameter(), aPext1);
      aCint2->D0(aCint2->LastParameter(), aPext2);

      if (aPext1.Distance(aCP2onArc.Point()) <= aTol3d
          || aPext2.Distance(aCP2onArc.Point()) <= aTol3d)
      {
        isFound = true;
      }
    }

    if (!isFound)
    {
      return;
    }
  }
  else
  {
    aTrCracc = new Geom_TrimmedCurve(aCracc, aCracc->FirstParameter(), aCracc->LastParameter());
  }

  // Storage of the data structure

  // calculation of the orientation of line of surfdata number
  // anIndPrev which contains aCP2onArc

  occ::handle<Geom2d_Curve> aPCraccS = GeomProjLib::Curve2d(aTrCracc, aSurf);

  if (is2ndCP1OnArc)
  {
    aFI      = aSurfData->InterferenceOnS1();
    indShape = aSurfData->IndexOfS1();
  }
  else
  {
    aFI      = aSurfData->InterferenceOnS2();
    indShape = aSurfData->IndexOfS2();
  }

  if (indShape <= 0)
  {
    return;
  }

  TopAbs_Orientation aCurOrient;

  aCurOrient = DStr.Shape(indShape).Orientation();
  aCurOrient = TopAbs::Compose(aCurOrient, aSurfData->Orientation());
  aCurOrient = TopAbs::Compose(TopAbs::Reverse(aFI.Transition()), aCurOrient);

  // Filling the data structure
  aSurfData = aSeqSurfData.Value(anInd);

  TopOpeBRepDS_Point aPtCP1(aCP1onArc.Point(), aCP1onArc.Tolerance());
  int                indCP1onArc = DStr.AddPoint(aPtCP1);
  int                indSurf1    = aSurfData->Surf();
  int                indArc1     = DStr.AddShape(aCP1onArc.Arc());
  int                indSol      = aStripe->SolidIndex();

  occ::handle<TopOpeBRepDS_CurvePointInterference> anInterfp1;
  occ::handle<TopOpeBRepDS_CurvePointInterference> anInterfp2;

  anInterfp1 = ChFi3d_FilPointInDS(aCP1onArc.TransitionOnArc(),
                                   indArc1,
                                   indCP1onArc,
                                   aCP1onArc.ParameterOnArc());
  DStr.ChangeShapeInterferences(aCP1onArc.Arc()).Append(anInterfp1);

  NCollection_List<occ::handle<TopOpeBRepDS_Interference>>& SolidInterfs =
    DStr.ChangeShapeInterferences(indSol);
  occ::handle<TopOpeBRepDS_SolidSurfaceInterference> SSI =
    new TopOpeBRepDS_SolidSurfaceInterference(TopOpeBRepDS_Transition(anOrSD1),
                                              TopOpeBRepDS_SOLID,
                                              indSol,
                                              TopOpeBRepDS_SURFACE,
                                              indSurf1);
  SolidInterfs.Append(SSI);

  // deletion of Surface Data.
  aSeqSurfData.Remove(anInd);

  if (!isFirst)
  {
    anInd--;
  }

  aSurfData = aSeqSurfData.Value(anInd);

  // definition of indices of common points in Data Structure

  int indCP2onArc;
  int indCP2NotonArc;

  if (is2ndCP1OnArc)
  {
    aStripe->SetIndexPoint(ChFi3d_IndexPointInDS(aCP2onArc, DStr), isFirst, 1);
    aStripe->SetIndexPoint(ChFi3d_IndexPointInDS(aCP2NotonArc, DStr), isFirst, 2);

    if (isFirst)
    {
      indCP2onArc    = aStripe->IndexFirstPointOnS1();
      indCP2NotonArc = aStripe->IndexFirstPointOnS2();
    }
    else
    {
      indCP2onArc    = aStripe->IndexLastPointOnS1();
      indCP2NotonArc = aStripe->IndexLastPointOnS2();
    }
  }
  else
  {
    aStripe->SetIndexPoint(ChFi3d_IndexPointInDS(aCP2onArc, DStr), isFirst, 2);
    aStripe->SetIndexPoint(ChFi3d_IndexPointInDS(aCP2NotonArc, DStr), isFirst, 1);

    if (isFirst)
    {
      indCP2onArc    = aStripe->IndexFirstPointOnS2();
      indCP2NotonArc = aStripe->IndexFirstPointOnS1();
    }
    else
    {
      indCP2onArc    = aStripe->IndexLastPointOnS2();
      indCP2NotonArc = aStripe->IndexLastPointOnS1();
    }
  }

  int    indPoint1;
  int    indPoint2;
  gp_Pnt aPoint1;
  gp_Pnt aPoint2;

  if (is2ndCP1OnArc)
  {
    aFI      = aSurfData->InterferenceOnS1();
    indShape = aSurfData->IndexOfS1();
  }
  else
  {
    aFI      = aSurfData->InterferenceOnS2();
    indShape = aSurfData->IndexOfS2();
  }

  gp_Pnt2d                                           aP2d;
  occ::handle<TopOpeBRepDS_SurfaceCurveInterference> anInterfc;
  TopAbs_Orientation                                 anOrSurf = aCurOrient;
  TopAbs_Orientation                                 anOrFace = aFace.Orientation();
  int                                                indaFace = DStr.AddShape(aFace);
  int                                                indPoint = indCP2onArc;
  int                                                indCurve;

  aFI.PCurveOnFace()->D0(aFI.LastParameter(), aP2d);
  occ::handle<Geom_Surface> Stemp2 = BRep_Tool::Surface(TopoDS::Face(DStr.Shape(indShape)));
  Stemp2->D0(aP2d.X(), aP2d.Y(), aPoint2);
  aFI.PCurveOnFace()->D0(aFI.FirstParameter(), aP2d);
  Stemp2->D0(aP2d.X(), aP2d.Y(), aPoint1);

  if (isDoSecondSection)
  {
    TopOpeBRepDS_Point tpoint(aPext, aTolex2);
    TopOpeBRepDS_Curve tcint2(aCint2, aTolex2);

    indPoint = DStr.AddPoint(tpoint);
    indCurve = DStr.AddCurve(tcint2);

    aCint2->D0(aCint2->FirstParameter(), aPext1);
    aCint2->D0(aCint2->LastParameter(), aPext2);

    if (aPext1.Distance(aPext) <= aTol3d)
    {
      indPoint1 = indPoint;
      indPoint2 = indCP2onArc;
    }
    else
    {
      indPoint1 = indCP2onArc;
      indPoint2 = indPoint;
    }

    // define the orientation of aCint2
    if (aPext1.Distance(aPoint2) > aTol3d && aPext2.Distance(aPoint1) > aTol3d)
    {
      anOrSurf = TopAbs::Reverse(anOrSurf);
    }

    // ---------------------------------------------------------------
    // storage of aCint2
    anInterfp1 = ChFi3d_FilPointInDS(TopAbs_FORWARD, indCurve, indPoint1, aCint2->FirstParameter());
    anInterfp2 = ChFi3d_FilPointInDS(TopAbs_REVERSED, indCurve, indPoint2, aCint2->LastParameter());
    DStr.ChangeCurveInterferences(indCurve).Append(anInterfp1);
    DStr.ChangeCurveInterferences(indCurve).Append(anInterfp2);

    // interference of aCint2 on the SurfData number anIndPrev
    anInterfc = ChFi3d_FilCurveInDS(indCurve, aSurfData->Surf(), aPCint22, anOrSurf);

    DStr.ChangeSurfaceInterferences(aSurfData->Surf()).Append(anInterfc);
    // interference of aCint2 on aFace

    if (anOrFace == anOrSD2)
    {
      anOrFace = TopAbs::Reverse(anOrSurf);
    }
    else
    {
      anOrFace = anOrSurf;
    }

    anInterfc = ChFi3d_FilCurveInDS(indCurve, indaFace, aPCint21, anOrFace);
    DStr.ChangeShapeInterferences(indaFace).Append(anInterfc);
  }

  aTrCracc->D0(aTrCracc->FirstParameter(), aPext1);
  aTrCracc->D0(aTrCracc->LastParameter(), aPext2);
  if (aPext1.Distance(aCP2NotonArc.Point()) <= aTol3d)
  {
    indPoint1 = indCP2NotonArc;
    indPoint2 = indPoint;
  }
  else
  {
    indPoint1 = indPoint;
    indPoint2 = indCP2NotonArc;
  }

  // Define the orientation of aTrCracc
  bool   isToReverse;
  gp_Pnt aP1;
  gp_Pnt aP2;
  gp_Pnt aP3;
  gp_Pnt aP4;

  if (isDoSecondSection)
  {
    aTrCracc->D0(aTrCracc->FirstParameter(), aP1);
    aTrCracc->D0(aTrCracc->LastParameter(), aP2);
    aCint2->D0(aCint2->FirstParameter(), aP3);
    aCint2->D0(aCint2->LastParameter(), aP4);
    isToReverse = (aP1.Distance(aP4) > aTol3d && aP2.Distance(aP3) > aTol3d);
  }
  else
  {
    isToReverse = (aPext1.Distance(aPoint2) > aTol3d && aPext2.Distance(aPoint1) > aTol3d);
  }

  if (isToReverse)
  {
    anOrSurf = TopAbs::Reverse(anOrSurf);
  }

  // ---------------------------------------------------------------
  // storage of aTrCracc
  TopOpeBRepDS_Curve tct2(aTrCracc, aTolReached);

  indCurve   = DStr.AddCurve(tct2);
  anInterfp1 = ChFi3d_FilPointInDS(TopAbs_FORWARD, indCurve, indPoint1, aTrCracc->FirstParameter());
  anInterfp2 = ChFi3d_FilPointInDS(TopAbs_REVERSED, indCurve, indPoint2, aTrCracc->LastParameter());
  DStr.ChangeCurveInterferences(indCurve).Append(anInterfp1);
  DStr.ChangeCurveInterferences(indCurve).Append(anInterfp2);

  // interference of aTrCracc on the SurfData number anIndPrev

  anInterfc = ChFi3d_FilCurveInDS(indCurve, aSurfData->Surf(), aPCraccS, anOrSurf);
  DStr.ChangeSurfaceInterferences(aSurfData->Surf()).Append(anInterfc);
  aStripe->InDS(isFirst);

  // interference of aTrCracc on the SurfData number anInd
  if (anOrSD1 == anOrSD2)
  {
    anOrSurf = TopAbs::Reverse(anOrSurf);
  }

  anInterfc = ChFi3d_FilCurveInDS(indCurve, indSurf1, aPCurv1, anOrSurf);
  DStr.ChangeSurfaceInterferences(indSurf1).Append(anInterfc);

  // ---------------------------------------------------------------
  // storage of aCint1

  aCint1->D0(aCint1->FirstParameter(), aPext1);
  if (aPext1.Distance(aCP1onArc.Point()) <= aTol3d)
  {
    indPoint1 = indCP1onArc;
    indPoint2 = indPoint;
  }
  else
  {
    indPoint1 = indPoint;
    indPoint2 = indCP1onArc;
  }

  //  definition of the orientation of aCint1

  aCint1->D0(aCint1->FirstParameter(), aP1);
  aCint1->D0(aCint1->LastParameter(), aP2);
  aTrCracc->D0(aTrCracc->FirstParameter(), aP3);
  aTrCracc->D0(aTrCracc->LastParameter(), aP4);

  if (aP1.Distance(aP4) > aTol3d && aP2.Distance(aP3) > aTol3d)
  {
    anOrSurf = TopAbs::Reverse(anOrSurf);
  }

  TopOpeBRepDS_Curve aTCint1(aCint1, aTolex1);
  indCurve   = DStr.AddCurve(aTCint1);
  anInterfp1 = ChFi3d_FilPointInDS(TopAbs_FORWARD, indCurve, indPoint1, aCint1->FirstParameter());
  anInterfp2 = ChFi3d_FilPointInDS(TopAbs_REVERSED, indCurve, indPoint2, aCint1->LastParameter());
  DStr.ChangeCurveInterferences(indCurve).Append(anInterfp1);
  DStr.ChangeCurveInterferences(indCurve).Append(anInterfp2);

  // interference of aCint1 on the SurfData number anInd

  anInterfc = ChFi3d_FilCurveInDS(indCurve, indSurf1, aPCint12, anOrSurf);
  DStr.ChangeSurfaceInterferences(indSurf1).Append(anInterfc);

  // interference of aCint1 on aFace

  anOrFace = aFace.Orientation();

  if (anOrFace == anOrSD1)
  {
    anOrFace = TopAbs::Reverse(anOrSurf);
  }
  else
  {
    anOrFace = anOrSurf;
  }

  anInterfc = ChFi3d_FilCurveInDS(indCurve, indaFace, aPCint11, anOrFace);
  DStr.ChangeShapeInterferences(indaFace).Append(anInterfc);
  // ---------------------------------------------------------------
  // storage of aCline passing through aCP1onArc and aCP2NotonArc

  occ::handle<Geom_Curve> aTrCline =
    new Geom_TrimmedCurve(aCline, aCline->FirstParameter(), aCline->LastParameter());
  double             aTolerance = DStr.Curve(indLine).Tolerance();
  TopOpeBRepDS_Curve aTct3(aTrCline, aTolerance);

  indCurve = DStr.AddCurve(aTct3);

  aTrCline->D0(aTrCline->FirstParameter(), aPext1);

  if (aPext1.Distance(aCP1onArc.Point()) < aTol3d)
  {
    indPoint1 = indCP1onArc;
    indPoint2 = indCP2NotonArc;
  }
  else
  {
    indPoint1 = indCP2NotonArc;
    indPoint2 = indCP1onArc;
  }
  //  definition of the orientation of aTrCline

  aTrCline->D0(aTrCline->FirstParameter(), aP1);
  aTrCline->D0(aTrCline->LastParameter(), aP2);
  aCint1->D0(aCint1->FirstParameter(), aP3);
  aCint1->D0(aCint1->LastParameter(), aP4);

  if (aP1.Distance(aP4) > aTol3d && aP2.Distance(aP3) > aTol3d)
  {
    anOrSurf = TopAbs::Reverse(anOrSurf);
  }

  anInterfp1 = ChFi3d_FilPointInDS(TopAbs_FORWARD, indCurve, indPoint1, aTrCline->FirstParameter());
  anInterfp2 = ChFi3d_FilPointInDS(TopAbs_REVERSED, indCurve, indPoint2, aTrCline->LastParameter());
  DStr.ChangeCurveInterferences(indCurve).Append(anInterfp1);
  DStr.ChangeCurveInterferences(indCurve).Append(anInterfp2);

  // interference of aTrCline on the SurfData number anInd

  anInterfc = ChFi3d_FilCurveInDS(indCurve, indSurf1, aPClineOnSurf, anOrSurf);
  DStr.ChangeSurfaceInterferences(indSurf1).Append(anInterfc);

  // interference de ctlin par rapport a Fvoisin
  indShape = DStr.AddShape(aNeighborFace);
  anOrFace = aNeighborFace.Orientation();

  if (anOrFace == anOrSD1)
  {
    anOrFace = TopAbs::Reverse(anOrSurf);
  }
  else
  {
    anOrFace = anOrSurf;
  }

  anInterfc = ChFi3d_FilCurveInDS(indCurve, indShape, aPClineOnFace, anOrFace);
  DStr.ChangeShapeInterferences(indShape).Append(anInterfc);
}

//  Modified by Sergey KHROMOV - Thu Apr 11 12:23:40 2002 End

//==============================================================
// function : FindFace
// purpose  : attention it works only if there is only one common face
//           between P1,P2,V
//===========================================================

bool ChFi3d_Builder::FindFace(const TopoDS_Vertex&      V,
                              const ChFiDS_CommonPoint& P1,
                              const ChFiDS_CommonPoint& P2,
                              TopoDS_Face&              Fv) const
{
  TopoDS_Face Favoid;
  return FindFace(V, P1, P2, Fv, Favoid);
}

bool ChFi3d_Builder::FindFace(const TopoDS_Vertex&      V,
                              const ChFiDS_CommonPoint& P1,
                              const ChFiDS_CommonPoint& P2,
                              TopoDS_Face&              Fv,
                              const TopoDS_Face&        Favoid) const
{
  if (P1.IsVertex() || P2.IsVertex())
  {
#ifdef OCCT_DEBUG
    std::cout << "change of face on vertex" << std::endl;
#endif
  }
  if (!(P1.IsOnArc() && P2.IsOnArc()))
  {
    return false;
  }
  NCollection_List<TopoDS_Shape>::Iterator It, Jt;
  bool                                     Found = false;
  for (It.Initialize(myEFMap(P1.Arc())); It.More() && !Found; It.Next())
  {
    Fv = TopoDS::Face(It.Value());
    if (!Fv.IsSame(Favoid))
    {
      for (Jt.Initialize(myEFMap(P2.Arc())); Jt.More() && !Found; Jt.Next())
      {
        if (TopoDS::Face(Jt.Value()).IsSame(Fv))
        {
          Found = true;
        }
      }
    }
  }
#ifdef OCCT_DEBUG
  bool ContainsV = false;
  if (Found)
  {
    for (It.Initialize(myVFMap(V)); It.More(); It.Next())
    {
      if (TopoDS::Face(It.Value()).IsSame(Fv))
      {
        ContainsV = true;
        break;
      }
    }
  }
  if (!ContainsV)
  {
    std::cout << "FindFace : the extremity of the spine is not in the end face" << std::endl;
  }
#else
  (void)V; // avoid compiler warning on unused variable
#endif
  return Found;
}

//=======================================================================
// function : MoreSurfdata
// purpose  : detects if the intersection at end concerns several Surfdata
//=======================================================================
bool ChFi3d_Builder::MoreSurfdata(const int Index) const
{
  // intersection at end is created on several surfdata if :
  // - the number of surfdata concerning the vertex is more than 1.
  // - and if the last but one surfdata has one of commonpoints on one of
  // two arcs, which constitute the intersections of the face at end and of the fillet

  NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator It;
  It.Initialize(myVDataMap(Index));
  occ::handle<ChFiDS_Stripe>&                         stripe = It.ChangeValue();
  NCollection_Sequence<occ::handle<ChFiDS_SurfData>>& SeqFil =
    stripe->ChangeSetOfSurfData()->ChangeSequence();
  const TopoDS_Vertex&          Vtx     = myVDataMap.FindKey(Index);
  int                           sens    = 0;
  int                           num     = ChFi3d_IndexOfSurfData(Vtx, stripe, sens);
  bool                          isfirst = (sens == 1);
  occ::handle<ChFiDS_SurfData>& Fd      = SeqFil.ChangeValue(num);
  ChFiDS_CommonPoint&           CV1     = Fd->ChangeVertex(isfirst, 1);
  ChFiDS_CommonPoint&           CV2     = Fd->ChangeVertex(isfirst, 2);

  int         num1, num2, nbsurf;
  TopoDS_Face Fv;
  bool        inters, oksurf;
  nbsurf = stripe->SetOfSurfData()->Length();
  // Fv is the face at end
  inters = FindFace(Vtx, CV1, CV2, Fv);
  if (sens == 1)
  {
    num1 = 1;
    num2 = num1 + 1;
  }
  else
  {
    num1 = nbsurf;
    num2 = num1 - 1;
  }

  oksurf = false;

  if (nbsurf != 1 && inters)
  {

    // determination of arc1 and arc2 intersection of the fillet and the face at end

    TopoDS_Edge                              arc1, arc2;
    NCollection_List<TopoDS_Shape>::Iterator ItE;
    bool                                     trouve = false;
    for (ItE.Initialize(myVEMap(Vtx)); ItE.More() && !trouve; ItE.Next())
    {
      arc1 = TopoDS::Edge(ItE.Value());
      if (containE(Fv, arc1))
      {
        trouve = true;
      }
    }
    trouve = false;
    for (ItE.Initialize(myVEMap(Vtx)); ItE.More() && !trouve; ItE.Next())
    {
      arc2 = TopoDS::Edge(ItE.Value());
      if (containE(Fv, arc2) && !arc2.IsSame(arc1))
      {
        trouve = true;
      }
    }

    occ::handle<ChFiDS_SurfData> Fd1 = SeqFil.ChangeValue(num2);
    ChFiDS_CommonPoint&          CV3 = Fd1->ChangeVertex(isfirst, 1);
    ChFiDS_CommonPoint&          CV4 = Fd1->ChangeVertex(isfirst, 2);

    if (CV3.IsOnArc())
    {
      if (CV3.Arc().IsSame(arc1))
      {
        if (CV1.Point().Distance(CV3.Point()) < 1.e-4)
        {
          oksurf = true;
        }
      }
      else if (CV3.Arc().IsSame(arc2))
      {
        if (CV2.Point().Distance(CV3.Point()) < 1.e-4)
        {
          oksurf = true;
        }
      }
    }

    if (CV4.IsOnArc())
    {
      if (CV1.Point().Distance(CV4.Point()) < 1.e-4)
      {
        oksurf = true;
      }
      else if (CV4.Arc().IsSame(arc2))
      {
        if (CV2.Point().Distance(CV4.Point()) < 1.e-4)
        {
          oksurf = true;
        }
      }
    }
  }
  return oksurf;
}

// Case of fillets on top with 4 edges, one of them is on the same geometry as the edgeof the fillet

void ChFi3d_Builder::IntersectMoreCorner(const int Index)
{
  TopOpeBRepDS_DataStructure& DStr = myDS->ChangeDS();

#ifdef OCCT_DEBUG
  OSD_Chronometer ch; // init perf pour PerformSetOfKPart
#endif
  // The fillet is returned,
  NCollection_List<occ::handle<ChFiDS_Stripe>>::Iterator StrIt;
  StrIt.Initialize(myVDataMap(Index));
  occ::handle<ChFiDS_Stripe>                          stripe = StrIt.Value();
  const occ::handle<ChFiDS_Spine>                     spine  = stripe->Spine();
  NCollection_Sequence<occ::handle<ChFiDS_SurfData>>& SeqFil =
    stripe->ChangeSetOfSurfData()->ChangeSequence();
  // the top,
  const TopoDS_Vertex& Vtx = myVDataMap.FindKey(Index);
  // the SurfData concerned and its CommonPoints,
  int sens = 0;

  // Choose the proper SurfData
  int  num     = ChFi3d_IndexOfSurfData(Vtx, stripe, sens);
  bool isfirst = (sens == 1);
  if (isfirst)
  {
    for (; num < SeqFil.Length()
           && ((SeqFil.Value(num)->IndexOfS1() == 0) || (SeqFil.Value(num)->IndexOfS2() == 0));)
    {
      SeqFil.Remove(num); // The surplus is removed
    }
  }
  else
  {
    for (; num > 1
           && ((SeqFil.Value(num)->IndexOfS1() == 0) || (SeqFil.Value(num)->IndexOfS2() == 0));)
    {
      SeqFil.Remove(num); // The surplus is removed
      num--;
    }
  }

  occ::handle<ChFiDS_SurfData>& Fd  = SeqFil.ChangeValue(num);
  ChFiDS_CommonPoint&           CV1 = Fd->ChangeVertex(isfirst, 1);
  ChFiDS_CommonPoint&           CV2 = Fd->ChangeVertex(isfirst, 2);
  // To evaluate the cloud of new points.
  Bnd_Box box1, box2;

  // The cases of cap are processed separately from intersection.
  // ----------------------------------------------------------

  TopoDS_Face Fv, Fad, Fop, Fopbis;
  TopoDS_Edge Arcpiv, Arcprol, Arcspine, Arcprolbis;
  if (isfirst)
  {
    Arcspine = spine->Edges(1);
  }
  else
  {
    Arcspine = spine->Edges(spine->NbEdges());
  }
  TopAbs_Orientation               OArcprolbis = TopAbs_FORWARD;
  TopAbs_Orientation               OArcprolv = TopAbs_FORWARD, OArcprolop = TopAbs_FORWARD;
  int                              ICurve;
  occ::handle<BRepAdaptor_Surface> HBs  = new BRepAdaptor_Surface();
  occ::handle<BRepAdaptor_Surface> HBad = new BRepAdaptor_Surface();
  occ::handle<BRepAdaptor_Surface> HBop = new BRepAdaptor_Surface();
  BRepAdaptor_Surface&             Bs   = *HBs;
  BRepAdaptor_Surface&             Bad  = *HBad;
  BRepAdaptor_Surface&             Bop  = *HBop;
  occ::handle<Geom_Curve>          Cc;
  occ::handle<Geom2d_Curve>        Pc, Ps;
  double                           Ubid, Vbid; //,mu,Mu,mv,Mv;
  double                           Udeb = 0., Ufin = 0.;
  // gp_Pnt2d UVf1,UVl1,UVf2,UVl2;
  // double Du,Dv,Step;
  bool inters  = true;
  int  IFadArc = 1, IFopArc = 2;
  Fop = TopoDS::Face(DStr.Shape(Fd->Index(IFopArc)));
  TopExp_Explorer ex;

#ifdef OCCT_DEBUG
  ChFi3d_InitChron(ch); // init perf condition
#endif
  {
    if (!CV1.IsOnArc() && !CV2.IsOnArc())
    {
      throw Standard_Failure("Corner intersmore : no point on arc");
    }
    else if (CV1.IsOnArc() && CV2.IsOnArc())
    {
      bool sur2 = false;
      for (ex.Init(CV1.Arc(), TopAbs_VERTEX); ex.More(); ex.Next())
      {
        if (Vtx.IsSame(ex.Current()))
        {
          break;
        }
      }
      for (ex.Init(CV2.Arc(), TopAbs_VERTEX); ex.More(); ex.Next())
      {
        if (Vtx.IsSame(ex.Current()))
        {
          sur2 = true;
          break;
        }
      }
      if (sur2)
      {
        IFadArc = 2;
      }
    }
    else if (CV2.IsOnArc())
    {
      IFadArc = 2;
    }
    IFopArc = 3 - IFadArc;

    Arcpiv = Fd->Vertex(isfirst, IFadArc).Arc();
    Fad    = TopoDS::Face(DStr.Shape(Fd->Index(IFadArc)));
    Fop    = TopoDS::Face(DStr.Shape(Fd->Index(IFopArc)));
    NCollection_List<TopoDS_Shape>::Iterator It;
    // The face at end is returned without control of its unicity.
    for (It.Initialize(myEFMap(Arcpiv)); It.More(); It.Next())
    {
      if (!Fad.IsSame(It.Value()))
      {
        Fv = TopoDS::Face(It.Value());
        break;
      }
    }

    // does the face at end contain the Vertex ?
    bool isinface = false;
    for (ex.Init(Fv, TopAbs_VERTEX); ex.More(); ex.Next())
    {
      if (ex.Current().IsSame(Vtx))
      {
        isinface = true;
        break;
      }
    }
    if (!isinface)
    {
      IFadArc = 3 - IFadArc;
      IFopArc = 3 - IFopArc;
      Arcpiv  = Fd->Vertex(isfirst, IFadArc).Arc();
      Fad     = TopoDS::Face(DStr.Shape(Fd->Index(IFadArc)));
      Fop     = TopoDS::Face(DStr.Shape(Fd->Index(IFopArc)));
      // NCollection_List<TopoDS_Shape>::Iterator It;
      // The face at end is returned without control of its unicity.
      for (It.Initialize(myEFMap(Arcpiv)); It.More(); It.Next())
      {
        if (!Fad.IsSame(It.Value()))
        {
          Fv = TopoDS::Face(It.Value());
          break;
        }
      }
    }

    if (Fv.IsNull())
    {
      throw StdFail_NotDone("OneCorner : face at end is not found");
    }

    Fv.Orientation(TopAbs_FORWARD);
    Fad.Orientation(TopAbs_FORWARD);

    // In the same way the edge to be extended is returned.
    for (It.Initialize(myVEMap(Vtx)); It.More() && Arcprol.IsNull(); It.Next())
    {
      if (!Arcpiv.IsSame(It.Value()))
      {
        for (ex.Init(Fv, TopAbs_EDGE); ex.More(); ex.Next())
        {
          if (It.Value().IsSame(ex.Current()))
          {
            Arcprol   = TopoDS::Edge(It.Value());
            OArcprolv = ex.Current().Orientation();
            break;
          }
        }
      }
    }

    // Guard: Arcprol may be null if the loop above found no matching edge.
    if (Arcprol.IsNull())
    {
      throw StdFail_NotDone("IntersectMoreCorner: edge to be extended is not found");
    }

    // Fopbis is the face containing the trace of fillet CP.Arc() which of does not contain Vtx.
    // Normally Fobis is either the same as Fop (cylinder), or Fobis is G1 with Fop.
    Fopbis.Orientation(TopAbs_FORWARD);

    // Fop calls the 4th face non-used for the vertex
    cherche_face(myVFMap(Vtx), Arcprol, Fad, Fv, Fv, Fopbis);
    Fop.Orientation(TopAbs_FORWARD);
    for (ex.Init(Fopbis, TopAbs_EDGE); ex.More(); ex.Next())
    {
      if (Arcprol.IsSame(ex.Current()))
      {
        OArcprolop = ex.Current().Orientation();
        break;
      }
    }
    TopoDS_Face               FFv;
    double                    tol;
    int                       prol = 0;
    BRep_Builder              BRE;
    occ::handle<Geom_Surface> Sface;
    Sface = BRep_Tool::Surface(Fv);
    ChFi3d_ExtendSurface(Sface, prol);
    tol = BRep_Tool::Tolerance(Fv);
    BRE.MakeFace(FFv, Sface, tol);
    if (prol)
    {
      Bs.Initialize(FFv, false);
      DStr.SetNewSurface(Fv, Sface);
    }
    else
    {
      Bs.Initialize(Fv, false);
    }
    Bad.Initialize(Fad);
    Bop.Initialize(Fop);
  }
  // it is necessary to modify the CommonPoint
  // in the space and its parameter in FaceInterference.
  // So both of them are returned in references
  // non const. Attention the modifications are done behind
  // CV1,CV2,Fi1,Fi2.
  ChFiDS_CommonPoint&      CPopArc = Fd->ChangeVertex(isfirst, IFopArc);
  ChFiDS_FaceInterference& FiopArc = Fd->ChangeInterference(IFopArc);
  ChFiDS_CommonPoint&      CPadArc = Fd->ChangeVertex(isfirst, IFadArc);
  ChFiDS_FaceInterference& FiadArc = Fd->ChangeInterference(IFadArc);
  // the parameter of the vertex is initialized with the value
  // of its opposing vertex (point on arc).
  double                           wop = Fd->ChangeInterference(IFadArc).Parameter(isfirst);
  occ::handle<Geom_Curve>          c3df;
  occ::handle<GeomAdaptor_Surface> HGs =
    new GeomAdaptor_Surface(DStr.Surface(Fd->Surf()).Surface());
  gp_Pnt2d p2dbout;
  {

    // add here more or less restrictive criteria to
    // decide if the intersection with face is done at the
    // extended end or if there will be a cap on sharp end.
    c3df = DStr.Curve(FiopArc.LineIndex()).Curve();
    if (c3df.IsNull())
    {
      throw Standard_ConstructionError("IntersectMoreCorner : no fillet line on the face");
    }
    double                         uf = FiopArc.FirstParameter();
    double                         ul = FiopArc.LastParameter();
    occ::handle<GeomAdaptor_Curve> Hc3df;
    if (c3df->IsPeriodic())
    {
      Hc3df = new GeomAdaptor_Curve(c3df);
    }
    else
    {
      Hc3df = new GeomAdaptor_Curve(c3df, uf, ul);
    }
    inters = Update(HBs, Hc3df, FiopArc, CPopArc, p2dbout, isfirst, wop);
    //  Modified by Sergey KHROMOV - Fri Dec 21 18:08:27 2001 Begin
    //  if(!inters && BRep_Tool::Continuity(Arcprol,Fv,Fop) != GeomAbs_C0){
    if (!inters && ChFi3d::IsTangentFaces(Arcprol, Fv, Fop))
    {
      //  Modified by Sergey KHROMOV - Fri Dec 21 18:08:29 2001 End
      // Arcprol is an edge of tangency, ultimate adjustment by an extrema curve/curve is attempted.
      double                    ff, ll;
      occ::handle<Geom2d_Curve> gpcprol = PCurveInFace(Arcprol, Fv, ff, ll);
      if (gpcprol.IsNull())
      {
        throw Standard_ConstructionError("Failed to get p-curve of edge");
      }
      occ::handle<Geom2dAdaptor_Curve> pcprol  = new Geom2dAdaptor_Curve(gpcprol);
      double                           partemp = BRep_Tool::Parameter(Vtx, Arcprol);
      inters =
        Update(HBs, pcprol, HGs, FiopArc, CPopArc, p2dbout, isfirst, partemp, wop, 10 * tolapp3d);
    }
    occ::handle<BRepAdaptor_Curve2d> pced = new BRepAdaptor_Curve2d();
    pced->Initialize(CPadArc.Arc(), Fv);
    Update(HBs, pced, HGs, FiadArc, CPadArc, isfirst);
  }
#ifdef OCCT_DEBUG
  ChFi3d_ResultChron(ch, t_same); // result perf condition if (same)
  ChFi3d_InitChron(ch);           // init perf condition if (inters)
#endif

  TopoDS_Edge                    edgecouture;
  bool                           couture, intcouture = false;
  double                         tolreached = tolapp3d;
  double                         par1 = 0., par2 = 0.;
  int                            indpt = 0, Icurv1 = 0, Icurv2 = 0;
  occ::handle<Geom_TrimmedCurve> curv1, curv2;
  occ::handle<Geom2d_Curve>      c2d1, c2d2;

  int Isurf = Fd->Surf();

  if (inters)
  {
    HGs                                = ChFi3d_BoundSurf(DStr, Fd, 1, 2);
    const ChFiDS_FaceInterference& Fi1 = Fd->InterferenceOnS1();
    const ChFiDS_FaceInterference& Fi2 = Fd->InterferenceOnS2();
    NCollection_Array1<double>     Pardeb(1, 4), Parfin(1, 4);
    gp_Pnt2d                       pfil1, pfac1, pfil2, pfac2;
    occ::handle<Geom2d_Curve>      Hc1, Hc2;
    if (IFopArc == 1)
    {
      pfac1 = p2dbout;
    }
    else
    {
      Hc1 = PCurveInFace(CV1.Arc(), Fv, Ubid, Ubid);
      if (Hc1.IsNull())
      {
        throw Standard_ConstructionError("Failed to get p-curve of edge");
      }
      pfac1 = Hc1->Value(CV1.ParameterOnArc());
    }
    if (IFopArc == 2)
    {
      pfac2 = p2dbout;
    }
    else
    {
      Hc2 = PCurveInFace(CV2.Arc(), Fv, Ubid, Ubid);
      if (Hc2.IsNull())
      {
        throw Standard_ConstructionError("Failed to get p-curve of edge");
      }
      pfac2 = Hc2->Value(CV2.ParameterOnArc());
    }
    if (Fi1.LineIndex() != 0)
    {
      pfil1 = Fi1.PCurveOnSurf()->Value(Fi1.Parameter(isfirst));
    }
    else
    {
      pfil1 = Fi1.PCurveOnSurf()->Value(Fi1.Parameter(!isfirst));
    }
    if (Fi2.LineIndex() != 0)
    {
      pfil2 = Fi2.PCurveOnSurf()->Value(Fi2.Parameter(isfirst));
    }
    else
    {
      pfil2 = Fi2.PCurveOnSurf()->Value(Fi2.Parameter(!isfirst));
    }
    ChFi3d_Recale(Bs, pfac1, pfac2, (IFadArc == 1));
    Pardeb(1) = pfil1.X();
    Pardeb(2) = pfil1.Y();
    Pardeb(3) = pfac1.X();
    Pardeb(4) = pfac1.Y();
    Parfin(1) = pfil2.X();
    Parfin(2) = pfil2.Y();
    Parfin(3) = pfac2.X();
    Parfin(4) = pfac2.Y();
    double uu1, uu2, vv1, vv2;
    ChFi3d_Boite(pfac1, pfac2, uu1, uu2, vv1, vv2);
    ChFi3d_BoundFac(Bs, uu1, uu2, vv1, vv2);

    if (!ChFi3d_ComputeCurves(HGs, HBs, Pardeb, Parfin, Cc, Ps, Pc, tolapp3d, tol2d, tolreached))
    {
      throw Standard_Failure("OneCorner : failed calculation intersection");
    }

    Udeb = Cc->FirstParameter();
    Ufin = Cc->LastParameter();

    // check if the curve has an intersection with sewing edge

    ChFi3d_Couture(Fv, couture, edgecouture);

    if (couture && !BRep_Tool::Degenerated(edgecouture))
    {

      // double Ubid,Vbid;
      occ::handle<Geom_Curve>        C     = BRep_Tool::Curve(edgecouture, Ubid, Vbid);
      occ::handle<Geom_TrimmedCurve> Ctrim = new Geom_TrimmedCurve(C, Ubid, Vbid);
      GeomAdaptor_Curve              cur1(Ctrim->BasisCurve());
      GeomAdaptor_Curve              cur2(Cc);
      Extrema_ExtCC                  extCC(cur1, cur2);
      if (extCC.IsDone() && extCC.NbExt() != 0)
      {
        int    imin     = 0;
        double dist2min = RealLast();
        for (int i = 1; i <= extCC.NbExt(); i++)
        {
          if (extCC.SquareDistance(i) < dist2min)
          {
            dist2min = extCC.SquareDistance(i);
            imin     = i;
          }
        }
        if (dist2min <= Precision::SquareConfusion())
        {
          Extrema_POnCurv ponc1, ponc2;
          extCC.Points(imin, ponc1, ponc2);
          par1       = ponc1.Parameter();
          par2       = ponc2.Parameter();
          double Tol = 1.e-4;
          if (std::abs(par2 - Udeb) > Tol && std::abs(Ufin - par2) > Tol)
          {
            gp_Pnt             P1 = ponc1.Value();
            TopOpeBRepDS_Point tpoint(P1, Tol);
            indpt      = DStr.AddPoint(tpoint);
            intcouture = true;
            curv1      = new Geom_TrimmedCurve(Cc, Udeb, par2);
            curv2      = new Geom_TrimmedCurve(Cc, par2, Ufin);
            TopOpeBRepDS_Curve tcurv1(curv1, tolreached);
            TopOpeBRepDS_Curve tcurv2(curv2, tolreached);
            Icurv1 = DStr.AddCurve(tcurv1);
            Icurv2 = DStr.AddCurve(tcurv2);
          }
        }
      }
    }
  }

  else
  {
    throw Standard_NotImplemented("OneCorner : cap not written");
  }
  int                IShape = DStr.AddShape(Fv);
  TopAbs_Orientation Et     = TopAbs_FORWARD;
  if (IFadArc == 1)
  {
    TopExp_Explorer Exp;
    for (Exp.Init(Fv.Oriented(TopAbs_FORWARD), TopAbs_EDGE); Exp.More(); Exp.Next())
    {
      if (Exp.Current().IsSame(CV1.Arc()))
      {
        Et = TopAbs::Reverse(TopAbs::Compose(Exp.Current().Orientation(), CV1.TransitionOnArc()));
        break;
      }
    }
  }
  else
  {
    TopExp_Explorer Exp;
    for (Exp.Init(Fv.Oriented(TopAbs_FORWARD), TopAbs_EDGE); Exp.More(); Exp.Next())
    {
      if (Exp.Current().IsSame(CV2.Arc()))
      {
        Et = TopAbs::Compose(Exp.Current().Orientation(), CV2.TransitionOnArc());
        break;
      }
    }

    //
  }

#ifdef OCCT_DEBUG
  ChFi3d_ResultChron(ch, t_inter); // result perf condition if (inter)
  ChFi3d_InitChron(ch);            // init perf condition  if ( inters)
#endif

  stripe->SetIndexPoint(ChFi3d_IndexPointInDS(CV1, DStr), isfirst, 1);
  stripe->SetIndexPoint(ChFi3d_IndexPointInDS(CV2, DStr), isfirst, 2);

  if (!intcouture)
  {
    // there is no intersection with edge of sewing
    // curve Cc is stored in the stripe
    // the storage in the DS is done by FILDS.

    TopOpeBRepDS_Curve Tc(Cc, tolreached);
    ICurve = DStr.AddCurve(Tc);
    occ::handle<TopOpeBRepDS_SurfaceCurveInterference> Interfc =
      ChFi3d_FilCurveInDS(ICurve, IShape, Pc, Et);
    DStr.ChangeShapeInterferences(IShape).Append(Interfc);
    stripe->ChangePCurve(isfirst) = Ps;
    stripe->SetCurve(ICurve, isfirst);
    stripe->SetParameters(isfirst, Udeb, Ufin);
  }
  else
  {
    // curves curv1 and curv2 are stored in the DS
    // these curves are not reconstructed by FILDS as
    // stripe->InDS(isfirst) is placed;

    // interferences of curv1 and curv2 on Fv
    ComputeCurve2d(curv1, Fv, c2d1);
    occ::handle<TopOpeBRepDS_SurfaceCurveInterference> InterFv;
    InterFv = ChFi3d_FilCurveInDS(Icurv1, IShape, c2d1, Et);
    DStr.ChangeShapeInterferences(IShape).Append(InterFv);
    ComputeCurve2d(curv2, Fv, c2d2);
    InterFv = ChFi3d_FilCurveInDS(Icurv2, IShape, c2d2, Et);
    DStr.ChangeShapeInterferences(IShape).Append(InterFv);
    // interferences of curv1 and curv2 on Isurf
    if (Fd->Orientation() == Fv.Orientation())
    {
      Et = TopAbs::Reverse(Et);
    }
    c2d1    = new Geom2d_TrimmedCurve(Ps, Udeb, par2);
    InterFv = ChFi3d_FilCurveInDS(Icurv1, Isurf, c2d1, Et);
    DStr.ChangeSurfaceInterferences(Isurf).Append(InterFv);
    c2d2    = new Geom2d_TrimmedCurve(Ps, par2, Ufin);
    InterFv = ChFi3d_FilCurveInDS(Icurv2, Isurf, c2d2, Et);
    DStr.ChangeSurfaceInterferences(Isurf).Append(InterFv);

    // limitation of the sewing edge
    int                                              Iarc = DStr.AddShape(edgecouture);
    occ::handle<TopOpeBRepDS_CurvePointInterference> Interfedge;
    TopAbs_Orientation                               ori;
    TopoDS_Vertex                                    Vdeb, Vfin;
    Vdeb = TopExp::FirstVertex(edgecouture);
    Vfin = TopExp::LastVertex(edgecouture);
    double pard, parf;
    pard = BRep_Tool::Parameter(Vdeb, edgecouture);
    parf = BRep_Tool::Parameter(Vfin, edgecouture);
    if (std::abs(par1 - pard) < std::abs(parf - par1))
    {
      ori = TopAbs_FORWARD;
    }
    else
    {
      ori = TopAbs_REVERSED;
    }
    Interfedge = ChFi3d_FilPointInDS(ori, Iarc, indpt, par1);
    DStr.ChangeShapeInterferences(Iarc).Append(Interfedge);

    // creation of CurveInterferences from Icurv1 and Icurv2
    stripe->InDS(isfirst);
    int                                              ind1 = stripe->IndexPoint(isfirst, 1);
    int                                              ind2 = stripe->IndexPoint(isfirst, 2);
    occ::handle<TopOpeBRepDS_CurvePointInterference> interfprol =
      ChFi3d_FilPointInDS(TopAbs_FORWARD, Icurv1, ind1, Udeb);
    DStr.ChangeCurveInterferences(Icurv1).Append(interfprol);
    interfprol = ChFi3d_FilPointInDS(TopAbs_REVERSED, Icurv1, indpt, par2);
    DStr.ChangeCurveInterferences(Icurv1).Append(interfprol);
    interfprol = ChFi3d_FilPointInDS(TopAbs_FORWARD, Icurv2, indpt, par2);
    DStr.ChangeCurveInterferences(Icurv2).Append(interfprol);
    interfprol = ChFi3d_FilPointInDS(TopAbs_REVERSED, Icurv2, ind2, Ufin);
    DStr.ChangeCurveInterferences(Icurv2).Append(interfprol);
  }

  ChFi3d_EnlargeBox(HBs, Pc, Udeb, Ufin, box1, box2);

  if (inters)
  {
    //

    // The small end of curve missing for the extension
    // of the face at end and the limitation of the opposing face is added.

    // Above all the points cut the points with the edge of the spine.
    int                IArcspine = DStr.AddShape(Arcspine);
    int                IVtx      = DStr.AddShape(Vtx);
    TopAbs_Orientation OVtx2     = TopAbs_FORWARD;
    TopAbs_Orientation OVtx      = TopAbs_FORWARD;
    for (ex.Init(Arcspine.Oriented(TopAbs_FORWARD), TopAbs_VERTEX); ex.More(); ex.Next())
    {
      if (Vtx.IsSame(ex.Current()))
      {
        OVtx = ex.Current().Orientation();
        break;
      }
    }
    OVtx                                                    = TopAbs::Reverse(OVtx);
    double                                           parVtx = BRep_Tool::Parameter(Vtx, Arcspine);
    occ::handle<TopOpeBRepDS_CurvePointInterference> interfv =
      ChFi3d_FilVertexInDS(OVtx, IArcspine, IVtx, parVtx);
    DStr.ChangeShapeInterferences(IArcspine).Append(interfv);

    // Modif of lvt to find the suite of Arcprol in the other face
    {
      NCollection_List<TopoDS_Shape>::Iterator It;
      for (It.Initialize(myVEMap(Vtx)); It.More(); It.Next())
      {
        if (!(Arcprol.IsSame(It.Value()) || Arcspine.IsSame(It.Value())
              || Arcpiv.IsSame(It.Value())))
        {
          Arcprolbis = TopoDS::Edge(It.Value());
          break;
        }
      }
    }
    // end of modif

    // Guard: Arcprolbis may be null if myVEMap(Vtx) contained only
    // Arcprol, Arcspine, and Arcpiv (e.g. when topology maps carry stale
    // data from shared TShape pointers after prior operations).
    if (Arcprolbis.IsNull())
    {
      throw Standard_ConstructionError(
        "IntersectMoreCorner: continuation edge (Arcprolbis) is not found");
    }

    // Now the missing curves are constructed.
    for (ex.Init(Arcprolbis.Oriented(TopAbs_FORWARD), TopAbs_VERTEX); ex.More(); ex.Next())
    {
      if (Vtx.IsSame(ex.Current()))
      {
        OVtx2 = ex.Current().Orientation();
        break;
      }
    }
    for (ex.Init(Arcprol.Oriented(TopAbs_FORWARD), TopAbs_VERTEX); ex.More(); ex.Next())
    {
      if (Vtx.IsSame(ex.Current()))
      {
        OVtx = ex.Current().Orientation();
        break;
      }
    }
    // it is checked if Fop has a sewing edge

    //     TopoDS_Edge edgecouture;
    //     bool couture;
    ChFi3d_Couture(Fop, couture, edgecouture);
    occ::handle<Geom2d_Curve> Hc;
    //    parVtx = BRep_Tool::Parameter(Vtx,Arcprol);
    const ChFiDS_FaceInterference& Fiop = Fd->Interference(IFopArc);
    gp_Pnt2d                       pop1, pop2, pv1, pv2;
    // deb modif
    parVtx = BRep_Tool::Parameter(Vtx, Arcprolbis);
    //  Modified by skv - Thu Aug 21 11:55:58 2008 OCC20222 Begin
    //    if(Fop.IsSame(Fopbis)) OArcprolbis = OArcprolop;
    //    else OArcprolbis = Arcprolbis.Orientation();
    if (Fop.IsSame(Fopbis))
    {
      OArcprolbis = OArcprolop;
    }
    else
    {
      for (ex.Init(Fop, TopAbs_EDGE); ex.More(); ex.Next())
      {
        if (Arcprolbis.IsSame(ex.Current()))
        {
          OArcprolbis = ex.Current().Orientation();
          break;
        }
      }
    }
    //  Modified by skv - Thu Aug 21 11:55:58 2008 OCC20222 End
    // fin modif
    Hc = PCurveInFace(Arcprolbis, Fop, Ubid, Ubid);
    if (Hc.IsNull())
    {
      throw Standard_ConstructionError("Failed to get p-curve of edge");
    }
    pop1 = Hc->Value(parVtx);
    pop2 = Fiop.PCurveOnFace()->Value(Fiop.Parameter(isfirst));
    Hc   = PCurveInFace(Arcprol, Fv, Ubid, Ubid);
    if (Hc.IsNull())
    {
      throw Standard_ConstructionError("Failed to get p-curve of edge");
    }
    // modif
    parVtx = BRep_Tool::Parameter(Vtx, Arcprol);
    // fin modif
    pv1 = Hc->Value(parVtx);
    pv2 = p2dbout;
    ChFi3d_Recale(Bs, pv1, pv2, true);
    NCollection_Array1<double> Pardeb(1, 4), Parfin(1, 4);
    Pardeb(1) = pop1.X();
    Pardeb(2) = pop1.Y();
    Pardeb(3) = pv1.X();
    Pardeb(4) = pv1.Y();
    Parfin(1) = pop2.X();
    Parfin(2) = pop2.Y();
    Parfin(3) = pv2.X();
    Parfin(4) = pv2.Y();
    double uu1, uu2, vv1, vv2;
    ChFi3d_Boite(pv1, pv2, uu1, uu2, vv1, vv2);
    ChFi3d_BoundFac(Bs, uu1, uu2, vv1, vv2);
    ChFi3d_Boite(pop1, pop2, uu1, uu2, vv1, vv2);
    ChFi3d_BoundFac(Bop, uu1, uu2, vv1, vv2);

    occ::handle<Geom_Curve>   zob3d;
    occ::handle<Geom2d_Curve> zob2dop, zob2dv;
    //    double tolreached;
    if (!ChFi3d_ComputeCurves(HBop,
                              HBs,
                              Pardeb,
                              Parfin,
                              zob3d,
                              zob2dop,
                              zob2dv,
                              tolapp3d,
                              tol2d,
                              tolreached))
    {
      throw Standard_Failure("OneCorner : echec calcul intersection");
    }

    Udeb = zob3d->FirstParameter();
    Ufin = zob3d->LastParameter();
    TopOpeBRepDS_Curve Zob(zob3d, tolreached);
    int                IZob = DStr.AddCurve(Zob);

    // it is not determined if the curve has an intersection with the sewing edge

    {
      Et = TopAbs::Reverse(TopAbs::Compose(OVtx, OArcprolv));
      int                                                Iop = DStr.AddShape(Fop);
      occ::handle<TopOpeBRepDS_SurfaceCurveInterference> InterFv =
        ChFi3d_FilCurveInDS(IZob, IShape, zob2dv, Et);
      DStr.ChangeShapeInterferences(IShape).Append(InterFv);
      // OVtx = TopAbs::Reverse(OVtx);
      //  Modified by skv - Thu Aug 21 11:55:58 2008 OCC20222 Begin
      //      Et = TopAbs::Reverse(TopAbs::Compose(OVtx,OArcprolbis));
      Et = TopAbs::Reverse(TopAbs::Compose(OVtx2, OArcprolbis));
      //  Modified by skv - Thu Aug 21 11:55:58 2008 OCC20222 End
      // OVtx = TopAbs::Reverse(OVtx);
      //      Et = TopAbs::Reverse(Et);
      occ::handle<TopOpeBRepDS_SurfaceCurveInterference> Interfop =
        ChFi3d_FilCurveInDS(IZob, Iop, zob2dop, Et);
      DStr.ChangeShapeInterferences(Iop).Append(Interfop);
      occ::handle<TopOpeBRepDS_CurvePointInterference> interfprol =
        ChFi3d_FilVertexInDS(TopAbs_FORWARD, IZob, IVtx, Udeb);
      DStr.ChangeCurveInterferences(IZob).Append(interfprol);
      int icc    = stripe->IndexPoint(isfirst, IFopArc);
      interfprol = ChFi3d_FilPointInDS(TopAbs_REVERSED, IZob, icc, Ufin);
      DStr.ChangeCurveInterferences(IZob).Append(interfprol);
    }
  }
  ChFi3d_EnlargeBox(DStr, stripe, Fd, box1, box2, isfirst);
  if (CV1.IsOnArc())
  {
    ChFi3d_EnlargeBox(CV1.Arc(), myEFMap(CV1.Arc()), CV1.ParameterOnArc(), box1);
  }
  if (CV2.IsOnArc())
  {
    ChFi3d_EnlargeBox(CV2.Arc(), myEFMap(CV2.Arc()), CV2.ParameterOnArc(), box2);
  }
  if (!CV1.IsVertex())
  {
    ChFi3d_SetPointTolerance(DStr, box1, stripe->IndexPoint(isfirst, 1));
  }
  if (!CV2.IsVertex())
  {
    ChFi3d_SetPointTolerance(DStr, box2, stripe->IndexPoint(isfirst, 2));
  }

#ifdef OCCT_DEBUG
  ChFi3d_ResultChron(ch, t_sameinter); // result perf condition if (same &&inter)
#endif
}
