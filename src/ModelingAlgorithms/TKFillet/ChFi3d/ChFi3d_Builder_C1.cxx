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
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepExtrema_ExtCC.hxx>
#include <BRepLib_MakeEdge.hxx>
#include <BRepTools.hxx>
#include <BRepTopAdaptor_FClass2d.hxx>
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
#include <ElSLib.hxx>
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
#include <Geom2dAPI_InterCurveCurve.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <Geom2dInt_GInter.hxx>
#include <Geom_BezierSurface.hxx>
#include <Geom_Conic.hxx>
#include <Geom_Line.hxx>
#include <Geom_Plane.hxx>
#include <GeomAPI_IntCS.hxx>
#include <GeomAPI_Interpolate.hxx>
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
    // the end is past the end -- but for one the last digits past it: the
    // line of a wall kept in pieces ends on a split that lies on the
    // surface (the face at end goes on up the split), and its end is the
    // crossing.
    if (isol == 0 && !isperiodic)
    {
      const double tolw = Precision::PConfusion();
      const double tolr = 1.e-12 * std::max(1., std::max(std::abs(uf), std::abs(ul)));
      for (i = 1; i <= nbp; i++)
      {
        const IntCurveSurface_IntersectionPoint& ip = Intersection.Point(i);
        w                                           = ip.W();
        if (uf - tolw > w || ul + tolw < w || std::abs(w - wop) >= dist)
        {
          continue;
        }
        if (uf - tolr <= w && ul + tolr >= w)
        {
          isol = i;
          dist = std::abs(w - wop);
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
  // at most once round the face
  int aNbEdges = 0;
  for (TopExp_Explorer ex(Fv, TopAbs_EDGE); ex.More(); ex.Next())
  {
    aNbEdges++;
  }
  for (int aStep = 0; aStep < aNbEdges; aStep++)
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
  // on the way there, on Fop: one crossing, an edge into a coplanar piece.
  // From a little before the end: the split can cross the line within the
  // tolerance of its end, on either side of it, where the walk ended at the
  // split without seeing it -- the line is straight, cut there all the same.
  NCollection_List<LineCrossing> aOnOp;
  if (!CrossingsOnFace(L, W0 - 2. * aSense * Tol, WP, Fop, Tol, aOnOp) || aOnOp.Extent() != 1)
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
  if (!ChFi3d_CoplanarPieces(Fop, Fn, Tol))
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
                               const occ::handle<Geom_Curve>&          L,
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
// function : ChFi3d_StoreLineOverSplit
// purpose  : StoreLineOverSplit, for the corners of two fillets.
//=======================================================================

void ChFi3d_StoreLineOverSplit(TopOpeBRepDS_DataStructure&             DStr,
                               const occ::handle<ChFiDS_Stripe>&       Stripe,
                               const occ::handle<ChFiDS_SurfData>&     Fd,
                               const bool                              isFirst,
                               const int                               Ons,
                               const int                               Isurf,
                               const occ::handle<GeomAdaptor_Surface>& HGs,
                               const occ::handle<Geom_Curve>&          L,
                               const TopoDS_Edge&                      Esplit,
                               const TopoDS_Face&                      Fn,
                               const double                            WQ,
                               const double                            WP,
                               const double                            ParQ,
                               const int                               IndP,
                               const double                            Tol,
                               NCollection_List<ChFiDS_Regul>&         Regul)
{
  StoreLineOverSplit(DStr, Stripe, Fd, isFirst, Ons, Isurf, HGs, L, Esplit, Fn, WQ, WP, ParQ, IndP,
                     Tol, Regul);
}

//=======================================================================
// function : ChFi3d_TransitionOnEdge
// purpose  : TransitionOnEdge, for the corners of two fillets.
//=======================================================================

TopAbs_Orientation ChFi3d_TransitionOnEdge(const ChFiDS_FaceInterference& Fi,
                                          const TopoDS_Face&             F,
                                          const TopoDS_Edge&             E,
                                          const bool                     isFirst)
{
  return TransitionOnEdge(Fi, F, E, isFirst);
}

//=======================================================================
// function : ChFi3d_PivotInTwoPieces
// purpose  : <Enear>, an edge of Vtx, and <Efar>, not one, are one edge in
//           two pieces: they meet at a vertex, the other end of <Enear>,
//           tangent to one another there and going on the same way. With
//           <Esplit> given, that vertex is an end of it as well -- the foot
//           of a split between two coplanar pieces of a face the edge
//           bounds. The pivot of a corner of two fillets, one of them wider
//           than the piece beside the vertex.
//=======================================================================

bool ChFi3d_PivotInTwoPieces(const TopoDS_Vertex& Vtx,
                             const TopoDS_Edge&   Enear,
                             const TopoDS_Edge&   Efar,
                             const TopoDS_Edge&   Esplit)
{
  if (Enear.IsSame(Efar) || BRep_Tool::Degenerated(Enear) || BRep_Tool::Degenerated(Efar))
  {
    return false;
  }
  TopoDS_Vertex aVn1, aVn2, aVf1, aVf2;
  TopExp::Vertices(Enear, aVn1, aVn2);
  TopExp::Vertices(Efar, aVf1, aVf2);
  if (aVn1.IsSame(aVn2) || aVf1.IsSame(aVf2) || aVf1.IsSame(Vtx) || aVf2.IsSame(Vtx)
      || aVn1.IsSame(Vtx) == aVn2.IsSame(Vtx))
  {
    return false;
  }
  const TopoDS_Vertex& aFoot = aVn1.IsSame(Vtx) ? aVn2 : aVn1;
  if (!aFoot.IsSame(aVf1) && !aFoot.IsSame(aVf2))
  {
    return false;
  }
  if (!Esplit.IsNull())
  {
    TopoDS_Vertex aVs1, aVs2;
    TopExp::Vertices(Esplit, aVs1, aVs2);
    if (!aFoot.IsSame(aVs1) && !aFoot.IsSame(aVs2))
    {
      return false;
    }
  }
  const BRepAdaptor_Curve aCn(Enear), aCf(Efar);
  gp_Pnt                  aP;
  gp_Vec                  aDn, aDf;
  aCn.D1(BRep_Tool::Parameter(aFoot, Enear), aP, aDn);
  aCf.D1(BRep_Tool::Parameter(aFoot, Efar), aP, aDf);
  if (aDn.Magnitude() < gp::Resolution() || aDf.Magnitude() < gp::Resolution())
  {
    return false;
  }
  // away from Vtx on the near piece, away from the foot on the far one
  if (aVn2.IsSame(Vtx))
  {
    aDn.Reverse();
  }
  if (aVf2.IsSame(aFoot))
  {
    aDf.Reverse();
  }
  return aDn.Angle(aDf) < 1.e-4;
}

//=======================================================================
// function : ExtensionOverSplit
// purpose  : The extension <C> of <Arcprol> from Vtx (at <First>) to the
//           end of the fillet's line on <Fop> (at <Last>), <Arcprol> an
//           edge of <Fprol>, a piece of the same plane beside <Fop>: it
//           starts in <Fprol> and crosses into <Fop> over one edge of the
//           two, <Esplit>, at <WQ> (<ParQ> on the edge), inside the edge,
//           and meets no other edge of either on the way. A wall kept in
//           coplanar pieces, the fillet wider than the piece beside its
//           spine. <AtEnd> when it meets the split within <TolEnd> of its
//           end -- the fillet as wide as the piece, its line along the
//           split: <WQ> is then <Last>, and all of the extension lies in
//           <Fprol>.
//=======================================================================

static bool ExtensionOverSplit(const occ::handle<Geom_Curve>& C,
                               const double                   First,
                               const double                   Last,
                               const TopoDS_Face&             Fprol,
                               const TopoDS_Face&             Fop,
                               const TopoDS_Edge&             Arcprol,
                               const TopoDS_Vertex&           Vtx,
                               const double                   Tol,
                               const double                   TolEnd,
                               TopoDS_Edge&                   Esplit,
                               double&                        WQ,
                               double&                        ParQ,
                               bool&                          AtEnd)
{
  // a crossing at the end is one: the curve a little past it
  occ::handle<Geom_Curve>              aBasis = C;
  const occ::handle<Geom_TrimmedCurve> aTrim  = occ::down_cast<Geom_TrimmedCurve>(C);
  if (!aTrim.IsNull())
  {
    aBasis = aTrim->BasisCurve();
  }
  const double            aDw = GeomAdaptor_Curve(aBasis).Resolution(Tol);
  const GeomAdaptor_Curve aSeg(aBasis, First, Last + aDw);
  const gp_Pnt            aPF = C->Value(First), aPL = C->Value(Last);
  AtEnd                       = false;
  for (int k = 0; k < 2; k++)
  {
    const TopoDS_Face& aF = k == 0 ? Fprol : Fop;
    for (TopExp_Explorer ex(aF, TopAbs_EDGE); ex.More(); ex.Next())
    {
      const TopoDS_Edge& anE = TopoDS::Edge(ex.Current());
      if (BRep_Tool::Degenerated(anE) || anE.IsSame(Arcprol))
      {
        continue;
      }
      const bool isSplit = containE(k == 0 ? Fop : Fprol, anE);
      if (isSplit && k == 1)
      {
        continue; // seen with Fprol
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
        // the extension leaves Vtx between the edges of Fprol there
        if (k == 0 && hasVertex(anE, Vtx) && aP1.Value().Distance(aPF) <= 2. * aTol)
        {
          continue;
        }
        if (!isSplit || !Esplit.IsNull() || aP1.Value().Distance(aPF) <= 2. * aTol)
        {
          return false;
        }
        double f, l;
        BRep_Tool::Range(anE, f, l);
        const double aTolP = aC.Resolution(aTol);
        if (aP2.Parameter() <= f + aTolP || aP2.Parameter() >= l - aTolP)
        {
          return false;
        }
        Esplit = anE;
        ParQ   = aP2.Parameter();
        AtEnd  = aP1.Parameter() >= Last || aP1.Value().Distance(aPL) <= TolEnd;
        WQ     = AtEnd ? Last : aP1.Parameter();
      }
    }
  }
  return !Esplit.IsNull();
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
//
//           <theAtFourEdges>: the vertex has four sharp edges, one of them
//           carrying the spine on past it, and the face at end is in two
//           coplanar pieces, the fillet wider than the one at the vertex
//           (see PerformIntersectionAtEnd). The corner is IntersectMoreCorner's
//           -- the face at end extended to the fillet, as at an end OnSame --
//           with the cut over the two pieces.
//=======================================================================

static bool ArcPastWallSplit(const TopoDS_Vertex& Vtx,
                             const TopoDS_Edge&   Arc,
                             const TopoDS_Face&   Fad,
                             const ChFiDS_Map&    VEMap,
                             const ChFiDS_Map&    EFMap,
                             const double         Tol,
                             TopoDS_Edge&         Eat);

//=======================================================================
// function : ChFi3d_PerformTangentEnd
// purpose  : The end at <Vtx> of a stripe whose two faces turn tangent there
//           (a tangency extremity), carried past it until both sides have
//           reached the face at end, tangent to the two at <Vtx> -- a chamfer
//           along the mitre of two fillets that both turn flush with a third
//           face (ChFi3d_ChBuilder::ExtentOneCorner). Each side ends on an
//           edge of that face, or in it (the side that reached it first ran
//           on over it). The last section lies in the face and is the cut:
//           the chamfer meets the face smoothly along it, and the corner of
//           the face beyond it goes. There is nothing to intersect, the
//           surfaces touching all along the section. False, and nothing
//           done, at any other end.
//=======================================================================

bool ChFi3d_PerformTangentEnd(const TopoDS_Vertex&              Vtx,
                              const occ::handle<ChFiDS_Stripe>& Stripe,
                              TopOpeBRepDS_DataStructure&       DStr,
                              const ChFiDS_Map&                 EFMap,
                              const double                      Tol3d)
{
  int        aSens   = 0;
  const int  aNum    = ChFi3d_IndexOfSurfData(Vtx, Stripe, aSens);
  const bool isfirst = (aSens == 1);
  if (!Stripe->Spine()->IsTangencyExtremity(isfirst))
  {
    return false;
  }
  const occ::handle<ChFiDS_SurfData>& Fd = Stripe->SetOfSurfData()->Value(aNum);
  if (Fd->IndexOfS1() == 0 || Fd->IndexOfS2() == 0)
  {
    return false;
  }
  const ChFiDS_CommonPoint& CV1 = Fd->Vertex(isfirst, 1);
  const ChFiDS_CommonPoint& CV2 = Fd->Vertex(isfirst, 2);
  if (!CV1.IsOnArc() && !CV2.IsOnArc())
  {
    return false;
  }
  // the face each side ends on: beyond its edge, or its own
  TopoDS_Face aFEnd[2];
  for (int k = 0; k < 2; k++)
  {
    const ChFiDS_CommonPoint& aCP = k == 0 ? CV1 : CV2;
    const TopoDS_Face         aF  = TopoDS::Face(DStr.Shape(Fd->Index(k + 1)));
    if (!aCP.IsOnArc())
    {
      aFEnd[k] = aF;
      continue;
    }
    for (NCollection_List<TopoDS_Shape>::Iterator itF(EFMap(aCP.Arc())); itF.More(); itF.Next())
    {
      if (!aF.IsSame(itF.Value()))
      {
        aFEnd[k] = TopoDS::Face(itF.Value());
      }
    }
  }
  if (aFEnd[0].IsNull() || !aFEnd[0].IsSame(aFEnd[1]) || !containV(aFEnd[0], Vtx))
  {
    return false;
  }
  TopoDS_Face Fv = aFEnd[0];
  Fv.Orientation(TopAbs_FORWARD);

  // The last section: an iso of the surface across the two lines (theirs
  // are the U isos at their two ends, ChFi3d_Builder::StoreData), from CV1's
  // side to CV2's.
  const occ::handle<Geom_Surface>& aSurf = DStr.Surface(Fd->Surf()).Surface();
  const ChFiDS_FaceInterference&   Fi1   = Fd->InterferenceOnS1();
  const ChFiDS_FaceInterference&   Fi2   = Fd->InterferenceOnS2();
  const gp_Pnt2d                   pfil1 = Fi1.PCurveOnSurf()->Value(Fi1.Parameter(isfirst));
  const gp_Pnt2d                   pfil2 = Fi2.PCurveOnSurf()->Value(Fi2.Parameter(isfirst));
  if (std::abs(pfil1.Y() - pfil2.Y()) > Precision::PConfusion()
      || std::abs(pfil1.X() - pfil2.X()) <= Precision::PConfusion())
  {
    return false;
  }
  const double                  aV   = pfil1.Y();
  const bool                    isUp = pfil1.X() < pfil2.X();
  const occ::handle<Geom_Curve> aIso = aSurf->VIso(aV);
  // the parameter on the iso, run backward when CV1's side is above
  const double                    aRev   = isUp ? 0. : aIso->ReversedParameter(0.);
  const occ::handle<Geom_Curve>   aBasis = isUp ? aIso : aIso->Reversed();
  const double                    Udeb   = isUp ? pfil1.X() : aRev - pfil1.X();
  const double                    Ufin   = isUp ? pfil2.X() : aRev - pfil2.X();
  const occ::handle<Geom_Curve>   Cc     = new Geom_TrimmedCurve(aBasis, Udeb, Ufin);
  const occ::handle<Geom2d_Curve> Ps =
    new Geom2d_Line(gp_Pnt2d(aRev, aV), isUp ? gp::DX2d() : gp_Dir2d(-1., 0.));
  const occ::handle<Geom_Surface> aSv      = BRep_Tool::Surface(Fv);
  double                          aTolProj = Tol3d;
  const occ::handle<Geom2d_Curve> Pc       = GeomProjLib::Curve2d(Cc, Udeb, Ufin, aSv, aTolProj);
  if (Pc.IsNull())
  {
    return false;
  }
  const occ::handle<GeomAdaptor_Surface> HGs = new GeomAdaptor_Surface(aSurf);
  const occ::handle<BRepAdaptor_Surface> HBs = new BRepAdaptor_Surface(Fv);
  const double tolreached                    = ChFi3d_EvalTolReached(HBs, Pc, HGs, Ps, Cc);
  // the section lies in the face, or this is not the end looked for
  if (tolreached > std::max(10. * Tol3d, std::max(CV1.Tolerance(), CV2.Tolerance())))
  {
    return false;
  }

  // the cut's orientation on the face, from a side's edge as at any corner
  // (ChFi3d_Builder::PerformOneCorner)
  const bool                on1 = CV1.IsOnArc();
  const ChFiDS_CommonPoint& aCP = on1 ? CV1 : CV2;
  TopAbs_Orientation        Et  = TopAbs_FORWARD;
  for (TopExp_Explorer ex(Fv, TopAbs_EDGE); ex.More(); ex.Next())
  {
    if (ex.Current().IsSame(aCP.Arc()))
    {
      Et = TopAbs::Compose(ex.Current().Orientation(), aCP.TransitionOnArc());
      Et = on1 ? TopAbs::Reverse(Et) : Et;
      break;
    }
  }

  const int IShape = DStr.AddShape(Fv);
  Stripe->SetIndexPoint(ChFi3d_IndexPointInDS(CV1, DStr), isfirst, 1);
  Stripe->SetIndexPoint(ChFi3d_IndexPointInDS(CV2, DStr), isfirst, 2);
  const int ICurve = DStr.AddCurve(TopOpeBRepDS_Curve(Cc, tolreached));
  DStr.ChangeShapeInterferences(IShape).Append(ChFi3d_FilCurveInDS(ICurve, IShape, Pc, Et));
  Stripe->ChangePCurve(isfirst) = Ps;
  Stripe->SetCurve(ICurve, isfirst);
  Stripe->SetParameters(isfirst, Udeb, Ufin);

  Bnd_Box box1, box2;
  ChFi3d_EnlargeBox(HBs, Pc, Udeb, Ufin, box1, box2);
  ChFi3d_EnlargeBox(DStr, Stripe, Fd, box1, box2, isfirst);
  if (CV1.IsOnArc())
  {
    ChFi3d_EnlargeBox(CV1.Arc(), EFMap(CV1.Arc()), CV1.ParameterOnArc(), box1);
  }
  if (CV2.IsOnArc())
  {
    ChFi3d_EnlargeBox(CV2.Arc(), EFMap(CV2.Arc()), CV2.ParameterOnArc(), box2);
  }
  if (!CV1.IsVertex())
  {
    ChFi3d_SetPointTolerance(DStr, box1, Stripe->IndexPoint(isfirst, 1));
  }
  if (!CV2.IsVertex())
  {
    ChFi3d_SetPointTolerance(DStr, box2, Stripe->IndexPoint(isfirst, 2));
  }
  return true;
}

// The distance of <theP> from <theS>, signed by the surface's normal.
static double SignedDistance(const gp_Pnt& theP, const occ::handle<Geom_Surface>& theS, bool& isOk)
{
  GeomAPI_ProjectPointOnSurf aProj(theP, theS);
  isOk = aProj.NbPoints() > 0;
  if (!isOk)
  {
    return 0.;
  }
  double u, v;
  aProj.LowerDistanceParameters(u, v);
  gp_Pnt aPs;
  gp_Vec aDu, aDv;
  theS->D1(u, v, aPs, aDu, aDv);
  const gp_Vec aN = aDu.Crossed(aDv);
  isOk            = aN.Magnitude() > gp::Resolution();
  return isOk ? gp_Vec(aPs, theP).Dot(aN) / aN.Magnitude() : 0.;
}

// The point of <theFil> on <theS1> and <theS2> both, by Newton's method from
// <theUV> on <theFil>: there in <theUV> and <theP>.
static bool PointOnThree(const occ::handle<Geom_Surface>& theFil,
                         const occ::handle<Geom_Surface>& theS1,
                         const occ::handle<Geom_Surface>& theS2,
                         gp_Pnt2d&                        theUV,
                         gp_Pnt&                          theP)
{
  double u = theUV.X(), v = theUV.Y();
  bool   isOk1 = true, isOk2 = true;
  auto   aF = [&](const double theU, const double theV, double& theF1, double& theF2) {
    const gp_Pnt aX = theFil->Value(theU, theV);
    theF1           = SignedDistance(aX, theS1, isOk1);
    theF2           = SignedDistance(aX, theS2, isOk2);
    return isOk1 && isOk2;
  };
  for (int it = 0; it < 50; it++)
  {
    double f1, f2;
    if (!aF(u, v, f1, f2))
    {
      return false;
    }
    if (std::max(std::abs(f1), std::abs(f2)) < 0.01 * Precision::Confusion())
    {
      theUV.SetCoord(u, v);
      theP = theFil->Value(u, v);
      return true;
    }
    const double hu = 1.e-7 * (1. + std::abs(u)), hv = 1.e-7 * (1. + std::abs(v));
    double       f1u, f2u, f1v, f2v;
    if (!aF(u + hu, v, f1u, f2u) || !aF(u, v + hv, f1v, f2v))
    {
      return false;
    }
    const double a = (f1u - f1) / hu, b = (f1v - f1) / hv;
    const double c = (f2u - f2) / hu, d = (f2v - f2) / hv;
    const double aDet = a * d - b * c;
    if (std::abs(aDet) <= gp::Resolution())
    {
      return false;
    }
    u -= (f1 * d - f2 * b) / aDet;
    v -= (a * f2 - c * f1) / aDet;
  }
  return false;
}

// The curve of <theFil> and <theS> from <theUVA> to <theUVB> on <theFil>,
// by its points: between the ends, each on the plane across the chord at it
// (PointOnThree, from the point between the ends' (u, v)), interpolated by
// chord length, its pcurve on <theFil> through the points' (u, v) on the
// same parameters (an end can lie past <theFil>'s bounds, where a projection
// finds nothing), on <theS> projected. For two surfaces tangent at an end,
// where the intersection finds no line. False where a point is not found or
// the curve is off the two by more than <theTolMax>.
static bool CurveOnTwoByPoints(const occ::handle<Geom_Surface>& theFil,
                               const occ::handle<Geom_Surface>& theS,
                               const gp_Pnt2d&                  theUVA,
                               const gp_Pnt2d&                  theUVB,
                               const double                     theTolMax,
                               occ::handle<Geom_Curve>&         theC,
                               occ::handle<Geom2d_Curve>&       thePcFil,
                               occ::handle<Geom2d_Curve>&       thePcS,
                               double&                          theTol)
{
  const int    aNb = 8;
  const gp_Pnt aA  = theFil->Value(theUVA.X(), theUVA.Y());
  const gp_Pnt aB  = theFil->Value(theUVB.X(), theUVB.Y());
  const gp_Vec aAB(aA, aB);
  if (aAB.Magnitude() <= Precision::Confusion())
  {
    return false;
  }
  occ::handle<NCollection_HArray1<gp_Pnt>>   aPnts = new NCollection_HArray1<gp_Pnt>(1, aNb + 1);
  occ::handle<NCollection_HArray1<gp_Pnt2d>> aUVs  = new NCollection_HArray1<gp_Pnt2d>(1, aNb + 1);
  occ::handle<NCollection_HArray1<double>>   aPars = new NCollection_HArray1<double>(1, aNb + 1);
  aPnts->SetValue(1, aA);
  aUVs->SetValue(1, theUVA);
  aPars->SetValue(1, 0.);
  for (int i = 1; i <= aNb; i++)
  {
    const double t = double(i) / aNb;
    gp_Pnt2d     aUV(theUVA.XY() + t * (theUVB.XY() - theUVA.XY()));
    gp_Pnt       aP = aB;
    if (i == aNb)
    {
      aUV = theUVB;
    }
    else if (!PointOnThree(theFil,
                           theS,
                           new Geom_Plane(aA.Translated(t * aAB), gp_Dir(aAB)),
                           aUV,
                           aP))
    {
      return false;
    }
    const double aStep = aP.Distance(aPnts->Value(i));
    if (aStep <= Precision::Confusion())
    {
      return false;
    }
    aPnts->SetValue(i + 1, aP);
    aUVs->SetValue(i + 1, aUV);
    aPars->SetValue(i + 1, aPars->Value(i) + aStep);
  }
  GeomAPI_Interpolate   anInterp(aPnts, aPars, false, Precision::Confusion());
  Geom2dAPI_Interpolate anInterpUV(aUVs, aPars, false, Precision::PConfusion());
  anInterp.Perform();
  anInterpUV.Perform();
  if (!anInterp.IsDone() || !anInterpUV.IsDone())
  {
    return false;
  }
  theC     = anInterp.Curve();
  thePcFil = anInterpUV.Curve();
  thePcS   = GeomProjLib::Curve2d(theC, aPars->Value(1), aPars->Value(aNb + 1), theS);
  if (thePcS.IsNull())
  {
    return false;
  }
  theTol = ChFi3d_EvalTolReached(new GeomAdaptor_Surface(theFil),
                                 thePcFil,
                                 new GeomAdaptor_Surface(theS),
                                 thePcS,
                                 theC);
  return theTol <= theTolMax;
}

// The natural normal of <theS> at <theUV>.
static gp_Vec NormalAt(const occ::handle<Geom_Surface>& theS, const gp_Pnt2d& theUV)
{
  gp_Pnt aP;
  gp_Vec aDu, aDv;
  theS->D1(theUV.X(), theUV.Y(), aP, aDu, aDv);
  return aDu.Crossed(aDv);
}

// Whether the face of surface <theS> lies, across the edge of curve <theC>
// and pcurve <thePc> on it, the way of <theKept> from the edge's middle:
// the edge is FORWARD in the face then (the face on its left).
static bool IsKeptLeft(const occ::handle<Geom_Surface>& theS,
                       const occ::handle<Geom2d_Curve>& thePc,
                       const occ::handle<Geom_Curve>&   theC,
                       const gp_Vec&                    theKept)
{
  const double w = 0.5 * (theC->FirstParameter() + theC->LastParameter());
  gp_Pnt       aP;
  gp_Vec       aT;
  theC->D1(w, aP, aT);
  return NormalAt(theS, thePc->Value(w)).Crossed(aT).Dot(theKept) > 0.;
}

// Whether the line or conic <theC>, carried on whole, never meets the
// surface of <theF>: a fillet's line that the face at its end, carried on,
// does not reach anywhere. False for any other curve, whose extension is
// not known.
static bool NeverMeets(const occ::handle<Geom_Curve>& theC, const TopoDS_Face& theF)
{
  occ::handle<Geom_Curve> aC = theC;
  while (!occ::down_cast<Geom_TrimmedCurve>(aC).IsNull())
  {
    aC = occ::down_cast<Geom_TrimmedCurve>(aC)->BasisCurve();
  }
  if (occ::down_cast<Geom_Line>(aC).IsNull() && occ::down_cast<Geom_Conic>(aC).IsNull())
  {
    return false;
  }
  occ::handle<Geom_Surface>                         aS = BRep_Tool::Surface(theF);
  const occ::handle<Geom_RectangularTrimmedSurface> aTrimmed =
    occ::down_cast<Geom_RectangularTrimmedSurface>(aS);
  if (!aTrimmed.IsNull())
  {
    aS = aTrimmed->BasisSurface();
  }
  GeomAPI_IntCS anInt(aC, aS);
  return anInt.IsDone() && anInt.NbPoints() == 0 && anInt.NbSegments() == 0;
}

//=======================================================================
// function : MakeFlatEnd
// purpose  : The end of a fillet at <Vtx> flat (ChFi3d_Builder::
//           SetFlatEndAtVertex). <Vtx> is a vertex of three faces; side
//           <Ip> of <Fd> ends on the edge of <Vtx> between its face and the
//           face at end (Epiv), and the other side's line runs on its face
//           (Fo). The face at end, carried on, would meet that line only
//           past <Vtx>, along its other edge at <Vtx> (Arcprol, which Fo
//           shares) carried on: the far side of a pocket's acute tip, where
//           the cap would cut a groove into the body beyond the tip. Or,
//           <theMissed>, the cap was not made (IntersUpdateOnSame, and the
//           cut off Fop, found nothing) and Fo's line, a line or a conic, the
//           face at end carried on meets nowhere (#962's ring). The
//           fillet ends instead in the plane through <Vtx> and the first
//           side's point on Epiv, as square to the spine as that allows. The
//           face it leaves there, flat, lies between the section, the piece
//           of Epiv from <Vtx> to that point, which the face at end keeps,
//           and the line from <Vtx> to the other side's point on Fo. All of
//           it is put in the DS here, and <EdgeSpine>, the edge filleted,
//           ended at <Vtx>, which stays. False, and nothing done, where this
//           is not the case or the face cannot be made: the end is then the
//           usual one.
//=======================================================================

static bool MakeFlatEnd(
  TopOpeBRepDS_DataStructure&                                                        DStr,
  const occ::handle<ChFiDS_Stripe>&                                                  Stripe,
  const occ::handle<ChFiDS_SurfData>&                                                Fd,
  const bool                                                                         isFirst,
  const TopoDS_Vertex&                                                               Vtx,
  const TopoDS_Edge&                                                                 EdgeSpine,
  const ChFiDS_Map&                                                                  VEMap,
  const ChFiDS_Map&                                                                  EFMap,
  const ChFiDS_Map&                                                                  VFMap,
  const double                                                                       Tol3d,
  const double                                                                       Tol2d,
  NCollection_DataMap<TopoDS_Shape, NCollection_List<int>, TopTools_ShapeMapHasher>& EVIMap,
  const bool                                                                         theMissed)
{
  if (Fd->Index(1) <= 0 || Fd->Index(2) <= 0 || Fd->Surf() <= 0
      || DStr.Shape(Fd->Index(1)).ShapeType() != TopAbs_FACE
      || DStr.Shape(Fd->Index(2)).ShapeType() != TopAbs_FACE || !VFMap.Contains(Vtx)
      || ChFi3d_nbface(VFMap(Vtx)) != 3)
  {
    return false;
  }
  // the side on an edge of Vtx, the only one
  int Ip = 0;
  for (int i = 1; i <= 2; i++)
  {
    const ChFiDS_CommonPoint& aCV = Fd->Vertex(isFirst, i);
    if (aCV.IsOnArc() && !aCV.IsVertex() && hasVertex(aCV.Arc(), Vtx))
    {
      Ip = Ip == 0 ? i : -1;
    }
  }
  if (Ip <= 0)
  {
    return false;
  }
  const int                 Io  = 3 - Ip;
  const ChFiDS_CommonPoint& CVp = Fd->Vertex(isFirst, Ip);
  ChFiDS_CommonPoint&       CVo = Fd->ChangeVertex(isFirst, Io);
  if (CVo.IsVertex())
  {
    return false;
  }
  const TopoDS_Face              Fp   = TopoDS::Face(DStr.Shape(Fd->Index(Ip)));
  const TopoDS_Face              Fo   = TopoDS::Face(DStr.Shape(Fd->Index(Io)));
  const TopoDS_Edge              Epiv = CVp.Arc();
  const ChFiDS_FaceInterference& Fip  = Fd->Interference(Ip);
  const ChFiDS_FaceInterference& Fio  = Fd->Interference(Io);
  if (Fip.LineIndex() == 0 || Fio.LineIndex() == 0 || Fip.PCurveOnSurf().IsNull()
      || Fio.PCurveOnSurf().IsNull() || Fio.PCurveOnFace().IsNull() || !EFMap.Contains(Epiv)
      || EFMap(Epiv).Extent() != 2)
  {
    return false;
  }
  // the face at end, across Epiv, and its other edge at Vtx, Fo's too
  TopoDS_Face Fv;
  for (NCollection_List<TopoDS_Shape>::Iterator It(EFMap(Epiv)); It.More(); It.Next())
  {
    if (!Fp.IsSame(It.Value()))
    {
      Fv = TopoDS::Face(It.Value());
    }
  }
  if (Fv.IsNull() || Fv.IsSame(Fo))
  {
    return false;
  }
  TopoDS_Edge Arcprol;
  for (NCollection_List<TopoDS_Shape>::Iterator It(VEMap(Vtx)); It.More(); It.Next())
  {
    const TopoDS_Edge& E = TopoDS::Edge(It.Value());
    if (!E.IsSame(Epiv) && !E.IsSame(EdgeSpine) && !BRep_Tool::Degenerated(E) && containE(Fv, E)
        && containE(Fo, E))
    {
      Arcprol = E;
    }
  }
  if (Arcprol.IsNull())
  {
    return false;
  }
  // the edge as it lies in the face, whatever way up the face is
  auto anOrIn = [](const TopoDS_Edge& theE, const TopoDS_Face& theF) {
    for (TopExp_Explorer ex(theF.Oriented(TopAbs_FORWARD), TopAbs_EDGE); ex.More(); ex.Next())
    {
      if (theE.IsSame(ex.Current()))
      {
        return ex.Current().Orientation();
      }
    }
    return theE.Orientation();
  };
  const gp_Pnt PV = BRep_Tool::Pnt(Vtx);
  const gp_Pnt PP = CVp.Point();

  // the spine's way from Vtx
  BRepAdaptor_Curve aCsp(EdgeSpine);
  gp_Pnt            aPbid;
  gp_Vec            aTS;
  aCsp.D1(BRep_Tool::Parameter(Vtx, EdgeSpine), aPbid, aTS);
  if (aTS.Magnitude() <= gp::Resolution())
  {
    return false;
  }
  if (gp_Vec(PV, aCsp.Value(0.5 * (aCsp.FirstParameter() + aCsp.LastParameter()))).Dot(aTS) < 0.)
  {
    aTS.Reverse();
  }
  aTS.Normalize();

  // The cap in the face at end lands past Vtx: Fo's line, along its tangent
  // at the end, meets Arcprol's line on the far side of Vtx, and past Vtx's
  // section along the spine. (Behind Vtx on Arcprol's line, yet short of the
  // section, the cap leans back over the fillet and cuts nothing past Vtx;
  // in the section, the face at end is square to the spine already.)
  const occ::handle<Geom_Curve>& aLo = DStr.Curve(Fio.LineIndex()).Curve();
  if (aLo.IsNull())
  {
    return false;
  }
  const double      wo0 = Fio.Parameter(isFirst);
  double            wo  = wo0;
  gp_Pnt            Po  = CVo.Point();
  // A line run on past the point (a KPart's, to its face's edge): the end is
  // the point, where it lies on the line
  if (aLo->Value(wo).Distance(Po) > Tol3d)
  {
    GeomAPI_ProjectPointOnCurve aProjO(Po, aLo, Fio.FirstParameter(), Fio.LastParameter());
    if (aProjO.NbPoints() == 0 || aProjO.LowerDistance() > 10 * Tol3d)
    {
      return false;
    }
    wo = aProjO.LowerDistanceParameter();
  }
  BRepAdaptor_Curve aCpr(Arcprol);
  const double      parVpr = BRep_Tool::Parameter(Vtx, Arcprol);
  gp_Vec            aDpr, aTo;
  aCpr.D1(parVpr, aPbid, aDpr);
  if (std::abs(parVpr - aCpr.LastParameter()) < std::abs(parVpr - aCpr.FirstParameter()))
  {
    aDpr.Reverse();
  }
  aLo->D1(wo, aPbid, aTo);
  const gp_Pnt2d aUVo = BRep_Tool::Parameters(Vtx, Fo);
  const gp_Vec   aNo  = NormalAt(BRep_Tool::Surface(Fo), aUVo);
  if (aDpr.Magnitude() <= gp::Resolution() || aTo.Magnitude() <= gp::Resolution()
      || aNo.Magnitude() <= gp::Resolution())
  {
    return false;
  }
  aDpr.Normalize();
  aTo.Normalize();
  const gp_Vec aNoU = aNo.Normalized();
  const double aDen = aTo.Crossed(aDpr).Dot(aNoU);
  if (std::abs(aDen) <= 1.e-6)
  {
    return false;
  }
  const gp_Pnt aPl  = aLo->Value(wo);
  const gp_Pnt aX   = aPl.Translated(aTo * (-gp_Vec(PV, aPl).Crossed(aDpr).Dot(aNoU) / aDen));
  const double aLam = gp_Vec(PV, aX).Dot(aDpr);
  const double aSig = gp_Vec(PV, aX).Dot(aTS);
  // Or, asked once the cap was not made, Fo's line, carried on whole, never
  // meets the face at end. #962's ring: its rim's fillet ends where an arm's
  // side leaves the rim 18.6 deg off its tangent, and the line on the floor,
  // a circle of radius R - r, misses the side's plane, carried back through
  // the ring, from r = R (1 - cos 18.6 deg), 1.03 -- the tangent above still
  // crosses Arcprol's line, behind Vtx.
  const bool isMissed = aLam >= -10 * Tol3d || aSig >= -10 * Tol3d;
  if (isMissed != theMissed || (isMissed && !NeverMeets(aLo, Fv)))
  {
    return false;
  }

  // the plane, its normal the way of the spine: out of the face it bounds
  gp_Vec aChord(PV, PP);
  if (aChord.Magnitude() <= 10 * Tol3d)
  {
    return false;
  }
  aChord.Normalize();
  gp_Vec aN = aTS - aChord * aTS.Dot(aChord);
  if (aN.Magnitude() <= 0.1)
  {
    return false;
  }
  const occ::handle<Geom_Plane> aGPl = new Geom_Plane(PV, gp_Dir(aN));
  const gp_Pln                  aPln = aGPl->Pln();
  // nor along the face at end, there
  {
    double                          aF0, aL0;
    const occ::handle<Geom2d_Curve> aPC = PCurveInFace(Epiv, Fv, aF0, aL0);
    if (aPC.IsNull())
    {
      return false;
    }
    const gp_Vec aNv = NormalAt(BRep_Tool::Surface(Fv), aPC->Value(CVp.ParameterOnArc()));
    if (aNv.Magnitude() <= gp::Resolution()
        || aNv.Normalized().Crossed(gp_Vec(aPln.Axis().Direction())).Magnitude() <= 1.e-3)
    {
      return false;
    }
  }
  // Epiv from Vtx to the point lies in it
  double                        aFE, aLE;
  const occ::handle<Geom_Curve> aCpiv = BRep_Tool::Curve(Epiv, aFE, aLE);
  const double                  parV  = BRep_Tool::Parameter(Vtx, Epiv);
  const double                  parP  = CVp.ParameterOnArc();
  if (aCpiv.IsNull() || std::abs(parP - parV) <= Precision::PConfusion())
  {
    return false;
  }
  for (int k = 1; k <= 3; k++)
  {
    if (aPln.Distance(aCpiv->Value(parV + 0.25 * k * (parP - parV))) > Tol3d)
    {
      return false;
    }
  }
  // Fo's line across the plane: at its end where the section is square to
  // the spine, else where it crosses
  if (aPln.Distance(Po) > Tol3d)
  {
    GeomAPI_IntCS anInt(new Geom_TrimmedCurve(aLo, Fio.FirstParameter(), Fio.LastParameter()),
                        aGPl);
    if (!anInt.IsDone() || anInt.NbPoints() == 0)
    {
      return false;
    }
    // Where no cap could be made, the point in the air is not yet the
    // line's (#962's ring: on the far side of the gap, its crossing nearest
    // a turn away): the crossing nearest Vtx
    const double aWs   = wo;
    double       aBest = RealLast();
    for (int i = 1; i <= anInt.NbPoints(); i++)
    {
      double aU, aV, aW;
      anInt.Parameters(i, aU, aV, aW);
      const double aD = isMissed ? PV.Distance(anInt.Point(i)) : std::abs(aW - aWs);
      if (aD < aBest)
      {
        aBest = aD;
        wo    = aW;
      }
    }
    Po = aLo->Value(wo);
  }
  if (aPln.Distance(Po) > 10 * Tol3d || PV.Distance(Po) <= 10 * Tol3d)
  {
    return false;
  }
  // the point within Fo, and the line to it from Vtx: Fo's, near enough
  const gp_Pnt2d          aUVPo = Fio.PCurveOnFace()->Value(wo);
  BRepTopAdaptor_FClass2d aClo(Fo, Precision::PConfusion());
  if (aClo.Perform(aUVPo) != TopAbs_IN)
  {
    return false;
  }
  const double                  aLt = PV.Distance(Po);
  const occ::handle<Geom_Curve> aCt =
    new Geom_TrimmedCurve(new Geom_Line(PV, gp_Dir(gp_Vec(PV, Po))), 0., aLt);
  const occ::handle<Geom_Surface> aSo    = BRep_Tool::Surface(Fo);
  double                          aTolPj = Tol3d;
  const occ::handle<Geom2d_Curve> aPctO  = GeomProjLib::Curve2d(aCt, 0., aLt, aSo, aTolPj);
  aTolPj                                 = Tol3d;
  const occ::handle<Geom2d_Curve> aPctL  = GeomProjLib::Curve2d(aCt, 0., aLt, aGPl, aTolPj);
  if (aPctO.IsNull() || aPctL.IsNull() || aClo.Perform(aPctO->Value(0.5 * aLt)) != TopAbs_IN)
  {
    return false;
  }
  const occ::handle<GeomAdaptor_Surface> aHSo = new GeomAdaptor_Surface(aSo);
  double                                 aUP0, aVP0, aUP1, aVP1, aUP2, aVP2;
  ElSLib::Parameters(aPln, PV, aUP0, aVP0);
  ElSLib::Parameters(aPln, PP, aUP1, aVP1);
  ElSLib::Parameters(aPln, Po, aUP2, aVP2);
  const double                           aMarg = PV.Distance(PP) + aLt;
  const occ::handle<GeomAdaptor_Surface> aHPl =
    new GeomAdaptor_Surface(aGPl,
                            std::min({aUP0, aUP1, aUP2}) - aMarg,
                            std::max({aUP0, aUP1, aUP2}) + aMarg,
                            std::min({aVP0, aVP1, aVP2}) - aMarg,
                            std::max({aVP0, aVP1, aVP2}) + aMarg);
  const double aTolT =
    std::max(ChFi3d_EvalTolReached(aHSo, aPctO, aHPl, aPctL, aCt), Precision::Confusion());
  if (aTolT > 10 * Tol3d)
  {
    return false;
  }

  // the section, from side 1 to side 2, as the corners put their curves
  const occ::handle<Geom_Surface>& aSfil = DStr.Surface(Fd->Surf()).Surface();
  occ::handle<GeomAdaptor_Surface> aHGs  = new GeomAdaptor_Surface(aSfil);
  const gp_Pnt2d                   aUVsp = Fip.PCurveOnSurf()->Value(Fip.Parameter(isFirst));
  const gp_Pnt2d                   aUVso = Fio.PCurveOnSurf()->Value(wo);
  if (isMissed)
  {
    // The fillet's surface about the section's ends only: a torus meets a
    // plane through its axis in whole circles, and the section ran the long
    // way round the tube (#962's ring, three quarters of it)
    const double aDu = std::abs(aUVsp.X() - aUVso.X()), aDv = std::abs(aUVsp.Y() - aUVso.Y());
    const double aMarg = std::max({aDu, aDv, 1.e-3});
    aHGs               = new GeomAdaptor_Surface(aSfil,
                                   std::min(aUVsp.X(), aUVso.X()) - aMarg,
                                   std::max(aUVsp.X(), aUVso.X()) + aMarg,
                                   std::min(aUVsp.Y(), aUVso.Y()) - 0.1 * aMarg,
                                   std::max(aUVsp.Y(), aUVso.Y()) + 0.1 * aMarg);
  }
  const gp_Pnt2d             aUVlp(aUP1, aVP1), aUVlo(aUP2, aVP2);
  const gp_Pnt2d&            aUVs1 = Ip == 1 ? aUVsp : aUVso;
  const gp_Pnt2d&            aUVs2 = Ip == 1 ? aUVso : aUVsp;
  const gp_Pnt2d&            aUVl1 = Ip == 1 ? aUVlp : aUVlo;
  const gp_Pnt2d&            aUVl2 = Ip == 1 ? aUVlo : aUVlp;
  NCollection_Array1<double> Pardeb(1, 4), Parfin(1, 4);
  Pardeb(1) = aUVs1.X();
  Pardeb(2) = aUVs1.Y();
  Pardeb(3) = aUVl1.X();
  Pardeb(4) = aUVl1.Y();
  Parfin(1) = aUVs2.X();
  Parfin(2) = aUVs2.Y();
  Parfin(3) = aUVl2.X();
  Parfin(4) = aUVl2.Y();
  occ::handle<Geom_Curve>   aCc;
  occ::handle<Geom2d_Curve> aPs, aPc;
  double                    aTolR = Tol3d;
  if (!ChFi3d_ComputeCurves(aHGs, aHPl, Pardeb, Parfin, aCc, aPs, aPc, Tol3d, Tol2d, aTolR, false)
      || aCc.IsNull() || aPs.IsNull() || aPc.IsNull())
  {
    return false;
  }
  aTolR =
    std::min(aTolR,
             std::max(ChFi3d_EvalTolReached(aHGs, aPs, aHPl, aPc, aCc), Precision::Confusion()));
  const double p1 = aCc->FirstParameter(), p2 = aCc->LastParameter();
  if (aTolR > 10 * Tol3d || aCc->Value(p1).Distance(Ip == 1 ? PP : Po) > 10 * Tol3d
      || aCc->Value(p2).Distance(Ip == 1 ? Po : PP) > 10 * Tol3d)
  {
    return false;
  }

  // Epiv's piece from Vtx to the point
  double                          aFV, aLV;
  const occ::handle<Geom2d_Curve> aPCpivV = PCurveInFace(Epiv, Fv, aFV, aLV);
  if (aPCpivV.IsNull())
  {
    return false;
  }
  const double                    aP0 = std::min(parV, parP), aP1 = std::max(parV, parP);
  const occ::handle<Geom_Curve>   aCp   = new Geom_TrimmedCurve(aCpiv, aP0, aP1);
  const occ::handle<Geom2d_Curve> aPcpV = new Geom2d_TrimmedCurve(aPCpivV, aP0, aP1);
  aTolPj                                = Tol3d;
  const occ::handle<Geom2d_Curve> aPcpL = GeomProjLib::Curve2d(aCp, aP0, aP1, aGPl, aTolPj);
  if (aPcpL.IsNull())
  {
    return false;
  }
  const double aTolCp = std::max(BRep_Tool::Tolerance(Epiv), Precision::Confusion());

  // The sides. The flat face lies, from the section, toward Vtx; from the
  // line on Fo, toward Epiv; from Epiv, toward the line. What Fo keeps lies
  // past the line, against the spine; what the face at end keeps, as it
  // lies on Epiv.
  const TopAbs_Orientation anOrSL =
    IsKeptLeft(aGPl, aPc, aCc, gp_Vec(aCc->Value(0.5 * (p1 + p2)), PV)) ? TopAbs_FORWARD
                                                                        : TopAbs_REVERSED;
  const TopAbs_Orientation anOrTL =
    IsKeptLeft(aGPl, aPctL, aCt, gp_Vec(PV, PP)) ? TopAbs_FORWARD : TopAbs_REVERSED;
  const TopAbs_Orientation anOrTO =
    IsKeptLeft(aSo, aPctO, aCt, -aTS) ? TopAbs_FORWARD : TopAbs_REVERSED;
  const TopAbs_Orientation anOrPL =
    IsKeptLeft(aGPl, aPcpL, aCp, gp_Vec(PV, Po)) ? TopAbs_FORWARD : TopAbs_REVERSED;
  const TopAbs_Orientation anOrPV = anOrIn(Epiv, Fv);
  // on the fillet, as the cut at an end goes
  TopAbs_Orientation aTrafil1 = TopAbs_FORWARD;
  if (Fd->IndexOfS1() > 0)
  {
    aTrafil1 = DStr.Shape(Fd->IndexOfS1()).Orientation();
  }
  aTrafil1 = TopAbs::Compose(aTrafil1, Fd->Orientation());
  aTrafil1 = TopAbs::Compose(TopAbs::Reverse(Fd->InterferenceOnS1().Transition()), aTrafil1);
  const TopAbs_Orientation anOrS = isFirst ? TopAbs::Reverse(aTrafil1) : aTrafil1;

  // -- the DS, from here on --
  const int Isurf  = Fd->Surf();
  const int IsurfC = DStr.AddSurface(TopOpeBRepDS_Surface(aGPl, Tol3d));
  if (!EVIMap.IsBound(Vtx))
  {
    EVIMap.Bind(Vtx, NCollection_List<int>());
  }
  EVIMap.ChangeFind(Vtx).Append(IsurfC);
  const int SolInd = Stripe->SolidIndex();
  DStr.ChangeShapeInterferences(SolInd).Append(
    new TopOpeBRepDS_SolidSurfaceInterference(TopOpeBRepDS_Transition(TopAbs_FORWARD),
                                              TopOpeBRepDS_SOLID,
                                              SolInd,
                                              TopOpeBRepDS_SURFACE,
                                              IsurfC));
  // the points
  if (wo != wo0)
  {
    Fd->ChangeInterference(Io).SetParameter(wo, isFirst);
  }
  if (Po.Distance(CVo.Point()) > Precision::Confusion())
  {
    CVo.Reset();
    CVo.SetPoint(Po);
    CVo.SetTolerance(std::max(aPln.Distance(Po), Precision::Confusion()));
  }
  Stripe->SetIndexPoint(ChFi3d_IndexPointInDS(CVp, DStr), isFirst, Ip);
  Stripe->SetIndexPoint(ChFi3d_IndexPointInDS(CVo, DStr), isFirst, Io);
  const int IndP = Stripe->IndexPoint(isFirst, Ip);
  const int IndO = Stripe->IndexPoint(isFirst, Io);
  const int IVtx = DStr.AddShape(Vtx);
  {
    TopOpeBRepDS_Point& aTP = DStr.ChangePoint(IndP);
    aTP.Tolerance(std::max({aTP.Tolerance(), aTolR, aTolCp}));
    TopOpeBRepDS_Point& aTO = DStr.ChangePoint(IndO);
    aTO.Tolerance(std::max({aTO.Tolerance(), aTolR, aTolT}));
  }
  // the section
  const int ICs = DStr.AddCurve(TopOpeBRepDS_Curve(aCc, aTolR));
  DStr.ChangeCurveInterferences(ICs).Append(
    ChFi3d_FilPointInDS(TopAbs_FORWARD, ICs, Ip == 1 ? IndP : IndO, p1));
  DStr.ChangeCurveInterferences(ICs).Append(
    ChFi3d_FilPointInDS(TopAbs_REVERSED, ICs, Ip == 1 ? IndO : IndP, p2));
  DStr.ChangeSurfaceInterferences(Isurf).Append(ChFi3d_FilCurveInDS(ICs, Isurf, aPs, anOrS));
  DStr.ChangeSurfaceInterferences(IsurfC).Append(ChFi3d_FilCurveInDS(ICs, IsurfC, aPc, anOrSL));
  // the line on Fo, from Vtx
  const int ICt = DStr.AddCurve(TopOpeBRepDS_Curve(aCt, aTolT));
  DStr.ChangeCurveInterferences(ICt).Append(
    ChFi3d_FilPointInDS(TopAbs_FORWARD, ICt, IVtx, 0., true));
  DStr.ChangeCurveInterferences(ICt).Append(ChFi3d_FilPointInDS(TopAbs_REVERSED, ICt, IndO, aLt));
  const int IFo = DStr.AddShape(Fo);
  DStr.ChangeShapeInterferences(IFo).Append(ChFi3d_FilCurveInDS(ICt, IFo, aPctO, anOrTO));
  DStr.ChangeSurfaceInterferences(IsurfC).Append(ChFi3d_FilCurveInDS(ICt, IsurfC, aPctL, anOrTL));
  // Epiv's piece, which the face at end keeps and the flat face shares
  const int  ICp      = DStr.AddCurve(TopOpeBRepDS_Curve(aCp, aTolCp));
  const bool isVtxLow = parV < parP;
  DStr.ChangeCurveInterferences(ICp).Append(
    ChFi3d_FilPointInDS(TopAbs_FORWARD, ICp, isVtxLow ? IVtx : IndP, aP0, isVtxLow));
  DStr.ChangeCurveInterferences(ICp).Append(
    ChFi3d_FilPointInDS(TopAbs_REVERSED, ICp, isVtxLow ? IndP : IVtx, aP1, !isVtxLow));
  const int IFv = DStr.AddShape(Fv);
  DStr.ChangeShapeInterferences(IFv).Append(ChFi3d_FilCurveInDS(ICp, IFv, aPcpV, anOrPV));
  DStr.ChangeSurfaceInterferences(IsurfC).Append(ChFi3d_FilCurveInDS(ICp, IsurfC, aPcpL, anOrPL));
  // Epiv cut at the point
  const int IEpiv = DStr.AddShape(Epiv);
  DStr.ChangeShapeInterferences(IEpiv).Append(
    ChFi3d_FilPointInDS(CVp.TransitionOnArc(), IEpiv, IndP, CVp.ParameterOnArc()));
  // Vtx stays, with Arcprol and the flat face's edges: the edge filleted
  // ends there with nothing kept of it
  TopAbs_Orientation anOrV = TopAbs_FORWARD;
  for (TopExp_Explorer ex(EdgeSpine.Oriented(TopAbs_FORWARD), TopAbs_VERTEX); ex.More(); ex.Next())
  {
    if (Vtx.IsSame(ex.Current()))
    {
      anOrV = ex.Current().Orientation();
      break;
    }
  }
  const int ISpine = DStr.AddShape(EdgeSpine);
  DStr.ChangeShapeInterferences(ISpine).Append(
    ChFi3d_FilVertexInDS(TopAbs::Reverse(anOrV),
                         ISpine,
                         IVtx,
                         BRep_Tool::Parameter(Vtx, EdgeSpine)));
  Stripe->InDS(isFirst);
  return true;
}

void ChFi3d_Builder::PerformOneCorner(const int Index, const bool thePrepareOnSame)
{
  PerformOneCorner(Index, thePrepareOnSame, false);
}

void ChFi3d_Builder::PerformOneCorner(const int  Index,
                                      const bool thePrepareOnSame,
                                      const bool theAtFourEdges)
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
  bool        onsame = (stat == ChFiDS_OnSame) || theAtFourEdges;
  TopoDS_Face Fv, Fad, Fop;
  TopoDS_Edge Arcpiv, Arcprol, Arcspine;
  // At four edges: the edge of Vtx on Fop that is not the spine's, which
  // the extension carries on in Fop as it does Arcprol in Fv, and Vtx's
  // orientation in it
  TopoDS_Edge        Arcprolbis;
  TopAbs_Orientation OVtxBis = TopAbs_FORWARD;
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
  // The same cut at an end that is not OnSame (see below): it crosses the
  // edge between the two faces itself, of this length from Vp
  bool   cutEtg = false;
  double lenEtg = 0.;
  // The cut crosses Etg's line carried on within this of Vp: no edge that
  // short is made of the line (BOP's check takes one under some 1e-6 for no
  // edge at all) -- Etg starts at the crossing point, which takes Vp's
  // place. So it does when the cut crosses Etg itself too near its end to
  // be cut there.
  const double aTolVp = 100 * Precision::Confusion();
  bool         atVp   = false;
  // The fillet's line over a wall in coplanar pieces (see LineOverSplit):
  // its side, the line, the edge it crosses and the one it ends on, the
  // piece beyond, and the parameters of the two points
  int                    onsOS = 0;
  occ::handle<Geom_Line> LinOS;
  TopoDS_Edge            EsplitOS, EarcOS;
  TopoDS_Face            FnOS;
  double                 wQOS = 0., wPOS = 0., parQOS = 0., parPOS = 0.;
  // The cut that runs off Fop (see below): the edge of Fop the fillet's line
  // there ends on (at Q, the stripe's point), the face beyond it, the point
  // P where the cut on Fv's surface crosses that face's, and R where the
  // edge crosses Fv's surface -- with their parameters on the fillet, Fv,
  // Fob and Fop, and R's on Eob
  bool                             offFop = false;
  TopoDS_Edge                      Eob;
  TopoDS_Face                      Fob;
  occ::handle<BRepAdaptor_Surface> HBob = new BRepAdaptor_Surface();
  gp_Pnt                           POff, ROff;
  gp_Pnt2d                         pfilP, pvP, pobP, pvR, pobR, popR, pobQ;
  double                           parR = 0.;
  // the orientations of Fob and Fv in the shell
  TopAbs_Orientation OFob = TopAbs_FORWARD, OFvSh = TopAbs_FORWARD;
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

  // the end as it was, for the flat end where no cap is made (below)
  const ChFiDS_CommonPoint      aCVIn1 = Fd->Vertex(isfirst, 1), aCVIn2 = Fd->Vertex(isfirst, 2);
  const ChFiDS_FaceInterference aFiIn1 = Fd->Interference(1), aFiIn2 = Fd->Interference(2);
  if (onsame)
  {
    // Where the cap in the face at end would cut into the body past Vtx,
    // the end flat, if asked for (SetFlatEndAtVertex)
    if (!thePrepareOnSame && !theAtFourEdges && stat == ChFiDS_OnSame
        && ChFi3d_Builder::FlatEndAtVertex()
        && MakeFlatEnd(DStr,
                       stripe,
                       Fd,
                       isfirst,
                       Vtx,
                       Arcspine,
                       myVEMap,
                       myEFMap,
                       myVFMap,
                       tolapp3d,
                       tol2d,
                       myEVIMap,
                       false))
    {
      return;
    }
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
      // At four edges: the side on an edge past the foot of the split of
      // the face at end is on Vtx as on that face in one piece.
      for (int ons = 1; ons <= 2 && theAtFourEdges && !sur1 && !sur2; ons++)
      {
        TopoDS_Edge anEat;
        if (ArcPastWallSplit(Vtx,
                             ons == 1 ? CV1.Arc() : CV2.Arc(),
                             TopoDS::Face(DStr.Shape(Fd->Index(ons))),
                             myVEMap,
                             myEFMap,
                             10 * tolapp3d,
                             anEat))
        {
          (ons == 1 ? sur1 : sur2) = true;
        }
      }
      for (int ons = 1; ons <= 2 && sur1 != sur2 && !theAtFourEdges; ons++)
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
        // (the wall in pieces is the fillet's face, or the one the arc bounds)
        TopoDS_Edge anEV = ChFi3d_EdgeOnSplitToVertex(anArc, aF, anOther, Vtx, myEFMap, myVEMap);
        if (anEV.IsNull())
        {
          anEV = ChFi3d_EdgeOnSplitToVertex(anArc, anOther, aF, Vtx, myEFMap, myVEMap);
        }
        if (!anEV.IsNull())
        {
          (ons == 1 ? EV1 : EV2)   = anEV;
          (ons == 1 ? sur1 : sur2) = true;
        }
      }
      if (sur1 && sur2 && !theAtFourEdges)
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

    // The fillet's line on Fad passes the end of the edge between Fad and
    // the face at Vtx (a drafted wall's rounded corner, a cone) onto an
    // edge of the face beyond it (the wall's plane), into which the face at
    // Vtx runs on, tangent, across a straight edge from Arcpiv's vertex:
    // one wall in two faces, and Fv does not hold Vtx. The cut runs over
    // both, from Arcpiv on the face beyond (FvT) across that edge carried
    // on past the vertex, to Fop on the face at Vtx, which is Fv from here.
    // (looked for before the other side's arc is tried for the face at end:
    // that one can lie far off, where the line ran on past the spine's end)
    auto aTwoFaces = [&]() {
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
    };

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
    if (!isinface)
    {
      aTwoFaces();
    }
    if (!isinface && FvT.IsNull() && Fd->Vertex(isfirst, 3 - IFadArc).IsOnArc())
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
      aTwoFaces();
    }

    if (Fv.IsNull())
    {
      throw StdFail_NotDone("OneCorner : face at end not found");
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
        if (Arcpiv.IsSame(It.Value()) || Eunder.IsSame(It.Value())
            || (theAtFourEdges && Arcspine.IsSame(It.Value())))
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
      // PerformIntersectionAtEnd is what sent the four-edge case here
      if (theAtFourEdges)
      {
        throw Standard_Failure("OneCorner : the cut over two pieces has no edge to extend");
      }
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
    if (theAtFourEdges)
    {
      if (FvT.IsNull())
      {
        throw Standard_Failure("OneCorner : the face at end is not in two pieces");
      }
      for (It.Initialize(myVEMap(Vtx)); It.More() && Arcprolbis.IsNull(); It.Next())
      {
        if (Arcprol.IsSame(It.Value()) || Arcspine.IsSame(It.Value())
            || Arcpiv.IsSame(It.Value()) || Eunder.IsSame(It.Value()))
        {
          continue;
        }
        for (ex.Init(Fop, TopAbs_EDGE); ex.More(); ex.Next())
        {
          if (It.Value().IsSame(ex.Current()))
          {
            Arcprolbis = TopoDS::Edge(It.Value());
            OArcprolop = ex.Current().Orientation();
            break;
          }
        }
      }
      if (Arcprolbis.IsNull())
      {
        throw Standard_Failure("OneCorner : no edge of the vertex carries the spine on");
      }
      for (ex.Init(Arcprolbis.Oriented(TopAbs_FORWARD), TopAbs_VERTEX); ex.More(); ex.Next())
      {
        if (Vtx.IsSame(ex.Current()))
        {
          OVtxBis = ex.Current().Orientation();
          break;
        }
      }
    }
    for (ex.Init(Fop, TopAbs_EDGE); ex.More() && Arcprolbis.IsNull(); ex.Next())
    {
      if (Arcprol.IsSame(ex.Current()))
      {
        OArcprolop = ex.Current().Orientation();
        break;
      }
    }
    if (!ex.More() && Arcprolbis.IsNull())
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

    // The stripe ended where its side on Fad met Arcpiv, its point on Fop
    // inside Fop, and the fillet's line there runs on to an edge of Fop
    // before it reaches Fv's surface: the line's end on that edge is where
    // the cut runs off Fop, as if the stripe had ended there (#273's
    // Fillet001, E40 -- the mirror of E39, whose side on Fop reached the
    // edge first). Taken back below when the cut off Fop is not made.
    const ChFiDS_CommonPoint aCPopIn     = CPopArc;
    bool                     isCPopMoved = false;
    if (!inters && !CPopArc.IsOnArc() && !c3df.IsNull() && FvT.IsNull() && Etan.IsNull()
        && Arcprolbis.IsNull())
    {
      const double aParEnd = FiopArc.Parameter(isfirst);
      const gp_Pnt aPEnd   = c3df->Value(aParEnd);
      TopoDS_Edge  anEnd;
      double       aParOnEnd = 0., aTolEnd = 0., aDMin = RealLast();
      for (ex.Init(Fop, TopAbs_EDGE); ex.More(); ex.Next())
      {
        const TopoDS_Edge& anE = TopoDS::Edge(ex.Current());
        if (BRep_Tool::Degenerated(anE) || hasVertex(anE, Vtx) || myEFMap(anE).Extent() != 2)
        {
          continue;
        }
        BRepAdaptor_Curve aCE(anE);
        Extrema_ExtPC     anExt(aPEnd, aCE);
        const double      aTol = std::max(BRep_Tool::Tolerance(anE), 10. * tolapp3d);
        for (int i = 1; anExt.IsDone() && i <= anExt.NbExt(); i++)
        {
          const double aD = std::sqrt(anExt.SquareDistance(i));
          if (aD <= aTol && aD < aDMin)
          {
            aDMin     = aD;
            anEnd     = anE;
            aParOnEnd = anExt.Point(i).Parameter();
            aTolEnd   = aTol;
          }
        }
      }
      if (!anEnd.IsNull())
      {
        // the transition the walk gives a line ending on an arc (anEnd is
        // as oriented in Fop)
        TopAbs_Orientation aTr = TopAbs::Compose(FiopArc.Transition(), anEnd.Orientation());
        if (isfirst)
        {
          aTr = TopAbs::Reverse(aTr);
        }
        CPopArc.Reset();
        CPopArc.SetPoint(aPEnd);
        CPopArc.SetArc(aTolEnd, anEnd, aParOnEnd, aTr);
        isCPopMoved = true;
      }
    }

    // The cut runs off Fop: the fillet's line there, carried on past Vtx,
    // ends on an edge of Fop (Eob) before it reaches Fv's surface, and the
    // cap in Fv's surface would run past Fop too (#273's pocket in the
    // handle's top, its apex 22 deg: a rim filleted at r, the cap in the
    // other wall r / sin(22 deg) past the apex, the top 1.30 -- r over
    // 0.486). The cut is carried on
    // into the face beyond Eob (Fob): on Fv's surface from the arc to P,
    // where it crosses Fob's, and on Fob from P to the line's end on Eob.
    // Fv's surface carried on past Vtx meets Fob from P to R on Eob, and
    // Arcprol is carried on to R; Eob between R and the line's end lies
    // under the fillet. P must lie on Fob, under Fop's surface and off Fv.
    if (!inters && CPopArc.IsOnArc() && FvT.IsNull() && Etan.IsNull() && Arcprolbis.IsNull()
        && !hasVertex(CPopArc.Arc(), Vtx) && myEFMap(CPopArc.Arc()).Extent() == 2)
    {
      Eob = CPopArc.Arc();
      for (NCollection_List<TopoDS_Shape>::Iterator itF(myEFMap(Eob)); itF.More(); itF.Next())
      {
        if (!itF.Value().IsSame(Fop))
        {
          Fob = TopoDS::Face(itF.Value());
        }
      }
      const TopoDS_Face             FopSh = TopoDS::Face(DStr.Shape(Fd->Index(IFopArc)));
      double                        f = 0., l = 0.;
      const occ::handle<Geom_Curve> aCob =
        Fob.IsNull() ? occ::handle<Geom_Curve>() : BRep_Tool::Curve(Eob, f, l);
      occ::handle<Geom_Surface> aSv = HBs->Surface().Surface();
      if (HBs->Trsf().Form() != gp_Identity)
      {
        aSv = occ::down_cast<Geom_Surface>(aSv->Transformed(HBs->Trsf()));
      }
      // R: of Eob's crossings with Fv's surface, the nearest the line's end
      const double aParQ = CPopArc.ParameterOnArc();
      bool         hasR  = false;
      if (!aCob.IsNull())
      {
        GeomAPI_IntCS anInt(aCob, aSv);
        double        aDMin = RealLast();
        for (int i = 1; anInt.IsDone() && i <= anInt.NbPoints(); i++)
        {
          double uS, vS, w;
          anInt.Parameters(i, uS, vS, w);
          if (aCob->IsPeriodic())
          {
            w = ElCLib::InPeriod(w, f, f + aCob->Period());
          }
          const double aD = anInt.Point(i).Distance(CPopArc.Point());
          if (w > f + Precision::PConfusion() && w < l - Precision::PConfusion() && aD < aDMin)
          {
            aDMin = aD;
            parR  = w;
            ROff  = anInt.Point(i);
            pvR.SetCoord(uS, vS);
            hasR = true;
          }
        }
      }
      const occ::handle<Geom_Surface>& aSfil = DStr.Surface(Fd->Surf()).Surface();
      pfilP = FiopArc.PCurveOnSurf()->Value(FiopArc.Parameter(isfirst));
      if (hasR && std::abs(parR - aParQ) > Precision::PConfusion()
          && PointOnThree(aSfil, aSv, BRep_Tool::Surface(Fob), pfilP, POff))
      {
        GeomAPI_ProjectPointOnSurf aPrV(POff, aSv), aPrOb(POff, BRep_Tool::Surface(Fob));
        bool                       isIn = aPrV.NbPoints() > 0 && aPrOb.NbPoints() > 0;
        if (isIn)
        {
          double u, v;
          aPrV.LowerDistanceParameters(u, v);
          pvP.SetCoord(u, v);
          aPrOb.LowerDistanceParameters(u, v);
          pobP.SetCoord(u, v);
          const double aSide  = SignedDistance(POff, BRep_Tool::Surface(FopSh), isIn);
          const double aDepth = FopSh.Orientation() == TopAbs_REVERSED ? aSide : -aSide;
          isIn = isIn && aDepth > 10. * Precision::Confusion()
                 && BRepTopAdaptor_FClass2d(Fob, Precision::PConfusion()).Perform(pobP) == TopAbs_IN
                 && BRepTopAdaptor_FClass2d(Fv, Precision::PConfusion()).Perform(pvP) != TopAbs_IN;
        }
        double                          fp, lp;
        const occ::handle<Geom2d_Curve> aPcOb =
          isIn ? PCurveInFace(Eob, Fob, fp, lp) : occ::handle<Geom2d_Curve>();
        const occ::handle<Geom2d_Curve> aPcOp =
          isIn ? PCurveInFace(Eob, Fop, fp, lp) : occ::handle<Geom2d_Curve>();
        if (!aPcOb.IsNull() && !aPcOp.IsNull())
        {
          pobQ = aPcOb->Value(aParQ);
          pobR = aPcOb->Value(parR);
          popR = aPcOp->Value(parR);
          OFob = Fob.Orientation();
          for (NCollection_List<TopoDS_Shape>::Iterator itF(myEFMap(Arcprol)); itF.More();
               itF.Next())
          {
            if (itF.Value().IsSame(Fv))
            {
              OFvSh = itF.Value().Orientation();
            }
          }
          Fob.Orientation(TopAbs_FORWARD);
          HBob->Initialize(Fob, false);
          offFop = true;
          inters = true;
        }
      }
    }
    if (isCPopMoved && !offFop)
    {
      CPopArc = aCPopIn;
    }

    // No cap: where Fo's line misses the face at end whole, the end flat,
    // if asked for, from the end as it was -- only once everything else
    // failed (ChFi3d_MissedFlatEnd): a fillet refused here can still be made
    // by the fallbacks after it (#273's D-hole, its flat's rim at the bottom,
    // E36: the line on the bottom face passes the hole's cylinder, and the
    // corner setback fallback, the end set back, takes it round into it)
    if (!inters && !theAtFourEdges && stat == ChFiDS_OnSame && ChFi3d_Builder::FlatEndAtVertex()
        && !ChFi3d_MissedFlatEnd())
    {
      ChFi3d_UncappedEnd() = true;
    }
    else if (!inters && !theAtFourEdges && stat == ChFiDS_OnSame
             && ChFi3d_Builder::FlatEndAtVertex())
    {
      const ChFiDS_CommonPoint      aCV1 = Fd->Vertex(isfirst, 1), aCV2 = Fd->Vertex(isfirst, 2);
      const ChFiDS_FaceInterference aFi1 = Fd->Interference(1), aFi2 = Fd->Interference(2);
      Fd->ChangeVertex(isfirst, 1) = aCVIn1;
      Fd->ChangeVertex(isfirst, 2) = aCVIn2;
      Fd->ChangeInterference(1)    = aFiIn1;
      Fd->ChangeInterference(2)    = aFiIn2;
      if (MakeFlatEnd(DStr,
                      stripe,
                      Fd,
                      isfirst,
                      Vtx,
                      Arcspine,
                      myVEMap,
                      myEFMap,
                      myVFMap,
                      tolapp3d,
                      tol2d,
                      myEVIMap,
                      true))
      {
        return;
      }
      Fd->ChangeVertex(isfirst, 1) = aCV1;
      Fd->ChangeVertex(isfirst, 2) = aCV2;
      Fd->ChangeInterference(1)    = aFi1;
      Fd->ChangeInterference(2)    = aFi2;
    }
  }
  else
  {
    inters = FindFace(Vtx, CV1, CV2, Fv, Fop);
    // A side ending on a vertex, by an edge that is not the corner's: its
    // line ran along a split of its wall, the fillet as wide as the piece
    // the spine is on, and the walk, on the piece beyond, ended on that
    // piece's own edge at the vertex. The edge of Vtx to that vertex -- the
    // near piece's, under the fillet whole -- is the corner's: the side
    // ends at its far end, and the face at end is the one between the two
    // edges of Vtx.
    if (!inters && CV1.IsOnArc() && CV2.IsOnArc())
    {
      const ChFiDS_CommonPoint aSaved1 = CV1, aSaved2 = CV2;
      for (int ons = 1; ons <= 2; ons++)
      {
        ChFiDS_CommonPoint&      aCP    = Fd->ChangeVertex(isfirst, ons);
        const ChFiDS_CommonPoint aSaved = aCP;
        if (!aCP.IsVertex() || hasVertex(aCP.Arc(), Vtx))
        {
          continue;
        }
        for (NCollection_List<TopoDS_Shape>::Iterator itE(myVEMap(Vtx)); itE.More(); itE.Next())
        {
          const TopoDS_Edge& anE = TopoDS::Edge(itE.Value());
          if (anE.IsSame(Arcspine) || BRep_Tool::Degenerated(anE) || myEFMap(anE).Extent() != 2
              || !hasVertex(anE, aSaved.Vertex()))
          {
            continue;
          }
          // nothing of the edge is kept: cut from its far end on
          TopAbs_Orientation aOV = TopAbs_FORWARD;
          for (ex.Init(anE.Oriented(TopAbs_FORWARD), TopAbs_VERTEX); ex.More(); ex.Next())
          {
            if (aSaved.Vertex().IsSame(ex.Current()))
            {
              aOV = ex.Current().Orientation();
              break;
            }
          }
          aCP.SetArc(aSaved.Tolerance(),
                     anE,
                     BRep_Tool::Parameter(aSaved.Vertex(), anE),
                     TopAbs::Reverse(aOV));
          break;
        }
      }
      inters = FindFace(Vtx, CV1, CV2, Fv, Fop);
      if (!inters)
      {
        CV1 = aSaved1;
        CV2 = aSaved2;
      }
    }
    // The walk stopped at the spine's end with one side on an arc of Vtx
    // and the other in its face, a plane: the line there, carried on,
    // crosses into a coplanar piece of the same wall and meets the face at
    // end beyond (a wall kept in pieces, the radius wider than the piece's
    // edge at Vtx). The line ends on that edge of the piece, as it would on
    // the wall in one face, and is carried over the split below.
    // Or that side stopped on the split itself, which ends at Vtx: the
    // piece the spine is on ends with the spine, and the line stands at the
    // split as one stopped in its face a hair before it.
    int onsOnSplit = 0;
    for (int ons = 1; ons <= 2 && !inters && CV1.IsOnArc() && CV2.IsOnArc(); ons++)
    {
      const TopoDS_Edge& anArc = Fd->Vertex(isfirst, ons).Arc();
      const TopoDS_Edge& anOpp = Fd->Vertex(isfirst, 3 - ons).Arc();
      const TopoDS_Face  aF    = TopoDS::Face(DStr.Shape(Fd->Index(ons)));
      if (!hasVertex(anArc, Vtx) || !hasVertex(anOpp, Vtx) || myEFMap(anArc).Extent() != 2
          || Fd->Vertex(isfirst, ons).IsVertex())
      {
        continue;
      }
      TopoDS_Face aFn;
      for (NCollection_List<TopoDS_Shape>::Iterator itF(myEFMap(anArc)); itF.More(); itF.Next())
      {
        if (!aF.IsSame(itF.Value()))
        {
          aFn = TopoDS::Face(itF.Value());
        }
      }
      if (!aFn.IsNull() && containE(aF, anArc) && ChFi3d_CoplanarPieces(aF, aFn, 10 * tolapp3d))
      {
        onsOnSplit = onsOnSplit == 0 ? ons : -1;
      }
    }
    if (!inters && (CV1.IsOnArc() != CV2.IsOnArc() || onsOnSplit > 0))
    {
      const int                 onsArc = onsOnSplit > 0 ? 3 - onsOnSplit : CV1.IsOnArc() ? 1 : 2;
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
    // The face at end kept in two tangent pieces, the fillet wider than the
    // one at Vtx: one side ends on an edge of Vtx, the other on an edge of
    // the piece beyond (FvT), past the vertex Vp where the edge between the
    // two pieces (Etg, straight) starts. The cut runs over both pieces, as
    // over a wall in two faces at an end OnSame, but crosses Etg itself:
    // Etg is cut there, and what lies between Vp and the crossing, with the
    // edge of Fv from Vp to Vtx, is under the fillet.
    if (!inters && CV1.IsOnArc() && CV2.IsOnArc())
    {
      for (int ad = 1; ad <= 2 && FvT.IsNull(); ad++)
      {
        const TopoDS_Edge& anArcAd = Fd->Vertex(isfirst, ad).Arc();
        const TopoDS_Edge& anArcOp = Fd->Vertex(isfirst, 3 - ad).Arc();
        if (hasVertex(anArcAd, Vtx) || !hasVertex(anArcOp, Vtx) || anArcOp.IsSame(Arcspine)
            || myEFMap(anArcAd).Extent() != 2 || myEFMap(anArcOp).Extent() != 2)
        {
          continue;
        }
        const TopoDS_Face aFOp = TopoDS::Face(DStr.Shape(Fd->Index(3 - ad)));
        TopoDS_Face       aFV;
        for (NCollection_List<TopoDS_Shape>::Iterator itF(myEFMap(anArcOp)); itF.More(); itF.Next())
        {
          if (!aFOp.IsSame(itF.Value()))
          {
            aFV = TopoDS::Face(itF.Value());
          }
        }
        if (aFV.IsNull() || !containV(aFV, Vtx) || containE(aFV, anArcAd))
        {
          continue;
        }
        for (ex.Init(anArcAd, TopAbs_VERTEX); ex.More() && FvT.IsNull(); ex.Next())
        {
          const TopoDS_Vertex& aVp = TopoDS::Vertex(ex.Current());
          TopoDS_Edge          anUnder;
          for (NCollection_List<TopoDS_Shape>::Iterator itE(myVEMap(aVp)); itE.More(); itE.Next())
          {
            const TopoDS_Edge& anE = TopoDS::Edge(itE.Value());
            if (!anE.IsSame(anArcAd) && !anE.IsSame(Arcspine) && !anE.IsSame(anArcOp)
                && hasVertex(anE, Vtx) && containE(aFV, anE))
            {
              anUnder = anE;
            }
          }
          if (anUnder.IsNull())
          {
            continue;
          }
          for (NCollection_List<TopoDS_Shape>::Iterator itF(myEFMap(anArcAd));
               itF.More() && FvT.IsNull();
               itF.Next())
          {
            const TopoDS_Face& aFT  = TopoDS::Face(itF.Value());
            TopAbs_Orientation anOT = TopAbs_FORWARD, anOV = TopAbs_FORWARD;
            const TopoDS_Edge  anEtg = TangentNeighbour(aVp, aFT, aFV, myVEMap, myEFMap, anOT);
            gp_Lin             aLin;
            if (anEtg.IsNull() || anEtg.IsSame(anUnder) || !StraightFromVertex(anEtg, aVp, aLin))
            {
              continue;
            }
            TangentNeighbour(aVp, aFV, aFT, myVEMap, myEFMap, anOV);
            TopoDS_Vertex aV1, aV2;
            TopExp::Vertices(anEtg, aV1, aV2);
            FvT     = aFT;
            Fv      = aFV;
            Etg     = anEtg;
            Eunder  = anUnder;
            Vp      = aVp;
            LinTg   = aLin.Reversed();
            OFvT    = anOT;
            OFv     = anOV;
            cutEtg  = true;
            lenEtg  = BRep_Tool::Pnt(aV1).Distance(BRep_Tool::Pnt(aV2));
            IFadArc = ad;
            IFopArc = 3 - ad;
          }
        }
      }
      if (cutEtg)
      {
        FvT.Orientation(TopAbs_FORWARD);
        Fv.Orientation(TopAbs_FORWARD);
        Swallowed.Append(Eunder);
        Swallowed.Append(Vtx);
        inters = true;
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
    for (int ons = 1; ons <= 2 && !cutEtg; ons++)
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
    if (cutEtg)
    {
      HBsT->Initialize(FvT);
    }
    occ::handle<BRepAdaptor_Curve2d> pced = new BRepAdaptor_Curve2d();
    // (each on the piece its arc bounds)
    const bool onT1 = cutEtg && IFadArc == 1, onT2 = cutEtg && IFadArc == 2;
    pced->Initialize(CV1.Arc(), onT1 ? FvT : Fv);
    Update(onT1 ? HBsT : HBs, pced, HGs, Fd->ChangeInterferenceOnS1(), CV1, isfirst);
    pced->Initialize(CV2.Arc(), onT2 ? FvT : Fv);
    Update(onT2 ? HBsT : HBs, pced, HGs, Fd->ChangeInterferenceOnS2(), CV2, isfirst);
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
  // The cut off Fop (offFop): its piece on Fob (Cc is the piece on Fv), Fv's
  // surface on Fob from P to R, the tolerances they reached, and the arc's
  // point on Fv
  occ::handle<Geom_Curve>   CcOb, C3;
  occ::handle<Geom2d_Curve> PsOb, PcOb, Pc3v, Pc3ob;
  double                    tolOb = tolapp3d, tol3 = tolapp3d;
  gp_Pnt2d                  pvA;

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
      Hc1 = PCurveInFace(CV1.Arc(), FvT.IsNull() || IFadArc != 1 ? Fv : FvT, Ubid, Ubid);
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
      Hc2 = PCurveInFace(CV2.Arc(), FvT.IsNull() || IFadArc != 2 ? Fv : FvT, Ubid, Ubid);
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
    if (offFop)
    {
      // The cut off Fop: on Fv from the arc to P and on Fob from P to the
      // line's end, each the way of the cut from CV1's side to CV2's; Fv's
      // surface on Fob from P to R.
      const gp_Pnt2d& pfilA = IFadArc == 1 ? pfil1 : pfil2;
      const gp_Pnt2d& pfilQ = IFopArc == 1 ? pfil1 : pfil2;
      pvA                   = IFadArc == 1 ? pfac1 : pfac2;
      auto aCurve           = [&](const occ::handle<Adaptor3d_Surface>& theS1,
                        const occ::handle<Adaptor3d_Surface>& theS2,
                        const gp_Pnt2d&                       theA1,
                        const gp_Pnt2d&                       theA2,
                        const gp_Pnt2d&                       theB1,
                        const gp_Pnt2d&                       theB2,
                        occ::handle<Geom_Curve>&              theC,
                        occ::handle<Geom2d_Curve>&            thePc1,
                        occ::handle<Geom2d_Curve>&            thePc2,
                        double&                               theTol) {
        Pardeb(1) = theA1.X();
        Pardeb(2) = theA1.Y();
        Pardeb(3) = theA2.X();
        Pardeb(4) = theA2.Y();
        Parfin(1) = theB1.X();
        Parfin(2) = theB1.Y();
        Parfin(3) = theB2.X();
        Parfin(4) = theB2.Y();
        theTol    = tolapp3d;
        return ChFi3d_ComputeCurves(theS1,
                                    theS2,
                                    Pardeb,
                                    Parfin,
                                    theC,
                                    thePc1,
                                    thePc2,
                                    tolapp3d,
                                    tol2d,
                                    theTol);
      };
      double uu1, uu2, vv1, vv2;
      ChFi3d_Boite(pvA, pvP, uu1, uu2, vv1, vv2);
      ChFi3d_BoundFac(Bs, uu1, uu2, vv1, vv2);
      bool isOk = IFadArc == 1 ? aCurve(HGs, HBs, pfilA, pvA, pfilP, pvP, Cc, Ps, Pc, tolV)
                               : aCurve(HGs, HBs, pfilP, pvP, pfilA, pvA, Cc, Ps, Pc, tolV);
      ChFi3d_Boite(pobQ, pobP, uu1, uu2, vv1, vv2);
      ChFi3d_BoundFac(*HBob, uu1, uu2, vv1, vv2);
      const bool isOkV = isOk;
      isOk =
        isOk
        && (IFopArc == 1 ? aCurve(HGs, HBob, pfilQ, pobQ, pfilP, pobP, CcOb, PsOb, PcOb, tolOb)
                         : aCurve(HGs, HBob, pfilP, pobP, pfilQ, pobQ, CcOb, PsOb, PcOb, tolOb));
      // (Fob a round across Eob, tangent to Fop there as the fillet is: the
      // two touch at the line's end Q, and the intersection finds no line --
      // #273's Fillet001, its handle's top rounded at the end. Then the
      // curve through their points from Q to P.)
      if (isOkV && !isOk)
      {
        const occ::handle<Geom_Surface>& aSfil = DStr.Surface(Fd->Surf()).Surface();
        isOk = CurveOnTwoByPoints(aSfil,
                                  BRep_Tool::Surface(Fob),
                                  IFopArc == 1 ? pfilQ : pfilP,
                                  IFopArc == 1 ? pfilP : pfilQ,
                                  10. * tolapp3d,
                                  CcOb,
                                  PsOb,
                                  PcOb,
                                  tolOb);
        // (failed, it is none: the segment below is tried only after it)
        if (!isOk)
        {
          CcOb.Nullify();
        }
      }
      ChFi3d_Boite(pvP, pvR, uu1, uu2, vv1, vv2);
      ChFi3d_BoundFac(Bs, uu1, uu2, vv1, vv2);
      ChFi3d_Boite(pobP, pobR, uu1, uu2, vv1, vv2);
      ChFi3d_BoundFac(*HBob, uu1, uu2, vv1, vv2);
      isOk = isOk && aCurve(HBs, HBob, pvP, pobP, pvR, pobR, C3, Pc3v, Pc3ob, tol3);
      // (short, the cut's end just past Fop's edge, it can fail: then the
      // segment from P to R, on the two surfaces within the tolerance it
      // reaches)
      if (!isOk && !Cc.IsNull() && !CcOb.IsNull() && POff.Distance(ROff) > Precision::Confusion())
      {
        occ::handle<Geom_Surface> aSv = HBs->Surface().Surface();
        if (HBs->Trsf().Form() != gp_Identity)
        {
          aSv = occ::down_cast<Geom_Surface>(aSv->Transformed(HBs->Trsf()));
        }
        const double aLen = POff.Distance(ROff);
        C3    = new Geom_TrimmedCurve(new Geom_Line(POff, gp_Vec(POff, ROff)), 0., aLen);
        Pc3v  = GeomProjLib::Curve2d(C3, 0., aLen, aSv);
        Pc3ob = GeomProjLib::Curve2d(C3, 0., aLen, BRep_Tool::Surface(Fob));
        isOk  = !Pc3v.IsNull() && !Pc3ob.IsNull();
        if (isOk)
        {
          tol3 = ChFi3d_EvalTolReached(HBs, Pc3v, HBob, Pc3ob, C3);
          isOk = tol3 <= 10. * tolapp3d;
        }
      }
      if (!isOk)
      {
        throw Standard_Failure("OneCorner : the cut off the face fails");
      }
      tolreached = std::max(tolV, tolOb);
      Udeb       = Cc->FirstParameter();
      Ufin       = Cc->LastParameter();
      couture    = false;
      aCutOnOp   = pvR;
    }
    else if (!FvT.IsNull())
    {
      // A face's surface where the face lies: the adaptor of a face that
      // carries a location holds the surface without it.
      auto aSurfOf = [](const occ::handle<BRepAdaptor_Surface>& theHB) {
        const occ::handle<Geom_Surface>& aS = theHB->Surface().Surface();
        return theHB->Trsf().Form() == gp_Identity
                 ? aS
                 : occ::down_cast<Geom_Surface>(aS->Transformed(theHB->Trsf()));
      };
      // The point where the cut crosses from FvT to Fv: on Etg's line
      // carried on past Vp (or on Etg, at an end not OnSame), the nearest
      // of its intersections with the fillet's surface within the
      // surface's bounds.
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
        if (w < -aTolVp || u < aU1 - aMargU || u > aU2 + aMargU || v < aV1 - aMargV
            || v > aV2 + aMargV || (isFound && w >= wTg)
            || (cutEtg && w >= lenEtg - Precision::Confusion()))
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
      atVp = wTg <= (cutEtg ? Precision::Confusion() : aTolVp);
      gp_Pnt2d pfacT, pfacV;
      if (atVp)
      {
        // the point on each face's surface, on the period of Vp there
        for (int k = 0; k < 2; k++)
        {
          occ::handle<Geom_Surface> aS = aSurfOf(k == 0 ? HBsT : HBs);
          const occ::handle<Geom_RectangularTrimmedSurface> aTr =
            occ::down_cast<Geom_RectangularTrimmedSurface>(aS);
          if (!aTr.IsNull())
          {
            aS = aTr->BasisSurface();
          }
          GeomAPI_ProjectPointOnSurf aProj(PTg, aS);
          if (aProj.NbPoints() == 0 || aProj.LowerDistance() > 10 * tolapp3d)
          {
            throw Standard_Failure("OneCorner : the wall's tangent edge does not carry on");
          }
          double u, v;
          aProj.LowerDistanceParameters(u, v);
          const gp_Pnt2d aUVv = BRep_Tool::Parameters(Vp, k == 0 ? FvT : Fv);
          if (aS->IsUPeriodic())
          {
            u += aS->UPeriod() * std::floor((aUVv.X() - u) / aS->UPeriod() + 0.5);
          }
          if (aS->IsVPeriodic())
          {
            v += aS->VPeriod() * std::floor((aUVv.Y() - v) / aS->VPeriod() + 0.5);
          }
          (k == 0 ? pfacT : pfacV).SetCoord(u, v);
        }
      }
      else
      {
        CTg  = new Geom_TrimmedCurve(new Geom_Line(LinTg), 0., wTg);
        PTgT = CurveOnFacePeriod(CTg, aSurfOf(HBsT), FvT, Vp, 10 * tolapp3d);
        PTgV = CurveOnFacePeriod(CTg, aSurfOf(HBs), Fv, Vp, 10 * tolapp3d);
        if (PTgT.IsNull() || PTgV.IsNull())
        {
          throw Standard_Failure("OneCorner : the wall's tangent edge does not carry on");
        }
        pfacT = PTgT->Value(wTg);
        pfacV = PTgV->Value(wTg);
      }
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
        if (!isBeside && !cutEtg)
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
                                         aSurfOf(HBs));
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
      // The same for the piece on FvT, when the fillet's line on Fop's side
      // ends on FvT's surface as well (the two faces pieces of one plane).
      gp_Pnt2d pfacTOp;
      bool     hasTOp = false;
      {
        const gp_Pnt2d&            pfilOp = IFopArc == 1 ? pfil1 : pfil2;
        GeomAPI_ProjectPointOnSurf aProj(HGs->Value(pfilOp.X(), pfilOp.Y()),
                                         aSurfOf(HBsT));
        if (aProj.NbPoints() > 0 && aProj.LowerDistance() <= tolapp3d)
        {
          double u, v;
          aProj.LowerDistanceParameters(u, v);
          pfacTOp.SetCoord(u, v);
          ChFi3d_Recale(*HBsT, pfacT, pfacTOp, true);
          hasTOp = true;
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
        const gp_Pnt2d&      pfacTV =
          theWhole ? (onT ? pfacTOp : pfacVAd) : (onT ? pfacT : pfacV);
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
        // Cut at the point: the whole curve, and a piece the intersection
        // ended elsewhere (it takes a line's own end for the one asked when
        // the two are close for the piece's length)
        const double aParEnd = k == 2 ? theC->FirstParameter() : theC->LastParameter();
        const double aGap    = theC->Value(aParEnd).Distance(PTg);
        if (theWhole || aGap > Precision::Confusion())
        {
          GeomAPI_ProjectPointOnCurve aProjC(PTg, theC);
          const bool isOff = aProjC.NbPoints() == 0 || aProjC.LowerDistance() > 10 * tolapp3d;
          const double aPar = isOff ? aParEnd : aProjC.LowerDistanceParameter();
          if (isOff || aPar <= theC->FirstParameter() + Precision::PConfusion()
              || aPar >= theC->LastParameter() - Precision::PConfusion())
          {
            // a piece that stops short of the point is as far off as that
            theTol = std::max(theTol, aGap);
            return !theWhole && aGap <= 10 * tolapp3d;
          }
          theC = k == 2 ? new Geom_TrimmedCurve(theC, aPar, theC->LastParameter())
                        : new Geom_TrimmedCurve(theC, theC->FirstParameter(), aPar);
          // the piece kept's own: past the point the curve runs on off the
          // parts of the surfaces it is computed on
          theTol = std::max(ChFi3d_EvalTolReached(HGs, thePs, onT ? HBsT : HBs, thePc, theC),
                            aProjC.LowerDistance());
        }
        return true;
      };
      double aTolT = 0., aTolV = 0.;
      // (two pieces of one surface: the whole curve is as good as a piece)
      const double aTolTry = hasTOp ? 10 * Precision::Confusion() : tolapp3d;
      const bool   isCutT  = aPiece(IFadArc, false, CcT, PsT, PcT, aTolT);
      if ((!isCutT || aTolT > aTolTry) && hasTOp)
      {
        occ::handle<Geom_Curve>   aC;
        occ::handle<Geom2d_Curve> aPs, aPc;
        double                    aTolW = 0.;
        if (aPiece(IFadArc, true, aC, aPs, aPc, aTolW) && (!isCutT || aTolW < aTolT))
        {
          CcT   = aC;
          PsT   = aPs;
          PcT   = aPc;
          aTolT = aTolW;
        }
        else if (!isCutT)
        {
          throw Standard_Failure("OneCorner : echec calcul intersection");
        }
      }
      else if (!isCutT)
      {
        throw Standard_Failure("OneCorner : echec calcul intersection");
      }
      const bool isCut = aPiece(IFopArc, false, Cc, Ps, Pc, aTolV);
      if (!isCut || aTolV > aTolTry)
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

    // a curve ending tangent to its face's boundary kept on the face (the
    // cut off Fop ends off Fv, on its surface carried on)
    const double aKeptV = offFop ? 0. : KeepCurveOnFace(Pc, Fv, Udeb, Ufin);
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
  // the cut off Fop's: P and R, and their indices
  Bnd_Box boxP, boxR;
  int     indP = 0, indR = 0;
  if (offFop)
  {
    // The cut off Fop, in two curves meeting at P, as the cut over two faces
    // below; Fv's surface on Fob from P to R; and Eob cut at R and at the
    // line's end, between them under the fillet. On the fillet each piece of
    // the cut goes the way the plain cut does, and the piece on Fv too (Et).
    indP            = DStr.AddPoint(TopOpeBRepDS_Point(POff, tolreached));
    indR            = DStr.AddPoint(TopOpeBRepDS_Point(ROff, tol3));
    const int ICV   = DStr.AddCurve(TopOpeBRepDS_Curve(Cc, tolV));
    const int ICOb  = DStr.AddCurve(TopOpeBRepDS_Curve(CcOb, tolOb));
    const int IC3   = DStr.AddCurve(TopOpeBRepDS_Curve(C3, tol3));
    const int IShOb = DStr.AddShape(Fob);
    // At the cut, the cap -- the end of the groove the fillet leaves -- lies
    // on Fv the way of the fillet's outward normal, what is kept of Fob the
    // other way (the notch in it is the groove's end). Which way the
    // fillet's normal runs, the piece on Fv says (Et). On Fv's surface on
    // Fob, what is kept of Fv lies within Fob, what is kept of Fob outside
    // Fv (the notch beyond).
    const occ::handle<Geom_Surface>& aSfil = DStr.Surface(Isurf).Surface();
    const occ::handle<Geom_Surface>  aSv   = BRep_Tool::Surface(Fv);
    const occ::handle<Geom_Surface>  aSob  = BRep_Tool::Surface(Fob);
    const double                     wV    = 0.5 * (Udeb + Ufin);
    const bool                       isOutN =
      IsKeptLeft(aSv, Pc, Cc, NormalAt(aSfil, Ps->Value(wV))) == (Et == TopAbs_FORWARD);
    const double             wOb  = 0.5 * (CcOb->FirstParameter() + CcOb->LastParameter());
    const gp_Vec             aNOb = NormalAt(aSfil, PsOb->Value(wOb));
    const TopAbs_Orientation EtOb =
      IsKeptLeft(aSob, PcOb, CcOb, isOutN ? -aNOb : aNOb) ? TopAbs_FORWARD : TopAbs_REVERSED;
    const double w3    = 0.5 * (C3->FirstParameter() + C3->LastParameter());
    const gp_Vec aNOb3 = NormalAt(aSob, Pc3ob->Value(w3)) * (OFob == TopAbs_REVERSED ? -1. : 1.);
    const gp_Vec aNV3  = NormalAt(aSv, Pc3v->Value(w3)) * (OFvSh == TopAbs_REVERSED ? -1. : 1.);
    const TopAbs_Orientation Et3V =
      IsKeptLeft(aSv, Pc3v, C3, -aNOb3) ? TopAbs_FORWARD : TopAbs_REVERSED;
    const TopAbs_Orientation Et3Ob =
      IsKeptLeft(aSob, Pc3ob, C3, -aNV3) ? TopAbs_FORWARD : TopAbs_REVERSED;
    DStr.ChangeShapeInterferences(IShape).Append(ChFi3d_FilCurveInDS(ICV, IShape, Pc, Et));
    DStr.ChangeShapeInterferences(IShOb).Append(ChFi3d_FilCurveInDS(ICOb, IShOb, PcOb, EtOb));
    TopAbs_Orientation aTrafil1 = TopAbs_FORWARD;
    if (Fd->IndexOfS1() > 0)
    {
      aTrafil1 = DStr.Shape(Fd->IndexOfS1()).Orientation();
    }
    aTrafil1 = TopAbs::Compose(aTrafil1, Fd->Orientation());
    aTrafil1 = TopAbs::Compose(TopAbs::Reverse(Fd->InterferenceOnS1().Transition()), aTrafil1);
    const TopAbs_Orientation EtS = isfirst ? TopAbs::Reverse(aTrafil1) : aTrafil1;
    DStr.ChangeSurfaceInterferences(Isurf).Append(ChFi3d_FilCurveInDS(ICOb, Isurf, PsOb, EtS));
    DStr.ChangeSurfaceInterferences(Isurf).Append(ChFi3d_FilCurveInDS(ICV, Isurf, Ps, EtS));
    stripe->InDS(isfirst);
    const int                      ind1 = stripe->IndexPoint(isfirst, 1);
    const int                      ind2 = stripe->IndexPoint(isfirst, 2);
    const int                      IC1  = IFadArc == 1 ? ICV : ICOb;
    const int                      IC2  = IFadArc == 1 ? ICOb : ICV;
    const occ::handle<Geom_Curve>& C1   = IFadArc == 1 ? Cc : CcOb;
    const occ::handle<Geom_Curve>& C2   = IFadArc == 1 ? CcOb : Cc;
    DStr.ChangeCurveInterferences(IC1).Append(
      ChFi3d_FilPointInDS(TopAbs_FORWARD, IC1, ind1, C1->FirstParameter()));
    DStr.ChangeCurveInterferences(IC1).Append(
      ChFi3d_FilPointInDS(TopAbs_REVERSED, IC1, indP, C1->LastParameter()));
    DStr.ChangeCurveInterferences(IC2).Append(
      ChFi3d_FilPointInDS(TopAbs_FORWARD, IC2, indP, C2->FirstParameter()));
    DStr.ChangeCurveInterferences(IC2).Append(
      ChFi3d_FilPointInDS(TopAbs_REVERSED, IC2, ind2, C2->LastParameter()));

    const double f3 = C3->FirstParameter(), l3 = C3->LastParameter();
    DStr.ChangeShapeInterferences(IShape).Append(ChFi3d_FilCurveInDS(IC3, IShape, Pc3v, Et3V));
    DStr.ChangeShapeInterferences(IShOb).Append(ChFi3d_FilCurveInDS(IC3, IShOb, Pc3ob, Et3Ob));
    DStr.ChangeCurveInterferences(IC3).Append(ChFi3d_FilPointInDS(TopAbs_FORWARD, IC3, indP, f3));
    DStr.ChangeCurveInterferences(IC3).Append(ChFi3d_FilPointInDS(TopAbs_REVERSED, IC3, indR, l3));

    // Eob ends at R and starts again at the line's end, or the other way
    const int    IEob = DStr.AddShape(Eob);
    const int    indQ = stripe->IndexPoint(isfirst, IFopArc);
    const double parQ = Fd->Vertex(isfirst, IFopArc).ParameterOnArc();
    DStr.ChangeShapeInterferences(IEob).Append(
      ChFi3d_FilPointInDS(parQ > parR ? TopAbs_REVERSED : TopAbs_FORWARD, IEob, indR, parR));
    DStr.ChangeShapeInterferences(IEob).Append(
      ChFi3d_FilPointInDS(parQ > parR ? TopAbs_FORWARD : TopAbs_REVERSED, IEob, indQ, parQ));

    Bnd_Box& boxA = IFadArc == 1 ? box1 : box2;
    Bnd_Box& boxQ = IFopArc == 1 ? box1 : box2;
    ChFi3d_EnlargeBox(HBs, Pc, Udeb, Ufin, IFadArc == 1 ? boxA : boxP, IFadArc == 1 ? boxP : boxA);
    ChFi3d_EnlargeBox(HBob,
                      PcOb,
                      CcOb->FirstParameter(),
                      CcOb->LastParameter(),
                      IFopArc == 1 ? boxQ : boxP,
                      IFopArc == 1 ? boxP : boxQ);
    ChFi3d_EnlargeBox(HBs, Pc3v, f3, l3, boxP, boxR);
    ChFi3d_EnlargeBox(HBob, Pc3ob, f3, l3, boxP, boxR);
    ChFi3d_EnlargeBox(Eob, myEFMap(Eob), parR, boxR);
    ChFi3d_EnlargeBox(Eob, myEFMap(Eob), parQ, boxQ);
  }
  else if (!FvT.IsNull())
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
    const bool noTgEdge = cutEtg || atVp;
    if (noTgEdge)
    {
      // Etg itself is crossed, or its line within reach of Vp: it starts
      // at the point, on Vp's side
      double parE = BRep_Tool::Parameter(Vp, Etg);
      if (!atVp)
      {
        double                        f, l;
        const occ::handle<Geom_Curve> aCE = BRep_Tool::Curve(Etg, f, l);
        GeomAPI_ProjectPointOnCurve   aProjE(PTg, aCE, f, l);
        if (aProjE.NbPoints() == 0 || aProjE.LowerDistance() > 10 * tolapp3d)
        {
          throw Standard_Failure(
            "OneCorner : the cut's crossing is off the edge between the pieces");
        }
        parE = aProjE.LowerDistanceParameter();
      }
      else
      {
        // (a straight edge is carried on to the point itself: its end left
        // at Vp, the faces' wires are open there by the distance between)
        double                  f, l;
        occ::handle<Geom_Curve> aCE = BRep_Tool::Curve(Etg, f, l);
        if (!aCE.IsNull() && aCE->IsKind(STANDARD_TYPE(Geom_TrimmedCurve)))
        {
          aCE = occ::down_cast<Geom_TrimmedCurve>(aCE)->BasisCurve();
        }
        const occ::handle<Geom_Line> aLE = occ::down_cast<Geom_Line>(aCE);
        if (!aLE.IsNull())
        {
          parE = ElCLib::Parameter(aLE->Lin(), PTg);
        }
      }
      const int IEtg = DStr.AddShape(Etg);
      DStr.ChangeShapeInterferences(IEtg).Append(
        ChFi3d_FilPointInDS(sameDir ? TopAbs_REVERSED : TopAbs_FORWARD, IEtg, indpt, parE));
      ChFi3d_EnlargeBox(Etg, myEFMap(Etg), parE, boxT);
    }
    const double tolG = noTgEdge ? 0.
                                 : std::max(ChFi3d_EvalTolReached(HBsT, PTgT, HBs, PTgV, CTg),
                                            Precision::Confusion());
    const int    IG   = noTgEdge ? 0 : DStr.AddCurve(TopOpeBRepDS_Curve(CTg, tolG));
    for (int k = 0; k < 2 && !noTgEdge; k++)
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
    if (!noTgEdge)
    {
      DStr.ChangeCurveInterferences(IG).Append(
        ChFi3d_FilVertexInDS(TopAbs_FORWARD, IG, DStr.AddShape(Vp), 0.));
      DStr.ChangeCurveInterferences(IG).Append(
        ChFi3d_FilPointInDS(TopAbs_REVERSED, IG, indpt, wTg));
    }
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

  if (FvT.IsNull() && onsOS == 0 && !offFop)
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
    Hc = PCurveInFace(Arcprolbis.IsNull() ? Arcprol : Arcprolbis, Fop, Ubid, Ubid);
    if (Hc.IsNull())
    {
      throw Standard_ConstructionError("Failed to get p-curve of edge");
    }
    pop1 = Hc->Value(Arcprolbis.IsNull() ? parVtx : BRep_Tool::Parameter(Vtx, Arcprolbis));
    // (to R, for the cut off Fop)
    pop2 = offFop ? popR : Fiop.PCurveOnFace()->Value(Fiop.Parameter(isfirst));
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
    // The extension over a split: Arcprol bounds a piece of Fop's plane
    // beside Fop (a wall kept in coplanar pieces, the fillet wider than the
    // piece beside its spine: its line lies on the next piece), and the
    // extension, starting in that piece, crosses the edge between the two
    // on its way to the line's end. It is stored in two, each part on its
    // piece, and the edge crossed is cut there: what lies on the spine's
    // side is under the fillet. All of it on Fop left that face with an open
    // wire, and the piece beside the spine without its end.
    TopoDS_Edge       EsplitZ;
    double            wQZ = 0., parQZ = 0.;
    bool              atEndZ = false;
    const TopoDS_Face FopSh  = TopoDS::Face(DStr.Shape(Fd->Index(IFopArc)));
    if (!zobOnEtan && Etan.IsNull() && !Fprol.IsNull() && !Fprol.IsSame(Fop)
        && ChFi3d_CoplanarPieces(FopSh, Fprol, 10 * tolapp3d)
        && !ExtensionOverSplit(zob3d,
                               Udeb,
                               Ufin,
                               Fprol,
                               Fop,
                               Arcprol,
                               Vtx,
                               10 * tolapp3d,
                               tolapp3d,
                               EsplitZ,
                               wQZ,
                               parQZ,
                               atEndZ))
    {
      EsplitZ.Nullify();
    }
    // The parts keep the tolerance the extension has, not the one it was
    // asked for: the cut of the corner leaves the extension's end tangent
    // to it (the fillet is tangent to the wall along its line), and between
    // there and the split runs within that tolerance of the first part, an
    // edge it shares no vertex with -- a self-intersection for BOP's check.
    // Where the cut passes the crossing itself within the tolerance, the
    // crossing is the line's end, as when the line runs along the split.
    double tolZ = tolreached;
    if (!EsplitZ.IsNull())
    {
      tolZ = std::min(tolreached,
                      std::max(ChFi3d_EvalTolReached(HBop, zob2dop, HBs, zob2dv, zob3d),
                               Precision::Confusion()));
      if (!atEndZ && !Cc.IsNull())
      {
        GeomAPI_ProjectPointOnCurve aProj(zob3d->Value(wQZ), Cc);
        if (aProj.NbPoints() > 0 && aProj.LowerDistance() <= 4. * tolZ)
        {
          atEndZ = true;
          wQZ    = Ufin;
        }
      }
    }
    TopOpeBRepDS_Curve Zob(zob3d, tolreached);
    int                IZob = EsplitZ.IsNull() ? DStr.AddCurve(Zob) : 0;

    // it is determined if Fop has an edge of sewing
    // it is determined if the curve has an intersection with the edge of sewing

    // TopoDS_Edge edgecouture;
    // bool couture;
    ChFi3d_Couture(Fop, couture, edgecouture);
    if (offFop && couture && !BRep_Tool::Degenerated(edgecouture))
    {
      throw Standard_Failure("OneCorner : the cut off the face crosses a seam");
    }

    if (EsplitZ.IsNull() && couture && !BRep_Tool::Degenerated(edgecouture))
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
    if (!EsplitZ.IsNull())
    {
      // the two parts, from Vtx to Q on the split and from Q to the line's
      // end, both on Fv; with the line along the split (atEndZ), Q is the
      // line's end and there is no second part
      const int    icc   = stripe->IndexPoint(isfirst, IFopArc);
      const gp_Pnt aPQ   = zob3d->Value(wQZ);
      const gp_Pnt aPQe  = BRepAdaptor_Curve(EsplitZ).Value(parQZ);
      const double aTolQ = std::max(tolZ, aPQ.Distance(aPQe));
      const int    indQ  = atEndZ ? icc : DStr.AddPoint(TopOpeBRepDS_Point(aPQ, aTolQ));
      const int    IZ1   = DStr.AddCurve(
        TopOpeBRepDS_Curve(new Geom_TrimmedCurve(zob3d, Udeb, wQZ), tolZ));
      Et = TopAbs::Reverse(TopAbs::Compose(OVtx, OArcprolv));
      DStr.ChangeShapeInterferences(IShape).Append(
        ChFi3d_FilCurveInDS(IZ1, IShape, new Geom2d_TrimmedCurve(zob2dv, Udeb, wQZ), Et));
      // the first on Arcprol's piece, as Arcprol lies in it; the second on
      // Fop, the same way round or, Fop the other way up on its surface,
      // the other
      TopoDS_Face aFprolF = Fprol;
      aFprolF.Orientation(TopAbs_FORWARD);
      TopAbs_Orientation OArcprolp = TopAbs_FORWARD;
      for (ex.Init(aFprolF, TopAbs_EDGE); ex.More(); ex.Next())
      {
        if (Arcprol.IsSame(ex.Current()))
        {
          OArcprolp = ex.Current().Orientation();
          break;
        }
      }
      const occ::handle<Geom2d_Curve> aPcProl =
        GeomProjLib::Curve2d(zob3d, Udeb, wQZ, BRep_Tool::Surface(aFprolF));
      if (aPcProl.IsNull())
      {
        throw Standard_Failure("OneCorner : the extension over a split has no pcurve");
      }
      const int Iprol = DStr.AddShape(Fprol);
      DStr.ChangeShapeInterferences(Iprol).Append(
        ChFi3d_FilCurveInDS(IZ1,
                            Iprol,
                            aPcProl,
                            TopAbs::Reverse(TopAbs::Compose(OVtx, OArcprolp))));
      DStr.ChangeCurveInterferences(IZ1).Append(
        ChFi3d_FilVertexInDS(TopAbs_FORWARD, IZ1, IVtx, Udeb));
      DStr.ChangeCurveInterferences(IZ1).Append(
        ChFi3d_FilPointInDS(TopAbs_REVERSED, IZ1, indQ, wQZ));
      if (!atEndZ)
      {
        const int IZ2 = DStr.AddCurve(
          TopOpeBRepDS_Curve(new Geom_TrimmedCurve(zob3d, wQZ, Ufin), tolZ));
        DStr.ChangeShapeInterferences(IShape).Append(
          ChFi3d_FilCurveInDS(IZ2, IShape, new Geom2d_TrimmedCurve(zob2dv, wQZ, Ufin), Et));
        if (Fprol.Orientation() != FopSh.Orientation())
        {
          OArcprolp = TopAbs::Reverse(OArcprolp);
        }
        const int Iop = DStr.AddShape(Fop);
        DStr.ChangeShapeInterferences(Iop).Append(
          ChFi3d_FilCurveInDS(IZ2,
                              Iop,
                              new Geom2d_TrimmedCurve(zob2dop, wQZ, Ufin),
                              TopAbs::Reverse(TopAbs::Compose(OVtx, OArcprolp))));
        DStr.ChangeCurveInterferences(IZ2).Append(
          ChFi3d_FilPointInDS(TopAbs_FORWARD, IZ2, indQ, wQZ));
        DStr.ChangeCurveInterferences(IZ2).Append(
          ChFi3d_FilPointInDS(TopAbs_REVERSED, IZ2, icc, Ufin));
      }
      // the split is cut at Q: of its two sides there, the one the spine
      // leaves Vtx toward, across the extension, is under the fillet
      BRepAdaptor_Curve aCsp(EsplitZ), aCspine(Arcspine);
      gp_Pnt            aPbid;
      gp_Vec            aDsp, aDz, aDspine;
      aCsp.D1(parQZ, aPbid, aDsp);
      zob3d->D1(wQZ, aPbid, aDz);
      const double aParV = BRep_Tool::Parameter(Vtx, Arcspine);
      aCspine.D1(aParV, aPbid, aDspine);
      if (std::abs(aParV - aCspine.LastParameter()) < std::abs(aParV - aCspine.FirstParameter()))
      {
        aDspine.Reverse();
      }
      const gp_Vec aNz = gp_Vec(BRepAdaptor_Surface(aFprolF, false).Plane().Axis().Direction())
                           .Crossed(aDz);
      const double aSideSp = aDsp.Dot(aNz), aSideSpine = aDspine.Dot(aNz);
      if (std::abs(aSideSp) <= gp::Resolution() || std::abs(aSideSpine) <= gp::Resolution())
      {
        throw Standard_Failure("OneCorner : the extension runs along the split it crosses");
      }
      const int IEs = DStr.AddShape(EsplitZ);
      DStr.ChangeShapeInterferences(IEs).Append(
        ChFi3d_FilPointInDS(aSideSp * aSideSpine > 0. ? TopAbs_REVERSED : TopAbs_FORWARD,
                            IEs,
                            indQ,
                            parQZ));
      (IFopArc == 1 ? box1 : box2).Add(zob3d->Value(Ufin));
      if (atEndZ)
      {
        (IFopArc == 1 ? box1 : box2).Add(aPQe);
      }
    }
    else if (intcouture)
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
      Et      = TopAbs::Reverse(TopAbs::Compose(Arcprolbis.IsNull() ? OVtx : OVtxBis, OArcprolop));
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
      Et = TopAbs::Reverse(TopAbs::Compose(Arcprolbis.IsNull() ? OVtx : OVtxBis, OArcprolop));
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
      int icc    = offFop ? indR : stripe->IndexPoint(isfirst, IFopArc);
      interfprol = ChFi3d_FilPointInDS(TopAbs_REVERSED, IZob, icc, Ufin);
      DStr.ChangeCurveInterferences(IZob).Append(interfprol);
#ifdef VARIANT1
      {
        if (offFop)
        {
          boxR.Add(zob3d->Value(Ufin));
        }
        else if (IFopArc == 1)
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
  if (offFop)
  {
    ChFi3d_SetPointTolerance(DStr, boxP, indP);
    ChFi3d_SetPointTolerance(DStr, boxR, indR);
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
// function : PointOnSurfaceAsItIs
// purpose  : The parameters on <S> of <P>, a point of an edge that lies on
//           the surface and is not an edge of its face <BS> -- no pcurve
//           to read them from: the nearest point of the surface, brought
//           into the face's range on a periodic surface (<BS> restricted
//           to the face). <UV> is left as it is when the projection fails.
//=======================================================================

static void PointOnSurfaceAsItIs(const gp_Pnt&                    P,
                                 const occ::handle<Geom_Surface>& S,
                                 const BRepAdaptor_Surface&       BS,
                                 gp_Pnt2d&                        UV)
{
  if (S.IsNull())
  {
    return;
  }
  GeomAPI_ProjectPointOnSurf aProj(P, S);
  if (!aProj.IsDone() || aProj.NbPoints() == 0)
  {
    return;
  }
  double aU = 0., aV = 0.;
  aProj.LowerDistanceParameters(aU, aV);
  if (S->IsUPeriodic())
  {
    const double aMid = 0.5 * (BS.FirstUParameter() + BS.LastUParameter());
    aU                = ElCLib::InPeriod(aU, aMid - 0.5 * S->UPeriod(), aMid + 0.5 * S->UPeriod());
  }
  if (S->IsVPeriodic())
  {
    const double aMid = 0.5 * (BS.FirstVParameter() + BS.LastVParameter());
    aV                = ElCLib::InPeriod(aV, aMid - 0.5 * S->VPeriod(), aMid + 0.5 * S->VPeriod());
  }
  UV.SetCoord(aU, aV);
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
    if (Fn.IsNull() || Fn.IsSame(Fs) || !ChFi3d_CoplanarPieces(Fs, Fn, Tol))
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

//=======================================================================
// function : ArcPastWallSplit
// purpose  : ArcPastSplit for the face at end in pieces: the side of the
//           fillet on <Fad>, a face of <Vtx>, ends on <Arc>, the edge of
//           <Fad> and the far piece of a wall kept in coplanar pieces, past
//           the foot of the split between that piece and the one at <Vtx>;
//           <Eat> is the edge of <Vtx> between <Fad> and the piece at
//           <Vtx>, which <Arc> carries on.
//=======================================================================

static bool ArcPastWallSplit(const TopoDS_Vertex& Vtx,
                             const TopoDS_Edge&   Arc,
                             const TopoDS_Face&   Fad,
                             const ChFiDS_Map&    VEMap,
                             const ChFiDS_Map&    EFMap,
                             const double         Tol,
                             TopoDS_Edge&         Eat)
{
  if (EFMap(Arc).Extent() != 2 || !containE(Fad, Arc) || !containV(Fad, Vtx))
  {
    return false;
  }
  TopoDS_Face aFw;
  for (NCollection_List<TopoDS_Shape>::Iterator It(EFMap(Arc)); It.More(); It.Next())
  {
    if (!Fad.IsSame(It.Value()))
    {
      aFw = TopoDS::Face(It.Value());
    }
  }
  return !aFw.IsNull() && ArcPastSplit(Vtx, Arc, aFw, VEMap, EFMap, Tol, Eat);
}

//=======================================================================
// function : CeilingAtEnd
// purpose  : The end of a fillet under a wall that stands on the face at
//           end (SectionCrossesEndFace). Side <Ip> of <Fd> ends on an edge
//           of <Vtx> between its face and the face at end, a plane; the
//           other side's line crosses that plane at a point of its own
//           face, which runs on past <Vtx> (or at a point of a split of
//           that face in coplanar pieces, an edge of <Vtx>); and the face
//           at end's other edge at <Vtx>, the wall's, runs between the two.
//           The fillet ends in the plane. Its section there cuts the face
//           at end from the one side's edge to X, where it crosses the
//           wall's edge. From X to the other side it runs under the wall,
//           and what closes the body there is a face in the same plane
//           turned the other way, a ceiling under the wall: between the
//           section, the piece of the wall's edge from <Vtx> to X, which
//           the wall keeps, and the line from the other side's point to
//           <Vtx> on that side's face. All of it is put in the DS here,
//           and <EdgeSpine>, the edge filleted, ended at <Vtx>, which stays.
//           False, and nothing done, where this is not the case.
//=======================================================================

static bool CeilingAtEnd(
  TopOpeBRepDS_DataStructure&                                                        DStr,
  const occ::handle<ChFiDS_Stripe>&                                                  Stripe,
  const occ::handle<ChFiDS_SurfData>&                                                Fd,
  const bool                                                                         isFirst,
  const int                                                                          Ip,
  const TopoDS_Vertex&                                                               Vtx,
  const TopoDS_Edge&                                                                 EdgeSpine,
  const occ::handle<GeomAdaptor_Surface>&                                            HGs,
  const TopAbs_Orientation                                                           OrCurve,
  const ChFiDS_Map&                                                                  VEMap,
  const ChFiDS_Map&                                                                  EFMap,
  const ChFiDS_Map&                                                                  VFMap,
  const double                                                                       Tol3d,
  const double                                                                       Tol2d,
  NCollection_DataMap<TopoDS_Shape, NCollection_List<int>, TopTools_ShapeMapHasher>& EVIMap)
{
  const int                 Io  = 3 - Ip;
  const ChFiDS_CommonPoint& CVp = Fd->Vertex(isFirst, Ip);
  ChFiDS_CommonPoint&       CVo = Fd->ChangeVertex(isFirst, Io);
  if (Fd->Index(1) <= 0 || Fd->Index(2) <= 0 || !CVp.IsOnArc() || CVp.IsVertex()
      || !hasVertex(CVp.Arc(), Vtx) || CVo.IsVertex()
      || DStr.Shape(Fd->Index(1)).ShapeType() != TopAbs_FACE
      || DStr.Shape(Fd->Index(2)).ShapeType() != TopAbs_FACE)
  {
    return false;
  }
  const TopoDS_Face              Fp   = TopoDS::Face(DStr.Shape(Fd->Index(Ip)));
  const TopoDS_Face              Fo   = TopoDS::Face(DStr.Shape(Fd->Index(Io)));
  const TopoDS_Edge              Epiv = CVp.Arc();
  const ChFiDS_FaceInterference& Fip  = Fd->Interference(Ip);
  const ChFiDS_FaceInterference& Fio  = Fd->Interference(Io);
  if (Fip.LineIndex() == 0 || Fio.LineIndex() == 0 || Fip.PCurveOnSurf().IsNull()
      || Fio.PCurveOnSurf().IsNull())
  {
    return false;
  }
  auto anAcross = [&](const TopoDS_Edge& theE, const TopoDS_Face& theF) {
    TopoDS_Face aF;
    if (EFMap(theE).Extent() != 2)
    {
      return aF;
    }
    for (NCollection_List<TopoDS_Shape>::Iterator It(EFMap(theE)); It.More(); It.Next())
    {
      if (!theF.IsSame(It.Value()))
      {
        aF = TopoDS::Face(It.Value());
      }
    }
    return aF;
  };
  auto aPlaneOf = [](const TopoDS_Face& theF) {
    occ::handle<Geom_Surface>                         aS = BRep_Tool::Surface(theF);
    const occ::handle<Geom_RectangularTrimmedSurface> aRT =
      occ::down_cast<Geom_RectangularTrimmedSurface>(aS);
    if (!aRT.IsNull())
    {
      aS = aRT->BasisSurface();
    }
    return occ::down_cast<Geom_Plane>(aS);
  };
  // the edge as it lies in the face, whatever way up the face is
  auto anOrIn = [](const TopoDS_Edge& theE, const TopoDS_Face& theF) {
    for (TopExp_Explorer ex(theF.Oriented(TopAbs_FORWARD), TopAbs_EDGE); ex.More(); ex.Next())
    {
      if (theE.IsSame(ex.Current()))
      {
        return ex.Current().Orientation();
      }
    }
    return theE.Orientation();
  };

  // the face at end, its other edge at Vtx, and the wall on that edge
  const TopoDS_Face Fv = anAcross(Epiv, Fp);
  if (Fv.IsNull())
  {
    return false;
  }
  TopoDS_Edge Arcprol;
  for (NCollection_List<TopoDS_Shape>::Iterator It(VEMap(Vtx)); It.More(); It.Next())
  {
    const TopoDS_Edge& E = TopoDS::Edge(It.Value());
    if (!E.IsSame(Epiv) && !BRep_Tool::Degenerated(E) && containE(Fv, E))
    {
      Arcprol = E;
      break;
    }
  }
  if (Arcprol.IsNull())
  {
    return false;
  }
  const TopoDS_Face Fwall = anAcross(Arcprol, Fv);
  if (Fwall.IsNull() || Fwall.IsSame(Fp) || Fwall.IsSame(Fo))
  {
    return false;
  }
  const occ::handle<Geom_Plane> aGPv = aPlaneOf(Fv);
  if (aGPv.IsNull())
  {
    return false;
  }
  const gp_Pln aPln = aGPv->Pln();
  const gp_Pnt PV   = BRep_Tool::Pnt(Vtx);

  // the other side: where its line crosses the plane
  const bool  onSplit = CVo.IsOnArc() && hasVertex(CVo.Arc(), Vtx);
  TopoDS_Edge Esplit;
  TopoDS_Face Fup;
  double      wo = Fio.Parameter(isFirst);
  gp_Pnt      Po = CVo.Point();
  if (onSplit)
  {
    Esplit = CVo.Arc();
    Fup    = anAcross(Esplit, Fo);
    if (Fup.IsNull() || Fup.IsSame(Fwall) || !ChFi3d_CoplanarPieces(Fo, Fup, 10 * Tol3d))
    {
      return false;
    }
  }
  else
  {
    const occ::handle<Geom_Curve>& aLo = DStr.Curve(Fio.LineIndex()).Curve();
    wo                                 = Fip.Parameter(isFirst);
    if (aLo.IsNull() || wo < Fio.FirstParameter() || wo > Fio.LastParameter())
    {
      return false;
    }
    Po = aLo->Value(wo);
    if (aPln.Distance(Po) > Tol3d)
    {
      // (the plane not square to the spine: next to the section through
      // the first side's end)
      const double  aWs = wo;
      GeomAPI_IntCS anInt(new Geom_TrimmedCurve(aLo, Fio.FirstParameter(), Fio.LastParameter()),
                          aGPv);
      if (!anInt.IsDone() || anInt.NbPoints() == 0)
      {
        return false;
      }
      double aBest = RealLast();
      for (int i = 1; i <= anInt.NbPoints(); i++)
      {
        double aU, aV, aW;
        anInt.Parameters(i, aU, aV, aW);
        if (std::abs(aW - aWs) < aBest)
        {
          aBest = std::abs(aW - aWs);
          wo    = aW;
        }
      }
      Po = aLo->Value(wo);
    }
  }
  if (aPln.Distance(Po) > 10 * Tol3d)
  {
    return false;
  }
  // no faces at Vtx but these, the wall next to the other side's face
  const TopoDS_Face& Fbeside = onSplit ? Fup : Fo;
  bool               isNext  = false;
  for (NCollection_List<TopoDS_Shape>::Iterator It(VEMap(Vtx)); It.More() && !isNext; It.Next())
  {
    const TopoDS_Edge& E = TopoDS::Edge(It.Value());
    isNext               = containE(Fwall, E) && containE(Fbeside, E);
  }
  if (!isNext || ChFi3d_nbface(VFMap(Vtx)) != (onSplit ? 5 : 4)
      || !SectionCrossesEndFace(Vtx, Epiv, Fp, CVp.Point(), Po, VEMap, EFMap))
  {
    return false;
  }

  // the section, from side 1 to side 2, as the corners put their curves
  double                          aF, aL;
  const occ::handle<Geom2d_Curve> aPCpiv = PCurveInFace(Epiv, Fv, aF, aL);
  if (aPCpiv.IsNull())
  {
    return false;
  }
  const gp_Pnt2d aUVsp = Fip.PCurveOnSurf()->Value(Fip.Parameter(isFirst));
  const gp_Pnt2d aUVso = Fio.PCurveOnSurf()->Value(wo);
  const gp_Pnt2d aUVvp = aPCpiv->Value(CVp.ParameterOnArc());
  double         aUo, aVo;
  ElSLib::Parameters(aPln, Po, aUo, aVo);
  const gp_Pnt2d             aUVvo(aUo, aVo);
  const gp_Pnt2d&            aUVs1 = Ip == 1 ? aUVsp : aUVso;
  const gp_Pnt2d&            aUVs2 = Ip == 1 ? aUVso : aUVsp;
  const gp_Pnt2d&            aUVv1 = Ip == 1 ? aUVvp : aUVvo;
  const gp_Pnt2d&            aUVv2 = Ip == 1 ? aUVvo : aUVvp;
  NCollection_Array1<double> Pardeb(1, 4), Parfin(1, 4);
  Pardeb(1) = aUVs1.X();
  Pardeb(2) = aUVs1.Y();
  Pardeb(3) = aUVv1.X();
  Pardeb(4) = aUVv1.Y();
  Parfin(1) = aUVs2.X();
  Parfin(2) = aUVs2.Y();
  Parfin(3) = aUVv2.X();
  Parfin(4) = aUVv2.Y();
  occ::handle<BRepAdaptor_Surface> aHBv = new BRepAdaptor_Surface(Fv);
  double                           uu1, uu2, vv1, vv2;
  ChFi3d_Boite(aUVv1, aUVv2, uu1, uu2, vv1, vv2);
  ChFi3d_BoundFac(*aHBv, uu1, uu2, vv1, vv2);
  occ::handle<Geom_Curve>   aCc;
  occ::handle<Geom2d_Curve> aPs, aPc;
  double                    aTolR = Tol3d;
  if (!ChFi3d_ComputeCurves(HGs, aHBv, Pardeb, Parfin, aCc, aPs, aPc, Tol3d, Tol2d, aTolR, false)
      || aCc.IsNull() || aPs.IsNull() || aPc.IsNull())
  {
    return false;
  }
  // (two planes are given no tolerance but the one asked for)
  aTolR =
    std::min(aTolR,
             std::max(ChFi3d_EvalTolReached(HGs, aPs, aHBv, aPc, aCc), Precision::Confusion()));
  const double p1 = aCc->FirstParameter(), p2 = aCc->LastParameter();
  if (aCc->Value(p1).Distance(Ip == 1 ? CVp.Point() : Po) > 10 * Tol3d
      || aCc->Value(p2).Distance(Ip == 1 ? Po : CVp.Point()) > 10 * Tol3d)
  {
    return false;
  }

  // X: the section across the wall's edge, once, between the two sides
  double                          aFA, aLA;
  const occ::handle<Geom2d_Curve> aPCprolV = PCurveInFace(Arcprol, Fv, aFA, aLA);
  const occ::handle<Geom2d_Curve> aPCprolW = PCurveInFace(Arcprol, Fwall, aFA, aLA);
  const occ::handle<Geom_Curve>   aCprol   = BRep_Tool::Curve(Arcprol, aFA, aLA);
  if (aPCprolV.IsNull() || aPCprolW.IsNull() || aCprol.IsNull()
      || p2 - p1 <= Precision::PConfusion() || aLA - aFA <= Precision::PConfusion())
  {
    return false;
  }
  Geom2dAPI_InterCurveCurve anInt2d(new Geom2d_TrimmedCurve(aPc, p1, p2),
                                    new Geom2d_TrimmedCurve(aPCprolV, aFA, aLA),
                                    Precision::Confusion());
  if (anInt2d.NbPoints() != 1 || anInt2d.NbSegments() != 0)
  {
    return false;
  }
  const double tX   = anInt2d.Intersector().Point(1).ParamOnFirst();
  const double parX = anInt2d.Intersector().Point(1).ParamOnSecond();
  const double parV = BRep_Tool::Parameter(Vtx, Arcprol);
  const gp_Pnt PX   = aCc->Value(tX);
  const double aTolX =
    std::max({aTolR, PX.Distance(aCprol->Value(parX)), Precision::Confusion()});
  const double aLim =
    4. * std::max({aTolX, CVp.Tolerance(), BRep_Tool::Tolerance(Vtx), Precision::Confusion()});
  const double parFar = std::abs(parV - aFA) < std::abs(parV - aLA) ? aLA : aFA;
  if (tX <= p1 || tX >= p2 || parX <= aFA || parX >= aLA || PX.Distance(CVp.Point()) <= aLim
      || PX.Distance(Po) <= aLim || PX.Distance(PV) <= aLim || PV.Distance(Po) <= aLim
      || PX.Distance(aCprol->Value(parFar)) <= aLim)
  {
    return false;
  }

  // the ceiling's side on the other face: the split from Vtx to the line's
  // end, the face above keeping it; or, the face in one, the line from Vtx
  // to the point on it
  occ::handle<Geom_Curve>   aC4;
  occ::handle<Geom2d_Curve> aPc4f;
  occ::handle<Geom_Surface> aS4;
  TopAbs_Orientation        anOr4f;
  double                    a4 = 0., b4 = 0.;
  bool                      isVtxFirst4 = true;
  if (onSplit)
  {
    double                        aFS, aLS;
    const occ::handle<Geom_Curve> aCs = BRep_Tool::Curve(Esplit, aFS, aLS);
    aPc4f                             = PCurveInFace(Esplit, Fup, aFS, aLS);
    if (aCs.IsNull() || aPc4f.IsNull())
    {
      return false;
    }
    const double parVs = BRep_Tool::Parameter(Vtx, Esplit), parO = CVo.ParameterOnArc();
    a4          = std::min(parVs, parO);
    b4          = std::max(parVs, parO);
    isVtxFirst4 = parVs < parO;
    if (b4 - a4 <= Precision::PConfusion())
    {
      return false;
    }
    aC4    = new Geom_TrimmedCurve(aCs, a4, b4);
    aS4    = BRep_Tool::Surface(Fup);
    anOr4f = anOrIn(Esplit, Fup);
  }
  else
  {
    const occ::handle<Geom_Plane> aGPo = aPlaneOf(Fo);
    if (aGPo.IsNull())
    {
      return false;
    }
    b4           = PV.Distance(Po);
    aC4          = new Geom_TrimmedCurve(new Geom_Line(PV, gp_Dir(gp_Vec(PV, Po))), a4, b4);
    aS4          = aGPo;
    double aTolP = Tol3d;
    aPc4f        = GeomProjLib::Curve2d(aC4, a4, b4, aS4, aTolP);
    // as the line of the fillet runs on that face, to its end and on to Vtx
    anOr4f = isFirst ? Fio.Transition() : TopAbs::Reverse(Fio.Transition());
  }
  double                          aTolP4 = Tol3d;
  const occ::handle<Geom2d_Curve> aPc4c  = GeomProjLib::Curve2d(aC4, a4, b4, aGPv, aTolP4);
  if (aPc4f.IsNull() || aPc4c.IsNull())
  {
    return false;
  }
  const occ::handle<GeomAdaptor_Surface> aHSv = new GeomAdaptor_Surface(aGPv);
  const double                           aTol4 =
    std::max(ChFi3d_EvalTolReached(new GeomAdaptor_Surface(aS4), aPc4f, aHSv, aPc4c, aC4),
             Precision::Confusion());
  if (aTol4 > 10 * Tol3d)
  {
    return false;
  }

  // -- the DS, from here on --
  // The ceiling lies in the plane of the face at end, turned the other
  // way and on the other side of the section and of the wall's edge: its
  // curves the other way round on the same surface.
  const TopAbs_Orientation anOrC = TopAbs::Reverse(Fv.Orientation());
  const TopAbs_Orientation anOrienV =
    Fv.Orientation() == Fd->Orientation() ? TopAbs::Reverse(OrCurve) : OrCurve;
  const TopAbs_Orientation anOrienC = TopAbs::Reverse(anOrienV);
  const int                Isurf    = Fd->Surf();
  const int                IsurfC   = DStr.AddSurface(
    TopOpeBRepDS_Surface(occ::down_cast<Geom_Surface>(aGPv->Copy()), BRep_Tool::Tolerance(Fv)));
  if (!EVIMap.IsBound(Vtx))
  {
    EVIMap.Bind(Vtx, NCollection_List<int>());
  }
  EVIMap.ChangeFind(Vtx).Append(IsurfC);
  const int SolInd = Stripe->SolidIndex();
  DStr.ChangeShapeInterferences(SolInd).Append(
    new TopOpeBRepDS_SolidSurfaceInterference(TopOpeBRepDS_Transition(anOrC),
                                              TopOpeBRepDS_SOLID,
                                              SolInd,
                                              TopOpeBRepDS_SURFACE,
                                              IsurfC));

  // the points
  if (!onSplit)
  {
    // the line is cut where it crosses the plane, as an edge splitting
    // its face there would cut it
    Fd->ChangeInterference(Io).SetParameter(wo, isFirst);
    CVo.Reset();
    CVo.SetPoint(Po);
    CVo.SetTolerance(std::max(aPln.Distance(Po), Precision::Confusion()));
    Stripe->SetIndexPoint(ChFi3d_IndexPointInDS(CVo, DStr), isFirst, Io);
  }
  const int IndP = Stripe->IndexPoint(isFirst, Ip);
  const int IndO = Stripe->IndexPoint(isFirst, Io);
  const int IndX = DStr.AddPoint(TopOpeBRepDS_Point(PX, aTolX));
  const int IVtx = DStr.AddShape(Vtx);
  for (int k = 0; k < 2; k++)
  {
    TopOpeBRepDS_Point& aTP = DStr.ChangePoint(k == 0 ? IndP : IndO);
    aTP.Tolerance(std::max(aTP.Tolerance(), aTolR));
  }

  // the section in two: on the face at end up to X, on the ceiling beyond
  const int IFv         = DStr.AddShape(Fv);
  auto      aPutSection = [&](const double theA,
                         const double theB,
                         const int    theInd1,
                         const int    theInd2,
                         const bool   theOnCeiling) {
    const int IC = DStr.AddCurve(TopOpeBRepDS_Curve(new Geom_TrimmedCurve(aCc, theA, theB), aTolR));
    DStr.ChangeCurveInterferences(IC).Append(
      ChFi3d_FilPointInDS(TopAbs_FORWARD, IC, theInd1, theA));
    DStr.ChangeCurveInterferences(IC).Append(
      ChFi3d_FilPointInDS(TopAbs_REVERSED, IC, theInd2, theB));
    DStr.ChangeSurfaceInterferences(Isurf).Append(
      ChFi3d_FilCurveInDS(IC, Isurf, new Geom2d_TrimmedCurve(aPs, theA, theB), OrCurve));
    const occ::handle<Geom2d_Curve> aPcAB = new Geom2d_TrimmedCurve(aPc, theA, theB);
    if (theOnCeiling)
    {
      DStr.ChangeSurfaceInterferences(IsurfC).Append(
        ChFi3d_FilCurveInDS(IC, IsurfC, aPcAB, anOrienC));
    }
    else
    {
      DStr.ChangeShapeInterferences(IFv).Append(ChFi3d_FilCurveInDS(IC, IFv, aPcAB, anOrienV));
    }
  };
  if (Ip == 1)
  {
    aPutSection(p1, tX, IndP, IndX, false);
    aPutSection(tX, p2, IndX, IndO, true);
  }
  else
  {
    aPutSection(p1, tX, IndO, IndX, true);
    aPutSection(tX, p2, IndX, IndP, false);
  }

  // an edge from Vtx between a face of the body and the ceiling
  auto aPutSide = [&](const occ::handle<Geom_Curve>&   theC,
                      const double                     theA,
                      const double                     theB,
                      const double                     theTol,
                      const bool                       theVtxFirst,
                      const int                        theInd,
                      const TopoDS_Face&               theF,
                      const occ::handle<Geom2d_Curve>& thePcF,
                      const TopAbs_Orientation         theOrF,
                      const occ::handle<Geom2d_Curve>& thePcC,
                      const TopAbs_Orientation         theOrC) {
    const int IC = DStr.AddCurve(TopOpeBRepDS_Curve(theC, theTol));
    DStr.ChangeCurveInterferences(IC).Append(
      ChFi3d_FilPointInDS(TopAbs_FORWARD, IC, theVtxFirst ? IVtx : theInd, theA, theVtxFirst));
    DStr.ChangeCurveInterferences(IC).Append(
      ChFi3d_FilPointInDS(TopAbs_REVERSED, IC, theVtxFirst ? theInd : IVtx, theB, !theVtxFirst));
    const int IF = DStr.AddShape(theF);
    DStr.ChangeShapeInterferences(IF).Append(ChFi3d_FilCurveInDS(IC, IF, thePcF, theOrF));
    DStr.ChangeSurfaceInterferences(IsurfC).Append(ChFi3d_FilCurveInDS(IC, IsurfC, thePcC, theOrC));
  };
  // the wall's edge from Vtx to X: the face at end lay on its one side,
  // the ceiling lies on the other
  aPutSide(new Geom_TrimmedCurve(aCprol, std::min(parV, parX), std::max(parV, parX)),
           std::min(parV, parX),
           std::max(parV, parX),
           std::max(BRep_Tool::Tolerance(Arcprol), Precision::Confusion()),
           parV < parX,
           IndX,
           Fwall,
           aPCprolW,
           anOrIn(Arcprol, Fwall),
           aPCprolV,
           TopAbs::Reverse(anOrIn(Arcprol, Fv)));
  // the other face's: the two faces as two faces of a shell meeting there
  const TopoDS_Face& F4 = onSplit ? Fup : Fo;
  aPutSide(aC4,
           a4,
           b4,
           aTol4,
           isVtxFirst4,
           IndO,
           F4,
           aPc4f,
           anOr4f,
           aPc4c,
           anOrC == F4.Orientation() ? TopAbs::Reverse(anOr4f) : anOr4f);

  // the edges cut
  const int IEpiv = DStr.AddShape(Epiv);
  DStr.ChangeShapeInterferences(IEpiv).Append(
    ChFi3d_FilPointInDS(CVp.TransitionOnArc(), IEpiv, IndP, CVp.ParameterOnArc()));
  const int IArcprol = DStr.AddShape(Arcprol);
  DStr.ChangeShapeInterferences(IArcprol).Append(
    ChFi3d_FilPointInDS(TopExp::FirstVertex(Arcprol).IsSame(Vtx) ? TopAbs_FORWARD : TopAbs_REVERSED,
                        IArcprol,
                        IndX,
                        parX));
  if (onSplit)
  {
    const int IEsplit = DStr.AddShape(Esplit);
    DStr.ChangeShapeInterferences(IEsplit).Append(
      ChFi3d_FilPointInDS(CVo.TransitionOnArc(), IEsplit, IndO, CVo.ParameterOnArc()));
  }
  // Vtx stays, with the wall's edges and the ceiling's: the edge filleted
  // ends there with nothing kept of it, or the faces beside it that keep
  // Vtx would keep the edge
  TopAbs_Orientation anOrV = TopAbs_FORWARD;
  for (TopExp_Explorer ex(EdgeSpine.Oriented(TopAbs_FORWARD), TopAbs_VERTEX); ex.More(); ex.Next())
  {
    if (Vtx.IsSame(ex.Current()))
    {
      anOrV = ex.Current().Orientation();
      break;
    }
  }
  const int ISpine = DStr.AddShape(EdgeSpine);
  DStr.ChangeShapeInterferences(ISpine).Append(
    ChFi3d_FilVertexInDS(TopAbs::Reverse(anOrV),
                         ISpine,
                         IVtx,
                         BRep_Tool::Parameter(Vtx, EdgeSpine)));
  Stripe->InDS(isFirst);
  return true;
}

//=======================================================================
// function : IsRimEnd
// purpose  : Whether the end of <Fd> at <Vtx> is a rim end (MakeRimFlatEnd):
//           one side's common point on an edge of <Vtx> (Epiv), short of
//           <Vtx>, the other's on no edge of <Vtx>, and the two faces those
//           of <EdgeSpine>. Sets <Ip> to that side.
//=======================================================================

static bool IsRimEnd(const TopOpeBRepDS_DataStructure&   DStr,
                     const occ::handle<ChFiDS_SurfData>& Fd,
                     const bool                          isFirst,
                     const TopoDS_Vertex&                Vtx,
                     const TopoDS_Edge&                  EdgeSpine,
                     const double                        Tol3d,
                     int&                                Ip)
{
  Ip = 0;
  if (Fd.IsNull() || Fd->Index(1) <= 0 || Fd->Index(2) <= 0 || Fd->Surf() <= 0 || EdgeSpine.IsNull()
      || !hasVertex(EdgeSpine, Vtx) || DStr.Shape(Fd->Index(1)).ShapeType() != TopAbs_FACE
      || DStr.Shape(Fd->Index(2)).ShapeType() != TopAbs_FACE)
  {
    return false;
  }
  const gp_Pnt PV = BRep_Tool::Pnt(Vtx);
  for (int i = 1; i <= 2; i++)
  {
    const ChFiDS_CommonPoint& aCV = Fd->Vertex(isFirst, i);
    if (aCV.IsOnArc() && hasVertex(aCV.Arc(), Vtx))
    {
      if (aCV.IsVertex() || aCV.Point().Distance(PV) <= 10 * Tol3d || aCV.Arc().IsSame(EdgeSpine))
      {
        return false;
      }
      Ip = Ip == 0 ? i : -1;
    }
  }
  if (Ip <= 0)
  {
    Ip = 0;
    return false;
  }
  const TopoDS_Face Fp = TopoDS::Face(DStr.Shape(Fd->Index(Ip)));
  const TopoDS_Face Fo = TopoDS::Face(DStr.Shape(Fd->Index(3 - Ip)));
  return containE(Fp, EdgeSpine) && containE(Fo, EdgeSpine)
         && containE(Fp, Fd->Vertex(isFirst, Ip).Arc());
}

//=======================================================================
// function : MakeRimFlatEnd
// purpose  : The end of a fillet at <Vtx> flat in its section where one
//           side's line runs off its face (Fp) over an edge of <Vtx> (Epiv),
//           short of <Vtx> (IsRimEnd): the face narrows to nothing at <Vtx>
//           between the edge filleted and Epiv, tangent there -- a post's
//           foot on a plate whose rim runs flush into the post's rounded
//           corner (#876). Nothing past Epiv could take the fillet on: the
//           face at end carried on misses the other side's line, or meets it
//           only tangent to it (a post standing flush). The fillet ends in
//           the plane square to the spine through the point on Epiv; the
//           face it leaves there, flat, lies between the section, the line
//           to the spine (Pf) on Fp and the line from Pf to the other side's
//           point on its face (Fo). The edge filleted is kept from Pf to
//           <Vtx>, Epiv whole: what lies between them, on Fp, stays. All of
//           it is put in the DS here. False, and nothing done, where the
//           face cannot be made.
//=======================================================================

static bool MakeRimFlatEnd(
  TopOpeBRepDS_DataStructure&                                                        DStr,
  const occ::handle<ChFiDS_Stripe>&                                                  Stripe,
  const occ::handle<ChFiDS_SurfData>&                                                Fd,
  const bool                                                                         isFirst,
  const TopoDS_Vertex&                                                               Vtx,
  const TopoDS_Edge&                                                                 EdgeSpine,
  const TopoDS_Shape&                                                                theSolid,
  const ChFiDS_Map&                                                                  EFMap,
  const double                                                                       Tol3d,
  const double                                                                       Tol2d,
  NCollection_DataMap<TopoDS_Shape, NCollection_List<int>, TopTools_ShapeMapHasher>& EVIMap)
{
  int Ip = 0;
  if (!IsRimEnd(DStr, Fd, isFirst, Vtx, EdgeSpine, Tol3d, Ip))
  {
    return false;
  }
  const int                      Io   = 3 - Ip;
  const ChFiDS_CommonPoint&      CVp  = Fd->Vertex(isFirst, Ip);
  ChFiDS_CommonPoint&            CVo  = Fd->ChangeVertex(isFirst, Io);
  const TopoDS_Face              Fp   = TopoDS::Face(DStr.Shape(Fd->Index(Ip)));
  const TopoDS_Face              Fo   = TopoDS::Face(DStr.Shape(Fd->Index(Io)));
  const TopoDS_Edge              Epiv = CVp.Arc();
  const ChFiDS_FaceInterference& Fip  = Fd->Interference(Ip);
  const ChFiDS_FaceInterference& Fio  = Fd->Interference(Io);
  if (Fip.LineIndex() == 0 || Fio.LineIndex() == 0 || Fip.PCurveOnSurf().IsNull()
      || Fio.PCurveOnSurf().IsNull() || Fio.PCurveOnFace().IsNull() || !EFMap.Contains(Epiv)
      || EFMap(Epiv).Extent() != 2)
  {
    return false;
  }
  TopoDS_Face Fv;
  for (NCollection_List<TopoDS_Shape>::Iterator It(EFMap(Epiv)); It.More(); It.Next())
  {
    if (!Fp.IsSame(It.Value()))
    {
      Fv = TopoDS::Face(It.Value());
    }
  }
  if (Fv.IsNull() || Fv.IsSame(Fo))
  {
    return false;
  }
  // Fv of a width beside Epiv, not a needle whose other side Epiv lies on
  // (#474's Fillet003, a plane of no area at the top of its edge 51: Epiv's
  // piece kept would lie on the needle's other edge, the plane's beyond)
  {
    double                        aFE, aLE;
    const occ::handle<Geom_Curve> aCE = BRep_Tool::Curve(Epiv, aFE, aLE);
    if (aCE.IsNull())
    {
      return false;
    }
    const gp_Pnt aMid = aCE->Value(0.5 * (aFE + aLE));
    for (TopExp_Explorer ex(Fv, TopAbs_EDGE); ex.More(); ex.Next())
    {
      const TopoDS_Edge& anE = TopoDS::Edge(ex.Current());
      double             aF2, aL2;
      if (anE.IsSame(Epiv) || BRep_Tool::Degenerated(anE))
      {
        continue;
      }
      const occ::handle<Geom_Curve> aC2 = BRep_Tool::Curve(anE, aF2, aL2);
      if (aC2.IsNull())
      {
        continue;
      }
      GeomAPI_ProjectPointOnCurve aPr(aMid, aC2, aF2, aL2);
      if (aPr.NbPoints() > 0 && aPr.LowerDistance() <= Tol3d)
      {
        return false;
      }
    }
  }
  auto anOrIn = [](const TopoDS_Edge& theE, const TopoDS_Face& theF) {
    for (TopExp_Explorer ex(theF.Oriented(TopAbs_FORWARD), TopAbs_EDGE); ex.More(); ex.Next())
    {
      if (theE.IsSame(ex.Current()))
      {
        return ex.Current().Orientation();
      }
    }
    return theE.Orientation();
  };
  const gp_Pnt PV = BRep_Tool::Pnt(Vtx);
  const gp_Pnt PP = CVp.Point();

  // Pf, the spine's point square to PP, inside the edge filleted
  double                        aFS, aLS;
  const occ::handle<Geom_Curve> aCS = BRep_Tool::Curve(EdgeSpine, aFS, aLS);
  if (aCS.IsNull())
  {
    return false;
  }
  GeomAPI_ProjectPointOnCurve aProjS(PP, aCS, aFS, aLS);
  if (aProjS.NbPoints() == 0)
  {
    return false;
  }
  const double parF = aProjS.LowerDistanceParameter();
  const double parV = BRep_Tool::Parameter(Vtx, EdgeSpine);
  const gp_Pnt Pf   = aCS->Value(parF);
  if (std::abs(parF - aFS) <= Precision::PConfusion()
      || std::abs(parF - aLS) <= Precision::PConfusion() || Pf.Distance(PV) <= 10 * Tol3d
      || Pf.Distance(PP) <= 10 * Tol3d)
  {
    return false;
  }
  // the spine's way to Vtx there
  gp_Pnt aPbid;
  gp_Vec aT;
  aCS->D1(parF, aPbid, aT);
  if (aT.Magnitude() <= gp::Resolution())
  {
    return false;
  }
  if (parV < parF)
  {
    aT.Reverse();
  }
  aT.Normalize();
  // the plane through Pf and PP, as square to the spine as that allows
  gp_Vec aChord(Pf, PP);
  aChord.Normalize();
  gp_Vec aN = aT - aChord * aT.Dot(aChord);
  if (aN.Magnitude() <= 0.9)
  {
    return false;
  }
  aN.Normalize();

  // the line on Fp at PP, and Fo's line across the plane nearest PP
  const occ::handle<Geom_Curve>& aLp = DStr.Curve(Fip.LineIndex()).Curve();
  const occ::handle<Geom_Curve>& aLo = DStr.Curve(Fio.LineIndex()).Curve();
  if (aLp.IsNull() || aLo.IsNull())
  {
    return false;
  }
  double wp = Fip.Parameter(isFirst);
  if (aLp->Value(wp).Distance(PP) > Tol3d)
  {
    GeomAPI_ProjectPointOnCurve aProjP(PP, aLp, Fip.FirstParameter(), Fip.LastParameter());
    if (aProjP.NbPoints() == 0 || aProjP.LowerDistance() > 10 * Tol3d)
    {
      return false;
    }
    wp = aProjP.LowerDistanceParameter();
  }
  // Fo's line at its end where that lies in the plane (the walk stopped at
  // the section), else where it crosses it nearest PP
  const occ::handle<Geom_Plane> aPl0  = new Geom_Plane(Pf, gp_Dir(aN));
  double                        wo    = Fio.Parameter(isFirst);
  double                        aBest = PP.Distance(aLo->Value(wo));
  if (aPl0->Pln().Distance(aLo->Value(wo)) > Tol3d)
  {
    GeomAPI_IntCS anInt(new Geom_TrimmedCurve(aLo, Fio.FirstParameter(), Fio.LastParameter()),
                        aPl0);
    if (!anInt.IsDone() || anInt.NbPoints() == 0)
    {
      return false;
    }
    aBest = RealLast();
    for (int i = 1; i <= anInt.NbPoints(); i++)
    {
      double aU, aV, aW;
      anInt.Parameters(i, aU, aV, aW);
      if (PP.Distance(anInt.Point(i)) < aBest)
      {
        aBest = PP.Distance(anInt.Point(i));
        wo    = aW;
      }
    }
  }
  const gp_Pnt Po = aLo->Value(wo);
  // a section, no further across than four times the way to the spine
  if (aBest > 4. * Pf.Distance(PP) || Po.Distance(Pf) <= 10 * Tol3d)
  {
    return false;
  }
  BRepTopAdaptor_FClass2d aClo(Fo, Precision::PConfusion());
  if (aClo.Perform(Fio.PCurveOnFace()->Value(wo)) != TopAbs_IN)
  {
    return false;
  }

  // The flat face's outer side: toward Vtx where the fillet adds to the
  // solid (the swept region, between the section and Pf, lies outside it),
  // away from Vtx where it takes away
  const gp_Pnt2d                   aUVsp = Fip.PCurveOnSurf()->Value(wp);
  const gp_Pnt2d                   aUVso = Fio.PCurveOnSurf()->Value(wo);
  const occ::handle<Geom_Surface>& aSfil = DStr.Surface(Fd->Surf()).Surface();
  {
    const gp_Pnt aMid = aSfil->Value(0.5 * (aUVsp.X() + aUVso.X()), 0.5 * (aUVsp.Y() + aUVso.Y()));
    const gp_Pnt aIn((Pf.XYZ() + aMid.XYZ()) * 0.5);
    BRepClass3d_SolidClassifier aCl(theSolid, aIn, Tol3d);
    if (aCl.State() == TopAbs_ON || aCl.State() == TopAbs_UNKNOWN)
    {
      return false;
    }
    if (aCl.State() == TopAbs_IN)
    {
      aN.Reverse();
    }
  }
  const occ::handle<Geom_Plane> aGPl = new Geom_Plane(Pf, gp_Dir(aN));
  double                        aUF, aVF, aUP, aVP, aUO, aVO;
  ElSLib::Parameters(aGPl->Pln(), Pf, aUF, aVF);
  ElSLib::Parameters(aGPl->Pln(), PP, aUP, aVP);
  ElSLib::Parameters(aGPl->Pln(), Po, aUO, aVO);
  const double                           aMarg = Pf.Distance(PP) + Pf.Distance(Po);
  const occ::handle<GeomAdaptor_Surface> aHPl =
    new GeomAdaptor_Surface(aGPl,
                            std::min({aUF, aUP, aUO}) - aMarg,
                            std::max({aUF, aUP, aUO}) + aMarg,
                            std::min({aVF, aVP, aVO}) - aMarg,
                            std::max({aVF, aVP, aVO}) + aMarg);

  // the section, from side 1 to side 2, on the fillet's surface about its
  // two ends
  const double aDu = std::abs(aUVsp.X() - aUVso.X()), aDv = std::abs(aUVsp.Y() - aUVso.Y());
  const double aM2 = std::max({aDu, aDv, 1.e-3});
  const occ::handle<GeomAdaptor_Surface> aHGs =
    new GeomAdaptor_Surface(aSfil,
                            std::min(aUVsp.X(), aUVso.X()) - aM2,
                            std::max(aUVsp.X(), aUVso.X()) + aM2,
                            std::min(aUVsp.Y(), aUVso.Y()) - 0.1 * aM2,
                            std::max(aUVsp.Y(), aUVso.Y()) + 0.1 * aM2);
  const gp_Pnt2d             aUVlp(aUP, aVP), aUVlo(aUO, aVO);
  const gp_Pnt2d&            aUVs1 = Ip == 1 ? aUVsp : aUVso;
  const gp_Pnt2d&            aUVs2 = Ip == 1 ? aUVso : aUVsp;
  const gp_Pnt2d&            aUVl1 = Ip == 1 ? aUVlp : aUVlo;
  const gp_Pnt2d&            aUVl2 = Ip == 1 ? aUVlo : aUVlp;
  NCollection_Array1<double> Pardeb(1, 4), Parfin(1, 4);
  Pardeb(1) = aUVs1.X();
  Pardeb(2) = aUVs1.Y();
  Pardeb(3) = aUVl1.X();
  Pardeb(4) = aUVl1.Y();
  Parfin(1) = aUVs2.X();
  Parfin(2) = aUVs2.Y();
  Parfin(3) = aUVl2.X();
  Parfin(4) = aUVl2.Y();
  occ::handle<Geom_Curve>   aCc;
  occ::handle<Geom2d_Curve> aPs, aPc;
  double                    aTolR = Tol3d;
  if (!ChFi3d_ComputeCurves(aHGs, aHPl, Pardeb, Parfin, aCc, aPs, aPc, Tol3d, Tol2d, aTolR, false)
      || aCc.IsNull() || aPs.IsNull() || aPc.IsNull())
  {
    return false;
  }
  aTolR =
    std::min(aTolR,
             std::max(ChFi3d_EvalTolReached(aHGs, aPs, aHPl, aPc, aCc), Precision::Confusion()));
  const double p1 = aCc->FirstParameter(), p2 = aCc->LastParameter();
  if (aTolR > 10 * Tol3d || aCc->Value(p1).Distance(Ip == 1 ? PP : Po) > 10 * Tol3d
      || aCc->Value(p2).Distance(Ip == 1 ? Po : PP) > 10 * Tol3d)
  {
    return false;
  }

  // the two lines from Pf, each in its face and the plane
  struct Side
  {
    occ::handle<Geom_Curve>   C;
    occ::handle<Geom2d_Curve> PcF, PcL;
    double                    L = 0., Tol = 0.;
  };

  auto aLineFrom = [&](const gp_Pnt& theTo, const TopoDS_Face& theF, Side& theS) {
    theS.L = Pf.Distance(theTo);
    theS.C = new Geom_TrimmedCurve(new Geom_Line(Pf, gp_Dir(gp_Vec(Pf, theTo))), 0., theS.L);
    const occ::handle<Geom_Surface> aS     = BRep_Tool::Surface(theF);
    double                          aTolPj = Tol3d;
    theS.PcF                               = GeomProjLib::Curve2d(theS.C, 0., theS.L, aS, aTolPj);
    aTolPj                                 = Tol3d;
    theS.PcL                               = GeomProjLib::Curve2d(theS.C, 0., theS.L, aGPl, aTolPj);
    if (theS.PcF.IsNull() || theS.PcL.IsNull())
    {
      return false;
    }
    BRepTopAdaptor_FClass2d aCl(theF, Precision::PConfusion());
    if (aCl.Perform(theS.PcF->Value(0.5 * theS.L)) != TopAbs_IN)
    {
      return false;
    }
    theS.Tol =
      std::max(ChFi3d_EvalTolReached(new GeomAdaptor_Surface(aS), theS.PcF, aHPl, theS.PcL, theS.C),
               Precision::Confusion());
    return theS.Tol <= 10 * Tol3d;
  };
  Side aSp, aSo;
  if (!aLineFrom(PP, Fp, aSp) || !aLineFrom(Po, Fo, aSo))
  {
    return false;
  }
  // Epiv's piece from the point to Vtx, which Fp and Fv keep
  double                          aFV, aLV, aFP, aLP;
  double                          aFE, aLE;
  const occ::handle<Geom_Curve>   aCpiv = BRep_Tool::Curve(Epiv, aFE, aLE);
  const occ::handle<Geom2d_Curve> aPCpV = PCurveInFace(Epiv, Fv, aFV, aLV);
  const occ::handle<Geom2d_Curve> aPCpP = PCurveInFace(Epiv, Fp, aFP, aLP);
  const double                    parEV = BRep_Tool::Parameter(Vtx, Epiv);
  const double                    parEP = CVp.ParameterOnArc();
  if (aCpiv.IsNull() || aPCpV.IsNull() || aPCpP.IsNull()
      || std::abs(parEP - parEV) <= Precision::PConfusion())
  {
    return false;
  }
  const double                    aE0 = std::min(parEV, parEP), aE1 = std::max(parEV, parEP);
  const occ::handle<Geom_Curve>   aCr   = new Geom_TrimmedCurve(aCpiv, aE0, aE1);
  const occ::handle<Geom2d_Curve> aPcrV = new Geom2d_TrimmedCurve(aPCpV, aE0, aE1);
  const occ::handle<Geom2d_Curve> aPcrP = new Geom2d_TrimmedCurve(aPCpP, aE0, aE1);
  const double aTolCr = std::max(BRep_Tool::Tolerance(Epiv), Precision::Confusion());

  // The sides. The flat face lies, from the section, toward Pf; from the
  // line on Fp, toward Po; from the line on Fo, toward PP. What Fp and Fo
  // keep past their lines lies toward Vtx.
  const TopAbs_Orientation anOrSL =
    IsKeptLeft(aGPl, aPc, aCc, gp_Vec(aCc->Value(0.5 * (p1 + p2)), Pf)) ? TopAbs_FORWARD
                                                                        : TopAbs_REVERSED;
  const TopAbs_Orientation anOrPL =
    IsKeptLeft(aGPl, aSp.PcL, aSp.C, gp_Vec(aSp.C->Value(0.5 * aSp.L), Po)) ? TopAbs_FORWARD
                                                                            : TopAbs_REVERSED;
  const TopAbs_Orientation anOrOL =
    IsKeptLeft(aGPl, aSo.PcL, aSo.C, gp_Vec(aSo.C->Value(0.5 * aSo.L), PP)) ? TopAbs_FORWARD
                                                                            : TopAbs_REVERSED;
  const TopAbs_Orientation anOrPF =
    IsKeptLeft(BRep_Tool::Surface(Fp), aSp.PcF, aSp.C, aT) ? TopAbs_FORWARD : TopAbs_REVERSED;
  const TopAbs_Orientation anOrOF =
    IsKeptLeft(BRep_Tool::Surface(Fo), aSo.PcF, aSo.C, aT) ? TopAbs_FORWARD : TopAbs_REVERSED;
  TopAbs_Orientation aTrafil1 = TopAbs_FORWARD;
  if (Fd->IndexOfS1() > 0)
  {
    aTrafil1 = DStr.Shape(Fd->IndexOfS1()).Orientation();
  }
  aTrafil1 = TopAbs::Compose(aTrafil1, Fd->Orientation());
  aTrafil1 = TopAbs::Compose(TopAbs::Reverse(Fd->InterferenceOnS1().Transition()), aTrafil1);
  const TopAbs_Orientation anOrS = isFirst ? TopAbs::Reverse(aTrafil1) : aTrafil1;

  // -- the DS, from here on --
  const int Isurf  = Fd->Surf();
  const int IsurfC = DStr.AddSurface(TopOpeBRepDS_Surface(aGPl, Tol3d));
  if (!EVIMap.IsBound(Vtx))
  {
    EVIMap.Bind(Vtx, NCollection_List<int>());
  }
  EVIMap.ChangeFind(Vtx).Append(IsurfC);
  const int SolInd = Stripe->SolidIndex();
  DStr.ChangeShapeInterferences(SolInd).Append(
    new TopOpeBRepDS_SolidSurfaceInterference(TopOpeBRepDS_Transition(TopAbs_FORWARD),
                                              TopOpeBRepDS_SOLID,
                                              SolInd,
                                              TopOpeBRepDS_SURFACE,
                                              IsurfC));
  // the points: the lines' ends at the section, and Pf
  if (std::abs(wp - Fip.Parameter(isFirst)) > Precision::PConfusion())
  {
    Fd->ChangeInterference(Ip).SetParameter(wp, isFirst);
  }
  Fd->ChangeInterference(Io).SetParameter(wo, isFirst);
  if (Po.Distance(CVo.Point()) > Precision::Confusion())
  {
    CVo.Reset();
    CVo.SetPoint(Po);
    CVo.SetTolerance(std::max(aGPl->Pln().Distance(Po), Precision::Confusion()));
  }
  Stripe->SetIndexPoint(ChFi3d_IndexPointInDS(CVp, DStr), isFirst, Ip);
  Stripe->SetIndexPoint(ChFi3d_IndexPointInDS(CVo, DStr), isFirst, Io);
  const int IndP = Stripe->IndexPoint(isFirst, Ip);
  const int IndO = Stripe->IndexPoint(isFirst, Io);
  const int IndF = DStr.AddPoint(TopOpeBRepDS_Point(Pf, std::max(aSp.Tol, aSo.Tol)));
  {
    TopOpeBRepDS_Point& aTP = DStr.ChangePoint(IndP);
    aTP.Tolerance(std::max({aTP.Tolerance(), aTolR, aSp.Tol, aTolCr}));
    TopOpeBRepDS_Point& aTO = DStr.ChangePoint(IndO);
    aTO.Tolerance(std::max({aTO.Tolerance(), aTolR, aSo.Tol}));
  }
  // the section
  const int ICs = DStr.AddCurve(TopOpeBRepDS_Curve(aCc, aTolR));
  DStr.ChangeCurveInterferences(ICs).Append(
    ChFi3d_FilPointInDS(TopAbs_FORWARD, ICs, Ip == 1 ? IndP : IndO, p1));
  DStr.ChangeCurveInterferences(ICs).Append(
    ChFi3d_FilPointInDS(TopAbs_REVERSED, ICs, Ip == 1 ? IndO : IndP, p2));
  DStr.ChangeSurfaceInterferences(Isurf).Append(ChFi3d_FilCurveInDS(ICs, Isurf, aPs, anOrS));
  DStr.ChangeSurfaceInterferences(IsurfC).Append(ChFi3d_FilCurveInDS(ICs, IsurfC, aPc, anOrSL));
  // the lines from Pf on Fp and on Fo
  const int IFp = DStr.AddShape(Fp);
  const int IFo = DStr.AddShape(Fo);
  for (int k = 0; k < 2; k++)
  {
    const Side& aS   = k == 0 ? aSp : aSo;
    const int   ICl  = DStr.AddCurve(TopOpeBRepDS_Curve(aS.C, aS.Tol));
    const int   IndE = k == 0 ? IndP : IndO;
    DStr.ChangeCurveInterferences(ICl).Append(ChFi3d_FilPointInDS(TopAbs_FORWARD, ICl, IndF, 0.));
    DStr.ChangeCurveInterferences(ICl).Append(
      ChFi3d_FilPointInDS(TopAbs_REVERSED, ICl, IndE, aS.L));
    DStr.ChangeShapeInterferences(k == 0 ? IFp : IFo)
      .Append(ChFi3d_FilCurveInDS(ICl, k == 0 ? IFp : IFo, aS.PcF, k == 0 ? anOrPF : anOrOF));
    DStr.ChangeSurfaceInterferences(IsurfC).Append(
      ChFi3d_FilCurveInDS(ICl, IsurfC, aS.PcL, k == 0 ? anOrPL : anOrOL));
  }
  // Epiv's piece from the point to Vtx, in Fp and Fv as Epiv is, and Epiv
  // cut at the point
  const int  IVtx     = DStr.AddShape(Vtx);
  const int  ICr      = DStr.AddCurve(TopOpeBRepDS_Curve(aCr, aTolCr));
  const bool isVtxLow = parEV < parEP;
  DStr.ChangeCurveInterferences(ICr).Append(
    ChFi3d_FilPointInDS(TopAbs_FORWARD, ICr, isVtxLow ? IVtx : IndP, aE0, isVtxLow));
  DStr.ChangeCurveInterferences(ICr).Append(
    ChFi3d_FilPointInDS(TopAbs_REVERSED, ICr, isVtxLow ? IndP : IVtx, aE1, !isVtxLow));
  const int IFv = DStr.AddShape(Fv);
  DStr.ChangeShapeInterferences(IFv).Append(ChFi3d_FilCurveInDS(ICr, IFv, aPcrV, anOrIn(Epiv, Fv)));
  DStr.ChangeShapeInterferences(IFp).Append(ChFi3d_FilCurveInDS(ICr, IFp, aPcrP, anOrIn(Epiv, Fp)));
  const int IEpiv = DStr.AddShape(Epiv);
  DStr.ChangeShapeInterferences(IEpiv).Append(
    ChFi3d_FilPointInDS(CVp.TransitionOnArc(), IEpiv, IndP, parEP));
  // the edge filleted cut at Pf, kept from there to Vtx
  TopAbs_Orientation anOrV = TopAbs_FORWARD;
  for (TopExp_Explorer ex(EdgeSpine.Oriented(TopAbs_FORWARD), TopAbs_VERTEX); ex.More(); ex.Next())
  {
    if (Vtx.IsSame(ex.Current()))
    {
      anOrV = ex.Current().Orientation();
      break;
    }
  }
  const int ISpine = DStr.AddShape(EdgeSpine);
  DStr.ChangeShapeInterferences(ISpine).Append(
    ChFi3d_FilPointInDS(TopAbs::Reverse(anOrV), ISpine, IndF, parF));
  Stripe->InDS(isFirst);
  return true;
}

//=======================================================================
// function : RemoveSurfDataTo
// purpose  : The SurfData of <theSeq> from <theFrom>, the end's, up to
//           <theTo> removed, <theTo> kept. True.
//=======================================================================

static bool RemoveSurfDataTo(NCollection_Sequence<occ::handle<ChFiDS_SurfData>>& theSeq,
                             const int                                           theFrom,
                             const int                                           theTo)
{
  for (int i = std::max(theFrom, theTo); i > std::min(theFrom, theTo); i--)
  {
    theSeq.Remove(theFrom < theTo ? theFrom : i);
  }
  return true;
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

  // A side's line run off its face over an edge of Vtx, short of it, the
  // face narrowing to nothing at Vtx (MakeRimFlatEnd): with the flat end at
  // a vertex asked for, the end flat there -- only once everything else
  // failed (ChFi3d_MissedFlatEnd), as at an end OnSame whose line misses the
  // face at end: marked here (ChFi3d_UncappedEnd), whether or not the corner
  // below fails at it (on #876's foot it throws; on a model of it drafted,
  // it is made, and the result refused)
  if (ChFi3d_Builder::FlatEndAtVertex() && myVDataMap(Index).Extent() == 1)
  {
    int        aSens    = 0;
    const int  aNum0    = ChFi3d_IndexOfSurfData(Vtx, stripe, aSens);
    const bool isFirstR = (aSens == 1);
    // the end's SurfData past those of no face, as the corners below take it;
    // those removed only where the end is made here
    int aNum = aNum0;
    for (; (isFirstR ? aNum < SeqFil.Length() : aNum > 1)
           && (SeqFil.Value(aNum)->IndexOfS1() == 0 || SeqFil.Value(aNum)->IndexOfS2() == 0);
         aNum += isFirstR ? 1 : -1)
    {
    }
    const occ::handle<ChFiDS_SurfData> aFd = SeqFil.Value(aNum);
    const TopoDS_Edge& aEs = isFirstR ? spine->Edges(1) : spine->Edges(spine->NbEdges());
    int                aIp = 0;
    if (IsRimEnd(DStr, aFd, isFirstR, Vtx, aEs, tolapp3d, aIp))
    {
      if (!ChFi3d_MissedFlatEnd())
      {
        ChFi3d_UncappedEnd() = true;
      }
      else if (RemoveSurfDataTo(SeqFil, aNum0, aNum)
               && MakeRimFlatEnd(DStr,
                                 stripe,
                                 aFd,
                                 isFirstR,
                                 Vtx,
                                 aEs,
                                 myShape,
                                 myEFMap,
                                 tolapp3d,
                                 tolapp2d,
                                 myEVIMap))
      {
        return;
      }
    }
  }
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
  // A fillet's last surface past Vtx on the next piece of its face, kept in
  // coplanar pieces (a body without Refine): the line on that face stops
  // on the split, an edge of Vtx, and the other side's runs on a little to
  // an edge of Vtx that ends above the split's level. On the face in one
  // piece the two are one surface; the other side's end is taken over from
  // the surface past Vtx, which goes -- its line on the next piece is what
  // the corner cuts off.
  if (nbsurf != nbedge && nbsurf != 1 && num1 >= 1 && num1 <= SeqFil.Length())
  {
    const occ::handle<ChFiDS_SurfData> aSDE = SeqFil(num);
    const occ::handle<ChFiDS_SurfData> aSDM = SeqFil(num1);
    for (int ons = 1; ons <= 2; ons++)
    {
      const int ono = 3 - ons;
      if (aSDE->Index(ons) <= 0 || aSDM->Index(ons) <= 0 || aSDE->Index(ono) <= 0
          || aSDE->Index(ono) != aSDM->Index(ono))
      {
        continue;
      }
      const TopoDS_Shape& aSa = DStr.Shape(aSDM->Index(ons));
      const TopoDS_Shape& aSb = DStr.Shape(aSDE->Index(ons));
      if (aSa.ShapeType() != TopAbs_FACE || aSb.ShapeType() != TopAbs_FACE || aSa.IsSame(aSb)
          || !ChFi3d_CoplanarPieces(TopoDS::Face(aSa), TopoDS::Face(aSb), 10 * tolapp3d))
      {
        continue;
      }
      const ChFiDS_CommonPoint& aCPs = aSDM->Vertex(isfirst, ons);
      const ChFiDS_CommonPoint& aCPo = aSDE->Vertex(isfirst, ono);
      if (!aCPs.IsOnArc() || !hasVertex(aCPs.Arc(), Vtx) || !aCPo.IsOnArc()
          || !hasVertex(aCPo.Arc(), Vtx) || !containE(TopoDS::Face(aSa), aCPs.Arc())
          || !containE(TopoDS::Face(aSb), aCPs.Arc()))
      {
        continue;
      }
      // one line on the other side, of one parameter, on one surface
      const ChFiDS_FaceInterference& aFiE = aSDE->Interference(ono);
      ChFiDS_FaceInterference&       aFiM = aSDM->ChangeInterference(ono);
      if (aFiE.LineIndex() == 0 || aFiM.LineIndex() == 0 || aFiM.PCurveOnSurf().IsNull()
          || aFiM.PCurveOnFace().IsNull())
      {
        continue;
      }
      const double                    aW  = aFiE.Parameter(isfirst);
      const occ::handle<Geom_Curve>&  aCM = DStr.Curve(aFiM.LineIndex()).Curve();
      const occ::handle<Geom_Surface> aSM = DStr.Surface(aSDM->Surf()).Surface();
      if (aCM.IsNull() || aSM.IsNull() || aCM->Value(aW).Distance(aCPo.Point()) > 10 * tolapp3d)
      {
        continue;
      }
      const gp_Pnt2d aUVs = aFiM.PCurveOnSurf()->Value(aW);
      const gp_Pnt2d aUVf = aFiM.PCurveOnFace()->Value(aW);
      const BRepAdaptor_Surface aBSo(TopoDS::Face(DStr.Shape(aSDM->Index(ono))), false);
      if (aSM->Value(aUVs.X(), aUVs.Y()).Distance(aCPo.Point()) > 10 * tolapp3d
          || aBSo.Value(aUVf.X(), aUVf.Y()).Distance(aCPo.Point()) > 10 * tolapp3d)
      {
        continue;
      }
      aSDM->ChangeVertex(isfirst, ono) = aCPo;
      aFiM.SetParameter(aW, isfirst);
      SeqFil.Remove(num);
      num    = ChFi3d_IndexOfSurfData(Vtx, stripe, sens);
      num1   = isfirst ? num + 1 : num - 1;
      nbsurf = SeqFil.Length();
      break;
    }
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
  // the pieces of the wall F1 (F2) is a piece of, when it is away from Vtx
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aWall1, aWall2;
  auto anOnWall = [&](const TopoDS_Face&                                            theF,
                      const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& theWall,
                      const TopoDS_Vertex&                                          theV) {
    if (containV(theF, theV))
    {
      return true;
    }
    for (NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>::Iterator itW(theWall);
         itW.More();
         itW.Next())
    {
      if (containV(TopoDS::Face(itW.Value()), theV))
      {
        return true;
      }
    }
    return false;
  };
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
    // One side on an edge of Vtx, the other run on past it, far beyond the
    // fillet's width: the walk went on as far as the spine was extended
    // (ExtentOneCorner, half the edge), and the line on that side met no
    // edge. The faces at Vtx then cut the strip along its length rather
    // than across it, and the corner runs off the surface's end rather than
    // to a line (#273's rim under the handle, 1.6 deg from tangent: BOP's
    // check found the corner crossing itself, and its vertex 0.029 off).
    // The strip's end would want a cap in its section at Vtx, a corner the
    // builder has no form for. The corner is refused, as upstream refuses
    // the edge (PerformElement takes its faces for tangent at 0.1 rad); the
    // setback fallback may make it again set back from Vtx.
    if (CV1.IsOnArc() != CV2.IsOnArc())
    {
      const int                 aFree = CV1.IsOnArc() ? 2 : 1;
      const ChFiDS_CommonPoint& aCF   = aFree == 1 ? CV1 : CV2;
      const ChFiDS_CommonPoint& aCA   = aFree == 1 ? CV2 : CV1;
      const int                 aLI   = Fd->Interference(aFree).LineIndex();
      if (aLI != 0 && !DStr.Curve(aLI).Curve().IsNull())
      {
        // the width: from the end on the edge across to the other line
        GeomAPI_ProjectPointOnCurve aPr(aCA.Point(), DStr.Curve(aLI).Curve());
        if (aPr.NbPoints() > 0 && aPr.LowerDistance() > Precision::Confusion()
            && aCF.Point().Distance(BRep_Tool::Pnt(Vtx)) > 10. * aPr.LowerDistance())
        {
          throw Standard_Failure("IntersectionAtEnd : the fillet runs on past the vertex");
        }
      }
    }
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

      // The fillet's face away from Vtx, a piece of the wall the spine is
      // on (a wall kept in coplanar pieces, the fillet wider than the piece
      // beside its spine): the way from Vtx to the arc is the wall's
      // outline, over its pieces, and the faces at end are those across it,
      // as on the wall in one face.
      auto aWallOf = [&](const TopoDS_Face&                                      theF,
                         NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& theWall) {
        if (containV(theF, Vtx))
        {
          return;
        }
        NCollection_List<TopoDS_Shape> aList;
        aList.Append(theF);
        theWall.Add(theF);
        bool hasSpine = false;
        for (NCollection_List<TopoDS_Shape>::Iterator itL(aList); itL.More(); itL.Next())
        {
          const TopoDS_Face aP = TopoDS::Face(itL.Value());
          hasSpine             = hasSpine || containE(aP, EdgeSpine);
          for (TopExp_Explorer exE(aP, TopAbs_EDGE); exE.More(); exE.Next())
          {
            if (!myEFMap.Contains(exE.Current()) || myEFMap(exE.Current()).Extent() != 2)
            {
              continue;
            }
            for (NCollection_List<TopoDS_Shape>::Iterator itF(myEFMap(exE.Current())); itF.More();
                 itF.Next())
            {
              const TopoDS_Face& aFn = TopoDS::Face(itF.Value());
              if (!theWall.Contains(aFn) && ChFi3d_CoplanarPieces(aP, aFn, 10 * tolapp3d))
              {
                theWall.Add(aFn);
                aList.Append(aFn);
              }
            }
          }
        }
        if (!hasSpine)
        {
          theWall.Clear();
        }
      };
      aWallOf(F1, aWall1);
      aWallOf(F2, aWall2);
      // the number of the wall's pieces an edge bounds
      auto aNbIn = [&](const TopoDS_Shape&                                           theE,
                       const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& theWall) {
        int n = 0;
        for (NCollection_List<TopoDS_Shape>::Iterator itF(myEFMap(theE)); itF.More(); itF.Next())
        {
          if (theWall.Contains(itF.Value()))
          {
            n++;
          }
        }
        return n;
      };
      // the next edge of the face's outline (or the wall's) from theV, after
      // theE, and its other vertex
      auto aNext = [&](const TopoDS_Vertex&                                          theV,
                       const TopoDS_Edge&                                            theE,
                       const TopoDS_Face&                                            theF,
                       const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& theWall,
                       TopoDS_Edge&                                                  theNext,
                       TopoDS_Vertex&                                                theVNext) {
        if (theWall.IsEmpty())
        {
          ChFi3d_cherche_element(theV, theE, theF, theNext, theVNext);
          return;
        }
        TopoDS_Edge aFound;
        for (NCollection_List<TopoDS_Shape>::Iterator itE(myVEMap(theV)); itE.More(); itE.Next())
        {
          const TopoDS_Edge& anE = TopoDS::Edge(itE.Value());
          if (anE.IsSame(theE) || anE.IsSame(aFound) || BRep_Tool::Degenerated(anE)
              || aNbIn(anE, theWall) != 1)
          {
            continue;
          }
          if (!aFound.IsNull())
          {
            throw Standard_Failure("IntersectionAtEnd : the wall's outline forks");
          }
          aFound = anE;
        }
        TopoDS_Vertex aVa, aVb;
        if (!aFound.IsNull())
        {
          TopExp::Vertices(aFound, aVa, aVb);
        }
        if (aFound.IsNull() || aVa.IsSame(aVb))
        {
          throw Standard_Failure("IntersectionAtEnd : the wall's outline ends");
        }
        theNext  = aFound;
        theVNext = theV.IsSame(aVa) ? aVb : aVa;
      };
      // the face across theE from the face (or from the wall)
      auto anAcross = [&](const TopoDS_Edge&                                            theE,
                          const TopoDS_Face&                                            theF,
                          const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& theWall,
                          TopoDS_Face&                                                  theAcross) {
        if (theWall.IsEmpty())
        {
          ChFi3d_cherche_face1(myEFMap(theE), theF, theAcross);
          return;
        }
        for (NCollection_List<TopoDS_Shape>::Iterator itF(myEFMap(theE)); itF.More(); itF.Next())
        {
          if (!theWall.Contains(itF.Value()))
          {
            theAcross = TopoDS::Face(itF.Value());
            return;
          }
        }
        throw Standard_Failure("IntersectionAtEnd : no face across the wall's outline");
      };

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
          aNext(V1, Ei, F1, aWall1, E, V2);
          V1 = V2;
          Ei = E;
          anAcross(E, F1, aWall1, Face[nb]);
          // (two edges of the outline in line, a split's foot between them)
          if (!aWall1.IsEmpty() && Face[nb].IsSame(Face[nb - 1]))
          {
            continue;
          }
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
        aNext(V1, EdgeSpine, F2, aWall2, E, V2);
        Ei = E;
        V1 = V2;
        while (!V1.IsSame(Vfin))
        {
          aNext(V1, Ei, F2, aWall2, E, V2);
          Ei = E;
          V1 = V2;
          anAcross(E, F2, aWall2, Face[nb]);
          if (!aWall2.IsEmpty() && Face[nb].IsSame(Face[nb - 1]))
          {
            continue;
          }
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
      // One side on an edge of Vtx, the other running on under a wall that
      // stands on the face at end: the fillet ends in that face's plane,
      // which it cuts as far as the wall, and under the wall the body is
      // closed by a ceiling in the same plane (CeilingAtEnd).
      if (!couture && !bordlibre && !isOnSame1 && !isOnSame2)
      {
        for (int ip = 1; ip <= 2; ip++)
        {
          if ((ip == 1 ? possible1 : possible2)
              && CeilingAtEnd(DStr,
                              stripe,
                              Fd,
                              isfirst,
                              ip,
                              Vtx,
                              EdgeSpine,
                              HGs,
                              orcourbe,
                              myVEMap,
                              myEFMap,
                              myVFMap,
                              tolapp3d,
                              tol2d,
                              myEVIMap))
          {
            return;
          }
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
      // Neither side on an edge of Vtx, one of them on an edge of its face
      // past the foot of a split of the face at end -- a wall kept in
      // coplanar pieces, the fillet wider than the piece at Vtx. On the wall
      // in one face that side ends on an edge of Vtx and the other runs on
      // in a face that goes on past Vtx (IntersectMoreCorner, below); here
      // the cut runs over the two pieces, which PerformOneCorner makes.
      if (!possible1 && !possible2 && nbarete == 4 && !couture && !bordlibre && onsOS == 0)
      {
        for (int ip = 1; ip <= 2; ip++)
        {
          const int                      io  = 3 - ip;
          const ChFiDS_FaceInterference& Fio = Fd->Interference(io);
          const double                   w   = Fd->Interference(ip).Parameter(isfirst);
          TopoDS_Edge                    anEat;
          if (!ArcPastWallSplit(Vtx,
                                Fd->Vertex(isfirst, ip).Arc(),
                                ip == 1 ? F1 : F2,
                                myVEMap,
                                myEFMap,
                                10 * tolapp3d,
                                anEat))
          {
            continue;
          }
          // (not where the section crosses the face at end, see below)
          if (w >= Fio.FirstParameter() && w <= Fio.LastParameter())
          {
            const gp_Pnt2d uv = Fio.PCurveOnSurf()->Value(w);
            const gp_Pnt   Po = DStr.Surface(Fd->Surf()).Surface()->Value(uv.X(), uv.Y());
            if (SectionCrossesEndFace(Vtx,
                                      anEat,
                                      ip == 1 ? F1 : F2,
                                      Fd->Vertex(isfirst, ip).Point(),
                                      Po,
                                      myVEMap,
                                      myEFMap))
            {
              continue;
            }
          }
          PerformOneCorner(Index, false, true);
          return;
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
            // (No ceiling was made, above.) The line is cut there, as an
            // edge splitting its face there would cut it, and the corner is
            // filled as that of the split face, with a plate.
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
        // (tangent faces, not the pieces of one: a split of a face kept in
        // coplanar pieces is no edge, and the vertex one of four edges --
        // the case of the pieces, below)
        bool isTangent0 = ChFi3d::IsTangentFaces(Edge[0], F1, Face[0])
                          && !ChFi3d_CoplanarPieces(F1, Face[0], 10 * tolapp3d);
        bool isTangent1 = ChFi3d::IsTangentFaces(Edge[1], Face[0], Face[1])
                          && !ChFi3d_CoplanarPieces(Face[0], Face[1], 10 * tolapp3d);
        bool isTangent2 = ChFi3d::IsTangentFaces(Edge[2], Face[1], Face[2])
                          && !ChFi3d_CoplanarPieces(Face[1], Face[2], 10 * tolapp3d);
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
      // The same two edges in line, with the face between them kept in
      // coplanar pieces (a body without Refine): the pieces are one face
      // of the four-edge case above, and the fillet does not cut them. The
      // piece the line's own arc borders takes the edge from Vtx to the
      // line's end; the pieces between, and the splits, stay as they are.
      // Each piece cut for itself gave the section of the fillet in their
      // plane twice, there and back.
      if (!oneintersection1 && !oneintersection2 && nbface >= 3 && !couture && !bordlibre)
      {
        auto arePieces = [&](const int theFrom, const int theTo) {
          for (int k = theFrom; k < theTo; k++)
          {
            if (!ChFi3d_CoplanarPieces(Face[k], Face[k + 1], 10 * tolapp3d))
            {
              return false;
            }
          }
          return true;
        };
        auto inLine = [&](const TopoDS_Edge& theE1, const TopoDS_Edge& theE2) {
          return hasVertex(theE1, Vtx) && hasVertex(theE2, Vtx)
                 && std::abs(ChFi3d_AngleEdge(Vtx, theE1, theE2) - M_PI) < 0.01;
        };
        if (arePieces(1, nbface - 1) && inLine(Edge[1], Edge[nbface]))
        {
          oneintersection2 = true;
          facesau          = Face[nbface - 1];
          edgesau          = Edge[1];
          Edge[1]          = Edge[nbface];
          nbface           = 1;
        }
        else if (arePieces(0, nbface - 2) && inLine(Edge[0], Edge[nbface - 1]))
        {
          oneintersection1 = true;
          facesau          = Face[0];
          edgesau          = Edge[nbface - 1];
          Face[0]          = Face[nbface - 1];
          Edge[1]          = Edge[nbface];
          nbface           = 1;
        }
      }
      // A side of the fillet that ends on a split of its own face, kept in
      // coplanar pieces (a body without Refine), the split an edge of Vtx
      // lying on the surface of the face at end beside the other side -- its
      // plane, or the straight end of a ramp: the fillet's end is its
      // section by that surface, from the one side's edge to the split, and
      // that face alone is cut, as on the wall in one face
      // (IntersectMoreCorner there, the line running on past Vtx in its
      // face). The face grows over the fillet's end, up to the split from
      // Vtx to the line's end, which it shares with the piece beyond from
      // then on; the faces between are left as they are.
      // Not where a wall stands on that face in the way of the section
      // (SectionCrossesEndFace): the ceiling's above, or, where that was
      // not made, the plate's, there as on the wall in one face.
      if (!oneintersection1 && !oneintersection2 && nbface >= 2 && !couture && !bordlibre)
      {
        auto onSurfaceOf = [&](const TopoDS_Edge& theE, const TopoDS_Face& theF) {
          const BRepAdaptor_Surface aS(theF, false);
          const BRepAdaptor_Curve   aC(theE);
          for (int k = 0; k <= 2; k++)
          {
            const double aW = aC.FirstParameter()
                              + k * (aC.LastParameter() - aC.FirstParameter()) / 2.;
            const gp_Pnt aP = aC.Value(aW);
            if (aS.GetType() == GeomAbs_Plane)
            {
              if (aS.Plane().Distance(aP) > 10 * tolapp3d)
              {
                return false;
              }
              continue;
            }
            GeomAPI_ProjectPointOnSurf aProj(aP, aS.Surface().Surface());
            if (!aProj.IsDone() || aProj.NbPoints() == 0 || aProj.LowerDistance() > 10 * tolapp3d)
            {
              return false;
            }
          }
          return true;
        };
        if (hasVertex(Edge[nbface], Vtx)
            && ChFi3d_CoplanarPieces(F2, Face[nbface - 1], 10 * tolapp3d)
            && onSurfaceOf(Edge[nbface], Face[0])
            && !SectionCrossesEndFace(Vtx, Edge[0], F1, CV1.Point(), CV2.Point(), myVEMap, myEFMap))
        {
          oneintersection2 = true;
          facesau          = Face[nbface - 1];
          edgesau          = Edge[nbface];
          Edge[1]          = Edge[nbface];
          nbface           = 1;
        }
        else if (hasVertex(Edge[0], Vtx) && ChFi3d_CoplanarPieces(F1, Face[0], 10 * tolapp3d)
                 && onSurfaceOf(Edge[0], Face[nbface - 1])
                 && !SectionCrossesEndFace(Vtx,
                                           Edge[nbface],
                                           F2,
                                           CV2.Point(),
                                           CV1.Point(),
                                           myVEMap,
                                           myEFMap))
        {
          oneintersection1 = true;
          facesau          = Face[0];
          edgesau          = Edge[0];
          Face[0]          = Face[nbface - 1];
          Edge[1]          = Edge[nbface];
          nbface           = 1;
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
          else
          {
            // a surface that is not extended (a cylinder, a plane): the
            // point is projected on it as it is
            PointOnSurfaceAsItIs(CV1.Point(), Sfacemoins1, BRepAdaptor_Surface(F, true), pfac1);
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
      // A long curve -- the edge of a helical ramp, turn upon turn of it --
      // is sampled too coarsely to show the fillet crossing it close to
      // Vtx: once more on the piece of it about Vtx. Without the point the
      // edge was taken for one the fillet misses, the faces beside it were
      // intersected anew and the piece from Vtx added on top of the edge.
      if ((!inters.IsDone() || inters.NbPoints() == 0) && hasVertex(E2, Vtx))
      {
        const double aLen  = 10. * CV1.Point().Distance(CV2.Point());
        const double aDw   = Cad.Resolution(aLen);
        const double aParV = BRep_Tool::Parameter(Vtx, E2);
        const double aW1 = std::max(Ubid, aParV - aDw), aW2 = std::min(Vbid, aParV + aDw);
        if (aLen > Precision::Confusion() && aW2 - aW1 > Precision::PConfusion()
            && aW2 - aW1 < Vbid - Ubid)
        {
          occ::handle<GeomAdaptor_Curve> aHCloc = new GeomAdaptor_Curve(C, aW1, aW2);
          inters.Perform(aHCloc, HGs);
        }
      }
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
          // On a closed curve the point comes on the curve's first turn,
          // whatever the edge's range: an arc of a circle from 3.95 to 9.10
          // cut at 9.097 was given the point at 2.814, off the edge, and the
          // faces beside it lost the edge whole.
          if (!extend && Cad.IsPeriodic())
          {
            double aFE, aLE;
            BRep_Tool::Range(E2, aFE, aLE);
            const double aPeriod = Cad.Period();
            paredge2 += aPeriod * std::floor((0.5 * (aFE + aLE) - paredge2) / aPeriod + 0.5);
          }
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
        else
        {
          // (see above; pfac2 was left as it stood, (0,0) on the first face)
          PointOnSurfaceAsItIs(CV2.Point(), Sfacemoins1, BRepAdaptor_Surface(F, true), pfac2);
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

    // A face of a B-spline or a Bezier -- a loft's plane among them -- ends
    // at its bounds, where a plane runs on: the fillet's cut past them (the
    // end beyond the face's edge, on its line) is on the surface extended,
    // as it is on a plane; bounded, it was not found and the end was a plate
    // (#631: a loft on a pad, the fillet of their joint at 4 bulging out of
    // the pad's side)
    if (!prolface[nb - 1] && !Sfacemoins1->IsUPeriodic() && !Sfacemoins1->IsVPeriodic()
        && (Sfacemoins1->IsKind(STANDARD_TYPE(Geom_BSplineSurface))
            || Sfacemoins1->IsKind(STANDARD_TYPE(Geom_BezierSurface))))
    {
      double aU1, aU2, aV1, aV2;
      Sfacemoins1->Bounds(aU1, aU2, aV1, aV2);
      const double aTolU = Precision::PConfusion() * std::max(1., aU2 - aU1);
      const double aTolV = Precision::PConfusion() * std::max(1., aV2 - aV1);
      const auto   isOut = [&](const gp_Pnt2d& theP) {
        return theP.X() < aU1 - aTolU || theP.X() > aU2 + aTolU || theP.Y() < aV1 - aTolV
               || theP.Y() > aV2 + aTolV;
      };
      if (isOut(pfac1) || isOut(pfac2))
      {
        ChFi3d_ExtendSurface(Sfacemoins1, prolface[nb - 1]);
        if (prolface[nb - 1])
        {
          BRep_Builder BRE;
          BRE.MakeFace(faceprol[nb - 1], Sfacemoins1, F.Location(), BRep_Tool::Tolerance(F));
        }
      }
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
        // (on F1 or F2, or on the outline of the wall it is a piece of)
        const bool on1V1 = anOnWall(F1, aWall1, V1), on1V2 = anOnWall(F1, aWall1, V2);
        const bool on2V1 = anOnWall(F2, aWall2, V1), on2V2 = anOnWall(F2, aWall2, V2);
        if (on1V1 || on2V1)
        {
          ori = TopAbs_FORWARD;
        }
        else if (on1V2 || on2V2)
        {
          ori = TopAbs_REVERSED;
        }
        else
        {
          throw Standard_Failure("IntersectionAtEnd : pb orientation");
        }

        // both ends on the fillet's faces -- the same one, or one each (a
        // face beside the vertex narrower than the fillet: #474's sliver of
        // a top, its edge on F2 0.155 long and the cut edge's other end on
        // F1 2 away) -- the end nearer the vertex is the one under the fillet
        if ((on1V1 || on2V1) && (on1V2 || on2V2))
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
        // The nearest of the curve's points to pext, not the first found: a
        // curve winding round an axis (the edge of a helical ramp) has one
        // on every turn, and the first can be a turn away from Vtx -- the
        // piece kept was then the rest of the curve.
        auto nearestOn = [&](const gp_Pnt& theP) {
          Extrema_ExtPC anExt(theP, cad, tolpt);
          if (!anExt.IsDone() || anExt.NbExt() == 0)
          {
            throw Standard_Failure("IntersectionAtEnd : no point of the edge near the cut");
          }
          int iMin = 1;
          for (int i = 2; i <= anExt.NbExt(); i++)
          {
            if (anExt.SquareDistance(i) < anExt.SquareDistance(iMin))
            {
              iMin = i;
            }
          }
          return anExt.Point(iMin).Parameter();
        };
        double par1, par2, par, ParVtx;
        bool   vtx1 = false;
        bool   vtx2 = false;
        par1        = nearestOn(pext);
        ParVtx      = par1;
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
          par2 = nearestOn(pext);
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
        // The piece follows the cut on from the end they share: by the
        // nearer pairing of ends, not the first within tolpt. The cut,
        // computed on a face carried on, can end farther than that from the
        // point on the edge (#273's handle at the cup's rim, 0.003 off), and
        // then neither test held and the piece was turned regardless.
        if (std::min(P2.Distance(P3), P1.Distance(P4))
            <= std::min(P1.Distance(P3), P2.Distance(P4)))
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
            // (the split of the fillet's own face, on Face[0]'s surface)
            if (C2dint1.IsNull())
            {
              ComputeCurve2d(Ct, Face[0], C2dint1);
            }
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
    indSurface = aStripe->SetOfSurfData()->Value(anInd)->IndexOfS1();
  }
  else
  {
    indSurface = aStripe->SetOfSurfData()->Value(anInd)->IndexOfS2();
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
