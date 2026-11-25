// Created on: 2004-01-09
// Created by: Sergey KUUL
// Copyright (c) 2004-2014 OPEN CASCADE SAS
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

#ifndef _XCAFDoc_DimTolTool_HeaderFile
#define _XCAFDoc_DimTolTool_HeaderFile

#include <Standard.hxx>
#include <Standard_Type.hxx>

#include <TDataStd_GenericEmpty.hxx>
#include <Standard_Boolean.hxx>
#include <TDF_Label.hxx>
#include <NCollection_Sequence.hxx>
#include <Standard_Integer.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_HArray1.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Ax3.hxx>
#include <TCollection_AsciiString.hxx>
#include <NCollection_IndexedDataMap.hxx>

class XCAFDoc_ShapeTool;
class TDF_Label;
class Standard_GUID;
class TCollection_HAsciiString;
class TCollection_AsciiString;
class Quantity_Color;

//! Attribute containing GD&T section of XCAF document.
//! Provide tools for GD&T section management.
class XCAFDoc_DimTolTool : public TDataStd_GenericEmpty
{

public:
  Standard_EXPORT XCAFDoc_DimTolTool();

  //! Creates (if not exist) DimTolTool attribute.
  Standard_EXPORT static occ::handle<XCAFDoc_DimTolTool> Set(const TDF_Label& L);

  //! Returns the standard GD&T tool GUID.
  Standard_EXPORT static const Standard_GUID& GetID();

  //! Returns the label under which GD&T table is stored.
  Standard_EXPORT TDF_Label BaseLabel() const;

  //! Returns internal XCAFDoc_ShapeTool tool
  Standard_EXPORT const occ::handle<XCAFDoc_ShapeTool>& ShapeTool();

  //! Returns True if the label belongs to a GD&T table and
  //! is a Dimension definition.
  Standard_EXPORT bool IsDimension(const TDF_Label& theLab) const;

  //! Returns a sequence of Dimension labels currently stored
  //! in the GD&T table.
  Standard_EXPORT void GetDimensionLabels(NCollection_Sequence<TDF_Label>& theLabels) const;

  //! Sets a dimension to sequences target labels.
  Standard_EXPORT void SetDimension(const NCollection_Sequence<TDF_Label>& theFirstLS,
                                    const NCollection_Sequence<TDF_Label>& theSecondLS,
                                    const TDF_Label&                       theDimL) const;

  //! Sets a dimension to target labels.
  Standard_EXPORT void SetDimension(const TDF_Label& theFirstL,
                                    const TDF_Label& theSecondL,
                                    const TDF_Label& theDimL) const;

  //! Sets a dimension to the target label.
  Standard_EXPORT void SetDimension(const TDF_Label& theL, const TDF_Label& theDimL) const;

  //! Returns all Dimension labels defined for theShapeL.
  Standard_EXPORT bool GetRefDimensionLabels(const TDF_Label&                 theShapeL,
                                             NCollection_Sequence<TDF_Label>& theDimensions) const;

  //! Adds a dimension definition to the GD&T table and returns its label.
  Standard_EXPORT TDF_Label AddDimension();

  //! Returns True if the label belongs to the GD&T table and is a dimension tolerance.
  Standard_EXPORT bool IsGeomTolerance(const TDF_Label& theLab) const;

  //! Returns a sequence of Tolerance labels currently stored in the GD&T table.
  Standard_EXPORT void GetGeomToleranceLabels(NCollection_Sequence<TDF_Label>& theLabels) const;

  //! Sets a geometry tolerance from theGeomTolL to theL label.
  //! Checks if theGeomTolL is a geometry tolerance definition first.
  Standard_EXPORT void SetGeomTolerance(const TDF_Label& theL, const TDF_Label& theGeomTolL) const;

  //! Sets a geometry tolerance from theGeomTolL to sequence of labels theL.
  //! Checks if theGeomTolL is a geometry tolerance definition first.
  Standard_EXPORT void SetGeomTolerance(const NCollection_Sequence<TDF_Label>& theL,
                                        const TDF_Label&                       theGeomTolL) const;

  //! Returns all GeomTolerance labels defined for theShapeL.
  Standard_EXPORT bool GetRefGeomToleranceLabels(const TDF_Label&                 theShapeL,
                                                 NCollection_Sequence<TDF_Label>& theDimTols) const;

  //! Adds a GeomTolerance definition to the GD&T table and returns its label.
  Standard_EXPORT TDF_Label AddGeomTolerance();

  //! Returns True if theLab belongs to the GD&T table and is a dmension tolerance.
  Standard_EXPORT bool IsDimTol(const TDF_Label& theLab) const;

  //! Returns a sequence of D&GTs currently stored in the GD&T table.
  Standard_EXPORT void GetDimTolLabels(NCollection_Sequence<TDF_Label>& Labels) const;

  //! Finds a dimension tolerance definition in the GD&T table
  //! satisfying the specified kind, values, name and description
  //! and returns its label if found.
  //! Returns False if dimension tolerance is not found in DGTtable.
  Standard_EXPORT bool FindDimTol(const int                                       theKind,
                                  const occ::handle<NCollection_HArray1<double>>& theVal,
                                  const occ::handle<TCollection_HAsciiString>&    theName,
                                  const occ::handle<TCollection_HAsciiString>&    theDescription,
                                  TDF_Label&                                      lab) const;

  //! Finds a dimension tolerance in the GD&T table
  //! satisfying the specified kind, values, name and description
  //! and returns its label if found (or Null label else).
  Standard_EXPORT TDF_Label
    FindDimTol(const int                                       theKind,
               const occ::handle<NCollection_HArray1<double>>& theVal,
               const occ::handle<TCollection_HAsciiString>&    theName,
               const occ::handle<TCollection_HAsciiString>&    theDescription) const;

  //! Adds a dimension tolerance definition with the specified
  //! kind, value, name and description to the GD&T table and returns its label.
  Standard_EXPORT TDF_Label
    AddDimTol(const int                                       theKind,
              const occ::handle<NCollection_HArray1<double>>& theVal,
              const occ::handle<TCollection_HAsciiString>&    theName,
              const occ::handle<TCollection_HAsciiString>&    theDescription) const;

  //! Sets existing dimension tolerance to theL label.
  Standard_EXPORT void SetDimTol(const TDF_Label& theL, const TDF_Label& theDimTolL) const;

  //! Creates a dimension tolerance and sets it to theL label.
  Standard_EXPORT TDF_Label
    SetDimTol(const TDF_Label&                                theL,
              const int                                       theKind,
              const occ::handle<NCollection_HArray1<double>>& theVal,
              const occ::handle<TCollection_HAsciiString>&    theName,
              const occ::handle<TCollection_HAsciiString>&    theDescription) const;

  //! Gets all shape labels referred by theL label of the GD&T table.
  //! Returns False if there are no shape labels added to the sequences.
  Standard_EXPORT static bool GetRefShapeLabel(const TDF_Label&                 theL,
                                               NCollection_Sequence<TDF_Label>& theShapeLFirst,
                                               NCollection_Sequence<TDF_Label>& theShapeLSecond);

  //! Returns dimension tolerance assigned to theDimTolL label.
  //! Returns False if no such dimension tolerance is assigned.
  Standard_EXPORT bool GetDimTol(const TDF_Label&                          theDimTolL,
                                 int&                                      theKind,
                                 occ::handle<NCollection_HArray1<double>>& theVal,
                                 occ::handle<TCollection_HAsciiString>&    theName,
                                 occ::handle<TCollection_HAsciiString>&    theDescription) const;

  //! Returns True if label belongs to the GD&T table and
  //! is a Datum definition.
  Standard_EXPORT bool IsDatum(const TDF_Label& lab) const;

  //! Returns a sequence of Datums currently stored
  //! in the GD&T table.
  Standard_EXPORT void GetDatumLabels(NCollection_Sequence<TDF_Label>& Labels) const;

  //! Finds a datum satisfying the specified name, description and
  //! identification and returns its label if found.
  Standard_EXPORT bool FindDatum(const occ::handle<TCollection_HAsciiString>& theName,
                                 const occ::handle<TCollection_HAsciiString>& theDescription,
                                 const occ::handle<TCollection_HAsciiString>& theIdentification,
                                 TDF_Label&                                   lab) const;

  //! Adds a datum definition to the GD&T table and returns its label.
  Standard_EXPORT TDF_Label
    AddDatum(const occ::handle<TCollection_HAsciiString>& theName,
             const occ::handle<TCollection_HAsciiString>& theDescription,
             const occ::handle<TCollection_HAsciiString>& theIdentification) const;

  //! Adds a datum definition to the GD&T table and returns its label.
  Standard_EXPORT TDF_Label AddDatum();

  //! Sets a datum to the sequence of shape labels.
  Standard_EXPORT void SetDatum(const NCollection_Sequence<TDF_Label>& theShapeLabels,
                                const TDF_Label&                       theDatumL) const;

  //! Sets a datum to theL label and binds it with theTolerL label.
  //! A datum with the specified name, description and identification
  //! is created if it isn't found in the GD&T table.
  Standard_EXPORT void SetDatum(
    const TDF_Label&                             theL,
    const TDF_Label&                             theTolerL,
    const occ::handle<TCollection_HAsciiString>& theName,
    const occ::handle<TCollection_HAsciiString>& theDescription,
    const occ::handle<TCollection_HAsciiString>& theIdentification) const;

  //! Sets a datum from theDatumL label to theToletL label.
  Standard_EXPORT void SetDatumToGeomTol(const TDF_Label& theDatumL,
                                         const TDF_Label& theTolerL) const;

  //! Returns datum assigned to theDatumL label.
  //! Returns False if no such datum is assigned.
  Standard_EXPORT bool GetDatum(const TDF_Label&                       theDatumL,
                                occ::handle<TCollection_HAsciiString>& theName,
                                occ::handle<TCollection_HAsciiString>& theDescription,
                                occ::handle<TCollection_HAsciiString>& theIdentification) const;

  //! Returns all Datum labels defined for theDimTolL label.
  Standard_EXPORT static bool GetDatumOfTolerLabels(const TDF_Label&                 theDimTolL,
                                                    NCollection_Sequence<TDF_Label>& theDatums);

  //! Returns all Datum labels with XCAFDimTolObjects_DatumObject defined for label theDimTolL.
  Standard_EXPORT static bool GetDatumWithObjectOfTolerLabels(
    const TDF_Label&                 theDimTolL,
    NCollection_Sequence<TDF_Label>& theDatums);

  //! Returns all GeomToleranses labels defined for theDatumL label.
  Standard_EXPORT bool GetTolerOfDatumLabels(const TDF_Label&                 theDatumL,
                                             NCollection_Sequence<TDF_Label>& theTols) const;

  //! Makes the tolerance theLowerL a lower segment of the composite frame
  //! whose upper segment is theUpperL (AP242's geometric_tolerance_relationship
  //! 'composite': one characteristic, the lower segments refining the upper's
  //! pattern-locating zone). A lower segment is of one frame: an earlier
  //! upper is replaced.
  Standard_EXPORT void SetCompositeUpper(const TDF_Label& theLowerL,
                                         const TDF_Label& theUpperL) const;

  //! The upper segment of the composite frame the tolerance theLowerL is a
  //! lower segment of; false when it is none.
  Standard_EXPORT static bool GetCompositeUpper(const TDF_Label& theLowerL, TDF_Label& theUpperL);

  //! Says the GD&T theGDTL, which holds no graphic of its own, is presented by
  //! the graphic the GD&T theHostL holds (a hole's depth in its callout, a
  //! pattern's sizes in one "32X R"): written to STEP, it is linked to that
  //! graphic as what it presents. One host each: an earlier one is replaced.
  Standard_EXPORT void SetPresentedWith(const TDF_Label& theGDTL, const TDF_Label& theHostL) const;

  //! The GD&T whose graphic presents theGDTL too; false when there is none.
  Standard_EXPORT static bool GetPresentedWith(const TDF_Label& theGDTL, TDF_Label& theHostL);

  //! Says the GD&T theGDTL has an annotation placeholder (AP242 ed2: where a
  //! reader that draws PMI itself puts its text) of the kind theKind, the name
  //! of its geometric set: 'linear dimension', 'general tolerance', 'note'.
  //! theName is the placeholder's own name, which need not be its graphic's
  //! (NIST STC-07: a dimension drawn as 'Oriented Linear Dimension.1' has the
  //! placeholder 'Oriented Linear Dimension.5'); empty, it is named as the graphic.
  //! Written to STEP it is given one, where its text is. An empty kind removes it.
  //! theShared says it is the placeholder of another PMI written before it, of its kind, name
  //! and place (NIST STC-06: a hole callout's one placeholder, its diameter's and its depth's):
  //! written, it is linked to that one, not given one of its own.
  Standard_EXPORT static void SetPlaceholder(const TDF_Label&               theGDTL,
                                             const TCollection_AsciiString& theKind,
                                             const TCollection_AsciiString& theName,
                                             const bool                     theShared = false);

  //! The kind and name of the placeholder theGDTL has; false when it has none.
  Standard_EXPORT static bool GetPlaceholder(const TDF_Label&         theGDTL,
                                             TCollection_AsciiString& theKind,
                                             TCollection_AsciiString& theName);

  //! As above, and whether it is another's (see SetPlaceholder).
  Standard_EXPORT static bool GetPlaceholder(const TDF_Label&         theGDTL,
                                             TCollection_AsciiString& theKind,
                                             TCollection_AsciiString& theName,
                                             bool&                    theShared);

  //! Adds a group of supplemental geometry of the part theShapeL (STEP AP242: a
  //! CONSTRUCTIVE_GEOMETRY_REPRESENTATION related to the part's shape representation):
  //! reference and construction geometry beside the part's shape -- planes, cylinders, axis
  //! systems, sets of curves -- that defines nothing of it, though GD&T may be on its curves
  //! (NIST CTC-04: a cosmetic thread). theName is the representation's name, theRelationName
  //! its relationship's. Its items are added with AddSupplementalItem.
  Standard_EXPORT TDF_Label AddSupplementalGeometry(
    const TDF_Label&               theShapeL,
    const TCollection_AsciiString& theName,
    const TCollection_AsciiString& theRelationName) const;

  //! The groups of supplemental geometry, in the order they were added.
  Standard_EXPORT void GetSupplementalGeometry(NCollection_Sequence<TDF_Label>& theGroups) const;

  //! The part, name and relationship name of the group theGroupL; false when it is none.
  Standard_EXPORT static bool GetSupplementalGeometry(const TDF_Label&         theGroupL,
                                                      TDF_Label&               theShapeL,
                                                      TCollection_AsciiString& theName,
                                                      TCollection_AsciiString& theRelationName);

  //! Adds an item to the group theGroupL. theKind is "plane", "cylinder" (thePosition its
  //! axis, theRadius its radius), "axes" (an axis system, thePosition) or "curves" (a set of
  //! curves and points, theShape; a shape GD&T is on). Returns its label, a null one when
  //! theGroupL is no group.
  Standard_EXPORT static TDF_Label AddSupplementalItem(const TDF_Label&               theGroupL,
                                                       const TCollection_AsciiString& theKind,
                                                       const TCollection_AsciiString& theName,
                                                       const gp_Ax3&       thePosition,
                                                       const double        theRadius,
                                                       const TopoDS_Shape& theShape);

  //! The items of the group theGroupL, in the order they were added.
  Standard_EXPORT static void GetSupplementalItems(const TDF_Label&                 theGroupL,
                                                   NCollection_Sequence<TDF_Label>& theItems);

  //! The kind, name, position, radius and shape of the item theItemL; false when it is none.
  Standard_EXPORT static bool GetSupplementalItem(const TDF_Label&         theItemL,
                                                  TCollection_AsciiString& theKind,
                                                  TCollection_AsciiString& theName,
                                                  gp_Ax3&                  thePosition,
                                                  double&                  theRadius,
                                                  TopoDS_Shape&            theShape);

  //! Whether theItemL is an item of supplemental geometry: GD&T on its curves refers to it.
  Standard_EXPORT static bool IsSupplementalItem(const TDF_Label& theItemL);

  //! Says the view theViewL shows the items of supplemental geometry theItems, in that order
  //! (STEP: a 'supplemental geometry subset' in its draughting model -- NIST STC-06's views each
  //! show some of the part's planes and curves). An empty list says nothing.
  Standard_EXPORT static void SetViewSupplemental(const TDF_Label&                       theViewL,
                                                  const NCollection_Sequence<TDF_Label>& theItems);

  //! The items of supplemental geometry the view theViewL shows; false when it says none.
  Standard_EXPORT static bool GetViewSupplemental(const TDF_Label&                 theViewL,
                                                  NCollection_Sequence<TDF_Label>& theItems);

  //! Says how the shapes of the GD&T theGDTL are grouped, its first (theIsSecond false) or
  //! its second ones: STEP's composite shape aspects of them, which a pattern of features and a
  //! feature of several faces are given as (NIST: a composite group shape aspect 'pattern of
  //! features' of a composite shape aspect 'unified geometric item' a hole). theStructure is
  //! a node, after the composite's identity, if it says one:
  //!   structure := ['@' digits] node
  //!   node      := kind '{' name '|' description '}' '(' member (',' member)* ')'
  //!   kind      := 'C' | 'G' | 'D' | 'F' | 'H'
  //!   member    := node | entry ('+' entry)*
  //! 'C' a composite shape aspect, 'G' a composite group shape aspect, 'D' -- the node said
  //! first only -- a datum feature, 'F' one that is a composite group shape aspect too (its own
  //! shapes a member, then its members): on a GD&T, the GD&T on the feature itself, written on
  //! the one written for the datum on those shapes; on a datum (theGDTL a datum's label), how
  //! its feature is written, and named. 'H' a member that is
  //! a hole occurrence, its faces (XCAFDoc_DimTolTool::AddHoleOccurrence); an entry is one of the
  //! GD&T's shapes, by its label's entry (TDF_Tool::Entry); several, one shape aspect of them.
  //! '%', '{', '|', '}', '(', ')', ',', '+' in a name or description are written %XX. GD&T of
  //! one structure, its identity too, are on one composite (STEP: the reader's identity is the
  //! entity's number, NIST CTC-04's position and size of one pattern of holes share theirs,
  //! CTC-02's two alike do not). An empty structure says nothing.
  Standard_EXPORT static void SetAspectStructure(const TDF_Label&               theGDTL,
                                                 const bool                     theIsSecond,
                                                 const TCollection_AsciiString& theStructure);

  //! How the shapes of the GD&T theGDTL are grouped (SetAspectStructure); false when it does
  //! not say.
  Standard_EXPORT static bool GetAspectStructure(const TDF_Label&         theGDTL,
                                                 const bool               theIsSecond,
                                                 TCollection_AsciiString& theStructure);

  //! Gives the item theItemL the colour it is drawn in (NIST's supplemental geometry is
  //! white or black).
  Standard_EXPORT static void SetSupplementalColor(const TDF_Label&      theItemL,
                                                   const Quantity_Color& theColor);

  //! The colour of the item theItemL; false when it has none.
  Standard_EXPORT static bool GetSupplementalColor(const TDF_Label& theItemL,
                                                   Quantity_Color&  theColor);

  //! Puts the item theItemL on the layer theName (NIST's are on layers, '-1' 'Layer-1', which a
  //! file may hide: FTC-07 hides 13 items so). Its being hidden is SetPresentationHidden's.
  Standard_EXPORT static void SetSupplementalLayer(const TDF_Label&               theItemL,
                                                   const TCollection_AsciiString& theName,
                                                   const TCollection_AsciiString& theDescription);

  //! The layer of the item theItemL, its name and description; false when it is on none.
  Standard_EXPORT static bool GetSupplementalLayer(const TDF_Label&         theItemL,
                                                   TCollection_AsciiString& theName,
                                                   TCollection_AsciiString& theDescription);

  //! Says the presentation of the GD&T theGDTL is shown or not (STEP: its graphic listed by an
  //! INVISIBILITY -- NIST STC-08 hides 29, STC-09 75).
  Standard_EXPORT static void SetPresentationHidden(const TDF_Label& theGDTL, const bool theHidden);

  //! Whether the presentation of the GD&T theGDTL is not shown.
  Standard_EXPORT static bool IsPresentationHidden(const TDF_Label& theGDTL);

  //! Says the hidden GD&T theGDTL had its callout listed by the INVISIBILITY, not its graphic
  //! (NIST FTC-09): written so again.
  Standard_EXPORT static void SetPresentationHiddenCallout(const TDF_Label& theGDTL,
                                                           const bool       theCallout);

  //! Whether the hidden GD&T theGDTL is hidden by its callout.
  Standard_EXPORT static bool IsPresentationHiddenCallout(const TDF_Label& theGDTL);

  //! Says the annotation placeholder of the GD&T theGDTL is not shown, while it is hidden (STEP:
  //! listed by an INVISIBILITY; NIST STC-06 hides 19, STC-08 none of its hidden PMI's).
  Standard_EXPORT static void SetPlaceholderHidden(const TDF_Label& theGDTL, const bool theHidden);

  //! Whether the annotation placeholder of the GD&T theGDTL is not shown.
  Standard_EXPORT static bool IsPlaceholderHidden(const TDF_Label& theGDTL);

  //! Says the presentation of the GD&T theGDTL is given as curves (STEP: an annotation curve
  //! occurrence of polylines and circles, NIST CTC-05's), not tessellated.
  Standard_EXPORT static void SetPresentationCurves(const TDF_Label& theGDTL, const bool theCurves);

  //! Whether the presentation of the GD&T theGDTL is given as curves.
  Standard_EXPORT static bool IsPresentationCurves(const TDF_Label& theGDTL);

  //! A hole as a STEP file defines it (AP242 ed3 hole features: NIST STC-09's
  //! basic_round_hole): its kind (the STEP type), name, diameter, depth (none when
  //! HasDepth is false), each with a tolerance or none, whether it goes through, and its
  //! axes in its own coordinates.
  struct HoleDefinition
  {
    TCollection_AsciiString Kind;
    TCollection_AsciiString Name;
    double                  Diameter            = 0.0;
    bool                    DiameterToleranced  = false;
    double                  DiameterLower       = 0.0;
    double                  DiameterUpper       = 0.0;
    bool                    HasDepth            = false;
    double                  Depth               = 0.0;
    bool                    DepthToleranced     = false;
    double                  DepthLower          = 0.0;
    double                  DepthUpper          = 0.0;
    bool                    Through             = false;
    gp_Ax3                  Axes;
  };

  //! Adds a hole of the part theShapeL; its occurrences are added with AddHoleOccurrence.
  Standard_EXPORT TDF_Label AddHole(const TDF_Label& theShapeL, const HoleDefinition& theHole) const;

  //! The holes, in the order they were added.
  Standard_EXPORT void GetHoles(NCollection_Sequence<TDF_Label>& theHoles) const;

  //! The part and definition of the hole theHoleL; false when it is none.
  Standard_EXPORT static bool GetHole(const TDF_Label& theHoleL,
                                      TDF_Label&       theShapeL,
                                      HoleDefinition&  theHole);

  //! Adds an occurrence of the hole theHoleL: its name, where it is in its part's
  //! coordinates, and the faces it is of.
  Standard_EXPORT static TDF_Label AddHoleOccurrence(const TDF_Label&               theHoleL,
                                                     const TCollection_AsciiString& theName,
                                                     const gp_Ax3&                  thePosition,
                                                     const NCollection_Sequence<TDF_Label>& theFaces);

  //! The occurrences of the hole theHoleL.
  Standard_EXPORT static void GetHoleOccurrences(const TDF_Label&                 theHoleL,
                                                 NCollection_Sequence<TDF_Label>& theOccurrences);

  //! The name, position and faces of the occurrence theOccurrenceL; false when it is none.
  Standard_EXPORT static bool GetHoleOccurrence(const TDF_Label&                 theOccurrenceL,
                                                TCollection_AsciiString&         theName,
                                                gp_Ax3&                          thePosition,
                                                NCollection_Sequence<TDF_Label>& theFaces);

  //! Gives the label theL -- a product's shape, a GD&T, a hole or an occurrence of one -- the
  //! UUID a STEP file identifies it by (AP242 ed4 uuid_attribute: NIST STC-09's
  //! v5_uuid_attribute), of the type theType ('V5_UUID_ATTRIBUTE'). An empty theUUID removes it.
  Standard_EXPORT static void SetUUID(const TDF_Label&               theL,
                                      const TCollection_AsciiString& theType,
                                      const TCollection_AsciiString& theUUID);

  //! The type and UUID of the label theL; false when it has none.
  Standard_EXPORT static bool GetUUID(const TDF_Label&         theL,
                                      TCollection_AsciiString& theType,
                                      TCollection_AsciiString& theUUID);

  //! Returns Datum label defined for theShapeL label.
  Standard_EXPORT bool GetRefDatumLabel(const TDF_Label&                 theShapeL,
                                        NCollection_Sequence<TDF_Label>& theDatum) const;

  //! Returns true if the given GDT is marked as locked.
  Standard_EXPORT bool IsLocked(const TDF_Label& theViewL) const;

  //! Mark the given GDT as locked.
  Standard_EXPORT void Lock(const TDF_Label& theViewL) const;

  //! fill the map GDT label -> shape presentation
  Standard_EXPORT void GetGDTPresentations(
    NCollection_IndexedDataMap<TDF_Label, TopoDS_Shape>& theGDTLabelToShape) const;

  //! Set shape presentation for GDT labels according to given map (theGDTLabelToPrs)
  //! theGDTLabelToPrsName map is an additional argument, can be used to set presentation names.
  //! If label is not in the theGDTLabelToPrsName map, the presentation name will be empty
  Standard_EXPORT void SetGDTPresentations(
    NCollection_IndexedDataMap<TDF_Label, TopoDS_Shape>& theGDTLabelToPrs);

  //! Unlock the given GDT.
  Standard_EXPORT void Unlock(const TDF_Label& theViewL) const;

  Standard_EXPORT const Standard_GUID& ID() const override;

  //! Dumps the content of me into the stream
  Standard_EXPORT void DumpJson(Standard_OStream& theOStream, int theDepth = -1) const override;

  DEFINE_DERIVED_ATTRIBUTE(XCAFDoc_DimTolTool, TDataStd_GenericEmpty)

private:
  occ::handle<XCAFDoc_ShapeTool> myShapeTool;
};

#endif // _XCAFDoc_DimTolTool_HeaderFile
