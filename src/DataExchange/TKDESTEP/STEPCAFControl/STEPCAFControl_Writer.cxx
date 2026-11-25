// Created on: 2000-08-15
// Created by: Andrey BETENEV
// Copyright (c) 2000-2014 OPEN CASCADE SAS
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

// CURRENT LIMITATIONS:
// when val props and names assigned to instance of
// component in assembly, it is in fact supposed that only one CDSR corresponds
// to such shape. This can be wrong and should be handled more carefully
// (analysis of SDRs which the CDSR links should be done)
// Names and validation props are supported for top-level shapes only

#include <STEPCAFControl_Writer.hxx>

#include <cstdio>

#include <BRep_Builder.hxx>
#include <GeomToStep_MakeAxis2Placement3d.hxx>
#include <GeomToStep_MakeCartesianPoint.hxx>
#include <GeomToStep_MakePlane.hxx>
#include <gp_Pln.hxx>
#include <HeaderSection_FileSchema.hxx>
#include <Interface_EntityIterator.hxx>
#include <Interface_Static.hxx>
#include <Message.hxx>
#include <Message_ProgressScope.hxx>
#include <TopoDS_Shape.hxx>
#include <Standard_Transient.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <NCollection_DataMap.hxx>
#include <OSD_Path.hxx>
#include <StepAP242_DraughtingModelItemAssociation.hxx>
#include <StepAP242_GeometricItemSpecificUsage.hxx>
#include <StepBasic_ConversionBasedUnitAndLengthUnit.hxx>
#include <StepBasic_ConversionBasedUnitAndPlaneAngleUnit.hxx>
#include <StepBasic_DerivedUnit.hxx>
#include <StepBasic_DerivedUnitElement.hxx>
#include <StepBasic_DimensionalExponents.hxx>
#include <StepBasic_GeneralProperty.hxx>
#include <StepBasic_GeneralPropertyAssociation.hxx>
#include <StepBasic_LengthMeasureWithUnit.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_HArray1.hxx>
#include <StepBasic_MeasureValueMember.hxx>
#include <StepBasic_PlaneAngleMeasureWithUnit.hxx>
#include <StepBasic_Product.hxx>
#include <StepBasic_ProductDefinition.hxx>
#include <StepBasic_ProductDefinitionFormation.hxx>
#include <StepBasic_SiUnitAndLengthUnit.hxx>
#include <StepBasic_SiUnitAndMassUnit.hxx>
#include <StepBasic_SiUnitAndPlaneAngleUnit.hxx>
#include <STEPCAFControl_ActorWrite.hxx>
#include <STEPCAFControl_Controller.hxx>
#include <STEPCAFControl_ExternFile.hxx>
#include <STEPConstruct.hxx>
#include <STEPConstruct_ExternRefs.hxx>
#include <STEPCAFControl_GDTProperty.hxx>
#include <STEPConstruct_Styles.hxx>
#include <STEPConstruct_ValidationProps.hxx>
#include <STEPControl_StepModelType.hxx>
#include <StepData_Factors.hxx>
#include <StepData_Logical.hxx>
#include <StepData_StepModel.hxx>
#include <StepDimTol_ProjectedZoneDefinition.hxx>
#include <StepDimTol_AngularityTolerance.hxx>
#include <StepDimTol_CircularRunoutTolerance.hxx>
#include <StepDimTol_CoaxialityTolerance.hxx>
#include <StepDimTol_ConcentricityTolerance.hxx>
#include <StepDimTol_Datum.hxx>
#include <StepDimTol_DatumFeature.hxx>
#include <StepDimTol_DatumReference.hxx>
#include <StepDimTol_DatumSystem.hxx>
#include <StepDimTol_DatumSystemOrReference.hxx>
#include <StepDimTol_GeometricToleranceRelationship.hxx>
#include <StepDimTol_GeometricToleranceType.hxx>
#include <StepDimTol_GeometricToleranceWithDatumReference.hxx>
#include <StepDimTol_GeometricToleranceWithModifiers.hxx>
#include <StepDimTol_GeoTolAndGeoTolWthDatRef.hxx>
#include <StepDimTol_GeoTolAndGeoTolWthDatRefAndGeoTolWthMaxTol.hxx>
#include <StepDimTol_GeoTolAndGeoTolWthDatRefAndGeoTolWthMod.hxx>
#include <StepDimTol_GeoTolAndGeoTolWthDatRefAndModGeoTolAndPosTol.hxx>
#include <StepDimTol_GeoTolAndGeoTolWthMaxTol.hxx>
#include <StepDimTol_GeoTolAndGeoTolWthDefUnit.hxx>
#include <StepDimTol_GeoTolAndGeoTolWthMod.hxx>
#include <StepDimTol_DatumReferenceElement.hxx>
#include <StepDimTol_DatumReferenceModifier.hxx>
#include <StepDimTol_ModifiedGeometricTolerance.hxx>
#include <StepDimTol_ParallelismTolerance.hxx>
#include <StepDimTol_PerpendicularityTolerance.hxx>
#include <StepDimTol_PlacedDatumTargetFeature.hxx>
#include <StepDimTol_RunoutZoneDefinition.hxx>
#include <StepDimTol_RunoutZoneOrientation.hxx>
#include <StepDimTol_SymmetryTolerance.hxx>
#include <StepDimTol_ToleranceZone.hxx>
#include <StepDimTol_ToleranceZoneForm.hxx>
#include <StepDimTol_TotalRunoutTolerance.hxx>
#include <StepGeom_Axis2Placement3d.hxx>
#include <StepGeom_Axis2Placement.hxx>
#include <StepGeom_CartesianPoint.hxx>
#include <StepGeom_Direction.hxx>
#include <StepGeom_GeometricRepresentationContextAndGlobalUnitAssignedContext.hxx>
#include <StepGeom_GeomRepContextAndGlobUnitAssCtxAndGlobUncertaintyAssCtx.hxx>
#include <StepGeom_Plane.hxx>
#include <StepRepr_CompGroupShAspAndCompShAspAndDatumFeatAndShAsp.hxx>
#include <StepRepr_AllAroundShapeAspect.hxx>
#include <StepRepr_CompositeShapeAspect.hxx>
#include <StepRepr_ConstructiveGeometryRepresentation.hxx>
#include <StepRepr_ConstructiveGeometryRepresentationRelationship.hxx>
#include <GeomToStep_MakeCylindricalSurface.hxx>
#include <StepGeom_CylindricalSurface.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <TopoDSToStep_MakeGeometricCurveSet.hxx>
#include <StepShape_GeometricCurveSet.hxx>
#include <StepVisual_SurfaceStyleBoundary.hxx>
#include <StepVisual_SurfaceSideStyle.hxx>
#include <StepVisual_SurfaceStyleUsage.hxx>
#include <StepVisual_CurveStyle.hxx>
#include <StepVisual_PresentationStyleAssignment.hxx>
#include <StepVisual_StyledItem.hxx>
#include <StepShape_GeometricSet.hxx>
#include <StepGeom_TrimmedCurve.hxx>
#include <STEPConstruct_RenderingProperties.hxx>
#include <Quantity_Color.hxx>
#include <StepVisual_Colour.hxx>
#include <Transfer_SimpleBinderOfTransient.hxx>
#include <Standard_ErrorHandler.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <StepRepr_DerivedShapeAspect.hxx>
#include <StepRepr_DescriptiveRepresentationItem.hxx>
#include <StepRepr_FeatureForDatumTargetRelationship.hxx>
#include <StepRepr_RepresentationItem.hxx>
#include <StepRepr_IntegerRepresentationItem.hxx>
#include <StepRepr_MeasureRepresentationItem.hxx>
#include <StepRepr_NextAssemblyUsageOccurrence.hxx>
#include <StepRepr_ProductDefinitionShape.hxx>
#include <StepRepr_CharacterizedDefinition.hxx>
#include <StepRepr_PropertyDefinition.hxx>
#include <StepRepr_PropertyDefinitionRepresentation.hxx>
#include <StepRepr_RepresentedDefinition.hxx>
#include <StepRepr_RealRepresentationItem.hxx>
#include <StepRepr_Representation.hxx>
#include <StepRepr_ReprItemAndLengthMeasureWithUnit.hxx>
#include <StepRepr_ReprItemAndLengthMeasureWithUnitAndQRI.hxx>
#include <StepRepr_ReprItemAndPlaneAngleMeasureWithUnit.hxx>
#include <StepRepr_ReprItemAndPlaneAngleMeasureWithUnitAndQRI.hxx>
#include <StepRepr_ShapeAspect.hxx>
#include <StepRepr_ShapeAspectDerivingRelationship.hxx>
#include <StepRepr_ShapeAspectRelationship.hxx>
#include <StepRepr_SpecifiedHigherUsageOccurrence.hxx>
#include <StepRepr_ValueRange.hxx>
#include <StepShape_AdvancedFace.hxx>
#include <StepShape_AngularLocation.hxx>
#include <StepShape_AngularSize.hxx>
#include <StepShape_ConnectedFaceSet.hxx>
#include <StepShape_ContextDependentShapeRepresentation.hxx>
#include <StepShape_DimensionalCharacteristicRepresentation.hxx>
#include <StepShape_DimensionalLocation.hxx>
#include <StepShape_DimensionalLocationWithPath.hxx>
#include <StepShape_DimensionalSize.hxx>
#include <StepShape_DimensionalSizeWithPath.hxx>
#include <StepShape_EdgeCurve.hxx>
#include <StepShape_EdgeLoop.hxx>
#include <StepShape_FaceBound.hxx>
#include <StepShape_LimitsAndFits.hxx>
#include <StepShape_OrientedEdge.hxx>
#include <StepShape_PlusMinusTolerance.hxx>
#include <StepShape_QualifiedRepresentationItem.hxx>
#include <StepShape_ShapeDefinitionRepresentation.hxx>
#include <StepShape_ShapeDimensionRepresentation.hxx>
#include <StepShape_ShapeRepresentation.hxx>
#include <StepShape_ShapeRepresentationWithParameters.hxx>
#include <StepShape_ToleranceValue.hxx>
#include <StepShape_TypeQualifier.hxx>
#include <StepShape_ValueFormatTypeQualifier.hxx>
#include <StepVisual_AnnotationPlane.hxx>
#include <StepVisual_CurveStyle.hxx>
#include <StepVisual_DraughtingCallout.hxx>
#include <StepVisual_DraughtingModel.hxx>
#include <StepVisual_PresentationStyleAssignment.hxx>
#include <StepVisual_PresentationStyleSelect.hxx>
#include <StepVisual_Invisibility.hxx>
#include <StepVisual_InvisibleItem.hxx>
#include <StepVisual_MechanicalDesignGeometricPresentationRepresentation.hxx>
#include <StepVisual_NullStyleMember.hxx>
#include <StepVisual_PointStyle.hxx>
#include <StepVisual_PresentationLayerAssignment.hxx>
#include <StepVisual_PresentationRepresentation.hxx>
#include <StepVisual_PresentationStyleByContext.hxx>
#include <StepVisual_StyledItem.hxx>
#include <StepVisual_SurfaceStyleUsage.hxx>
#include <StepVisual_TessellatedAnnotationOccurrence.hxx>
#include <StepVisual_TessellatedGeometricSet.hxx>
#include <TCollection_AsciiString.hxx>
#include <TCollection_HAsciiString.hxx>
#include <NCollection_Map.hxx>
#include <TDataStd_Name.hxx>
#include <TDataStd_TreeNode.hxx>
#include <TDataStd_UAttribute.hxx>
#include <TDF_ChildIterator.hxx>
#include <Interface_UndefinedContent.hxx>
#include <StepAP242_DraughtingModelItemAssociationWithPlaceholder.hxx>
#include <StepData_UndefinedEntity.hxx>
#include <StepShape_GeometricSet.hxx>
#include <StepShape_GeometricSetSelect.hxx>
#include <TDF_Label.hxx>
#include <TDF_Tool.hxx>
#include <NCollection_Sequence.hxx>
#include <TDocStd_Document.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Iterator.hxx>
#include <Transfer_Binder.hxx>
#include <Transfer_FinderProcess.hxx>
#include <Transfer_TransientListBinder.hxx>
#include <StepRepr_MappedItem.hxx>
#include <StepRepr_MechanicalDesignAndDraughtingRelationship.hxx>
#include <StepRepr_RepresentationMap.hxx>
#include <StepVisual_CameraModelD3.hxx>
#include <StepVisual_CameraModelD3MultiClipping.hxx>
#include <StepVisual_CameraModelD3MultiClippingInterectionSelect.hxx>
#include <StepVisual_PlanarBox.hxx>
#include <StepVisual_ViewVolume.hxx>
#include <Precision.hxx>
#include <TransferBRep.hxx>
#include <TransferBRep_ShapeMapper.hxx>
#include <XCAFDimTolObjects_DatumObject.hxx>
#include <XCAFDimTolObjects_DimensionFormVariance.hxx>
#include <XCAFDimTolObjects_DimensionGrade.hxx>
#include <XCAFDimTolObjects_DimensionObject.hxx>
#include <XCAFDimTolObjects_DimensionModif.hxx>
#include <XCAFDimTolObjects_DimensionQualifier.hxx>
#include <XCAFDimTolObjects_GeomToleranceObject.hxx>
#include <XCAFDoc.hxx>
#include <XCAFDoc_Area.hxx>
#include <XCAFDoc_Centroid.hxx>
#include <XCAFDoc_ColorTool.hxx>
#include <XCAFDoc_Datum.hxx>
#include <XCAFDoc_Dimension.hxx>
#include <XCAFDoc_DimTol.hxx>
#include <XCAFDoc_DimTolTool.hxx>
#include <XCAFDoc_ClippingPlaneTool.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_GeomTolerance.hxx>
#include <XCAFDoc_GraphNode.hxx>
#include <XCAFDoc_LengthUnit.hxx>
#include <XCAFDoc_LayerTool.hxx>
#include <XCAFDoc_MaterialTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <XCAFDoc_Volume.hxx>
#include <XCAFDoc_VisMaterialTool.hxx>
#include <XCAFDoc_View.hxx>
#include <XCAFDoc_ViewTool.hxx>
#include <XCAFView_Object.hxx>
#include <XCAFPrs_Style.hxx>
#include <NCollection_IndexedDataMap.hxx>
#include <XSAlgo.hxx>
#include <XSAlgo_ShapeProcessor.hxx>
#include <XSControl_TransferWriter.hxx>
#include <XSControl_WorkSession.hxx>
#include <UnitsMethods.hxx>

//=======================================================================
// function : GetLabelName
// purpose  : auxiliary function: take name of label and append it to str
//=======================================================================
static bool GetLabelName(const TDF_Label& theL, occ::handle<TCollection_HAsciiString>& theStr)
{
  occ::handle<TDataStd_Name> anAttrName;
  if (!theL.FindAttribute(TDataStd_Name::GetID(), anAttrName))
  {
    return false;
  }
  TCollection_ExtendedString aName = anAttrName->Get();
  if (aName.IsEmpty())
  {
    return false;
  }

  // set name, removing spaces around it
  TCollection_AsciiString aBuffer(aName);
  aBuffer.LeftAdjust();
  aBuffer.RightAdjust();
  theStr->AssignCat(aBuffer.ToCString());
  return true;
}

//=================================================================================================

//! The draughting model of the GD&T presentations, with the writer's state
//! that has no member in it: STEPCAFControl_Writer is constructed on its
//! user's stack, its layout is its ABI. Of no type of its own -- it is
//! written as the draughting model it is.
namespace
{
class GDTPresentationModel : public StepVisual_DraughtingModel
{
public:
  //! The colour GD&T presentations are drawn in
  Quantity_Color Color = Quantity_Color(Quantity_NOC_WHITE);
  //! The presentation of each GD&T label written: its draughting callout,
  //! and its annotation plane when it has one -- what a view of it lists
  NCollection_DataMap<TDF_Label, NCollection_Sequence<occ::handle<StepRepr_RepresentationItem>>>
    Items;
  //! A PMI with no graphic of its own, presented by the one another PMI holds
  //! (XCAFDoc_DimTolTool::SetPresentedWith; NIST STC-07's 32 radii by one '32X R', a
  //! hole's depth by its callout): linked to it once all are written, in whatever
  //! order they were
  struct Link
  {
    TDF_Label                                               Label;
    TDF_Label                                               Host;
    occ::handle<Standard_Transient>                         Dimension;
    NCollection_Sequence<occ::handle<StepRepr_ShapeAspect>> Geometry;
  };
  NCollection_Sequence<Link> Links;
  //! The placeholders written, by kind, name and where they are: two PMI with one (NIST
  //! STC-06's hole callout, of its diameter and its depth) are linked to it both, not given
  //! one each. Its callout first, then the placeholder
  NCollection_DataMap<TCollection_AsciiString, NCollection_Array1<occ::handle<Standard_Transient>>>
    Placeholders;
  //! Supplemental geometry written (XCAFDoc_DimTolTool::AddSupplementalGeometry): each set of
  //! curves, the item GD&T on it identifies, with its constructive geometry representation --
  //! the representation the usage is of -- and its part's shape
  struct Supplemental
  {
    occ::handle<StepRepr_Representation>         Representation;
    occ::handle<StepRepr_ProductDefinitionShape> PDS;
  };
  NCollection_DataMap<occ::handle<Standard_Transient>, Supplemental> Supplementals;
};

//! theModel is the writer's myGDTPresentationDM, made in Init() and nowhere else
GDTPresentationModel& gdtState(const occ::handle<StepVisual_DraughtingModel>& theModel)
{
  return *static_cast<GDTPresentationModel*>(theModel.get());
}
} // namespace

//=================================================================================================

void STEPCAFControl_Writer::SetGDTPresentationColor(const Quantity_Color& theColor)
{
  gdtState(myGDTPresentationDM).Color = theColor;
}

//=================================================================================================

const Quantity_Color& STEPCAFControl_Writer::GetGDTPresentationColor() const
{
  return gdtState(myGDTPresentationDM).Color;
}

//=================================================================================================

STEPCAFControl_Writer::STEPCAFControl_Writer()
    : myColorMode(true),
      myNameMode(true),
      myLayerMode(true),
      myPropsMode(true),
      myMetadataMode(true),
      mySHUOMode(true),
      myGDTMode(true),
      myMatMode(true),
      myVisMatMode(false),
      myIsCleanDuplicates(false)
{
  STEPCAFControl_Controller::Init();
  occ::handle<XSControl_WorkSession> aWS = new XSControl_WorkSession;
  Init(aWS);
}

//=================================================================================================

STEPCAFControl_Writer::STEPCAFControl_Writer(const occ::handle<XSControl_WorkSession>& theWS,
                                             const bool                                theScratch)
    : myColorMode(true),
      myNameMode(true),
      myLayerMode(true),
      myPropsMode(true),
      myMetadataMode(true),
      mySHUOMode(true),
      myGDTMode(true),
      myMatMode(true),
      myVisMatMode(false),
      myIsCleanDuplicates(false)
{
  STEPCAFControl_Controller::Init();
  Init(theWS, theScratch);
}

//=================================================================================================

void STEPCAFControl_Writer::Init(const occ::handle<XSControl_WorkSession>& theWS,
                                 const bool                                theScratch)
{
  theWS->SelectNorm("STEP");
  myWriter.SetWS(theWS, theScratch);
  myFiles.Clear();
  myLabEF.Clear();
  myLabels.Clear();
  myPureRefLabels.Clear();
  myRootLabels.Clear();
  // The colour set is kept through an Init, as a member was
  Quantity_Color aColor(Quantity_NOC_WHITE);
  if (!myGDTPresentationDM.IsNull())
  {
    aColor = gdtState(myGDTPresentationDM).Color;
  }
  occ::handle<GDTPresentationModel> aDM = new GDTPresentationModel();
  aDM->Color                            = aColor;
  myGDTPresentationDM                   = aDM;
  myGDTPrsCurveStyle =
    new NCollection_HArray1<occ::handle<StepVisual_PresentationStyleAssignment>>(1, 1);
}

//=================================================================================================

IFSelect_ReturnStatus STEPCAFControl_Writer::Write(const char* const theFileName)
{
  if (myIsCleanDuplicates)
  {
    // remove duplicates
    myWriter.CleanDuplicateEntities();
  }

  IFSelect_ReturnStatus aStatus = myWriter.Write(theFileName);
  if (aStatus != IFSelect_RetDone)
  {
    return aStatus;
  }

  // get directory name of the main file
  TCollection_AsciiString aDirPath;
  {
    OSD_Path aMainFile(theFileName);
    aMainFile.SetName("");
    aMainFile.SetExtension("");
    aMainFile.SystemName(aDirPath);
  }

  for (NCollection_DataMap<TCollection_AsciiString,
                           occ::handle<STEPCAFControl_ExternFile>>::Iterator anExtFileIter(myFiles);
       anExtFileIter.More();
       anExtFileIter.Next())
  {
    const occ::handle<STEPCAFControl_ExternFile>& anExtFile = anExtFileIter.Value();
    if (anExtFile->GetWriteStatus() != IFSelect_RetVoid)
    {
      continue;
    }

    // construct extern file name
    TCollection_AsciiString aFileName =
      OSD_Path::AbsolutePath(aDirPath, anExtFile->GetName()->String());
    if (aFileName.Length() <= 0)
    {
      aFileName = anExtFile->GetName()->String();
    }
    Message::SendTrace() << "Writing external file: " << aFileName << "\n";

    const IFSelect_ReturnStatus anExtStatus = anExtFile->GetWS()->SendAll(aFileName.ToCString());
    anExtFile->SetWriteStatus(anExtStatus);
    if (anExtStatus != IFSelect_RetDone)
    {
      aStatus = anExtStatus;
    }
  }

  return aStatus;
}

//=================================================================================================

void STEPCAFControl_Writer::prepareUnit(const TDF_Label&                       theLabel,
                                        const occ::handle<StepData_StepModel>& theModel,
                                        StepData_Factors&                      theLocalFactors)
{
  occ::handle<XCAFDoc_LengthUnit> aLengthAttr;
  if (!theLabel.IsNull() && theLabel.Root().FindAttribute(XCAFDoc_LengthUnit::GetID(), aLengthAttr))
  {
    theModel->SetLocalLengthUnit(aLengthAttr->GetUnitValue() * 1000); // convert to mm
    theLocalFactors.SetCascadeUnit(aLengthAttr->GetUnitValue() * 1000);
  }
  else
  {
    XSAlgo_ShapeProcessor::PrepareForTransfer(); // update unit info
    theModel->SetLocalLengthUnit(UnitsMethods::GetCasCadeLengthUnit());
    theLocalFactors.SetCascadeUnit(UnitsMethods::GetCasCadeLengthUnit());
  }
}

//=================================================================================================

IFSelect_ReturnStatus STEPCAFControl_Writer::WriteStream(std::ostream& theStream)
{
  if (!myFiles.IsEmpty())
  {
    // writing external files is unsupported via stream interface
    return IFSelect_RetError;
  }

  return myWriter.WriteStream(theStream);
}

//=================================================================================================

bool STEPCAFControl_Writer::Transfer(const occ::handle<TDocStd_Document>& theDoc,
                                     const STEPControl_StepModelType      theMode,
                                     const char* const                    theMulti,
                                     const Message_ProgressRange&         theProgress)
{
  const occ::handle<StepData_StepModel> aModel =
    occ::down_cast<StepData_StepModel>(myWriter.WS()->Model());
  aModel->InternalParameters.InitFromStatic();
  return Transfer(theDoc, aModel->InternalParameters, theMode, theMulti, theProgress);
}

//=================================================================================================

bool STEPCAFControl_Writer::Transfer(const occ::handle<TDocStd_Document>& theDoc,
                                     const DESTEP_Parameters&             theParams,
                                     const STEPControl_StepModelType      theMode,
                                     const char* const                    theMulti,
                                     const Message_ProgressRange&         theProgress)
{
  occ::handle<XCAFDoc_ShapeTool> aShTool = XCAFDoc_DocumentTool::ShapeTool(theDoc->Main());
  if (aShTool.IsNull())
  {
    return false;
  }

  myRootLabels.Clear();
  myRootLabels.Add(theDoc->Main().Root());
  NCollection_Sequence<TDF_Label> aLabels;
  aShTool->GetFreeShapes(aLabels);
  const occ::handle<StepData_StepModel> aModel =
    occ::down_cast<StepData_StepModel>(myWriter.WS()->Model());
  aModel->InternalParameters = theParams;
  return transfer(myWriter, aLabels, theMode, theMulti, false, theProgress);
}

//=================================================================================================

bool STEPCAFControl_Writer::Transfer(const TDF_Label&                theLabel,
                                     const STEPControl_StepModelType theMode,
                                     const char* const               theIsMulti,
                                     const Message_ProgressRange&    theProgress)
{
  const occ::handle<StepData_StepModel> aModel =
    occ::down_cast<StepData_StepModel>(myWriter.WS()->Model());
  aModel->InternalParameters.InitFromStatic();
  return Transfer(theLabel, aModel->InternalParameters, theMode, theIsMulti, theProgress);
}

//=================================================================================================

bool STEPCAFControl_Writer::Transfer(const TDF_Label&                theLabel,
                                     const DESTEP_Parameters&        theParams,
                                     const STEPControl_StepModelType theMode,
                                     const char* const               theIsMulti,
                                     const Message_ProgressRange&    theProgress)
{
  if (theLabel.IsNull())
  {
    return false;
  }
  NCollection_Sequence<TDF_Label> aLabels;
  aLabels.Append(theLabel);
  myRootLabels.Clear();
  myRootLabels.Add(theLabel.Root());
  const occ::handle<StepData_StepModel> aModel =
    occ::down_cast<StepData_StepModel>(myWriter.WS()->Model());
  aModel->InternalParameters = theParams;
  return transfer(myWriter, aLabels, theMode, theIsMulti, false, theProgress);
}

//=================================================================================================

bool STEPCAFControl_Writer::Transfer(const NCollection_Sequence<TDF_Label>& theLabels,
                                     const STEPControl_StepModelType        theMode,
                                     const char* const                      theIsMulti,
                                     const Message_ProgressRange&           theProgress)
{
  const occ::handle<StepData_StepModel> aModel =
    occ::down_cast<StepData_StepModel>(myWriter.WS()->Model());
  aModel->InternalParameters.InitFromStatic();
  return Transfer(theLabels, aModel->InternalParameters, theMode, theIsMulti, theProgress);
}

//=================================================================================================

bool STEPCAFControl_Writer::Transfer(const NCollection_Sequence<TDF_Label>& theLabels,
                                     const DESTEP_Parameters&               theParams,
                                     const STEPControl_StepModelType        theMode,
                                     const char* const                      theIsMulti,
                                     const Message_ProgressRange&           theProgress)
{
  myRootLabels.Clear();
  for (NCollection_Sequence<TDF_Label>::Iterator aLabelIter(theLabels); aLabelIter.More();
       aLabelIter.Next())
  {
    const TDF_Label& aLabel = aLabelIter.Value();
    if (!aLabel.IsNull())
    {
      myRootLabels.Add(aLabel.Root());
    }
  }
  const occ::handle<StepData_StepModel> aModel =
    occ::down_cast<StepData_StepModel>(myWriter.WS()->Model());
  aModel->InternalParameters = theParams;
  return transfer(myWriter, theLabels, theMode, theIsMulti, false, theProgress);
}

//=================================================================================================

bool STEPCAFControl_Writer::Perform(const occ::handle<TDocStd_Document>& theDoc,
                                    const char* const                    theFileName,
                                    const Message_ProgressRange&         theProgress)
{
  if (!Transfer(theDoc, STEPControl_AsIs, nullptr, theProgress))
  {
    return false;
  }
  return Write(theFileName) == IFSelect_RetDone;
}

//=================================================================================================

bool STEPCAFControl_Writer::Perform(const occ::handle<TDocStd_Document>& theDoc,
                                    const char* const                    theFileName,
                                    const DESTEP_Parameters&             theParams,
                                    const Message_ProgressRange&         theProgress)
{
  if (!Transfer(theDoc, theParams, STEPControl_AsIs, nullptr, theProgress))
  {
    return false;
  }
  return Write(theFileName) == IFSelect_RetDone;
}

//=================================================================================================

bool STEPCAFControl_Writer::Perform(const occ::handle<TDocStd_Document>& theDoc,
                                    const TCollection_AsciiString&       theFileName,
                                    const Message_ProgressRange&         theProgress)
{
  if (!Transfer(theDoc, STEPControl_AsIs, nullptr, theProgress))
  {
    return false;
  }
  return Write(theFileName.ToCString()) == IFSelect_RetDone;
}

//=================================================================================================

bool STEPCAFControl_Writer::ExternFile(const TDF_Label&                        theLabel,
                                       occ::handle<STEPCAFControl_ExternFile>& theExtFile) const
{
  theExtFile.Nullify();
  if (!myLabEF.IsBound(theLabel))
  {
    return false;
  }
  theExtFile = myLabEF.Find(theLabel);
  return true;
}

//=================================================================================================

bool STEPCAFControl_Writer::ExternFile(const char* const                       theName,
                                       occ::handle<STEPCAFControl_ExternFile>& theExtFile) const
{
  theExtFile.Nullify();
  if (!myFiles.IsEmpty() || !myFiles.IsBound(theName))
  {
    return false;
  }
  theExtFile = myFiles.Find(theName);
  return true;
}

//=============================================================================

void STEPCAFControl_Writer::SetShapeFixParameters(
  const XSAlgo_ShapeProcessor::ParameterMap& theParameters)
{
  myWriter.SetShapeFixParameters(theParameters);
}

//=============================================================================

void STEPCAFControl_Writer::SetShapeFixParameters(
  XSAlgo_ShapeProcessor::ParameterMap&& theParameters)
{
  myWriter.SetShapeFixParameters(std::move(theParameters));
}

//=============================================================================

void STEPCAFControl_Writer::SetShapeFixParameters(
  const DE_ShapeFixParameters&               theParameters,
  const XSAlgo_ShapeProcessor::ParameterMap& theAdditionalParameters)
{
  myWriter.SetShapeFixParameters(theParameters, theAdditionalParameters);
}

//=============================================================================

const XSAlgo_ShapeProcessor::ParameterMap& STEPCAFControl_Writer::GetShapeFixParameters() const
{
  return myWriter.GetShapeFixParameters();
}

//=============================================================================

void STEPCAFControl_Writer::SetShapeProcessFlags(const ShapeProcess::OperationsFlags& theFlags)
{
  myWriter.SetShapeProcessFlags(theFlags);
}

//=============================================================================

const XSAlgo_ShapeProcessor::ProcessingFlags& STEPCAFControl_Writer::GetShapeProcessFlags() const
{
  return myWriter.GetShapeProcessFlags();
}

//=================================================================================================

bool STEPCAFControl_Writer::transfer(STEPControl_Writer&                    theWriter,
                                     const NCollection_Sequence<TDF_Label>& theLabels,
                                     const STEPControl_StepModelType        theMode,
                                     const char* const                      theIsMulti,
                                     const bool                             theIsExternFile,
                                     const Message_ProgressRange&           theProgress)
{
  if (theLabels.IsEmpty())
  {
    return false;
  }

  occ::handle<STEPCAFControl_ActorWrite> anActor =
    occ::down_cast<STEPCAFControl_ActorWrite>(theWriter.WS()->NormAdaptor()->ActorWrite());

  StepData_Factors                      aLocalFactors;
  const occ::handle<StepData_StepModel> aModel =
    occ::down_cast<StepData_StepModel>(theWriter.WS()->Model());
  prepareUnit(theLabels.First(), aModel, aLocalFactors); // set local length unit to the model
  // translate free top-level shapes of the DECAF document
  NCollection_Sequence<TDF_Label> aSubLabels;
  Message_ProgressScope           aPS(theProgress, "Labels", theLabels.Length());
  // Iterate on requested shapes
  for (NCollection_Sequence<TDF_Label>::Iterator aLabelIter(theLabels);
       aLabelIter.More() && aPS.More();
       aLabelIter.Next())
  {
    Message_ProgressRange aRange = aPS.Next();
    TDF_Label             aCurL  = aLabelIter.Value();
    if (myLabels.IsBound(aCurL))
    {
      continue; // already processed
    }

    TopoDS_Shape aCurShape = XCAFDoc_ShapeTool::GetShape(aCurL);
    if (aCurShape.IsNull())
    {
      continue;
    }
    TopoDS_Shape aShapeForBind = aCurShape;
    // write shape either as a whole, or as multifile (with extern refs)
    if (!theIsMulti)
    {
      anActor->SetStdMode(false);

      NCollection_Sequence<TDF_Label> aCompLabels;

      // For case when only part of assembly structure should be written in the document
      // if specified label is component of the assembly then
      // in order to save location of this component in the high-level assembly
      // and save name of high-level assembly it is necessary to represent structure of high-level
      // assembly as assembly with one component specified by current label. For that compound
      // containing only specified component is binded to the label of the high-level assembly. The
      // such way full structure of high-level assembly was replaced on the assembly containing one
      // component. For case when free shape reference is (located root) also create an auxiliary
      // assembly.
      if (XCAFDoc_ShapeTool::IsReference(aCurL))
      {
        TopoDS_Compound aComp;
        BRep_Builder    aB;
        aB.MakeCompound(aComp);
        aB.Add(aComp, aCurShape);
        aCurShape = aComp;
        myPureRefLabels.Bind(aCurL, aComp);
        aCompLabels.Append(aCurL);
        TDF_Label aRefL;
        if (XCAFDoc_ShapeTool::GetReferredShape(aCurL, aRefL))
        {
          if (XCAFDoc_ShapeTool::IsAssembly(aRefL))
          {
            XCAFDoc_ShapeTool::GetComponents(aRefL, aCompLabels, true);
          }
        }
        if (!XCAFDoc_ShapeTool::IsFree(aCurL))
        {
          aCurL = aCurL.Father();
        }
      }
      else
      {
        // fill sequence of (sub) shapes for which attributes should be written
        // and set actor to handle assemblies in a proper way
        if (XCAFDoc_ShapeTool::IsAssembly(aCurL))
        {
          XCAFDoc_ShapeTool::GetComponents(aCurL, aCompLabels, true);
        }
      }

      for (NCollection_Sequence<TDF_Label>::Iterator aCompIter(aCompLabels); aCompIter.More();
           aCompIter.Next())
      {
        const TDF_Label aCurSubShL = aCompIter.Value();
        if (myLabels.IsBound(aCurSubShL))
        {
          continue;
        }
        const TopoDS_Shape aCurSh = XCAFDoc_ShapeTool::GetShape(aCurSubShL);
        if (aCurSh.IsNull())
        {
          continue;
        }
        myLabels.Bind(aCurSubShL, aCurSh);
        aSubLabels.Append(aCurSubShL);
        TDF_Label aRefL;
        if (!XCAFDoc_ShapeTool::GetReferredShape(aCurSubShL, aRefL))
        {
          continue;
        }
        if (!myLabels.IsBound(aRefL))
        {
          TopoDS_Shape aRefSh = XCAFDoc_ShapeTool::GetShape(aRefL);
          myLabels.Bind(aRefL, aRefSh);
          aSubLabels.Append(aRefL);
          if (XCAFDoc_ShapeTool::IsAssembly(aRefL))
          {
            anActor->RegisterAssembly(aRefSh);
          }
        }
      }
      myLabels.Bind(aCurL, aShapeForBind);
      aSubLabels.Append(aCurL);

      if (XCAFDoc_ShapeTool::IsAssembly(aCurL) || XCAFDoc_ShapeTool::IsReference(aCurL))
      {
        anActor->RegisterAssembly(aCurShape);
      }

      theWriter.Transfer(aCurShape, theMode, aModel->InternalParameters, false, aRange);
      anActor->SetStdMode(true); // restore default behaviour
    }
    else
    {
      // translate final solids
      Message_ProgressScope aPS1(aRange, nullptr, 2);
      TopoDS_Shape          aSass =
        transferExternFiles(aCurL, theMode, aSubLabels, aLocalFactors, theIsMulti, aPS1.Next());
      if (aPS1.UserBreak())
      {
        return false;
      }

      DESTEP_Parameters::WriteMode_Assembly assemblymode = aModel->InternalParameters.WriteAssembly;
      aModel->InternalParameters.WriteAssembly           = DESTEP_Parameters::WriteMode_Assembly_On;
      theWriter.Transfer(aSass, STEPControl_AsIs, aModel->InternalParameters, true, aPS1.Next());
      aModel->InternalParameters.WriteAssembly = assemblymode;
    }
  }
  if (aPS.UserBreak())
  {
    return false;
  }

  theWriter.WS()->ComputeGraph(true); // Setting to clear and regenerate graph

  // write names
  if (GetNameMode())
  {
    writeNames(theWriter.WS(), aSubLabels);
  }

  if (!theIsMulti)
  {
    // write colors
    if (GetColorMode())
    {
      writeColors(theWriter.WS(), aSubLabels);
    }

    // write layers
    if (GetLayerMode())
    {
      writeLayers(theWriter.WS(), aSubLabels);
    }

    // write SHUO entities
    if (GetSHUOMode() && !theIsExternFile)
    {
      // do not store SHUO for extern reference for the moment
      writeSHUOs(theWriter.WS(), aSubLabels);
    }

    // write G&DTs
    if (GetDimTolMode())
    {
      if (aModel->InternalParameters.WriteSchema == 5)
      {
        // In the file's length unit, as the shapes are (STEPControl_ActorWrite): the
        // factors' length factor is otherwise 1, the GD&T written in mm in an inch file
        StepData_Factors aGDTFactors = aLocalFactors;
        aGDTFactors.InitializeFactors(aModel->WriteLengthUnit() / aModel->LocalLengthUnit(),
                                      aLocalFactors.PlaneAngleFactor(),
                                      aLocalFactors.SolidAngleFactor());
        writeDGTsAP242(theWriter.WS(), aSubLabels, aGDTFactors);
      }
      else
      {
        writeDGTs(theWriter.WS(), aSubLabels);
      }
    }

    // write Materials
    if (GetMaterialMode())
    {
      writeMaterials(theWriter.WS(), aSubLabels);
    }

    // register all MDGPRs in model
    for (NCollection_DataMap<TopoDS_Shape,
                             occ::handle<Standard_Transient>,
                             TopTools_ShapeMapHasher>::Iterator anItr(myMapCompMDGPR);
         anItr.More();
         anItr.Next())
    {
      aModel->AddWithRefs(anItr.Value());
    }
  }

  if (theIsMulti)
  { // external refs
    writeExternRefs(theWriter.WS(), aSubLabels);
  }

  // write validation props
  if (GetPropsMode())
  {
    writeValProps(theWriter.WS(), aSubLabels, theIsMulti);
  }

  // write metadata
  if (GetMetadataMode())
  {
    writeMetadata(theWriter.WS(), aSubLabels);
  }

  // refresh graph
  theWriter.WS()->ComputeGraph(true);

  // Write names for the sub-shapes
  if (aModel->InternalParameters.WriteSubshapeNames != 0)
  {
    const occ::handle<XSControl_TransferWriter>& TW  = this->ChangeWriter().WS()->TransferWriter();
    const occ::handle<Transfer_FinderProcess>&   aFP = TW->FinderProcess();

    // Iterate on requested sub shapes
    for (NCollection_Sequence<TDF_Label>::Iterator aLabelIter(aSubLabels); aLabelIter.More();
         aLabelIter.Next())
    {
      const TDF_Label& aCurL = aLabelIter.Value();

      for (TDF_ChildIterator aChildIter(aCurL, true); aChildIter.More(); aChildIter.Next())
      {
        const TDF_Label& aSubL = aChildIter.Value();

        // Access name recorded in OCAF TDataStd_Name attribute
        occ::handle<TCollection_HAsciiString> aHSubName = new TCollection_HAsciiString;
        if (!GetLabelName(aSubL, aHSubName))
        {
          continue;
        }

        // Access topological data
        TopoDS_Shape aSubS = XCAFDoc_ShapeTool::GetShape(aSubL);
        if (aSubS.IsNull())
        {
          continue;
        }

        // Access the correspondent STEP Representation Item
        occ::handle<StepRepr_RepresentationItem> aRI;
        occ::handle<TransferBRep_ShapeMapper>    aShMapper = TransferBRep::ShapeMapper(aFP, aSubS);
        if (!aFP->FindTypedTransient(aShMapper, STANDARD_TYPE(StepRepr_RepresentationItem), aRI))
        {
          continue;
        }

        // Record the name
        aRI->SetName(aHSubName);
      }
    }
  }

  return true;
}

//=================================================================================================

TopoDS_Shape STEPCAFControl_Writer::transferExternFiles(const TDF_Label&                 theLabel,
                                                        const STEPControl_StepModelType  theMode,
                                                        NCollection_Sequence<TDF_Label>& theLabels,
                                                        const StepData_Factors& theLocalFactors,
                                                        const char* const       thePrefix,
                                                        const Message_ProgressRange& theProgress)
{
  // if label already translated, just return the shape
  if (myLabels.IsBound(theLabel))
  {
    return myLabels.Find(theLabel);
  }

  TopoDS_Compound aComp;
  BRep_Builder    aBuilder;
  aBuilder.MakeCompound(aComp);
  // if not assembly, write to separate file
  if (!XCAFDoc_ShapeTool::IsAssembly(theLabel) && !XCAFDoc_ShapeTool::IsComponent(theLabel))
  {
    theLabels.Append(theLabel);
    // prepare for transfer
    occ::handle<XSControl_WorkSession> aNewWS = new XSControl_WorkSession;
    aNewWS->SelectNorm("STEP");
    STEPControl_Writer              aStepWriter(aNewWS, true);
    NCollection_Sequence<TDF_Label> aLabelSeq;
    aLabelSeq.Append(theLabel);

    // construct the name for extern file
    occ::handle<TCollection_HAsciiString> aBaseName = new TCollection_HAsciiString;
    if (thePrefix && thePrefix[0])
    {
      aBaseName->AssignCat(thePrefix);
    }
    GetLabelName(theLabel, aBaseName);
    occ::handle<TCollection_HAsciiString> aNewName = new TCollection_HAsciiString(aBaseName);
    aNewName->AssignCat(".stp");
    if (myFiles.IsBound(aNewName->ToCString()))
    { // avoid confusions
      for (int k = 1; k < 32000; k++)
      {
        aNewName = new TCollection_HAsciiString(aBaseName);
        aNewName->AssignCat("_");
        aNewName->AssignCat(TCollection_AsciiString(k).ToCString());
        aNewName->AssignCat(".stp");
        if (!myFiles.IsBound(aNewName->ToCString()))
        {
          break;
        }
      }
    }

    // translate and record extern file
    occ::handle<STEPCAFControl_ExternFile> anExtFile = new STEPCAFControl_ExternFile;
    anExtFile->SetWS(aNewWS);
    anExtFile->SetName(aNewName);
    anExtFile->SetLabel(theLabel);
    DESTEP_Parameters::WriteMode_Assembly anAssemblymode =
      aStepWriter.Model()->InternalParameters.WriteAssembly;
    aStepWriter.Model()->InternalParameters.WriteAssembly =
      DESTEP_Parameters::WriteMode_Assembly_Off;
    const char* const anIsMulti = nullptr;
    anExtFile->SetTransferStatus(
      transfer(aStepWriter, aLabelSeq, theMode, anIsMulti, true, theProgress));
    aStepWriter.Model()->InternalParameters.WriteAssembly = anAssemblymode;
    myLabEF.Bind(theLabel, anExtFile);
    myFiles.Bind(aNewName->ToCString(), anExtFile);

    // return empty compound as replacement for the shape
    myLabels.Bind(theLabel, aComp);
    return aComp;
  }
  NCollection_Sequence<TDF_Label> aCompLabels;
  TDF_Label                       aLabel = theLabel;
  // if specified shape is component then high-level assembly is considered
  // to get valid structure with location
  if (XCAFDoc_ShapeTool::IsComponent(theLabel))
  {
    aCompLabels.Append(theLabel);
    aLabel = theLabel.Father();
  }
  // else iterate on components add create structure of empty compounds
  // representing the assembly
  else if (XCAFDoc_ShapeTool::IsAssembly(theLabel))
  {
    XCAFDoc_ShapeTool::GetComponents(theLabel, aCompLabels, false);
  }

  theLabels.Append(aLabel);
  Message_ProgressScope aPS(theProgress, nullptr, aCompLabels.Length());
  // Iterate on requested shapes
  for (NCollection_Sequence<TDF_Label>::Iterator aLabelIter(aCompLabels);
       aLabelIter.More() && aPS.More();
       aLabelIter.Next())
  {
    const TDF_Label& aCurL = aLabelIter.Value();
    TDF_Label        aRefL;
    if (!XCAFDoc_ShapeTool::GetReferredShape(aCurL, aRefL))
    {
      continue;
    }
    TopoDS_Shape aShComp =
      transferExternFiles(aRefL, theMode, theLabels, theLocalFactors, thePrefix, aPS.Next());
    aShComp.Location(XCAFDoc_ShapeTool::GetLocation(aCurL));
    aBuilder.Add(aComp, aShComp);
  }
  myLabels.Bind(aLabel, aComp);
  return aComp;
}

//=================================================================================================

bool STEPCAFControl_Writer::writeExternRefs(const occ::handle<XSControl_WorkSession>& theWS,
                                            const NCollection_Sequence<TDF_Label>& theLabels) const
{
  if (theLabels.IsEmpty())
  {
    return false;
  }

  const occ::handle<XSControl_TransferWriter>& aTW = theWS->TransferWriter();
  const occ::handle<Transfer_FinderProcess>&   aFP = aTW->FinderProcess();
  STEPConstruct_ExternRefs                     anEFTool(theWS);
  occ::handle<StepData_StepModel> aStepModel  = occ::down_cast<StepData_StepModel>(theWS->Model());
  int                             aStepSchema = aStepModel->InternalParameters.WriteSchema;
  // Iterate on requested shapes
  for (NCollection_Sequence<TDF_Label>::Iterator aLabelIter(theLabels); aLabelIter.More();
       aLabelIter.Next())
  {
    const TDF_Label& aLab = aLabelIter.Value();
    if (XCAFDoc_ShapeTool::IsAssembly(aLab))
    {
      continue; // skip assemblies
    }

    // get extern file
    occ::handle<STEPCAFControl_ExternFile> anExtFile;
    if (!ExternFile(aLab, anExtFile))
    {
      continue; // should never be
    }

    // find SDR
    if (!myLabels.IsBound(aLab))
    {
      continue; // not recorded as translated, skip
    }
    TopoDS_Shape aShape = myLabels.Find(aLab);

    occ::handle<StepShape_ShapeDefinitionRepresentation> aSDR;
    occ::handle<TransferBRep_ShapeMapper> mapper = TransferBRep::ShapeMapper(aFP, aShape);
    if (!aFP->FindTypedTransient(mapper,
                                 STANDARD_TYPE(StepShape_ShapeDefinitionRepresentation),
                                 aSDR))
    {
      Message::SendTrace() << "Warning: Cannot find SDR for "
                           << aShape.TShape()->DynamicType()->Name() << "\n";
      continue;
    }

    // add extern ref
    const char* aStepFormat = (aStepSchema == 3 ? "STEP AP203" : "STEP AP214");
    // try to get PD from SDR
    StepRepr_RepresentedDefinition           aRD      = aSDR->Definition();
    occ::handle<StepRepr_PropertyDefinition> aPropDef = aRD.PropertyDefinition();
    if (aPropDef.IsNull())
    {
      Message::SendTrace() << "Warning: STEPCAFControl_Writer::writeExternRefs "
                              "StepRepr_PropertyDefinition is null for "
                           << aShape.TShape()->DynamicType()->Name() << "\n";
      continue;
    }
    StepRepr_CharacterizedDefinition         aCharDef = aPropDef->Definition();
    occ::handle<StepBasic_ProductDefinition> aPD      = aCharDef.ProductDefinition();
    if (aPD.IsNull())
    {
      Message::SendTrace() << "Warning: STEPCAFControl_Writer::writeExternRefs "
                              "StepBasic_ProductDefinition is null for "
                           << aShape.TShape()->DynamicType()->Name() << "\n";
      continue;
    }
    anEFTool.AddExternRef(anExtFile->GetName()->ToCString(), aPD, aStepFormat);
  }
  anEFTool.WriteExternRefs(aStepSchema);
  return true;
}

//=================================================================================================

static int FindEntities(const occ::handle<Transfer_FinderProcess>&             theFP,
                        const TopoDS_Shape&                                    theShape,
                        TopLoc_Location&                                       theLocation,
                        NCollection_Sequence<occ::handle<Standard_Transient>>& theSeqRI)
{
  occ::handle<StepRepr_RepresentationItem> anItem =
    STEPConstruct::FindEntity(theFP, theShape, theLocation);

  if (!anItem.IsNull())
  {
    theSeqRI.Append(anItem);
    return 1;
  }

  // may be S was split during shape processing
  occ::handle<TransferBRep_ShapeMapper> aMapper = TransferBRep::ShapeMapper(theFP, theShape);
  occ::handle<Transfer_Binder>          aBinder = theFP->Find(aMapper);
  if (aBinder.IsNull())
  {
    return 0;
  }

  occ::handle<Transfer_TransientListBinder> aTransientListBinder =
    // occ::down_cast<Transfer_TransientListBinder>( bnd->Next(true) );
    occ::down_cast<Transfer_TransientListBinder>(aBinder);
  int aResCount = 0;
  if (aTransientListBinder.IsNull())
  {
    for (TopoDS_Iterator anIter(theShape); anIter.More(); anIter.Next())
    {
      occ::handle<StepRepr_RepresentationItem> aLocalItem =
        STEPConstruct::FindEntity(theFP, anIter.Value(), theLocation);
      if (aLocalItem.IsNull())
      {
        continue;
      }
      aResCount++;
      theSeqRI.Append(aLocalItem);
    }
  }
  else if (!aTransientListBinder.IsNull())
  {
    const int aNbTransient = aTransientListBinder->NbTransients();
    for (int anInd = 1; anInd <= aNbTransient; anInd++)
    {
      occ::handle<Standard_Transient> anEntity = aTransientListBinder->Transient(anInd);
      anItem = occ::down_cast<StepRepr_RepresentationItem>(anEntity);
      if (anItem.IsNull())
      {
        continue;
      }
      aResCount++;
      theSeqRI.Append(anItem);
    }
  }
  return aResCount;
}

//=================================================================================================

static bool getStyledItem(
  const TopoDS_Shape&                   theShape,
  const occ::handle<XCAFDoc_ShapeTool>& theShapeTool,
  const STEPConstruct_Styles&           theStyles,
  occ::handle<StepVisual_StyledItem>&   theResSelItem,
  const NCollection_DataMap<TopoDS_Shape, occ::handle<Standard_Transient>, TopTools_ShapeMapHasher>&
    theMapCompMDGPR)
{
  const TDF_Label aTopShL   = theShapeTool->FindShape(theShape, false);
  TopoDS_Shape    aTopLevSh = theShapeTool->GetShape(aTopShL);
  bool            anIsFound = false;
  if (aTopLevSh.IsNull() || !theMapCompMDGPR.IsBound(aTopLevSh))
  {
    return false;
  }
  occ::handle<StepVisual_PresentationRepresentation> aMDGPR =
    occ::down_cast<StepVisual_PresentationRepresentation>(theMapCompMDGPR.Find(aTopLevSh));
  occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> anSelItmHArr =
    aMDGPR->Items();
  if (anSelItmHArr.IsNull())
  {
    return false;
  }
  // Search for PSA of Manifold solid
  NCollection_Sequence<occ::handle<Standard_Transient>> aNewSeqRI;
  bool                                                  isFilled = false;
  for (NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>::Iterator anIter(
         anSelItmHArr->Array1());
       anIter.More() && !anIsFound;
       anIter.Next())
  {
    const occ::handle<StepVisual_StyledItem>& aSelItm =
      occ::down_cast<StepVisual_StyledItem>(anIter.Value());

    if (aSelItm.IsNull())
    {
      continue;
    }
    // Check that it is a styled item for manifold solid brep
    if (!isFilled)
    {
      TopLoc_Location aLoc;
      FindEntities(theStyles.FinderProcess(), aTopLevSh, aLoc, aNewSeqRI);
      isFilled = true;
    }
    if (!aNewSeqRI.IsEmpty())
    {
      const occ::handle<StepRepr_RepresentationItem>& anItem         = aSelItm->Item();
      bool                                            isSameMonSolBR = false;
      for (NCollection_Sequence<occ::handle<Standard_Transient>>::Iterator aIterRI(aNewSeqRI);
           aIterRI.More();
           aIterRI.Next())
      {
        if (!anItem.IsNull() && anItem == aIterRI.Value())
        {
          isSameMonSolBR = true;
          break;
        }
      }
      if (!isSameMonSolBR)
      {
        continue;
      }
    }
    for (NCollection_HArray1<occ::handle<StepVisual_PresentationStyleAssignment>>::Iterator
           aStyleIter(aSelItm->Styles()->Array1());
         aStyleIter.More() && !anIsFound;
         aStyleIter.Next())
    {
      const occ::handle<StepVisual_PresentationStyleAssignment>& aFatherPSA = aStyleIter.Value();
      // check for PSA for top-level (not Presentation style by context for NAUO)
      if (aFatherPSA.IsNull()
          || aFatherPSA->IsKind(STANDARD_TYPE(StepVisual_PresentationStyleByContext)))
      {
        continue;
      }
      theResSelItem = aSelItm;
      anIsFound     = true;
    }
  }
  return anIsFound;
}

//=================================================================================================

static bool setDefaultInstanceColor(const occ::handle<StepVisual_StyledItem>& theStyleItem,
                                    occ::handle<StepVisual_PresentationStyleAssignment>& thePSA)
{
  // No styled item of the referred shape to copy from -- getStyledItem() finds
  // none for the components of an assembly written in one pass, which come
  // before the shape they refer to: keep the default PSA as it is.
  if (theStyleItem.IsNull() || theStyleItem->Styles().IsNull())
  {
    return false;
  }
  bool anIsFound = false;
  for (NCollection_HArray1<occ::handle<StepVisual_PresentationStyleAssignment>>::Iterator
         aStyleIter(theStyleItem->Styles()->Array1());
       aStyleIter.More() && !anIsFound;
       aStyleIter.Next())
  {
    const occ::handle<StepVisual_PresentationStyleAssignment>& aFatherPSA = aStyleIter.Value();
    // check for PSA for top-level (not Presentation style by context for NAUO)
    if (aFatherPSA.IsNull()
        || aFatherPSA->IsKind(STANDARD_TYPE(StepVisual_PresentationStyleByContext)))
    {
      return false;
    }

    // get style select from father PSA
    if (aFatherPSA->NbStyles() > 0)
    {
      occ::handle<NCollection_HArray1<StepVisual_PresentationStyleSelect>> aFatherStyles =
        new NCollection_HArray1<StepVisual_PresentationStyleSelect>(1, aFatherPSA->NbStyles());
      int aSettingInd = 1;
      for (NCollection_HArray1<StepVisual_PresentationStyleSelect>::Iterator aFatherStyleIter(
             aFatherPSA->Styles()->Array1());
           aFatherStyleIter.More();
           aFatherStyleIter.Next())
      {
        StepVisual_PresentationStyleSelect        aPSS;
        const StepVisual_PresentationStyleSelect& anOlDPSS = aFatherStyleIter.Value();
        if (!anOlDPSS.PointStyle().IsNull())
        {
          aPSS.SetValue(anOlDPSS.PointStyle());
        }
        else if (!anOlDPSS.CurveStyle().IsNull())
        {
          aPSS.SetValue(anOlDPSS.CurveStyle());
        }
        else if (!anOlDPSS.SurfaceStyleUsage().IsNull())
        {
          aPSS.SetValue(anOlDPSS.SurfaceStyleUsage());
        }
        else
        {
          anIsFound = false;
          break;
        }
        aFatherStyles->SetValue(aSettingInd++, anOlDPSS);
        anIsFound = true;
      }
      // init PSA of NAUO
      if (anIsFound)
      {
        thePSA->Init(aFatherStyles);
      }
    }
  }
  return anIsFound;
}

//=================================================================================================

static void MakeSTEPStyles(
  STEPConstruct_Styles& theStyles,
  const TopoDS_Shape&   theShape,
  const NCollection_IndexedDataMap<TopoDS_Shape, XCAFPrs_Style, TopTools_ShapeMapHasher>&
                                                          theSettings,
  occ::handle<StepVisual_StyledItem>&                     theOverride,
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& theMap,
  const NCollection_DataMap<TopoDS_Shape, occ::handle<Standard_Transient>, TopTools_ShapeMapHasher>&
                                                                                 theMapCompMDGPR,
  NCollection_DataMap<TCollection_AsciiString, occ::handle<Standard_Transient>>& theDPDCs,
  NCollection_DataMap<gp_Pnt, occ::handle<Standard_Transient>>&                  theColRGBs,
  const occ::handle<XCAFDoc_ShapeTool>&                                          theShTool,
  const XCAFPrs_Style*                                                           theInherit,
  const bool                                                                     theIsComponent,
  const bool                                                                     theVisMaterialMode,
  NCollection_Sequence<occ::handle<StepVisual_StyledItem>>*                      theHidden = nullptr)
{
  // skip already processed shapes
  if (!theMap.Add(theShape))
  {
    return;
  }

  // check if shape has its own style (r inherits from ancestor)
  XCAFPrs_Style aStyle;
  if (theInherit)
  {
    aStyle = *theInherit;
  }
  if (theSettings.Contains(theShape))
  {
    const XCAFPrs_Style& anOwnStyle = theSettings.FindFromKey(theShape);
    if (!anOwnStyle.IsVisible())
    {
      aStyle.SetVisibility(false);
    }
    if (anOwnStyle.IsSetColorCurv())
    {
      aStyle.SetColorCurv(anOwnStyle.GetColorCurv());
    }
    if (anOwnStyle.IsSetColorSurf())
    {
      aStyle.SetColorSurf(anOwnStyle.GetColorSurfRGBA());
    }
    if (!anOwnStyle.Material().IsNull())
    {
      aStyle.SetMaterial(anOwnStyle.Material());
    }
  }

  // translate colors to STEP
  occ::handle<StepVisual_Colour> aSurfColor, aCurvColor;
  double                         aRenderTransp = 0.0;
  if (aStyle.IsSetColorSurf())
  {
    Quantity_ColorRGBA aSurfCol = aStyle.GetColorSurfRGBA();
    aRenderTransp               = 1.0 - aSurfCol.Alpha();
    aSurfColor = STEPConstruct_Styles::EncodeColor(aSurfCol.GetRGB(), theDPDCs, theColRGBs);
  }
  if (aStyle.IsSetColorCurv())
  {
    aCurvColor = STEPConstruct_Styles::EncodeColor(aStyle.GetColorCurv(), theDPDCs, theColRGBs);
  }

  bool aHasOwn = (!aSurfColor.IsNull() || !aCurvColor.IsNull() || !aStyle.IsVisible());

  // find target item and assign style to it
  occ::handle<StepVisual_StyledItem> aSTEPstyle = theOverride;
  if (aHasOwn)
  {
    if (theShape.ShapeType() != TopAbs_COMPOUND || theIsComponent)
    { // skip compounds, let subshapes inherit its colors
      TopLoc_Location                                       aLocation;
      NCollection_Sequence<occ::handle<Standard_Transient>> aSeqRI;
      int aNbEntities = FindEntities(theStyles.FinderProcess(), theShape, aLocation, aSeqRI);
      if (aNbEntities <= 0)
      {
        Message::SendTrace() << "Warning: Cannot find RI for "
                             << theShape.TShape()->DynamicType()->Name() << "\n";
      }
      // Get overridden style gka 10.06.03
      if (theIsComponent && aNbEntities > 0)
      {
        getStyledItem(theShape, theShTool, theStyles, theOverride, theMapCompMDGPR);
      }

      for (NCollection_Sequence<occ::handle<Standard_Transient>>::Iterator anEntIter(aSeqRI);
           anEntIter.More();
           anEntIter.Next())
      {
        const occ::handle<StepRepr_RepresentationItem>& anItem =
          occ::down_cast<StepRepr_RepresentationItem>(anEntIter.Value());
        occ::handle<StepVisual_PresentationStyleAssignment> aPSA;
        if (aStyle.IsVisible() || !aSurfColor.IsNull() || !aCurvColor.IsNull()
            || (theVisMaterialMode && !aStyle.Material().IsNull() && !aStyle.Material()->IsEmpty()))
        {
          STEPConstruct_RenderingProperties aRenderProps;
          if (theVisMaterialMode && !aStyle.Material().IsNull() && !aStyle.Material()->IsEmpty())
          {
            aRenderProps.Init(aStyle.Material());
          }
          else if (aRenderTransp > 0.0)
          {
            aRenderProps.Init(aStyle.GetColorSurfRGBA());
          }
          aPSA =
            theStyles.MakeColorPSA(anItem, aSurfColor, aCurvColor, aRenderProps, theIsComponent);
        }
        else
        {
          // default white color
          aSurfColor = STEPConstruct_Styles::EncodeColor(Quantity_Color(Quantity_NOC_WHITE),
                                                         theDPDCs,
                                                         theColRGBs);
          aPSA       = theStyles.MakeColorPSA(anItem,
                                        aSurfColor,
                                        aCurvColor,
                                        STEPConstruct_RenderingProperties(),
                                        theIsComponent);
          if (theIsComponent)
          {
            setDefaultInstanceColor(theOverride, aPSA);
          }

        } // end of component case

        aSTEPstyle = theStyles.AddStyle(anItem, aPSA, theOverride);
        // A part of the shape hidden alone: its style listed by an invisibility
        if (theHidden && !aStyle.IsVisible() && !aSTEPstyle.IsNull())
        {
          theHidden->Append(aSTEPstyle);
        }
        aHasOwn = false;
      }
    }
  }

  // iterate on subshapes (except vertices :)
  if (theShape.ShapeType() == TopAbs_EDGE)
  {
    return;
  }
  if (theIsComponent)
  {
    return;
  }
  for (TopoDS_Iterator anIter(theShape); anIter.More(); anIter.Next())
  {
    MakeSTEPStyles(theStyles,
                   anIter.Value(),
                   theSettings,
                   aSTEPstyle,
                   theMap,
                   theMapCompMDGPR,
                   theDPDCs,
                   theColRGBs,
                   theShTool,
                   (aHasOwn ? &aStyle : nullptr),
                   false,
                   theVisMaterialMode,
                   theHidden);
  }
}

//=================================================================================================

bool STEPCAFControl_Writer::writeColors(const occ::handle<XSControl_WorkSession>& theWS,
                                        const NCollection_Sequence<TDF_Label>&    theLabels)
{
  if (theLabels.IsEmpty())
  {
    return false;
  }

  STEPConstruct_Styles                                                          Styles(theWS);
  NCollection_DataMap<TCollection_AsciiString, occ::handle<Standard_Transient>> DPDCs;
  NCollection_DataMap<gp_Pnt, occ::handle<Standard_Transient>>                  ColRGBs;
  // Iterate on requested shapes
  for (NCollection_Sequence<TDF_Label>::Iterator aLabelIter(theLabels); aLabelIter.More();
       aLabelIter.Next())
  {
    const TDF_Label aLabel = aLabelIter.Value();
    // Iterate on shapes in the document
    occ::handle<XCAFDoc_ShapeTool> aSTool = XCAFDoc_DocumentTool::ShapeTool(aLabel);
    // Skip assemblies: colors assigned to assemblies and their instances
    // are not supported (it is not clear how to encode that in STEP)
    if (XCAFDoc_ShapeTool::IsAssembly(aLabel))
    {
      Message::SendTrace() << "Warning: Cannot write color  for Assembly" << "\n";
      Message::SendTrace() << "Info: Check for colors assigned to components in assembly" << "\n";
      // PTV 22.01.2003 Write color for instances.
      NCollection_Sequence<TDF_Label> compLabels;
      if (!aSTool->GetComponents(aLabel, compLabels))
      {
        continue;
      }
      writeColors(theWS, compLabels);
      continue;
    }
    Styles.ClearStyles();

    // get a target shape and try to find corresponding context
    // (all the colors set under that label will be put into that context)
    TopoDS_Shape aShape;
    if (!XCAFDoc_ShapeTool::GetShape(aLabel, aShape))
    {
      continue;
    }
    bool         anIsComponent = aSTool->IsComponent(aLabel) || myPureRefLabels.IsBound(aLabel);
    TopoDS_Shape aTopSh        = aShape;
    occ::handle<StepRepr_RepresentationContext> aContext = Styles.FindContext(aShape);
    if (anIsComponent)
    {
      TDF_Label aTopShL = aSTool->FindShape(aShape, false);
      if (aTopShL.IsNull())
      {
        continue;
      }
      aTopSh   = aSTool->GetShape(aTopShL);
      aContext = Styles.FindContext(aTopSh);
    }
    if (aContext.IsNull())
    {
      continue;
    }

    // collect settings set on that label
    NCollection_IndexedDataMap<TopoDS_Shape, XCAFPrs_Style, TopTools_ShapeMapHasher> aSettings;
    NCollection_Sequence<TDF_Label>                                                  aSeq;
    aSeq.Append(aLabel);
    XCAFDoc_ShapeTool::GetSubShapes(aLabel, aSeq);
    bool anIsVisible = true;
    for (NCollection_Sequence<TDF_Label>::Iterator aSeqIter(aSeq); aSeqIter.More(); aSeqIter.Next())
    {
      const TDF_Label&   aSeqValue = aSeqIter.Value();
      XCAFPrs_Style      aStyle;
      Quantity_ColorRGBA aColor;
      // check for invisible status of object on label, and of a part of it hidden alone
      const bool isSubVisible = XCAFDoc_ColorTool::IsVisible(aSeqValue);
      if (!isSubVisible)
      {
        aStyle.SetVisibility(false);
        if (aSeqValue == aLabel)
        {
          anIsVisible = false;
        }
      }
      if (XCAFDoc_ColorTool::GetColor(aSeqValue, XCAFDoc_ColorGen, aColor))
      {
        aStyle.SetColorCurv(aColor.GetRGB());
        aStyle.SetColorSurf(aColor);
      }
      if (XCAFDoc_ColorTool::GetColor(aSeqValue, XCAFDoc_ColorSurf, aColor))
      {
        aStyle.SetColorSurf(aColor);
      }
      if (XCAFDoc_ColorTool::GetColor(aSeqValue, XCAFDoc_ColorCurv, aColor))
      {
        aStyle.SetColorCurv(aColor.GetRGB());
      }
      {
        occ::handle<XCAFDoc_VisMaterial> aVisMat =
          XCAFDoc_VisMaterialTool::GetShapeMaterial(aSeqValue);
        if (!aVisMat.IsNull() && !aVisMat->IsEmpty())
        {
          // the material rides the style: written out as rendering
          // properties when the visual material mode asks for it, and
          // its base colour keeps serving readers that only take colours
          aStyle.SetMaterial(aVisMat);
          if (!aStyle.IsSetColorSurf())
          {
            aStyle.SetColorSurf(aVisMat->BaseColor());
          }
        }
      }
      if (!aStyle.IsSetColorCurv() && !aStyle.IsSetColorSurf() && aStyle.Material().IsNull()
          && isSubVisible)
      {
        continue;
      }

      TopoDS_Shape   aSub      = XCAFDoc_ShapeTool::GetShape(aSeqValue);
      XCAFPrs_Style* aMapStyle = aSettings.ChangeSeek(aSub);
      if (aMapStyle == nullptr)
      {
        aSettings.Add(aSub, aStyle);
      }
      else
      {
        *aMapStyle = aStyle;
      }
    }

    if (aSettings.Extent() <= 0)
    {
      continue;
    }

    // iterate on subshapes and create STEP styles
    occ::handle<StepVisual_StyledItem>                     anOverride;
    NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aMap;

    NCollection_Sequence<occ::handle<StepVisual_StyledItem>> aHidden;
    MakeSTEPStyles(Styles,
                   aShape,
                   aSettings,
                   anOverride,
                   aMap,
                   myMapCompMDGPR,
                   DPDCs,
                   ColRGBs,
                   aSTool,
                   nullptr,
                   anIsComponent,
                   GetVisualMaterialMode(),
                   &aHidden);

    const occ::handle<XSControl_TransferWriter>& aTW = theWS->TransferWriter();
    const occ::handle<Transfer_FinderProcess>&   aFP = aTW->FinderProcess();
    occ::handle<StepData_StepModel> aStepModel = occ::down_cast<StepData_StepModel>(aFP->Model());

    // create MDGPR and record it in model
    occ::handle<StepVisual_MechanicalDesignGeometricPresentationRepresentation> aMDGPR;

    if (!anIsComponent)
    {
      if (myMapCompMDGPR.IsBound(aTopSh))
      {
        Message::SendTrace() << "Error: Current Top-Level shape have MDGPR already " << "\n";
      }
      Styles.CreateMDGPR(aContext, aMDGPR, aStepModel);
      if (!aMDGPR.IsNull())
      {
        myMapCompMDGPR.Bind(aTopSh, aMDGPR);
      }
    }
    else
    {
      // create SDR and add to model.
      occ::handle<TransferBRep_ShapeMapper> aMapper = TransferBRep::ShapeMapper(aFP, aShape);
      occ::handle<StepShape_ContextDependentShapeRepresentation> aCDSR;
      if (aFP->FindTypedTransient(aMapper,
                                  STANDARD_TYPE(StepShape_ContextDependentShapeRepresentation),
                                  aCDSR))
      {
        // create SDR for NAUO
        occ::handle<StepRepr_ProductDefinitionShape> aNullPDS; // important to be NULL
        Styles.CreateNAUOSRD(aContext, aCDSR, aNullPDS);

        // search for MDGPR of the component top-level shape
        if (myMapCompMDGPR.IsBound(aTopSh))
        {
          aMDGPR = occ::down_cast<StepVisual_MechanicalDesignGeometricPresentationRepresentation>(
            myMapCompMDGPR.Find(aTopSh));
        }
        else
        {
          aMDGPR = new StepVisual_MechanicalDesignGeometricPresentationRepresentation;
          occ::handle<TCollection_HAsciiString> aReprName = new TCollection_HAsciiString("");
          aMDGPR->SetName(aReprName);
          aMDGPR->SetContextOfItems(aContext);
          myMapCompMDGPR.Bind(aTopSh, aMDGPR);
        }
        occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> anOldItems =
          aMDGPR->Items();
        int oldLengthlen = 0;
        if (!anOldItems.IsNull())
        {
          oldLengthlen = anOldItems->Length();
        }
        const int aNbIt = oldLengthlen + Styles.NbStyles();
        if (!aNbIt)
        {
          continue;
        }
        occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> aNewItems =
          new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, aNbIt);
        int anElemInd = 1;
        for (int aStyleInd = 1; aStyleInd <= oldLengthlen; aStyleInd++)
        {
          aNewItems->SetValue(anElemInd++, anOldItems->Value(aStyleInd));
        }
        for (int aStyleInd = 1; aStyleInd <= Styles.NbStyles(); aStyleInd++)
        {
          aNewItems->SetValue(anElemInd++, Styles.Style(aStyleInd));
        }

        if (aNewItems->Length() > 0)
        {
          aMDGPR->SetItems(aNewItems);
        }
      } // end of work with CDSR
    }
    if (anIsVisible && !aHidden.IsEmpty())
    {
      // The parts hidden alone (NIST STC-08's surface beside its solid): their styles
      occ::handle<NCollection_HArray1<StepVisual_InvisibleItem>> anItems =
        new NCollection_HArray1<StepVisual_InvisibleItem>(1, aHidden.Length());
      for (int i = 1; i <= aHidden.Length(); ++i)
      {
        StepVisual_InvisibleItem anInvItem;
        anInvItem.SetValue(aHidden.Value(i));
        anItems->SetValue(i, anInvItem);
      }
      occ::handle<StepVisual_Invisibility> anInvisibility = new StepVisual_Invisibility();
      anInvisibility->Init(anItems);
      theWS->Model()->AddWithRefs(anInvisibility);
    }
    if (!anIsVisible)
    {
      // create invisibility item and refer for styledItem
      occ::handle<StepVisual_Invisibility> anInvisibility = new StepVisual_Invisibility();
      occ::handle<NCollection_HArray1<StepVisual_InvisibleItem>> anInvisibilitySeq =
        new NCollection_HArray1<StepVisual_InvisibleItem>(1, Styles.NbStyles());
      // put all style item into the harray
      for (int aStyleInd = 1; aStyleInd <= Styles.NbStyles(); aStyleInd++)
      {
        occ::handle<StepRepr_RepresentationItem> aStyledItm = Styles.Style(aStyleInd);
        StepVisual_InvisibleItem                 anInvItem;
        anInvItem.SetValue(aStyledItm);
        anInvisibilitySeq->SetValue(aStyleInd, anInvItem);
      }
      // set the invisibility of items
      anInvisibility->Init(anInvisibilitySeq);
      theWS->Model()->AddWithRefs(anInvisibility);
    }
  }
  return true;
}

//=================================================================================================

bool STEPCAFControl_Writer::writeNames(const occ::handle<XSControl_WorkSession>& theWS,
                                       const NCollection_Sequence<TDF_Label>&    theLabels) const
{
  if (theLabels.IsEmpty())
  {
    return false;
  }

  // get working data
  const occ::handle<XSControl_TransferWriter>& aTW = theWS->TransferWriter();
  const occ::handle<Transfer_FinderProcess>&   aFP = aTW->FinderProcess();

  // Iterate on requested shapes
  for (NCollection_Sequence<TDF_Label>::Iterator aLabelIter(theLabels); aLabelIter.More();
       aLabelIter.Next())
  {
    const TDF_Label& aLabel = aLabelIter.Value();
    // find target STEP entity for the current shape
    if (!myLabels.IsBound(aLabel))
    {
      continue; // not recorded as translated, skip
    }
    // get name
    occ::handle<TCollection_HAsciiString> aHName = new TCollection_HAsciiString;
    if (!GetLabelName(aLabel, aHName))
    {
      continue;
    }
    const TopoDS_Shape&                                        aShape = myLabels.Find(aLabel);
    occ::handle<StepShape_ShapeDefinitionRepresentation>       aSDR;
    occ::handle<StepShape_ContextDependentShapeRepresentation> aCDSR;
    bool isComponent = XCAFDoc_ShapeTool::IsComponent(aLabel) || myPureRefLabels.IsBound(aLabel);
    occ::handle<TransferBRep_ShapeMapper> aMapper = TransferBRep::ShapeMapper(aFP, aShape);
    if (isComponent
        && aFP->FindTypedTransient(aMapper,
                                   STANDARD_TYPE(StepShape_ContextDependentShapeRepresentation),
                                   aCDSR))
    {
      occ::handle<StepRepr_ProductDefinitionShape> aPDS = aCDSR->RepresentedProductRelation();
      occ::handle<StepBasic_ProductDefinitionRelationship> aNAUO =
        aPDS->Definition().ProductDefinitionRelationship();
      if (!aNAUO.IsNull())
      {
        aNAUO->SetName(aHName);
      }
      TopoDS_Shape anInternalAssembly;
      if (myPureRefLabels.Find(aLabel, anInternalAssembly))
      {
        occ::handle<TransferBRep_ShapeMapper> aMapperOfInternalShape =
          TransferBRep::ShapeMapper(aFP, anInternalAssembly);
        aFP->FindTypedTransient(aMapperOfInternalShape,
                                STANDARD_TYPE(StepShape_ShapeDefinitionRepresentation),
                                aSDR);
      }
    }
    if (!aSDR.IsNull()
        || (aCDSR.IsNull()
            && aFP->FindTypedTransient(aMapper,
                                       STANDARD_TYPE(StepShape_ShapeDefinitionRepresentation),
                                       aSDR)))
    {
      // set the name to the PRODUCT
      occ::handle<StepRepr_PropertyDefinition> aPropD = aSDR->Definition().PropertyDefinition();
      if (aPropD.IsNull())
      {
        continue;
      }
      occ::handle<StepBasic_ProductDefinition> aPD = aPropD->Definition().ProductDefinition();
      if (aPD.IsNull())
      {
        continue;
      }
      occ::handle<StepBasic_Product> aProd = aPD->Formation()->OfProduct();

      aProd->SetId(aHName);
      aProd->SetName(aHName);
    }
    else
    {
      Message::SendTrace() << "Warning: Cannot find RI for "
                           << aShape.TShape()->DynamicType()->Name() << "\n";
      continue;
    }
  }
  return true;
}

//=================================================================================================

bool STEPCAFControl_Writer::writeMetadata(const occ::handle<XSControl_WorkSession>& theWS,
                                          const NCollection_Sequence<TDF_Label>&    theLabels) const
{
  if (theLabels.IsEmpty())
  {
    return false;
  }

  // Iterate on requested shapes.
  for (NCollection_Sequence<TDF_Label>::Iterator aLabelIter(theLabels); aLabelIter.More();
       aLabelIter.Next())
  {
    writeMetadataForLabel(theWS, aLabelIter.Value());
  }

  return true;
}

//=================================================================================================

bool STEPCAFControl_Writer::writeMetadataForLabel(const occ::handle<XSControl_WorkSession>& theWS,
                                                  const TDF_Label& theLabel) const
{
  if (theLabel.IsNull())
  {
    return false;
  }

  // Get working data:
  const occ::handle<XSControl_TransferWriter>& aTW = theWS->TransferWriter();
  const occ::handle<Transfer_FinderProcess>&   aFP = aTW->FinderProcess();
  const occ::handle<StepData_StepModel>&       aModel =
    occ::down_cast<StepData_StepModel>(theWS->Model());

  // Check if label has metadata (NamedData):
  occ::handle<TDataStd_NamedData> aNamedData;
  if (!theLabel.FindAttribute(TDataStd_NamedData::GetID(), aNamedData))
  {
    return false; // No metadata on this label
  }

  // Find target STEP entity for the current shape:
  if (!myLabels.IsBound(theLabel))
  {
    return false; // Not recorded as translated, skip
  }

  const TopoDS_Shape&                                  aShape = myLabels.Find(theLabel);
  occ::handle<StepShape_ShapeDefinitionRepresentation> aSDR;
  const occ::handle<TransferBRep_ShapeMapper> aMapper = TransferBRep::ShapeMapper(aFP, aShape);
  if (!aFP->FindTypedTransient(aMapper,
                               STANDARD_TYPE(StepShape_ShapeDefinitionRepresentation),
                               aSDR))
  {
    return false; // Cannot find STEP representation
  }

  // Get the product definition from the shape definition representation:
  const occ::handle<StepRepr_PropertyDefinition> aPropDef = aSDR->Definition().PropertyDefinition();
  if (aPropDef.IsNull())
  {
    return false;
  }
  const occ::handle<StepBasic_ProductDefinition> aProdDef =
    aPropDef->Definition().ProductDefinition();
  if (aProdDef.IsNull())
  {
    return false;
  }

  // Export string metadata.
  const NCollection_DataMap<TCollection_ExtendedString, TCollection_ExtendedString>& aStringMap =
    aNamedData->GetStringsContainer();
  for (NCollection_DataMap<TCollection_ExtendedString, TCollection_ExtendedString>::Iterator anIter(
         aStringMap);
       anIter.More();
       anIter.Next())
  {
    const TCollection_ExtendedString& aKey   = anIter.Key();
    const TCollection_ExtendedString& aValue = anIter.Value();

    // Create descriptive representation item for the value.
    const occ::handle<StepRepr_DescriptiveRepresentationItem> aDescrItem =
      new StepRepr_DescriptiveRepresentationItem();
    const occ::handle<TCollection_HAsciiString> aItemName  = new TCollection_HAsciiString(aKey);
    const occ::handle<TCollection_HAsciiString> aItemValue = new TCollection_HAsciiString(aValue);
    aDescrItem->SetName(aItemName);
    aDescrItem->SetDescription(aItemValue);

    writeMetadataRepresentationItem(aKey, aModel, aSDR, aProdDef, aDescrItem);
  }

  // Export integer metadata.
  const NCollection_DataMap<TCollection_ExtendedString, int>& aIntMap =
    aNamedData->GetIntegersContainer();
  for (NCollection_DataMap<TCollection_ExtendedString, int>::Iterator anIter(aIntMap);
       anIter.More();
       anIter.Next())
  {
    const TCollection_ExtendedString& aKey   = anIter.Key();
    const int                         aValue = anIter.Value();
    // Create integer representation item for the value.
    const occ::handle<StepRepr_IntegerRepresentationItem> aIntItem =
      new StepRepr_IntegerRepresentationItem();
    const occ::handle<TCollection_HAsciiString> aItemName = new TCollection_HAsciiString(aKey);
    aIntItem->Init(aItemName, aValue);

    writeMetadataRepresentationItem(aKey, aModel, aSDR, aProdDef, aIntItem);
  }

  // Export real metadata.
  const NCollection_DataMap<TCollection_ExtendedString, double>& aRealMap =
    aNamedData->GetRealsContainer();
  for (NCollection_DataMap<TCollection_ExtendedString, double>::Iterator anIter(aRealMap);
       anIter.More();
       anIter.Next())
  {
    const TCollection_ExtendedString& aKey   = anIter.Key();
    const double                      aValue = anIter.Value();
    // Create real representation item for the value.
    const occ::handle<StepRepr_RealRepresentationItem> aRealItem =
      new StepRepr_RealRepresentationItem();
    const occ::handle<TCollection_HAsciiString> aItemName = new TCollection_HAsciiString(aKey);
    aRealItem->Init(aItemName, aValue);

    writeMetadataRepresentationItem(aKey, aModel, aSDR, aProdDef, aRealItem);
  }

  // Process label children recursively:
  if (theLabel.HasChild())
  {
    for (int aChildInd = 1; aChildInd <= theLabel.NbChildren(); aChildInd++)
    {
      const TDF_Label& aChildLabel = theLabel.FindChild(aChildInd);
      writeMetadataForLabel(theWS, aChildLabel);
    }
  }

  return true;
}

//=================================================================================================

void STEPCAFControl_Writer::writeMetadataRepresentationItem(
  const TCollection_AsciiString&                              theKey,
  const occ::handle<StepData_StepModel>&                      theModel,
  const occ::handle<StepShape_ShapeDefinitionRepresentation>& theShapeDefRep,
  const occ::handle<StepBasic_ProductDefinition>&             theProdDef,
  const occ::handle<StepRepr_RepresentationItem>&             theItem) const
{
  // Empty string to use for empty values:
  const occ::handle<TCollection_HAsciiString> anEmptyStr = new TCollection_HAsciiString("");

  // Create property_definition:
  const occ::handle<StepRepr_PropertyDefinition> aMetaPropDef = new StepRepr_PropertyDefinition();
  const occ::handle<TCollection_HAsciiString>    aPropName = new TCollection_HAsciiString(theKey);
  const occ::handle<TCollection_HAsciiString>    aPropDesc =
    new TCollection_HAsciiString("user defined attribute");
  aMetaPropDef->SetName(aPropName);
  aMetaPropDef->SetDescription(aPropDesc);

  // Create a general_property:
  const occ::handle<StepBasic_GeneralProperty> aGeneralProp = new StepBasic_GeneralProperty();
  aGeneralProp->SetId(anEmptyStr);
  aGeneralProp->SetName(aPropName);
  aGeneralProp->SetDescription(anEmptyStr);

  // Create a general_property_association:
  const occ::handle<StepBasic_GeneralPropertyAssociation> aGeneralPropAssoc =
    new StepBasic_GeneralPropertyAssociation();
  aGeneralPropAssoc->SetName(anEmptyStr);
  aGeneralPropAssoc->SetDescription(anEmptyStr);
  aGeneralPropAssoc->SetPropertyDefinition(aMetaPropDef);
  aGeneralPropAssoc->SetGeneralProperty(aGeneralProp);

  // Set the definition to point to the product definition:
  StepRepr_CharacterizedDefinition aCharDef;
  aCharDef.SetValue(theProdDef);
  aMetaPropDef->SetDefinition(aCharDef);

  // Create property_definition_representation:
  const occ::handle<StepRepr_PropertyDefinitionRepresentation> aPropDefRepr =
    new StepRepr_PropertyDefinitionRepresentation();
  StepRepr_RepresentedDefinition aRepDef;
  aRepDef.SetValue(aMetaPropDef);
  aPropDefRepr->SetDefinition(aRepDef);

  // Create representation with descriptive_representation_item:
  const occ::handle<StepRepr_Representation>  aRepr     = new StepRepr_Representation();
  const occ::handle<TCollection_HAsciiString> aReprName = new TCollection_HAsciiString(theKey);
  aRepr->SetName(aReprName);

  // Add item to representation:
  const occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> aItems =
    new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, 1);
  aItems->SetValue(1, theItem);
  aRepr->SetItems(aItems);

  // Set representation context (reuse from shape representation):
  const occ::handle<StepRepr_RepresentationContext> aRC =
    theShapeDefRep->UsedRepresentation()->ContextOfItems();
  aRepr->SetContextOfItems(aRC);

  aPropDefRepr->SetUsedRepresentation(aRepr);

  // Add entities to the model:
  theModel->AddWithRefs(aMetaPropDef);
  theModel->AddWithRefs(aGeneralProp);
  theModel->AddWithRefs(aGeneralPropAssoc);
  theModel->AddWithRefs(aPropDefRepr);
  theModel->AddWithRefs(aRepr);
  theModel->AddWithRefs(theItem);
}

//=================================================================================================

static bool WritePropsForLabel(const occ::handle<XSControl_WorkSession>&           theWS,
                               const NCollection_DataMap<TDF_Label, TopoDS_Shape>& theLabels,
                               const TDF_Label&                                    theLabel,
                               const char* const                                   theIsMulti)
{
  if (theLabel.IsNull())
  {
    return false;
  }

  STEPConstruct_ValidationProps aProps(theWS);

  TopoDS_Shape aShape = XCAFDoc_ShapeTool::GetShape(theLabel);
  if (aShape.IsNull())
  {
    return false;
  }

  if (!theIsMulti || XCAFDoc_ShapeTool::IsAssembly(theLabel))
  {
    // write area
    occ::handle<XCAFDoc_Area> anArea;
    theLabel.FindAttribute(XCAFDoc_Area::GetID(), anArea);
    if (!anArea.IsNull())
    {
      aProps.AddArea(aShape, anArea->Get());
    }
    // write volume
    occ::handle<XCAFDoc_Volume> aVolume;
    theLabel.FindAttribute(XCAFDoc_Volume::GetID(), aVolume);
    if (!aVolume.IsNull())
    {
      aProps.AddVolume(aShape, aVolume->Get());
    }
  }
  // write centroid
  occ::handle<XCAFDoc_Centroid> aCentroid;
  theLabel.FindAttribute(XCAFDoc_Centroid::GetID(), aCentroid);
  if (!aCentroid.IsNull())
  {
    aProps.AddCentroid(aShape, aCentroid->Get());
  }

  if (XCAFDoc_ShapeTool::IsCompound(theLabel) || XCAFDoc_ShapeTool::IsAssembly(theLabel))
  {
    if (theLabel.HasChild())
    {
      for (int aChildInd = 1; aChildInd <= theLabel.NbChildren(); aChildInd++)
      {
        WritePropsForLabel(theWS, theLabels, theLabel.FindChild(aChildInd), theIsMulti);
      }
    }
  }

  return true;
}

//=================================================================================================

bool STEPCAFControl_Writer::writeValProps(const occ::handle<XSControl_WorkSession>& theWS,
                                          const NCollection_Sequence<TDF_Label>&    theLabels,
                                          const char* const theIsMulti) const
{
  if (theLabels.IsEmpty())
  {
    return false;
  }

  // Iterate on requested shapes
  for (NCollection_Sequence<TDF_Label>::Iterator aLabelIter(theLabels); aLabelIter.More();
       aLabelIter.Next())
  {
    const TDF_Label& aLabel = aLabelIter.Value();

    WritePropsForLabel(theWS, myLabels, aLabel, theIsMulti);
  }

  return true;
}

//=================================================================================================

bool STEPCAFControl_Writer::writeLayers(const occ::handle<XSControl_WorkSession>& theWS,
                                        const NCollection_Sequence<TDF_Label>&    theLabels) const
{

  if (theLabels.IsEmpty())
  {
    return false;
  }

  // get working data
  const occ::handle<Interface_InterfaceModel>& aModel = theWS->Model();
  const occ::handle<XSControl_TransferWriter>& aTW    = theWS->TransferWriter();
  const occ::handle<Transfer_FinderProcess>&   aFP    = aTW->FinderProcess();
  NCollection_Sequence<TDF_Label>              aGlobalLayerLS;
  // Iterate on requested shapes collect Tools
  for (NCollection_Map<TDF_Label>::Iterator aLabelIter(myRootLabels); aLabelIter.More();
       aLabelIter.Next())
  {
    const TDF_Label&                aLabel = aLabelIter.Value();
    occ::handle<XCAFDoc_LayerTool>  aLTool = XCAFDoc_DocumentTool::LayerTool(aLabel);
    NCollection_Sequence<TDF_Label> aLayerLS;
    aLTool->GetLayerLabels(aLayerLS);
    aGlobalLayerLS.Append(aLayerLS);
  }
  for (NCollection_Sequence<TDF_Label>::Iterator aLayerIter(aGlobalLayerLS); aLayerIter.More();
       aLayerIter.Next())
  {
    const TDF_Label& aLayerL = aLayerIter.Value();
    // get labels of shapes in that layer
    NCollection_Sequence<TDF_Label> aShapeLabels;
    XCAFDoc_LayerTool::GetShapesOfLayer(aLayerL, aShapeLabels);
    if (aShapeLabels.IsEmpty())
    {
      continue;
    }

    // name of layer: if not set, is considered as being empty
    occ::handle<TCollection_HAsciiString> aHName = new TCollection_HAsciiString;
    GetLabelName(aLayerL, aHName);

    // Find target STEP entity for each shape and add to StepVisual_PresentationLayerAssignment
    // items.
    NCollection_Sequence<occ::handle<Standard_Transient>> aSeqRI;
    for (NCollection_Sequence<TDF_Label>::Iterator aShapeIter(aShapeLabels); aShapeIter.More();
         aShapeIter.Next())
    {
      const TDF_Label& aShLabel = aShapeIter.Value();
      if (aShLabel.IsNull())
      {
        continue;
      }

      // there is no way to assign layer to instance in STEP
      if (XCAFDoc_ShapeTool::IsAssembly(aShLabel) || XCAFDoc_ShapeTool::IsReference(aShLabel))
      {
        continue;
      }
      // check that the shape is one of (uub)labels written during current transfer
      bool anIsWritten = false;
      for (NCollection_Sequence<TDF_Label>::Iterator aLabelIter(theLabels); aLabelIter.More();
           aLabelIter.Next())
      {
        if (aShLabel.IsDescendant(aLabelIter.Value()))
        {
          anIsWritten = true;
          break;
        }
      }
      if (!anIsWritten)
      {
        continue;
      }

      // get target STEP entity
      TopoDS_Shape anOneShape = XCAFDoc_ShapeTool::GetShape(aShLabel);

      TopLoc_Location aLoc;
      int             aNb = FindEntities(aFP, anOneShape, aLoc, aSeqRI);
      if (aNb <= 0)
      {
        Message::SendTrace() << "Warning: Cannot find RI for "
                             << anOneShape.TShape()->DynamicType()->Name() << "\n";
      }
    }
    if (aSeqRI.IsEmpty())
    {
      continue;
    }

    // analyze visibility
    occ::handle<StepVisual_PresentationLayerAssignment> aStepLayerAs =
      new StepVisual_PresentationLayerAssignment;
    occ::handle<TCollection_HAsciiString> aDescr;
    occ::handle<TDataStd_UAttribute>      aUAttr;
    bool                                  isLinv = false;
    if (aLayerL.FindAttribute(XCAFDoc::InvisibleGUID(), aUAttr))
    {
      aDescr = new TCollection_HAsciiString("invisible");
      Message::SendTrace() << "\tLayer \"" << aHName->String().ToCString() << "\" is invisible"
                           << "\n";
      isLinv = true;
    }
    else
    {
      aDescr = new TCollection_HAsciiString("visible");
    }

    // create layer entity
    int                                                      anSetStyleInd = 1;
    occ::handle<NCollection_HArray1<StepVisual_LayeredItem>> aHArrayOfLItem =
      new NCollection_HArray1<StepVisual_LayeredItem>(1, aSeqRI.Length());
    for (NCollection_Sequence<occ::handle<Standard_Transient>>::Iterator aRIIter(aSeqRI);
         aRIIter.More();
         aRIIter.Next())
    {
      StepVisual_LayeredItem aLI;
      aLI.SetValue(aRIIter.Value());
      aHArrayOfLItem->SetValue(anSetStyleInd++, aLI);
    }
    aStepLayerAs->Init(aHName, aDescr, aHArrayOfLItem);
    aModel->AddWithRefs(aStepLayerAs);
    // PTV 23.01.2003 add the invisibility AFTER adding layer into the model.
    // add the invisibility for the layer
    if (isLinv)
    {
      // Invisibility Item for containing invisible layers.
      occ::handle<NCollection_HArray1<StepVisual_InvisibleItem>> aHInvsblItm =
        new NCollection_HArray1<StepVisual_InvisibleItem>(1, 1);
      StepVisual_InvisibleItem aInvIt;
      aInvIt.SetValue(aStepLayerAs);
      aHInvsblItm->SetValue(1, aInvIt);

      occ::handle<StepVisual_Invisibility> aInvsblt = new StepVisual_Invisibility();
      aInvsblt->Init(aHInvsblItm);
      aModel->AddWithRefs(aInvsblt);
    }
  }
  return true;
}

//=================================================================================================

static bool getSHUOstyle(const TDF_Label& theSHUOlab, XCAFPrs_Style& theSHUOstyle)
{
  Quantity_Color aColor;
  if (!XCAFDoc_ColorTool::IsVisible(theSHUOlab))
  {
    theSHUOstyle.SetVisibility(false);
  }
  else
  {
    if (XCAFDoc_ColorTool::GetColor(theSHUOlab, XCAFDoc_ColorGen, aColor))
    {
      theSHUOstyle.SetColorCurv(aColor);
      theSHUOstyle.SetColorSurf(aColor);
    }
    if (XCAFDoc_ColorTool::GetColor(theSHUOlab, XCAFDoc_ColorSurf, aColor))
    {
      theSHUOstyle.SetColorSurf(aColor);
    }
    if (XCAFDoc_ColorTool::GetColor(theSHUOlab, XCAFDoc_ColorCurv, aColor))
    {
      theSHUOstyle.SetColorCurv(aColor);
    }
    if (!theSHUOstyle.IsSetColorSurf())
    {
      occ::handle<XCAFDoc_VisMaterial> aVisMat =
        XCAFDoc_VisMaterialTool::GetShapeMaterial(theSHUOlab);
      if (!aVisMat.IsNull() && !aVisMat->IsEmpty())
      {
        // only color can be stored in STEP
        theSHUOstyle.SetColorSurf(aVisMat->BaseColor());
      }
    }
  }
  if (!theSHUOstyle.IsSetColorCurv() && !theSHUOstyle.IsSetColorSurf() && theSHUOstyle.IsVisible())
  {
    return false;
  }
  return true;
}

//=================================================================================================

static bool getProDefinitionOfNAUO(const occ::handle<XSControl_WorkSession>&          theWS,
                                   const TopoDS_Shape&                                theShape,
                                   occ::handle<StepBasic_ProductDefinition>&          thePD,
                                   occ::handle<StepRepr_NextAssemblyUsageOccurrence>& theNAUO,
                                   bool                                               theIsRelating)
{
  if (theShape.IsNull())
  {
    return false;
  }
  // get CDSR
  const occ::handle<XSControl_TransferWriter>&               aTW = theWS->TransferWriter();
  const occ::handle<Transfer_FinderProcess>&                 aFP = aTW->FinderProcess();
  occ::handle<StepShape_ContextDependentShapeRepresentation> aCDSR;
  occ::handle<TransferBRep_ShapeMapper> aMapper = TransferBRep::ShapeMapper(aFP, theShape);
  if (!aFP->FindTypedTransient(aMapper,
                               STANDARD_TYPE(StepShape_ContextDependentShapeRepresentation),
                               aCDSR))
  {
    return false;
  }
  // get PDS of NAUO
  occ::handle<StepRepr_ProductDefinitionShape> aPDS = aCDSR->RepresentedProductRelation();
  if (aPDS.IsNull())
  {
    return false;
  }
  // get the NAUO entity
  Interface_Graph aGraph = theWS->HGraph()->Graph();
  for (Interface_EntityIterator aShareds = aGraph.Shareds(aPDS); aShareds.More(); aShareds.Next())
  {
    theNAUO = occ::down_cast<StepRepr_NextAssemblyUsageOccurrence>(aShareds.Value());
    if (!theNAUO.IsNull())
    {
      break;
    }
  }
  if (theNAUO.IsNull())
  {
    return false;
  }
  // get Relatinf or Related product definition
  if (!theIsRelating)
  {
    thePD = theNAUO->RelatedProductDefinition();
  }
  else
  {
    thePD = theNAUO->RelatingProductDefinition();
  }
  if (thePD.IsNull())
  {
    return false;
  }
  return true;
}

//=================================================================================================

static bool writeSHUO(const occ::handle<XCAFDoc_GraphNode>&                 theSHUO,
                      const occ::handle<XSControl_WorkSession>&             theWS,
                      occ::handle<StepRepr_SpecifiedHigherUsageOccurrence>& theTopSHUO,
                      TopoDS_Shape&                                         theNAUOShape,
                      occ::handle<StepBasic_ProductDefinition>&             theRelatingPD,
                      bool&                                                 theIsDeepest)
{
  // set the ProductDefinitionRelationship descriptin information as empty strings.
  occ::handle<TCollection_HAsciiString> anEmptyString = new TCollection_HAsciiString("");

  NCollection_Sequence<TDF_Label> aNextUsageLabs;
  XCAFDoc_ShapeTool::GetSHUONextUsage(theSHUO->Label(), aNextUsageLabs);
  occ::handle<XCAFDoc_GraphNode> aNuSHUO;
  if (theTopSHUO.IsNull())
  {
    // the top SHUO
    if (aNextUsageLabs.Length() < 1)
    {
      return false;
    }
    XCAFDoc_ShapeTool::GetSHUO(aNextUsageLabs.Value(1), aNuSHUO);
    if (aNuSHUO.IsNull())
    {
      return false;
    }
    // get relating product definition
    TopoDS_Shape aTopCompShape = XCAFDoc_ShapeTool::GetShape(theSHUO->Label().Father());
    occ::handle<StepRepr_NextAssemblyUsageOccurrence> aRelatingNAUO;
    if (!getProDefinitionOfNAUO(theWS, aTopCompShape, theRelatingPD, aRelatingNAUO, true))
    {
      return false;
    }
    // get related product definition
    TopoDS_Shape aNUShape = XCAFDoc_ShapeTool::GetShape(aNuSHUO->Label().Father());
    occ::handle<StepBasic_ProductDefinition>          aRelatedPD;
    occ::handle<StepRepr_NextAssemblyUsageOccurrence> aRelatedNAUO;
    if (!getProDefinitionOfNAUO(theWS, aNUShape, aRelatedPD, aRelatedNAUO, false))
    {
      return false;
    }

    theTopSHUO = new StepRepr_SpecifiedHigherUsageOccurrence;
    // create deepest shuo EmptyString
    theTopSHUO->Init(/*id*/ anEmptyString,
                     /*name*/ anEmptyString,
                     /*no description*/ false,
                     /*description*/ anEmptyString,
                     theRelatingPD,
                     aRelatedPD,
                     /*no ACURefDesignator*/ false,
                     /*ACURefDesignator*/ anEmptyString,
                     /*upper_usage*/ aRelatingNAUO,
                     /*next_usage*/ aRelatedNAUO);
    // write the other SHUO.
    if (!writeSHUO(aNuSHUO, theWS, theTopSHUO, theNAUOShape, theRelatingPD, theIsDeepest))
    {
      theTopSHUO.Nullify();
      return false;
    }

    return true;
  }
  //   occ::handle<XCAFDoc_GraphNode> NuSHUO;
  if (aNextUsageLabs.Length() > 0)
  {
    // store SHUO recursive
    if (aNextUsageLabs.Length() > 1)
    {
      Message::SendTrace() << "Warning: store only one next_usage of current SHUO" << "\n";
    }
    XCAFDoc_ShapeTool::GetSHUO(aNextUsageLabs.Value(1), aNuSHUO);
    occ::handle<StepRepr_SpecifiedHigherUsageOccurrence> aNUEntSHUO =
      new StepRepr_SpecifiedHigherUsageOccurrence;
    if (!writeSHUO(aNuSHUO, theWS, aNUEntSHUO, theNAUOShape, theRelatingPD, theIsDeepest))
    {
      return false;
    }

    // store the deepest SHUO to the dociment
    TopoDS_Shape aNUSh, aUUSh;
    aNUSh = XCAFDoc_ShapeTool::GetShape(aNuSHUO->Label().Father());
    aUUSh = XCAFDoc_ShapeTool::GetShape(theSHUO->Label().Father());
    // get relating PD with upper_usage and related PD with next_usage
    occ::handle<StepBasic_ProductDefinition>
      nullPD; // no need to use,case have shared <theRelatingPD>
    occ::handle<StepBasic_ProductDefinition>          aRelatedPD;
    occ::handle<StepRepr_NextAssemblyUsageOccurrence> UUNAUO, NUNAUO;
    if (!getProDefinitionOfNAUO(theWS, aUUSh, nullPD, UUNAUO, true)
        || !getProDefinitionOfNAUO(theWS, aNUSh, aRelatedPD, NUNAUO, false))
    {
      Message::SendTrace() << "Warning: cannot get related or relating PD" << "\n";
      return false;
    }
    aNUEntSHUO->Init(/*id*/ anEmptyString,
                     /*name*/ anEmptyString,
                     /*no description*/ false,
                     /*description*/ anEmptyString,
                     theRelatingPD,
                     aRelatedPD,
                     /*no ACURefDesignator*/ false,
                     /*ACURefDesignator*/ anEmptyString,
                     /*upper_usage*/ theTopSHUO,
                     /*next_usage*/ NUNAUO);
    if (theIsDeepest)
    {
      theIsDeepest = false;
    }
    theWS->Model()->AddWithRefs(aNUEntSHUO);
    return true;
  } // end of recurse storing

  // get shape
  TDF_Label aShapeL = theSHUO->Label().Father();
  theNAUOShape      = XCAFDoc_ShapeTool::GetShape(aShapeL);
  // return to the deepest level from SHUO shape level
  // it is because SHUO is attribute on deep level and shape level.
  theIsDeepest = true;
  return true;
}

//=================================================================================================

static bool createSHUOStyledItem(
  const XCAFPrs_Style&                                theStyle,
  const occ::handle<StepRepr_ProductDefinitionShape>& thePDS,
  const occ::handle<XSControl_WorkSession>&           theWS,
  const TopoDS_Shape&                                 theShape,
  const occ::handle<XCAFDoc_ShapeTool>&               theSTool,
  NCollection_DataMap<TopoDS_Shape, occ::handle<Standard_Transient>, TopTools_ShapeMapHasher>&
    theMapCompMDGPR)
{
  // create styled item for the indicated SHUO and store to the model
  STEPConstruct_Styles aStyles(theWS);
  // translate colors to STEP
  occ::handle<StepVisual_Colour> aSurfColor, aCurvColor;
  double                         aRenderTransp = 0.0;
  if (theStyle.IsSetColorSurf())
  {
    Quantity_ColorRGBA aSurfCol = theStyle.GetColorSurfRGBA();
    aRenderTransp               = 1.0 - aSurfCol.Alpha();
    aSurfColor                  = STEPConstruct_Styles::EncodeColor(aSurfCol.GetRGB());
  }
  if (theStyle.IsSetColorCurv())
  {
    aCurvColor = STEPConstruct_Styles::EncodeColor(theStyle.GetColorCurv());
  }
  bool                                     isComponent = true; // cause need to get PSBC
  occ::handle<StepRepr_RepresentationItem> anItem;
  // set default color for invisible SHUO.
  bool isSetDefaultColor = false;
  if (aSurfColor.IsNull() && aCurvColor.IsNull() && !theStyle.IsVisible())
  {
    aSurfColor        = STEPConstruct_Styles::EncodeColor(Quantity_Color(Quantity_NOC_WHITE));
    isSetDefaultColor = true;
  }
  STEPConstruct_RenderingProperties aRenderProps;
  if (aRenderTransp > 0.0)
  {
    aRenderProps.Init(theStyle.GetColorSurfRGBA());
  }
  occ::handle<StepVisual_PresentationStyleAssignment> aPSA =
    aStyles.MakeColorPSA(anItem, aSurfColor, aCurvColor, aRenderProps, isComponent);
  occ::handle<StepVisual_StyledItem> anOverride; // null styled item

  // find the repr item of the shape
  const occ::handle<XSControl_TransferWriter>& aTW = theWS->TransferWriter();
  const occ::handle<Transfer_FinderProcess>&   aFP = aTW->FinderProcess();
  occ::handle<StepData_StepModel> aStepModel    = occ::down_cast<StepData_StepModel>(aFP->Model());
  occ::handle<TransferBRep_ShapeMapper> aMapper = TransferBRep::ShapeMapper(aFP, theShape);
  occ::handle<StepShape_ContextDependentShapeRepresentation> aCDSR;
  aFP->FindTypedTransient(aMapper,
                          STANDARD_TYPE(StepShape_ContextDependentShapeRepresentation),
                          aCDSR);
  if (aCDSR.IsNull())
  {
    return false;
  }
  // find context
  occ::handle<StepRepr_RepresentationContext> aContext = aStyles.FindContext(theShape);
  TopoDS_Shape                                aTopSh   = theShape;
  if (aContext.IsNull())
  {
    TDF_Label aTopShL = theSTool->FindShape(theShape, false);
    if (aTopShL.IsNull())
    {
      return false;
    }
    aTopSh   = XCAFDoc_ShapeTool::GetShape(aTopShL);
    aContext = aStyles.FindContext(aTopSh);
  }
  if (aContext.IsNull())
  {
    return false;
  }
  // get representation item of the shape
  TopLoc_Location                                       aLocation;
  NCollection_Sequence<occ::handle<Standard_Transient>> aSeqRI;
  FindEntities(aFP, theShape, aLocation, aSeqRI);
  if (aSeqRI.Length() <= 0)
  {
    Message::SendTrace() << "Warning: Cannot find RI for "
                         << theShape.TShape()->DynamicType()->Name() << "\n";
  }
  anItem = occ::down_cast<StepRepr_RepresentationItem>(aSeqRI(1));
  // get overridden styled item
  getStyledItem(theShape, theSTool, aStyles, anOverride, theMapCompMDGPR);

  // get STEP STYLED ITEM
  occ::handle<StepVisual_StyledItem> aSTEPstyle = aStyles.AddStyle(anItem, aPSA, anOverride);
  // create SR, SDR and all necessary references between them and ST, PDS, PSBC, GRC
  aStyles.CreateNAUOSRD(aContext, aCDSR, thePDS);

  // add step styled item of SHUO to the model
  // do it by additing styled item to the MDGPR
  if (!aTopSh.IsNull() && !theMapCompMDGPR.IsBound(aTopSh))
  {
    // create MDGPR and record it in model
    Message::SendTrace() << "Warning: " << __FILE__ << ": Create new MDGPR for SHUO instance"
                         << "\n";
    occ::handle<StepVisual_MechanicalDesignGeometricPresentationRepresentation> aMDGPR;
    aStyles.CreateMDGPR(aContext, aMDGPR, aStepModel);
    if (!aMDGPR.IsNull())
    {
      theMapCompMDGPR.Bind(aTopSh, aMDGPR);
    }
  }
  else if (!aTopSh.IsNull() && theMapCompMDGPR.IsBound(aTopSh))
  {
    // get MDGPR of the top-level shape
    occ::handle<StepVisual_PresentationRepresentation> aMDGPR =
      occ::down_cast<StepVisual_PresentationRepresentation>(theMapCompMDGPR.Find(aTopSh));
    // get old styled items to not lose it
    occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> anOldItems =
      aMDGPR->Items();
    int anOldLengthlen = 0;
    if (!anOldItems.IsNull())
    {
      anOldLengthlen = anOldItems->Length();
    }
    // create new array of styled items by an olds and new one
    occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> aNewItems =
      new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, anOldLengthlen + 1);
    int aSetStyleInd = 1;
    for (NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>::Iterator aStyleIter(
           anOldItems->Array1());
         aStyleIter.More();
         aStyleIter.Next())
    {
      aNewItems->SetValue(aSetStyleInd++, aStyleIter.Value());
    }
    aNewItems->SetValue(aSetStyleInd++, aSTEPstyle);
    // init MDGPR be new array of styled items
    if (aNewItems->Length() > 0)
    {
      aMDGPR->SetItems(aNewItems);
    }
  }
  else
  {
    theWS->Model()->AddWithRefs(aSTEPstyle); // add as root to the model, but it is not good
    Message::SendTrace() << "Warning: " << __FILE__
                         << ": adds styled item of SHUO as root, cause cannot find MDGPR" << "\n";
  }
  // create invisibility item for the styled item
  if (!theStyle.IsVisible())
  {
    if (isSetDefaultColor)
    {
      // try to set default color from top-level shape

      setDefaultInstanceColor(anOverride, aPSA);
    }
    // create invisibility item and refer for styledItem
    occ::handle<StepVisual_Invisibility> aInvsblt = new StepVisual_Invisibility();
    occ::handle<NCollection_HArray1<StepVisual_InvisibleItem>> aHInvsblItm =
      new NCollection_HArray1<StepVisual_InvisibleItem>(1, 1);
    // put all style item into the harray
    StepVisual_InvisibleItem anInvItem;
    anInvItem.SetValue(aSTEPstyle);
    aHInvsblItm->SetValue(1, anInvItem);
    aInvsblt->Init(aHInvsblItm);
    theWS->Model()->AddWithRefs(aInvsblt);
  }

  return true;
}

//=================================================================================================

bool STEPCAFControl_Writer::writeSHUOs(const occ::handle<XSControl_WorkSession>& theWS,
                                       const NCollection_Sequence<TDF_Label>&    theLabels)
{
  if (theLabels.IsEmpty())
  {
    return false;
  }

  // map of transferred SHUO
  NCollection_Map<occ::handle<Standard_Transient>> aMapOfMainSHUO;
  // Iterate on requested shapes
  for (NCollection_Sequence<TDF_Label>::Iterator aLabelIter(theLabels); aLabelIter.More();
       aLabelIter.Next())
  {
    const TDF_Label& aLabel = aLabelIter.Value();
    if (!myLabels.IsBound(aLabel))
    {
      continue; // not recorded as translated, skip
    }
    if (!XCAFDoc_ShapeTool::IsAssembly(aLabel))
    {
      continue;
    }
    NCollection_Sequence<TDF_Label> aLabelSeq;
    XCAFDoc_ShapeTool::GetComponents(aLabel, aLabelSeq);
    // iterates on components of assembly
    for (NCollection_Sequence<TDF_Label>::Iterator aCompIter(aLabelSeq); aCompIter.More();
         aCompIter.Next())
    {
      const TDF_Label&                                 aCompL = aCompIter.Value();
      NCollection_Sequence<occ::handle<TDF_Attribute>> anAttrSeq;
      XCAFDoc_ShapeTool::GetAllComponentSHUO(aCompL, anAttrSeq);
      // work with SHUO
      for (NCollection_Sequence<occ::handle<TDF_Attribute>>::Iterator anAttrIter(anAttrSeq);
           anAttrIter.More();
           anAttrIter.Next())
      {
        occ::handle<XCAFDoc_GraphNode> aSHUO =
          occ::down_cast<XCAFDoc_GraphNode>(anAttrIter.Value());
        // take label of SHUO
        TDF_Label                       aSHUOlab = aSHUO->Label();
        NCollection_Sequence<TDF_Label> aUpLabels;
        // check is it SHUO of upper_usage
        XCAFDoc_ShapeTool::GetSHUOUpperUsage(aSHUOlab, aUpLabels);
        if (aUpLabels.Length() > 0)
        {
          continue; // transfer only main SHUO
        }
        if (aMapOfMainSHUO.Contains(aSHUO))
        {
          continue; // do not try to transfer SHUO twice
        }
        aMapOfMainSHUO.Add(aSHUO);
        // check if it is styled SHUO
        XCAFPrs_Style aSHUOstyle;
        if (!getSHUOstyle(aSHUOlab, aSHUOstyle))
        {
          Message::SendTrace() << "Warning: " << __FILE__
                               << ": do not store SHUO without any style to the STEP model" << "\n";
          continue;
        }
        // write SHUO to the model amd then add structure type.
        TopoDS_Shape aNAUOShape; // shape of the deepest NAUO in the SHUO structure
        bool         isDeepest = false;
        occ::handle<StepRepr_SpecifiedHigherUsageOccurrence> anEntOfSHUO;
        occ::handle<StepBasic_ProductDefinition>             aRelatingPD;
        // create the top SHUO and all other.
        writeSHUO(aSHUO, theWS, anEntOfSHUO, aNAUOShape, aRelatingPD, isDeepest);
        if (anEntOfSHUO.IsNull() || aNAUOShape.IsNull())
        {
          Message::SendTrace() << "Warning: " << __FILE__ << ": Cannot store SHUO" << "\n";
          continue;
        }
        // create new Product Definition Shape for TOP SHUO
        Message::SendTrace() << "Info: " << __FILE__ << ": Create NEW PDS for current SHUO "
                             << "\n";
        occ::handle<StepRepr_ProductDefinitionShape> aPDS = new StepRepr_ProductDefinitionShape;
        occ::handle<TCollection_HAsciiString> aPDSname    = new TCollection_HAsciiString("SHUO");
        occ::handle<TCollection_HAsciiString> aDescrStr   = new TCollection_HAsciiString("");
        StepRepr_CharacterizedDefinition      aCharDef;
        aCharDef.SetValue(anEntOfSHUO);
        aPDS->Init(aPDSname, false, aDescrStr, aCharDef);

        // create styled item for SHUO and add to the model
        createSHUOStyledItem(aSHUOstyle,
                             aPDS,
                             theWS,
                             aNAUOShape,
                             XCAFDoc_DocumentTool::ShapeTool(aLabel),
                             myMapCompMDGPR);

      } // end work with SHUO
    } // end of an assembly components
  } // end of iterates on indicated labels
  return true;
}

//=======================================================================
// function : FindPDSforDGT
// purpose  : auxiliary: find PDS for AdvancedFace or EdgeCurve for creation
//                     needed ShapeAspect in D&GT structure
//=======================================================================
static bool FindPDSforDGT(const Interface_Graph&                        theGraph,
                          const occ::handle<Standard_Transient>&        theEnt,
                          occ::handle<StepRepr_ProductDefinitionShape>& thePDS,
                          occ::handle<StepRepr_RepresentationContext>&  theRC,
                          occ::handle<StepShape_AdvancedFace>&          theAF,
                          occ::handle<StepShape_EdgeCurve>&             theEC)
{
  if (theEnt.IsNull())
  {
    return false;
  }
  if (!theEnt->IsKind(STANDARD_TYPE(StepShape_EdgeCurve))
      && !theEnt->IsKind(STANDARD_TYPE(StepShape_AdvancedFace)))
  {
    return false;
  }

  theAF = occ::down_cast<StepShape_AdvancedFace>(theEnt);
  if (theAF.IsNull())
  {
    theEC = occ::down_cast<StepShape_EdgeCurve>(theEnt);
    for (Interface_EntityIterator aSharingsIter = theGraph.Sharings(theEC);
         aSharingsIter.More() && theAF.IsNull();
         aSharingsIter.Next())
    {
      occ::handle<StepShape_OrientedEdge> anOE =
        occ::down_cast<StepShape_OrientedEdge>(aSharingsIter.Value());
      if (anOE.IsNull())
      {
        continue;
      }
      for (Interface_EntityIterator aSubs1 = theGraph.Sharings(anOE);
           aSubs1.More() && theAF.IsNull();
           aSubs1.Next())
      {
        occ::handle<StepShape_EdgeLoop> aEL = occ::down_cast<StepShape_EdgeLoop>(aSubs1.Value());
        if (aEL.IsNull())
        {
          continue;
        }
        for (Interface_EntityIterator aSubs2 = theGraph.Sharings(aEL);
             aSubs2.More() && theAF.IsNull();
             aSubs2.Next())
        {
          occ::handle<StepShape_FaceBound> aFB =
            occ::down_cast<StepShape_FaceBound>(aSubs2.Value());
          if (aFB.IsNull())
          {
            continue;
          }
          for (Interface_EntityIterator aSubs3 = theGraph.Sharings(aFB);
               aSubs3.More() && theAF.IsNull();
               aSubs3.Next())
          {
            theAF = occ::down_cast<StepShape_AdvancedFace>(aSubs3.Value());
          }
        }
      }
    }
  }
  if (theAF.IsNull())
  {
    return false;
  }

  for (Interface_EntityIterator aSharingsIter = theGraph.Sharings(theAF);
       aSharingsIter.More() && thePDS.IsNull();
       aSharingsIter.Next())
  {
    occ::handle<StepShape_ConnectedFaceSet> aCFS =
      occ::down_cast<StepShape_ConnectedFaceSet>(aSharingsIter.Value());
    if (aCFS.IsNull())
    {
      continue;
    }
    for (Interface_EntityIterator aSubs1 = theGraph.Sharings(aCFS);
         aSubs1.More() && thePDS.IsNull();
         aSubs1.Next())
    {
      occ::handle<StepRepr_RepresentationItem> aRI =
        occ::down_cast<StepRepr_RepresentationItem>(aSubs1.Value());
      if (aRI.IsNull())
      {
        continue;
      }
      for (Interface_EntityIterator aSubs2 = theGraph.Sharings(aRI);
           aSubs2.More() && thePDS.IsNull();
           aSubs2.Next())
      {
        occ::handle<StepShape_ShapeRepresentation> aSR =
          occ::down_cast<StepShape_ShapeRepresentation>(aSubs2.Value());
        if (aSR.IsNull())
        {
          continue;
        }
        theRC = aSR->ContextOfItems();
        for (Interface_EntityIterator aSubs3 = theGraph.Sharings(aSR);
             aSubs3.More() && thePDS.IsNull();
             aSubs3.Next())
        {
          occ::handle<StepShape_ShapeDefinitionRepresentation> aSDR =
            occ::down_cast<StepShape_ShapeDefinitionRepresentation>(aSubs3.Value());
          if (aSDR.IsNull())
          {
            continue;
          }
          occ::handle<StepRepr_PropertyDefinition> aPropD = aSDR->Definition().PropertyDefinition();
          if (aPropD.IsNull())
          {
            continue;
          }
          thePDS = occ::down_cast<StepRepr_ProductDefinitionShape>(aPropD);
        }
      }
    }
  }
  return true;
}

//=======================================================================
// function : FindPDS
// purpose  : auxiliary: find Product_definition_shape entity for given entity
//=======================================================================
static occ::handle<StepRepr_ProductDefinitionShape> FindPDS(
  const Interface_Graph&                       theGraph,
  const occ::handle<Standard_Transient>&       theEnt,
  occ::handle<StepRepr_RepresentationContext>& theRC)
{
  if (theEnt.IsNull())
  {
    return nullptr;
  }
  occ::handle<StepRepr_ProductDefinitionShape> aPDS;

  // try to find shape_representation in sharings
  for (Interface_EntityIterator anIter = theGraph.Sharings(theEnt); anIter.More() && aPDS.IsNull();
       anIter.Next())
  {
    occ::handle<StepShape_ShapeRepresentation> aSR =
      occ::down_cast<StepShape_ShapeRepresentation>(anIter.Value());
    if (aSR.IsNull())
    {
      continue;
    }
    theRC                           = aSR->ContextOfItems();
    Interface_EntityIterator aSDRIt = theGraph.Sharings(aSR);
    for (aSDRIt.Start(); aSDRIt.More() && aPDS.IsNull(); aSDRIt.Next())
    {
      occ::handle<StepShape_ShapeDefinitionRepresentation> aSDR =
        occ::down_cast<StepShape_ShapeDefinitionRepresentation>(aSDRIt.Value());
      if (aSDR.IsNull())
      {
        continue;
      }
      occ::handle<StepRepr_PropertyDefinition> aPropD = aSDR->Definition().PropertyDefinition();
      if (aPropD.IsNull())
      {
        continue;
      }
      aPDS = occ::down_cast<StepRepr_ProductDefinitionShape>(aPropD);
    }
  }
  if (!aPDS.IsNull())
  {
    return aPDS;
  }

  for (Interface_EntityIterator anIter = theGraph.Sharings(theEnt); anIter.More(); anIter.Next())
  {
    if (anIter.Value()->IsKind(STANDARD_TYPE(StepShape_TopologicalRepresentationItem))
        || anIter.Value()->IsKind(STANDARD_TYPE(StepGeom_GeometricRepresentationItem)))
    {
      aPDS = FindPDS(theGraph, anIter.Value(), theRC);
      if (!aPDS.IsNull())
      {
        return aPDS;
      }
    }
  }
  return aPDS;
}

//=================================================================================================

static StepBasic_Unit GetUnit(const occ::handle<StepRepr_RepresentationContext>& theRC,
                              const bool                                         theIsAngle = false)
{
  StepBasic_Unit                   aUnit;
  occ::handle<StepBasic_NamedUnit> aCurrentUnit;
  if (theIsAngle)
  {
    occ::handle<StepGeom_GeometricRepresentationContextAndGlobalUnitAssignedContext> aCtx =
      occ::down_cast<StepGeom_GeometricRepresentationContextAndGlobalUnitAssignedContext>(theRC);
    if (!aCtx.IsNull())
    {
      for (NCollection_HArray1<occ::handle<StepBasic_NamedUnit>>::Iterator aUnitIter(
             aCtx->Units()->Array1());
           aUnitIter.More();
           aUnitIter.Next())
      {
        const occ::handle<StepBasic_NamedUnit>& aCurUnitValue = aUnitIter.Value();
        if (aCurUnitValue->IsKind(STANDARD_TYPE(StepBasic_SiUnitAndPlaneAngleUnit))
            || aCurUnitValue->IsKind(STANDARD_TYPE(StepBasic_ConversionBasedUnitAndPlaneAngleUnit)))
        {
          aCurrentUnit = aCurUnitValue;
          break;
        }
      }
    }
    if (aCurrentUnit.IsNull())
    {
      occ::handle<StepGeom_GeomRepContextAndGlobUnitAssCtxAndGlobUncertaintyAssCtx> aCtx1 =
        occ::down_cast<StepGeom_GeomRepContextAndGlobUnitAssCtxAndGlobUncertaintyAssCtx>(theRC);
      if (!aCtx1.IsNull())
      {
        for (NCollection_HArray1<occ::handle<StepBasic_NamedUnit>>::Iterator aUnitIter(
               aCtx1->Units()->Array1());
             aUnitIter.More();
             aUnitIter.Next())
        {
          const occ::handle<StepBasic_NamedUnit>& aCurUnitValue = aUnitIter.Value();
          if (aCurUnitValue->IsKind(STANDARD_TYPE(StepBasic_SiUnitAndPlaneAngleUnit))
              || aCurUnitValue->IsKind(
                STANDARD_TYPE(StepBasic_ConversionBasedUnitAndPlaneAngleUnit)))
          {
            aCurrentUnit = aCurUnitValue;
            break;
          }
        }
      }
    }
    if (aCurrentUnit.IsNull())
    {
      aCurrentUnit = new StepBasic_SiUnitAndPlaneAngleUnit;
    }
    // XCAF holds angles in degrees, the context's unit is the radian: the values go in a
    // DEGREE converted from it, CONVERSION_BASED_UNIT('DEGREE', 0.0174532925199433 rad), as
    // the context of a file in degrees has it
    if (aCurrentUnit->IsKind(STANDARD_TYPE(StepBasic_SiUnitAndPlaneAngleUnit)))
    {
      StepBasic_Unit aRadian;
      aRadian.SetValue(aCurrentUnit);
      occ::handle<StepBasic_MeasureValueMember> aFactor = new StepBasic_MeasureValueMember();
      aFactor->SetName("PLANE_ANGLE_MEASURE");
      aFactor->SetReal(M_PI / 180.);
      occ::handle<StepBasic_PlaneAngleMeasureWithUnit> aConversion =
        new StepBasic_PlaneAngleMeasureWithUnit();
      aConversion->Init(aFactor, aRadian);
      occ::handle<StepBasic_DimensionalExponents> anExponents =
        new StepBasic_DimensionalExponents();
      anExponents->Init(0., 0., 0., 0., 0., 0., 0.);
      occ::handle<StepBasic_ConversionBasedUnitAndPlaneAngleUnit> aDegree =
        new StepBasic_ConversionBasedUnitAndPlaneAngleUnit();
      aDegree->Init(anExponents, new TCollection_HAsciiString("DEGREE"), aConversion);
      aCurrentUnit = aDegree;
    }
  }
  else
  {
    occ::handle<StepGeom_GeometricRepresentationContextAndGlobalUnitAssignedContext> aCtx =
      occ::down_cast<StepGeom_GeometricRepresentationContextAndGlobalUnitAssignedContext>(theRC);
    if (!aCtx.IsNull())
    {
      for (NCollection_HArray1<occ::handle<StepBasic_NamedUnit>>::Iterator aUnitIter(
             aCtx->Units()->Array1());
           aUnitIter.More();
           aUnitIter.Next())
      {
        const occ::handle<StepBasic_NamedUnit>& aCurUnitValue = aUnitIter.Value();
        if (aCurUnitValue->IsKind(STANDARD_TYPE(StepBasic_SiUnitAndLengthUnit))
            || aCurUnitValue->IsKind(STANDARD_TYPE(StepBasic_ConversionBasedUnitAndLengthUnit)))
        {
          aCurrentUnit = aCurUnitValue;
          break;
        }
      }
    }
    if (aCurrentUnit.IsNull())
    {
      occ::handle<StepGeom_GeomRepContextAndGlobUnitAssCtxAndGlobUncertaintyAssCtx> aCtx1 =
        occ::down_cast<StepGeom_GeomRepContextAndGlobUnitAssCtxAndGlobUncertaintyAssCtx>(theRC);
      if (!aCtx1.IsNull())
      {
        for (NCollection_HArray1<occ::handle<StepBasic_NamedUnit>>::Iterator aUnitIter(
               aCtx1->Units()->Array1());
             aUnitIter.More();
             aUnitIter.Next())
        {
          const occ::handle<StepBasic_NamedUnit>& aCurUnitValue = aUnitIter.Value();
          if (aCurUnitValue->IsKind(STANDARD_TYPE(StepBasic_SiUnitAndLengthUnit))
              || aCurUnitValue->IsKind(STANDARD_TYPE(StepBasic_ConversionBasedUnitAndLengthUnit)))
          {
            aCurrentUnit = aCurUnitValue;
            break;
          }
        }
      }
    }
    if (aCurrentUnit.IsNull())
    {
      aCurrentUnit = new StepBasic_SiUnitAndLengthUnit;
    }
  }

  aUnit.SetValue(aCurrentUnit);
  return aUnit;
}

//=================================================================================================

static occ::handle<StepRepr_ReprItemAndMeasureWithUnit> CreateDimValue(
  const double                                              theValue,
  const StepBasic_Unit&                                     theUnit,
  const occ::handle<TCollection_HAsciiString>&              theName,
  const char* const                                         theMeasureName,
  const bool                                                theIsAngle,
  const bool                                                theIsQualified = false,
  const occ::handle<StepShape_QualifiedRepresentationItem>& theQRI         = nullptr)
{
  occ::handle<StepRepr_RepresentationItem> aReprItem = new StepRepr_RepresentationItem();
  aReprItem->Init(new TCollection_HAsciiString(theName));
  occ::handle<StepBasic_MeasureWithUnit>    aMWU         = new StepBasic_MeasureWithUnit();
  occ::handle<StepBasic_MeasureValueMember> aValueMember = new StepBasic_MeasureValueMember();
  aValueMember->SetName(theMeasureName);
  aValueMember->SetReal(theValue);
  aMWU->Init(aValueMember, theUnit);
  if (theIsQualified)
  {
    if (theIsAngle)
    {
      // Angle & with qualifiers
      occ::handle<StepRepr_ReprItemAndPlaneAngleMeasureWithUnitAndQRI> anItem =
        new StepRepr_ReprItemAndPlaneAngleMeasureWithUnitAndQRI();
      anItem->Init(aMWU, aReprItem, theQRI);
      return anItem;
    }
    else
    {
      // Length & with qualifiers
      occ::handle<StepRepr_ReprItemAndLengthMeasureWithUnitAndQRI> anItem =
        new StepRepr_ReprItemAndLengthMeasureWithUnitAndQRI();
      anItem->Init(aMWU, aReprItem, theQRI);
      return anItem;
    }
  }
  else
  {
    if (theIsAngle)
    {
      // Angle & without qualifiers
      occ::handle<StepRepr_ReprItemAndPlaneAngleMeasureWithUnit> anItem =
        new StepRepr_ReprItemAndPlaneAngleMeasureWithUnit();
      anItem->Init(aMWU, aReprItem);
      return anItem;
    }
    else
    {
      // Length & without qualifiers
      occ::handle<StepRepr_ReprItemAndLengthMeasureWithUnit> anItem =
        new StepRepr_ReprItemAndLengthMeasureWithUnit();
      anItem->Init(aMWU, aReprItem);
      return anItem;
    }
  }
}

//=======================================================================
// function : writeShapeAspect
// purpose  : auxiliary (write Shape_Aspect entity for given shape)
//=======================================================================
occ::handle<StepRepr_ShapeAspect> STEPCAFControl_Writer::writeShapeAspect(
  const occ::handle<XSControl_WorkSession>&          theWS,
  const TDF_Label                                    theLabel,
  const TopoDS_Shape&                                theShape,
  occ::handle<StepRepr_RepresentationContext>&       theRC,
  occ::handle<StepAP242_GeometricItemSpecificUsage>& theGISU)
{
  // Get working data
  const occ::handle<Interface_InterfaceModel>& aModel  = theWS->Model();
  const occ::handle<XSControl_TransferWriter>& aTW     = theWS->TransferWriter();
  const occ::handle<Transfer_FinderProcess>&   aFP     = aTW->FinderProcess();
  const occ::handle<Interface_HGraph>          aHGraph = theWS->HGraph();
  if (aHGraph.IsNull())
  {
    return nullptr;
  }
  Interface_Graph aGraph = aHGraph->Graph();

  TopLoc_Location                                       aLoc;
  NCollection_Sequence<occ::handle<Standard_Transient>> aSeqRI;
  FindEntities(aFP, theShape, aLoc, aSeqRI);
  if (aSeqRI.Length() <= 0)
  {
    Message::SendTrace() << "Warning: Cannot find RI for "
                         << theShape.TShape()->DynamicType()->Name() << "\n";
    return nullptr;
  }

  occ::handle<StepRepr_ProductDefinitionShape> aPDS;
  occ::handle<StepRepr_RepresentationContext>  aRC;
  occ::handle<Standard_Transient>              anEnt = aSeqRI.Value(1);
  // Supplemental geometry's: of its part, in its constructive geometry representation
  occ::handle<StepRepr_Representation> aUsedRep;
  if (const GDTPresentationModel::Supplemental* aSupplemental =
        gdtState(myGDTPresentationDM).Supplementals.Seek(anEnt))
  {
    aPDS     = aSupplemental->PDS;
    aUsedRep = aSupplemental->Representation;
    aRC      = aUsedRep->ContextOfItems();
  }
  else
  {
    aPDS = FindPDS(aGraph, anEnt, aRC);
  }
  if (aPDS.IsNull())
  {
    return nullptr;
  }

  theRC = aRC;
  // Shape_Aspect
  occ::handle<TCollection_HAsciiString> aName = new TCollection_HAsciiString();
  occ::handle<TDataStd_Name>            aNameAttr;
  if (theLabel.FindAttribute(TDataStd_Name::GetID(), aNameAttr))
  {
    aName           = new TCollection_HAsciiString(TCollection_AsciiString(aNameAttr->Get()));
    int aFirstSpace = aName->Search(" ");
    if (aFirstSpace != -1)
    {
      aName = aName->SubString(aFirstSpace + 1, aName->Length());
    }
    else
    {
      aName = new TCollection_HAsciiString();
    }
  }
  occ::handle<TCollection_HAsciiString> aDescription = new TCollection_HAsciiString();
  occ::handle<StepRepr_ShapeAspect>     aSA          = new StepRepr_ShapeAspect;
  aSA->Init(aName, aDescription, aPDS, StepData_LTrue);

  // Geometric_Item_Specific_Usage
  occ::handle<StepAP242_GeometricItemSpecificUsage> aGISU =
    new StepAP242_GeometricItemSpecificUsage();
  StepAP242_ItemIdentifiedRepresentationUsageDefinition aDefinition;
  aDefinition.SetValue(aSA);
  occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> anReprItems =
    new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, 1);
  occ::handle<StepRepr_RepresentationItem> anIdentifiedItem =
    occ::down_cast<StepRepr_RepresentationItem>(anEnt);
  if (anIdentifiedItem.IsNull())
  {
    // What the shape was written as is no item to identify (a tessellated
    // shape's face): a usage of a null item leaves a model that cannot be
    // written
    return nullptr;
  }
  anReprItems->SetValue(1, anIdentifiedItem);
  if (aUsedRep.IsNull())
  {
    occ::handle<StepShape_ShapeDefinitionRepresentation> aSDR;
    for (Interface_EntityIterator aSharingIter = aGraph.Sharings(aPDS);
         aSharingIter.More() && aSDR.IsNull();
         aSharingIter.Next())
    {
      const occ::handle<Standard_Transient>& anEntity = aSharingIter.Value();
      aSDR = occ::down_cast<StepShape_ShapeDefinitionRepresentation>(anEntity);
    }
    if (aSDR.IsNull())
    {
      return nullptr;
    }
    aUsedRep = aSDR->UsedRepresentation();
  }

  // Set entities to model
  aGISU->Init(aName, aDescription, aDefinition, aUsedRep, anReprItems);
  aModel->AddWithRefs(aSA);
  aModel->AddWithRefs(aGISU);
  theGISU = aGISU;
  return aSA;
}

//=======================================================================
// function : shapeAspects
// purpose  : a sequence of the one shape aspect
//=======================================================================
static NCollection_Sequence<occ::handle<StepRepr_ShapeAspect>> shapeAspects(
  const occ::handle<StepRepr_ShapeAspect>& theSA)
{
  NCollection_Sequence<occ::handle<StepRepr_ShapeAspect>> aSeq;
  aSeq.Append(theSA);
  return aSeq;
}

//! The kind of graphic a file named it (its set's name), kept to be written again; else the one
//! its semantics give
static const char* keptOr(const occ::handle<TCollection_HAsciiString>& theKept,
                          const char*                                  theType)
{
  return !theKept.IsNull() && theKept->Length() > 0 ? theKept->ToCString() : theType;
}

//=======================================================================
// function : setEdition2, isEdition2
// purpose  : the file identified as of AP242's second edition, where the
//            annotation placeholder and the model_geometric_view are from:
//            NIST's own are { 1 0 10303 442 3 1 4 } (CTC-02) and the third's
//=======================================================================
static const char* THE_AP242_ED2 =
  "AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 3 1 4 }";

static occ::handle<HeaderSection_FileSchema> fileSchema(
  const occ::handle<Interface_InterfaceModel>& theModel)
{
  occ::handle<StepData_StepModel> aStepModel = occ::down_cast<StepData_StepModel>(theModel);
  if (aStepModel.IsNull())
  {
    return nullptr;
  }
  return occ::down_cast<HeaderSection_FileSchema>(
    aStepModel->HeaderEntity(STANDARD_TYPE(HeaderSection_FileSchema)));
}

static void setEdition2(const occ::handle<Interface_InterfaceModel>& theModel)
{
  occ::handle<HeaderSection_FileSchema> aSchema = fileSchema(theModel);
  for (int i = 1; !aSchema.IsNull() && i <= aSchema->NbSchemaIdentifiers(); i++)
  {
    occ::handle<TCollection_HAsciiString> anId = aSchema->SchemaIdentifiersValue(i);
    if (!anId.IsNull() && anId->Search("AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF") == 1
        && anId->Search(" 442 1 1 4") > 0)
    {
      aSchema->SchemaIdentifiers()->SetValue(i, new TCollection_HAsciiString(THE_AP242_ED2));
    }
  }
}

static bool isEdition2(const occ::handle<Interface_InterfaceModel>& theModel)
{
  occ::handle<HeaderSection_FileSchema> aSchema = fileSchema(theModel);
  for (int i = 1; !aSchema.IsNull() && i <= aSchema->NbSchemaIdentifiers(); i++)
  {
    occ::handle<TCollection_HAsciiString> anId = aSchema->SchemaIdentifiersValue(i);
    if (!anId.IsNull() && anId->IsSameString(new TCollection_HAsciiString(THE_AP242_ED2)))
    {
      return true;
    }
  }
  return false;
}

//=======================================================================
// function : writePlaceholder
// purpose  : a PMI's annotation placeholder (AP242 ed2), for the one that has
//            one (XCAFDoc_DimTolTool::GetPlaceholder): where a reader that
//            draws PMI itself puts its text, nothing drawn. As NIST's STC-07
//            has it: a callout of the placeholder alone, named as the PMI's
//            graphic, linked to the PMI by an association with a
//            placeholder; the placeholder's geometric set named by its kind
//            ('linear dimension', 'note': what SFA counts it by), holding its
//            placement. No class is for the placeholder: it is written as an
//            entity of parameters
//=======================================================================
static void writePlaceholder(
  const occ::handle<XSControl_WorkSession>&                    theWS,
  NCollection_DynamicArray<occ::handle<StepRepr_RepresentationItem>>& theAnnotations,
                             const occ::handle<StepVisual_DraughtingModel>& theDM,
                             const TDF_Label&                               theLabel,
                             const occ::handle<TCollection_HAsciiString>&   thePrsName,
                             const occ::handle<Standard_Transient>&         theDefinition,
                             const bool                                     theHasPlane,
                             const gp_Ax2&                                  thePlane,
                             const gp_Pnt&                                  theTextPosition,
                             const StepData_Factors&                        theLocalFactors)
{
  TCollection_AsciiString aKind, aName;
  bool                    isShared = false;
  if (theLabel.IsNull() || theDefinition.IsNull() || theDM.IsNull()
      || !XCAFDoc_DimTolTool::GetPlaceholder(theLabel, aKind, aName, isShared))
  {
    return;
  }
  const occ::handle<Interface_InterfaceModel>& aModel = theWS->Model();
  // Named as it was, else as the PMI's graphic
  if (aName.IsEmpty() && !thePrsName.IsNull())
  {
    aName = thePrsName->String();
  }
  // Where its text is: the PMI's plane, at its text
  gp_Ax2 aPlace = theHasPlane ? thePlane : gp_Ax2();
  aPlace.SetLocation(theTextPosition);
  // Another's, written already of this kind and name, there: this PMI's too
  char aWhere[128];
  std::snprintf(aWhere,
                sizeof(aWhere),
                "\n%.6g %.6g %.6g",
                theTextPosition.X(),
                theTextPosition.Y(),
                theTextPosition.Z());
  const TCollection_AsciiString aKey = aKind + "\n" + aName + aWhere;
  const NCollection_Array1<occ::handle<Standard_Transient>>* aWritten =
    isShared ? gdtState(theDM).Placeholders.Seek(aKey) : nullptr;
  if (aWritten != nullptr)
  {
    occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> anItems =
      new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, 1);
    anItems->SetValue(1, occ::down_cast<StepRepr_RepresentationItem>(aWritten->Value(1)));
    StepAP242_ItemIdentifiedRepresentationUsageDefinition aDefinition;
    aDefinition.SetValue(theDefinition);
    occ::handle<StepAP242_DraughtingModelItemAssociationWithPlaceholder> aLink =
      new StepAP242_DraughtingModelItemAssociationWithPlaceholder();
    aLink->Init(new TCollection_HAsciiString(TCollection_AsciiString("PPH - ") + aName),
                new TCollection_HAsciiString(),
                aDefinition,
                theDM,
                anItems);
    aLink->SetAnnotationPlaceholder(aWritten->Value(2));
    aModel->AddWithRefs(aLink);
    return;
  }
  GeomToStep_MakeAxis2Placement3d               anAxisMaker(aPlace, theLocalFactors);
  const occ::handle<StepGeom_Axis2Placement3d>& anAxis = anAxisMaker.Value();
  anAxis->SetName(new TCollection_HAsciiString(TCollection_AsciiString("PPH - ") + aName));
  occ::handle<NCollection_HArray1<StepShape_GeometricSetSelect>> anElements =
    new NCollection_HArray1<StepShape_GeometricSetSelect>(1, 1);
  StepShape_GeometricSetSelect anElement;
  anElement.SetValue(anAxis);
  anElements->SetValue(1, anElement);
  occ::handle<StepShape_GeometricSet> aSet = new StepShape_GeometricSet();
  aSet->Init(new TCollection_HAsciiString(aKind), anElements);
  // Its style: none
  occ::handle<StepVisual_NullStyleMember> aNullStyle = new StepVisual_NullStyleMember();
  aNullStyle->SetEnumText(0, ".NULL.");
  StepVisual_PresentationStyleSelect aStyleItem;
  aStyleItem.SetValue(aNullStyle);
  occ::handle<NCollection_HArray1<StepVisual_PresentationStyleSelect>> aStyles =
    new NCollection_HArray1<StepVisual_PresentationStyleSelect>(1, 1);
  aStyles->SetValue(1, aStyleItem);
  occ::handle<StepVisual_PresentationStyleAssignment> aPrsStyle =
    new StepVisual_PresentationStyleAssignment();
  aPrsStyle->Init(aStyles);
  // The placeholder: its name, styles, item, role and line spacing
  TCollection_AsciiString aQuoted("'");
  for (int i = 1; i <= aName.Length(); i++)
  {
    aQuoted += aName.Value(i);
    if (aName.Value(i) == '\'')
    {
      aQuoted += '\'';
    }
  }
  aQuoted += "'";
  occ::handle<StepData_UndefinedEntity> aStyleList = new StepData_UndefinedEntity(true);
  aStyleList->UndefinedContent()->AddEntity(Interface_ParamIdent, aPrsStyle);
  occ::handle<StepData_UndefinedEntity> aPlaceholder = new StepData_UndefinedEntity();
  aPlaceholder->SetStepType("ANNOTATION_PLACEHOLDER_OCCURRENCE");
  occ::handle<Interface_UndefinedContent> aParams = aPlaceholder->UndefinedContent();
  aParams->AddLiteral(Interface_ParamText, new TCollection_HAsciiString(aQuoted));
  aParams->AddEntity(Interface_ParamSub, aStyleList);
  aParams->AddEntity(Interface_ParamIdent, aSet);
  aParams->AddLiteral(Interface_ParamEnum, new TCollection_HAsciiString(".GPS_DATA."));
  aParams->AddLiteral(Interface_ParamReal, new TCollection_HAsciiString("1."));
  // Its callout, and the link to what it is the placeholder of
  StepVisual_DraughtingCalloutElement aContent;
  aContent.SetValue(aPlaceholder);
  occ::handle<NCollection_HArray1<StepVisual_DraughtingCalloutElement>> aContents =
    new NCollection_HArray1<StepVisual_DraughtingCalloutElement>(1, 1);
  aContents->SetValue(1, aContent);
  occ::handle<StepVisual_DraughtingCallout> aCallout = new StepVisual_DraughtingCallout();
  aCallout->Init(new TCollection_HAsciiString(aName), aContents);
  occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> anItems =
    new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, 1);
  anItems->SetValue(1, aCallout);
  StepAP242_ItemIdentifiedRepresentationUsageDefinition aDefinition;
  aDefinition.SetValue(theDefinition);
  occ::handle<StepAP242_DraughtingModelItemAssociationWithPlaceholder> aLink =
    new StepAP242_DraughtingModelItemAssociationWithPlaceholder();
  aLink->Init(new TCollection_HAsciiString(TCollection_AsciiString("PPH - ") + aName),
              new TCollection_HAsciiString(),
              aDefinition,
              theDM,
              anItems);
  aLink->SetAnnotationPlaceholder(aPlaceholder);
  aModel->AddWithRefs(aLink);
  NCollection_Array1<occ::handle<Standard_Transient>> aThis(1, 2);
  aThis.SetValue(1, aCallout);
  aThis.SetValue(2, aPlaceholder);
  if (!gdtState(theDM).Placeholders.IsBound(aKey))
  {
    gdtState(theDM).Placeholders.Bind(aKey, aThis);
  }
  // A placeholder and its association are of AP242's second edition: the file says so (one
  // identified as the first edition and holding them is told off by SFA)
  setEdition2(aModel);
  // In an annotation plane, as every callout is (one in none is an error to SFA), and in the
  // draughting model with it
  occ::handle<NCollection_HArray1<occ::handle<StepVisual_PresentationStyleAssignment>>> aPrsStyles =
    new NCollection_HArray1<occ::handle<StepVisual_PresentationStyleAssignment>>(1, 1);
  aPrsStyles->SetValue(1, aPrsStyle);
  GeomToStep_MakeAxis2Placement3d aPlaneAxis(aPlace, theLocalFactors);
  occ::handle<StepGeom_Plane>     aPlane = new StepGeom_Plane();
  aPlane->Init(new TCollection_HAsciiString(), aPlaneAxis.Value());
  StepVisual_AnnotationPlaneElement aPlaneElement;
  aPlaneElement.SetValue(aCallout);
  occ::handle<NCollection_HArray1<StepVisual_AnnotationPlaneElement>> aPlaneElements =
    new NCollection_HArray1<StepVisual_AnnotationPlaneElement>(1, 1);
  aPlaneElements->SetValue(1, aPlaneElement);
  occ::handle<StepVisual_AnnotationPlane> anAnnPlane = new StepVisual_AnnotationPlane();
  anAnnPlane->Init(new TCollection_HAsciiString(), aPrsStyles, aPlane, aPlaneElements);
  aModel->AddWithRefs(anAnnPlane);
  theAnnotations.Append(aCallout);
  theAnnotations.Append(anAnnPlane);
}

//=======================================================================
// function : writePresentation
// purpose  : auxiliary (write annotation plane and presentation)
//======================================================================
void STEPCAFControl_Writer::writePresentation(
  const occ::handle<XSControl_WorkSession>&    theWS,
  const TDF_Label&                             theLabel,
  const TopoDS_Shape&                          thePresentation,
  const occ::handle<TCollection_HAsciiString>& thePrsName,
  const bool                                   theHasSemantic,
  const bool                                   theHasPlane,
  const gp_Ax2&                                theAnnotationPlane,
  const gp_Pnt&                                theTextPosition,
  const occ::handle<Standard_Transient>&       theDimension,
  const char*                                  thePrsType,
  const NCollection_Sequence<occ::handle<StepRepr_ShapeAspect>>& theGeometry,
  const StepData_Factors&                      theLocalFactors,
  const bool                                   theAlone)
{
  // A graphic of nothing, or of no item: a draughting model item association
  // without its definition leaves a model the writer cannot write. Alone -- a
  // graphic-only note on nothing, NIST CTC-05's rectangular regions -- the
  // callout is in the draughting model with no association: one to an empty
  // shape aspect a reader (SFA) takes for geometry it cannot find
  writePlaceholder(theWS,
                   myGDTAnnotations,
                   myGDTPresentationDM,
                   theLabel,
                   thePrsName,
                   theDimension,
                   theHasPlane,
                   theAnnotationPlane,
                   theTextPosition,
                   theLocalFactors);
  TDF_Label aHostL;
  if (thePresentation.IsNull() && !theDimension.IsNull() && theHasSemantic
      && XCAFDoc_DimTolTool::GetPresentedWith(theLabel, aHostL))
  {
    // Not drawn: another's graphic presents it too
    GDTPresentationModel::Link aLink;
    aLink.Label     = theLabel;
    aLink.Host      = aHostL;
    aLink.Dimension = theDimension;
    aLink.Geometry  = theGeometry;
    gdtState(myGDTPresentationDM).Links.Append(aLink);
    return;
  }
  if (thePresentation.IsNull() || (theDimension.IsNull() && !theAlone))
  {
    return;
  }
  // Get working data
  const occ::handle<Interface_InterfaceModel>& aModel = theWS->Model();

  // Presentation
  occ::handle<StepVisual_TessellatedGeometricSet> aGeomSet =
    STEPCAFControl_GDTProperty::GetTessellation(thePresentation, theLocalFactors.LengthFactor());
  if (aGeomSet.IsNull())
  {
    return;
  }
  occ::handle<TCollection_HAsciiString> aPrsName =
    thePrsName.IsNull() ? new TCollection_HAsciiString() : thePrsName;
  // The occurrence named as the callout, the set by what it presents: one of
  // the names of AP242 PMI recommended practice 8.4 ("datum", "position",
  // "diameter dimension" ...), by which a reader (SFA) tells what it is; a
  // name missing is an error to it (8.1.1)
  aGeomSet->SetName(new TCollection_HAsciiString(thePrsType ? thePrsType : "general note"));
  occ::handle<StepVisual_TessellatedAnnotationOccurrence> aTAO =
    new StepVisual_TessellatedAnnotationOccurrence();
  aTAO->Init(aPrsName, myGDTPrsCurveStyle, aGeomSet);
  StepVisual_DraughtingCalloutElement aDCElement;
  aDCElement.SetValue(aTAO);
  occ::handle<NCollection_HArray1<StepVisual_DraughtingCalloutElement>> aTAOs =
    new NCollection_HArray1<StepVisual_DraughtingCalloutElement>(1, 1);
  aTAOs->SetValue(1, aDCElement);
  occ::handle<StepVisual_DraughtingCallout> aDCallout = new StepVisual_DraughtingCallout();
  aDCallout->Init(aPrsName, aTAOs);
  occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> aDCsForDMIA =
    new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, 1);
  aDCsForDMIA->SetValue(1, aDCallout);
  myGDTAnnotations.Append(aDCallout);
  if (!theLabel.IsNull())
  {
    if (!gdtState(myGDTPresentationDM).Items.IsBound(theLabel))
    {
      gdtState(myGDTPresentationDM).Items.Bind(theLabel,
                         NCollection_Sequence<occ::handle<StepRepr_RepresentationItem>>());
    }
    gdtState(myGDTPresentationDM).Items.ChangeFind(theLabel).Append(aDCallout);
  }
  // Its associations, to what it presents and to its geometry: none for a
  // graphic alone, which still has its plane (an annotation in no
  // annotation plane is an error to SFA)
  if (!theDimension.IsNull())
  {
    StepAP242_ItemIdentifiedRepresentationUsageDefinition aDimension;
    aDimension.SetValue(theDimension);
    occ::handle<TCollection_HAsciiString> aDMIAName;
    if (theHasSemantic)
    {
      aDMIAName = new TCollection_HAsciiString("PMI representation to presentation link");
    }
    else
    {
      aDMIAName = new TCollection_HAsciiString();
    }
    occ::handle<StepAP242_DraughtingModelItemAssociation> aDMIA =
      new StepAP242_DraughtingModelItemAssociation();
    aDMIA->Init(aDMIAName,
                new TCollection_HAsciiString(),
                aDimension,
                myGDTPresentationDM,
                aDCsForDMIA);
    aModel->AddWithRefs(aDMIA);
    // And to what it is about: the shape aspects of the faces, each by an
    // association of its own (as the NIST files write it), for a reader to
    // find the geometry a graphic presents without going through the semantics
    for (NCollection_Sequence<occ::handle<StepRepr_ShapeAspect>>::Iterator anIt(theGeometry);
         anIt.More();
         anIt.Next())
    {
      if (anIt.Value().IsNull())
      {
        continue;
      }
      StepAP242_ItemIdentifiedRepresentationUsageDefinition aGeometry;
      aGeometry.SetValue(anIt.Value());
      occ::handle<StepAP242_DraughtingModelItemAssociation> aGeomDMIA =
        new StepAP242_DraughtingModelItemAssociation();
      aGeomDMIA->Init(new TCollection_HAsciiString(),
                      new TCollection_HAsciiString(),
                      aGeometry,
                      myGDTPresentationDM,
                      aDCsForDMIA);
      aModel->AddWithRefs(aGeomDMIA);
    }
  }

  if (!theHasPlane)
  {
    return;
  }

  // Annotation plane
  // Presentation Style
  occ::handle<StepVisual_NullStyleMember> aNullStyle = new StepVisual_NullStyleMember();
  aNullStyle->SetEnumText(0, ".NULL.");
  StepVisual_PresentationStyleSelect aStyleItem;
  aStyleItem.SetValue(aNullStyle);
  occ::handle<NCollection_HArray1<StepVisual_PresentationStyleSelect>> aStyles =
    new NCollection_HArray1<StepVisual_PresentationStyleSelect>(1, 1);
  aStyles->SetValue(1, aStyleItem);
  occ::handle<StepVisual_PresentationStyleAssignment> aPrsStyle =
    new StepVisual_PresentationStyleAssignment();
  aPrsStyle->Init(aStyles);
  occ::handle<NCollection_HArray1<occ::handle<StepVisual_PresentationStyleAssignment>>> aPrsStyles =
    new NCollection_HArray1<occ::handle<StepVisual_PresentationStyleAssignment>>(1, 1);
  aPrsStyles->SetValue(1, aPrsStyle);
  // Plane
  occ::handle<StepGeom_Plane>                   aPlane = new StepGeom_Plane();
  GeomToStep_MakeAxis2Placement3d               anAxisMaker(theAnnotationPlane, theLocalFactors);
  const occ::handle<StepGeom_Axis2Placement3d>& anAxis = anAxisMaker.Value();
  // Set text position to plane origin
  occ::handle<StepGeom_CartesianPoint>     aTextPos = new StepGeom_CartesianPoint();
  occ::handle<NCollection_HArray1<double>> aCoords  = new NCollection_HArray1<double>(1, 3);
  for (int i = 1; i <= 3; i++)
  {
    // In the file's length unit, as the plane and the graphic are
    aCoords->SetValue(i, theTextPosition.Coord(i) / theLocalFactors.LengthFactor());
  }
  aTextPos->Init(new TCollection_HAsciiString(), aCoords);
  anAxis->SetLocation(aTextPos);
  aPlane->Init(new TCollection_HAsciiString(), anAxis);
  // Annotation plane element
  StepVisual_AnnotationPlaneElement aPlaneElement;
  aPlaneElement.SetValue(aDCallout);
  occ::handle<NCollection_HArray1<StepVisual_AnnotationPlaneElement>> aDCsForAnnPln =
    new NCollection_HArray1<StepVisual_AnnotationPlaneElement>(1, 1);
  aDCsForAnnPln->SetValue(1, aPlaneElement);
  // Init AnnotationPlane entity
  occ::handle<StepVisual_AnnotationPlane> anAnnPlane = new StepVisual_AnnotationPlane();
  anAnnPlane->Init(new TCollection_HAsciiString(), aPrsStyles, aPlane, aDCsForAnnPln);
  myGDTAnnotations.Append(anAnnPlane);
  if (!theLabel.IsNull())
  {
    gdtState(myGDTPresentationDM).Items.ChangeFind(theLabel).Append(anAnnPlane);
  }
  aModel->AddWithRefs(anAnnPlane);
}

//=======================================================================
// function : writeDatumAP242
// purpose  : auxiliary (write Datum entity for given shape or write all
//           necessary entities and link them to already written datum
//           in case of multiple features association)
//=======================================================================
occ::handle<StepDimTol_Datum> STEPCAFControl_Writer::writeDatumAP242(
  const occ::handle<XSControl_WorkSession>& theWS,
  const NCollection_Sequence<TDF_Label>&    theShapeL,
  const TDF_Label&                          theDatumL,
  const bool                                theIsFirstDTarget,
  const occ::handle<StepDimTol_Datum>&      theWrittenDatum,
  const StepData_Factors&                   theLocalFactors)
{
  // Get working data
  const occ::handle<Interface_InterfaceModel>& aModel  = theWS->Model();
  const occ::handle<XSControl_TransferWriter>& aTW     = theWS->TransferWriter();
  const occ::handle<Transfer_FinderProcess>&   aFP     = aTW->FinderProcess();
  const occ::handle<Interface_HGraph>          aHGraph = theWS->HGraph();
  if (aHGraph.IsNull())
  {
    return nullptr;
  }
  Interface_Graph aGraph = aHGraph->Graph();

  occ::handle<StepRepr_ShapeAspect>                       aSA;
  occ::handle<StepRepr_RepresentationContext>             aRC;
  occ::handle<StepRepr_ProductDefinitionShape>            aPDS;
  NCollection_Sequence<occ::handle<StepRepr_ShapeAspect>> aSASeq;
  occ::handle<StepAP242_GeometricItemSpecificUsage>       aGISU;
  int                                                     aSANum = 0, aGISUNum = 0;
  // Link with datum feature
  for (NCollection_Sequence<TDF_Label>::Iterator aLabelIter(theShapeL); aLabelIter.More();
       aLabelIter.Next())
  {
    occ::handle<Standard_Transient>                       anEnt;
    TopLoc_Location                                       aLoc;
    NCollection_Sequence<occ::handle<Standard_Transient>> aSeqRI;

    TopoDS_Shape aShape = XCAFDoc_ShapeTool::GetShape(aLabelIter.Value());
    FindEntities(aFP, aShape, aLoc, aSeqRI);
    if (aSeqRI.Length() <= 0)
    {
      Message::SendTrace() << "Warning: Cannot find RI for "
                           << aShape.TShape()->DynamicType()->Name() << "\n";
      continue;
    }
    anEnt = aSeqRI.Value(1);
    aPDS  = FindPDS(aGraph, anEnt, aRC);
    if (aPDS.IsNull())
    {
      continue;
    }

    occ::handle<StepRepr_ShapeAspect> aCurrentSA =
      writeShapeAspect(theWS, theDatumL, aShape, aRC, aGISU);
    if (aCurrentSA.IsNull())
    {
      continue;
    }
    aSASeq.Append(aCurrentSA);
    aSANum   = aModel->Number(aCurrentSA);
    aGISUNum = aModel->Number(aGISU);
  }
  if (aPDS.IsNull())
  {
    // Workaround for datums without shape
    aPDS                                                        = myGDTCommonPDS;
    Interface_EntityIterator                             aSDRIt = aGraph.Sharings(aPDS);
    occ::handle<StepShape_ShapeDefinitionRepresentation> aSDR;
    for (aSDRIt.Start(); aSDRIt.More() && aSDR.IsNull(); aSDRIt.Next())
    {
      aSDR = occ::down_cast<StepShape_ShapeDefinitionRepresentation>(aSDRIt.Value());
    }
    if (!aSDR.IsNull())
    {
      occ::handle<StepRepr_Representation> aRepr = aSDR->UsedRepresentation();
      if (!aRepr.IsNull())
      {
        aRC = aRepr->ContextOfItems();
      }
    }
  }

  // Find if datum has datum targets and get common datum attributes
  occ::handle<XCAFDoc_Datum> aDatumAttr;
  if (!theDatumL.FindAttribute(XCAFDoc_Datum::GetID(), aDatumAttr))
  {
    return nullptr;
  }
  occ::handle<XCAFDimTolObjects_DatumObject> anObject = aDatumAttr->GetObject();
  if (anObject.IsNull())
  {
    return nullptr;
  }
  bool                                  isSimpleDatum = !anObject->IsDatumTarget();
  occ::handle<TCollection_HAsciiString> anIdentifier  = anObject->GetName();
  occ::handle<TCollection_HAsciiString> aTargetId =
    (anObject->GetDatumTargetNumber() == 0
       ? new TCollection_HAsciiString()
       : new TCollection_HAsciiString(anObject->GetDatumTargetNumber()));

  // If datum type is area, but there is no area in object, write as simple datum
  if (anObject->GetDatumTargetType() == XCAFDimTolObjects_DatumTargetType_Area
      && anObject->GetDatumTarget().IsNull())
  {
    isSimpleDatum = true;
  }

  // Simple datum
  if (isSimpleDatum)
  {
    if (aSASeq.Length() == 0)
    {
      // Create empty datum with name and presentation only
      occ::handle<StepDimTol_DatumFeature> aDF = new StepDimTol_DatumFeature();
      aDF->Init(new TCollection_HAsciiString(),
                new TCollection_HAsciiString(),
                aPDS,
                StepData_LTrue);
      aSA = aDF;
      aModel->AddWithRefs(aDF);
    }
    else if (aSASeq.Length() == 1)
    {
      occ::handle<StepDimTol_DatumFeature> aDF = new StepDimTol_DatumFeature();
      aDF->Init(new TCollection_HAsciiString(),
                new TCollection_HAsciiString(),
                aPDS,
                StepData_LTrue);
      aModel->ReplaceEntity(aSANum, aDF);
      aSA = aDF;
      StepAP242_ItemIdentifiedRepresentationUsageDefinition aDefinition;
      aDefinition.SetValue(aDF);
      aGISU->SetDefinition(aDefinition);
      aModel->ReplaceEntity(aGISUNum, aGISU);
    }
    else if (aSASeq.Length() > 1)
    {
      occ::handle<StepRepr_CompShAspAndDatumFeatAndShAsp> aDF =
        new StepRepr_CompShAspAndDatumFeatAndShAsp();
      aDF->Init(new TCollection_HAsciiString(),
                new TCollection_HAsciiString(),
                aPDS,
                StepData_LTrue);
      for (int i = 1; i <= aSASeq.Length(); i++)
      {
        occ::handle<StepRepr_ShapeAspectRelationship> aSAR = new StepRepr_ShapeAspectRelationship();
        aSAR->Init(new TCollection_HAsciiString(),
                   false,
                   new TCollection_HAsciiString(),
                   aDF,
                   aSASeq.Value(i));
        aModel->AddWithRefs(aSAR);
      }
      aSA = aDF;
      aModel->AddWithRefs(aDF);
    }
  }
  // Datum with datum targets
  else
  {
    XCAFDimTolObjects_DatumTargetType   aDatumType = anObject->GetDatumTargetType();
    occ::handle<StepDimTol_DatumTarget> aDatumTarget;
    // Note: the given way to write such datum type may be incorrect (too little information)
    if (aDatumType == XCAFDimTolObjects_DatumTargetType_Area)
    {
      TopoDS_Shape                                      aDTShape = anObject->GetDatumTarget();
      occ::handle<StepAP242_GeometricItemSpecificUsage> anAreaGISU;
      occ::handle<StepRepr_ShapeAspect>                 anAreaSA =
        writeShapeAspect(theWS, theDatumL, aDTShape, aRC, anAreaGISU);
      aSANum                                  = aModel->Number(anAreaSA);
      aGISUNum                                = aModel->Number(anAreaGISU);
      occ::handle<StepDimTol_DatumTarget> aDT = new StepDimTol_DatumTarget();
      aDT->Init(new TCollection_HAsciiString(),
                new TCollection_HAsciiString("area"),
                anAreaSA->OfShape(),
                StepData_LTrue,
                aTargetId);
      aModel->ReplaceEntity(aSANum, aDT);
      StepAP242_ItemIdentifiedRepresentationUsageDefinition aDefinition;
      aDefinition.SetValue(aDT);
      anAreaGISU->SetDefinition(aDefinition);
      aModel->ReplaceEntity(aGISUNum, anAreaGISU);
      // The target it is: related to its datum below (it was not, and the area was lost)
      aDatumTarget = aDT;
    }
    else
    {
      occ::handle<StepDimTol_PlacedDatumTargetFeature> aPDTF =
        new StepDimTol_PlacedDatumTargetFeature();
      aPDTF->Init(new TCollection_HAsciiString(),
                  STEPCAFControl_GDTProperty::GetDatumTargetName(aDatumType),
                  aPDS,
                  StepData_LTrue,
                  aTargetId);
      aModel->AddWithRefs(aPDTF);
      aDatumTarget = aPDTF;
      // Datum targets
      occ::handle<StepRepr_PropertyDefinition> aPD = new StepRepr_PropertyDefinition();
      StepRepr_CharacterizedDefinition         aCDefinition;
      aCDefinition.SetValue(aPDTF);
      aPD->Init(new TCollection_HAsciiString(), false, nullptr, aCDefinition);
      if (anObject->HasDatumTargetParams())
      {
        // write all parameters of datum target
        occ::handle<StepShape_ShapeRepresentationWithParameters> aSRWP =
          new StepShape_ShapeRepresentationWithParameters();
        // Common for all datum targets
        StepBasic_Unit                         aUnit   = GetUnit(aRC);
        gp_Ax2                                 aDTAxis = anObject->GetDatumTargetAxis();
        GeomToStep_MakeAxis2Placement3d        anAxisMaker(aDTAxis, theLocalFactors);
        occ::handle<StepGeom_Axis2Placement3d> anA2P3D = anAxisMaker.Value();
        anA2P3D->SetName(new TCollection_HAsciiString("orientation"));
        occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> anItems;
        // Process each datum target type
        if (aDatumType == XCAFDimTolObjects_DatumTargetType_Point)
        {
          anItems = new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, 1);
        }
        else
        {
          occ::handle<TCollection_HAsciiString> aTargetValueName;
          if (aDatumType == XCAFDimTolObjects_DatumTargetType_Line)
          {
            anItems = new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, 2);
            aTargetValueName = new TCollection_HAsciiString("target length");
          }
          else if (aDatumType == XCAFDimTolObjects_DatumTargetType_Rectangle)
          {
            anItems = new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, 3);
            aTargetValueName = new TCollection_HAsciiString("target length");
            // Additional value
            occ::handle<StepRepr_ReprItemAndMeasureWithUnit> aTargetValue =
              CreateDimValue(anObject->GetDatumTargetWidth() / theLocalFactors.LengthFactor(),
                             aUnit,
                             new TCollection_HAsciiString("target width"),
                             "POSITIVE_LENGTH_MEASURE",
                             false);
            anItems->SetValue(2, aTargetValue);
            aModel->AddWithRefs(aTargetValue);
          }
          else if (aDatumType == XCAFDimTolObjects_DatumTargetType_Circle
                   || aDatumType == XCAFDimTolObjects_DatumTargetType_CircularCurve)
          {
            anItems = new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, 2);
            aTargetValueName = new TCollection_HAsciiString("target diameter");
          }
          // Value
          occ::handle<StepRepr_ReprItemAndMeasureWithUnit> aTargetValue =
            CreateDimValue(anObject->GetDatumTargetLength() / theLocalFactors.LengthFactor(),
                           aUnit,
                           aTargetValueName,
                           "POSITIVE_LENGTH_MEASURE",
                           false);
          anItems->SetValue(1, aTargetValue);
          aModel->AddWithRefs(aTargetValue);
        }
        anItems->SetValue(anItems->Length(), anA2P3D);
        aSRWP->Init(new TCollection_HAsciiString(), anItems, aRC);
        // Create and write auxiliary entities
        occ::handle<StepShape_ShapeDefinitionRepresentation> aSDR =
          new StepShape_ShapeDefinitionRepresentation();
        StepRepr_RepresentedDefinition aRDefinition;
        aRDefinition.SetValue(aPD);
        aSDR->Init(aRDefinition, aSRWP);
        aModel->AddWithRefs(aPD);
        aModel->AddWithRefs(aSRWP);
        aModel->AddWithRefs(aSDR);
      }
    }
    // Link datum target to datum feature
    // if aSASeq.Length() == 0 nothing to do
    if (aSASeq.Length() == 1)
    {
      occ::handle<StepRepr_FeatureForDatumTargetRelationship> aFFDTR =
        new StepRepr_FeatureForDatumTargetRelationship();
      aFFDTR->Init(new TCollection_HAsciiString(), false, nullptr, aSASeq.Value(1), aDatumTarget);
      aModel->AddWithRefs(aFFDTR);
    }
    else if (aSASeq.Length() > 1)
    {
      occ::handle<StepRepr_CompositeShapeAspect> aCompSA = new StepRepr_CompositeShapeAspect();
      aCompSA->Init(new TCollection_HAsciiString(),
                    new TCollection_HAsciiString(),
                    aPDS,
                    aSASeq.Value(1)->ProductDefinitional());
      for (int i = 1; i <= aSASeq.Length(); i++)
      {
        occ::handle<StepRepr_ShapeAspectRelationship> aSAR = new StepRepr_ShapeAspectRelationship();
        aSAR->Init(new TCollection_HAsciiString(),
                   false,
                   new TCollection_HAsciiString(),
                   aCompSA,
                   aSASeq.Value(i));
        aModel->AddWithRefs(aSAR);
      }
      occ::handle<StepRepr_FeatureForDatumTargetRelationship> aFFDTR =
        new StepRepr_FeatureForDatumTargetRelationship();
      aFFDTR->Init(new TCollection_HAsciiString(), false, nullptr, aCompSA, aDatumTarget);
      aModel->AddWithRefs(aFFDTR);
    }
    aSA = aDatumTarget;
  }

  // Datum
  occ::handle<StepDimTol_Datum> aDatum = theWrittenDatum;
  if (theIsFirstDTarget)
  {
    aDatum = new StepDimTol_Datum();
    aDatum->Init(new TCollection_HAsciiString(),
                 new TCollection_HAsciiString(),
                 aPDS,
                 StepData_LFalse,
                 anIdentifier);
    aModel->AddWithRefs(aDatum);
  }

  // Shape_Aspect_Relationship; a target's described as one, 'datum target', as NIST's files
  // have it
  if (!aSA.IsNull())
  {
    occ::handle<StepRepr_ShapeAspectRelationship> aSAR = new StepRepr_ShapeAspectRelationship();
    if (anObject->IsDatumTarget())
    {
      aSAR->Init(new TCollection_HAsciiString(),
                 true,
                 new TCollection_HAsciiString("datum target"),
                 aSA,
                 aDatum);
    }
    else
    {
      aSAR->Init(new TCollection_HAsciiString(), false, nullptr, aSA, aDatum);
    }
    aModel->AddWithRefs(aSAR);
  }

  // Annotation plane and Presentation
  writePresentation(theWS,
                    theDatumL,
                    anObject->GetPresentation(),
                    anObject->GetPresentationName(),
                    true,
                    anObject->HasPlane(),
                    anObject->GetPlane(),
                    anObject->GetPointTextAttach(),
                    aSA,
                    keptOr(anObject->GetPresentationType(),
                           anObject->IsDatumTarget() ? "datum target" : "datum"),
                    // The geometry a target is on: its feature's shape aspects, which the
                    // faces are of, not the target, which has none
                    anObject->IsDatumTarget() && aSASeq.Length() > 0 ? aSASeq
                                                                     : shapeAspects(aSA),
                    theLocalFactors);

  return aDatum;
}

//=======================================================================
// function : WriteDimValues
// purpose  : auxiliary (write all data for given dimension: values,
//           qualifiers, modifiers, orientation and tolerance class)
//======================================================================
static void WriteDimValues(const occ::handle<XSControl_WorkSession>&             theWS,
                           const occ::handle<XCAFDimTolObjects_DimensionObject>& theObject,
                           const occ::handle<StepRepr_RepresentationContext>&    theRC,
                           const StepShape_DimensionalCharacteristic&            theDimension,
                           const StepData_Factors&                               theLocalFactors)
{
  // Get working data
  const occ::handle<Interface_InterfaceModel>&           aModel     = theWS->Model();
  NCollection_Sequence<XCAFDimTolObjects_DimensionModif> aModifiers = theObject->GetModifiers();
  const occ::handle<Standard_Transient>&                 aDim       = theDimension.Value();
  bool isAngle = aDim->IsKind(STANDARD_TYPE(StepShape_AngularLocation))
                 || aDim->IsKind(STANDARD_TYPE(StepShape_AngularSize));

  // Unit; a length in the file's (an inch file's in inches), an angle in degrees as it is
  StepBasic_Unit aUnit = GetUnit(theRC, isAngle);
  const double   aLF   = isAngle ? 1. : theLocalFactors.LengthFactor();
  const char*    aMeasureName;
  if (isAngle)
  {
    aMeasureName = "POSITIVE_PLANE_ANGLE_MEASURE";
  }
  else
  {
    aMeasureName = "POSITIVE_LENGTH_MEASURE";
  }

  // Values
  occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> aValues;
  int aNbItems = 1, aValIt = 1;
  if (theObject->IsDimWithRange())
  {
    aNbItems += 2;
  }
  if (aModifiers.Length() > 0)
  {
    aNbItems++;
  }
  if (theObject->GetType() == XCAFDimTolObjects_DimensionType_Location_Oriented)
  {
    aNbItems++;
  }
  aNbItems += theObject->NbDescriptions();
  aValues = new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, aNbItems);

  // Nominal value
  double aNominal = theObject->GetValue() / aLF;
  int    aLeftNbDigits, aRightNbDigits;
  theObject->GetNbOfDecimalPlaces(aLeftNbDigits, aRightNbDigits);
  if (aLeftNbDigits > 0 && aLF != 1.)
  {
    // The digits before the point are those of the value as written, in the file's unit: 12 in
    // has two, the 304.8 mm it is kept as three (an 'NR2 3.2' of 12 is more than it needs)
    aLeftNbDigits = 1;
    for (double aWhole = std::abs(aNominal); aWhole >= 10. && aLeftNbDigits < 15; aWhole /= 10.)
    {
      aLeftNbDigits++;
    }
  }
  int aNbQualifiers = 0;
  if (theObject->HasQualifier() && !isAngle)
  {
    aNbQualifiers++;
  }
  if (aLeftNbDigits > 0 || aRightNbDigits > 0)
  {
    aNbQualifiers++;
  }
  // With qualifiers
  if (aNbQualifiers > 0)
  {
    occ::handle<StepShape_QualifiedRepresentationItem> aQRI =
      new StepShape_QualifiedRepresentationItem();
    occ::handle<NCollection_HArray1<StepShape_ValueQualifier>> aQualifiers =
      new NCollection_HArray1<StepShape_ValueQualifier>(1, aNbQualifiers);
    // Type qualifier
    if (theObject->HasQualifier() && !isAngle)
    {
      StepShape_ValueQualifier             anItem;
      occ::handle<StepShape_TypeQualifier> aType      = new StepShape_TypeQualifier();
      XCAFDimTolObjects_DimensionQualifier aQualifier = theObject->GetQualifier();
      aType->Init(STEPCAFControl_GDTProperty::GetDimQualifierName(aQualifier));
      aModel->AddWithRefs(aType);
      anItem.SetValue(aType);
      aQualifiers->SetValue(1, anItem);
    }
    // Number of decimal places
    if (aLeftNbDigits > 0 || aRightNbDigits > 0)
    {
      StepShape_ValueQualifier                        anItem;
      occ::handle<StepShape_ValueFormatTypeQualifier> aType =
        new StepShape_ValueFormatTypeQualifier();
      occ::handle<TCollection_HAsciiString> aFormatType = new TCollection_HAsciiString("NR2 ");
      aFormatType->AssignCat(new TCollection_HAsciiString(aLeftNbDigits));
      aFormatType->AssignCat(new TCollection_HAsciiString("."));
      aFormatType->AssignCat(new TCollection_HAsciiString(aRightNbDigits));
      aType->Init(aFormatType);
      aModel->AddWithRefs(aType);
      anItem.SetValue(aType);
      aQualifiers->SetValue(aNbQualifiers, anItem);
    }
    // Set qualifiers
    aQRI->SetQualifiers(aQualifiers);
    occ::handle<StepRepr_ReprItemAndMeasureWithUnit> anItem =
      CreateDimValue(aNominal,
                     aUnit,
                     new TCollection_HAsciiString("nominal value"),
                     aMeasureName,
                     isAngle,
                     true,
                     aQRI);
    aValues->SetValue(aValIt, anItem);
    aValIt++;
  }
  // Without qualifiers
  else
  {
    occ::handle<StepRepr_ReprItemAndMeasureWithUnit> anItem =
      CreateDimValue(aNominal,
                     aUnit,
                     new TCollection_HAsciiString("nominal value"),
                     aMeasureName,
                     isAngle);
    aValues->SetValue(aValIt, anItem);
    aValIt++;
  }

  // Ranges
  if (theObject->IsDimWithRange())
  {
    occ::handle<StepRepr_ReprItemAndMeasureWithUnit> aLowerItem =
      CreateDimValue(theObject->GetLowerBound() / aLF,
                     aUnit,
                     new TCollection_HAsciiString("lower limit"),
                     aMeasureName,
                     isAngle);
    occ::handle<StepRepr_ReprItemAndMeasureWithUnit> anUpperItem =
      CreateDimValue(theObject->GetUpperBound() / aLF,
                     aUnit,
                     new TCollection_HAsciiString("upper limit"),
                     aMeasureName,
                     isAngle);
    aValues->SetValue(aValIt, aLowerItem);
    aValIt++;
    aValues->SetValue(aValIt, anUpperItem);
    aValIt++;
  }

  // Modifiers
  if (aModifiers.Length() > 0)
  {
    occ::handle<StepRepr_CompoundRepresentationItem> aCompoundRI =
      new StepRepr_CompoundRepresentationItem();
    occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> aModifItems =
      new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, aModifiers.Length());
    for (int i = 1; i <= aModifiers.Length(); i++)
    {
      XCAFDimTolObjects_DimensionModif                    aModif = aModifiers.Value(i);
      occ::handle<StepRepr_DescriptiveRepresentationItem> aModifItem =
        new StepRepr_DescriptiveRepresentationItem();
      aModifItem->Init(new TCollection_HAsciiString(),
                       STEPCAFControl_GDTProperty::GetDimModifierName(aModif));
      aModel->AddWithRefs(aModifItem);
      aModifItems->SetValue(i, aModifItem);
    }
    aCompoundRI->Init(new TCollection_HAsciiString(), aModifItems);
    aValues->SetValue(aValIt, aCompoundRI);
    aValIt++;
  }

  // Orientation
  if (theObject->GetType() == XCAFDimTolObjects_DimensionType_Location_Oriented)
  {
    occ::handle<StepGeom_Axis2Placement3d> anOrientation = new StepGeom_Axis2Placement3d();
    gp_Dir                                 aDir;
    theObject->GetDirection(aDir);
    GeomToStep_MakeCartesianPoint MkPoint(gp_Pnt(0, 0, 0), theLocalFactors.LengthFactor());
    const occ::handle<StepGeom_CartesianPoint>& aLoc    = MkPoint.Value();
    occ::handle<StepGeom_Direction>             anAxis  = new StepGeom_Direction();
    occ::handle<NCollection_HArray1<double>>    aCoords = new NCollection_HArray1<double>(1, 3);
    aCoords->SetValue(1, aDir.X());
    aCoords->SetValue(2, aDir.Y());
    aCoords->SetValue(3, aDir.Z());
    anAxis->Init(new TCollection_HAsciiString(), aCoords);
    anOrientation
      ->Init(new TCollection_HAsciiString("orientation"), aLoc, true, anAxis, false, nullptr);
    aValues->SetValue(aValIt, anOrientation);
    aValIt++;
  }

  // Descriptions
  if (theObject->HasDescriptions())
  {
    for (int i = 0; i < theObject->NbDescriptions(); i++)
    {
      occ::handle<StepRepr_DescriptiveRepresentationItem> aDRI =
        new StepRepr_DescriptiveRepresentationItem();
      aDRI->Init(theObject->GetDescriptionName(i), theObject->GetDescription(i));
      aValues->SetValue(aValIt, aDRI);
      aValIt++;
    }
  }

  for (int i = 1; i <= aValues->Length(); i++)
  {
    aModel->AddWithRefs(aValues->Value(i));
  }

  // Create resulting Shape_Dimension_Representation
  occ::handle<StepShape_ShapeDimensionRepresentation> aSDR =
    new StepShape_ShapeDimensionRepresentation();
  aSDR->Init(new TCollection_HAsciiString(), aValues, theRC);
  aModel->AddWithRefs(aSDR);
  occ::handle<StepShape_DimensionalCharacteristicRepresentation> aDCR =
    new StepShape_DimensionalCharacteristicRepresentation();
  aDCR->Init(theDimension, aSDR);
  aModel->AddWithRefs(aDCR);

  // Plus_Minus_Tolerance
  if (theObject->IsDimWithPlusMinusTolerance())
  {
    occ::handle<TCollection_HAsciiString> aDummyName = new TCollection_HAsciiString(aMeasureName);
    // Delete "POSITIVE_", 9 chars
    aDummyName             = aDummyName->SubString(10, aDummyName->Length());
    aMeasureName           = aDummyName->ToCString();
    double aLowerTolValue  = -theObject->GetLowerTolValue() / aLF,
           anUpperTolValue = theObject->GetUpperTolValue() / aLF;
    // Upper
    // Typed, as the measure is: PLANE_ANGLE_ or LENGTH_MEASURE_WITH_UNIT
    occ::handle<StepBasic_MeasureWithUnit> anUpperMWU =
      isAngle ? occ::handle<StepBasic_MeasureWithUnit>(new StepBasic_PlaneAngleMeasureWithUnit())
              : occ::handle<StepBasic_MeasureWithUnit>(new StepBasic_LengthMeasureWithUnit());
    occ::handle<StepBasic_MeasureValueMember> anUpperValue = new StepBasic_MeasureValueMember();
    anUpperValue->SetName(aMeasureName);
    anUpperValue->SetReal(anUpperTolValue);
    anUpperMWU->Init(anUpperValue, aUnit);
    aModel->AddWithRefs(anUpperMWU);
    // Lower
    // Typed, as the measure is: PLANE_ANGLE_ or LENGTH_MEASURE_WITH_UNIT
    occ::handle<StepBasic_MeasureWithUnit> aLowerMWU =
      isAngle ? occ::handle<StepBasic_MeasureWithUnit>(new StepBasic_PlaneAngleMeasureWithUnit())
              : occ::handle<StepBasic_MeasureWithUnit>(new StepBasic_LengthMeasureWithUnit());
    occ::handle<StepBasic_MeasureValueMember> aLowerValue = new StepBasic_MeasureValueMember();
    aLowerValue->SetName(aMeasureName);
    aLowerValue->SetReal(aLowerTolValue);
    aLowerMWU->Init(aLowerValue, aUnit);
    aModel->AddWithRefs(aLowerMWU);
    // Tolerance
    occ::handle<StepShape_ToleranceValue> aTolValue = new StepShape_ToleranceValue();
    aTolValue->Init(aLowerMWU, anUpperMWU);
    aModel->AddWithRefs(aTolValue);
    StepShape_ToleranceMethodDefinition aMethod;
    aMethod.SetValue(aTolValue);
    occ::handle<StepShape_PlusMinusTolerance> aPlusMinusTol = new StepShape_PlusMinusTolerance();
    aPlusMinusTol->Init(aMethod, theDimension);
    aModel->AddWithRefs(aPlusMinusTol);
  }
  // Tolerance class
  if (theObject->IsDimWithClassOfTolerance())
  {
    bool                                    isHole;
    XCAFDimTolObjects_DimensionFormVariance aFormVariance;
    XCAFDimTolObjects_DimensionGrade        aGrade;
    if (!theObject->GetClassOfTolerance(isHole, aFormVariance, aGrade))
    {
      return;
    }
    occ::handle<StepShape_LimitsAndFits> aLAF =
      STEPCAFControl_GDTProperty::GetLimitsAndFits(isHole, aFormVariance, aGrade);
    aModel->AddWithRefs(aLAF);
    StepShape_ToleranceMethodDefinition aMethod;
    aMethod.SetValue(aLAF);
    occ::handle<StepShape_PlusMinusTolerance> aPlusMinusTol = new StepShape_PlusMinusTolerance();
    aPlusMinusTol->Init(aMethod, theDimension);
    aModel->AddWithRefs(aPlusMinusTol);
  }
}

//=======================================================================
// function : WriteDerivedGeometry
// purpose  : auxiliary (write connection point for dimensions)
//======================================================================
static void WriteDerivedGeometry(
  const occ::handle<XSControl_WorkSession>&                       theWS,
  const occ::handle<XCAFDimTolObjects_DimensionObject>&           theObject,
  const occ::handle<StepRepr_ConstructiveGeometryRepresentation>& theRepr,
  occ::handle<StepRepr_ShapeAspect>&                              theFirstSA,
  occ::handle<StepRepr_ShapeAspect>&                              theSecondSA,
  NCollection_DynamicArray<occ::handle<StepGeom_CartesianPoint>>& thePnts,
  const StepData_Factors&                                         theLocalFactors)
{
  const occ::handle<Interface_InterfaceModel>& aModel = theWS->Model();
  // First point
  if (theObject->HasPoint())
  {
    GeomToStep_MakeCartesianPoint        aPointMaker(theObject->GetPoint(),
                                              theLocalFactors.LengthFactor());
    occ::handle<StepGeom_CartesianPoint> aPoint = aPointMaker.Value();
    thePnts.Append(aPoint);
    occ::handle<StepRepr_DerivedShapeAspect> aDSA = new StepRepr_DerivedShapeAspect();
    aDSA->Init(new TCollection_HAsciiString(),
               new TCollection_HAsciiString(),
               theFirstSA->OfShape(),
               StepData_LFalse);
    occ::handle<StepAP242_GeometricItemSpecificUsage> aGISU =
      new StepAP242_GeometricItemSpecificUsage();
    StepAP242_ItemIdentifiedRepresentationUsageDefinition aDefinition;
    aDefinition.SetValue(aDSA);
    occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> anItem =
      new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, 1);
    anItem->SetValue(1, aPoint);
    aGISU->Init(new TCollection_HAsciiString(),
                new TCollection_HAsciiString(),
                aDefinition,
                theRepr,
                anItem);
    occ::handle<StepRepr_ShapeAspectDerivingRelationship> aSADR =
      new StepRepr_ShapeAspectDerivingRelationship();
    aSADR->Init(new TCollection_HAsciiString(),
                false,
                new TCollection_HAsciiString(),
                aDSA,
                theFirstSA);
    theFirstSA = aDSA;
    aModel->AddWithRefs(aGISU);
    aModel->AddWithRefs(aSADR);
  }

  // Second point (for locations)
  if (theObject->HasPoint2())
  {
    GeomToStep_MakeCartesianPoint        aPointMaker(theObject->GetPoint2(),
                                              theLocalFactors.LengthFactor());
    occ::handle<StepGeom_CartesianPoint> aPoint = aPointMaker.Value();
    thePnts.Append(aPoint);
    occ::handle<StepRepr_DerivedShapeAspect> aDSA = new StepRepr_DerivedShapeAspect();
    aDSA->Init(new TCollection_HAsciiString(),
               new TCollection_HAsciiString(),
               theFirstSA->OfShape(),
               StepData_LFalse);
    occ::handle<StepAP242_GeometricItemSpecificUsage> aGISU =
      new StepAP242_GeometricItemSpecificUsage();
    StepAP242_ItemIdentifiedRepresentationUsageDefinition aDefinition;
    aDefinition.SetValue(aDSA);
    occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> anItem =
      new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, 1);
    anItem->SetValue(1, aPoint);
    aGISU->Init(new TCollection_HAsciiString(),
                new TCollection_HAsciiString(),
                aDefinition,
                theRepr,
                anItem);
    occ::handle<StepRepr_ShapeAspectDerivingRelationship> aSADR =
      new StepRepr_ShapeAspectDerivingRelationship();
    aSADR->Init(new TCollection_HAsciiString(),
                false,
                new TCollection_HAsciiString(),
                aDSA,
                theSecondSA);
    theSecondSA = aDSA;
    aModel->AddWithRefs(aGISU);
    aModel->AddWithRefs(aSADR);
  }
}

//=======================================================================
// function : datumPartKey
// purpose  : the part a datum or a tolerance is on, as a key: the entry of the shape its first
//           shape is of, empty when it is on none
//=======================================================================
static TCollection_AsciiString datumPartKey(const TDF_Label& theGDTL)
{
  TCollection_AsciiString         aKey;
  NCollection_Sequence<TDF_Label> aFirst, aSecond;
  if (!XCAFDoc_DimTolTool::GetRefShapeLabel(theGDTL, aFirst, aSecond) || aFirst.IsEmpty())
  {
    return aKey;
  }
  TDF_Label aPartL = aFirst.First();
  if (XCAFDoc_ShapeTool::IsSubShape(aPartL))
  {
    aPartL = aPartL.Father();
  }
  TDF_Tool::Entry(aPartL, aKey);
  return aKey;
}

//=======================================================================
// function : WriteDatumSystem
// purpose  : auxiliary (write Write datum system for given
//           geometric_tolerance)
//======================================================================
static occ::handle<NCollection_HArray1<StepDimTol_DatumSystemOrReference>> WriteDatumSystem(
  const occ::handle<XSControl_WorkSession>&                                            theWS,
  const TDF_Label                                                                      theGeomTolL,
  const NCollection_Sequence<TDF_Label>&                                               theDatumSeq,
  const NCollection_DataMap<TCollection_AsciiString, occ::handle<Standard_Transient>>& theDatumMap,
  const occ::handle<StepRepr_RepresentationContext>&                                   theRC,
  const StepData_Factors& theLocalFactors)
{
  // Get working data
  const occ::handle<Interface_InterfaceModel>& aModel  = theWS->Model();
  const occ::handle<Interface_HGraph>          aHGraph = theWS->HGraph();
  if (aHGraph.IsNull())
  {
    return nullptr;
  }
  Interface_Graph                    aGraph = aHGraph->Graph();
  occ::handle<XCAFDoc_GeomTolerance> aGTAttr;
  if (!theGeomTolL.FindAttribute(XCAFDoc_GeomTolerance::GetID(), aGTAttr))
  {
    return nullptr;
  }
  occ::handle<XCAFDimTolObjects_GeomToleranceObject> anObject = aGTAttr->GetObject();
  if (anObject.IsNull())
  {
    return nullptr;
  }

  // Unit
  StepBasic_Unit aUnit = GetUnit(theRC);

  NCollection_Sequence<occ::handle<XCAFDimTolObjects_DatumObject>> aDatums;
  int                                                              aMaxDatumNum = 0;

  for (NCollection_Sequence<TDF_Label>::Iterator aDatumIter(theDatumSeq); aDatumIter.More();
       aDatumIter.Next())
  {
    occ::handle<XCAFDoc_Datum> aDatumAttr;
    if (!aDatumIter.Value().FindAttribute(XCAFDoc_Datum::GetID(), aDatumAttr))
    {
      continue;
    }
    occ::handle<XCAFDimTolObjects_DatumObject> aDatumObj = aDatumAttr->GetObject();
    if (aDatumObj.IsNull())
    {
      continue;
    }
    aDatums.Append(aDatumObj);
    aMaxDatumNum = std::max(aMaxDatumNum, aDatumObj->GetPosition());
  }
  if (aMaxDatumNum == 0)
  {
    return nullptr;
  }

  // The datums of the tolerance's own part: a letter is one part's
  const TCollection_AsciiString aPartKey = datumPartKey(theGeomTolL);

  occ::handle<NCollection_HArray1<occ::handle<StepDimTol_DatumReferenceCompartment>>>
    aConstituents =
      new NCollection_HArray1<occ::handle<StepDimTol_DatumReferenceCompartment>>(1, aMaxDatumNum);
  // Auxiliary datum to initialize attributes in Datum_System
  occ::handle<StepDimTol_Datum> aFirstDatum;
  for (int aConstituentsNum = 1; aConstituentsNum <= aMaxDatumNum; aConstituentsNum++)
  {
    // Collect datums with i-th position
    NCollection_Sequence<occ::handle<XCAFDimTolObjects_DatumObject>> aDatumSeqPos;
    for (NCollection_Sequence<occ::handle<XCAFDimTolObjects_DatumObject>>::Iterator aDatumIter(
           aDatums);
         aDatumIter.More();
         aDatumIter.Next())
    {
      if (aDatumIter.Value()->GetPosition() == aConstituentsNum)
      {
        aDatumSeqPos.Append(aDatumIter.Value());
      }
    }
    if (aDatumSeqPos.Length() < 1)
    {
      continue;
    }
    // Initialize Datum_Reference_Compartment
    StepDimTol_DatumOrCommonDatum                     aDatumRef;
    occ::handle<StepDimTol_DatumReferenceCompartment> aCompartment =
      new StepDimTol_DatumReferenceCompartment();
    occ::handle<NCollection_HArray1<StepDimTol_DatumReferenceModifier>> aModifiers;
    if (aDatumSeqPos.Length() == 1)
    {
      // Datum entity
      occ::handle<Standard_Transient> aFDValue;
      if ((theDatumMap.Find(aDatumSeqPos.Value(1)->GetName()->String().Cat("@").Cat(aPartKey),
                            aFDValue)
           || theDatumMap.Find(aDatumSeqPos.Value(1)->GetName()->String(), aFDValue))
          && !aFDValue.IsNull())
      {
        aFirstDatum = occ::down_cast<StepDimTol_Datum>(aFDValue);
      }
      aDatumRef.SetValue(aFirstDatum);
      // Modifiers
      NCollection_Sequence<XCAFDimTolObjects_DatumSingleModif> aSimpleModifiers =
        aDatumSeqPos.Value(1)->GetModifiers();
      XCAFDimTolObjects_DatumModifWithValue aModifWithVal;
      double                                aValue = 0;
      aDatumSeqPos.Value(1)->GetModifierWithValue(aModifWithVal, aValue);
      aValue /= theLocalFactors.LengthFactor();
      aModifiers = STEPCAFControl_GDTProperty::GetDatumRefModifiers(aSimpleModifiers,
                                                                    aModifWithVal,
                                                                    aValue,
                                                                    aUnit);
      // Add Datum_Reference_Modifier_With_Value
      if (!aModifiers.IsNull())
      {
        occ::handle<StepDimTol_DatumReferenceModifierWithValue> aDRMWV =
          aModifiers->Value(aModifiers->Length()).DatumReferenceModifierWithValue();
        if (!aDRMWV.IsNull())
        {
          aModel->AddWithRefs(aDRMWV);
        }
      }
    }
    else
    {
      int aSetInd = 1;
      occ::handle<NCollection_HArray1<occ::handle<StepDimTol_DatumReferenceElement>>>
        aCommonDatumList = new NCollection_HArray1<occ::handle<StepDimTol_DatumReferenceElement>>(
          1,
          aDatumSeqPos.Length());
      for (NCollection_Sequence<occ::handle<XCAFDimTolObjects_DatumObject>>::Iterator aDatumIter(
             aDatumSeqPos);
           aDatumIter.More();
           aDatumIter.Next())
      {
        // Datum entity
        const occ::handle<XCAFDimTolObjects_DatumObject>& aDatumObj = aDatumIter.Value();
        occ::handle<StepDimTol_Datum>                     aDatum;
        occ::handle<Standard_Transient>                   aDValue;
        if (theDatumMap.Find(aDatumObj->GetName()->String().Cat("@").Cat(aPartKey), aDValue)
            || theDatumMap.Find(aDatumObj->GetName()->String(), aDValue))
        {
          aDatum = occ::down_cast<StepDimTol_Datum>(aDValue);
        }
        StepDimTol_DatumOrCommonDatum anElemDatumRef;
        anElemDatumRef.SetValue(aDatum);
        if (aFirstDatum.IsNull())
        {
          aFirstDatum = aDatum;
        }
        // Modifiers
        NCollection_Sequence<XCAFDimTolObjects_DatumSingleModif> aSimpleModifiers =
          aDatumObj->GetModifiers();
        XCAFDimTolObjects_DatumModifWithValue aModifWithVal;
        double                                aValue = 0;
        aDatumObj->GetModifierWithValue(aModifWithVal, aValue);
        aValue /= theLocalFactors.LengthFactor();
        occ::handle<NCollection_HArray1<StepDimTol_DatumReferenceModifier>> anElemModifiers =
          STEPCAFControl_GDTProperty::GetDatumRefModifiers(aSimpleModifiers,
                                                           aModifWithVal,
                                                           aValue,
                                                           aUnit);
        // Add Datum_Reference_Modifier_With_Value
        if (!anElemModifiers.IsNull())
        {
          occ::handle<StepDimTol_DatumReferenceModifierWithValue> aDRMWV =
            anElemModifiers->Value(anElemModifiers->Length()).DatumReferenceModifierWithValue();
          if (!aDRMWV.IsNull())
          {
            aModel->AddWithRefs(aDRMWV);
          }
        }
        // Datum_Reference_Element
        occ::handle<StepDimTol_DatumReferenceElement> anElement =
          new StepDimTol_DatumReferenceElement();
        anElement->Init(new TCollection_HAsciiString(),
                        new TCollection_HAsciiString(),
                        aDatum->OfShape(),
                        aDatum->ProductDefinitional(),
                        anElemDatumRef,
                        !anElemModifiers.IsNull(),
                        anElemModifiers);
        aModel->AddWithRefs(anElement);
        aCommonDatumList->SetValue(aSetInd++, anElement);
      }
      aDatumRef.SetValue(aCommonDatumList);
    }
    aCompartment->Init(new TCollection_HAsciiString(),
                       new TCollection_HAsciiString(),
                       aFirstDatum->OfShape(),
                       aFirstDatum->ProductDefinitional(),
                       aDatumRef,
                       !aModifiers.IsNull(),
                       aModifiers);
    aModel->AddWithRefs(aCompartment);
    aConstituents->SetValue(aConstituentsNum, aCompartment);
  }
  // Remove null elements from aConstituents
  int aNbConstituents = 0;
  for (NCollection_HArray1<occ::handle<StepDimTol_DatumReferenceCompartment>>::Iterator
         aConstituentIter(aConstituents->Array1());
       aConstituentIter.More();
       aConstituentIter.Next())
  {
    if (!aConstituentIter.Value().IsNull())
    {
      aNbConstituents++;
    }
  }
  occ::handle<NCollection_HArray1<occ::handle<StepDimTol_DatumReferenceCompartment>>>
    aResConstituents =
      new NCollection_HArray1<occ::handle<StepDimTol_DatumReferenceCompartment>>(1,
                                                                                 aNbConstituents);
  int aConstituentsIt = 0;
  for (NCollection_HArray1<occ::handle<StepDimTol_DatumReferenceCompartment>>::Iterator
         aConstituentIter(aConstituents->Array1());
       aConstituentIter.More();
       aConstituentIter.Next())
  {
    if (!aConstituentIter.Value().IsNull())
    {
      aConstituentsIt++;
      aResConstituents->SetValue(aConstituentsIt, aConstituentIter.Value());
    }
  }

  occ::handle<NCollection_HArray1<StepDimTol_DatumSystemOrReference>> aDatumSystem;
  occ::handle<StepDimTol_DatumSystem> aDS = new StepDimTol_DatumSystem();
  aDS->Init(new TCollection_HAsciiString(),
            new TCollection_HAsciiString(),
            aFirstDatum->OfShape(),
            aFirstDatum->ProductDefinitional(),
            aResConstituents);
  aModel->AddWithRefs(aDS);
  StepDimTol_DatumSystemOrReference anArrayValue;
  anArrayValue.SetValue(aDS);
  aDatumSystem = new NCollection_HArray1<StepDimTol_DatumSystemOrReference>(1, 1);
  aDatumSystem->SetValue(1, anArrayValue);

  // Axis
  if (anObject->HasAxis())
  {
    GeomToStep_MakeAxis2Placement3d        anAxisMaker(anObject->GetAxis(), theLocalFactors);
    occ::handle<StepGeom_Axis2Placement3d> anAxis = anAxisMaker.Value();
    anAxis->SetName(new TCollection_HAsciiString("orientation"));
    occ::handle<StepAP242_GeometricItemSpecificUsage> aGISU =
      new StepAP242_GeometricItemSpecificUsage();
    StepAP242_ItemIdentifiedRepresentationUsageDefinition aDefinition;
    aDefinition.SetValue(aDS);
    occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> anReprItems =
      new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, 1);
    occ::handle<StepRepr_RepresentationItem> anIdentifiedItem = anAxis;
    anReprItems->SetValue(1, anIdentifiedItem);
    occ::handle<StepShape_ShapeDefinitionRepresentation> aSDR;
    for (Interface_EntityIterator aSharingsIter = aGraph.Sharings(aFirstDatum->OfShape());
         aSharingsIter.More() && aSDR.IsNull();
         aSharingsIter.Next())
    {
      aSDR = occ::down_cast<StepShape_ShapeDefinitionRepresentation>(aSharingsIter.Value());
    }
    if (aSDR.IsNull())
    {
      return aDatumSystem;
    }

    aGISU->Init(new TCollection_HAsciiString(),
                new TCollection_HAsciiString(),
                aDefinition,
                aSDR->UsedRepresentation(),
                anReprItems);
    aModel->AddWithRefs(anAxis);
    aModel->AddWithRefs(aGISU);
  }

  return aDatumSystem;
}

//=======================================================================
// function : writeToleranceZone
// purpose  : auxiliary (write tolerace zones)
//=======================================================================
void STEPCAFControl_Writer::writeToleranceZone(
  const occ::handle<XSControl_WorkSession>&                 theWS,
  const occ::handle<XCAFDimTolObjects_GeomToleranceObject>& theObject,
  const occ::handle<StepDimTol_GeometricTolerance>&         theEntity,
  const occ::handle<StepRepr_RepresentationContext>&        theRC,
  const StepData_Factors&                                   theLocalFactors)
{
  // Get working data
  const occ::handle<Interface_InterfaceModel>& aModel = theWS->Model();
  if (theEntity.IsNull() || theObject.IsNull())
  {
    return;
  }

  // Return if there is no tolerance zones
  const XCAFDimTolObjects_GeomToleranceZoneModif aZoneModif = theObject->GetZoneModifier();
  if (theObject->GetTypeOfValue() == XCAFDimTolObjects_GeomToleranceTypeValue_None
      && aZoneModif != XCAFDimTolObjects_GeomToleranceZoneModif_Runout
      && aZoneModif != XCAFDimTolObjects_GeomToleranceZoneModif_Projected)
  {
    return;
  }

  // Create Tolerance_Zone
  occ::handle<StepDimTol_ToleranceZoneForm> aForm = new StepDimTol_ToleranceZoneForm();
  aModel->AddWithRefs(aForm);
  aForm->Init(STEPCAFControl_GDTProperty::GetTolValueType(theObject->GetTypeOfValue()));
  occ::handle<NCollection_HArray1<StepDimTol_ToleranceZoneTarget>> aZoneTargetArray =
    new NCollection_HArray1<StepDimTol_ToleranceZoneTarget>(1, 1);
  StepDimTol_ToleranceZoneTarget aTarget;
  aTarget.SetValue(theEntity);
  aZoneTargetArray->SetValue(1, aTarget);
  occ::handle<StepDimTol_ToleranceZone> aZone = new StepDimTol_ToleranceZone();
  aZone->Init(new TCollection_HAsciiString(),
              new TCollection_HAsciiString(),
              theEntity->TolerancedShapeAspect().ShapeAspect().IsNull()
                ? theEntity->TolerancedShapeAspect().ProductDefinitionShape()
                : theEntity->TolerancedShapeAspect().ShapeAspect()->OfShape(),
              StepData_LFalse,
              aZoneTargetArray,
              aForm);
  aModel->AddWithRefs(aZone);

  if (aZoneModif == XCAFDimTolObjects_GeomToleranceZoneModif_Projected)
  {
    // Projected_Zone_Definition: projected from the part, by the length the
    // modifier has (as the NIST files write it: no boundaries, the end a
    // shape aspect of the part)
    occ::handle<StepBasic_LengthMeasureWithUnit> aLMWU = new StepBasic_LengthMeasureWithUnit();
    occ::handle<StepBasic_MeasureValueMember>    aLength = new StepBasic_MeasureValueMember();
    aLength->SetName("LENGTH_MEASURE");
    aLength->SetReal(theObject->GetValueOfZoneModifier() / theLocalFactors.LengthFactor());
    aLMWU->Init(aLength, GetUnit(theRC));
    occ::handle<StepRepr_ShapeAspect> anEnd = new StepRepr_ShapeAspect();
    anEnd->Init(new TCollection_HAsciiString(),
                new TCollection_HAsciiString(),
                aZone->OfShape(),
                StepData_LFalse);
    occ::handle<StepDimTol_ProjectedZoneDefinition> aDefinition =
      new StepDimTol_ProjectedZoneDefinition();
    aDefinition->Init(aZone, nullptr, anEnd, aLMWU);
    aModel->AddWithRefs(aDefinition);
    return;
  }
  if (aZoneModif != XCAFDimTolObjects_GeomToleranceZoneModif_Runout)
  {
    // The zone alone: no definition of it but a runout's or a projected one's
    return;
  }

  // Runout_Tolerance_Zone
  occ::handle<StepBasic_PlaneAngleMeasureWithUnit> aPAMWU =
    new StepBasic_PlaneAngleMeasureWithUnit();
  occ::handle<StepBasic_MeasureValueMember> aValueMember = new StepBasic_MeasureValueMember();
  aValueMember->SetName("PLANE_ANGLE_MEASURE");
  aValueMember->SetReal(theObject->GetValueOfZoneModifier());
  aPAMWU->Init(aValueMember, GetUnit(theRC, true));
  occ::handle<StepDimTol_RunoutZoneOrientation> anOrientation =
    new StepDimTol_RunoutZoneOrientation();
  anOrientation->Init(aPAMWU);
  occ::handle<StepDimTol_RunoutZoneDefinition> aDefinition = new StepDimTol_RunoutZoneDefinition();
  aDefinition->Init(aZone, nullptr, anOrientation);
  aModel->AddWithRefs(aDefinition);
  aModel->AddWithRefs(anOrientation);
  aModel->AddWithRefs(aPAMWU);
}

//=======================================================================
// function : presentedType
// purpose  : the name AP242 PMI recommended practice 8.4 gives the graphic
//           of a tolerance or a dimension of the type
//=======================================================================
static const char* presentedType(const XCAFDimTolObjects_GeomToleranceType theType)
{
  switch (theType)
  {
    case XCAFDimTolObjects_GeomToleranceType_Angularity:
      return "angularity";
    case XCAFDimTolObjects_GeomToleranceType_CircularRunout:
      return "circular runout";
    case XCAFDimTolObjects_GeomToleranceType_CircularityOrRoundness:
      return "circularity";
    case XCAFDimTolObjects_GeomToleranceType_Coaxiality:
      return "coaxiality";
    case XCAFDimTolObjects_GeomToleranceType_Concentricity:
      return "concentricity";
    case XCAFDimTolObjects_GeomToleranceType_Cylindricity:
      return "cylindricity";
    case XCAFDimTolObjects_GeomToleranceType_Flatness:
      return "flatness";
    case XCAFDimTolObjects_GeomToleranceType_Parallelism:
      return "parallelism";
    case XCAFDimTolObjects_GeomToleranceType_Perpendicularity:
      return "perpendicularity";
    case XCAFDimTolObjects_GeomToleranceType_Position:
      return "position";
    case XCAFDimTolObjects_GeomToleranceType_ProfileOfLine:
      return "profile of line";
    case XCAFDimTolObjects_GeomToleranceType_ProfileOfSurface:
      return "profile of surface";
    case XCAFDimTolObjects_GeomToleranceType_Straightness:
      return "straightness";
    case XCAFDimTolObjects_GeomToleranceType_Symmetry:
      return "symmetry";
    case XCAFDimTolObjects_GeomToleranceType_TotalRunout:
      return "total runout";
    default:
      return "general tolerance";
  }
}

static const char* presentedType(const XCAFDimTolObjects_DimensionType theType)
{
  switch (theType)
  {
    case XCAFDimTolObjects_DimensionType_Location_Angular:
    case XCAFDimTolObjects_DimensionType_Size_Angular:
      return "angular dimension";
    case XCAFDimTolObjects_DimensionType_Size_Diameter:
    case XCAFDimTolObjects_DimensionType_Size_SphericalDiameter:
    case XCAFDimTolObjects_DimensionType_Size_ToroidalHighMajorDiameter:
    case XCAFDimTolObjects_DimensionType_Size_ToroidalLowMajorDiameter:
    case XCAFDimTolObjects_DimensionType_Size_ToroidalMajorDiameter:
    case XCAFDimTolObjects_DimensionType_Size_ToroidalMinorDiameter:
      return "diameter dimension";
    case XCAFDimTolObjects_DimensionType_Size_Radius:
    case XCAFDimTolObjects_DimensionType_Size_SphericalRadius:
    case XCAFDimTolObjects_DimensionType_Size_ToroidalHighMajorRadius:
    case XCAFDimTolObjects_DimensionType_Size_ToroidalLowMajorRadius:
    case XCAFDimTolObjects_DimensionType_Size_ToroidalMajorRadius:
    case XCAFDimTolObjects_DimensionType_Size_ToroidalMinorRadius:
      return "radial dimension";
    case XCAFDimTolObjects_DimensionType_Location_CurvedDistance:
    case XCAFDimTolObjects_DimensionType_Location_WithPath:
    case XCAFDimTolObjects_DimensionType_Size_CurveLength:
    case XCAFDimTolObjects_DimensionType_Size_WithPath:
      return "curve dimension";
    case XCAFDimTolObjects_DimensionType_CommonLabel:
      return "note";
    case XCAFDimTolObjects_DimensionType_DimensionPresentation:
      return "general dimension";
    default:
      return "linear dimension";
  }
}

//=======================================================================
// function : toleranceWithDatumReference
// purpose  : the entity of a tolerance type that is a
//           geometric_tolerance_with_datum_reference, or null
//=======================================================================
static occ::handle<StepDimTol_GeometricToleranceWithDatumReference> toleranceWithDatumReference(
  const XCAFDimTolObjects_GeomToleranceType theType)
{
  switch (theType)
  {
    case XCAFDimTolObjects_GeomToleranceType_Angularity:
      return new StepDimTol_AngularityTolerance();
    case XCAFDimTolObjects_GeomToleranceType_CircularRunout:
      return new StepDimTol_CircularRunoutTolerance();
    case XCAFDimTolObjects_GeomToleranceType_Coaxiality:
      return new StepDimTol_CoaxialityTolerance();
    case XCAFDimTolObjects_GeomToleranceType_Concentricity:
      return new StepDimTol_ConcentricityTolerance();
    case XCAFDimTolObjects_GeomToleranceType_Parallelism:
      return new StepDimTol_ParallelismTolerance();
    case XCAFDimTolObjects_GeomToleranceType_Perpendicularity:
      return new StepDimTol_PerpendicularityTolerance();
    case XCAFDimTolObjects_GeomToleranceType_Symmetry:
      return new StepDimTol_SymmetryTolerance();
    case XCAFDimTolObjects_GeomToleranceType_TotalRunout:
      return new StepDimTol_TotalRunoutTolerance();
    default:
      return nullptr;
  }
}

//=======================================================================
// function : writeGeomTolerance
// purpose  : auxiliary (write Geometric_Tolerance entity for given shapes,
//           label and datum system)
//======================================================================
occ::handle<StepDimTol_GeometricTolerance> STEPCAFControl_Writer::writeGeomTolerance(
  const occ::handle<XSControl_WorkSession>&                                  theWS,
  const NCollection_Sequence<TDF_Label>&                                     theShapeSeqL,
  const TDF_Label&                                                           theGeomTolL,
  const occ::handle<NCollection_HArray1<StepDimTol_DatumSystemOrReference>>& theDatumSystem,
  const occ::handle<StepRepr_RepresentationContext>&                         theRC,
  const StepData_Factors&                                                    theLocalFactors)
{
  // Get working data
  const occ::handle<Interface_InterfaceModel>& aModel = theWS->Model();
  occ::handle<XCAFDoc_GeomTolerance>           aGTAttr;
  if (!theGeomTolL.FindAttribute(XCAFDoc_GeomTolerance::GetID(), aGTAttr))
  {
    return nullptr;
  }
  occ::handle<XCAFDimTolObjects_GeomToleranceObject> anObject = aGTAttr->GetObject();
  if (anObject.IsNull())
  {
    return nullptr;
  }

  // Value
  occ::handle<StepBasic_LengthMeasureWithUnit> aLMWU        = new StepBasic_LengthMeasureWithUnit();
  StepBasic_Unit                               aUnit        = GetUnit(theRC);
  occ::handle<StepBasic_MeasureValueMember>    aValueMember = new StepBasic_MeasureValueMember();
  aValueMember->SetName("LENGTH_MEASURE");
  aValueMember->SetReal(anObject->GetValue() / theLocalFactors.LengthFactor());
  aLMWU->Init(aValueMember, aUnit);
  aModel->AddWithRefs(aLMWU);

  // Geometric_Tolerance target
  // All around: the target is an all_around_shape_aspect of the faces' aspects, as it is
  // read. The modifier was read and never written: every profile all around (seven of
  // NIST's files) came back one of its faces alone
  bool                                                       isAllAround = false;
  bool                                                       isAllOver   = false;
  const NCollection_Sequence<XCAFDimTolObjects_GeomToleranceModif> anAllModifs =
    anObject->GetModifiers();
  for (NCollection_Sequence<XCAFDimTolObjects_GeomToleranceModif>::Iterator aModIt(anAllModifs);
       aModIt.More();
       aModIt.Next())
  {
    isAllAround = isAllAround || aModIt.Value() == XCAFDimTolObjects_GeomToleranceModif_All_Around;
    isAllOver   = isAllOver || aModIt.Value() == XCAFDimTolObjects_GeomToleranceModif_All_Over;
  }
  occ::handle<StepRepr_ShapeAspect>                 aMainSA;
  occ::handle<StepRepr_RepresentationContext>       dummyRC;
  occ::handle<StepAP242_GeometricItemSpecificUsage> dummyGISU;
  // All over, on the part's own shape: the target is the product's shape itself, no aspect of
  // it and no geometry (as NIST's general profile is, FTC-07's and FTC-10's)
  occ::handle<StepRepr_ProductDefinitionShape> aWholePDS;
  if (isAllOver && theShapeSeqL.Length() == 1
      && !XCAFDoc_ShapeTool::IsSubShape(theShapeSeqL.Value(1)))
  {
    const occ::handle<Interface_HGraph> aHGraph = theWS->HGraph();
    TopLoc_Location                                       aLoc;
    NCollection_Sequence<occ::handle<Standard_Transient>> aSeqRI;
    FindEntities(theWS->TransferWriter()->FinderProcess(),
                 XCAFDoc_ShapeTool::GetShape(theShapeSeqL.Value(1)),
                 aLoc,
                 aSeqRI);
    if (!aHGraph.IsNull() && !aSeqRI.IsEmpty())
    {
      occ::handle<StepRepr_RepresentationContext> aRC;
      aWholePDS = FindPDS(aHGraph->Graph(), aSeqRI.Value(1), aRC);
    }
  }
  if (!aWholePDS.IsNull())
  {
  }
  else if (theShapeSeqL.Length() == 1 && !isAllAround)
  {
    TopoDS_Shape aShape = XCAFDoc_ShapeTool::GetShape(theShapeSeqL.Value(1));
    aMainSA             = writeShapeAspect(theWS, theGeomTolL, aShape, dummyRC, dummyGISU);
    aModel->AddWithRefs(aMainSA);
  }
  else
  {
    occ::handle<StepRepr_CompositeShapeAspect> aCSA;
    for (NCollection_Sequence<TDF_Label>::Iterator aShIter(theShapeSeqL); aShIter.More();
         aShIter.Next())
    {
      TopoDS_Shape                      aShape = XCAFDoc_ShapeTool::GetShape(aShIter.Value());
      occ::handle<StepRepr_ShapeAspect> aSA =
        writeShapeAspect(theWS, theGeomTolL, aShape, dummyRC, dummyGISU);
      if (aSA.IsNull())
      {
        continue;
      }
      if (aCSA.IsNull())
      {
        if (isAllAround)
        {
          aCSA = new StepRepr_AllAroundShapeAspect();
        }
        else
        {
          aCSA = new StepRepr_CompositeShapeAspect();
        }
        aCSA->Init(aSA->Name(), aSA->Description(), aSA->OfShape(), aSA->ProductDefinitional());
        aModel->AddWithRefs(aCSA);
      }
      occ::handle<StepRepr_ShapeAspectRelationship> aSAR = new StepRepr_ShapeAspectRelationship();
      aSAR->Init(new TCollection_HAsciiString(), false, nullptr, aCSA, aSA);
      aModel->AddWithRefs(aSAR);
    }
    aMainSA = aCSA;
  }
  StepDimTol_GeometricToleranceTarget aGTTarget;
  if (!aWholePDS.IsNull())
  {
    aGTTarget.SetValue(aWholePDS);
  }
  else
  {
    aGTTarget.SetValue(aMainSA);
  }

  bool isWithModif = false, isWithDatRef = false, isWithMaxTol = false;
  // Modifiers
  // Simple modifiers
  NCollection_Sequence<XCAFDimTolObjects_GeomToleranceModif> aModifiers = anObject->GetModifiers();
  occ::handle<NCollection_HArray1<StepDimTol_GeometricToleranceModifier>> aModifArray;
  occ::handle<StepBasic_LengthMeasureWithUnit>                            aMaxLMWU;
  int aModifNb = aModifiers.Length();
  if (anObject->GetMaterialRequirementModifier() != XCAFDimTolObjects_GeomToleranceMatReqModif_None)
  {
    aModifNb++;
  }
  for (int i = 1; i <= aModifiers.Length(); i++)
  {
    if (aModifiers.Value(i) == XCAFDimTolObjects_GeomToleranceModif_All_Around
        || aModifiers.Value(i) == XCAFDimTolObjects_GeomToleranceModif_All_Over)
    {
      aModifNb--;
    }
  }
  if (aModifNb > 0)
  {
    isWithModif = true;
    aModifArray = new NCollection_HArray1<StepDimTol_GeometricToleranceModifier>(1, aModifNb);
    int k       = 1;
    for (int i = 1; i <= aModifiers.Length(); i++)
    {
      if (aModifiers.Value(i) == XCAFDimTolObjects_GeomToleranceModif_All_Around
          || aModifiers.Value(i) == XCAFDimTolObjects_GeomToleranceModif_All_Over)
      {
        continue;
      }
      StepDimTol_GeometricToleranceModifier aModif =
        STEPCAFControl_GDTProperty::GetGeomToleranceModifier(aModifiers.Value(i));
      aModifArray->SetValue(k, aModif);
      k++;
    }
    if (anObject->GetMaterialRequirementModifier() == XCAFDimTolObjects_GeomToleranceMatReqModif_L)
    {
      aModifArray->SetValue(aModifNb, StepDimTol_GTMLeastMaterialRequirement);
    }
    else if (anObject->GetMaterialRequirementModifier()
             == XCAFDimTolObjects_GeomToleranceMatReqModif_M)
    {
      aModifArray->SetValue(aModifNb, StepDimTol_GTMMaximumMaterialRequirement);
    }
    // Modifier with value
    if (anObject->GetMaxValueModifier() != 0)
    {
      isWithMaxTol = true;
      aMaxLMWU     = new StepBasic_LengthMeasureWithUnit();
      occ::handle<StepBasic_MeasureValueMember> aModifierValueMember =
        new StepBasic_MeasureValueMember();
      aModifierValueMember->SetName("LENGTH_MEASURE");
      aModifierValueMember->SetReal(anObject->GetMaxValueModifier()
                                    / theLocalFactors.LengthFactor());
      aMaxLMWU->Init(aModifierValueMember, aUnit);
      aModel->AddWithRefs(aMaxLMWU);
    }
  }

  // Datum Reference
  isWithDatRef = !theDatumSystem.IsNull();

  // Collect all attributes
  occ::handle<TCollection_HAsciiString> aName        = new TCollection_HAsciiString(),
                                        aDescription = new TCollection_HAsciiString();
  occ::handle<StepDimTol_GeometricToleranceWithDatumReference> aGTWDR =
    new StepDimTol_GeometricToleranceWithDatumReference();
  aGTWDR->SetDatumSystem(theDatumSystem);
  occ::handle<StepDimTol_GeometricToleranceWithModifiers> aGTWM =
    new StepDimTol_GeometricToleranceWithModifiers();
  aGTWM->SetModifiers(aModifArray);
  StepDimTol_GeometricToleranceType aType =
    STEPCAFControl_GDTProperty::GetGeomToleranceType(anObject->GetType());

  // Init and write necessary subtype of Geometric_Tolerance entity
  occ::handle<StepDimTol_GeometricTolerance> aGeomTol;
  const XCAFDimTolObjects_GeomToleranceUnitBasis aUnitBasis = anObject->GetUnitBasis();
  if (aUnitBasis != XCAFDimTolObjects_GeomToleranceUnitBasis_None && anObject->GetUnitSize() > 0.
      && !isWithModif && !isWithDatRef)
  {
    // Per unit length or area: Geometric_Tolerance & Geometric_Tolerance_With_Defined_Unit
    // (& Geometric_Tolerance_With_Defined_Area_Unit) & its type. With modifiers or datums it is
    // a complex no reader knows: the value is written for the whole feature, as before.
    const auto aLength = [&](const double theValue) {
      occ::handle<StepBasic_LengthMeasureWithUnit> aMeasure = new StepBasic_LengthMeasureWithUnit();
      occ::handle<StepBasic_MeasureValueMember>    aMember  = new StepBasic_MeasureValueMember();
      aMember->SetName("LENGTH_MEASURE");
      aMember->SetReal(theValue / theLocalFactors.LengthFactor());
      aMeasure->Init(aMember, aUnit);
      aModel->AddWithRefs(aMeasure);
      return aMeasure;
    };
    occ::handle<StepDimTol_GeoTolAndGeoTolWthDefUnit> aResult =
      new StepDimTol_GeoTolAndGeoTolWthDefUnit();
    aResult->Init(aName, aDescription, aLMWU, aGTTarget);
    aResult->SetToleranceType(aType);
    aResult->SetUnitSize(aLength(anObject->GetUnitSize()));
    if (aUnitBasis != XCAFDimTolObjects_GeomToleranceUnitBasis_Length)
    {
      aResult->SetPerArea(true);
      aResult->SetAreaType(aUnitBasis == XCAFDimTolObjects_GeomToleranceUnitBasis_Circle
                             ? StepDimTol_Circular
                           : aUnitBasis == XCAFDimTolObjects_GeomToleranceUnitBasis_Rectangle
                             ? StepDimTol_Rectangular
                             : StepDimTol_Square);
      if (aUnitBasis == XCAFDimTolObjects_GeomToleranceUnitBasis_Rectangle
          && anObject->GetSecondUnitSize() > 0.)
      {
        aResult->SetSecondUnitSize(aLength(anObject->GetSecondUnitSize()));
      }
    }
    aGeomTol = aResult;
  }
  else if (isWithModif)
  {
    if (isWithMaxTol)
    {
      if (isWithDatRef)
      {
        // Geometric_Tolerance & Geometric_Tolerance_With_Datum_Reference &
        // Geometric_Tolerance_With_Maximum_Tolerance & Geometric_Tolerance_With_Modifiers
        occ::handle<StepDimTol_GeoTolAndGeoTolWthDatRefAndGeoTolWthMaxTol> aResult =
          new StepDimTol_GeoTolAndGeoTolWthDatRefAndGeoTolWthMaxTol();
        aResult->Init(aName, aDescription, aLMWU, aGTTarget, aGTWDR, aGTWM, aMaxLMWU, aType);
        aGeomTol = aResult;
      }
      else
      {
        // Geometric_Tolerance & Geometric_Tolerance_With_Maximum_Tolerance &
        // Geometric_Tolerance_With_Modifiers
        occ::handle<StepDimTol_GeoTolAndGeoTolWthMaxTol> aResult =
          new StepDimTol_GeoTolAndGeoTolWthMaxTol();
        aResult->Init(aName, aDescription, aLMWU, aGTTarget, aGTWM, aMaxLMWU, aType);
        aGeomTol = aResult;
      }
    }
    else
    {
      if (isWithDatRef)
      {
        // Geometric_Tolerance & Geometric_Tolerance_With_Datum_Reference &
        // Geometric_Tolerance_With_Modifiers
        occ::handle<StepDimTol_GeoTolAndGeoTolWthDatRefAndGeoTolWthMod> aResult =
          new StepDimTol_GeoTolAndGeoTolWthDatRefAndGeoTolWthMod();
        aResult->Init(aName, aDescription, aLMWU, aGTTarget, aGTWDR, aGTWM, aType);
        aGeomTol = aResult;
      }
      else
      {
        // Geometric_Tolerance & Geometric_Tolerance_With_Modifiers
        occ::handle<StepDimTol_GeoTolAndGeoTolWthMod> aResult =
          new StepDimTol_GeoTolAndGeoTolWthMod();
        aResult->Init(aName, aDescription, aLMWU, aGTTarget, aGTWM, aType);
        aGeomTol = aResult;
      }
    }
  }
  else
  {
    occ::handle<StepDimTol_GeometricToleranceWithDatumReference> aWithDatRef =
      isWithDatRef ? toleranceWithDatumReference(anObject->GetType()) : nullptr;
    if (!aWithDatRef.IsNull())
    {
      // A type that is a Geometric_Tolerance_With_Datum_Reference itself
      // (Perpendicularity_Tolerance ...): the entity alone, as the schema
      // has it, not a complex with its own supertype
      aWithDatRef->Init(aName, aDescription, aLMWU, aGTTarget, theDatumSystem);
      aGeomTol = aWithDatRef;
    }
    else if (isWithDatRef)
    {
      // Geometric_Tolerance & Geometric_Tolerance_With_Datum_Reference
      occ::handle<StepDimTol_GeoTolAndGeoTolWthDatRef> aResult =
        new StepDimTol_GeoTolAndGeoTolWthDatRef();
      aResult->Init(aName, aDescription, aLMWU, aGTTarget, aGTWDR, aType);
      aGeomTol = aResult;
    }
    else
    {
      // Geometric_Tolerance
      occ::handle<StepDimTol_GeometricTolerance> aResult =
        STEPCAFControl_GDTProperty::GetGeomTolerance(anObject->GetType());
      if (!aResult.IsNull())
      {
        aResult->Init(aName, aDescription, aLMWU, aGTTarget);
        aGeomTol = aResult;
      }
    }
  }
  aModel->AddWithRefs(aGeomTol);
  writeToleranceZone(theWS, anObject, aGeomTol, theRC, theLocalFactors);
  // Annotation plane and Presentation
  writePresentation(theWS,
                    theGeomTolL,
                    anObject->GetPresentation(),
                    anObject->GetPresentationName(),
                    true,
                    anObject->HasPlane(),
                    anObject->GetPlane(),
                    anObject->GetPointTextAttach(),
                    aGeomTol,
                    keptOr(anObject->GetPresentationType(), presentedType(anObject->GetType())),
                    shapeAspects(aMainSA),
                    theLocalFactors);
  return aGeomTol;
}

//=================================================================================================

bool STEPCAFControl_Writer::writeDGTs(const occ::handle<XSControl_WorkSession>& theWS,
                                      const NCollection_Sequence<TDF_Label>&    theLabels) const
{

  if (theLabels.IsEmpty())
  {
    return false;
  }

  // get working data
  const occ::handle<Interface_InterfaceModel>& aModel = theWS->Model();
  const occ::handle<XSControl_TransferWriter>& aTW    = theWS->TransferWriter();
  const occ::handle<Transfer_FinderProcess>&   aFP    = aTW->FinderProcess();

  const occ::handle<Interface_HGraph> aHGraph = theWS->HGraph();
  if (aHGraph.IsNull())
  {
    return false;
  }

  Interface_Graph aGraph = aHGraph->Graph();
  NCollection_DataMap<TCollection_AsciiString, occ::handle<Standard_Transient>> aDatumMap;

  NCollection_Sequence<TDF_Label> aDGTLabels;
  // Iterate on requested shapes collect Tools
  for (NCollection_Map<TDF_Label>::Iterator aLabelIter(myRootLabels); aLabelIter.More();
       aLabelIter.Next())
  {
    const TDF_Label&                aLabel   = aLabelIter.Value();
    occ::handle<XCAFDoc_DimTolTool> aDGTTool = XCAFDoc_DocumentTool::DimTolTool(aLabel);
    NCollection_Sequence<TDF_Label> aDGTLS;
    aDGTTool->GetDatumLabels(aDGTLS);
    aDGTLabels.Append(aDGTLS);
  }

  if (aDGTLabels.IsEmpty())
  {
    return false;
  }

  for (NCollection_Sequence<TDF_Label>::Iterator aDGTIter(aDGTLabels); aDGTIter.More();
       aDGTIter.Next())
  {
    const TDF_Label&                aDatumL = aDGTIter.Value();
    NCollection_Sequence<TDF_Label> aShapeL;
    NCollection_Sequence<TDF_Label> aNullSeq;
    if (!XCAFDoc_DimTolTool::GetRefShapeLabel(aDatumL, aShapeL, aNullSeq))
    {
      continue;
    }
    // find target shape
    TopoDS_Shape    aShape = XCAFDoc_ShapeTool::GetShape(aShapeL.Value(1));
    TopLoc_Location aLoc;
    NCollection_Sequence<occ::handle<Standard_Transient>> aSeqRI;
    FindEntities(aFP, aShape, aLoc, aSeqRI);
    if (aSeqRI.IsEmpty())
    {
      Message::SendTrace() << "Warning: Cannot find RI for "
                           << aShape.TShape()->DynamicType()->Name() << "\n";
      continue;
    }
    occ::handle<StepRepr_ProductDefinitionShape> aPDS;
    occ::handle<StepRepr_RepresentationContext>  aRC;
    occ::handle<Standard_Transient>              anEnt = aSeqRI.Value(1);
    occ::handle<StepShape_AdvancedFace>          anAF;
    occ::handle<StepShape_EdgeCurve>             anEC;
    FindPDSforDGT(aGraph, anEnt, aPDS, aRC, anAF, anEC);
    if (aPDS.IsNull())
    {
      continue;
    }
    occ::handle<XCAFDoc_Datum> aDatumAttr;
    if (!aDatumL.FindAttribute(XCAFDoc_Datum::GetID(), aDatumAttr))
    {
      continue;
    }
    occ::handle<TCollection_HAsciiString> aName            = aDatumAttr->GetName();
    occ::handle<TCollection_HAsciiString> anIdentification = aDatumAttr->GetIdentification();
    occ::handle<TCollection_HAsciiString> aDescription     = aDatumAttr->GetDescription();
    if (aDescription.IsNull())
    {
      aDescription = new TCollection_HAsciiString();
    }
    occ::handle<StepDimTol_DatumFeature> aDF    = new StepDimTol_DatumFeature;
    occ::handle<StepDimTol_Datum>        aDatum = new StepDimTol_Datum;
    aDF->Init(aName, new TCollection_HAsciiString, aPDS, StepData_LTrue);
    aModel->AddWithRefs(aDF);
    aDatum->Init(aName, new TCollection_HAsciiString, aPDS, StepData_LFalse, anIdentification);
    aModel->AddWithRefs(aDatum);
    occ::handle<StepRepr_ShapeAspectRelationship> aSAR = new StepRepr_ShapeAspectRelationship;
    aSAR->SetName(aName);
    aSAR->SetRelatingShapeAspect(aDF);
    aSAR->SetRelatedShapeAspect(aDatum);
    aModel->AddWithRefs(aSAR);
    // write chain for DatumFeature
    StepRepr_CharacterizedDefinition aCD;
    aCD.SetValue(aDF);
    occ::handle<StepRepr_PropertyDefinition> aPropD = new StepRepr_PropertyDefinition;
    aPropD->Init(aName, true, aDescription, aCD);
    aModel->AddWithRefs(aPropD);
    StepRepr_RepresentedDefinition aRD;
    aRD.SetValue(aPropD);
    occ::handle<StepShape_ShapeRepresentation> aSR = new StepShape_ShapeRepresentation;
    occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> aHARI =
      new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, 1);
    aHARI->SetValue(1, anAF);
    aSR->Init(aName, aHARI, aRC);
    occ::handle<StepShape_ShapeDefinitionRepresentation> aSDR =
      new StepShape_ShapeDefinitionRepresentation;
    aSDR->Init(aRD, aSR);
    aModel->AddWithRefs(aSDR);
    // write chain for Datum
    StepRepr_CharacterizedDefinition aCD1;
    aCD1.SetValue(aDatum);
    occ::handle<StepRepr_PropertyDefinition> aPropD1 = new StepRepr_PropertyDefinition;
    aPropD1->Init(aName, true, aDescription, aCD1);
    aModel->AddWithRefs(aPropD1);
    StepRepr_RepresentedDefinition aRD1;
    aRD1.SetValue(aPropD1);
    occ::handle<StepShape_ShapeRepresentation> aSR1 = new StepShape_ShapeRepresentation;
    occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> aHARI1 =
      new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, 1);
    aHARI1->SetValue(1, anAF->FaceGeometry());
    aSR1->Init(aName, aHARI1, aRC);
    aModel->AddWithRefs(aSR1);
    occ::handle<StepShape_ShapeDefinitionRepresentation> aSDR1 =
      new StepShape_ShapeDefinitionRepresentation;
    aSDR1->Init(aRD1, aSR1);
    aModel->AddWithRefs(aSDR1);
    // add created Datum into Map
    TCollection_AsciiString aStmp(aName->ToCString());
    aStmp.AssignCat(aDescription->ToCString());
    aStmp.AssignCat(anIdentification->ToCString());
    aDatumMap.Bind(aStmp, aDatum);
  }

  // write Tolerances and Dimensions
  aDGTLabels.Clear();
  // Iterate on requested shapes collect Tools
  for (NCollection_Map<TDF_Label>::Iterator aLabelIter(myRootLabels); aLabelIter.More();
       aLabelIter.Next())
  {
    const TDF_Label&                aLabel   = aLabelIter.Value();
    occ::handle<XCAFDoc_DimTolTool> aDGTTool = XCAFDoc_DocumentTool::DimTolTool(aLabel);
    NCollection_Sequence<TDF_Label> aaDGTLS;
    aDGTTool->GetDimTolLabels(aDGTLabels);
    aDGTLabels.Append(aaDGTLS);
  }

  if (aDGTLabels.IsEmpty())
  {
    return false;
  }
  for (NCollection_Sequence<TDF_Label>::Iterator aDGTIter(aDGTLabels); aDGTIter.More();
       aDGTIter.Next())
  {
    const TDF_Label&                aDimTolL = aDGTIter.Value();
    NCollection_Sequence<TDF_Label> aShapeL;
    NCollection_Sequence<TDF_Label> aNullSeq;
    if (!XCAFDoc_DimTolTool::GetRefShapeLabel(aDimTolL, aShapeL, aNullSeq))
    {
      continue;
    }
    // find target shape
    TopoDS_Shape    aShape = XCAFDoc_ShapeTool::GetShape(aShapeL.Value(1));
    TopLoc_Location aLoc;
    NCollection_Sequence<occ::handle<Standard_Transient>> seqRI;
    FindEntities(aFP, aShape, aLoc, seqRI);
    if (seqRI.IsEmpty())
    {
      Message::SendTrace() << "Warning: Cannot find RI for "
                           << aShape.TShape()->DynamicType()->Name() << "\n";
      continue;
    }
    occ::handle<StepRepr_ProductDefinitionShape> aPDS;
    occ::handle<StepRepr_RepresentationContext>  aRC;
    occ::handle<Standard_Transient>              anEnt = seqRI.Value(1);
    occ::handle<StepShape_AdvancedFace>          anAF;
    occ::handle<StepShape_EdgeCurve>             anEC;
    FindPDSforDGT(aGraph, anEnt, aPDS, aRC, anAF, anEC);
    if (aPDS.IsNull())
    {
      continue;
    }

    occ::handle<XCAFDoc_DimTol> aDimTolAttr;
    if (!aDimTolL.FindAttribute(XCAFDoc_DimTol::GetID(), aDimTolAttr))
    {
      continue;
    }
    int                                      aKind        = aDimTolAttr->GetKind();
    occ::handle<NCollection_HArray1<double>> aVal         = aDimTolAttr->GetVal();
    occ::handle<TCollection_HAsciiString>    aName        = aDimTolAttr->GetName();
    occ::handle<TCollection_HAsciiString>    aDescription = aDimTolAttr->GetDescription();

    // common part of writing D&GT entities
    StepRepr_CharacterizedDefinition  aCD;
    occ::handle<StepRepr_ShapeAspect> aSA = new StepRepr_ShapeAspect;
    aSA->Init(aName, new TCollection_HAsciiString, aPDS, StepData_LTrue);
    aModel->AddWithRefs(aSA);
    aCD.SetValue(aSA);
    occ::handle<StepRepr_PropertyDefinition> aPropD = new StepRepr_PropertyDefinition;
    aPropD->Init(aName, true, aDescription, aCD);
    aModel->AddWithRefs(aPropD);
    StepRepr_RepresentedDefinition aRD;
    aRD.SetValue(aPropD);
    occ::handle<StepShape_ShapeRepresentation> aSR = new StepShape_ShapeRepresentation;
    occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> aHARI =
      new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, 1);
    if (aKind < 20)
    {
      aHARI->SetValue(1, anEC);
    }
    else
    {
      aHARI->SetValue(1, anAF);
    }
    aSR->Init(aName, aHARI, aRC);
    occ::handle<StepShape_ShapeDefinitionRepresentation> aSDR =
      new StepShape_ShapeDefinitionRepresentation;
    aSDR->Init(aRD, aSR);
    aModel->AddWithRefs(aSDR);
    // define aUnit for creation LengthMeasureWithUnit (common for all)
    StepBasic_Unit aUnit;
    aUnit = GetUnit(aRC);

    // specific part of writing D&GT entities
    if (aKind < 20)
    { // dimension
      occ::handle<StepShape_DimensionalSize> aDimSize = new StepShape_DimensionalSize;
      aDimSize->Init(aSA, aDescription);
      aModel->AddWithRefs(aDimSize);
      if (aVal->Length() > 1)
      {
        // create MeasureWithUnits
        occ::handle<StepBasic_MeasureValueMember> aMVM1 = new StepBasic_MeasureValueMember;
        aMVM1->SetName("POSITIVE_LENGTH_MEASURE");
        aMVM1->SetReal(aVal->Value(1));
        occ::handle<StepBasic_MeasureWithUnit> aMWU1 = new StepBasic_MeasureWithUnit;
        aMWU1->Init(aMVM1, aUnit);
        occ::handle<StepBasic_MeasureValueMember> aMVM2 = new StepBasic_MeasureValueMember;
        aMVM2->SetName("POSITIVE_LENGTH_MEASURE");
        aMVM2->SetReal(aVal->Value(2));
        occ::handle<StepBasic_MeasureWithUnit> aMWU2 = new StepBasic_MeasureWithUnit;
        aMWU2->Init(aMVM2, aUnit);
        occ::handle<StepRepr_RepresentationItem> aRI1 = new StepRepr_RepresentationItem;
        aRI1->Init(new TCollection_HAsciiString("lower limit"));
        occ::handle<StepRepr_RepresentationItem> aRI2 = new StepRepr_RepresentationItem;
        aRI2->Init(new TCollection_HAsciiString("upper limit"));
        occ::handle<StepRepr_ReprItemAndLengthMeasureWithUnit> aRILMU1 =
          new StepRepr_ReprItemAndLengthMeasureWithUnit;
        aRILMU1->Init(aMWU1, aRI1);
        occ::handle<StepRepr_ReprItemAndLengthMeasureWithUnit> aRILMU2 =
          new StepRepr_ReprItemAndLengthMeasureWithUnit;
        aRILMU2->Init(aMWU2, aRI2);
        aModel->AddWithRefs(aRILMU1);
        aModel->AddWithRefs(aRILMU2);
        occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> aHARIVR =
          new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, 2);
        aHARIVR->SetValue(1, aRILMU1);
        aHARIVR->SetValue(2, aRILMU2);
        occ::handle<StepRepr_ValueRange> aVR = new StepRepr_ValueRange;
        // VR->Init(aName,CID);
        aVR->Init(aName, aHARIVR);
        aModel->AddEntity(aVR);
        occ::handle<StepShape_ShapeDimensionRepresentation> aSDimR =
          new StepShape_ShapeDimensionRepresentation;
        occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> aHARI1 =
          new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, 1);
        aHARI1->SetValue(1, aVR);
        aSDimR->Init(aName, aHARI1, aRC);
        aModel->AddWithRefs(aSDimR);
        occ::handle<StepShape_DimensionalCharacteristicRepresentation> aDimCharR =
          new StepShape_DimensionalCharacteristicRepresentation;
        StepShape_DimensionalCharacteristic aDimChar;
        aDimChar.SetValue(aDimSize);
        aDimCharR->Init(aDimChar, aSDimR);
        aModel->AddEntity(aDimCharR);
      }
    }
    else if (aKind < 50)
    { // tolerance
      if (aKind < 35)
      { // tolerance with datum system
        NCollection_Sequence<TDF_Label> aDatumLabels;
        XCAFDoc_DimTolTool::GetDatumOfTolerLabels(aDimTolL, aDatumLabels);
        int                                                                      aSetDatumInd = 1;
        occ::handle<NCollection_HArray1<occ::handle<StepDimTol_DatumReference>>> aHADR =
          new NCollection_HArray1<occ::handle<StepDimTol_DatumReference>>(1, aDatumLabels.Length());
        for (NCollection_Sequence<TDF_Label>::Iterator aDatumIter(aDatumLabels); aDatumIter.More();
             aDatumIter.Next(), aSetDatumInd++)
        {
          occ::handle<XCAFDoc_Datum> aDatumAttr;
          const TDF_Label&           aDatumL = aDatumIter.Value();
          if (!aDatumL.FindAttribute(XCAFDoc_Datum::GetID(), aDatumAttr))
          {
            continue;
          }
          occ::handle<TCollection_HAsciiString> aNameD            = aDatumAttr->GetName();
          occ::handle<TCollection_HAsciiString> aDescriptionD     = aDatumAttr->GetDescription();
          occ::handle<TCollection_HAsciiString> anIdentificationD = aDatumAttr->GetIdentification();
          TCollection_AsciiString               aStmp(aNameD->ToCString());
          aStmp.AssignCat(aDescriptionD->ToCString());
          aStmp.AssignCat(anIdentificationD->ToCString());
          if (aDatumMap.IsBound(aStmp))
          {
            occ::handle<StepDimTol_Datum> aDatum =
              occ::down_cast<StepDimTol_Datum>(aDatumMap.Find(aStmp));
            occ::handle<StepDimTol_DatumReference> aDR = new StepDimTol_DatumReference;
            aDR->Init(aSetDatumInd, aDatum);
            aModel->AddWithRefs(aDR);
            aHADR->SetValue(aSetDatumInd, aDR);
          }
        }
        // create LengthMeasureWithUnit
        occ::handle<StepBasic_MeasureValueMember> aMVM = new StepBasic_MeasureValueMember;
        aMVM->SetName("LENGTH_MEASURE");
        aMVM->SetReal(aVal->Value(1));
        occ::handle<StepBasic_LengthMeasureWithUnit> aLMWU = new StepBasic_LengthMeasureWithUnit;
        aLMWU->Init(aMVM, aUnit);
        // create tolerance by it's type
        if (aKind < 24)
        {
          occ::handle<StepDimTol_GeometricToleranceWithDatumReference> aGTWDR =
            new StepDimTol_GeometricToleranceWithDatumReference;
          aGTWDR->SetDatumSystem(aHADR);
          occ::handle<StepDimTol_ModifiedGeometricTolerance> aMGT =
            new StepDimTol_ModifiedGeometricTolerance;
          if (aKind == 21)
          {
            aMGT->SetModifier(StepDimTol_MaximumMaterialCondition);
          }
          else if (aKind == 22)
          {
            aMGT->SetModifier(StepDimTol_LeastMaterialCondition);
          }
          else if (aKind == 23)
          {
            aMGT->SetModifier(StepDimTol_RegardlessOfFeatureSize);
          }
          occ::handle<StepDimTol_GeoTolAndGeoTolWthDatRefAndModGeoTolAndPosTol> aGTComplex =
            new StepDimTol_GeoTolAndGeoTolWthDatRefAndModGeoTolAndPosTol;
          aGTComplex->Init(aName, aDescription, aLMWU, aSA, aGTWDR, aMGT);
          aModel->AddWithRefs(aGTComplex);
        }
        else if (aKind == 24)
        {
          occ::handle<StepDimTol_AngularityTolerance> aToler = new StepDimTol_AngularityTolerance;
          aToler->Init(aName, aDescription, aLMWU, aSA, aHADR);
          aModel->AddWithRefs(aToler);
        }
        else if (aKind == 25)
        {
          occ::handle<StepDimTol_CircularRunoutTolerance> aToler =
            new StepDimTol_CircularRunoutTolerance;
          aToler->Init(aName, aDescription, aLMWU, aSA, aHADR);
          aModel->AddWithRefs(aToler);
        }
        else if (aKind == 26)
        {
          occ::handle<StepDimTol_CoaxialityTolerance> aToler = new StepDimTol_CoaxialityTolerance;
          aToler->Init(aName, aDescription, aLMWU, aSA, aHADR);
          aModel->AddWithRefs(aToler);
        }
        else if (aKind == 27)
        {
          occ::handle<StepDimTol_ConcentricityTolerance> aToler =
            new StepDimTol_ConcentricityTolerance;
          aToler->Init(aName, aDescription, aLMWU, aSA, aHADR);
          aModel->AddWithRefs(aToler);
        }
        else if (aKind == 28)
        {
          occ::handle<StepDimTol_ParallelismTolerance> aToler = new StepDimTol_ParallelismTolerance;
          aToler->Init(aName, aDescription, aLMWU, aSA, aHADR);
          aModel->AddWithRefs(aToler);
        }
        else if (aKind == 29)
        {
          occ::handle<StepDimTol_PerpendicularityTolerance> aToler =
            new StepDimTol_PerpendicularityTolerance;
          aToler->Init(aName, aDescription, aLMWU, aSA, aHADR);
          aModel->AddWithRefs(aToler);
        }
        else if (aKind == 30)
        {
          occ::handle<StepDimTol_SymmetryTolerance> aToler = new StepDimTol_SymmetryTolerance;
          aToler->Init(aName, aDescription, aLMWU, aSA, aHADR);
          aModel->AddWithRefs(aToler);
        }
        else if (aKind == 31)
        {
          occ::handle<StepDimTol_TotalRunoutTolerance> aToler = new StepDimTol_TotalRunoutTolerance;
          aToler->Init(aName, aDescription, aLMWU, aSA, aHADR);
          aModel->AddWithRefs(aToler);
        }
      }
    }
  }

  return true;
}

//! Text as a STEP file has it: ASCII, the rest in \X\ and \X2\ directives (ISO 10303-21), as
//! the reader decodes them. The XCAF strings are UTF-8
static occ::handle<TCollection_HAsciiString> stepText(
  const occ::handle<TCollection_HAsciiString>& theText)
{
  if (theText.IsNull())
  {
    return new TCollection_HAsciiString();
  }
  const TCollection_ExtendedString aWide(theText->ToCString(), true);
  TCollection_AsciiString          anOut;
  bool                             isWide = false;
  for (int i = 1; i <= aWide.Length(); i++)
  {
    const char16_t aChar = aWide.Value(i);
    // Printable ASCII as it is; the rest -- a note's line breaks too, as
    // NIST's files give them -- in UTF-16
    if (aChar >= 0x20 && aChar < 0x7F)
    {
      if (isWide)
      {
        anOut += "\\X0\\";
        isWide = false;
      }
      anOut += static_cast<char>(aChar);
      continue;
    }
    // A character of ISO 8859-1 -- a line break, a degree sign -- in the 8-bit directive, as
    // NIST's files give a note's line breaks ('\X\0A'): in UTF-16 a reader (SFA) counts the
    // text as Unicode
    if (aChar < 0x100)
    {
      if (isWide)
      {
        anOut += "\\X0\\";
        isWide = false;
      }
      char aHex[3];
      std::snprintf(aHex, sizeof(aHex), "%02X", static_cast<unsigned>(aChar));
      anOut += "\\X\\";
      anOut += aHex;
      continue;
    }
    if (!isWide)
    {
      anOut += "\\X2\\";
      isWide = true;
    }
    char aHex[5];
    std::snprintf(aHex, sizeof(aHex), "%04X", static_cast<unsigned>(aChar));
    anOut += aHex;
  }
  if (isWide)
  {
    anOut += "\\X0\\";
  }
  return new TCollection_HAsciiString(anOut);
}

//! A note's words, its descriptions, as AP242's recommended practice writes them: a 'semantic
//! text', a property of what the note is on, its words a representation's descriptive items.
//! The representation is given its context when the draughting model has one (theTexts).
//! Returns the property; null if the note has no words
static occ::handle<StepRepr_PropertyDefinition> writeSemanticText(
  const occ::handle<Interface_InterfaceModel>&                theModel,
  const occ::handle<XCAFDimTolObjects_DimensionObject>&       theObject,
  const occ::handle<StepRepr_ShapeAspect>&                    theOn,
  NCollection_Sequence<occ::handle<StepRepr_Representation>>& theTexts)
{
  occ::handle<StepRepr_PropertyDefinition> aText;
  if (!theObject->HasDescriptions() || theOn.IsNull())
  {
    return aText;
  }
  occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> anItems =
    new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(
      1,
      theObject->NbDescriptions());
  for (int i = 0; i < theObject->NbDescriptions(); i++)
  {
    occ::handle<StepRepr_DescriptiveRepresentationItem> aWord =
      new StepRepr_DescriptiveRepresentationItem();
    occ::handle<TCollection_HAsciiString> aName = theObject->GetDescriptionName(i);
    occ::handle<TCollection_HAsciiString> aWords = theObject->GetDescription(i);
    aWord->Init(stepText(aName), stepText(aWords));
    anItems->SetValue(i + 1, aWord);
  }
  occ::handle<StepRepr_Representation> aRepr = new StepRepr_Representation();
  aRepr->Init(new TCollection_HAsciiString(), anItems, nullptr);
  StepRepr_CharacterizedDefinition anOn;
  anOn.SetValue(theOn);
  aText = new StepRepr_PropertyDefinition();
  aText->Init(new TCollection_HAsciiString("semantic text"),
              true,
              new TCollection_HAsciiString(),
              anOn);
  StepRepr_RepresentedDefinition aDefinition;
  aDefinition.SetValue(aText);
  occ::handle<StepRepr_PropertyDefinitionRepresentation> aPDR =
    new StepRepr_PropertyDefinitionRepresentation();
  aPDR->Init(aDefinition, aRepr);
  theModel->AddWithRefs(aPDR);
  // Which property it is, as NIST's files say it
  occ::handle<StepBasic_GeneralProperty> aKind = new StepBasic_GeneralProperty();
  aKind->Init(new TCollection_HAsciiString(),
              new TCollection_HAsciiString("semantic text"),
              false,
              nullptr);
  occ::handle<StepBasic_GeneralPropertyAssociation> anAssociation =
    new StepBasic_GeneralPropertyAssociation();
  anAssociation->Init(new TCollection_HAsciiString(),
                      new TCollection_HAsciiString(),
                      aKind,
                      aText);
  theModel->AddWithRefs(anAssociation);
  theTexts.Append(aRepr);
  return aText;
}

//=======================================================================
// function : writeSupplementalGeometry
// purpose  : each group of supplemental geometry as a CONSTRUCTIVE_GEOMETRY_REPRESENTATION
//            related to its part's shape representation; each set of curves bound to the
//            shape it is, for GD&T on it to find
//=======================================================================
static void writeSupplementalGeometry(const occ::handle<XSControl_WorkSession>& theWS,
                                      const Interface_Graph&                    theGraph,
                                      const occ::handle<XCAFDoc_DimTolTool>&    theDGTTool,
                                      GDTPresentationModel&                     theState,
                                      const StepData_Factors&                   theLocalFactors)
{
  const occ::handle<Interface_InterfaceModel>& aModel = theWS->Model();
  const occ::handle<Transfer_FinderProcess>&   aFP    = theWS->TransferWriter()->FinderProcess();
  NCollection_Sequence<TDF_Label>              aGroups;
  theDGTTool->GetSupplementalGeometry(aGroups);
  // The items' colours, styles of a presentation of their own
  STEPConstruct_Styles                        aStyles(theWS);
  occ::handle<StepRepr_RepresentationContext> aStyleContext;
  for (NCollection_Sequence<TDF_Label>::Iterator aGroupIter(aGroups); aGroupIter.More();
       aGroupIter.Next())
  {
    TDF_Label               aPartL;
    TCollection_AsciiString aName, aRelationName;
    if (!XCAFDoc_DimTolTool::GetSupplementalGeometry(aGroupIter.Value(),
                                                     aPartL,
                                                     aName,
                                                     aRelationName)
        || aPartL.IsNull())
    {
      continue;
    }
    TDF_Label aProtoL = aPartL;
    if (XCAFDoc_ShapeTool::IsReference(aPartL))
    {
      XCAFDoc_ShapeTool::GetReferredShape(aPartL, aProtoL);
    }
    // The part's shape representation: the one an element of its shape was written in
    const TopoDS_Shape                           aPart = XCAFDoc_ShapeTool::GetShape(aProtoL);
    occ::handle<StepRepr_ProductDefinitionShape> aPDS;
    occ::handle<StepRepr_RepresentationContext>  aRC;
    for (TopAbs_ShapeEnum aType : {TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX})
    {
      for (TopExp_Explorer anExp(aPart, aType); anExp.More() && aPDS.IsNull(); anExp.Next())
      {
        TopLoc_Location                                       aLoc;
        NCollection_Sequence<occ::handle<Standard_Transient>> aSeqRI;
        if (FindEntities(aFP, anExp.Current(), aLoc, aSeqRI) > 0)
        {
          aPDS = FindPDS(theGraph, aSeqRI.First(), aRC);
        }
      }
      if (!aPDS.IsNull())
      {
        break;
      }
    }
    // A tessellated part (NIST FTC-08 -e4-tg) has no item for its faces: its shapes' own
    for (TopoDS_Iterator anIter(aPart); anIter.More() && aPDS.IsNull(); anIter.Next())
    {
      TopLoc_Location                                       aLoc;
      NCollection_Sequence<occ::handle<Standard_Transient>> aSeqRI;
      if (FindEntities(aFP, anIter.Value(), aLoc, aSeqRI) > 0)
      {
        aPDS = FindPDS(theGraph, aSeqRI.First(), aRC);
      }
    }
    if (aPDS.IsNull())
    {
      TopLoc_Location                                       aLoc;
      NCollection_Sequence<occ::handle<Standard_Transient>> aSeqRI;
      if (FindEntities(aFP, aPart, aLoc, aSeqRI) > 0)
      {
        aPDS = FindPDS(theGraph, aSeqRI.First(), aRC);
      }
    }
    occ::handle<StepShape_ShapeDefinitionRepresentation> aSDR;
    for (Interface_EntityIterator anIter = aPDS.IsNull() ? Interface_EntityIterator()
                                                         : theGraph.Sharings(aPDS);
         anIter.More() && aSDR.IsNull();
         anIter.Next())
    {
      aSDR = occ::down_cast<StepShape_ShapeDefinitionRepresentation>(anIter.Value());
    }
    if (aSDR.IsNull() || aSDR->UsedRepresentation().IsNull())
    {
      continue;
    }
    occ::handle<StepRepr_Representation> aPartRep = aSDR->UsedRepresentation();

    NCollection_Sequence<TDF_Label> anItemLs;
    XCAFDoc_DimTolTool::GetSupplementalItems(aGroupIter.Value(), anItemLs);
    NCollection_Sequence<occ::handle<StepRepr_RepresentationItem>> anItems;
    NCollection_Sequence<std::pair<occ::handle<Standard_Transient>, TopoDS_Shape>> aSets;
    for (NCollection_Sequence<TDF_Label>::Iterator anItemIter(anItemLs); anItemIter.More();
         anItemIter.Next())
    {
      TCollection_AsciiString aKind, anItemName;
      gp_Ax3                  aPos;
      double                  aRadius = 0.0;
      TopoDS_Shape            aShape;
      if (!XCAFDoc_DimTolTool::GetSupplementalItem(anItemIter.Value(),
                                                   aKind,
                                                   anItemName,
                                                   aPos,
                                                   aRadius,
                                                   aShape))
      {
        continue;
      }
      occ::handle<StepRepr_RepresentationItem> anItem;
      try
      {
        OCC_CATCH_SIGNALS
        if (aKind == "plane")
        {
          GeomToStep_MakePlane aMaker(gp_Pln(aPos), theLocalFactors);
          anItem = aMaker.Value();
        }
        else if (aKind == "cylinder" && aRadius > 0.0)
        {
          GeomToStep_MakeCylindricalSurface aMaker(new Geom_CylindricalSurface(aPos, aRadius),
                                                   theLocalFactors);
          anItem = aMaker.Value();
        }
        else if (aKind == "axes")
        {
          GeomToStep_MakeAxis2Placement3d aMaker(aPos, theLocalFactors);
          anItem = aMaker.Value();
        }
        else if (aKind == "curves" && !aShape.IsNull())
        {
          // Its curves, and its points that are on none
          NCollection_Sequence<StepShape_GeometricSetSelect> anElements;
          BRep_Builder                                       aBuilder;
          TopoDS_Compound                                    anEdges;
          aBuilder.MakeCompound(anEdges);
          bool hasEdges = false;
          for (TopExp_Explorer anExp(aShape, TopAbs_EDGE); anExp.More(); anExp.Next())
          {
            aBuilder.Add(anEdges, anExp.Current());
            hasEdges = true;
          }
          if (hasEdges)
          {
            TopoDSToStep_MakeGeometricCurveSet aMaker(anEdges, aFP, theLocalFactors);
            if (aMaker.IsDone())
            {
              for (int i = 1; i <= aMaker.Value()->NbElements(); ++i)
              {
                anElements.Append(aMaker.Value()->ElementsValue(i));
              }
            }
          }
          for (TopExp_Explorer anExp(aShape, TopAbs_VERTEX, TopAbs_EDGE); anExp.More();
               anExp.Next())
          {
            GeomToStep_MakeCartesianPoint aMaker(BRep_Tool::Pnt(TopoDS::Vertex(anExp.Current())),
                                                 theLocalFactors.LengthFactor());
            StepShape_GeometricSetSelect  aSelect;
            aSelect.SetValue(aMaker.Value());
            anElements.Append(aSelect);
          }
          if (!anElements.IsEmpty())
          {
            occ::handle<NCollection_HArray1<StepShape_GeometricSetSelect>> anArray =
              new NCollection_HArray1<StepShape_GeometricSetSelect>(1, anElements.Length());
            for (int i = 1; i <= anElements.Length(); ++i)
            {
              anArray->SetValue(i, anElements.Value(i));
            }
            occ::handle<StepShape_GeometricCurveSet> aSet = new StepShape_GeometricCurveSet();
            aSet->Init(new TCollection_HAsciiString(), anArray);
            anItem = aSet;
            aSets.Append(std::make_pair(occ::handle<Standard_Transient>(aSet), aShape));
          }
        }
      }
      catch (const Standard_Failure&)
      {
        anItem.Nullify();
      }
      if (anItem.IsNull())
      {
        continue;
      }
      anItem->SetName(new TCollection_HAsciiString(anItemName));
      anItems.Append(anItem);
      Quantity_Color aColor;
      if (XCAFDoc_DimTolTool::GetSupplementalColor(anItemIter.Value(), aColor))
      {
        // As NIST's: one style each. A surface's is its fill and its boundary in the curve's
        // colour (OCCT cannot write the parameter lines NIST's have too) (SFA reads the one style of an assignment, and stops at
        // two: "Error getting color"); a curve set's, its curves' -- and each of its curves
        // and points is drawn so, where SFA takes the colour from
        occ::handle<StepVisual_Colour> aColour = STEPConstruct_Styles::EncodeColor(aColor);
        const STEPConstruct_RenderingProperties aProps;
        occ::handle<StepVisual_PresentationStyleAssignment> aCurvePSA =
          aStyles.MakeColorPSA(anItem, nullptr, aColour, aProps);
        occ::handle<StepVisual_PresentationStyleAssignment> aPSA = aCurvePSA;
        if (aKind == "plane" || aKind == "cylinder")
        {
          aPSA = aStyles.MakeColorPSA(anItem, aColour, nullptr, aProps);
          occ::handle<StepVisual_CurveStyle> aCurve =
            aCurvePSA.IsNull() || aCurvePSA->NbStyles() < 1
              ? nullptr
              : aCurvePSA->StylesValue(1).CurveStyle();
          occ::handle<StepVisual_SurfaceStyleUsage> aUsage =
            aPSA.IsNull() || aPSA->NbStyles() < 1 ? nullptr
                                                  : aPSA->StylesValue(1).SurfaceStyleUsage();
          occ::handle<StepVisual_SurfaceSideStyle> aSide =
            aUsage.IsNull() ? nullptr : aUsage->Style();
          if (!aCurve.IsNull() && !aSide.IsNull() && !aSide->Styles().IsNull())
          {
            occ::handle<StepVisual_SurfaceStyleBoundary> aBoundary =
              new StepVisual_SurfaceStyleBoundary();
            aBoundary->Init(aCurve);
            const occ::handle<NCollection_HArray1<StepVisual_SurfaceStyleElementSelect>>& anOld =
              aSide->Styles();
            occ::handle<NCollection_HArray1<StepVisual_SurfaceStyleElementSelect>> aNew =
              new NCollection_HArray1<StepVisual_SurfaceStyleElementSelect>(1,
                                                                            anOld->Length() + 1);
            int j = 1;
            for (int i = anOld->Lower(); i <= anOld->Upper(); ++i)
            {
              aNew->SetValue(j++, anOld->Value(i));
            }
            StepVisual_SurfaceStyleElementSelect aBoundarySelect;
            aBoundarySelect.SetValue(aBoundary);
            aNew->SetValue(j, aBoundarySelect);
            aSide->SetStyles(aNew);
          }
        }
        auto aStyle = [&aStyles](const occ::handle<StepRepr_RepresentationItem>&           theItem,
                                 const occ::handle<StepVisual_PresentationStyleAssignment>& thePSA) {
          if (theItem.IsNull() || thePSA.IsNull())
          {
            return;
          }
          occ::handle<NCollection_HArray1<occ::handle<StepVisual_PresentationStyleAssignment>>>
            aPSAs = new NCollection_HArray1<occ::handle<StepVisual_PresentationStyleAssignment>>(1,
                                                                                                1);
          aPSAs->SetValue(1, thePSA);
          occ::handle<StepVisual_StyledItem> aStyled = new StepVisual_StyledItem();
          aStyled->Init(new TCollection_HAsciiString(), aPSAs, theItem);
          aStyles.AddStyle(aStyled);
        };
        aStyle(anItem, aPSA);
        occ::handle<StepShape_GeometricSet> aSet = occ::down_cast<StepShape_GeometricSet>(anItem);
        // SFA colours a trimmed curve by its basis curve, any other member by itself
        // (NIST FTC-10's lines)
        for (int i = 1; !aSet.IsNull() && i <= aSet->NbElements(); ++i)
        {
          occ::handle<StepRepr_RepresentationItem> aMember =
            occ::down_cast<StepRepr_RepresentationItem>(aSet->ElementsValue(i).Value());
          if (occ::handle<StepGeom_TrimmedCurve> aTrimmed =
                occ::down_cast<StepGeom_TrimmedCurve>(aMember))
          {
            aMember = aTrimmed->BasisCurve();
          }
          aStyle(aMember, aCurvePSA);
        }
      }
    }
    if (anItems.IsEmpty())
    {
      continue;
    }
    occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> anArray =
      new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, anItems.Length());
    for (int i = 1; i <= anItems.Length(); ++i)
    {
      anArray->SetValue(i, anItems.Value(i));
    }
    occ::handle<StepRepr_ConstructiveGeometryRepresentation> aCGR =
      new StepRepr_ConstructiveGeometryRepresentation();
    aCGR->Init(new TCollection_HAsciiString(aName), anArray, aPartRep->ContextOfItems());
    occ::handle<StepRepr_ConstructiveGeometryRepresentationRelationship> aCGRR =
      new StepRepr_ConstructiveGeometryRepresentationRelationship();
    aCGRR->Init(new TCollection_HAsciiString(aRelationName),
                new TCollection_HAsciiString(),
                aPartRep,
                aCGR);
    aModel->AddWithRefs(aCGRR);
    if (aStyleContext.IsNull())
    {
      aStyleContext = aPartRep->ContextOfItems();
    }
    // GD&T on a set of curves: the set, in this representation
    for (const auto& aSet : aSets)
    {
      GDTPresentationModel::Supplemental anEntry;
      anEntry.Representation = aCGR;
      anEntry.PDS            = aPDS;
      theState.Supplementals.Bind(aSet.first, anEntry);
      occ::handle<Transfer_SimpleBinderOfTransient> aBinder =
        new Transfer_SimpleBinderOfTransient();
      aBinder->SetResult(aSet.first);
      aFP->Bind(TransferBRep::ShapeMapper(aFP, aSet.second), aBinder);
    }
  }
  if (aStyles.NbStyles() > 0 && !aStyleContext.IsNull())
  {
    occ::handle<StepVisual_MechanicalDesignGeometricPresentationRepresentation> aMDGPR;
    occ::handle<StepData_StepModel> aStepModel = occ::down_cast<StepData_StepModel>(aModel);
    if (aStyles.CreateMDGPR(aStyleContext, aMDGPR, aStepModel) && !aMDGPR.IsNull())
    {
      aModel->AddWithRefs(aMDGPR);
    }
  }
}

//=================================================================================================

bool STEPCAFControl_Writer::writeDGTsAP242(const occ::handle<XSControl_WorkSession>& theWS,
                                           const NCollection_Sequence<TDF_Label>&    theLabels,
                                           const StepData_Factors& theLocalFactors)
{
  (void)theLabels;
  // Get working data
  const occ::handle<Interface_InterfaceModel>& aModel = theWS->Model();

  const occ::handle<Interface_HGraph> aHGraph = theWS->HGraph();
  if (aHGraph.IsNull())
  {
    return false;
  }

  Interface_Graph aGraph = aHGraph->Graph();

  // Common entities for presentation
  STEPConstruct_Styles           aStyles(theWS);
  occ::handle<StepVisual_Colour> aCurvColor = STEPConstruct_Styles::EncodeColor(gdtState(myGDTPresentationDM).Color);
  occ::handle<StepRepr_RepresentationItem> anItem = nullptr;
  myGDTPrsCurveStyle->SetValue(
    1,
    aStyles.MakeColorPSA(anItem, aCurvColor, aCurvColor, STEPConstruct_RenderingProperties()));
  for (Interface_EntityIterator aModelIter = aModel->Entities();
       aModelIter.More() && myGDTCommonPDS.IsNull();
       aModelIter.Next())
  {
    myGDTCommonPDS = occ::down_cast<StepRepr_ProductDefinitionShape>(aModelIter.Value());
  }

  // Supplemental geometry first: GD&T on its curves identifies its items
  {
    NCollection_Map<occ::handle<Standard_Transient>> aTools;
    for (NCollection_Map<TDF_Label>::Iterator aLabelIter(myRootLabels); aLabelIter.More();
         aLabelIter.Next())
    {
      occ::handle<XCAFDoc_DimTolTool> aDGTTool =
        XCAFDoc_DocumentTool::DimTolTool(aLabelIter.Value());
      if (!aDGTTool.IsNull() && aTools.Add(aDGTTool))
      {
        writeSupplementalGeometry(theWS,
                                  aGraph,
                                  aDGTTool,
                                  gdtState(myGDTPresentationDM),
                                  theLocalFactors);
      }
    }
  }

  NCollection_DataMap<TCollection_AsciiString, occ::handle<Standard_Transient>> aDatumMap;
  occ::handle<StepRepr_RepresentationContext>                                   aRC;

  NCollection_Sequence<TDF_Label> aDGTLabels;
  // Iterate on requested shapes collect Tools
  for (NCollection_Map<TDF_Label>::Iterator aLabelIter(myRootLabels); aLabelIter.More();
       aLabelIter.Next())
  {
    const TDF_Label&                aLabel   = aLabelIter.Value();
    occ::handle<XCAFDoc_DimTolTool> aDGTTool = XCAFDoc_DocumentTool::DimTolTool(aLabel);
    NCollection_Sequence<TDF_Label> aaDGTLS;
    aDGTTool->GetDatumLabels(aDGTLabels);
    aDGTLabels.Append(aaDGTLS);
  }

  // Find all shapes with datums
  NCollection_Map<TCollection_AsciiString> aNameIdMap, aWrittenNameIds, aOfNoPart;
  for (NCollection_Sequence<TDF_Label>::Iterator aDGTIter(aDGTLabels); aDGTIter.More();
       aDGTIter.Next())
  {
    const TDF_Label&                aDatumL = aDGTIter.Value();
    NCollection_Sequence<TDF_Label> aShapeL, aNullSeq;
    XCAFDoc_DimTolTool::GetRefShapeLabel(aDatumL, aShapeL, aNullSeq);
    occ::handle<XCAFDoc_Datum> aDatumAttr;
    aDatumL.FindAttribute(XCAFDoc_Datum::GetID(), aDatumAttr);
    occ::handle<XCAFDimTolObjects_DatumObject> anObject   = aDatumAttr->GetObject();
    TCollection_AsciiString                    aDatumName = anObject->GetName()->String();
    TCollection_AsciiString                    aDatumTargetId =
      TCollection_AsciiString(anObject->GetDatumTargetNumber());
    // A letter is one part's: two parts written to one file may each have a datum A, each
    // written with its own feature (by name alone the second part's was taken for the first's
    // given again, and dropped). A label of no shape -- a frame citing the datum -- is of the
    // letter as first written
    const TCollection_AsciiString aPart     = datumPartKey(aDatumL);
    const TCollection_AsciiString aNameId   = aDatumName.Cat(aDatumTargetId);
    const TCollection_AsciiString aDatumKey = aPart.IsEmpty() ? aDatumName
                                                              : aDatumName.Cat("@").Cat(aPart);
    if ((aPart.IsEmpty() && aWrittenNameIds.Contains(aNameId))
        || !aNameIdMap.Add(aPart.IsEmpty() ? aNameId : aNameId.Cat("@").Cat(aPart)))
    {
      continue;
    }
    aWrittenNameIds.Add(aNameId);
    // A datum written of no part -- from a frame citing it, before its own label came: one
    // established by targets alone, cited from another set (NIST CTC-02's K) -- is the one its
    // targets are of, whatever part they are on
    occ::handle<Standard_Transient> aWrittenDatum;
    bool                            isFirstDT = !aDatumMap.Find(aDatumKey, aWrittenDatum);
    // A target is of its letter's datum on its own part, else of its letter's datum wherever
    // it is: a datum's targets may be on faces written as a part of their own (NIST FTC-10's K,
    // its feature on the solid, its areas on surfaces beside it)
    if (isFirstDT && !aPart.IsEmpty()
        && (aOfNoPart.Contains(aDatumName) || anObject->IsDatumTarget()))
    {
      isFirstDT = !aDatumMap.Find(aDatumName, aWrittenDatum);
    }
    if (aPart.IsEmpty())
    {
      aOfNoPart.Add(aDatumName);
    }
    occ::handle<StepDimTol_Datum>   aDatum =
      writeDatumAP242(theWS,
                      aShapeL,
                      aDatumL,
                      isFirstDT,
                      occ::down_cast<StepDimTol_Datum>(aWrittenDatum),
                      theLocalFactors);
    // Add created Datum into Map: as its part's, and as the letter's when it is the first
    aDatumMap.Bind(aDatumKey, aDatum);
    if (!aDatumMap.IsBound(aDatumName))
    {
      aDatumMap.Bind(aDatumName, aDatum);
    }
  }

  // write Dimensions
  aDGTLabels.Clear();
  for (NCollection_Map<TDF_Label>::Iterator aLabelIter(myRootLabels); aLabelIter.More();
       aLabelIter.Next())
  {
    const TDF_Label&                aLabel   = aLabelIter.Value();
    occ::handle<XCAFDoc_DimTolTool> aDGTTool = XCAFDoc_DocumentTool::DimTolTool(aLabel);
    NCollection_Sequence<TDF_Label> aaDGTLS;
    aDGTTool->GetDimensionLabels(aDGTLabels);
    aDGTLabels.Append(aaDGTLS);
  }

  // Auxiliary entities for derived geometry
  occ::handle<StepRepr_ConstructiveGeometryRepresentation> aCGRepr =
    new StepRepr_ConstructiveGeometryRepresentation();
  occ::handle<StepRepr_ConstructiveGeometryRepresentationRelationship> aCGReprRel =
    new StepRepr_ConstructiveGeometryRepresentationRelationship();
  NCollection_DynamicArray<occ::handle<StepGeom_CartesianPoint>> aConnectionPnts;
  occ::handle<StepRepr_RepresentationContext>                    dummyRC;
  occ::handle<StepAP242_GeometricItemSpecificUsage>              dummyGISU;
  // The notes' words, given the draughting model's context once it has one
  NCollection_Sequence<occ::handle<StepRepr_Representation>> aTexts;
  for (NCollection_Sequence<TDF_Label>::Iterator aDGTIter(aDGTLabels); aDGTIter.More();
       aDGTIter.Next())
  {
    const TDF_Label&                aDimensionL = aDGTIter.Value();
    NCollection_Sequence<TDF_Label> aFirstShapeL, aSecondShapeL;
    occ::handle<XCAFDoc_Dimension>  aDimAttr;
    if (!aDimensionL.FindAttribute(XCAFDoc_Dimension::GetID(), aDimAttr))
    {
      continue;
    }
    occ::handle<XCAFDimTolObjects_DimensionObject> anObject = aDimAttr->GetObject();
    if (anObject.IsNull())
    {
      continue;
    }
    if (anObject->GetType() == XCAFDimTolObjects_DimensionType_CommonLabel)
    {
      // The product's: its words a property of the product's aspect. With
      // none, on nothing: its graphic associated with nothing, as a graphic
      // alone is (NIST CTC-05's regions)
      const bool                        onNothing = anObject->NbDescriptions() == 0;
      occ::handle<StepRepr_ShapeAspect> aSA;
      if (!onNothing)
      {
        aSA = new StepRepr_ShapeAspect();
        aSA->Init(new TCollection_HAsciiString(),
                  new TCollection_HAsciiString(),
                  myGDTCommonPDS,
                  StepData_LTrue);
        aModel->AddWithRefs(aSA);
      }
      // No semantics but its words, if it has any, and its plane if it has one (the two were
      // swapped, which left every graphic-only item without its annotation plane)
      occ::handle<StepRepr_PropertyDefinition> aText =
        aSA.IsNull() ? occ::handle<StepRepr_PropertyDefinition>()
                     : writeSemanticText(aModel, anObject, aSA, aTexts);
      writePresentation(theWS,
                        aDimensionL,
                        anObject->GetPresentation(),
                        anObject->GetPresentationName(),
                        !aText.IsNull(),
                        anObject->HasPlane(),
                        anObject->GetPlane(),
                        anObject->GetPointTextAttach(),
                        aText.IsNull() ? occ::handle<Standard_Transient>(aSA) : aText,
                        keptOr(anObject->GetPresentationType(), "note"),
                        NCollection_Sequence<occ::handle<StepRepr_ShapeAspect>>(),
                        theLocalFactors,
                        onNothing);
    }

    if (!XCAFDoc_DimTolTool::GetRefShapeLabel(aDimensionL, aFirstShapeL, aSecondShapeL))
    {
      continue;
    }

    // Write links with shapes
    occ::handle<StepRepr_ShapeAspect> aFirstSA, aSecondSA;
    if (aFirstShapeL.Length() == 1)
    {
      TopoDS_Shape aShape = XCAFDoc_ShapeTool::GetShape(aFirstShapeL.Value(1));
      aFirstSA            = writeShapeAspect(theWS, aDimensionL, aShape, dummyRC, dummyGISU);
      if (aRC.IsNull() && !dummyRC.IsNull())
      {
        aRC = dummyRC;
      }
    }
    else if (aFirstShapeL.Length() > 1)
    {
      occ::handle<StepRepr_CompositeShapeAspect> aCSA;
      for (int shIt = 1; shIt <= aFirstShapeL.Length(); shIt++)
      {
        TopoDS_Shape aShape = XCAFDoc_ShapeTool::GetShape(aFirstShapeL.Value(shIt));
        occ::handle<StepRepr_ShapeAspect> aSA =
          writeShapeAspect(theWS, aDimensionL, aShape, dummyRC, dummyGISU);
        if (aSA.IsNull())
        {
          continue;
        }
        if (aCSA.IsNull())
        {
          aCSA = new StepRepr_CompositeShapeAspect();
          aCSA->Init(aSA->Name(), aSA->Description(), aSA->OfShape(), aSA->ProductDefinitional());
          aModel->AddWithRefs(aCSA);
        }
        occ::handle<StepRepr_ShapeAspectRelationship> aSAR = new StepRepr_ShapeAspectRelationship();
        aSAR->Init(new TCollection_HAsciiString(),
                   false,
                   new TCollection_HAsciiString(),
                   aCSA,
                   aSA);
        aModel->AddWithRefs(aSAR);
        if (aRC.IsNull() && !dummyRC.IsNull())
        {
          aRC = dummyRC;
        }
      }
      aFirstSA = aCSA;
    }
    if (aSecondShapeL.Length() == 1)
    {
      TopoDS_Shape aShape = XCAFDoc_ShapeTool::GetShape(aSecondShapeL.Value(1));
      aSecondSA           = writeShapeAspect(theWS, aDimensionL, aShape, dummyRC, dummyGISU);
      if (aRC.IsNull() && !dummyRC.IsNull())
      {
        aRC = dummyRC;
      }
    }
    else if (aSecondShapeL.Length() > 1)
    {
      occ::handle<StepRepr_CompositeShapeAspect> aCSA;
      for (int shIt = 1; shIt <= aSecondShapeL.Length(); shIt++)
      {
        TopoDS_Shape aShape = XCAFDoc_ShapeTool::GetShape(aSecondShapeL.Value(shIt));
        occ::handle<StepRepr_ShapeAspect> aSA =
          writeShapeAspect(theWS, aDimensionL, aShape, dummyRC, dummyGISU);
        if (aCSA.IsNull() && !aSA.IsNull())
        {
          aCSA = new StepRepr_CompositeShapeAspect();
        }
        aCSA->Init(aSA->Name(), aSA->Description(), aSA->OfShape(), aSA->ProductDefinitional());
        if (!aSA.IsNull())
        {
          occ::handle<StepRepr_ShapeAspectRelationship> aSAR =
            new StepRepr_ShapeAspectRelationship();
          aSAR->Init(new TCollection_HAsciiString(),
                     false,
                     new TCollection_HAsciiString(),
                     aCSA,
                     aSA);
          aModel->AddWithRefs(aSAR);
        }
        if (aRC.IsNull() && !dummyRC.IsNull())
        {
          aRC = dummyRC;
        }
      }
      aSecondSA = aCSA;
    }

    if (anObject->GetType() == XCAFDimTolObjects_DimensionType_DimensionPresentation)
    {
      // On nothing, with no words: its graphic associated with nothing
      const bool onNothing = aFirstShapeL.IsEmpty() && anObject->NbDescriptions() == 0;
      if (aFirstSA.IsNull() && !onNothing)
      {
        // Its shapes have no entity to refer to (a tessellated part's faces):
        // the graphic is the product's, as a common label's is
        aFirstSA = new StepRepr_ShapeAspect();
        aFirstSA->Init(new TCollection_HAsciiString(),
                       new TCollection_HAsciiString(),
                       myGDTCommonPDS,
                       StepData_LTrue);
        aModel->AddWithRefs(aFirstSA);
      }
      // No semantics but its words, if it has any -- the graphic their presentation, associated
      // with the faces they are on -- and its plane if it has one (the two were swapped, which
      // left every graphic-only item without its annotation plane)
      occ::handle<StepRepr_PropertyDefinition> aText =
        aFirstSA.IsNull() ? occ::handle<StepRepr_PropertyDefinition>()
                          : writeSemanticText(aModel, anObject, aFirstSA, aTexts);
      writePresentation(theWS,
                        aDimensionL,
                        anObject->GetPresentation(),
                        anObject->GetPresentationName(),
                        !aText.IsNull(),
                        anObject->HasPlane(),
                        anObject->GetPlane(),
                        anObject->GetPointTextAttach(),
                        aText.IsNull() ? occ::handle<Standard_Transient>(aFirstSA) : aText,
                        keptOr(anObject->GetPresentationType(), "note"),
                        aText.IsNull() ? NCollection_Sequence<occ::handle<StepRepr_ShapeAspect>>()
                                       : shapeAspects(aFirstSA),
                        theLocalFactors,
                        onNothing);
      continue;
    }

    // Write dimensions
    StepShape_DimensionalCharacteristic aDimension;
    if (anObject->HasPoint() || anObject->HasPoint2())
    {
      WriteDerivedGeometry(theWS,
                           anObject,
                           aCGRepr,
                           aFirstSA,
                           aSecondSA,
                           aConnectionPnts,
                           theLocalFactors);
    }
    XCAFDimTolObjects_DimensionType aDimType = anObject->GetType();
    if (XCAFDimTolObjects_DimensionObject::IsDimensionalLocation(aDimType))
    {
      // Dimensional_Location
      occ::handle<StepShape_DimensionalLocation> aDim = new StepShape_DimensionalLocation();
      aDim->Init(STEPCAFControl_GDTProperty::GetDimTypeName(aDimType),
                 false,
                 nullptr,
                 aFirstSA,
                 aSecondSA);
      aDimension.SetValue(aDim);
    }
    else if (aDimType == XCAFDimTolObjects_DimensionType_Location_Angular)
    {
      // Angular_Location
      occ::handle<StepShape_AngularLocation> aDim     = new StepShape_AngularLocation();
      StepShape_AngleRelator                 aRelator = StepShape_Equal;
      if (anObject->HasQualifier())
      {
        XCAFDimTolObjects_AngularQualifier aQualifier = anObject->GetAngularQualifier();
        switch (aQualifier)
        {
          case XCAFDimTolObjects_AngularQualifier_Small:
            aRelator = StepShape_Small;
            break;
          case XCAFDimTolObjects_AngularQualifier_Large:
            aRelator = StepShape_Large;
            break;
          default:
            aRelator = StepShape_Equal;
        }
      }
      // Named 'angle', as AP242's recommended practice names it (SFA's dimension name)
      aDim->Init(new TCollection_HAsciiString("angle"), false, nullptr, aFirstSA, aSecondSA,
                 aRelator);
      aDimension.SetValue(aDim);
    }
    else if (aDimType == XCAFDimTolObjects_DimensionType_Location_WithPath)
    {
      // Dimensional_Location_With_Path
      occ::handle<StepShape_DimensionalLocationWithPath> aDim =
        new StepShape_DimensionalLocationWithPath();
      occ::handle<StepRepr_ShapeAspect> aPathSA =
        writeShapeAspect(theWS, aDimensionL, anObject->GetPath(), dummyRC, dummyGISU);
      aDim->Init(new TCollection_HAsciiString(), false, nullptr, aFirstSA, aSecondSA, aPathSA);
      aDimension.SetValue(aDim);
    }
    else if (XCAFDimTolObjects_DimensionObject::IsDimensionalSize(aDimType))
    {
      // Dimensional_Size
      occ::handle<StepShape_DimensionalSize> aDim = new StepShape_DimensionalSize();
      aDim->Init(aFirstSA, STEPCAFControl_GDTProperty::GetDimTypeName(aDimType));
      aDimension.SetValue(aDim);
    }
    else if (aDimType == XCAFDimTolObjects_DimensionType_Size_Angular)
    {
      // Angular_Size
      occ::handle<StepShape_AngularSize> aDim     = new StepShape_AngularSize();
      StepShape_AngleRelator             aRelator = StepShape_Equal;
      if (anObject->HasQualifier())
      {
        XCAFDimTolObjects_AngularQualifier aQualifier = anObject->GetAngularQualifier();
        switch (aQualifier)
        {
          case XCAFDimTolObjects_AngularQualifier_Small:
            aRelator = StepShape_Small;
            break;
          case XCAFDimTolObjects_AngularQualifier_Large:
            aRelator = StepShape_Large;
            break;
          default:
            aRelator = StepShape_Equal;
        }
      }
      aDim->Init(aFirstSA, new TCollection_HAsciiString("angle"), aRelator);
      aDimension.SetValue(aDim);
    }
    else if (aDimType == XCAFDimTolObjects_DimensionType_Size_WithPath)
    {
      // Dimensional_Size_With_Path
      occ::handle<StepShape_DimensionalSizeWithPath> aDim = new StepShape_DimensionalSizeWithPath();
      occ::handle<StepRepr_ShapeAspect>              aPathSA =
        writeShapeAspect(theWS, aDimensionL, anObject->GetPath(), dummyRC, dummyGISU);
      aDim->Init(aFirstSA, new TCollection_HAsciiString(), aPathSA);
      aDimension.SetValue(aDim);
    }

    // Write values
    WriteDimValues(theWS, anObject, aRC, aDimension, theLocalFactors);
    // Annotation plane and Presentation, and what it measures
    NCollection_Sequence<occ::handle<StepRepr_ShapeAspect>> aDimGeometry;
    aDimGeometry.Append(aFirstSA);
    if (!aSecondSA.IsNull())
    {
      aDimGeometry.Append(aSecondSA);
    }
    writePresentation(theWS,
                      aDimensionL,
                      anObject->GetPresentation(),
                      anObject->GetPresentationName(),
                      true,
                      anObject->HasPlane(),
                      anObject->GetPlane(),
                      anObject->GetPointTextAttach(),
                      aDimension.Value(),
                      keptOr(anObject->GetPresentationType(),
                             presentedType(anObject->GetType())),
                      aDimGeometry,
                      theLocalFactors);
  }
  // Write Derived geometry
  if (aConnectionPnts.Length() > 0)
  {
    occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> anItems =
      new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1,
                                                                        aConnectionPnts.Length());
    for (int i = 0; i < aConnectionPnts.Length(); i++)
    {
      anItems->SetValue(i + 1, aConnectionPnts(i));
    }
    aCGRepr->Init(new TCollection_HAsciiString(), anItems, dummyRC);
    aCGReprRel->Init(new TCollection_HAsciiString(),
                     new TCollection_HAsciiString(),
                     dummyGISU->UsedRepresentation(),
                     aCGRepr);
    aModel->AddWithRefs(aCGReprRel);
  }

  // write Geometric Tolerances
  aDGTLabels.Clear();
  for (NCollection_Map<TDF_Label>::Iterator aLabelIter(myRootLabels); aLabelIter.More();
       aLabelIter.Next())
  {
    const TDF_Label&                aLabel   = aLabelIter.Value();
    occ::handle<XCAFDoc_DimTolTool> aDGTTool = XCAFDoc_DocumentTool::DimTolTool(aLabel);
    NCollection_Sequence<TDF_Label> aaDGTLS;
    aDGTTool->GetGeomToleranceLabels(aDGTLabels);
    aDGTLabels.Append(aaDGTLS);
  }
  NCollection_DataMap<TDF_Label, occ::handle<StepDimTol_GeometricTolerance>> aWrittenTols;
  for (NCollection_Sequence<TDF_Label>::Iterator aDGTIter(aDGTLabels); aDGTIter.More();
       aDGTIter.Next())
  {
    const TDF_Label                 aGeomTolL = aDGTIter.Value();
    NCollection_Sequence<TDF_Label> aFirstShapeL, aNullSeqL;
    if (!XCAFDoc_DimTolTool::GetRefShapeLabel(aGeomTolL, aFirstShapeL, aNullSeqL))
    {
      continue;
    }
    NCollection_Sequence<TDF_Label> aDatumSeq;
    XCAFDoc_DimTolTool::GetDatumWithObjectOfTolerLabels(aGeomTolL, aDatumSeq);
    occ::handle<NCollection_HArray1<StepDimTol_DatumSystemOrReference>> aDatumSystem;
    if (aDatumSeq.Length() > 0)
    {
      aDatumSystem = WriteDatumSystem(theWS, aGeomTolL, aDatumSeq, aDatumMap, aRC, theLocalFactors);
    }
    occ::handle<StepDimTol_GeometricTolerance> aTol =
      writeGeomTolerance(theWS, aFirstShapeL, aGeomTolL, aDatumSystem, aRC, theLocalFactors);
    if (!aTol.IsNull())
    {
      aWrittenTols.Bind(aGeomTolL, aTol);
    }
  }
  // Composite frames: each lower segment related to its upper, the upper
  // relating (AP242's geometric_tolerance_relationship 'composite')
  // In the labels' order: the same file each time
  for (NCollection_Sequence<TDF_Label>::Iterator aTolIter(aDGTLabels); aTolIter.More();
       aTolIter.Next())
  {
    TDF_Label                                  anUpperL;
    occ::handle<StepDimTol_GeometricTolerance> anUpper, aLower;
    if (!aWrittenTols.Find(aTolIter.Value(), aLower)
        || !XCAFDoc_DimTolTool::GetCompositeUpper(aTolIter.Value(), anUpperL)
        || !aWrittenTols.Find(anUpperL, anUpper))
    {
      continue;
    }
    occ::handle<StepDimTol_GeometricToleranceRelationship> aRel =
      new StepDimTol_GeometricToleranceRelationship();
    aRel->Init(new TCollection_HAsciiString("composite"),
               new TCollection_HAsciiString(),
               anUpper,
               aLower);
    aModel->AddWithRefs(aRel);
  }

  // The context: a shape aspect's, else -- none written, every item a
  // graphic of a part whose faces have no entity (a tessellated part) -- the
  // product's own. A draughting model without one cannot be written.
  for (Interface_EntityIterator aSharingIter = aGraph.Sharings(myGDTCommonPDS);
       aRC.IsNull() && !myGDTCommonPDS.IsNull() && aSharingIter.More();
       aSharingIter.Next())
  {
    occ::handle<StepShape_ShapeDefinitionRepresentation> aSDR =
      occ::down_cast<StepShape_ShapeDefinitionRepresentation>(aSharingIter.Value());
    if (!aSDR.IsNull() && !aSDR->UsedRepresentation().IsNull())
    {
      aRC = aSDR->UsedRepresentation()->ContextOfItems();
    }
  }
  if (aRC.IsNull())
  {
    return true;
  }
  for (NCollection_Sequence<occ::handle<StepRepr_Representation>>::Iterator aTextIt(aTexts);
       aTextIt.More();
       aTextIt.Next())
  {
    aTextIt.Value()->SetContextOfItems(aRC);
  }

  // The PMI presented by another's graphic: linked to it, and to its geometry, as one drawn
  // is; a view of it shows that graphic
  GDTPresentationModel& aState = gdtState(myGDTPresentationDM);
  for (NCollection_Sequence<GDTPresentationModel::Link>::Iterator aLinkIt(aState.Links);
       aLinkIt.More();
       aLinkIt.Next())
  {
    const GDTPresentationModel::Link&         aLink = aLinkIt.Value();
    // The host's callout: the first of what was written for it (its plane follows)
    occ::handle<StepVisual_DraughtingCallout> aCallout;
    if (const NCollection_Sequence<occ::handle<StepRepr_RepresentationItem>>* aShownOfHost =
          aState.Items.Seek(aLink.Host))
    {
      for (NCollection_Sequence<occ::handle<StepRepr_RepresentationItem>>::Iterator anItemIt(
             *aShownOfHost);
           anItemIt.More() && aCallout.IsNull();
           anItemIt.Next())
      {
        aCallout = occ::down_cast<StepVisual_DraughtingCallout>(anItemIt.Value());
      }
    }
    if (aCallout.IsNull())
    {
      continue;
    }
    occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> aCallouts =
      new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, 1);
    aCallouts->SetValue(1, aCallout);
    StepAP242_ItemIdentifiedRepresentationUsageDefinition aDefinition;
    aDefinition.SetValue(aLink.Dimension);
    occ::handle<StepAP242_DraughtingModelItemAssociation> aDMIA =
      new StepAP242_DraughtingModelItemAssociation();
    aDMIA->Init(new TCollection_HAsciiString("PMI representation to presentation link"),
                new TCollection_HAsciiString(),
                aDefinition,
                myGDTPresentationDM,
                aCallouts);
    aModel->AddWithRefs(aDMIA);
    for (NCollection_Sequence<occ::handle<StepRepr_ShapeAspect>>::Iterator aGeomIt(aLink.Geometry);
         aGeomIt.More();
         aGeomIt.Next())
    {
      if (aGeomIt.Value().IsNull())
      {
        continue;
      }
      StepAP242_ItemIdentifiedRepresentationUsageDefinition aGeometry;
      aGeometry.SetValue(aGeomIt.Value());
      occ::handle<StepAP242_DraughtingModelItemAssociation> aGeomDMIA =
        new StepAP242_DraughtingModelItemAssociation();
      aGeomDMIA->Init(new TCollection_HAsciiString(),
                      new TCollection_HAsciiString(),
                      aGeometry,
                      myGDTPresentationDM,
                      aCallouts);
      aModel->AddWithRefs(aGeomDMIA);
    }
    if (!aLink.Label.IsNull() && !aState.Items.IsBound(aLink.Label))
    {
      NCollection_Sequence<occ::handle<StepRepr_RepresentationItem>> aShown;
      aShown.Append(aCallout);
      aState.Items.Bind(aLink.Label, aShown);
    }
  }
  aState.Links.Clear();
  aState.Placeholders.Clear();
  // Write Draughting model for Annotation Planes
  if (myGDTAnnotations.Length() > 0)
  {
    occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> aItems =
      new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(
        1,
        myGDTAnnotations.Length());
    for (int i = 1; i <= aItems->Length(); i++)
    {
      aItems->SetValue(i, myGDTAnnotations.Value(i - 1));
    }
    myGDTPresentationDM->Init(new TCollection_HAsciiString(), aItems, aRC);
    aModel->AddWithRefs(myGDTPresentationDM);
  }

  // The presentations not shown (XCAFDoc_DimTolTool::SetPresentationHidden): their callouts'
  // occurrences listed by an invisibility, as NIST's are. One drawn in another's graphic hides
  // nothing of it
  {
    NCollection_Sequence<occ::handle<Standard_Transient>> aHidden;
    using ItemMap =
      NCollection_DataMap<TDF_Label, NCollection_Sequence<occ::handle<StepRepr_RepresentationItem>>>;
    for (ItemMap::Iterator anIter(aState.Items); anIter.More(); anIter.Next())
    {
      TDF_Label aHost;
      if (!XCAFDoc_DimTolTool::IsPresentationHidden(anIter.Key())
          || XCAFDoc_DimTolTool::GetPresentedWith(anIter.Key(), aHost))
      {
        continue;
      }
      for (const auto& anItem : anIter.Value())
      {
        occ::handle<StepVisual_DraughtingCallout> aCallout =
          occ::down_cast<StepVisual_DraughtingCallout>(anItem);
        for (int i = 1; !aCallout.IsNull() && i <= aCallout->NbContents(); ++i)
        {
          occ::handle<StepVisual_StyledItem> anOccurrence =
            occ::down_cast<StepVisual_StyledItem>(aCallout->ContentsValue(i).Value());
          if (!anOccurrence.IsNull())
          {
            aHidden.Append(anOccurrence);
          }
        }
      }
    }
    if (!aHidden.IsEmpty())
    {
      occ::handle<NCollection_HArray1<StepVisual_InvisibleItem>> anItems =
        new NCollection_HArray1<StepVisual_InvisibleItem>(1, aHidden.Length());
      for (int i = 1; i <= aHidden.Length(); ++i)
      {
        StepVisual_InvisibleItem anInvItem;
        anInvItem.SetValue(aHidden.Value(i));
        anItems->SetValue(i, anInvItem);
      }
      occ::handle<StepVisual_Invisibility> anInvisibility = new StepVisual_Invisibility();
      anInvisibility->Init(anItems);
      aModel->AddWithRefs(anInvisibility);
    }
  }

  writeViews(theWS, aRC, theLocalFactors);
  return true;
}

//=================================================================================================

void STEPCAFControl_Writer::writeViews(const occ::handle<XSControl_WorkSession>&          theWS,
                                       const occ::handle<StepRepr_RepresentationContext>& theRC,
                                       const StepData_Factors& theLocalFactors)
{
  const occ::handle<Interface_InterfaceModel>& aModel  = theWS->Model();
  const occ::handle<Transfer_FinderProcess>&   aFP     = theWS->TransferWriter()->FinderProcess();
  const double                                 aFactor = theLocalFactors.LengthFactor();
  NCollection_Map<TDF_Label>                   aTools;
  for (NCollection_Map<TDF_Label>::Iterator aRootIter(myRootLabels); aRootIter.More();
       aRootIter.Next())
  {
    occ::handle<XCAFDoc_ViewTool> aViewTool = XCAFDoc_DocumentTool::ViewTool(aRootIter.Value());
    if (aViewTool.IsNull() || !aTools.Add(aViewTool->Label()))
    {
      continue;
    }
    NCollection_Sequence<TDF_Label> aViewLabels;
    aViewTool->GetViewLabels(aViewLabels);
    for (NCollection_Sequence<TDF_Label>::Iterator aViewIter(aViewLabels); aViewIter.More();
         aViewIter.Next())
    {
      const TDF_Label&          aViewL = aViewIter.Value();
      occ::handle<XCAFDoc_View> aViewAttr;
      if (!aViewL.FindAttribute(XCAFDoc_View::GetID(), aViewAttr))
      {
        continue;
      }
      occ::handle<XCAFView_Object> anObj = aViewAttr->GetObject();
      if (anObj.IsNull() || anObj->Type() == XCAFView_ProjectionType_NoCamera)
      {
        continue;
      }
      const gp_Dir aSight = anObj->ViewDirection();
      const gp_Dir anUp   = anObj->UpDirection();
      if (anUp.IsParallel(aSight, Precision::Angular()))
      {
        continue;
      }
      occ::handle<TCollection_HAsciiString> aName =
        anObj->Name().IsNull() ? new TCollection_HAsciiString()
                               : new TCollection_HAsciiString(anObj->Name()->String());

      // The camera: its view reference system at the eye, z along the line
      // of sight, x such that z ^ x is up -- as the reader takes it back,
      // and as the NIST test files have it; the projection point the
      // system's origin
      const gp_Ax2 aFrame(anObj->ProjectionPoint(), aSight, anUp.Crossed(aSight));
      occ::handle<StepGeom_Axis2Placement3d> aSystem =
        GeomToStep_MakeAxis2Placement3d(aFrame, theLocalFactors).Value();
      // The window: in the view plane, centred on the line of sight
      const double aDistance = anObj->ZoomFactor();
      const double aSizeX    = anObj->WindowHorizontalSize();
      const double aSizeY    = anObj->WindowVerticalSize();
      const gp_Ax2 aWindowFrame(gp_Pnt(-aSizeX / 2, -aSizeY / 2, aDistance), gp::DZ(), gp::DX());
      StepGeom_Axis2Placement aWindowPlacement;
      aWindowPlacement.SetValue(
        GeomToStep_MakeAxis2Placement3d(aWindowFrame, theLocalFactors).Value());
      occ::handle<StepVisual_PlanarBox> aWindow = new StepVisual_PlanarBox();
      aWindow->Init(aName, aSizeX / aFactor, aSizeY / aFactor, aWindowPlacement);
      occ::handle<StepVisual_ViewVolume> aVolume = new StepVisual_ViewVolume();
      aVolume->Init(anObj->Type() == XCAFView_ProjectionType_Central ? StepVisual_copCentral
                                                                     : StepVisual_copParallel,
                    GeomToStep_MakeCartesianPoint(gp::Origin(), aFactor).Value(),
                    aDistance / aFactor,
                    anObj->HasFrontPlaneClipping() ? anObj->FrontPlaneDistance() / aFactor : 0.,
                    anObj->HasFrontPlaneClipping(),
                    anObj->HasBackPlaneClipping() ? anObj->BackPlaneDistance() / aFactor : 0.,
                    anObj->HasBackPlaneClipping(),
                    anObj->HasViewVolumeSidesClipping(),
                    aWindow);
      // Its clipping planes, cutting all together (the intersection of
      // what each keeps): a plain camera without, one clipping with
      NCollection_Sequence<TDF_Label> aPlaneLabels;
      aViewTool->GetRefClippingPlaneLabel(aViewL, aPlaneLabels);
      occ::handle<XCAFDoc_ClippingPlaneTool> aPlaneTool =
        XCAFDoc_DocumentTool::ClippingPlaneTool(aRootIter.Value());
      NCollection_Sequence<occ::handle<StepGeom_Plane>> aPlanes;
      for (NCollection_Sequence<TDF_Label>::Iterator aPlaneIter(aPlaneLabels); aPlaneIter.More();
           aPlaneIter.Next())
      {
        gp_Pln                                aPln;
        occ::handle<TCollection_HAsciiString> aPlaneName;
        bool                                  aCapping = false;
        if (aPlaneTool.IsNull()
            || !aPlaneTool->GetClippingPlane(aPlaneIter.Value(), aPln, aPlaneName, aCapping))
        {
          continue;
        }
        GeomToStep_MakePlane aMaker(aPln, theLocalFactors);
        if (!aMaker.IsDone())
        {
          continue;
        }
        occ::handle<StepGeom_Plane> aPlane = aMaker.Value();
        aPlane->SetName(aPlaneName.IsNull() ? new TCollection_HAsciiString()
                                            : new TCollection_HAsciiString(aPlaneName));
        aPlanes.Append(aPlane);
      }
      occ::handle<StepVisual_CameraModelD3> aCamera;
      if (aPlanes.IsEmpty())
      {
        aCamera = new StepVisual_CameraModelD3();
        aCamera->Init(aName, aSystem, aVolume);
      }
      else
      {
        occ::handle<NCollection_HArray1<StepVisual_CameraModelD3MultiClippingInterectionSelect>>
          aClipping =
            new NCollection_HArray1<StepVisual_CameraModelD3MultiClippingInterectionSelect>(
              1,
              aPlanes.Length());
        for (int i = 1; i <= aPlanes.Length(); i++)
        {
          StepVisual_CameraModelD3MultiClippingInterectionSelect aSelect;
          aSelect.SetValue(aPlanes.Value(i));
          aClipping->SetValue(i, aSelect);
        }
        occ::handle<StepVisual_CameraModelD3MultiClipping> aClippingCamera =
          new StepVisual_CameraModelD3MultiClipping();
        aClippingCamera->Init(aName, aSystem, aVolume, aClipping);
        aCamera = aClippingCamera;
      }

      // Its items: the camera, the shapes it shows -- each one's
      // representation mapped through the camera -- and the presentations
      // of the GD&T it shows
      NCollection_DynamicArray<occ::handle<StepRepr_RepresentationItem>> anItems;
      NCollection_Map<occ::handle<Standard_Transient>>                    aListed;
      anItems.Append(aCamera);
      NCollection_Sequence<TDF_Label> aShapes;
      aViewTool->GetRefShapeLabel(aViewL, aShapes);
      for (NCollection_Sequence<TDF_Label>::Iterator aShapeIter(aShapes); aShapeIter.More();
           aShapeIter.Next())
      {
        TDF_Label aShapeL = aShapeIter.Value();
        if (XCAFDoc_ShapeTool::IsReference(aShapeL))
        {
          XCAFDoc_ShapeTool::GetReferredShape(aShapeIter.Value(), aShapeL);
        }
        if (!myLabels.IsBound(aShapeL))
        {
          continue;
        }
        occ::handle<TransferBRep_ShapeMapper> aMapper =
          TransferBRep::ShapeMapper(aFP, myLabels.Find(aShapeL));
        occ::handle<StepShape_ShapeDefinitionRepresentation> aSDR;
        if (!aFP->FindTypedTransient(aMapper,
                                     STANDARD_TYPE(StepShape_ShapeDefinitionRepresentation),
                                     aSDR)
            || aSDR->UsedRepresentation().IsNull()
            || !aListed.Add(aSDR->UsedRepresentation()))
        {
          continue;
        }
        occ::handle<StepRepr_RepresentationMap> aMap = new StepRepr_RepresentationMap();
        aMap->Init(aCamera, aSDR->UsedRepresentation());
        occ::handle<StepRepr_MappedItem> aMapped = new StepRepr_MappedItem();
        aMapped->Init(new TCollection_HAsciiString(),
                      aMap,
                      GeomToStep_MakeAxis2Placement3d(theLocalFactors).Value());
        anItems.Append(aMapped);
      }
      NCollection_Sequence<TDF_Label> aGDTs;
      aViewTool->GetRefGDTLabel(aViewL, aGDTs);
      for (NCollection_Sequence<TDF_Label>::Iterator aGDTIter(aGDTs); aGDTIter.More();
           aGDTIter.Next())
      {
        const NCollection_Sequence<occ::handle<StepRepr_RepresentationItem>>* aPrs =
          gdtState(myGDTPresentationDM).Items.Seek(aGDTIter.Value());
        if (aPrs == nullptr)
        {
          continue;
        }
        for (NCollection_Sequence<occ::handle<StepRepr_RepresentationItem>>::Iterator anIter(
               *aPrs);
             anIter.More();
             anIter.Next())
        {
          if (aListed.Add(anIter.Value()))
          {
            anItems.Append(anIter.Value());
          }
        }
      }

      occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> anArray =
        new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, anItems.Length());
      for (int i = 1; i <= anArray->Length(); i++)
      {
        anArray->SetValue(i, anItems.Value(i - 1));
      }
      occ::handle<StepVisual_DraughtingModel> aModelOfView = new StepVisual_DraughtingModel();
      aModelOfView->Init(aName, anArray, theRC);
      aModel->AddWithRefs(aModelOfView);
      // In a file of AP242's second edition (one written with placeholders), the camera and
      // its draughting model are connected by a model_geometric_view, as the recommended
      // practice has it (9.4.2) and NIST's files do: SFA takes its absence for a syntax error
      // there. No class is for it: an entity of its four parameters
      if (isEdition2(aModel))
      {
        occ::handle<StepData_UndefinedEntity> aView = new StepData_UndefinedEntity();
        aView->SetStepType("MODEL_GEOMETRIC_VIEW");
        occ::handle<Interface_UndefinedContent> aParams = aView->UndefinedContent();
        // Named as its view: SFA takes a model_geometric_view of no name for a syntax error
        // (recommended practice 9.4.2, figure 96), which NIST's own files have
        TCollection_AsciiString aViewName =
          aName.IsNull() || aName->IsEmpty() ? TCollection_AsciiString("view") : aName->String();
        TCollection_AsciiString aQuoted("'");
        for (int i = 1; i <= aViewName.Length(); ++i)
        {
          const char aChar = aViewName.Value(i);
          if (aChar == '\'' || aChar == '\\')
          {
            aQuoted += aChar;
          }
          aQuoted += aChar;
        }
        aQuoted += "'";
        aParams->AddLiteral(Interface_ParamText, new TCollection_HAsciiString(aQuoted));
        aParams->AddLiteral(Interface_ParamText, new TCollection_HAsciiString("''"));
        aParams->AddEntity(Interface_ParamIdent, aCamera);
        aParams->AddEntity(Interface_ParamIdent, aModelOfView);
        aModel->AddWithRefs(aView);
      }
      // Related to the draughting model of all the GD&T's presentations
      if (myGDTAnnotations.Length() > 0)
      {
        occ::handle<StepRepr_MechanicalDesignAndDraughtingRelationship> aRelation =
          new StepRepr_MechanicalDesignAndDraughtingRelationship();
        aRelation->Init(new TCollection_HAsciiString(),
                        new TCollection_HAsciiString(),
                        aModelOfView,
                        myGDTPresentationDM);
        aModel->AddWithRefs(aRelation);
      }
    }
  }
}

//=================================================================================================

static bool FindPDSforRI(const Interface_Graph&                        theGraph,
                         const occ::handle<Standard_Transient>&        theEnt,
                         occ::handle<StepRepr_ProductDefinitionShape>& thePDS,
                         occ::handle<StepRepr_RepresentationContext>&  theRC)
{
  if (theEnt.IsNull() || !theEnt->IsKind(STANDARD_TYPE(StepRepr_RepresentationItem)))
  {
    return false;
  }
  for (Interface_EntityIterator aSharingIter = theGraph.Sharings(theEnt);
       aSharingIter.More() && thePDS.IsNull();
       aSharingIter.Next())
  {
    occ::handle<StepShape_ShapeRepresentation> aSR =
      occ::down_cast<StepShape_ShapeRepresentation>(aSharingIter.Value());
    if (aSR.IsNull())
    {
      continue;
    }
    theRC = aSR->ContextOfItems();
    for (Interface_EntityIterator aSubs1 = theGraph.Sharings(aSR); aSubs1.More() && thePDS.IsNull();
         aSubs1.Next())
    {
      occ::handle<StepShape_ShapeDefinitionRepresentation> aSDR =
        occ::down_cast<StepShape_ShapeDefinitionRepresentation>(aSubs1.Value());
      if (aSDR.IsNull())
      {
        continue;
      }
      occ::handle<StepRepr_PropertyDefinition> aPropD = aSDR->Definition().PropertyDefinition();
      if (aPropD.IsNull())
      {
        continue;
      }
      thePDS = occ::down_cast<StepRepr_ProductDefinitionShape>(aPropD);
    }
  }
  return true;
}

//=================================================================================================

bool STEPCAFControl_Writer::writeMaterials(const occ::handle<XSControl_WorkSession>& theWS,
                                           const NCollection_Sequence<TDF_Label>& theLabels) const
{

  if (theLabels.IsEmpty())
  {
    return false;
  }

  // get working data
  const occ::handle<Interface_InterfaceModel>& aModel = theWS->Model();
  const occ::handle<XSControl_TransferWriter>& aTW    = theWS->TransferWriter();
  const occ::handle<Transfer_FinderProcess>&   aFP    = aTW->FinderProcess();

  const occ::handle<Interface_HGraph> aHGraph = theWS->HGraph();
  if (aHGraph.IsNull())
  {
    return false;
  }

  Interface_Graph aGraph = theWS->HGraph()->Graph();

  NCollection_DataMap<TCollection_AsciiString, occ::handle<Standard_Transient>> aMapDRI, aMapMRI;
  NCollection_Sequence<TDF_Label>                                               aTopLabels;
  // Iterate on requested shapes collect Tools
  for (NCollection_Map<TDF_Label>::Iterator aLabelIter(myRootLabels); aLabelIter.More();
       aLabelIter.Next())
  {
    const TDF_Label&                aLabel  = aLabelIter.Value();
    occ::handle<XCAFDoc_ShapeTool>  aShTool = XCAFDoc_DocumentTool::ShapeTool(aLabel);
    NCollection_Sequence<TDF_Label> aTopInterLabels;
    aShTool->GetShapes(aTopInterLabels);
    aTopLabels.Append(aTopInterLabels);
  }
  for (NCollection_Sequence<TDF_Label>::Iterator aTopLIter(aTopLabels); aTopLIter.More();
       aTopLIter.Next())
  {
    const TDF_Label&               aShL = aTopLIter.Value();
    occ::handle<TDataStd_TreeNode> aNode;
    if (!aShL.FindAttribute(XCAFDoc::MaterialRefGUID(), aNode) || !aNode->HasFather())
    {
      continue;
    }
    // find PDS for current shape
    TopoDS_Shape    aShape = XCAFDoc_ShapeTool::GetShape(aShL);
    TopLoc_Location aLocation;
    NCollection_Sequence<occ::handle<Standard_Transient>> aSeqRI;
    FindEntities(aFP, aShape, aLocation, aSeqRI);
    if (aSeqRI.Length() <= 0)
    {
      continue;
    }
    occ::handle<StepRepr_ProductDefinitionShape> aPDS;
    occ::handle<StepRepr_RepresentationContext>  aRC;
    occ::handle<Standard_Transient>              anEnt = aSeqRI.Value(1);
    FindPDSforRI(aGraph, anEnt, aPDS, aRC);
    if (aPDS.IsNull())
    {
      continue;
    }
    occ::handle<StepBasic_ProductDefinition> aProdDef = aPDS->Definition().ProductDefinition();
    if (aProdDef.IsNull())
    {
      continue;
    }
    // write material entities
    TDF_Label                             aMatL = aNode->Father()->Label();
    occ::handle<TCollection_HAsciiString> aName;
    occ::handle<TCollection_HAsciiString> aDescription;
    double                                aDensity;
    occ::handle<TCollection_HAsciiString> aDensName;
    occ::handle<TCollection_HAsciiString> aDensValType;
    occ::handle<StepRepr_Representation>  aRepDRI;
    occ::handle<StepRepr_Representation>  aRepMRI;
    if (XCAFDoc_MaterialTool::GetMaterial(aMatL,
                                          aName,
                                          aDescription,
                                          aDensity,
                                          aDensName,
                                          aDensValType))
    {
      if (aName->Length() == 0)
      {
        continue;
      }
      TCollection_AsciiString aKey(aName->ToCString());
      if (aMapDRI.IsBound(aKey))
      {
        aRepDRI = occ::down_cast<StepRepr_Representation>(aMapDRI.Find(aKey));
        if (aMapMRI.IsBound(aKey))
        {
          aRepMRI = occ::down_cast<StepRepr_Representation>(aMapMRI.Find(aKey));
        }
      }
      else
      {
        // write DRI
        occ::handle<StepRepr_DescriptiveRepresentationItem> aDRI =
          new StepRepr_DescriptiveRepresentationItem;
        aDRI->Init(aName, aDescription);
        occ::handle<NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>> aHARI =
          new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, 1);
        aHARI->SetValue(1, aDRI);
        aRepDRI = new StepRepr_Representation();
        aRepDRI->Init(new TCollection_HAsciiString("material name"), aHARI, aRC);
        aModel->AddWithRefs(aRepDRI);
        // write MRI
        if (aDensity > 0)
        {
          // mass
          occ::handle<StepBasic_SiUnitAndMassUnit> aSMU = new StepBasic_SiUnitAndMassUnit;
          aSMU->SetName(StepBasic_sunGram);
          occ::handle<StepBasic_DerivedUnitElement> aDUE1 = new StepBasic_DerivedUnitElement;
          aDUE1->Init(aSMU, 3.0);
          // length
          occ::handle<StepBasic_SiUnitAndLengthUnit> aSLU = new StepBasic_SiUnitAndLengthUnit;
          aSLU->Init(true, StepBasic_spCenti, StepBasic_sunMetre);
          occ::handle<StepBasic_DerivedUnitElement> aDUE2 = new StepBasic_DerivedUnitElement;
          aDUE2->Init(aSLU, 2.0);
          // other
          occ::handle<NCollection_HArray1<occ::handle<StepBasic_DerivedUnitElement>>> aHADUE =
            new NCollection_HArray1<occ::handle<StepBasic_DerivedUnitElement>>(1, 2);
          aHADUE->SetValue(1, aDUE1);
          aHADUE->SetValue(2, aDUE2);
          occ::handle<StepBasic_DerivedUnit> aDU = new StepBasic_DerivedUnit;
          aDU->Init(aHADUE);
          aModel->AddWithRefs(aDU);
          StepBasic_Unit aUnit;
          aUnit.SetValue(aDU);
          occ::handle<StepBasic_MeasureValueMember> aMVM = new StepBasic_MeasureValueMember;
          aMVM->SetName(aDensValType->ToCString());
          aMVM->SetReal(aDensity);
          occ::handle<StepRepr_MeasureRepresentationItem> aMRI =
            new StepRepr_MeasureRepresentationItem;
          aMRI->Init(aDensName, aMVM, aUnit);
          aHARI = new NCollection_HArray1<occ::handle<StepRepr_RepresentationItem>>(1, 1);
          aHARI->SetValue(1, aMRI);
          aRepMRI = new StepRepr_Representation();
          aRepMRI->Init(new TCollection_HAsciiString("density"), aHARI, aRC);
          aModel->AddWithRefs(aRepMRI);
        }
        // WriteNewMaterial(Model,aName,aDescription,aDensity,aDensName,aDensValType,RC,RepDRI,RepMRI);
        aMapDRI.Bind(aKey, aRepDRI);
        if (!aRepMRI.IsNull())
        {
          aMapMRI.Bind(aKey, aRepMRI);
        }
      }
    }

    if (!aRepDRI.IsNull())
    {
      StepRepr_CharacterizedDefinition aCD1;
      aCD1.SetValue(aProdDef);
      occ::handle<StepRepr_PropertyDefinition> aPropD1 = new StepRepr_PropertyDefinition;
      aPropD1->Init(new TCollection_HAsciiString("material property"),
                    true,
                    new TCollection_HAsciiString("material name"),
                    aCD1);
      aModel->AddWithRefs(aPropD1);
      StepRepr_RepresentedDefinition aRD1;
      aRD1.SetValue(aPropD1);
      occ::handle<StepRepr_PropertyDefinitionRepresentation> aPDR1 =
        new StepRepr_PropertyDefinitionRepresentation;
      aPDR1->Init(aRD1, aRepDRI);
      aModel->AddWithRefs(aPDR1);

      if (!aRepMRI.IsNull())
      {
        StepRepr_CharacterizedDefinition aCD2;
        aCD2.SetValue(aProdDef);
        occ::handle<StepRepr_PropertyDefinition> aPropD2 = new StepRepr_PropertyDefinition;
        aPropD2->Init(new TCollection_HAsciiString("material property"),
                      true,
                      new TCollection_HAsciiString("density"),
                      aCD2);
        aModel->AddWithRefs(aPropD2);
        StepRepr_RepresentedDefinition aRD2;
        aRD2.SetValue(aPropD2);
        occ::handle<StepRepr_PropertyDefinitionRepresentation> aPDR2 =
          new StepRepr_PropertyDefinitionRepresentation;
        aPDR2->Init(aRD2, aRepMRI);
        aModel->AddWithRefs(aPDR2);
      }
    }
  }

  return true;
}
