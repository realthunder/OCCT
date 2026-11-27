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

#ifndef _STEPCAFControl_Writer_HeaderFile
#define _STEPCAFControl_Writer_HeaderFile

#include <TopoDS_Shape.hxx>
#include <Quantity_Color.hxx>
#include <Standard_Transient.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <NCollection_DataMap.hxx>
#include <NCollection_DynamicArray.hxx>
#include <Standard_DefineAlloc.hxx>
#include <Standard_Handle.hxx>
#include <TDF_Label.hxx>
#include <STEPCAFControl_ExternFile.hxx>
#include <STEPControl_Writer.hxx>
#include <StepAP242_GeometricItemSpecificUsage.hxx>
#include <DESTEP_Parameters.hxx>
#include <StepData_Factors.hxx>
#include <StepDimTol_Datum.hxx>
#include <StepDimTol_GeometricTolerance.hxx>
#include <StepDimTol_DatumSystemOrReference.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_HArray1.hxx>
#include <StepRepr_ProductDefinitionShape.hxx>
#include <StepVisual_DraughtingModel.hxx>
#include <StepVisual_PresentationStyleAssignment.hxx>
#include <NCollection_Sequence.hxx>
#include <NCollection_Map.hxx>
#include <XCAFDimTolObjects_GeomToleranceObject.hxx>

class XSControl_WorkSession;
class TDocStd_Document;
class STEPCAFControl_ExternFile;
class TopoDS_Shape;
class StepShape_ShapeDefinitionRepresentation;

//! Provides a tool to write DECAF document to the
//! STEP file. Besides transfer of shapes (including
//! assemblies) provided by STEPControl, supports also
//! colors and part names
//!
//! Also supports multifile writing
class STEPCAFControl_Writer
{
public:
  DEFINE_STANDARD_ALLOC

public:
  //! Creates a writer with an empty
  //! STEP model and sets ColorMode, LayerMode, NameMode and
  //! PropsMode to true.
  Standard_EXPORT STEPCAFControl_Writer();

  //! Creates a reader tool and attaches it to an already existing Session
  //! Clears the session if it was not yet set for STEP
  //! Clears the internal data structures
  Standard_EXPORT STEPCAFControl_Writer(const occ::handle<XSControl_WorkSession>& theWS,
                                        const bool theScratch = true);

  //! Clears the internal data structures and attaches to a new session
  //! Clears the session if it was not yet set for STEP
  Standard_EXPORT void Init(const occ::handle<XSControl_WorkSession>& theWS,
                            const bool                                theScratch = true);

  //! Writes all the produced models into file
  //! In case of multimodel with extern references,
  //! filename will be a name of root file, all other files
  //! have names of corresponding parts
  //! Provided for use like single-file writer
  Standard_EXPORT IFSelect_ReturnStatus Write(const char* const theFileName);

  //! Writes all the produced models into the stream.
  //! Provided for use like single-file writer
  Standard_EXPORT IFSelect_ReturnStatus WriteStream(std::ostream& theStream);

  //! Transfers a document (or single label) to a STEP model
  //! The mode of translation of shape is AsIs
  //! If multi is not null pointer, it switches to multifile
  //! mode (with external refs), and string pointed by <multi>
  //! gives prefix for names of extern files (can be empty string)
  //! Returns True if translation is OK
  Standard_EXPORT bool Transfer(const occ::handle<TDocStd_Document>& theDoc,
                                const STEPControl_StepModelType      theMode    = STEPControl_AsIs,
                                const char* const                    theIsMulti = nullptr,
                                const Message_ProgressRange& theProgress = Message_ProgressRange());

  //! Transfers a document (or single label) to a STEP model
  //! This method uses if need to set parameters avoiding
  //! initialization from Interface_Static
  //! @param theParams  configuration parameters
  //! @param theMode    mode of translation of shape is AsIs
  //! @param theIsMulti if multi is not null pointer, it switches to multifile
  //!                   mode (with external refs), and string pointed by <multi>
  //!                   gives prefix for names of extern files (can be empty string)
  //! @param theProgress progress indicator
  //! Returns True if translation is OK
  Standard_EXPORT bool Transfer(const occ::handle<TDocStd_Document>& theDoc,
                                const DESTEP_Parameters&             theParams,
                                const STEPControl_StepModelType      theMode    = STEPControl_AsIs,
                                const char* const                    theIsMulti = nullptr,
                                const Message_ProgressRange& theProgress = Message_ProgressRange());

  //! Method to transfer part of the document specified by label
  Standard_EXPORT bool Transfer(const TDF_Label&                theLabel,
                                const STEPControl_StepModelType theMode    = STEPControl_AsIs,
                                const char* const               theIsMulti = nullptr,
                                const Message_ProgressRange& theProgress = Message_ProgressRange());

  //! Method to transfer part of the document specified by label
  //! This method uses if need to set parameters avoiding
  //! initialization from Interface_Static
  Standard_EXPORT bool Transfer(const TDF_Label&                theLabel,
                                const DESTEP_Parameters&        theParams,
                                const STEPControl_StepModelType theMode    = STEPControl_AsIs,
                                const char* const               theIsMulti = nullptr,
                                const Message_ProgressRange& theProgress = Message_ProgressRange());

  //! Method to writing sequence of root assemblies
  //! or part of the file specified by use by one label
  Standard_EXPORT bool Transfer(const NCollection_Sequence<TDF_Label>& theLabelSeq,
                                const STEPControl_StepModelType        theMode = STEPControl_AsIs,
                                const char* const                      theIsMulti = nullptr,
                                const Message_ProgressRange& theProgress = Message_ProgressRange());

  //! Method to writing sequence of root assemblies
  //! or part of the file specified by use by one label.
  //! This method is utilized if there's a need to set parameters avoiding
  //! initialization from Interface_Static
  Standard_EXPORT bool Transfer(const NCollection_Sequence<TDF_Label>& theLabelSeq,
                                const DESTEP_Parameters&               theParams,
                                const STEPControl_StepModelType        theMode = STEPControl_AsIs,
                                const char* const                      theIsMulti = nullptr,
                                const Message_ProgressRange& theProgress = Message_ProgressRange());

  Standard_EXPORT bool Perform(const occ::handle<TDocStd_Document>& theDoc,
                               const TCollection_AsciiString&       theFileName,
                               const Message_ProgressRange& theProgress = Message_ProgressRange());

  //! Transfers a document and writes it to a STEP file
  //! Returns True if translation is OK
  Standard_EXPORT bool Perform(const occ::handle<TDocStd_Document>& theDoc,
                               const char* const                    theFileName,
                               const Message_ProgressRange& theProgress = Message_ProgressRange());

  //! Transfers a document and writes it to a STEP file
  //! This method is utilized if there's a need to set parameters avoiding
  //! initialization from Interface_Static
  //! Returns True if translation is OK
  Standard_EXPORT bool Perform(const occ::handle<TDocStd_Document>& theDoc,
                               const char* const                    theFileName,
                               const DESTEP_Parameters&             theParams,
                               const Message_ProgressRange& theProgress = Message_ProgressRange());

  //! Returns data on external files
  //! Returns Null handle if no external files are read
  const NCollection_DataMap<TCollection_AsciiString, occ::handle<STEPCAFControl_ExternFile>>&
    ExternFiles() const
  {
    return myFiles;
  };

  //! Returns data on external file by its original label
  //! Returns False if no external file with given name is read
  Standard_EXPORT bool ExternFile(const TDF_Label&                        theLabel,
                                  occ::handle<STEPCAFControl_ExternFile>& theExtFile) const;

  //! Returns data on external file by its name
  //! Returns False if no external file with given name is read
  Standard_EXPORT bool ExternFile(const char* const                       theName,
                                  occ::handle<STEPCAFControl_ExternFile>& theExtFile) const;

  //! Returns basic reader for root file
  STEPControl_Writer& ChangeWriter() { return myWriter; }

  //! Returns basic reader as const
  const STEPControl_Writer& Writer() const { return myWriter; }

  //! Set ColorMode for indicate write Colors or not.
  void SetColorMode(const bool theColorMode) { myColorMode = theColorMode; }

  bool GetColorMode() const { return myColorMode; }

  //! Set NameMode for indicate write Name or not.
  void SetNameMode(const bool theNameMode) { myNameMode = theNameMode; }

  bool GetNameMode() const { return myNameMode; }

  //! Set LayerMode for indicate write Layers or not.
  void SetLayerMode(const bool theLayerMode) { myLayerMode = theLayerMode; }

  bool GetLayerMode() const { return myLayerMode; }

  //! PropsMode for indicate write Validation properties or not.
  void SetPropsMode(const bool thePropsMode) { myPropsMode = thePropsMode; }

  bool GetPropsMode() const { return myPropsMode; }

  //! Set MetadataMode for indicate write metadata or not.
  void SetMetadataMode(const bool theMetadataMode) { myMetadataMode = theMetadataMode; }

  bool GetMetadataMode() const { return myMetadataMode; }

  //! Set SHUO mode for indicate write SHUO or not.
  void SetSHUOMode(const bool theSHUOMode) { mySHUOMode = theSHUOMode; }

  bool GetSHUOMode() const { return mySHUOMode; }

  //! Set dimtolmode for indicate write D&GTs or not.
  void SetDimTolMode(const bool theDimTolMode) { myGDTMode = theDimTolMode; };

  bool GetDimTolMode() const { return myGDTMode; }

  //! Set the colour GD&T presentations are drawn in (AP242): white unless
  //! set, as it always was.
  //! Not inline, and kept with no data member of its own: this class is
  //! constructed on its user's stack, so its size and layout are its ABI --
  //! a binary built against an earlier header must run with this library.
  Standard_EXPORT void SetGDTPresentationColor(const Quantity_Color& theColor);

  Standard_EXPORT const Quantity_Color& GetGDTPresentationColor() const;

  //! Set flag for indicate write material or not.
  void SetMaterialMode(const bool theMaterialMode) { myMatMode = theMaterialMode; }

  bool GetMaterialMode() const { return myMatMode; }

  //! Set flag for indicate write visual material or not.
  void SetVisualMaterialMode(const bool theVisualMaterialMode)
  {
    myVisMatMode = theVisualMaterialMode;
  }

  bool GetVisualMaterialMode() const { return myVisMatMode; }

  //! Set clean duplicates flag.
  //! If set to True, duplicates will be removed from the model.
  //! @param theCleanDuplicates the flag to set.
  void SetCleanDuplicates(const bool theCleanDuplicates)
  {
    myIsCleanDuplicates = theCleanDuplicates;
  }

  //! Returns the flag indicating whether duplicates should be removed from the model.
  //! @return the flag indicating whether duplicates should be removed from the model.
  bool GetCleanDuplicates() const { return myIsCleanDuplicates; }

  //! Sets parameters for shape processing.
  //! @param theParameters the parameters for shape processing.
  Standard_EXPORT void SetShapeFixParameters(
    const XSAlgo_ShapeProcessor::ParameterMap& theParameters);

  //! Sets parameters for shape processing.
  //! Parameters are moved from the input map.
  //! @param theParameters the parameters for shape processing.
  Standard_EXPORT void SetShapeFixParameters(XSAlgo_ShapeProcessor::ParameterMap&& theParameters);

  //! Sets parameters for shape processing.
  //! Parameters from @p theParameters are copied to the internal map.
  //! Parameters from @p theAdditionalParameters are copied to the internal map
  //! if they are not present in @p theParameters.
  //! @param theParameters the parameters for shape processing.
  //! @param theAdditionalParameters the additional parameters for shape processing.
  Standard_EXPORT void SetShapeFixParameters(
    const DE_ShapeFixParameters&               theParameters,
    const XSAlgo_ShapeProcessor::ParameterMap& theAdditionalParameters = {});

  //! Returns parameters for shape processing that was set by SetParameters() method.
  //! @return the parameters for shape processing. Empty map if no parameters were set.
  Standard_EXPORT const XSAlgo_ShapeProcessor::ParameterMap& GetShapeFixParameters() const;

  //! Sets flags defining operations to be performed on shapes.
  //! @param theFlags The flags defining operations to be performed on shapes.
  Standard_EXPORT void SetShapeProcessFlags(const ShapeProcess::OperationsFlags& theFlags);

  //! Returns flags defining operations to be performed on shapes.
  //! @return Pair of values defining operations to be performed on shapes and a boolean value
  //!         that indicates whether the flags were set.
  Standard_EXPORT const XSAlgo_ShapeProcessor::ProcessingFlags& GetShapeProcessFlags() const;

protected:
  //! Transfers labels to a STEP model
  //! Returns True if translation is OK
  //! isExternFile setting from transferExternFiles method
  bool transfer(STEPControl_Writer&                    theWriter,
                const NCollection_Sequence<TDF_Label>& theLabels,
                const STEPControl_StepModelType        theMode      = STEPControl_AsIs,
                const char* const                      theIsMulti   = nullptr,
                const bool                             isExternFile = false,
                const Message_ProgressRange&           theProgress  = Message_ProgressRange());

  //! Parses assembly structure of label L, writes all the simple
  //! shapes each to its own file named by name of its label plus
  //! prefix
  //! Returns shape representing that assembly structure
  //! in the form of nested empty compounds (and a sequence of
  //! labels which are newly written nodes of this assembly)
  TopoDS_Shape transferExternFiles(
    const TDF_Label&                 theLabel,
    const STEPControl_StepModelType  theMode,
    NCollection_Sequence<TDF_Label>& theLabelSeq,
    const StepData_Factors&          theLocalFactors = StepData_Factors(),
    const char* const                thePrefix       = "",
    const Message_ProgressRange&     theProgress     = Message_ProgressRange());

  //! Write external references to STEP
  bool writeExternRefs(const occ::handle<XSControl_WorkSession>& theWS,
                       const NCollection_Sequence<TDF_Label>&    theLabels) const;

  //! Write colors assigned to specified labels, to STEP model
  bool writeColors(const occ::handle<XSControl_WorkSession>& theWS,
                   const NCollection_Sequence<TDF_Label>&    theLabels);

  //! Write names assigned to specified labels, to STEP model
  bool writeNames(const occ::handle<XSControl_WorkSession>& theWS,
                  const NCollection_Sequence<TDF_Label>&    theLabels) const;

  //! Write metadata assigned to specified labels, to STEP model
  bool writeMetadata(const occ::handle<XSControl_WorkSession>& theWS,
                     const NCollection_Sequence<TDF_Label>&    theLabels) const;

  //! Write metadata assigned to specified label, to STEP model.
  // Also recursively writes metadata for children labels.
  bool writeMetadataForLabel(const occ::handle<XSControl_WorkSession>& theWS,
                             const TDF_Label&                          theLabel) const;

  //! Write metadata representation item to STEP model.
  //! @param theKey The key for metadata item.
  //! @param theModel The STEP model to write to.
  //! @param theShapeDefRep The shape definition representation.
  //! @param theProdDef The product definition.
  //! @param theItem The representation item to write.
  void writeMetadataRepresentationItem(
    const TCollection_AsciiString&                              theKey,
    const occ::handle<StepData_StepModel>&                      theModel,
    const occ::handle<StepShape_ShapeDefinitionRepresentation>& theShapeDefRep,
    const occ::handle<StepBasic_ProductDefinition>&             theProdDef,
    const occ::handle<StepRepr_RepresentationItem>&             theItem) const;

  //! Write D&GTs assigned to specified labels, to STEP model
  bool writeDGTs(const occ::handle<XSControl_WorkSession>& theWS,
                 const NCollection_Sequence<TDF_Label>&    theLabels) const;

  //! Write D&GTs assigned to specified labels, to STEP model, according AP242
  bool writeDGTsAP242(const occ::handle<XSControl_WorkSession>& theWS,
                      const NCollection_Sequence<TDF_Label>&    theLabels,
                      const StepData_Factors& theLocalFactors = StepData_Factors());

  //! Write materials assigned to specified labels, to STEP model
  bool writeMaterials(const occ::handle<XSControl_WorkSession>& theWS,
                      const NCollection_Sequence<TDF_Label>&    theLabels) const;

  //! Write validation properties assigned to specified labels,
  //! to STEP model
  bool writeValProps(const occ::handle<XSControl_WorkSession>& theWS,
                     const NCollection_Sequence<TDF_Label>&    theLabels,
                     const char* const                         theIsMulti) const;

  //! Write the attribute validation properties a label and its children were read with
  //! (XCAFDoc::ValidationItemsGUID): the counts of the attributes written of each
  void writeAttributeValidation(const occ::handle<XSControl_WorkSession>& theWS,
                                const TDF_Label&                          theLabel) const;

  //! Write the pmi validation properties the GD&T, its graphics, the views and the parts were
  //! read with (XCAFDoc::ValidationItemsGUID), computed of what was written
  void writePMIValidation(const occ::handle<XSControl_WorkSession>&          theWS,
                          const occ::handle<StepRepr_RepresentationContext>& theRC,
                          const StepData_Factors&                            theLocalFactors) const;

  //! Write layers assigned to specified labels, to STEP model
  bool writeLayers(const occ::handle<XSControl_WorkSession>& theWS,
                   const NCollection_Sequence<TDF_Label>&    theLabels) const;

  //! Write SHUO assigned to specified component, to STEP model
  bool writeSHUOs(const occ::handle<XSControl_WorkSession>& theWS,
                  const NCollection_Sequence<TDF_Label>&    theLabels);

  //! Finds length units located in root of label
  //! If it exists, initializes local length unit from it
  //! Else initializes according to Cascade length unit
  void prepareUnit(const TDF_Label&                       theLabel,
                   const occ::handle<StepData_StepModel>& theModel,
                   StepData_Factors&                      theLocalFactors);

  occ::handle<StepRepr_ShapeAspect> writeShapeAspect(
    const occ::handle<XSControl_WorkSession>&          theWS,
    const TDF_Label                                    theLabel,
    const TopoDS_Shape&                                theShape,
    occ::handle<StepRepr_RepresentationContext>&       theRC,
    occ::handle<StepAP242_GeometricItemSpecificUsage>& theGISU);

  void writePresentation(const occ::handle<XSControl_WorkSession>&    theWS,
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
                         const StepData_Factors& theLocalFactors = StepData_Factors(),
                         const bool              theAlone        = false);

  //! The saved views of the documents written (XCAFDoc_ViewTool): a
  //! draughting model each, with its camera, the presentations of the GD&T
  //! it shows and the shapes it shows mapped into it, as the CAx-IF practice
  //! and the reader have them. Called once the GD&T is written.
  void writeViews(const occ::handle<XSControl_WorkSession>&          theWS,
                  const occ::handle<StepRepr_RepresentationContext>& theRC,
                  const StepData_Factors&                            theLocalFactors);

  occ::handle<StepDimTol_Datum> writeDatumAP242(
    const occ::handle<XSControl_WorkSession>& theWS,
    const NCollection_Sequence<TDF_Label>&    theShapeL,
    const TDF_Label&                          theDatumL,
    const bool                                isFirstDTarget,
    const occ::handle<StepDimTol_Datum>&      theWrittenDatum,
    const StepData_Factors&                   theLocalFactors = StepData_Factors());

  void writeToleranceZone(const occ::handle<XSControl_WorkSession>&                 theWS,
                          const occ::handle<XCAFDimTolObjects_GeomToleranceObject>& theObject,
                          const occ::handle<StepDimTol_GeometricTolerance>&         theEntity,
                          const occ::handle<StepRepr_RepresentationContext>&        theRC,
                          const StepData_Factors& theLocalFactors = StepData_Factors());

  //! The tolerance written, null if none
  occ::handle<StepDimTol_GeometricTolerance> writeGeomTolerance(
    const occ::handle<XSControl_WorkSession>&                                  theWS,
    const NCollection_Sequence<TDF_Label>&                                     theShapeSeqL,
    const TDF_Label&                                                           theGeomTolL,
    const occ::handle<NCollection_HArray1<StepDimTol_DatumSystemOrReference>>& theDatumSystem,
    const occ::handle<StepRepr_RepresentationContext>&                         theRC,
    const StepData_Factors& theLocalFactors = StepData_Factors());

private:
  STEPControl_Writer                                                                   myWriter;
  NCollection_DataMap<TCollection_AsciiString, occ::handle<STEPCAFControl_ExternFile>> myFiles;
  NCollection_Map<TDF_Label>                                                           myRootLabels;
  NCollection_DataMap<TDF_Label, TopoDS_Shape>                                         myLabels;
  NCollection_DataMap<TDF_Label, occ::handle<STEPCAFControl_ExternFile>>               myLabEF;
  NCollection_DataMap<TDF_Label, TopoDS_Shape> myPureRefLabels;
  bool                                         myColorMode;
  bool                                         myNameMode;
  bool                                         myLayerMode;
  bool                                         myPropsMode;
  bool                                         myMetadataMode;
  bool                                         mySHUOMode;
  NCollection_DataMap<TopoDS_Shape, occ::handle<Standard_Transient>, TopTools_ShapeMapHasher>
                                                                     myMapCompMDGPR;
  bool                                                               myGDTMode;
  bool                                                               myMatMode;
  bool                                                               myVisMatMode;
  bool                                                               myIsCleanDuplicates;
  NCollection_DynamicArray<occ::handle<StepRepr_RepresentationItem>> myGDTAnnotations;
  //! The draughting model of the GD&T presentations. It carries, out of
  //! this class's layout, what was added to the writer since (see
  //! SetGDTPresentationColor): the colour, and the presentation of each
  //! GD&T label written. Add no data member here: hold new state there.
  occ::handle<StepVisual_DraughtingModel>                            myGDTPresentationDM;
  occ::handle<NCollection_HArray1<occ::handle<StepVisual_PresentationStyleAssignment>>>
                                               myGDTPrsCurveStyle;
  occ::handle<StepRepr_ProductDefinitionShape> myGDTCommonPDS;
};

//! SetGDTPresentationColor as a C function, for a client to look up by name at run time and
//! do without (GD&T drawn white) where the library is upstream's, which has none: a client
//! calling it by its C++ name does not load there.
extern "C" Standard_EXPORT void STEPCAFControl_Writer_SetGDTPresentationColor(
  STEPCAFControl_Writer& theWriter,
  const Quantity_Color&  theColor);

#endif // _STEPCAFControl_Writer_HeaderFile
