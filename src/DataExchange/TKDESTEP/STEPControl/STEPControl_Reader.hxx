// Created on: 1996-07-08
// Created by: Christian CAILLET
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

#ifndef _STEPControl_Reader_HeaderFile
#define _STEPControl_Reader_HeaderFile

#include <Standard.hxx>
#include <Standard_DefineAlloc.hxx>
#include <Standard_Handle.hxx>

#include <XSControl_Reader.hxx>
#include <Standard_Integer.hxx>
#include <DESTEP_Parameters.hxx>
#include <TCollection_AsciiString.hxx>
#include <NCollection_HSequence.hxx>
#include <NCollection_Sequence.hxx>
#include <NCollection_Array1.hxx>
class XSControl_WorkSession;
class StepData_StepModel;
class StepRepr_RepresentationContext;

//! Reads STEP files, checks them and translates their contents
//! into Open CASCADE models. The STEP data can be that of
//! a whole model or that of a specific list of entities in the model.
//! As in XSControl_Reader, you specify the list using a selection.
//! For the translation of iges files it is possible to use next sequence:
//! To change translation parameters
//! class Interface_Static should be used before beginning of
//! translation (see STEP Parameters and General Parameters)
//! Creation of reader - STEPControl_Reader reader;
//! To load s file in a model use method reader.ReadFile("filename.stp")
//! To print load results reader.PrintCheckLoad(failsonly,mode)
//! where mode is equal to the value of enumeration IFSelect_PrintCount
//! For definition number of candidates :
//! int nbroots = reader. NbRootsForTransfer();
//! To transfer entities from a model the following methods can be used:
//! for the whole model - reader.TransferRoots();
//! to transfer a list of entities: reader.TransferList(list);
//! to transfer one entity occ::handle<Standard_Transient>
//! ent = reader.RootForTransfer(num);
//! reader.TransferEntity(ent), or
//! reader.TransferOneRoot(num), or
//! reader.TransferOne(num), or
//! reader.TransferRoot(num)
//! To obtain the result the following method can be used:
//! reader.NbShapes() and reader.Shape(num); or reader.OneShape();
//! To print the results of transfer use method:
//! reader.PrintCheckTransfer(failwarn,mode);
//! where printfail is equal to the value of enumeration
//! IFSelect_PrintFail, mode see above; or reader.PrintStatsTransfer();
//! Gets correspondence between a STEP entity and a result
//! shape obtained from it.
//! occ::handle<XSControl_WorkSession>
//! WS = reader.WS();
//! if ( WS->TransferReader()->HasResult(ent) )
//! TopoDS_Shape shape = WS->TransferReader()->ShapeResult(ent);
class STEPControl_Reader : public XSControl_Reader
{
public:
  DEFINE_STANDARD_ALLOC

  //! Creates a reader object with an empty STEP model.
  Standard_EXPORT STEPControl_Reader();

  //! Creates a Reader for STEP from an already existing Session
  //! Clears the session if it was not yet set for STEP
  Standard_EXPORT STEPControl_Reader(const occ::handle<XSControl_WorkSession>& WS,
                                     const bool                                scratch = true);

  //! Returns the model as a StepModel.
  //! It can then be consulted (header, product)
  Standard_EXPORT occ::handle<StepData_StepModel> StepModel() const;

  //! Loads a file and returns the read status
  //! Zero for a Model which compies with the Controller
  Standard_EXPORT IFSelect_ReturnStatus ReadFile(const char* const filename) override;

  //! Loads a file from stream and returns the read status
  Standard_EXPORT IFSelect_ReturnStatus ReadStream(const char* const theName,
                                                   std::istream&     theIStream) override;

  //! Loads a file and returns the read status
  //! Zero for a Model which compies with the Controller
  Standard_EXPORT IFSelect_ReturnStatus ReadFile(const char* const        filename,
                                                 const DESTEP_Parameters& theParams);

  //! Loads a file from stream and returns the read status
  Standard_EXPORT IFSelect_ReturnStatus ReadStream(const char* const        theName,
                                                   const DESTEP_Parameters& theParams,
                                                   std::istream&            theIStream);

  //! Transfers a root given its rank in the list of candidate roots
  //! Default is the first one
  //! Returns True if a shape has resulted, false else
  //! Same as inherited TransferOneRoot, kept for compatibility
  Standard_EXPORT bool TransferRoot(
    const int                    num         = 1,
    const Message_ProgressRange& theProgress = Message_ProgressRange());

  //! Determines the list of root entities from Model which are candidate for
  //! a transfer to a Shape (type of entities is PRODUCT)
  Standard_EXPORT int NbRootsForTransfer() override;

  //! Returns sequence of all unit names for shape representations
  //! found in file
  Standard_EXPORT void FileUnits(
    NCollection_Sequence<TCollection_AsciiString>& theUnitLengthNames,
    NCollection_Sequence<TCollection_AsciiString>& theUnitAngleNames,
    NCollection_Sequence<TCollection_AsciiString>& theUnitSolidAngleNames);

  //! Sets system length unit used by transfer process.
  //! Performs only if a model is not NULL
  Standard_EXPORT void SetSystemLengthUnit(const double theLengthUnit);

  //! Returns system length unit used by transfer process.
  //! Performs only if a model is not NULL
  Standard_EXPORT double SystemLengthUnit() const;

  //! Same as TransferRoots(), but runs the whole batch with the actor's deferred
  //! post-processing enabled: per-shape work such as shape healing is accumulated
  //! during translation and flushed once at the end (possibly in parallel), and
  //! result shapes are extracted only after the flush. Falls back to plain
  //! sequential behavior with an actor that does not support deferral.
  //! @param theFirst rank of the first root to transfer (1-based); values below 1
  //!        are clamped to 1.
  //! @param theLast rank of the last root to transfer; 0 (default) or values
  //!        beyond the root count mean "up to the last root".
  //! Warning - This function clears existing output shapes first.
  Standard_EXPORT int TransferRootsDeferred(
    const Message_ProgressRange& theProgress = Message_ProgressRange(),
    const int                    theFirst    = 1,
    const int                    theLast     = 0);

  //! Same as TransferRootsDeferred(), but for an arbitrary list of entities
  //! instead of a range of roots - a progressive import can hand over the
  //! components of one root this way, since a component entity translated on
  //! its own binds its result in the transfer process and the owning root
  //! reuses that binder instead of translating the component again.
  //! Warning - This function clears existing output shapes first.
  Standard_EXPORT int TransferListDeferred(
    const occ::handle<NCollection_HSequence<occ::handle<Standard_Transient>>>& theList,
    const Message_ProgressRange& theProgress = Message_ProgressRange());

protected:
  //! Returns default parameters for shape fixing.
  //! This method is used by the base class to get default parameters for shape fixing.
  //! @return default parameters for shape fixing.
  Standard_EXPORT DE_ShapeFixParameters GetDefaultShapeFixParameters() const override;

  //! Returns default flags for shape processing.
  //! @return Default flags for shape processing.
  Standard_EXPORT ShapeProcess::OperationsFlags GetDefaultShapeProcessFlags() const override;

private:
  //! Translates the given entities with the actor's deferred post-processing
  //! enabled, flushes it once, and only then records and collects the results.
  //! Shared implementation of TransferRootsDeferred() and TransferListDeferred().
  //! Lives here rather than on XSControl_Reader because only the STEP actor
  //! supports deferral, which keeps the generic base classes untouched.
  //! @param theScopeName progress scope name of the translation half, may be null.
  int transferDeferred(const NCollection_Sequence<occ::handle<Standard_Transient>>& theEntities,
                       const char*                                                 theScopeName,
                       const Message_ProgressRange&                                theProgress);

  //! Returns units for length, angle and solidangle for shape representations
  Standard_EXPORT bool findUnits(const occ::handle<StepRepr_RepresentationContext>& theReprContext,
                                 NCollection_Array1<TCollection_AsciiString>&       theNameUnits,
                                 NCollection_Array1<double>&                        theFactorUnits);
};

#endif // _STEPControl_Reader_HeaderFile
