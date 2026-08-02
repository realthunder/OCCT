// Created on: 1991-08-30
// Created by: Christian CAILLET
// Copyright (c) 1991-1999 Matra Datavision
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

#include <StepFile_Read.hxx>

#include <StepFile_ReadData.hxx>

#include <Interface_Check.hxx>
#include <Interface_InterfaceError.hxx>
#include <Interface_ParamType.hxx>
#include <Interface_Protocol.hxx>

#include <StepData_FileRecognizer.hxx>
#include <StepData_Protocol.hxx>
#include <StepData_StepModel.hxx>
#include <StepData_StepReaderData.hxx>
#include <StepData_StepReaderTool.hxx>

#include <Standard_ErrorHandler.hxx>
#include <Standard_Failure.hxx>

#include <Interface_Static.hxx>

#include <Message.hxx>
#include <Message_Messenger.hxx>

#include <OSD_FileSystem.hxx>
#include <OSD_Parallel.hxx>
#include <OSD_Timer.hxx>

#include "step.tab.hxx"

#include <algorithm>
#include <cstdio>
#include <memory>
#include <mutex>
#include <streambuf>
#include <vector>

namespace
{
static std::mutex& GetGlobalReadMutex()
{
  static std::mutex THE_GLOBAL_READ_MUTEX;
  return THE_GLOBAL_READ_MUTEX;
}

//! Chunking a small file costs more in thread hand-off than it saves, and past
//! a point more chunks only add merge work.
static const int64_t THE_MIN_CHUNKED_SIZE = 4 * 1024 * 1024;
static const int     THE_MAX_CHUNKS       = 16;

//! Stands in for the part of the file a chunk does not contain: the grammar
//! accepts an empty header (rule stepf2), so this is a complete STEP file
//! around any run of DATA records.
static const char THE_CHUNK_PREAMBLE[]  = "ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n";
static const char THE_CHUNK_POSTAMBLE[] = "\nENDSEC;\nEND-ISO-10303-21;\n";

//! A slice of the DATA section, in bytes, half-open: [Start, End).
struct StepFile_ChunkRange
{
  int64_t Start;
  int64_t End;
};

//! Wraps a byte range of the file so that it reads as a complete STEP file:
//! a synthetic "empty header" preamble, then the range, then a synthetic end
//! of the DATA section. The grammar accepts an empty header (rule stepf2), so
//! every chunk parses with the unmodified scanner and parser.
//! The first chunk carries the real header and needs no preamble; the last one
//! runs to the end of the file and so needs no postamble.
class StepFile_ChunkStreamBuf : public std::streambuf
{
public:
  StepFile_ChunkStreamBuf(const char* theName, const StepFile_ChunkRange& theRange, bool theIsFirst,
                          bool theIsLast)
      : myRemaining(theRange.End - theRange.Start),
        myBuffer(THE_BLOCK_SIZE),
        myStage(theIsFirst ? Stage_Body : Stage_Preamble),
        myIsLast(theIsLast)
  {
    const occ::handle<OSD_FileSystem>& aFileSystem = OSD_FileSystem::DefaultFileSystem();
    myStream = aFileSystem->OpenIStream(theName, std::ios::in | std::ios::binary);
    if (myStream != nullptr && !myStream->fail())
    {
      myStream->seekg((std::streamoff)theRange.Start, std::ios::beg);
    }
  }

  //! False when the file could not be opened or positioned.
  bool IsOpen() const { return myStream != nullptr && !myStream->fail(); }

protected:
  int_type underflow() override
  {
    if (gptr() != nullptr && gptr() < egptr())
    {
      return traits_type::to_int_type(*gptr());
    }
    switch (myStage)
    {
      case Stage_Preamble: {
        myStage = Stage_Body;
        setBuffer(THE_CHUNK_PREAMBLE, sizeof(THE_CHUNK_PREAMBLE) - 1);
        return traits_type::to_int_type(*gptr());
      }
      case Stage_Body: {
        if (myRemaining > 0 && IsOpen())
        {
          const std::streamsize aWanted =
            (std::streamsize)std::min<int64_t>(myRemaining, (int64_t)THE_BLOCK_SIZE);
          myStream->read(myBuffer.data(), aWanted);
          const std::streamsize aRead = myStream->gcount();
          if (aRead > 0)
          {
            myRemaining -= (int64_t)aRead;
            setg(myBuffer.data(), myBuffer.data(), myBuffer.data() + aRead);
            return traits_type::to_int_type(*gptr());
          }
        }
        myStage = Stage_Postamble;
        if (!myIsLast)
        {
          setBuffer(THE_CHUNK_POSTAMBLE, sizeof(THE_CHUNK_POSTAMBLE) - 1);
          return traits_type::to_int_type(*gptr());
        }
      }
        Standard_FALLTHROUGH
      case Stage_Postamble:
      default:
        return traits_type::eof();
    }
  }

private:
  void setBuffer(const char* theText, size_t theLength)
  {
    myFixed.assign(theText, theText + theLength);
    setg(myFixed.data(), myFixed.data(), myFixed.data() + myFixed.size());
  }

private:
  enum Stage
  {
    Stage_Preamble,
    Stage_Body,
    Stage_Postamble
  };

  static const size_t THE_BLOCK_SIZE = 1024 * 1024;

  std::shared_ptr<std::istream> myStream;
  int64_t                       myRemaining;
  std::vector<char>             myBuffer;
  std::vector<char>             myFixed;
  Stage                         myStage;
  bool                          myIsLast;
};

//! Finds the byte offsets at which the DATA section can be cut into chunks that
//! each parse on their own: only a ';' seen outside a comment and outside a
//! quoted string terminates a record.
//!
//! Mirrors the scanner's own view of the text (step.lex): comments are
//! "/* ... */", and a quoted string ends at a "'" followed - across blanks and
//! newlines - by ')' or ',', which is what makes doubled quotes and embedded
//! ';' safe.
//!
//! Returns false whenever the file is not safely divisible, in which case the
//! caller parses it in one piece: no DATA section, more than one of them, a
//! "&SCOPE" block (which spans several records), or a blank run so long that
//! the end of a string cannot be decided from the look-ahead window.
class StepFile_Prescanner
{
public:
  StepFile_Prescanner(std::istream& theStream, int64_t theSize)
      : myStream(theStream),
        mySize(theSize),
        myBuffer(THE_WINDOW + THE_LOOKAHEAD),
        myBase(0),
        myPos(0),
        myLength(0),
        myIsEOF(false)
  {
  }

  bool Perform(int theNbChunks, std::vector<StepFile_ChunkRange>& theChunks)
  {
    int64_t aDataStart = -1;
    int     aNbData    = 0;
    // Next offset at which a cut is wanted; the cut lands on the first record
    // boundary at or after it, so chunks are near equal but never mid-record.
    int     aNextChunk = 1;
    int64_t aTarget    = -1;
    std::vector<int64_t> aCuts;

    State aState = State_Normal;
    while (ensure(THE_LOOKAHEAD))
    {
      const char aChar = myBuffer[myPos];
      switch (aState)
      {
        case State_Normal: {
          if (aChar == '/' && peek(1) == '*')
          {
            aState = State_Comment;
            advance(2);
            continue;
          }
          if (aChar == '\'')
          {
            aState = State_Text;
            advance(1);
            continue;
          }
          if (aChar == ';')
          {
            advance(1);
            if (aDataStart >= 0 && offset() >= aTarget && aNextChunk < theNbChunks
                && offset() < mySize)
            {
              aCuts.push_back(offset());
              ++aNextChunk;
              aTarget = mySize * aNextChunk / theNbChunks;
            }
            continue;
          }
          if (isWordStart(aChar))
          {
            const Keyword aKeyword = readKeyword();
            if (aKeyword == Keyword_Scope)
            {
              return false; // a scope spans records: not divisible
            }
            if (aKeyword == Keyword_Data)
            {
              if (++aNbData > 1)
              {
                return false; // several DATA sections: not divisible
              }
              aDataStart = offset();
              aTarget    = mySize / theNbChunks;
            }
            continue;
          }
          advance(1);
          continue;
        }
        case State_Comment: {
          if (aChar == '*' && peek(1) == '/')
          {
            aState = State_Normal;
            advance(2);
            continue;
          }
          advance(1);
          continue;
        }
        case State_Text: {
          if (aChar != '\'')
          {
            advance(1);
            continue;
          }
          // The scanner ends a string only on a quote followed by ')' or ','.
          size_t aLook = 1;
          while (true)
          {
            if (!ensure(aLook + 1))
            {
              return false;
            }
            if (myPos + aLook >= myLength)
            {
              return false;
            }
            const char aNext = myBuffer[myPos + aLook];
            if (aNext == ' ' || aNext == '\n' || aNext == '\r' || aNext == '\t')
            {
              if (++aLook > THE_LOOKAHEAD)
              {
                return false; // cannot decide: leave the file undivided
              }
              continue;
            }
            if (aNext == ')' || aNext == ',')
            {
              aState = State_Normal;
            }
            break;
          }
          advance(1);
          continue;
        }
      }
    }

    if (aDataStart < 0 || aCuts.empty())
    {
      return false;
    }
    theChunks.clear();
    int64_t aStart = 0;
    for (size_t i = 0; i < aCuts.size(); ++i)
    {
      theChunks.push_back(StepFile_ChunkRange{aStart, aCuts[i]});
      aStart = aCuts[i];
    }
    theChunks.push_back(StepFile_ChunkRange{aStart, mySize});
    return true;
  }

private:
  enum State
  {
    State_Normal,
    State_Comment,
    State_Text
  };

  enum Keyword
  {
    Keyword_None,
    Keyword_Data,
    Keyword_Scope
  };

  static bool isWordStart(char theChar)
  {
    return theChar == '&' || (theChar >= 'A' && theChar <= 'Z') || (theChar >= 'a' && theChar <= 'z');
  }

  static bool isWordChar(char theChar)
  {
    return theChar == '&' || theChar == '-' || theChar == '_' || (theChar >= '0' && theChar <= '9')
           || (theChar >= 'A' && theChar <= 'Z') || (theChar >= 'a' && theChar <= 'z');
  }

  //! Consumes one word and reports whether it is a keyword the split must know
  //! about. Consuming it whole also stops "ENDSCOPE" being seen as "SCOPE".
  Keyword readKeyword()
  {
    char   aWord[16];
    size_t aLength = 0;
    while (myPos < myLength && isWordChar(myBuffer[myPos]))
    {
      if (aLength + 1 < sizeof(aWord))
      {
        const char aChar = myBuffer[myPos];
        aWord[aLength++] = (aChar >= 'a' && aChar <= 'z') ? (char)(aChar - 'a' + 'A') : aChar;
      }
      else
      {
        aLength = sizeof(aWord); // too long to be one of ours
      }
      advance(1);
      if (myPos >= myLength && !ensure(1))
      {
        break;
      }
    }
    if (aLength >= sizeof(aWord))
    {
      return Keyword_None;
    }
    aWord[aLength] = '\0';
    if (aLength == 6 && strcmp(aWord, "&SCOPE") == 0)
    {
      return Keyword_Scope;
    }
    if (aLength == 4 && strcmp(aWord, "DATA") == 0)
    {
      return Keyword_Data;
    }
    return Keyword_None;
  }

  int64_t offset() const { return myBase + (int64_t)myPos; }

  char peek(size_t theAhead) const
  {
    return (myPos + theAhead < myLength) ? myBuffer[myPos + theAhead] : '\0';
  }

  void advance(size_t theCount) { myPos += theCount; }

  //! Guarantees theCount readable bytes ahead of the cursor unless the file ends.
  bool ensure(size_t theCount)
  {
    if (myPos + theCount <= myLength)
    {
      return myPos < myLength;
    }
    if (myIsEOF)
    {
      return myPos < myLength;
    }
    const size_t aRest = myLength - myPos;
    memmove(myBuffer.data(), myBuffer.data() + myPos, aRest);
    myBase += (int64_t)myPos;
    myPos   = 0;
    myLength = aRest;
    myStream.read(myBuffer.data() + myLength, (std::streamsize)(myBuffer.size() - myLength));
    const std::streamsize aRead = myStream.gcount();
    if (aRead <= 0)
    {
      myIsEOF = true;
    }
    myLength += (size_t)aRead;
    return myPos < myLength;
  }

private:
  static const size_t THE_WINDOW    = 1024 * 1024;
  static const size_t THE_LOOKAHEAD = 4096;

  std::istream&     myStream;
  int64_t           mySize;
  std::vector<char> myBuffer;
  int64_t           myBase;
  size_t            myPos;
  size_t            myLength;
  bool              myIsEOF;
};
} // namespace

void StepFile_Interrupt(const char* theErrorMessage, const bool theIsFail)
{
  if (theErrorMessage == nullptr)
  {
    return;
  }

  Message_Messenger::StreamBuffer sout = theIsFail ? Message::SendFail() : Message::SendTrace();
  sout << "**** ERR StepFile : " << theErrorMessage << "    ****" << '\n';
}

static int StepFile_Read(const char*                                 theName,
                         std::istream*                               theIStream,
                         const occ::handle<StepData_StepModel>&      theStepModel,
                         const occ::handle<StepData_Protocol>&       theProtocol,
                         const occ::handle<StepData_FileRecognizer>& theRecogHeader,
                         const occ::handle<StepData_FileRecognizer>& theRecogData)
{
  // if stream is not provided, open file stream here
  std::istream*                 aStreamPtr = theIStream;
  std::shared_ptr<std::istream> aFileStream;
  if (aStreamPtr == nullptr)
  {
    const occ::handle<OSD_FileSystem>& aFileSystem = OSD_FileSystem::DefaultFileSystem();
    aFileStream = aFileSystem->OpenIStream(theName, std::ios::in | std::ios::binary);
    aStreamPtr  = aFileStream.get();
  }
  if (aStreamPtr == nullptr || aStreamPtr->fail())
  {
    return -1;
  }

  // Phase timings follow the trace stream this function already reports to,
  // so a slow read can be attributed (scan, records, prepare, load) without a
  // special build.
  OSD_Timer c;
  c.Reset();
  c.Start();

  Message_Messenger::StreamBuffer sout = Message::SendTrace();
  sout << "      ...    Step File Reading : '" << theName << "'";

  // The scanner and the parser only build a flat table of records - every #id
  // cross reference is resolved later, by Prepare() and LoadModel() - so the
  // DATA section can be cut into record-aligned chunks and scanned in
  // parallel, then transcribed in file order for an identical result.
  std::vector<std::unique_ptr<StepFile_ReadData>> aChunkData;
  std::vector<StepFile_ChunkRange>                aChunks;
  if (theIStream == nullptr && theName != nullptr
      && Interface_Static::IVal("read.step.parallel.parse") != 0)
  {
    int aNbChunks = std::min(OSD_Parallel::NbLogicalProcessors(), THE_MAX_CHUNKS);
    aStreamPtr->seekg(0, std::ios::end);
    const int64_t aSize = (int64_t)aStreamPtr->tellg();
    aStreamPtr->seekg(0, std::ios::beg);
    if (aSize < THE_MIN_CHUNKED_SIZE || aNbChunks < 2 || aStreamPtr->fail())
    {
      aNbChunks = 1;
    }
    if (aNbChunks > 1)
    {
      StepFile_Prescanner aPrescanner(*aStreamPtr, aSize);
      if (!aPrescanner.Perform(aNbChunks, aChunks))
      {
        aChunks.clear();
      }
      aStreamPtr->clear();
      aStreamPtr->seekg(0, std::ios::beg);
    }
  }

  if (aChunks.size() > 1)
  {
    aChunkData.resize(aChunks.size());
    std::vector<int> aFailed(aChunks.size(), 0);
    OSD_Parallel::For(
      0,
      (int)aChunks.size(),
      [&](const int theIndex) {
        aChunkData[theIndex].reset(new StepFile_ReadData());
        StepFile_ChunkStreamBuf aChunkBuf(theName,
                                          aChunks[theIndex],
                                          theIndex == 0,
                                          theIndex + 1 == (int)aChunks.size());
        if (!aChunkBuf.IsOpen())
        {
          aFailed[theIndex] = 1;
          return;
        }
        std::istream aChunkStream(&aChunkBuf);
        try
        {
          OCC_CATCH_SIGNALS
          step::scanner aScanner(aChunkData[theIndex].get(), &aChunkStream);
          aScanner.yyrestart(&aChunkStream);
          step::parser aParser(&aScanner);
          if (aParser.parse() != 0)
          {
            aFailed[theIndex] = 1;
          }
        }
        catch (Standard_Failure const&)
        {
          aFailed[theIndex] = 1;
        }
        // The grammar recovers from a malformed record instead of failing, so
        // a recorded error is the only sign that a cut landed badly. Whatever
        // the reason, give the file back to the undivided parser rather than
        // trust a partial result.
        if (aChunkData[theIndex]->GetLastError() != nullptr)
        {
          aFailed[theIndex] = 1;
        }
      },
      false);
    if (std::find(aFailed.begin(), aFailed.end(), 1) != aFailed.end())
    {
      aChunkData.clear();
      aStreamPtr->clear();
      aStreamPtr->seekg(0, std::ios::beg);
    }
  }

  if (aChunkData.empty())
  {
    aChunkData.emplace_back(new StepFile_ReadData());
    try
    {
      OCC_CATCH_SIGNALS
      int           aLetat = 0;
      step::scanner aScanner(aChunkData.front().get(), aStreamPtr);
      aScanner.yyrestart(aStreamPtr);
      step::parser aParser(&aScanner);
      aLetat = aParser.parse();
      if (aLetat != 0)
      {
        StepFile_Interrupt(aChunkData.front()->GetLastError(), true);
        return 1;
      }
    }
    catch (Standard_Failure const& anException)
    {
      Message::SendFail() << " ...  Exception Raised while reading Step File : '" << theName
                          << "':\n"
                          << anException << "    ...";
      return 1;
    }
  }

  c.Show(sout);

  sout << "      ...    STEP File   Read    ...\n";

  std::lock_guard<std::mutex> aLock(GetGlobalReadMutex());

  // Only the first chunk holds the header; the rest were parsed behind a
  // synthetic empty one. Records keep their file order across chunks, and
  // sub-list idents ($1, $2 ...) are numbered per entity rather than per file,
  // so concatenating the chunks reproduces the undivided parse exactly.
  int nbhead = 0, nbrec = 0, nbpar = 0;
  for (size_t aChunk = 0; aChunk < aChunkData.size(); ++aChunk)
  {
    int aNbHead = 0, aNbRec = 0, aNbPar = 0;
    aChunkData[aChunk]->GetFileNbR(&aNbHead, &aNbRec, &aNbPar); // renvoi par lex/yacc
    if (aChunk == 0)
    {
      nbhead = aNbHead;
    }
    nbrec += aNbRec;
    nbpar += aNbPar;
  }
  occ::handle<StepData_StepReaderData> undirec =
    // clang-format off
    new StepData_StepReaderData(nbhead,nbrec,nbpar, theStepModel->SourceCodePage());  // creation tableau de records
  // clang-format on
  int nr = 0;
  for (size_t aChunk = 0; aChunk < aChunkData.size(); ++aChunk)
  {
    StepFile_ReadData& aFileDataModel = *aChunkData[aChunk];
    int                aNbHead = 0, aNbRec = 0, aNbPar = 0;
    aFileDataModel.GetFileNbR(&aNbHead, &aNbRec, &aNbPar);
    for (int aRec = 1; aRec <= aNbRec; aRec++)
    {
      int   nbarg;
      char* ident;
      char* typrec = nullptr;
      aFileDataModel.GetRecordDescription(&ident, &typrec, &nbarg);
      undirec->SetRecord(++nr, ident, typrec, nbarg);

      if (nbarg > 0)
      {
        Interface_ParamType typa;
        char*               val;
        while (aFileDataModel.GetArgDescription(&typa, &val) == 1)
        {
          undirec->AddStepParam(nr, val, typa);
        }
      }
      undirec->InitParams(nr);
      aFileDataModel.NextRecord();
    }
    aFileDataModel.ErrorHandle(undirec->GlobalCheck());
  }

  int anFailsCount = undirec->GlobalCheck()->NbFails();
  if (anFailsCount > 0)
  {
    Message::SendInfo() << "**** ERR StepFile : Incorrect Syntax : Fails Count : " << anFailsCount
                        << " ****";
  }

  for (size_t aChunk = 0; aChunk < aChunkData.size(); ++aChunk)
  {
    aChunkData[aChunk]->ClearRecorder(1);
  }

  sout << "      ... Step File loaded  ...\n";
  sout << "   " << undirec->NbRecords() << " records (entities,sub-lists,scopes), " << nbpar
       << " parameters";

  c.Show(sout);

  //   Analyse : par StepReaderTool

  StepData_StepReaderTool readtool(undirec, theProtocol);
  readtool.SetErrorHandle(true);

  readtool.PrepareHeader(theRecogHeader); // Header. reco nul -> pour Protocol
  readtool.Prepare(theRecogData);         // Data.   reco nul -> pour Protocol

  sout << "      ... Parameters prepared ...\n";

  c.Show(sout);

  readtool.LoadModel(theStepModel);
  if (theStepModel->Protocol().IsNull())
  {
    theStepModel->SetProtocol(theProtocol);
  }
  // Parameter values point into each chunk's own text storage, so it may only
  // be released once the model has been loaded from them.
  for (size_t aChunk = 0; aChunk < aChunkData.size(); ++aChunk)
  {
    aChunkData[aChunk]->ClearRecorder(2);
  }
  anFailsCount = undirec->GlobalCheck()->NbFails() - anFailsCount;
  if (anFailsCount > 0)
  {
    Message::SendInfo() << "*** ERR StepReaderData : Unresolved Reference : Fails Count : "
                        << anFailsCount << " ***";
  }

  readtool.Clear();
  undirec.Nullify();

  sout << "      ...   Objects analysed  ...\n";
  int n = theStepModel->NbEntities();
  sout << "  STEP Loading done : " << n << " Entities";

  c.Show(sout);

  return 0;
}

int StepFile_Read(const char*                            theName,
                  std::istream*                          theIStream,
                  const occ::handle<StepData_StepModel>& theStepModel,
                  const occ::handle<StepData_Protocol>&  theProtocol)
{
  occ::handle<StepData_FileRecognizer> aNulRecog;
  return StepFile_Read(theName, theIStream, theStepModel, theProtocol, aNulRecog, aNulRecog);
}
