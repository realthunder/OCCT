// Created on: 1993-01-20
// Created by: Remi LEQUETTE
// Copyright (c) 1993-1999 Matra Datavision
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

#include <TopTools.hxx>

#include <TopTools_ShapeSet.hxx>

//=======================================================================
//function : Dump
//purpose  : 
//=======================================================================
void  TopTools::Dump(const TopoDS_Shape& Sh, Standard_OStream& S)
{
  TopTools_ShapeSet SSet;
  SSet.Add(Sh);
  SSet.Dump(Sh,S);
  SSet.Dump(S);
}


void TopTools::Dummy(const Standard_Integer)
{
}

static FuncShowTopoShape *_FuncShowTopoShape;

static thread_local char _ShapeName[256];
#define OCCT_EXT_VERSION 1
Standard_Integer SetFuncShowTopoShape(FuncShowTopoShape *func)
{
  _FuncShowTopoShape = func;
  return OCCT_EXT_VERSION;
}

static inline Standard_Boolean _ShowTopoShape (const char *Key, int line, const TopoDS_Shape& S, const char *Name, Standard_Boolean Oriented)
{
  if (!_FuncShowTopoShape)
    return Standard_False;

  if (Oriented) {
    const char *Postfix = nullptr;
    if (S.Orientation() == TopAbs_FORWARD)
      Postfix = "_F";
    else if (S.Orientation() == TopAbs_REVERSED)
      Postfix = "_R";
    else if (S.Orientation() == TopAbs_INTERNAL)
      Postfix = "_I";
    else if (S.Orientation() == TopAbs_EXTERNAL)
      Postfix = "_X";
    
    if (Postfix) {
      if (Name != _ShapeName) {
        snprintf(_ShapeName, sizeof(_ShapeName), "%s%s", Name, Postfix);
        Name = _ShapeName;
      }
      else {
        size_t len = strlen(Name);
        if (len+1 < sizeof(_ShapeName) - strlen(Postfix))
          sprintf(_ShapeName+len, "%s", Postfix);
      }
    }
  }

  return _FuncShowTopoShape(Key, line, S, Name);
}

void ShowTopoShape (const char *Key, int line, const TopoDS_Shape& S, const char *Name, Standard_Boolean Oriented)
{
  _ShowTopoShape(Key, line, S, Name, Oriented);
}


void ShowTopoShape (const char *Key, int line, const TopoDS_Shape& S, const char *Name, const TopoDS_Shape &S2, Standard_Boolean Oriented)
{
  if (!_FuncShowTopoShape)
    return;
  snprintf(_ShapeName, sizeof(_ShapeName), "%s1", Name);
  if (!_ShowTopoShape(Key, line, S, _ShapeName, Oriented))
    return;
  snprintf(_ShapeName, sizeof(_ShapeName), "%s2", Name);
  _ShowTopoShape(nullptr, line, S2, _ShapeName, Oriented);
}

void ShowTopoShape (const char *Key, int line, const TopoDS_Shape& S, const char *Name, const TopTools_ListOfShape &Shapes, Standard_Boolean Oriented)
{
  if (!_FuncShowTopoShape)
    return;
  if (!_ShowTopoShape(Key, line, S, Name, Oriented))
    return;
  TopTools_ListIteratorOfListOfShape it(Shapes);
  for (; it.More(); it.Next()) {
    snprintf(_ShapeName, sizeof(_ShapeName), "%s_item", Name);
    _ShowTopoShape(nullptr, 0, it.Value(), _ShapeName, Oriented);
  }
}

void ShowTopoShape (const char *Key, int line, const TopoDS_Shape& S, const char *Name, const TopTools_SequenceOfShape &Shapes, Standard_Boolean Oriented)
{
  if (!_FuncShowTopoShape)
    return;
  if (!_ShowTopoShape(Key, line, S, Name, Oriented))
    return;
  for (Standard_Integer ii = 1; ii <= Shapes.Length(); ++ii) {
    snprintf(_ShapeName, sizeof(_ShapeName), "%s_item", Name);
    _ShowTopoShape(nullptr, 0, Shapes(ii), _ShapeName, Oriented);
  }
}

void ShowTopoShape (const char *Key, int line, const TopoDS_Shape& S, const char *Name, const TopTools_IndexedMapOfShape &Shapes, Standard_Boolean Oriented)
{
  if (!_FuncShowTopoShape)
    return;
  if (!_ShowTopoShape(Key, line, S, Name, Oriented))
    return;
  for (Standard_Integer ii = 1; ii <= Shapes.Extent(); ++ii) {
    snprintf(_ShapeName, sizeof(_ShapeName), "%s_item", Name);
    _ShowTopoShape(nullptr, 0, Shapes(ii), _ShapeName, Oriented);
  }
}

void ShowTopoShape (const char *Key, int line, const TopoDS_Shape& S, const char *Name, const TopTools_MapOfShape &Shapes, Standard_Boolean Oriented)
{
  if (!_FuncShowTopoShape)
    return;
  if (!_ShowTopoShape(Key, line, S, Name, Oriented))
    return;
  TopTools_MapIteratorOfMapOfShape it(Shapes);
  for (; it.More(); it.Next()) {
    snprintf(_ShapeName, sizeof(_ShapeName), "%s_item", Name);
    _ShowTopoShape(nullptr, 0, it.Value(), _ShapeName, Oriented);
  }
}
