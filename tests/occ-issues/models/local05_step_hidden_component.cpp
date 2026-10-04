// local05: STEPCAFControl_Writer dereferences a null override item when an
// assembly component is INVISIBLE and carries no colour of its own. See
// ../README.md. The writer styles a component from its own label only; with
// no colour there MakeSTEPStyles takes its "default white" branch and calls
// setDefaultInstanceColor(theOverride, ...), and theOverride is null whenever
// getStyledItem() found nothing -- always, for the components of an assembly
// written in one pass, since writeColors() reaches them before the MDGPR of
// the shape they refer to exists.
//
// Build (Windows, OCCT install prefix P):
//   cl /EHa /MD /O2 /std:c++17 /IP\inc local05_step_hidden_component.cpp
//      /link /LIBPATH:P\win64\vc14\libi TK*.lib   (list the libs: link does not glob)
// Run: local05_step_hidden_component <out dir>
//   expected while open: "none" and "ref" CRASH (access violation), "own" OK.
//   FreeCAD hit "ref": every GUI export colours the shape, not the instance.
#include <BRepPrimAPI_MakeBox.hxx>
#include <OSD.hxx>
#include <Quantity_ColorRGBA.hxx>
#include <STEPCAFControl_Writer.hxx>
#include <Standard_ErrorHandler.hxx>
#include <Standard_Failure.hxx>
#include <TDF_Label.hxx>
#include <TDocStd_Document.hxx>
#include <XCAFApp_Application.hxx>
#include <XCAFDoc_ColorTool.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <gp_Trsf.hxx>
#include <cstdio>
#include <string>
#include <typeinfo>

// colour: "none" -- no colour anywhere; "ref" -- on the referred shape only;
// "own" -- on the component itself (FreeCAD's workaround)
static const char* writeOne(const std::string& dir, const char* colour)
{
    Handle(XCAFApp_Application) app = XCAFApp_Application::GetApplication();
    Handle(TDocStd_Document) doc;
    app->NewDocument("MDTV-XCAF", doc);
    Handle(XCAFDoc_ShapeTool) shapes = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
    Handle(XCAFDoc_ColorTool) colours = XCAFDoc_DocumentTool::ColorTool(doc->Main());

    // Two boxes, as FreeCAD writes an App::Part of Box1 + Box2: each box a
    // shape of its own, each placed by a component of the assembly.
    TDF_Label box1 = shapes->AddShape(BRepPrimAPI_MakeBox(10, 10, 10).Shape(), Standard_False);
    TDF_Label box2 = shapes->AddShape(BRepPrimAPI_MakeBox(10, 10, 10).Shape(), Standard_False);
    TDF_Label assembly = shapes->NewShape();
    gp_Trsf move;
    shapes->AddComponent(assembly, box1, TopLoc_Location());
    move.SetTranslation(gp_Vec(20, 0, 0));
    TDF_Label hidden = shapes->AddComponent(assembly, box2, TopLoc_Location(move));
    shapes->UpdateAssemblies();

    Quantity_ColorRGBA red(1, 0, 0, 1);
    if (std::string(colour) == "ref") {
        colours->SetColor(box1, red, XCAFDoc_ColorSurf);
        colours->SetColor(box2, red, XCAFDoc_ColorSurf);
    }
    else if (std::string(colour) == "own") {
        colours->SetColor(hidden, red, XCAFDoc_ColorSurf);
    }
    colours->SetVisibility(hidden, Standard_False);

    std::string path = dir + "/local05_" + colour + ".step";
    const char* result = "OK";
    try {
        OCC_CATCH_SIGNALS
        STEPCAFControl_Writer writer;
        writer.SetColorMode(Standard_True);
        if (!writer.Transfer(doc, STEPControl_AsIs)) {
            result = "TRANSFER FAILED";
        }
        else if (writer.Write(path.c_str()) != IFSelect_RetDone) {
            result = "WRITE FAILED";
        }
    }
    catch (const Standard_Failure& e) {
        static std::string msg;
        msg = std::string("CRASH ") + typeid(e).name() + ": " + e.what();
        result = msg.c_str();
    }
    app->Close(doc);
    return result;
}

int main(int argc, char** argv)
{
    OSD::SetSignal(Standard_False);
    std::string dir = argc > 1 ? argv[1] : ".";
    for (const char* colour : {"own", "ref", "none"}) {
        std::printf("%-5s %s\n", colour, writeOne(dir, colour));
        std::fflush(stdout);
    }
    return 0;
}
