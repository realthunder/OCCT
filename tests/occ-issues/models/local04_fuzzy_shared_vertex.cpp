// local04: non-destructive fuzzy Booleans on local04_fuzzy_shared_vertex.brep,
// reporting every INPUT vertex, edge or face whose tolerance an operation
// changed. See ../README.md. FreeCAD's Part booleans copy their tools when a
// fuzzy value is given, which hides this, so it is driven through OCCT here.
//
// Build (Windows, OCCT install prefix P):
//   cl /EHsc /MD /O2 /std:c++17 /IP\inc local04_fuzzy_shared_vertex.cpp
//      /link /LIBPATH:P\win64\vc14\libi TK*.lib   (list the libs: link does not glob)
// Run: local04_fuzzy_shared_vertex local04_fuzzy_shared_vertex.brep 2.512e-7 3.162e-7
//   expected while open: 3.162e-7 changes V (1.94e-6 -> 1e-5) in all four ops,
//   2.512e-7 changes nothing.
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <cstdio>
#include <cstdlib>
#include <vector>

static double tolOf(const TopoDS_Shape& s)
{
    switch (s.ShapeType()) {
    case TopAbs_VERTEX: return BRep_Tool::Tolerance(TopoDS::Vertex(s));
    case TopAbs_EDGE: return BRep_Tool::Tolerance(TopoDS::Edge(s));
    case TopAbs_FACE: return BRep_Tool::Tolerance(TopoDS::Face(s));
    default: return 0.;
    }
}

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::fprintf(stderr, "usage: local04_fuzzy_shared_vertex file.brep fuzzy...\n");
        return 2;
    }
    for (int k = 2; k < argc; ++k) {
        double fuzzy = std::atof(argv[k]);
        const char* ops[] = {"fuse", "cut", "common", "cut21"};
        for (int op = 0; op < 4; ++op) {
            TopoDS_Shape comp;
            BRep_Builder bb;
            if (!BRepTools::Read(comp, argv[1], bb)) {
                std::fprintf(stderr, "cannot read %s\n", argv[1]);
                return 2;
            }
            std::vector<TopoDS_Shape> args;
            for (TopoDS_Iterator it(comp); it.More(); it.Next())
                args.push_back(it.Value());
            TopTools_IndexedMapOfShape all;
            TopExp::MapShapes(comp, TopAbs_VERTEX, all);
            TopExp::MapShapes(comp, TopAbs_EDGE, all);
            TopExp::MapShapes(comp, TopAbs_FACE, all);
            std::vector<double> before(all.Extent());
            for (int i = 1; i <= all.Extent(); ++i)
                before[i - 1] = tolOf(all(i));
            TopTools_ListOfShape la, lt;
            const TopoDS_Shape& a = op == 3 ? args[1] : args[0];
            const TopoDS_Shape& t = op == 3 ? args[0] : args[1];
            la.Append(a);
            lt.Append(t);
            BRepAlgoAPI_BooleanOperation* mk = nullptr;
            if (op == 0) mk = new BRepAlgoAPI_Fuse;
            else if (op == 2) mk = new BRepAlgoAPI_Common;
            else mk = new BRepAlgoAPI_Cut;
            mk->SetArguments(la);
            mk->SetTools(lt);
            mk->SetNonDestructive(Standard_True);
            if (fuzzy > 0.)
                mk->SetFuzzyValue(fuzzy);
            std::fprintf(stderr, "CASE %s fuzzy=%g op=%s\n", argv[1], fuzzy, ops[op]);
            mk->Build();
            int changed = 0;
            for (int i = 1; i <= all.Extent(); ++i) {
                double now = tolOf(all(i));
                if (now != before[i - 1]) {
                    if (changed < 4) {
                        const char* tn = all(i).ShapeType() == TopAbs_VERTEX ? "V"
                                         : all(i).ShapeType() == TopAbs_EDGE ? "E" : "F";
                        std::printf("  input %s%d %.3g -> %.3g\n", tn, i, before[i - 1], now);
                    }
                    ++changed;
                }
            }
            double vol = -1.;
            bool valid = false;
            if (mk->IsDone()) {
                GProp_GProps g;
                BRepGProp::VolumeProperties(mk->Shape(), g);
                vol = g.Mass();
                valid = BRepCheck_Analyzer(mk->Shape()).IsValid();
            }
            std::printf("ROW %s fuzzy=%g op=%s done=%d vol=%.6f valid=%d changed=%d\n", argv[1], fuzzy,
                        ops[op], (int)mk->IsDone(), vol, (int)valid, changed);
            std::fflush(stdout);
            delete mk;
        }
    }
    return 0;
}
