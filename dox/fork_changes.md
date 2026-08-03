# What this fork changes

This is the **realthunder fork of Open CASCADE Technology**, maintained as the geometry
kernel of the [realthunder fork of FreeCAD](https://github.com/realthunder/FreeCAD). It
carries a small number of changes on top of a stock OCCT release: modelling fixes that
FreeCAD users hit in practice, a STEP reader that can be driven progressively and in
parallel, and the packaging used to publish the result.

Everything here is a delta over upstream. This document lists that delta; it does not
document OCCT itself.

## Branches

| Branch | Base | Role |
| --- | --- | --- |
| `LinkVibe-801` | `V8.0.1` | **Primary.** What FreeCAD builds against, and where new work lands. |
| `LinkVibe` | `V7_7_2` | Kept for the version-guarded fallback paths, and for what the released packages still link. |

`LinkVibe-801` is not a merge of `LinkVibe`: the 8.0 tree was restructured, so the older
branch's fixes were ported (`1789444318`) or re-derived rather than merged. Both branches
are maintained; a fix that matters to released packages may need doing twice.

To regenerate the raw lists behind this document:

```
git log --oneline V8.0.1..LinkVibe-801     # 30 commits
git log --oneline V7_7_2..LinkVibe         # 28 commits
git diff --stat V8.0.1..LinkVibe-801 -- src/
```

## The rule that governs all of it: stay ABI-compatible with upstream

A binary compiled against upstream OCCT headers must keep working against this library.
That has been broken four times and repaired each time; treat it as a hard constraint.

What it forbids in practice:

- **no new data members** on a public class — they change its size;
- **no new virtual methods** — they move every later vtable slot;
- **no appending parameters to an existing method** — that renames the symbol.

What to do instead — every one of these is used somewhere in this fork:

- put the state in a **`thread_local` in the `.cxx`**, when the thread is genuinely its
  owner (`STEPControl_ActorRead`'s deferral state, `XSAlgo_ShapeProcessor`'s messenger,
  `RWGltf_CafWriter`'s node-to-mesh map, the transfer reader's encoded-shape map);
- add a **new overload** and keep the old signature forwarding to it
  (`STEPCAFControl_Reader::Transfer`);
- add **non-virtual** methods, or **static** ones — neither affects layout or mangling;
- move a feature **down** to the concrete class that actually implements it, rather than
  adding a virtual to a generic base (`TransferRootsDeferred` moved from
  `XSControl_Reader` to `STEPControl_Reader`).

The check that settles it: compile a `sizeof` probe against upstream headers and against
these, and compare. The repairs — `db4c5ac51c`, `4eeaad075c` — record the expected values
(`296/136/136`, `64`, `328`).

## 1. STEP import: progressive and parallel

The largest workstream, and the reason most of the new public API exists. A stock STEP
import is one indivisible serial operation: nothing can be shown until the whole file is
through, and it uses one core throughout. This fork breaks it into pieces that can be
handed over as they are ready, and spreads the expensive parts across cores.

FreeCAD drives this from `src/Mod/Import/App/ReaderStep.cpp`.

**Reading the file** — `7d20eb5521` cuts the DATA section into record-aligned chunks and
scans them in parallel. The scanner and parser are used unmodified: they only build a flat
record table, `#id` references being resolved afterwards, so a serial pre-pass finding safe
cut points (respecting comments and quoted strings) is enough. Records keep file order and
sub-list idents are numbered per entity, so the chunked parse reproduces the undivided one
exactly. Anything unusual — no DATA section, several of them, a `&SCOPE` block spanning
records, a chunk reporting a syntax error — falls back to the undivided parser.
*FGC-9 (70 MB, 780935 entities), debug: `ReadFile` 13.2 s → 8.8 s; release 2.74 s → 1.76 s.*

**Translating in pieces** — three commits build up the streaming API:

- `c513c6a702` gives `TransferRootsDeferred()` a `[first, last]` range and exposes it as
  `STEPCAFControl_Reader::TransferRootRange()`.
- `a20a712ccb` adds `RootComponents()` / `TransferComponents()`, so the single root that a
  typical STEP file has can be handed over component by component. Two costs had to be
  paid off to make finer granularity worthwhile: only the styles and their usage overrides
  are now read per batch (the rest, names included, come with the root), and each deferred
  flush now starts rewriting binders where the previous one stopped.
- `44b20b893b` adds `RootAssemblyTree()`, walking the product structure all the way down —
  placements and product names read *without translating anything*, so an importer can
  create the containers of an assembly before any of its geometry exists. Occurrences whose
  interpretation needs the whole child list (a product used twice, a component sharing its
  representation with a sibling, an assembly resolving to fewer than two representations)
  are reported apart, because a reducing importer must decide about those itself.

**Healing in parallel** — shape healing is roughly half of a STEP import and ran inline per
solid, deep inside the serial transfer recursion.

- `b1060d3eb5` defers each manifold solid's healing request and heals the batch through
  `OSD_Parallel`, then merges the modifications in deterministic record order and applies
  them in one location-aware `BRepTools_ReShape` pass over every binder. Parent compounds
  are assembled from raw shapes before healing runs, so whole subtrees — not just edges —
  must be rewritten. One `Message_Report`-backed messenger per task, since the default
  messenger is not thread-safe.
  *Debug, 24 threads: MAN200 80 s → 41 s, FGC-9 196 s → 77 s, 92_razbor 196 s → 98 s.*
- `941d0c88e9` replaces the end-of-batch barrier with a worker pool fed as translation
  produces shapes, so translating and healing overlap and a long shape starts healing as
  soon as it exists. This is what stops a one-component batch from costing a fixed
  serial minimum. *FGC-9 streamed 130.5 s → 103.0 s, first object 12.3 s → 9.8 s.*

**Encoding regularity once, and off the reading thread** — every result a reader hands back
goes through `ShapeFix::EncodeRegularity`, which walks each face and edge of it.

- `e90e6283ec`: a progressive import asks for its results piece by piece *and then for the
  whole*, so a file was encoded twice over. `EncodeRegularity` already keeps a map of what
  it has processed; it is now handed to the caller, and the transfer reader keeps one per
  model. *FGC-9 root transfer 22.1 s → 0.65 s; streamed transfer 107 s → 80 s, below the
  82 s of the same file read in one go.*
- `cfa9bb392f`: the heal-ahead worker that healed a shape encodes it too, and the flush
  tells the reader what has been done. Workers are kept off one another's shapes by
  claiming the parts first — encoding writes continuity onto edges that an assembly and its
  parts, or one component used twice, share. *FGC-9 streamed transfer 66.6 s → 53.7 s; the
  reading thread's own encoding 17.5 s → 0 s.*

**Crash and correctness fixes found along the way** — `18fd4542b0`, a style naming a
component of a later batch produced an empty label that `GetShape()` threw on, killing a
streamed import of any coloured assembly on its first batch; `77a6bfb46a`, an
`oriented_closed_shell` whose element is not a closed shell segfaults, because
`RWStepShape_RWOrientedClosedShell::ReadStep` discards the status of its `ReadEntity` call
(the comment saying so has been there since 1999) and the redefined accessors are left
holding a null element.

**Instrumentation** — `e3b3681932` and `e8262e31a7` time the phases of a transfer, the
healing flush profile (how many shapes waited for, summed and longest — a batch is only as
fast as its slowest shape) and each healing operator, reporting on the trace stream. Both
guard message *preparation* behind the new `Message::IsAccepted()`, since counting a
shape's solids, faces and edges walks it three times. Nothing is printed unless a printer
accepts `Message_Trace`, so this costs nothing when nobody is listening.

### Runtime switches

All are `Interface_Static` enums in the `step` group, settable before a read. Each keeps
the previous behaviour one setting away, which is what makes them usable for bisecting a
suspect result.

| Parameter | Values | Default | Effect |
| --- | --- | --- | --- |
| `read.step.parallel.parse` | `Off`, `On` | `On` | Chunked parallel scan of the DATA section. |
| `read.step.parallel.healing` | `Off`, `Serial`, `On`, `Pipeline` | `Pipeline` | `Off` = classic inline per shape; `Serial` = deferred to one batch at the end; `On` = that batch in parallel; `Pipeline` = healed by workers as translation hands each shape over. |
| `read.step.parallel.encoding` | `Off`, `On` | `On` | Regularity encoded by the healing worker instead of by whoever asks for the result. Only bites when healing pipelines. |

### New public API

```cpp
// STEPCAFControl_Reader — progressive transfer
struct AssemblyNode { handle<Standard_Transient> Component; int Parent; gp_Trsf Location;
                      bool IsAssembly; TCollection_AsciiString Name; };
int          RootComponents(num, theUnique, theShared);        // shared ones reported apart
int          RootAssemblyTree(num, doc, theNodes);             // NCollection_Sequence<AssemblyNode>
bool         TransferComponents(theEntities, doc, progress = {});
bool         TransferRootRange(num, theLastNum, doc, progress = {});
TopoDS_Shape ComponentShape(theComponent) const;
bool         Transfer(doc, theNum, theLastNum, ...);   // new overload; old signature forwards

// STEPControl_Reader — deferred transfer (moved here from XSControl_Reader for ABI)
int  TransferRootsDeferred(progress = {}, theFirst = 1, theLast = 0);
int  TransferListDeferred(theList, progress = {});

// STEPControl_ActorRead — non-virtual, no new state
void SetDeferredProcessing(bool);
void FlushDeferredProcessing(tp, progress);

// Supporting
static void XSAlgo_ShapeProcessor::SetContextMessenger(messenger);   // thread-local
static void XSControl_TransferReader::NoteEncodedRegularity(model, shapes);
static void BRepLib::EncodeRegularity(shape, tolAng, processed);     // caller-owned map
static void ShapeFix::EncodeRegularity(shape, tolAng, processed);    //   "
static bool Message::IsAccepted(gravity, messenger);
```

## 2. Offset, thickness and the `BRepAlgo_Loop` chain

`BRepOffset_MakeOffset::MakeThickSolid()` on concave faces is the fork's oldest sore point,
tracked as realthunder/OCCT issues #1–#4. The 7.7.2 branch accumulated a fix chain over
several commits; `1789444318` ports the net result onto the restructured 8.0 tree, adapted
to the 8.0 API (`std::hash` instead of `HashCode()`, the new NCollection hasher concept,
`occ::handle` spelling).

The port then needed three follow-ups, all caused by 8.0's switch from `HashCode()` to
`std::hash` re-rolling the iteration order of every shape-keyed container:

- `933803af41` — first-wins decisions in `BRepAlgo_Loop` had come to depend on heap
  addresses: which near-coincident vertex survives canonicalization, which orientation of a
  wire reaches an orientation-blind dedup first, the order wires reach
  `BRepAlgo_FaceRestrictor`. The decision-feeding containers now iterate in insertion
  order. The same audit fixed six latent defects, including `BRepAlgo_AsDes::Replace(S, S)`
  erasing its own records, `Bubble` reading parameters before computing them, and a
  `static int` debug leftover that dereferenced an unchecked pcurve on the 11th call
  process-wide.
- `756262b692` — on periodic surfaces two wires touching at a point can each be classified
  as inside the other, so `MakeThickSolid` failed about one run in four. Mutual-containment
  pairs are resolved by sampled UV area.
- `72aae3644f` — a degenerated pole edge came out of the vertex-topology DFS as a bogus
  standalone one-edge wire, making both offset sphere faces unorientable. Degenerated edges
  are now diverted from the loop search and added back into the wire crossing their vertex.
  The #4 result is a valid single-shell solid — better than the 7.7.2 branch, which
  returned an invalid shape.

Earlier fixes in the same area: `dfce6bbb3b` (edge splitting in `BRepAlgo_Loop::Build()` —
`CutEdge()` relied on vertex orientation to choose which segment to keep; all segments are
kept and the loop-finding logic decides), `3afddec1dd` and `420b93158d`
(`Geom_RectangularTrimmedSurface` trimming only one parametric dimension, and `Copy()`
returning nullptr once that is possible), `ee10c96cec` (crash on negative offset of a
radius-one circle).

## 3. Other modelling fixes

| Commit | Fix |
| --- | --- |
| `d0ad46e75a` | `BRepFeat_MakePrism` empty result on perform-until-face (FreeCAD #953). |
| `e8e2758d6f` | `BRepFeat_MakePrism` empty result on `JustFeat` (FreeCAD #985). |
| `f6e6ff7879` | Thrusection, by reverting upstream `4607bd0747f`. |
| `b096b91e74` | Accumulated fixes by blobfish across `BRepFill`, `ChFi3d`, `ChFiDS`, `MAT`, `Extrema`. |

## 4. glTF export

`dc4f7a6547` — `RWGltf_CafWriter` emitted one glTF mesh per scene node, so several
instances of the same shape duplicated the mesh in the JSON and instancing was lost on
re-import, even though the accessors were shared. `writeMeshes()` now skips duplicates and
records each node's mesh index in a thread-local map that `writeNodes()` consults. Draco
compression keeps the old per-node behaviour, its buffer bookkeeping assuming it.

## 5. Debugging aid

`TopTools::SetFuncShowTopoShape` / `ShowTopoShape` (in the `1789444318` port, originally
`0e728753b5`/`68724cb966` on the 7.7.2 branch) let a debug build hand intermediate shapes
to the host application. FreeCAD resolves it by `dlsym` in `AppPartPy.cpp` and exposes it
as `Part.showShapeOCCT`, which turns the internals of a failing recompute into document
objects. It is declared `extern "C"` and exported unmangled from `libTKBRep`.

## 6. Packaging

`b41443f967` carries the conda and snap recipes and the GitHub workflows over from the
7.7.2 branch, with the conda version bumped to 8.0.1: `conda/` (recipe, build scripts,
build config), `snap/snapcraft.yaml`, and `.github/workflows/` for the builder, conda and
snap jobs. These feed the distribution images built from `~/works/sw/*-feedstock`.

## 7. Tests and issue corpus

Neither is wired into `ctest` — both are driven by FreeCAD, since the failures they capture
are FreeCAD recomputes.

- **`tests/thickness/`** (`578cec5f17`) — regression suite for the thickness fix chain: the
  four reported issue models as document recomputes with reference volumes, plus
  programmatic plain-solid cases (cylinder, cylinder with hole, hollowed one face at a time
  in both directions) because this chain has a history of breaking *ordinary* thickness.
  Run with `FreeCADCmd tests/thickness/run_tests.py`; exit 0 means no expected-pass case
  regressed. Eight configurations are currently `XFAIL` and documented with audit leads;
  one that starts passing is reported as `UNEXPECTED-PASS` for promotion. The suite has
  stood at PASS=10 XFAIL=8 across the whole STEP workstream.
- **`tests/occ-issues/`** (`ccaa86c63f`) — all 26 issues labelled `occ` on the FreeCAD
  tracker with every recoverable attachment (30 models), each summarized with its analysis
  and a crash-isolated recompute scan: 17 reproduce on plain recompute, 3 of those crashing
  outright, 8 need a scripted step, 1 has no repro. Raw material for growing the suite.

## Appendix: commit list

### `LinkVibe-801` (over `V8.0.1`), oldest first

```
d0ad46e75a  Fix empty output on BRepFeat_MakePrism perform until face
e8e2758d6f  Fix BRepFeat_MakePrism empty result on JustFeat
ee10c96cec  Fix crash on negative offset a circle with radius one
3afddec1dd  Fix Geom_RectangularTrimmedSurface that only trim one parametric dimension
b096b91e74  Accumulated fix by blobfish
dfce6bbb3b  Fixes edge splitting in BRepAlgo_Loop::Build()
420b93158d  Fix Geom_RectangularTrimmedSurface::Copy() returning nullptr
f6e6ff7879  Fix thrusection by reverting 4607bd0747f
dc4f7a6547  RWGltf_CafWriter: share one glTF mesh among instance nodes
b41443f967  Add conda/snap packaging and github workflows
1789444318  Port MakeThickSolid fix chain to OCCT 8.0.1
756262b692  Fix nondeterministic face loss in BRepAlgo_FaceRestrictor
933803af41  Fix hash-order nondeterminism and latent bugs in the ThickSolid chain
72aae3644f  Keep degenerated pole edges in the wire that crosses the singularity
578cec5f17  Add MakeThickSolid regression test suite
ccaa86c63f  Add occ-labeled issue corpus from the FreeCAD tracker
b1060d3eb5  Heal STEP shapes in parallel after the roots translate
c513c6a702  Let STEP readers transfer a contiguous batch of roots
a20a712ccb  Let STEP readers transfer the components of one root
44b20b893b  Let STEP readers list the assembly tree of a root
7d20eb5521  Read STEP files in parallel chunks
941d0c88e9  Heal STEP shapes while the reader goes on translating
db4c5ac51c  Keep the deferred-transfer work off the public ABI
4eeaad075c  Keep the shape processor and the CAF reader off the public ABI too
e3b3681932  Report what the phases of a STEP transfer cost
e90e6283ec  Encode the regularity of a streamed shape once
e8262e31a7  Report what a transfer's healing spends, when anyone is listening
18fd4542b0  STEP: let a style name a component that has not arrived yet
77a6bfb46a  STEP: survive an oriented_closed_shell with no element
cfa9bb392f  Encode a streamed shape's regularity where it is healed
```

### `LinkVibe` (over `V7_7_2`), oldest first

The modelling fixes here are the ancestors of what `1789444318` ported; the packaging
commits are what `b41443f967` squashed.

```
99f4565de2  Accumulated fix by blobfish
70f20cce3c  Fix empty output on BRepFeat_MakePrism perform until face
efdf87cddf  Fix crash on negative offset a circle with radius one
c39b7fbe1d  Fix Geom_RectangularTrimmedSurface that only trim one parametric dimension
1314a2fca3  Fixes edge splitting in BRepAlgo_Loop::Build()
d5fed06b5e  Add github workflow
39b09f2605  Fix github workflow
dee3be32fb  Fix Geom_RectangularTrimmedSurface::Copy() returning nullptr
1cb7c037af  Fix BRepFeat_MakePrism empty result on JustFeat
216cf40c1c  Change github workflow
dddf528dc9  Fix conda build recipe
88f32e4d37  Fix conda build number
4600acf0f2  Try workaround conda error on aarch64
f4846e163b  Fix conda build number with git describe
6b3929872c  Fix MakeThickSolid
68724cb966  Report extension version in SetFuncShowTopoShape
a284bd6c44  Fix thrusection by reverting 4607bd0747f
0e728753b5  Change shape debugging API
e89c971d2e  Fixes handling of periodic curve in MakeThickSolid
d82f5cfde5  Fix missing export
d0a3a771ab  Fix github workflow
5892247209  Handle null curve
ba9e8cadff  Fixing conda build
3a915fd6fa  Update snap build
2d22d37f96  Fixed hanlding of seam edges in MakeThickSolid
bdae59e515  Fix build with newer FreeType
d997dd3465  RWGltf_CafWriter: share one glTF mesh among instance nodes
002e0f0b8e  Fix nondeterministic face loss in BRepAlgo_FaceRestrictor
```
