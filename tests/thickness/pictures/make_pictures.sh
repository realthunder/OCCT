#!/bin/bash
# make_pictures.sh -- regenerate models/pictures (README.md, "Pictures").
#
# Builds scratch TKBool/TKOffset libraries for upstream's chain files and for
# the fork just before each fix (mkold.sh), runs every case of cases.py on
# each and on the fork as built (compute.py, the libraries preloaded), renders
# the panels under Xvfb (mkjobs.py, render.py) and lays them out (compose.py).
#
# env: THICK_WORK  scratch directory (default /tmp/thickness-pictures; the
#                  libraries in it are reused when present)
#      FCBUILD     a FreeCAD build linked against this OCCT tree
#      RUN         the wrapper that runs a command in the build's environment
#      OUTDIR      where the pictures go (default ../models/pictures)
set -e
HERE=$(cd $(dirname $0) && pwd)
export THICK_WORK=${THICK_WORK:-/tmp/thickness-pictures}
export RUN=${RUN:-$HOME/works/sw/fcad/.conda/run.sh}
FCBUILD=${FCBUILD:-$HOME/works/sw/fcad/build/conda-relwithdebinfo-801}
W=$THICK_WORK
mkdir -p $W/r $W/old $W/home
tables() { python3 -c "import sys; sys.path.insert(0, '$HERE'); import cases; $1"; }
LIBS="libTKBool.so.8.0.1 libTKOffset.so.8.0.1"

# The libraries: "up" and one per stage.
set -- $(tables "print(cases.UPSTREAM[0], ' '.join(cases.UPSTREAM[1]))")
[ -e $W/old/up/lib/libTKOffset.so.8.0.1 ] || $HERE/mkold.sh $1 $W/old/up "${@:2}" > $W/old/up.log
STAGES=$(tables "print(' '.join(cases.STAGES))")
for st in $STAGES; do
    commit=$(tables "print(cases.STAGES['$st'][0])")
    [ -e $W/old/$st/lib/libTKOffset.so.8.0.1 ] || $HERE/mkold.sh $commit $W/old/$st > $W/old/$st.log
done

# The cases: compute <libdir or -> <outdir> [VAR=value...]
compute() {
    local lib=$1 out=$2 pre=""
    shift 2
    [ "$lib" != "-" ] && pre=$(for l in $LIBS; do printf '%s ' $lib/$l; done)
    rm -rf $out
    env "$@" THICK_TOOLS=$HERE OUT=$out LD_PRELOAD="$pre" FREECAD_USER_HOME=$W/home \
        $RUN $FCBUILD/bin/FreeCADCmd $HERE/compute.py > $out.txt 2> $out.err
    grep -q '^DONE' $out.txt || { echo "compute $out failed, see $out.err"; exit 1; }
}
compute - $W/r/after STAGES=all
compute $W/old/up/lib $W/r/up STAGES=all
for st in $STAGES; do
    compute $W/old/$st/lib $W/r/$st STAGES=$st
done

# The pictures.
rm -rf $W/png
python3 $HERE/mkjobs.py
JOBS=$W/jobs.json FREECAD_USER_HOME=$W/home xvfb-run -a -s "-screen 0 1600x1200x24" \
    $RUN $FCBUILD/bin/FreeCAD $HERE/render.py > $W/render.log 2>&1
grep -q '^DONE' $W/render.log || { echo "render failed, see $W/render.log"; exit 1; }
$RUN python $HERE/compose.py ${OUTDIR:-$HERE/../models/pictures}
