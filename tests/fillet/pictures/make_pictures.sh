#!/bin/bash
# make_pictures.sh -- regenerate models/pictures (README.md, "Pictures").
#
# Builds scratch TKFillet libraries for upstream's files and for the fork just
# before each fix (mkold.sh), runs every case of cases.py on each and on the
# fork as built (compute.py, the libraries loaded ahead of the installed
# ones), renders the panels in the FreeCAD GUI (render.py) and lays them out
# (compose.py). Linux (Xvfb, LD_PRELOAD) and macOS (the real display,
# DYLD_LIBRARY_PATH) alike.
#
# env: FILLET_WORK scratch directory (default /tmp/fillet-pictures; the
#                  libraries in it are reused when present)
#      FCBUILD     a FreeCAD build linked against this OCCT tree
#      OCCT_BUILD  this OCCT's build tree
#      RUN         the wrapper that runs a command in the build's environment
#      OUTDIR      where the pictures go (default ../models/pictures)
set -e
HERE=$(cd $(dirname $0) && pwd)
O=$(cd $HERE/../../.. && pwd)
export FILLET_WORK=${FILLET_WORK:-/tmp/fillet-pictures}
export RUN=${RUN:-$HOME/works/sw/fcad/.conda/run.sh}
MAC=$([ "$(uname)" = Darwin ] && echo 1 || true)
if [ -n "$MAC" ]; then
    FCBUILD=${FCBUILD:-$HOME/works/sw/fcad/build/mac-relwithdebinfo-801}
    export OCCT_BUILD=${OCCT_BUILD:-$O/build_conda_rwdi_801}
else
    FCBUILD=${FCBUILD:-$HOME/works/sw/fcad/build/conda-relwithdebinfo-801}
    export OCCT_BUILD=${OCCT_BUILD:-$O/build_conda_relwithdebinfo_801}
fi
W=$FILLET_WORK
mkdir -p $W/r $W/old $W/home
tables() { python3 -c "import sys; sys.path.insert(0, '$HERE'); import cases; $1"; }

# The libraries: "up" and one per stage. A directory made earlier is reused
# when it holds every toolkit its column needs.
have_libs() {  # have_libs <dir> <toolkit...>
    local d=$1 tk
    shift
    for tk in "$@"; do
        ls $d/lib/lib$tk.* > /dev/null 2>&1 || return 1
    done
}
UP=$(tables "print(cases.UPSTREAM)")
UPTKS=$(tables "print(' '.join(cases.UPSTREAM_TOOLKITS))")
if ! have_libs $W/old/up $UPTKS; then
    # upstream's TKFillet, and the files of other toolkits a fix touched
    UPFILES="$(cd $O && git diff --name-only $UP HEAD -- src/*/TKFillet)
             $(tables "print(' '.join(cases.UPSTREAM_EXTRA))")"
    TOOLKITS="$UPTKS" $HERE/mkold.sh $UP $W/old/up $UPFILES > $W/old/up.log
fi
STAGES=$(tables "print(' '.join(cases.STAGES))")
for st in $STAGES; do
    commit=$(tables "print(cases.STAGES['$st'][0])")
    tks=$(tables "print(cases.STAGE_TOOLKITS.get('$st', 'TKFillet'))")
    have_libs $W/old/$st $tks || TOOLKITS="$tks" $HERE/mkold.sh $commit $W/old/$st > $W/old/$st.log
done

# Run a command with the libraries of <libdir> (or - for the fork as built).
with_libs() {
    local lib=$1
    shift
    if [ "$lib" = "-" ]; then
        $RUN "$@"
    elif [ -n "$MAC" ]; then
        # SIP strips DYLD_* from /bin/bash, which RUN is: set it inside.
        $RUN env DYLD_LIBRARY_PATH=$lib "$@"
    else
        # every toolkit the directory holds
        LD_PRELOAD="$(ls $lib/*.so.8.0.1 | tr '\n' ' ')" $RUN "$@"
    fi
}

# The cases: compute <libdir or -> <outdir> [VAR=value...]
compute() {
    local lib=$1 out=$2
    shift 2
    rm -rf $out
    env "$@" FILLET_TOOLS=$HERE OUT=$out FREECAD_USER_HOME=$W/home \
        bash -c "$(declare -f with_libs); MAC=$MAC RUN=$RUN with_libs $lib $FCBUILD/bin/FreeCADCmd $HERE/compute.py" \
        > $out.txt 2> $out.err
    grep -q '^DONE' $out.txt || { echo "compute $out failed, see $out.err"; exit 1; }
}
compute - $W/r/after STAGES=all
compute $W/old/up/lib $W/r/up STAGES=all
for st in $STAGES; do
    compute $W/old/$st/lib $W/r/$st STAGES=$st
done

# The pictures.
rm -rf $W/png
$RUN python $HERE/compose.py jobs
if [ -n "$MAC" ]; then
    # a .py given to the GUI opens in the editor on macOS; a macro runs
    cp $HERE/render.py $W/render.FCMacro
    JOBS=$W/jobs.json FREECAD_USER_HOME=$W/home $RUN $FCBUILD/bin/FreeCAD $W/render.FCMacro \
        > $W/render.log 2>&1
else
    JOBS=$W/jobs.json FREECAD_USER_HOME=$W/home xvfb-run -a -s "-screen 0 1600x1200x24" \
        $RUN $FCBUILD/bin/FreeCAD $HERE/render.py > $W/render.log 2>&1
fi
grep -q '^DONE' $W/render.log || { echo "render failed, see $W/render.log"; exit 1; }
$RUN python $HERE/compose.py ${OUTDIR:-$HERE/../models/pictures}
