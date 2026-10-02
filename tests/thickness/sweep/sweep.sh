#!/bin/bash
# sweep.sh <outfile> [libdir]  -- run sweep.py to the end, restarting after a
# crash (the run it died in is written as CRASH). With libdir, its TKBool and
# TKOffset are preloaded: the libraries of another build, upstream's or an
# older fork's (pictures/mkold.sh makes them).
# env: FCBUILD  a FreeCAD build linked against this OCCT tree
#      RUN      the wrapper that runs a command in the build's environment
#      SWEEP_HOME  a scratch FREECAD_USER_HOME (default /tmp/thickness-sweep-home)
#      MODES, THICK_FREEZE  as sweep.py
HERE=$(cd $(dirname $0) && pwd)
RUN=${RUN:-$HOME/works/sw/fcad/.conda/run.sh}
FCBUILD=${FCBUILD:-$HOME/works/sw/fcad/build/conda-relwithdebinfo-801}
H=${SWEEP_HOME:-/tmp/thickness-sweep-home}
mkdir -p $H
OUT=$(realpath $1)
PRE=""
[ -n "$2" ] && PRE="$2/libTKBool.so.8.0.1 $2/libTKOffset.so.8.0.1"
: > $OUT
START=0
while :; do
    START=$START LD_PRELOAD="$PRE" FREECAD_USER_HOME=$H \
        $RUN $FCBUILD/bin/FreeCADCmd $HERE/sweep.py >> $OUT 2> $OUT.err
    grep -q '^DONE' $OUT && break
    line=$(grep '^CASE' $OUT.err | tail -1)
    [ -z "$line" ] && { echo "sweep did not start, see $OUT.err"; exit 2; }
    last=$(echo "$line" | awk '{print $2}')
    echo "R $last ${line#CASE $last } CRASH" >> $OUT
    START=$((last + 1))
done
