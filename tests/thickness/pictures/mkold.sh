#!/bin/bash
# mkold.sh <commit> <dir> [path...]
#
# Scratch TKBool/TKOffset libraries with the files that differ between
# <commit> and HEAD taken as of <commit> -- or only the given paths -- and
# everything else as the work tree has it now. Compiled from the build tree's
# own compile commands (oldbuild.py) into <dir>/lib, for LD_PRELOAD.
set -e
HERE=$(cd $(dirname $0) && pwd)
C=$1; D=$2; shift 2
case "$D" in /*) ;; *) D=$PWD/$D ;; esac
O=$(cd $HERE/../../.. && pwd)
if [ -e "$D" ]; then
    # only ever remove a directory this script made
    [ -d "$D/lib" ] && [ -d "$D/obj" ] || { echo "$D is not a library directory of mkold.sh"; exit 1; }
    rm -rf "$D"
fi
mkdir -p $D/src $D/inc $D/lib $D/obj
cd $O
if [ $# -gt 0 ]; then
  FILES="$*"
else
  FILES=$(git diff --name-only $C HEAD -- src/ModelingAlgorithms/TKBool src/ModelingAlgorithms/TKOffset)
fi
for f in $FILES; do
  rel=${f#src/}
  mkdir -p $D/src/$(dirname $rel)
  git show $C:$f > $D/src/$rel
  case $f in *.hxx) cp $D/src/$rel $D/inc/ ;; esac
done
U=$D OCCT=$O ${RUN:-$HOME/works/sw/fcad/.conda/run.sh} python $HERE/oldbuild.py
