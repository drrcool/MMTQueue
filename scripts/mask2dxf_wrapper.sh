#!/bin/bash
# Wrapper to run the MMIRS task mask2dxf on all .msk file in the 
# current working directory.

# If mask2dxf is located somewhere else, you may need to change the
# direct path below.
MMTQUEUE_PATH=${MMTQUEUE_PATH:="$HOME/MMTQueue/"}
binfile=$MMTQUEUE_PATH/scripts/mask2dxf

# Loop through each of the msk files
for filename in *.msk
do 
    fileroot="${filename%.*}"
    $binfile < $fileroot.msk > $fileroot.dxf
done
