#!/bin/sh
# get_mmirs_masks.sh trimestername

# Fetches  mask files from fields and places them in the proper working directories.

# Check to ensure that the trimester name was passed.
if test "$#" -ne 1; then
    echo "You must provide a trimester name"
    exit 1
fi
TRIMESTER=$1
TRIMESTER_UPPERCASE=`echo $TRIMESTER | awk '{print toupper($0)}'`

# Check to see if the MMTQUEUE_PATH environment variable is set. Otherwise, assume a location
MMTQUEUE_PATH=${MMTQUEUE_PATH:="$HOME/MMTQueue/"}
CATALOG_PATH="$MMTQUEUE_PATH/catalogs/$TRIMESTER/"
REMOTE_PATH="/data/mmti/mmirs/masks/$TRIMESTER_UPPERCASE/*"

# Check that the target directory exists
if [ ! -d "$CATALOG_PATH" ]; then
    mkdir $CATALOG_PATH
fi

MASK_PATH="$CATALOG_PATH/masks/"
# Check to see that the mask path directory exists
if [ ! -d "$MASK_PATH" ]; then
    mkdir $MASK_PATH
fi

echo $REMOTE_PATH
echo "Starting rsync job.  You will likely need to enter the password for mmirs on fields"
rsync -e ssh -Cuav --progress mmirs@fields.mmto.arizona.edu:$REMOTE_PATH $MASK_PATH
