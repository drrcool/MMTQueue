#!/bin/sh
# get_mmirs_catalogs.sh trimestername

# Fetches catalogs from hopper and mask files from fields and places them in the proper working directories.

# Check to ensure that the trimester name was passed.
if test "$#" -ne 1; then
    echo "You must provide a trimester name"
    exit 1
fi


# Check to see if the MMTQUEUE_PATH environment variable is set. Otherwise, assume a location
MMTQUEUE_PATH=${MMTQUEUE_PATH:="$HOME/MMTQueue/"}
CATALOG_PATH="$MMTQUEUE_PATH/catalogs/$1/"
REMOTE_PATH="/var/www/mmtobs/hecto/catalogs/$1/mmirs/"

# Check that the target directory exists
if [ ! -d "$CATALOG_PATH" ]; then
    mkdir $CATALOG_PATH
fi

echo "Starting rsync job.  You will likely need to enter the password for mmtobs on hopper"
rsync -e ssh -Cuav --progress mmtobs@hopper.si.edu:$REMOTE_PATH $CATALOG_PATH
