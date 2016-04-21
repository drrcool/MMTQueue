#!/bin/sh
# get_all.sh trimestername

# Wrapper to get both catalogs and masks

# Check to ensure that the trimester name was passed.
if test "$#" -ne 1; then
    echo "You must provide a trimester name"
fi
TRIMESTER=$1

sh get_mmirs_catalogs.sh $TRIMESTER
sh get_mmirs_masks.sh $TRIMESTER
