#!/usr/bin/env bash
#
# Barcode
# Copyright E.G.P. Bos and F.S. Kitaura
#
# Distributed under the terms of the MIT License.
# The full license is in the file LICENSE, distributed with this software.
#

PROJECT_BINARY_DIR=$1
OUTNAME=$2

#echo "$0"
#echo "$1"
#echo "$2"

OLD_PATH="${PROJECT_BINARY_DIR}/${OUTNAME}"

LATEST_PATH="${PROJECT_BINARY_DIR}/${OUTNAME}_latest_build"

HG_REV=$(hg parents --template '{rev}')

if [[ $(hg status 2>/dev/null) ]]
then
    NEW_PATH="${PROJECT_BINARY_DIR}/${OUTNAME}_r${HG_REV}"_uncommitted
else
    NEW_PATH="${PROJECT_BINARY_DIR}/${OUTNAME}_r${HG_REV}"
fi

mv "$OLD_PATH" "$NEW_PATH"
rm "$LATEST_PATH"
ln -s "$NEW_PATH" "$LATEST_PATH"

#echo "Renamed output binary ${OLD_PATH} to ${NEW_PATH} and created symbolic link barcode_latest_build to new path."
