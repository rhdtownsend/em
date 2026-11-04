#!/bin/bash
#
# Configure and check pkgconf file for ForUM

IN_FILE=$1
OUT_FILE=$2

PKG_DIR=`dirname $OUT_FILE`

# Check that GYRE_DIR is set

if [[ -z "$GYRE_DIR" ]]; then
    echo "GYRE_DIR is not set"
    exit 1
fi

# Customize the	pkgconfig file

awk -v prefix=$GYRE_DIR '{sub(/%PREFIX%/, prefix); print}' \
    $IN_FILE > $OUT_FILE

# Check that ForUM has the correct version

version=`awk '/VERSION = .*/ { gsub(/'\''/, "", $NF); print $NF }' $GYRE_DIR/src/forum/src/lib/version_m.fypp`

pkgconf --with-path=$PKG_DIR --exists "forum = $version"

if [ $? -ne 0 ]; then

    req_version=`pkgconf --with-path=$PKG_DIR --modversion forum`

    echo "The version ($version) of ForUM at $GYRE_DIR does not match the required version ($req_version)"

    rm $OUT_FILE

    exit 1

fi
