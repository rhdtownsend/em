#!/bin/bash
#
# Configure and check pkgconf file for MESA

IN_FILE=$1
OUT_FILE=$2

PKG_DIR=`dirname $OUT_FILE`

# Check that MESA_DIR is set

if [[ -z "$MESA_DIR" ]]; then
    echo "MESA_DIR is not set"
    exit 1
fi

# Customize the	pkgconfig file

awk -v prefix=$MESA_DIR '{sub(/%PREFIX%/, prefix); print}' \
    $IN_FILE > $OUT_FILE

# Check that MESA has the correct version

version=`awk 'NR==FNR {sub(/^r/, ""); print}' $MESA_DIR/data/version_number`

pkgconf --with-path=$PKG_DIR --exists "mesa = $version"

if [ $? -ne 0 ]; then

    req_version=`pkgconf --with-path=$PKG_DIR --modversion mesa`

    echo "The version ($version) of MESA at $MESA_DIR does not match the required version ($req_version)"

    rm $OUT_FILE

    exit 1

fi
