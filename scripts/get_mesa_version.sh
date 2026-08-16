#!/bin/sh

awk -f- ${MESA_DIR}/data/version_number <<EOF
# Read MESA version number
NR==FNR {
    sub(/^r/, "", \$0)
    print
}
