#!/bin/bash

# loading modules
echo "...loading modules."
module load anaconda

# creating and activating conda environment
echo "...activating conda environment."
conda activate pmed

# [] Requries reference map files to be in the order:
# col1: chromosome (formatted "chr1")
# col2: variant identifier
# col3: position in centimorgans
# col4: base-pair coordinate
# If you have a map file whose ordering is different,
# this script provides a suggested way to re-order,
#
# NOTE: col1: variant identifier is not used,
# so you can pad with '.' or any data.
#
# NOTE: if your first column does not have "chr"
# before the chromosome number, you can add "chr" in front of $1.
#
# Assume the unordered map file you have comes as:
# col1: chromosome
# col2: position in centimorgans
# col3: base-pair coordinate
# col4: variant identifier

# TODO: fill in file paths appropriately.
unordered_reference_map=""
orderd_reference_map=""

awk '{ print $1, $4, $2, $3}' $unordered_reference_map > $orderd_reference_map
awk '{ print "chr"$1, $4, $2, $3}' $unordered_reference_map > $orderd_reference_map
