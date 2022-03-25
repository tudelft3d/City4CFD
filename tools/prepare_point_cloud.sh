#!/bin/bash

cd ${0%/*} || exit 1
WORKDIR="$PWD"

#!** Script to prepare point cloud for City4CFD reconstruction using LAStools **!#
#! What this script does:
    #- Merges multiple .las and .las tiles
    #- Removes points outside a bounding box
    #- Filters out ground and water points (LAS classes 2 and 9) in 'ground.xyz'
    #- Filters out buildings (LAS class 6) in 'buildings.xyz' 

#! The script will run for some time depending on the size of the point cloud

# USAGE:
#   ./prepare_point_cloud.sh filename, e.g. ./prepare_point_cloud.sh *.laz

#! PLEASE SET THE INPUT VARIABLES BEFORE RUNNING THE SCRIPT !#
#!---------------------------- INPUT VARIALBES --------------------------------#!
# Full path of the City4CFD source folder
CITY4CFD_DIR=${PWD%/*}

# Eastimate the bounding box of the dataset
# Min_x Min_y Max_x Max_y
BBOX="89973.860 435350.211 90632.581 435937.505"
#===============================================================================#

echo "Preparing point cloud using LAStools for City4CFD reconstruction"
echo "Working directory: $WORKDIR"
echo ""

lasmerge="$CITY4CFD_DIR/thirdparty/LAStools/bin64/lasmerge64"
las2las="$CITY4CFD_DIR/thirdparty/LAStools/bin64/las2las64"

mkdir pctempfiles

echo "Merging tiles and cropping area..."
echo "Files to be merged: $@"
echo "Bounding box: $BBOX"
$lasmerge -i $@ -o pctempfiles/temp.laz -keep_xy $BBOX -keep_class 2 6 9

cd pctempfiles
echo "Extracting buildings..."
$las2las -i temp.laz -keep_class 6 -oparse xyz -o buildings.xyz

echo "Extracting ground and water..."
$las2las -i temp.laz -keep_class 2 9 -oparse xyz -o ground.xyz

mv buildings.xyz ..; mv ground.xyz ..
cd ..; rm -r pctempfiles
echo "End"
