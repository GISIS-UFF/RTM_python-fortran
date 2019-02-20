#!/bin/bash
# Smoothing velocity model -  smooth_velocity.sh:
# This script uses Seismic Unix programs to Smooth Velocity Model.
# The velocity model is converted to vagarosity and then
# smoothed and finally converted again to velocity.
# This scripr is used during update velocity field to prepare
# to next iteration
#
# EXEMPLE OF USAGE:
# $> ./smooth_velocity.sh dh Nz  Nx    s_factor
# $> ./smooth_velocity.sh 5  501 1201  50
#
# INPUT:  
# file     = ../Models_initial/velocitymodel.bin 
# dh       = Incremental Space
# Nz       = Number of elements in z Direction
# Nx       = Number of elements in x Direction
# s_factor = smooth factor
#
# OUTPUT: 
# file     =  ../Models_initial/velocitymodel_smth.bin
#
# Code Written by Felipe Timoteo
#                 Last update: August 07, 2017
#
# Copyright (C) 2017 Grupo de Imageamento Sísmico e Inversão Sísmica (GISIS)
#                    Departamento de Geologia e Geofísica
#                    Universidade Federal Fluminense
#
#Set messages on
#set -x
#------------------------------------------------

dh=10 #$1     # dt in seismic unix is in microsecond 
dh=`python -c "print( float('$dh')*float('1.0e3'))"`

# number of points of model
Nz=141 #$2       #501
Nx=383 #$3       #1201

# Smooth Factor
s_factor=20 #$4 #100

# Check Input Parameters
echo $dh
echo $Nx
echo $Nz
echo $s_factor

# Path Folder of Model
pathfolder=../modelos_utilizados 
filename_input=marmousi_vp_383x141
sufixname=_SmoothedModel$s_factor
filename_output=$filename_input$sufixname

filename_output_brute=velocitymodel_inputmodel
indata=$pathfolder/$filename_input.bin
outdata=$pathfolder/$filename_output.bin
outdatabrute=$pathfolder/$filename_output_brute.bin

Wplot=800       # Width of plot (pixels)
Hplot=600        # Height of plot (pixels)

# Plot Input Velocity Model
suaddhead < $indata ns=$Nz | suchw key1=dt a=$dh |\
    suximage d2=$dh title="Input Velocity Model"  \
    xbox=$Wplot+100 ybox=10 wbox=$Wplot hbox=$Hplot cmap=hsv2\
    legend=1 label1="depth (km)" label2="distance (m)" &

# Obtain vagarosity Model
suaddhead < $indata ns=$Nz | suchw key1=dt a=$dh | suop op=inv |\
    sustrip   > tmp2

# Smoothing vagarosity model and convert to velocity model
smooth2 < tmp2 n1=$Nz n2=$Nx r1=$s_factor r2=$s_factor |\
    suaddhead ns=$Nz |suchw  key1=dt a=$dh | suop op=inv |\
    sustrip  > $outdata

# Plot Final Result
smooth2 < tmp2 n1=$Nz n2=$Nx r1=$s_factor r2=$s_factor |\
    suaddhead ns=$Nz |suchw  key1=dt a=$dh | suop op=inv |\
   suximage d2=$dh title="Smoothed Velocity Model"  \
    xbox=1200 ybox=10 wbox=$Wplot hbox=$Hplot cmap=hsv2 \
    legend=1 label1="depth (km)" label2="distance (m)" &

echo ""
echo ""
echo "Saving Smoothed model on: $outdata"
echo ""
rm -f tmp*

exit
