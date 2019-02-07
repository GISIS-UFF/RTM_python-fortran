#!/bin/bash
# LowPass Filter in Migrated Section with Seismic Unix
# file: filter_StackSection.sh
#----------------------------------------------------------
# Uses Seismic Unix programs to pass low bandpass filter in 
# Seismic Section. 
# The output is  Filtered Section with different
# cut frequencies. This cut frequencies goes from $f_cutmin
# to $f_cutmax incremented by $df.
# 1 - Convert binaries to SU format
# 2 - Set ns and dt
# 3 - Low Pass Filter in Seismic Section with f_cutmax cut frequency
# 4 - Remove SU header
# 5 - Save Seismogram with cut frequency name
#==========================================================
#example: ./filter_StacktSection.sh Nz    h   fcut
#         ./filter_StacktSection.sh 501  10    20
#
# INPUT:  
# 
# samples       = Total Number of Samples in Time
# sampling      = Sample rate of Seismogram
# f_cutmax      = Max Frequency permittedl receive Seismogram
 
# OUTPUT: ../Migrateddata/prefix_MigrateddataStack_f_cutHz.bin
 
# Code Written by Felipe Timoteo
#                 Last update: Jan 27th, 2018
 
# Copyright (C) 2018 Grupo de Imageamento Sísmico e Inversão Sísmica (GISIS)
#                    Departamento de Geologia e Geofísica
#                    Universidade Federal Fluminense
# Not Finished yet!!!
#------------------------------------------------
#set -x


#Plot Parameters
Wplot=1800       # Width of plot (pixels)
Hplot=900       # Height of plot (pixels)

plotymax=50 # y axis limit (frequency spectrum)

# Data Parameters
indata=../Imagem/Marmousi_FinalImage.bin

#indata=../Models_true/velocitymodel.bin
Nz=141
dh=10

# sample rate in micro in SU # 
sampling=`python -c "print(float('$dh')*float('1.0e3'))"`

# Cut Frequency
f_cut=6
f_cut_max=30

# End of filter taper   
d_fcut=5

sufix=_SUfilter_
Hz=Hz
outdata=$indata$sufix$f_cut-$f_cut_max$Hz.bin

#plot Migrated Seismic Section
suaddhead < $indata ns=$Nz | suchw key1=dt a=$sampling |\
            suximage  perc=98 title="Migrated Seismic Section" legend=1 units=" Amplitude  "\
            xbox=10 ybox=10 wbox=$Wplot hbox=$Hplot  legend=1 label1="Depth (km)" label2="Lenth (m)" &

#plot Filtered Frequency Spectrum of Migrated Seismic Section
suaddhead < $indata ns=$Nz | suchw key1=dt a=$sampling |           
           sufilter f=$f_cut,$(($f_cut+$d_fcut)),$(($f_cut_max-$d_fcut)),$f_cut_max amps=0.,1.,1.,0.| #BANDPASS FILTER           
           #sufilter f=$(($f_cut-$d_fcut)),$f_cut amps=1.,0.| #LOWPASS FILTER             
           suspecfx | #suop op=norm | #SHOW AMPLITUDE SPECTRUM  
           suximage  perc=98 title="Frequency Spectrum of Migrated Seismic Section with band filter with $f_cut and $f_cut_max Hz" \
           xbox=10 ybox=10 wbox=$Wplot hbox=$Hplot  \
           legend=1 label1="Depth (km)" label2="Lenth (m)"&# cmap=hsv2 &


echo "***********************"
echo "Normal end of Execution"
echo "***********************"

