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

plotymax=150


# Data Parameters
indata=../Migrateddata/RTM_Mrm2invrtd_MigrateddataStack.bin

#indata=../Models_true/velocitymodel.bin
Nz=351
dh=10.0e3

# Cut Frequency
f_cut=5
f_cut_max=60

# End of filter taper   
d_fcut=2        
sufix=_SUfilter_
Hz=Hz
outdata=$indata$sufix$f_cut-$f_cut_max$Hz.bin

#plot Migrated Seismic Section
suaddhead < $indata ns=$Nz | suchw key1=dt a=$dh |\
            suximage  perc=98 title="Migrated Seismic Section" legend=1 units=" Amplitude  "\
            xbox=10 ybox=10 wbox=$Wplot hbox=$Hplot  legend=1 label1="Depth (km)" label2="Lenth (m)" &

#plot Frequency Spectrum of Migrated Seismic Section
suaddhead < $indata ns=$Nz | suchw key1=dt a=$dh |\
            suspecfx |  suop op=norm | suwind tmin=0 tmax=$plotymax|\
            suximage  perc=98 title="Frequency Spectrum of Migrated Seismic Section" \
            xbox=10 ybox=10 wbox=$Wplot hbox=$Hplot  \
            legend=1 label1="Depth (km)" label2="Lenth (m)" cmap=hsv2 &




# #plot Filtered Frequency Spectrum of Migrated Seismic Section   - LOW PASS
# suaddhead < $indata ns=$Nz | suchw key1=dt a=$dh |\
#             sufilter f=$(($f_cut-$d_fcut)),$f_cut amps=1.,0.|\
#             suspecfx |  suop op=norm | suwind tmin=0 tmax=$plotymax|\
#             suximage  perc=98 title="Frequency Spectrum of Migrated Seismic Section with lowpass filter of $f_cut Hz" \
#             xbox=10 ybox=10 wbox=$Wplot hbox=$Hplot  \
#             legend=1 label1="Depth (km)" label2="Lenth (m)" cmap=hsv2 &

# #plot Migrated Seismic Section    - LOW PASS
# suaddhead < $indata ns=$Nz | suchw key1=dt a=$dh |\
#             sufilter f=$(($f_cut-$d_fcut)),$f_cut amps=1.,0.|\
#             suximage  perc=98 title="Migrated Seismic Section with lowpass filter of $f_cut Hz" \
#             xbox=10 ybox=10 wbox=$Wplot hbox=$Hplot  legend=1 label1="Depth (km)" label2="Lenth (m)" &




#plot Filtered Frequency Spectrum of Migrated Seismic Section    - BAND PASS
suaddhead < $indata ns=$Nz | suchw key1=dt a=$dh |\
           sufilter f=$(($f_cut-$d_fcut)),$f_cut,$f_cut_max,$(($f_cut_max+$d_fcut)) amps=0.,1.,1.,0.|\
           suspecfx |  suop op=norm | suwind tmin=0 tmax=$plotymax|\
           suximage  perc=98 title="Frequency Spectrum of Migrated Seismic Section with band filter with $f_cut and $f_cut_max Hz" \
           xbox=10 ybox=10 wbox=$Wplot hbox=$Hplot  \
           legend=1 label1="Depth (km)" label2="Lenth (m)" cmap=hsv2 &

#plot Migrated Seismic Section    - BAND PASS
suaddhead < $indata ns=$Nz | suchw key1=dt a=$dh |\
            sufilter f=$(($f_cut-$d_fcut)),$f_cut,$f_cut_max,$(($f_cut_max+$d_fcut)) amps=0.,1.,1.,0.|\
            suximage  perc=98 title="Migrated Seismic Section with lowpass filter of $f_cut Hz" \
            xbox=10 ybox=10 wbox=$Wplot hbox=$Hplot  legend=1 label1="Depth (km)" label2="Lenth (m)" &

# #plot Migrated Seismic Section
# echo ""
# echo "Saving Filtered Seismic as: $outdata"
# echo ""
# suaddhead < $indata ns=$Nz | suchw key1=dt a=$dh |\
#             sufilter f=$(($f_cut-$d_fcut)),$f_cut,$f_cut_max,$(($f_cut_max+$d_fcut)) amps=0.,1.,1.,0.|\
#             #sustrip > $outdatab
#               suximage  perc=95 title="Frequency Spectrum of Migrated Seismic Section with band filter with $f_cut and $f_cut_max Hz" \
#              xbox=10 ybox=10 wbox=$Wplot hbox=$Hplot  legend=1 label1="Depth (km)" labe2l2="Lenth (m)" &
echo "***********************"
echo "Normal end of Execution"
echo "***********************"

