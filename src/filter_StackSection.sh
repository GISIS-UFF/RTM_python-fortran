#!/bin/bash
# LowPass Filter in Migrated Section with Seismic Unix
# file: filter_StackSection.sh
#----------------------------------------------------------
# Uses Seismic Unix programs to pass low bandpass filter in 
# Seismic Section. Plot the T-X spectrum.
# Code Written by Felipe Timoteo
#                 Last update: Jun 18th, 2019
 
# Copyright (C) 2018 Grupo de Imageamento Sísmico e Inversão Sísmica (GISIS)
#                    Departamento de Geologia e Geofísica
#                    Universidade Federal Fluminense
#--------------------------------------------------------------------------#

#######    USER AREA     ##########

# Data Input
indata=../Imagem/Marmousi_v2_FinalImage_Laplacian.bin

Nz=281     # number of samples of 1st dimension (z direction or time)
dh=10.0e3  # sample rate (dt or dh). SU in miliseconds or milimeters

# Cut Frequency
#             ___________
#           /|  
#          / | 
#_________/  | 
#         | fcut
#         |                
#     fcut-d_fcut           
f_cut=7
d_fcut=3    # shift at begininng and end of the filter

# Data Output
sufix=_highpass_
Hz=Hz
outdata=$indata$sufix$((f_cut-d_fcut))-$f_cut$Hz.bin

#Plot Parameters
Wplot=1800       # Width of plot (pixels)
Hplot=900        # Height of plot (pixels)
plotymax=150     # max y direction

####### END OF USER AREA ##########
MENSAGEM_USO="
$(basename $0) [-p | -s ]
-p only plot and check the results
-s save the results on the disk

First, check the parameters in 
user area inside the 
$0 script.

Actual parameters:
Cut frequency                 = $f_cut
slope frequency               = $d_fcut
number of samples in vertical = $Nz
sample rate                   = $dh
output file                   = $outdata

Usage:     
	bash$ $0 -p # Plot and check the results
	bash$ $0 -s # Save the results on the disk
    
"
# Check if user provide some parameter
[ -z "$1" ] && {
    echo "" 
	echo -e "$MENSAGEM_USO"
    echo "" 
	exit 1 
}

case "$1" in
-p | --plot) ## Plot the results
MESSAGE="
        Plotting the Migrated data, 
        Filterer Migrated 
        data and spectrum
"
    echo -e $MESSAGE
	#########################################################
    #plot Migrated Seismic Section
    suaddhead < $indata ns=$Nz |                     # add to binary SU header
                suchw key1=dt a=$dh |                # add to header sample rate
                suximage  perc=98 \
                title="Migrated Seismic Section" \
                legend=1 units=" Amplitude  "\
                xbox=10 ybox=10 wbox=$Wplot\
                hbox=$Hplot  legend=1 \
                label1="Grid points"\
                label2="Grid points" &                 # generate plot
    #########################################################
    #plot Frequency Spectrum of Migrated Seismic Section
    suaddhead < $indata ns=$Nz |               # add to binary SU header
                suchw key1=dt a=$dh |          # add to header sample rate
                suspecfx |                     # calculates amplitude spectrum 
                suop op=norm |                 # improve visualization (normalize)
                suwind tmin=0 tmax=$plotymax| # limit the data window visualizations
                suximage  perc=98 \
                title="Frequency Spectrum of Migrated Seismic Section" \
                xbox=10 ybox=10 wbox=$Wplot \
                hbox=$Hplot  \
                legend=1 label1="Frequency (Hz)" \
                label2="Grid points" \
                cmap=hsv2 &                    # generate plot                
    #########################################################
    #plot Migrated Seismic Section - HIGH PASS
    suaddhead < $indata ns=$Nz |               # add to binary SU header
                suchw key1=dt a=$dh |          # add to header sample rate
                sufilter f=$(($f_cut-$d_fcut)),$f_cut amps=0.,1.| # apply filter 
                suximage  perc=98 \
                title="Migrated Seismic Section with lowpass filter of $((f_cut-d_fcut))-$f_cut Hz" \
                xbox=10 ybox=10 wbox=$Wplot\
                hbox=$Hplot  legend=1\
                label1="grid points" \
                label2="grid points" &           # generate plot
    #########################################################
    #plot Filtered Frequency Spectrum of Migrated Seismic Section   - HIGH PASS
    suaddhead < $indata ns=$Nz |               # add to binary SU header
                suchw key1=dt a=$dh |          # add to header sample rate
                sufilter f=$(($f_cut-$d_fcut)),$f_cut amps=0.,1.| # apply filter 
                suspecfx |                     # calculates amplitude spectrum 
                suop op=norm |                 # improve visualization (normalize)
                suwind tmin=0 tmax=$plotymax|  # limit the data window visualizations
                suximage  perc=98 title="Frequency Spectrum of Migrated Seismic Section with lowpass filter of $((f_cut - d_fcut))-$f_cut Hz" \
                xbox=10 ybox=10 wbox=$Wplot hbox=$Hplot  \
                legend=1 label1="Depth (km)" label2="Lenth (m)" cmap=hsv2 & # generate plot
;;

-s | --save) # save the filtered data
MESSAGE="
save filtered data in : $outdata
"
    echo -e $MESSAGE
    #plot Migrated Seismic Section    - BAND PASS
    suaddhead < $indata ns=$Nz |               # add to binary SU header
                suchw key1=dt a=$dh |          # add to header sample rate
                sufilter f=$(($f_cut-$d_fcut)),$f_cut amps=0.,1.| # apply filter                 
                sustrip > $outdata
;;
esac

echo ""
echo "***********************"
echo "Normal end of Execution"
echo "***********************"
echo ""
exit 0