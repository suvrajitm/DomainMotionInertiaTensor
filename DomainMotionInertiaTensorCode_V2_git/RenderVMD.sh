#!/bin/sh
#function RenderVMDT () {
echo "Running tachyon with AO lighting"

filename=$1
renderoption=$2
# TCB lab install location
# setenv TACHYONBIN /Projects/vmd/pub/linux64/lib/vmd186b13/tachyon_LINUXAMD64

# typical install location
#setenv TACHYONBIN /usr/local/lib/vmd/tachyon_LINUXAMD64 
#setenv TACHYONBIN /guam.raid.home/sbgrid/programs/x86_64-linux/vmd/1.9.2/lib/tachyon_LINUXAMD64

export TACHYONBIN='/guam.raid.home/sbgrid/programs/x86_64-linux/vmd/1.9.2/lib/tachyon_LINUXAMD64' 

if [ $renderoption = 1 ] ; then

    $TACHYONBIN -trans_vmd -res 1024 1024 -aasamples 4 -rescale_lights 0.4 -add_skylight 1.4 -format BMP $filename -o $filename.bmp 

elif [ $renderoption = 2 ] ; then
#### for axes only 

    $TACHYONBIN -trans_vmd -res 1024 1024 -aasamples 4 -rescale_lights 0.5 -add_skylight 0.8 -format BMP $1 -o $1.bmp 

    #$TACHYONBIN -trans_vmd -trans_max_surfaces 3 -res 1024 1024 -aasamples 4 -rescale_lights 0.3 -add_skylight 1.2 -format BMP $1 -o $1.bmp 

    #exec convert $1.tga $1.bmp

fi

#}