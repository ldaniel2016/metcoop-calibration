#!/bin/bash

# This program finds the topological variables such as
# slope/laplacian/variance/# roughness/tpi/tri for different
# resolutions from Digital Elevation Maps(DEM) file.  The topological
# variables are derived based on the HARP_dem_new.tif, a Digital
# Elevation Map (DEM) in the WGS-84 projection and with 100m
# resolution. It is then converted to the CMEPS projection which is
# Lambert azimuthal equal-area and polar stereographic map projection
# and it is proj=lcc +lat_0=63.3 +lon_0=15 +lat_1=63.3 +lat_2=63.3
# +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs.

# The area data in geotiff format use the lower left corner of the
# cells to specify coordinates, while NWP grids are centered. To
# correct this we need # to subtract half the pixel size for the lower
# x and y coordinates and add half# the pixel size for the upper x The
# DEM file converted to the MEPS projection has the coordinates xmin =
# -1060084, xmax = 1309916, ymin = -1332518 and ymax = 1337482.  As
# the extents of the CMEPS area are xmin = -1162860.0408, xmax =
# 1196095.6182, ymin= -1023334.1885 and ymax = 1207415.0381, we lose a
# few stations in the North East area as the DEM does not cover it.

# The terrain types are derived from the DEM file,  HARP_globcovernew.tif.
# This file also converted to MEPS projection

# The DEM file in the MEPS domain is wrapped to different resolutions
# (or pixel sizes) with the interpolation method 'average' using
# gdalwrap commands


# Usage  PATH/calcDerivedVars_meps.sh sourcefile
sourcefile=$1
fbase=`echo "${sourcefile%.*}"`
ext1=`echo "${sourcefile##*.}"`
ext2="nc"

path="$(dirname $sourcefile)"
scriptsdir=$(echo $path/scripts)

mepsourcefile=$(echo ${fbase}_meps.${ext1})
mepsfbase=`echo "${mepsourcefile%.*}"`

gcfile=$(echo $path/HARP_globcover_new.tif)
mepsgcfile=$(echo $path/HARP_globcover_new_meps.tif)
mepsgcfbase=`echo "${mepsgcfile%.*}"`

#Converting HARP_dem_new.tif and HARP_globcover_new.tif  from WGS84 to LCC 

gdalwarp -s_srs "+proj=longlat +datum=WGS84 +no_defs" -t_srs "+proj=lcc +lat_0=63.3 +lon_0=15 +lat_1=63.3 +lat_2=63.3 +no_defs +R=6.371e+06" -r average   $sourcefile $mepsourcefile

gdalwarp -s_srs "+proj=longlat +datum=WGS84 +no_defs" -t_srs "+proj=lcc +lat_0=63.3 +lon_0=15 +lat_1=63.3 +lat_2=63.3 +no_defs +R=6.371e+06" -r average   $gcfile $mepsgcfile


# Creating directories for slope and roughness and laplacian, variance, tpi (Topographic Position Index) and tri (Terrain Ruggedness Index)

for derived  in "slope" "roughness" "tpi" "tri" "laplacian" "variance";do
    mkdir -p $path/$derived
done	   


# Associative array for different resolutions. NOTE: always check the resolutions are dividable to the extents
declare -A arr

arr["100m"]=100
arr+=(["200m"]=200 ["500m"]=500 ["1km"]=1000 ["2.5km"]=2500 ["5km"]=5000 ["7.5km"]=7500 ["10km"]=10000 ["15km"]=15000 ["30km"]=30000 )


#Since we find the derived vars in lcc coordinates we use the value in the arr as the resolution.


# The projected coordinates of the DEM file  in the MEPS projection
xmin=-1060084
ymin=-1332518
xmax=1309916
ymax=1337482

for key in ${!arr[@]}; do
    val=${arr[${key}]}
    resolution=$val
    ((shift=resolution/2))
    echo $resolution $shift
    infile=$(echo ${mepsfbase}_${key}.${ext1})
    echo $infile
    gcinfile=$(echo ${mepsgcfbase}_${key}.${ext1})
    ncfile=$(echo ${mepsfbase}_${key}.${ext2})
    gcncfile=$(echo ${mepsgcfbase}_${key}.${ext2})
    outfileslope=$(echo $path/slope/gdalslope_$(basename $infile))
    outfilerough=$(echo $path/roughness/roughness_$(basename $infile))
    outfiletpi=$(echo $path/tpi/tpi_$(basename $infile))
    outfiletri=$(echo $path/tri/tri_$(basename $infile))

    
    
    # using gdalwarp to generate a file with a new resolution
    # using gdal_translate to convert the GeoTIFF to netCDF format
    if [ ! -e "$infile" ] ; then
	gdalwarp -co COMPRESS=DEFLATE -r average  -te $((xmin-shift)) $((ymin-shift)) $((xmax+shift)) $((ymax+shift))  -tr $resolution $resolution $mepsourcefile $infile
    fi

    if [ ! -e "$gcinfile" ] ; then
	gdalwarp -co COMPRESS=DEFLATE -r average  -te $((xmin - shift)) $((ymin - shift)) $((xmax + shift)) $((ymax + shift))  -tr $resolution $resolution $mepsgcfile $gcinfile
    fi
    
    if [ ! -e "$ncfile" ] ; then
	gdal_translate -of netCDF $infile $ncfile
    fi

        if [ ! -e "$gcncfile" ] ; then
	gdal_translate -of netCDF $gcinfile $gcncfile
    fi
  
   
    python2 $scriptsdir/calc_centeredSlope_scale.py $infile 1
    python2 $scriptsdir/calcLaplacian.py $infile
    python2 $scriptsdir/calcVariance.py $infile
    gdaldem roughness  $infile $outfilerough
    gdaldem TPI  $infile $outfiletpi
    gdaldem TRI  $infile $outfiletri
    gdaldem slope -p $infile $outfileslope
   
done

# Converting to netCDF files
for derived  in "slope" "roughness" "tpi" "tri" "laplacian" "variance";do
    dir=$path/$derived
    for file in $dir/*; do
	tifbase=`echo "${file%.*}"`
	ncfile=$(echo ${tifbase}.${ext2})
	if [ ! -e "$ncfile" ] ; then
	    gdal_translate -of netCDF $infile $ncfile
	fi
    done
done
