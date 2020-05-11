#!/bin/bash

# This program finds the slope/laplacian/variance/roughness for
# different resolutions for point data,
#calcDerivedVars is using the wgs84 coordinate system
# Here the sourcefile is HARP_dem_new.tif  with 100m resolutions.
#The interpolation method is average

# ./calcDerivedVars.sh sourcefile

sourcefile=$1
fbase=`echo "${sourcefile%.*}"`
ext1=`echo "${sourcefile##*.}"`
ext2="nc"

#This will not correctly give the filebase and extension in the case of HARP_dem_new_2.5km.tif
#fbase=`echo "$sourcefile" | cut -d'.' -f1`
#ext1=`echo "$sourcefile" | cut -d'.' -f2`

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


# Creating directories for slope and roughness and laplacian, variance not derived now

for derived  in "slope" "roughness" "tpi" "tri"  "laplacian" "variance";do
    mkdir -p $path/$derived
done	   


# Associative array for different resolutions. NOTE: always check the resolutions are dividable to the extents
declare -A arr

arr["100m"]=100
arr+=(["200m"]=200 ["500m"]=500 ["1km"]=1000 ["2.5km"]=2500 ["5km"]=5000 ["7.5km"]=7500 ["10km"]=10000 ["15km"]=15000 ["30km"]=30000 )


#Since we find the derived vars in lcc coordinates we use the value in the arr as the resolution.



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
