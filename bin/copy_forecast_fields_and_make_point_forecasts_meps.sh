#!/bin/bash

# This program copies the fields that are used in generating point
# forecasts to a common folder and changes their name so that they are
# compatible with the extract_values.py -script and the .dbf -format
# produced by it.

# Usage: $PATH/copy_forecast_fields_and_make_point_forecasts_meps.sh

outdir="/home/daniel/metcoop_calibration/derived_terrain_variables/DEM_fields/fields_interpolated_to_points/"
mkdir -p $outdir


cd /home/daniel/metcoop_calibration/derived_terrain_variables/elevation
for f in *.tif
do
    f2=${f/HARP_dem_new_meps_/elev}
    cp $f $outdir$f2
done


cd /home/daniel/metcoop_calibration/derived_terrain_variables/globcover
for f in *.tif
do
    f2=${f/HARP_globcover_new_meps_/gc}
    cp $f $outdir$f2
done


cd /home/daniel/metcoop_calibration/derived_terrain_variables/laplacian
for f in *.tif
do
    f2=${f/laplacian_HARP_dem_new_meps_/lap}
    cp $f $outdir$f2
done

cd /home/daniel/metcoop_calibration/derived_terrain_variables/roughness
for f in *.tif
do
    f2=${f/roughness_HARP_dem_new_meps_/rough}
    cp $f $outdir$f2
done

cd /home/daniel/metcoop_calibration/derived_terrain_variables/tpi
for f in *.tif
do
    f2=${f/tpi_HARP_dem_new_meps_/tpi}
    cp $f $outdir$f2
done

cd /home/daniel/metcoop_calibration/derived_terrain_variables/tri
for f in *.tif
do
    f2=${f/tri_HARP_dem_new_meps_/tri}
    cp $f $outdir$f2
done

cd /home/daniel/metcoop_calibration/derived_terrain_variables/slope

for f in *.tif
do
    f2=${f/slope_NS_ABS_HARP_dem_new_meps_/NSABS}
    f2=${f2/slope_WE_ABS_HARP_dem_new_meps_/WEABS}
    f2=${f2/slope_Z_ABS_HARP_dem_new_meps_/ZABS}
    f2=${f2/gdalslope_HARP_dem_new_meps_/ZPRCT}
    cp $f $outdir$f2
done

cd /home/daniel/metcoop_calibration/derived_terrain_variables/variance

for f in *.tif
do
    f2=${f/variance_HARP_dem_new_meps_/var}
    cp $f $outdir$f2
done

# The program extract_values.py extracts the values of different
# spatial variables to the station points given in the shape file. The
# shape file and the output directory are to be given as arguments.

python2 /home/daniel/metcoop_calibration/derived_terrain_variables/scripts/extract_values.py -c  /home/daniel/metcoop_calibration/derived_terrain_variables/metcoop_station_list_new_MEPS_translated/meps_stations_mepsCRS_new.shp -d $outdir -e tif

