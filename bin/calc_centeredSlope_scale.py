#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This program calculates the slope in X and Y directions of a DEM

# Usage: python2 $PATH/calc_centeredSlope_scale.py $inputDEMfile 1


import numpy as np
import sys
import os
import math
import gc

from osgeo import gdal
from osgeo.gdalconst import *



# register all of the GDAL drivers
gdal.AllRegister()



# The function 'assignBCs' is to assign the boundary values . This 
# function takes a numpy array and returns a numpy array that has an 
# additional row and column before and after those specified in the input 
# array. This function repeats the values on the edges of the array. 

def assignBCs(elevGrid):
    # Pads the boundaries of a grid
    # Boundary condition pads the boundaries with equivalent values 
    # to the data margins, e.g. x[-1,1] = x[1,1]
    # This creates a grid 2 rows and 2 columns larger than the input

    ny, nx = elevGrid.shape  # Size of array
    Zbc = np.zeros((ny + 2, nx + 2))  # Create boundary condition array
    Zbc[1:-1,1:-1] = elevGrid  # Insert old grid in center

    #Assign boundary conditions - sides
    Zbc[0, 1:-1] = elevGrid[0, :]
    Zbc[-1, 1:-1] = elevGrid[-1, :]
    Zbc[1:-1, 0] = elevGrid[:, 0]
    Zbc[1:-1, -1] = elevGrid[:,-1]

    #Assign boundary conditions - corners
    Zbc[0, 0] = elevGrid[0, 0]
    Zbc[0, -1] = elevGrid[0, -1]
    Zbc[-1, 0] = elevGrid[-1, 0]
    Zbc[-1, -1] = elevGrid[-1, 0]
    return Zbc


def calcFiniteSlope(elevGrid, weres, nsres, yoff, scale):
    # Sx,Sy = calcFiniteDiffs(elevGrid)
    # Sx, Sy - slopes in x and y directions
    # calculates finite differences in X and Y direction using the 
    # 2nd order/centered difference method.
    # Applies a boundary condition such that the size and location 
    # of the grids in is the same as that out. 

    # Assign boundary conditions
    # Compute finite differences and slope
    # We consider the change in dx in the latitudes. 
    # dx remains  the same for longitudes
    # for latitudes dx*cos(latitude)
    # latitude of any row = rotation * col + nsres *  col + yoff  
    # xoff, ewres, b, yoff, d, nsres = dem.GetGeoTransformw()
    # we calculate the latitude for each pixel using the formula
    # latitude = d* col + nsres * row + yoff

    Zbc = assignBCs(elevGrid)
    ny, nx = elevGrid.shape  # Size of array
    
    SxAbs = (Zbc[1:-1, 2:] - Zbc[1:-1, :-2])/2.0 
    SyAbs = (Zbc[2:,1:-1] - Zbc[:-2, 1:-1])/2.0
    
    if (scale == 1):
        dx = weres
        dy = nsres * -1.0
        Sx = SxAbs/dx
        Sy = SyAbs/dy

    else:            # scale 111120
        dxGrid = np.zeros((ny + 2, nx + 2))
        dx = weres * 111111.11
        dy = nsres * 111111.11
        for row in range(1,rows+1):
	    for col in range(1,cols+1):
	        dxGrid[row, col] = 2.0 * dx * math.cos(( nsres * row + yoff) * math.pi/180)
        Sx = (Zbc[1:-1, 2:] - Zbc[1:-1, :-2]) /dxGrid[1:-1,1:-1]
        Sy = (Zbc[2:,1:-1] - Zbc[:-2, 1:-1]) / (2.0 * dy)

    return SxAbs,Sx,SyAbs,Sy


          

# main program here

input_file = sys.argv[1]
scale = int(sys.argv[2])

path, filename = os.path.split(input_file)

slope_nsAbs =  path+"/slope/slope_NS_ABS_" + filename
slope_weAbs = path+"/slope/slope_WE_ABS_" + filename
slope_ns =  path+"/slope/slope_NS_PRCT_" + filename
slope_we = path+"/slope/slope_WE_PRCT_" + filename

slopeZAbs = path+"/slope/slope_Z_ABS_" + filename
slopeZ = path+"/slope/slope_Z_PRCT_" + filename
slopeZdeg = path+"/slope/slope_Z_DEG_" + filename

# open the dem
dem = gdal.Open(input_file)

# get info about it
band1 = dem.GetRasterBand(1)
rows = dem.RasterYSize
cols = dem.RasterXSize



                                        

# GDAL affine transform parameters, According to gdal documentation 
# xoff/yoff are image left corner, a/e are pixel wight/height and b/d 
# is rotation and is zero if image is north up.
#    adfGeoTransform[0]  /* top left x (origin X) xoff */
#    adfGeoTransform[1]  /* w-e pixel resolution (pixelWidth) */
#    adfGeoTransform[2]  /* rotation, 0 if image is north up */ 
#    adfGeoTransform[3]  /* top left y  (origin Y) yoff */
#    adfGeoTransform[4]  /* rotation, 0 if image is north up */
#    adfGeoTransform[5]  /* n-s pixel resolution  (pixelHeight) */

#    longitude of any column = pixelWidth * col + rotation * row + xoff
#    latitude of any row = rotation * col + pixelHeight * row + yoff  


xoff, weres, b, yoff, d, nsres = dem.GetGeoTransform()
elevGrid = band1.ReadAsArray(0,0,cols,rows)

for row in range(0,rows):
   for col in range(0,cols):
      if elevGrid[row,col]<-99999:
          elevGrid[row, col]=0
                                          

SxAbs,Sx, SyAbs, Sy = calcFiniteSlope(elevGrid, weres, nsres, yoff, scale)

SmagAbs =np.sqrt(SxAbs**2 + SyAbs**2)

Smag = 100* np.sqrt(Sx**2 + Sy**2)
SmagDeg = np.arctan(np.sqrt(Sx**2 + Sy**2)) * (180/math.pi)

# write the Abs slope in x direction  

# create the output image
driver = dem.GetDriver()
outDs = driver.Create(slope_weAbs, cols, rows, 1, GDT_Float32)
outBand = outDs.GetRasterBand(1)
outData = np.zeros((rows,cols), np.float32)

outData = SxAbs

# write the data
outBand.WriteArray(outData, 0, 0)
outBand.SetNoDataValue(-99999)

# georeference the image and set the projection
outDs.SetGeoTransform(dem.GetGeoTransform())
outDs.SetProjection(dem.GetProjection())

# del outData
outBand = None
outDs = None

# write the percentage slope in x direction 
# create the output image
driver = dem.GetDriver()
outDs = driver.Create(slope_we, cols, rows, 1, GDT_Float32)
outBand = outDs.GetRasterBand(1)
outData = np.zeros((rows,cols), np.float32)

outData = Sx

# write the data
outBand.WriteArray(outData, 0, 0)
outBand.SetNoDataValue(-99999)

# georeference the image and set the projection
outDs.SetGeoTransform(dem.GetGeoTransform())
outDs.SetProjection(dem.GetProjection())

# del outData
outBand = None
outDs = None



# write the Abs slope in y direction 

# create the output image
driver = dem.GetDriver()
outDs = driver.Create(slope_nsAbs, cols, rows, 1, GDT_Float32)
outBand = outDs.GetRasterBand(1)
outData = np.zeros((rows,cols), np.float32)

outData = SyAbs

# write the data
outBand.WriteArray(outData, 0, 0)
outBand.SetNoDataValue(-99999)

# georeference the image and set the projection
outDs.SetGeoTransform(dem.GetGeoTransform())
outDs.SetProjection(dem.GetProjection())

# del outData
outBand = None
outDs = None


# write the percentage slope in y direction 

# create the output image
driver = dem.GetDriver()
outDs = driver.Create(slope_ns, cols, rows, 1, GDT_Float32)
outBand = outDs.GetRasterBand(1)
outData = np.zeros((rows,cols), np.float32)

outData = Sy

# write the data
outBand.WriteArray(outData, 0, 0)
outBand.SetNoDataValue(-99999)

# georeference the image and set the projection
outDs.SetGeoTransform(dem.GetGeoTransform())
outDs.SetProjection(dem.GetProjection())

# del outData
outBand = None
outDs = None


# write the Abs slope in the z direction
# create the output image
driver = dem.GetDriver()
outDs = driver.Create(slopeZAbs, cols, rows, 1, GDT_Float32)
outBand = outDs.GetRasterBand(1)
outData = np.zeros((rows,cols), np.float32)

outData = SmagAbs

# write the data
outBand.WriteArray(outData, 0, 0)
outBand.SetNoDataValue(-99999)

# georeference the image and set the projection
outDs.SetGeoTransform(dem.GetGeoTransform())
outDs.SetProjection(dem.GetProjection())

# del outData
outBand = None
outDs = None

# write the percentage slope in the z direction
# create the output image
driver = dem.GetDriver()
outDs = driver.Create(slopeZ, cols, rows, 1, GDT_Float32)
outBand = outDs.GetRasterBand(1)
outData = np.zeros((rows,cols), np.float32)

outData = Smag

# write the data
outBand.WriteArray(outData, 0, 0)
outBand.SetNoDataValue(-99999)

# georeference the image and set the projection
outDs.SetGeoTransform(dem.GetGeoTransform())
outDs.SetProjection(dem.GetProjection())

# del outData
outBand = None
outDs = None

# write the radiance slope in the z direction
# create the output image
driver = dem.GetDriver()
outDs = driver.Create(slopeZdeg, cols, rows, 1, GDT_Float32)
outBand = outDs.GetRasterBand(1)
outData = np.zeros((rows,cols), np.float32)

outData = SmagDeg

# write the data
outBand.WriteArray(outData, 0, 0)
outBand.SetNoDataValue(-99999)

# georeference the image and set the projection
outDs.SetGeoTransform(dem.GetGeoTransform())
outDs.SetProjection(dem.GetProjection())

# del outData
outBand = None
outDs = None
dem = None
