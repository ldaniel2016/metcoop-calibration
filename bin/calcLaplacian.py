#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This program calculates the laplacian in X and Y directions
#laplacian is the second partial derivative of a function and
# it is calculated using finite differences here
##       1    / d^2            d^2         \
## Laplacian  = --- * | ---  M(x,y) +  ---  M(x,y) | 
##       4    \ dx^2           dx^2        /
# In 1-dimension laplacian is f(x+1,y) -2f(x,y) +f(x-1,y)
# In 2D a 3x3 kernel is [0 -1 0, -1 4 -1, 0 -1 0]
# adapted from  https://stackoverflow.com/questions/4692196/
# discrete-laplacian-del2-equivalent-in-python
# and https://searchcode.com/codesearch/view/9565671/


# Usage: python2 $PATH/calcLaplacian.py $inputDEMfile

import numpy as np
import sys
import os
import math

from osgeo import gdal
from osgeo.gdalconst import *



# register all of the GDAL drivers
gdal.AllRegister()

# del2(M) calculates the laplacian using finite difference  method
# on the matrix/image array using a 3x3 grid

def del2(M, dx, dy):

    # Here M is the elevGrid, dy is the distance in y direction
    # As dy is a constance we pass it as a float and generate the actual dy
    # matrix using np.kron function with np.ones
    # dx is a row matrix  M and has the  distance in x direction
    # multiplied by cosine of latitude of each pixel. 
    
    mr, mc = M.shape  # rows and columns
    
    dy = dy * np.ones ((mr-1, 1))

   
    # Laplacian is the D array of mrxmc initialized to zero
    D = np.zeros((mr,mc))
    
    if (mr >= 3):
        ## x direction
        ## left and right boundary
        D[:,0] = (M[:, 0] - 2 * M[:, 1] + M[:, 2])/(dx[:,0] * dx[:,1])
        D[:,mc-1] = (M[:, mc - 3] - 2 * M[:, mc - 2] + M[:, mc - 1]) \
                    / (dx[:, mc - 3] * dx[:, mc - 2])

        ## interior points
        ## python array starts from index 0 and
        ## M[:, 1:mc] is all rows from second column to the last but one
        ## Kronecker product of two arrays. Computes the Kronecker product,
        ## a composite array made of blocks of the second array scaled
        ## by the first.
        
        tmp1 = D[:, 1:mc - 1] 
        tmp2 = (M[:, 2:mc] - 2 * M[:, 1:mc - 1] + M[:, 0:mc - 2])
        tmp3 = dx[:,0:mc -2] * dx[:,1:mc - 1]
        D[:, 1:mc - 1] = tmp1 + tmp2 / tmp3



    if (mc >= 3):
        ## y direction
        ## top and bottom boundary
        ### CHK (dy[mr-3,:] * dx[:,mr-2]) OR (dy[mr-3,:] * dy[mr-2, :])
        
        D[0, :] = D[0,:] + \
                  (M[0, :] - 2 * M[1, :] + M[2, :] ) / (dy[0,:] * dy[1,:])
        
        
        
        D[mr-1, :] = D[mr-1, :] \
                     + (M[mr-3,:] - 2 * M[mr-2, :] + M[mr-1, :]) \
                     / (dy[mr-3,:] * dy[mr-2,:])

        ## interior points
        tmp1 = D[1:mr-1, :] 
        tmp2 = (M[2:mr, :] - 2 * M[1:mr - 1, :] + M[0:mr-2, :])
        tmp3 = np.kron (dy[0:mr-2,:] * dy[1:mr-1,:], np.ones ((1, mc)))
        D[1:mr-1, :] = tmp1 + tmp2 / tmp3

    return D / 4





# main program here

input_file = sys.argv[1]

path, filename = os.path.split(input_file)

laplacian =  path+"/laplacian/laplacian_" + filename

# Laplacian using scipy

# open the dem
# example  input_file = "/home/daniel/digital_elevation_maps/Germany_100m_DEM_da# ta/saksa_dem_2km.tif"

dem = gdal.Open(input_file)
band1 = dem.GetRasterBand(1)
band1 = dem.GetRasterBand(1)
rows = dem.RasterYSize
cols = dem.RasterXSize

elevGrid = band1.ReadAsArray(0,0,cols,rows)


for row in range(0,rows):
   for col in range(0,cols):
      if elevGrid[row,col]<-99999:
          elevGrid[row, col]=0

# This is for wgs84 coordinate
# find the dx and dy matrix, the difference in the x direction and y direction
# dx for a pixel is nsres*111111.11 * cos latitude
# dy is ewres*111111.11
# We consider the change in dx in the latitudes. 
# dy remains  the same for longitudes
# for latitudes dx*cos(latitude)
# latitude of any row = rotation * col + nsres *  col + yoff  
# xoff, weres, b, yoff, d, nsres = dem.GetGeoTransform()
# we calculate the latitude for each pixel using the formula
# latitude = d* col + nsres * row + yoff
# b/d is rotation and is zero if image is north up
# one degree at the centre of Earth is equivalent to 111111.11 mts
# sometimes this value is rounded as 111120 (gdal)


xoff, weres, b, yoff, d, nsres = dem.GetGeoTransform()


dx = np.zeros((rows, cols),dtype=np.float32)
distance_x = abs(nsres)
distance_y = weres
print(distance_x, distance_y)
for row in range(0,rows):
    for col in range(0,cols):
        dx[row, col] = distance_x
        #dx[row, col] = distance_x * math.cos((d* col + nsres * row + yoff) * math.pi/180)

                
D = del2(elevGrid, dx, distance_y)

# write the Laplacian in a file  

# create the output image
driver =  dem.GetDriver()
outDs =   driver.Create(laplacian, cols, rows, 1, GDT_Float32)
outBand = outDs.GetRasterBand(1)
outData = np.zeros((rows,cols), np.float32)

outData = D

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
