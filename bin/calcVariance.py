#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This program calculates the variance of a DEM

# Usage: python2 $PATH/calcVariance.py $inputDEMfile

import cv2
import numpy as np
import sys
import os
import math

from osgeo import gdal
from osgeo.gdalconst import *



# register all of the GDAL drivers
gdal.AllRegister()




def sliding_window(a, window, axis=-1):
    shape = list(a.shape) + [window]
    shape[axis] -= window - 1
    if shape[axis] < 0:
        raise ValueError("Array too small")
    strides = a.strides + (a.strides[axis],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def sliding_img_var(img, window):
    if window <= 0:
        raise ValueError("invalid window size")
    buf = sliding_window(img, 2*window, 0)
    buf = sliding_window(buf, 2*window, 1)
    out = np.zeros(img.shape, dtype=np.float32)
    np.var(buf[:-1,:-1], axis=(-1,-2), out=out[window:-window,window:-window])
    return out
        

def looping_img_var(img, w):
    nx, ny = img.shape
    varianceMatrix = np.zeros(img.shape, np.float32)
    for i in range(w,nx-w):
        for j in range(w,ny-w):
            sampleframe = img[i-w:i+w, j-w:j+w]
            variance    = np.var(sampleframe)
            varianceMatrix[i][j] = variance
    return varianceMatrix


# main program here

input_file = sys.argv[1]

path, filename = os.path.split(input_file)

variance =  path+"/variance/variance_" + filename
#variance_sw =  path+"/variance/variance_sw_" + filename
#variance_loop =  path+"/variance/variance_loop_" + filename

# open the dem
dem = gdal.Open(input_file)

# get info about it
band1 = dem.GetRasterBand(1)
rows = dem.RasterYSize
cols = dem.RasterXSize

elevGrid = band1.ReadAsArray(0,0,cols,rows)

for row in range(0,rows):
   for col in range(0,cols):
      if elevGrid[row,col]<-99999:
          elevGrid[row, col]=0


# Calculating the variance using sliding window and loop method
# Loop method takes more time compared to sliding window
# https://stackoverflow.com/questions/28931265/calculating-variance-of-an-image-python-efficiently

var_sw = sliding_img_var(elevGrid,1)
#var_loop = looping_img_var(elevGrid, 1)


# write the variance  

# create the output image - variance by sliding window
driver = dem.GetDriver()
outDs = driver.Create(variance, cols, rows, 1, GDT_Float32)
outBand = outDs.GetRasterBand(1)
outData = np.zeros((rows,cols), np.float32)

outData = var_sw

# write the data
outBand.WriteArray(outData, 0, 0)
outBand.SetNoDataValue(-99999)

# georeference the image and set the projection
outDs.SetGeoTransform(dem.GetGeoTransform())
outDs.SetProjection(dem.GetProjection())


# delete outData
outBand = None
outDs = None
dem = None
