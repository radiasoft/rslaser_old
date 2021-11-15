# -*- coding: utf-8 -*-
u"""Methods for plotting electromagnetic fields.
Copyright (c) 2021 RadiaSoft LLC. All rights reserved
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import rslaser.rspulse.gauss_hermite as rsgh
import rslaser.utils.plot_tools as rspt

def plot_zx(_zArr, _xArr, _pulse, _ax):
    numX = np.size(_xArr)
    minX = np.min(_xArr)
    maxX = np.max(_xArr)
    delX = (maxX-minX) / (numX-1)
        
    numZ = np.size(_zArr)
    minZ = np.min(_zArr)
    maxZ = np.max(_zArr)
    delZ = (maxZ-minZ) / (numZ-1)

    for iLoop in range(numZ):
        _zArr[iLoop] = minZ + iLoop * delZ

    for jLoop in range(numX):
        _xArr[jLoop] = minX + jLoop * delX

    # specify y position for plot
    yValue = 0.

    # Calculate Ex at the 2D array of x,y values
    zxEData = np.zeros((numX, numZ))
    for iLoop in range(numZ):
        for jLoop in range(numX):
            contour_value = np.real(_pulse.evaluate_envelope_ex(_xArr[jLoop], yValue, _zArr[iLoop]))
            zxEData[jLoop, iLoop] = contour_value

    # generate the contour plot
    n_levels = 10
    levels = rspt.generate_contour_levels(zxEData, n_levels)
    _ax.axis([minZ, maxZ, minX, maxX])
    _ax.set_xlabel('z [m]')
    _ax.set_ylabel('x [m]')
    _ax.set_title('ZX slice, at  y={0:4.2f} [{1}]'.format(yValue, '[m]'))

    rspt.scatter_contour('contour', 'linear', _xArr, _zArr, _ax, 10, n_levels)
    
    
def plot_zy(_pulse, _ax):

    # Specify the desired grid size
    zyNumH = 256
    zyNumV = 256
    zyNumCells = zyNumH * zyNumV

    # specify the min's and max's
    zyMinH = -20. * _pulse.lambda0
    zyMaxH =  20. * _pulse.lambda0

    zyMinV = -2. * _pulse.waist_y
    zyMaxV =  2. * _pulse.waist_y

    zArr = np.zeros(zyNumH)
    yArr = np.zeros(zyNumV)

    for iLoop in range(zyNumH):
        zArr[iLoop] = zyMinH + iLoop * (zyMaxH-zyMinH) / (zyNumH-1)

    for jLoop in range(zyNumV):
        yArr[jLoop] = zyMinV + jLoop * (zyMaxV-zyMinV) / (zyNumV-1)

    # Choose values of x,t for plot
    xValue = 0.
    tValue = 0.

    # Calculate Ex at the 2D array of x,y values
    zyEData = np.zeros((zyNumV, zyNumH))
    for iLoop in range(zyNumH):
        for jLoop in range(zyNumV):
            zyEData[jLoop, iLoop] = _pulse.evaluate_ex(xValue, yArr[jLoop], zArr[iLoop], tValue)

    # generate the contour plot
    _ax.clear()

    n_levels = 10
    levels = rspt.generate_contour_levels(zyEData, n_levels)
    _ax.contourf(zArr, yArr, zyEData, levels, extent='none')

    _ax.axis([zyMinH, zyMaxH, zyMinV, zyMaxV])
    _ax.set_xlabel('z [m]')
    _ax.set_ylabel('y [m]')
    _ax.set_title('ZY slice, at  x={0:4.2f} [m]')

    rspt.scatter_contour('contour', 'linear', yArr, zArr, _ax, 10, n_levels)

    
def plot_xy(_z_waist, _pulse, _ax):

    # Specify the desired grid size
    xyNumH = 64
    xyNumV = 64
    xyNumCells = xyNumH * xyNumV

    # specify the min's and max's
    xyMinH = -2. * _pulse.waist_x
    xyMaxH =  2. * _pulse.waist_x

    xyMinV = -2. * _pulse.waist_y
    xyMaxV =  2. * _pulse.waist_y

    xArr = np.zeros(xyNumH)
    yArr = np.zeros(xyNumV)
    xTmp = np.zeros(xyNumV)

    for iLoop in range(xyNumH):
        xArr[iLoop] = xyMinH + iLoop * (xyMaxH-xyMinH) / (xyNumH-1)

    for jLoop in range(xyNumV):
        yArr[jLoop] = xyMinV + jLoop * (xyMaxV-xyMinV) / (xyNumV-1)

    # Calculate Ex at the 2D array of x,y values
    xyEData = np.zeros((xyNumV, xyNumH))
    for iLoop in range(xyNumH):
        for jLoop in range(xyNumV):
            xTmp[jLoop] = xArr[iLoop]
        xyEData[0:xyNumV, iLoop] = _pulse.evaluate_envelope_ex(xTmp, yArr, _z_waist)

    # generate the contour plot
    _ax.clear()

    n_levels = 20
    levels = rspt.generate_contour_levels(xyEData, n_levels)
    _ax.contourf(xArr, yArr, xyEData, levels, extent='none')

    _ax.axis([xyMinH, xyMaxH, xyMinV, xyMaxV])
    _ax.set_xlabel('x [m]')
    _ax.set_ylabel('y [m]')
    _ax.set_title('XY slice, at waist location')

    rspt.scatter_contour('contour', 'linear', xArr, yArr, _ax, 10, n_levels)
