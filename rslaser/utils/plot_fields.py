# -*- coding: utf-8 -*-
u"""Methods for plotting electromagnetic fields.
Copyright (c) 2021 RadiaSoft LLC. All rights reserved
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import rslaser.rspulse.gauss_hermite as rsgh
import rslaser.utils.plot_tools as rspt

def plot_1d_x(_xArr, _pulse, _ax, _y=0., _z=0., _t=0., _time_explicit=False):
    """Plot a 1D (longitudinal, along z) lineout of a laser pulse field.

    For now, we are assuming the field is Ex

    Args:
        _xArr (1D numpy array): x positions where the field is to be evaluated
        _pulse (object): instance of the GaussHermite() class.
        _ax (object): matplotlib 'axis', used to generate the plot
        _y (float): [m] vertical location of the lineout
        _z (float): [m] longitudinal location of the lineout
        _t (float): [m] time of the lineout
        _time_explicit (bool): envelope only (False) or time explicit (True)
    """
    numX = np.size(_xArr)
#    print('\n _xArr = ', _xArr)

    # Calculate Ex at the 2D array of x,y values
    em_field = np.zeros(numX)
    if (_time_explicit):
#        for iLoop in range(numX):
#            em_field[iLoop] = _pulse.evaluate_ex(_xArr[iLoop], _y, _z, _t)
        em_field = _pulse.evaluate_ex(_xArr, _y, _z, _t)
    else:
#        for iLoop in range(numX):
#            em_field[iLoop] = np.real(_pulse.evaluate_envelope_ex(_xArr[iLoop], _y, _z))
        em_field = np.real(_pulse.evaluate_envelope_ex(_xArr, _y, _z))

    _ax.plot(_xArr, em_field)
    _ax.set_xlabel('x [m]')
    _ax.set_ylabel('Ex [V/m]')
    if (_time_explicit):
        _ax.set_title('Ex [V/m], at (y,z)=({0:4.2f},{0:4.2f}) [m] and t={0:4.2f} [s]'.format(_y, _z, _t))
    else:
        _ax.set_title('Ex (envelope) [V/m], at (y,z)=({0:4.2f},{0:4.2f}) [m]'.format(_y, _z))

        
def plot_1d_z(_zArr, _pulse, _ax, _x=0., _y=0., _t=0., _time_explicit=False):
    """Plot a 1D (longitudinal, along z) lineout of a laser pulse field.

    For now, we are assuming the field is Ex

    Args:
        _zArr (1D numpy array): z positions where the field is to be evaluated
        _pulse (object): instance of the GaussHermite() class.
        _ax (matplotlib axis): used to generate plot
        _x (float): [m] horizontal location of the lineout
        _y (float): [m] vertical location of the lineout
        _t (float): [m] time of the lineout
        _time_explicit (bool): envelope only (False) or time explicit (True)
    """
    numZ = np.size(_zArr)

    # Calculate Ex at the 2D array of x,y values
    em_field = np.zeros(numZ)
    if (_time_explicit):
        for iLoop in range(numZ):
            em_field[iLoop] = _pulse.evaluate_ex(_x, _y, _zArr[iLoop], _t)
    else:
        for iLoop in range(numZ):
            em_field[iLoop] = np.real(_pulse.evaluate_envelope_ex(_x, _y, _zArr[iLoop]))

    _ax.plot(_zArr, em_field)
    _ax.set_xlabel('z [m]')
    _ax.set_ylabel('Ex [V/m]')
    if (_time_explicit):
        _ax.set_title('Ex [V/m], at (x,y)=({0:4.2f},{0:4.2f}) [m] and t={0:4.2f} [s]'.format(_x, _y, _t))
    else:
        _ax.set_title('Ex (envelope) [V/m], at (x,y)=({0:4.2f},{0:4.2f}) [m]'.format(_x, _y))
    
    
def plot_zy(_pulse, _ax):
    
    # longitudinal resolution
    hres = 16
    dh = _pulse.lambda0 / hres

    # specify the horizontal min's and max's
    zyMinH = -16. * _pulse.lambda0
    zyMaxH =  16. * _pulse.lambda0

    # specify the number of horizontal mesh points
    zyNumH = int((zyMaxH - zyMinH) / dh)
    
    # vertical resolution
    vres = 64
    dv = _pulse.waist_y / vres

    # specify the number of vertical mesh points
    zyNumV = zyNumH

    # specify the vertical min's and max's
    zyMaxV =  zyNumV * dv / 2
    zyMinV = -zyMaxV

    # specify the number of cells in the 2D mesh
    zyNumCells = zyNumH * zyNumV
    
    print('nx, ny, n_cells = ', zyNumH, zyNumV, zyNumCells)

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
        zyEData[:, iLoop] = _pulse.evaluate_ex(xValue, yArr, zArr[iLoop], 0.)

    # generate the contour plot
    _ax.clear()

    n_levels = 10
    levels = rspt.generate_contour_levels(zyEData, n_levels)
    _ax.contourf(zArr, yArr, zyEData, levels, extent='none')

    _ax.axis([zyMinH, zyMaxH, zyMinV, zyMaxV])
    _ax.set_xlabel('z [m]')
    _ax.set_ylabel('y [m]')
    _ax.set_title('ZY slice, at  x={0:4.2f} [m]'.format(xValue))

    rspt.scatter_contour('contour', 'linear', yArr, zArr, _ax, 10, n_levels)


def plot_zx(_zArr, _xArr, _pulse, _ax):
    numX = np.size(_xArr)
    minX = np.min(_xArr)
    maxX = np.max(_xArr)
        
    numZ = np.size(_zArr)
    minZ = np.min(_zArr)
    maxZ = np.max(_zArr)

    # specify y position for plot
    yValue = 0.

    # Calculate Ex at the 2D array of x,y values
    zxEData = np.zeros((numX, numZ))
    for iLoop in range(numZ):
        zxEData[:, iLoop] = np.real(_pulse.evaluate_envelope_ex(_xArr, yValue, _zArr[iLoop]))

    # generate the contour plot
    n_levels = 10
    levels = rspt.generate_contour_levels(zxEData, n_levels)
    _ax.axis([minZ, maxZ, minX, maxX])
    _ax.set_xlabel('z [m]')
    _ax.set_ylabel('x [m]')
    _ax.set_title('ZX slice, at  y={0:4.2f} [{1}]'.format(yValue, '[m]'))

    rspt.scatter_contour('contour', 'linear', _xArr, _zArr, _ax, 10, n_levels)

    
def plot_xy(_z_waist, _pulse, _ax):

    # Specify the desired grid size
    xyNumH = 128
    xyNumV = 128

    # specify the min's and max's
    xyMinH = -2. * _pulse.waist_x
    xyMaxH =  2. * _pulse.waist_x

    xyMinV = -2. * _pulse.waist_y
    xyMaxV =  2. * _pulse.waist_y

    xArr = np.zeros(xyNumH)
    yArr = np.zeros(xyNumV)

    for iLoop in range(xyNumH):
        xArr[iLoop] = xyMinH + iLoop * (xyMaxH-xyMinH) / (xyNumH-1)

    for jLoop in range(xyNumV):
        yArr[jLoop] = xyMinV + jLoop * (xyMaxV-xyMinV) / (xyNumV-1)

    # Calculate Ex at the 2D array of x,y values
    xyEData = np.zeros((xyNumH, xyNumV))
    for iLoop in range(xyNumV):
        xyEData[:, iLoop] = np.real(_pulse.evaluate_envelope_ex(xArr, yArr, _z_waist))

    # generate the contour plot
    _ax.clear()
    _ax.axis([xyMinH, xyMaxH, xyMinV, xyMaxV])
    _ax.set_xlabel('x [m]')
    _ax.set_ylabel('y [m]')
    _ax.set_title('XY slice, at waist location')

    rspt.scatter_contour('contour', 'linear', xArr, yArr, _ax)
