# -*- coding: utf-8 -*-
u"""Merit functions used for nonlinear optimization
Copyright (c) 2022 RadiaSoft LLC. All rights reserved
"""
import numpy as np

# (x,y) are 2D numpy arrays, specifying the Cartesian domain of the function
# (xc, yc) is the center of the laser pulse
# (phi_min, phi_max) are the min/max values of the wavefront phase
# r_max is approximately the maximum radius of the image
# r_mid ~ 0.2*r_max is user-specified, with a value that might need tweaking
# r0 is the radial curvature of phi(r), for r < r_mid
def spline_wfs(x, y, xc, yc, r_mid, r_max, phi_min, phi_max, r0):
    rsq = (x-xc)**2 + (y-yc)**2
    r = np.sqrt(rsq)
    phi_mid = phi_max * (1. - (r_mid/r0)**2)
    c0 = (phi_mid - phi_min) / (1. - r_mid/r_max)
    
    f_1 = phi_max * (1. - rsq/r0**2)
    f_2 = phi_min + c0 * (1. - r/r_max)
    return np.where(r<=r_mid, f_1, f_2)

# r0 and phi_max are the fitting parameters
def spline_wfs_fit(params, x, y, data, xc, yc, r_mid, r_max, phi_min):
    r0 = params[0]
    phi_max = params[1]
    g = spline_wfs(x, y, xc, yc, r_mid, r_max, phi_min, phi_max, r0)
    r = np.sqrt((x-xc)**2 + (y-yc)**2)
    weighted_diff = np.where(r<=r_mid, (g - data.reshape(g.shape))/r, 0.)
    return np.linalg.norm(weighted_diff)

# (x,y) are 2D numpy arrays, specifying the Cartesian domain of the function
# (xc, yc) is the center of the laser pulse
# n0_max is the maximum # of photon counts on a pixel, near the pulse center
# r_rms, the RMS radius of the 2D cylindrically symmetric gaussian
def gaussian_ccd(x, y, xc, yc, n0_max, r_rms):
    return n0_max * np.exp(-(((x-xc)/r_rms)**2 + ((y-yc)/r_rms)**2))

# r_rms and n0_max are the fitting parameters
def gaussian_ccd_fit(params, x, y, xc, yc, data, r_mid):
    r_rms = params[0]
    n0_max = params[1]
    g = gaussian_ccd(x, y, xc, yc, n0_max, r_rms)
    rsq = (x-xc)**2 + (y-yc)**2
    r = np.sqrt(rsq)
    weighted_diff = np.where(r<r_mid, (g - data.reshape(g.shape))/r, 0.)
    return np.linalg.norm(weighted_diff)

# compute an azimuthally averaged radial profile
# Code adapted from https://github.com/vicbonj/radialprofile/blob/master/radialProfile.py

def azimuthalAverage(image, centerx, centery, type='mean'):
    '''
    Compute spherically symetric profiles around a center
    Returns
    -------
    profiles, errors, distance to the center in pixels
    '''

    y, x = np.indices(image.shape)
    r = np.hypot(x - centerx, y - centery)

    ind = np.argsort(r.flat)
    r_sorted = r.flat[ind]
    i_sorted = image.flat[ind]

    r_int = r_sorted.astype(int)

    deltar = r_int[1:] - r_int[:-1]
    rind = np.where(deltar)[0]
    rind2 = rind+1
    rind3 = np.zeros(len(rind2)+1)
    rind3[1:] = rind2
    rind3 = rind3.astype('int')

    if type == 'mean':
        aaa = [np.nanmean(i_sorted[rind3[i]:rind3[i+1]]) for i in range(len(rind3)-1)]
    elif type == 'median':
        aaa = [np.nanmedian(i_sorted[rind3[i]:rind3[i+1]]) for i in range(len(rind3)-1)]
    elif type == 'mode':
        aaa_list = [i_sorted[rind3[i]:rind3[i+1]] for i in range(len(rind3)-1)]
        aaa = []
        for part in aaa_list:
            if len(part) == 1:
                counts, xed = np.histogram(part, bins=len(part))
            elif (len(part) > 1) & (len(part) < 40):
                counts, xed = np.histogram(part, bins=int(len(part)/2))
            else:
                counts, xed = np.histogram(part, bins=20)
            if len(xed) != 2:
                aaa.append(0.5*(xed[1:]+xed[:-1])[np.where(counts == np.max(counts))[0]][0])
            elif len(xed) == 2:
                aaa.append(0.5*(xed[1:]+xed[:-1])[np.where(counts == np.max(counts))[0]][0])
    else:
        raise ValueError('Nope')
    aaa_std = [np.nanstd(i_sorted[rind3[i]:rind3[i+1]]) for i in range(len(rind3)-1)]
    dist_r = [np.mean(r_sorted[rind3[i]:rind3[i+1]]) for i in range(len(rind3)-1)]
    return np.array(aaa), np.array(aaa_std), np.array(dist_r)
