import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import sys


def super_gaussian(params, amplitude, xo, yo, sigma, expon):
    expon = 0.0001 * expon
    xo = float(xo)
    yo = float(yo)
    r = np.sqrt((params[0] - xo)**2 + (params[1] - yo)**2)
    g = amplitude * np.exp(-(r/sigma)**expon)
    return g.ravel()


img = Image.open(sys.argv[1])
img = np.array(img)

x = np.linspace(0, img.shape[1] - 1, img.shape[1])
y = np.linspace(0, img.shape[0] - 1, img.shape[0])
x, y = np.meshgrid(x, y)

# plot twoD_Gaussian data generated above
plt.figure()
plt.imshow(img, origin='lower')
plt.colorbar()

initial_guess = (3, 300, 300, 10, 2)

popt, pcov = opt.curve_fit(super_gaussian, (x, y), img.flatten(), p0=initial_guess, maxfev=10000)

print(popt)

data_fitted = super_gaussian((x, y), *popt).reshape(img.shape)

fig, ax = plt.subplots(1, 1)
im = ax.imshow(data_fitted, origin='lower', extent=(x.min(), x.max(), y.min(), y.max()))
fig.colorbar(im)
ax.contour(x, y, data_fitted, origin='lower')
plt.show()
