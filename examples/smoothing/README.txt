***
gauss_smoother.py 

This script applies a simple gaussian filter to the provided images and produces plots with
various gaussian sigmas.

Execute on the command line as follows, or see the smooth_*.ipynb notebooks in this directory

> python gauss_smoother.py ../../rslaser/package_data/raw_beam_profile_532nm_pump.bmp


***
supergauss_fit.py 

This script fits a super-gaussian function to a provided image, and plots the result along with the
contours of the image. It has not been tested on a large variety of images.

Execute on the command line as follows, or see the supergauss_*.ipynb notebooks in this directory

> python supergauss_fit.py ../../rslaser/package_data/raw_beam_profile_800nm_seed.bmp
