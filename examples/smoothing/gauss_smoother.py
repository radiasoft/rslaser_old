from PIL import Image
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage.filters import gaussian_filter
import sys

blurs = [5, 10, 15, 20]
img = Image.open(sys.argv[1])
img = np.array(img)
f, axarr = plt.subplots(5, 1)
axarr[0].imshow(img, cmap="gray", vmin=0, vmax=255)
axarr[0].set_title("Original Image")
for i, blur in enumerate(blurs):
    axarr[i + 1].imshow(gaussian_filter(img, sigma=blur), cmap="gray", vmin=0, vmax=255)
    axarr[i + 1].set_title("Image with guassian blue of sigma " + str(blur))
f.tight_layout()
plt.show()
