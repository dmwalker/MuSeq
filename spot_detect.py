from math import sqrt
from skimage import data, util
from skimage.feature import blob_log
from skimage.color import rgb2gray
import numpy as np
from tkinter import Tk
from tkinter.filedialog import askopenfilename, askdirectory
import math
import gc
import matplotlib.pyplot as plt

# This script will ask you to choose a directory file. Your spot detection figure will go there
# To get the CFP spots file, use thresh hold in photoshop or GIMP2 until only the foci spots are shown.
# Then you must load in your CFP, CFP spots, GFP, and finally GFP spots
# Don't be afraid to mess around with the values in blob_log to get better detection

gc.collect()

# This will give you an error, but makes the font better looking
font = {'family' : 'normal',
        'weight' : 'normal',
        'size' : 28}
plt.rc('font', **font)

def open_image():
    Tk().withdraw()
    return data.load(askopenfilename())

def bLength( old, new ):
    return math.sqrt(np.dot((old-new), (old-new)))


blobs={}
Tk().withdraw()
direc = askdirectory()


# Selecting files
print('Open the brightfield image....')
brightField = open_image()
print('Open the CFP image...')
cfp=open_image()
print('Open the CFP thresh hold worked image...')
cfp_thresh = open_image()
print('Open the GFP image...')
gfp = open_image()
print('Open the GFP thresh hold worked image...')
gfp_thresh = open_image()

# This will color your blobs on your final figure
blobs['CFP']=[cfp_thresh, cfp, 1, 'c']
blobs['GFP']=[gfp_thresh, gfp, 2, 'g']

fig, axes=plt.subplots(1, 3, figsize=(36,9), sharex=True, sharey=True)
ax = axes.ravel()
ax[0].set_title('Bright Field')
ax[0].imshow(brightField, cmap='gray')
ax[0].set_axis_off()

spotLoc={}
for key in blobs:
    image=blobs[key][0]
    r_image=blobs[key][1]
    idx = blobs[key][2]
    colour=blobs[key][3]
    image_gray = rgb2gray(image) # Makes sure the image is gray. It should already be in grayscale anyways.

    # Detect the blobs
    blobs_log = blob_log(image_gray, min_sigma=1, max_sigma=5, num_sigma=10, threshold=.2, overlap=0.4)
    
    
# Compute radii in the 3rd column.
    ax[idx].set_title(key)
    ax[idx].imshow(r_image)
    ax[idx].set_axis_off()
    
    blobs_log[:, 2] = blobs_log[:, 2] * sqrt(2)
    spotLoc[key]=blobs_log
    for blob in blobs_log:
        y,x,r=blob
        f=plt.Circle((x,y), r+0.2, color=colour, linewidth=0.5, fill=False, alpha=0.8) #estimates where the focus is on the BF
        ax[idx].add_patch(f)
        c=plt.Circle((x,y), r, color=colour, linewidth=2, fill=True, alpha=0.8) #estimates where the other focus is on the BF
        ax[0].add_patch(c)


plt.tight_layout()
plt.savefig('%s/%s' % (direc, 'full_spotDetect.png')) #This will save the file in the directory you selected earlier.
plt.close()
    
# Calculating distances between nearest GFP and CFP neighbors
distances=[]
if len(spotLoc['GFP']) > len(spotLoc['CFP']):
    k = 'GFP'
    z = 'CFP'
else:
    k = 'CFP'
    z = 'GFP'
for g in spotLoc[z]:
	d_close = 1920
	x,y,r = g
	for c in spotLoc[k]:
		u,v,w=c
		d = bLength(np.array([x,y]),np.array([u,v]))
		if d < d_close:
			d_close = d
			c_coord = [u,v]
	distances.append([[x,y], c_coord, d_close])

distances = np.array(distances); g,c,d = distances.transpose()
counts, bins = np.histogram(d, bins=1000, range=(0,156))
np.savetxt('%s/distances.txt' % direc, (bins[:-1], counts), delimiter=' ')
fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(d, bins=300, range=(0,25), color='k') # This will exclude anything greather than 750 uM apart.

# Estimated as distances in nm
ax.set_xticklabels([0, 150, 300, 450, 600, 750])

# This will show you the plot and will not save. You have to save manually.
plt.show()

