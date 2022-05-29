from skimage import io
import numpy as np
import matplotlib.pyplot as plt
import os
from os.path import join, dirname, abspath, exists
from PIL import Image
from PIL.ImageOps import expand
import seaborn as sns
import sys
import time
import argparse
from pandas import read_csv
import numpy as np
import pandas as pd
import json
#from math import *
import math
from numpy.ma import exp,sin,cos,abs
from pylab import *
import matplotlib.mlab as mlab  
from skimage import data,color,morphology,feature
from math import exp
from scipy.optimize import curve_fit

####################Function#######################
def gauss(x,mu,sigma,A):
    return A*exp(-(x-mu)**2/2/sigma**2)

def bimodal(x,mu1,sigma1,A1,mu2,sigma2,A2):
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)
def _solve(m1,std1,s1,m2,std2,s2):
    """solve equation: s1*N(m1,std1)=s2*N(m2,std2), return the intersection points of two weighted Gaussian"""
    a = 1/(2*std1**2) - 1/(2*std2**2)
    b = m2/(std2**2) - m1/(std1**2)
    c = m1**2 /(2*std1**2) - m2**2 / (2*std2**2) - np.log((std2*s1)/(std1*s2))
    return np.roots([a,b,c])

def image_init(image_path,barcode_path,
               spot_path,
               save_path,
               tissue_hires_scalef = 1,
               dx = 20,dy = 20,
               rgb_score = -1):
    img=io.imread(image_path)
    io.imshow(img)
    #img_gray=io.imread(image_path,as_grey=True)
    barcode = []
    tissue = []
    row = []
    col = []
    imagerow = []
    imagecol = []
    hires_scale_row = []
    hires_scale_col = []
    with open(spot_path,'r') as file:
        lines = file.readlines()
        for i in range(1,len(lines)):
            tm = lines[i].strip('\n').split(',')
            barcode = barcode+[tm[0]]
            tissue = tissue +[tm[1]]
            row = row + [int(tm[2])]
            col = col + [int(tm[3])]
            imagerow = imagerow + [int(float(tm[4]))]
            imagecol = imagecol + [int(float(tm[5]))]
            hires_scale_row = hires_scale_row +[int(int(tm[4])*tissue_hires_scalef)]
            hires_scale_col = hires_scale_col +[int(int(tm[5])*tissue_hires_scalef)]
    
    spot_loc= pd.DataFrame({'barcode':barcode, 'tissue':tissue, 'row':row,'col':col,
                            'imagerow' : imagerow ,'imagecol' : imagecol,
                            'hires_scale_row':hires_scale_row,'hires_scale_col':hires_scale_col})
    barcode_filter = []
    with open(barcode_path,'r') as file:
        lines = file.readlines()
        for i in range(1,len(lines)):
            tm = lines[i].strip('\n').split(',')
            barcode_filter = barcode_filter + [tm[0]]
    
    spot_loc_filter = spot_loc.loc[spot_loc['barcode'].isin(list(barcode_filter))]
    #filter spots
    searchID = []
    for j in range(0,int(dx/2+5)):
        for k in range(0,int(dy/2+5)):
            searchID = searchID + [(j,k),(-j,k),(j,-k),(-j,-k)]
    if(rgb_score==-1):
        show_r = img[:,:,0].reshape(-1,1)
        n, bins, patches = plt.hist(x = show_r, bins = 500, facecolor='blue', alpha=0.5)  
        expected=(np.percentile(show_r,35),0.99,1,np.percentile(show_r,90),0.1,1)
        params,cov=curve_fit(bimodal,bins[0:500],n,expected)
        result = _solve(*abs(params))
        rgb_score = min(result)
    mask = np.zeros(shape = (img.shape[0],img.shape[1]))
    for i in range(len(spot_loc_filter)):
        y = spot_loc_filter.hires_scale_col[spot_loc_filter.index[i]]
        x = spot_loc_filter.hires_scale_row[spot_loc_filter.index[i]]
        mask[x,y] = 1
        for (j,k) in searchID:
            tx = x + j
            ty = y + k
            if img[tx,ty,0]<=rgb_score:
                mask[tx,ty] = 1
    mask = mask.astype(bool)
    mask1=morphology.remove_small_objects(mask,min_size=180,connectivity=1)
    mask1 = 1-mask1
    mask1 = mask1.astype(bool)
    ##filter use spot info
    mask1=morphology.remove_small_objects(mask1,min_size=180,connectivity=1)
    mask1 = 1-mask1
    mask1 = mask1.astype(bool)
    #io.imshow(mask1)
    np.save(save_path+"mask1.npy",mask1)
    np.savetxt(save_path+"mask1.txt", mask1,fmt="%d", delimiter=",") 
    # draw red line on org img
    fig, ax_arr = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(10, 10))
    ax1, ax2= ax_arr.ravel()
    
    ax1.imshow(mask1)
    #m_slic = segmentation.slic(img, n_segments=100, mask=mask, start_label=1)
    ax2.imshow(img)
    ax2.contour(mask1, colors='red', linewidths=1)
    
    for ax in ax_arr.ravel():
        ax.set_axis_off()
    
    plt.tight_layout()
    plt.savefig(save_path+"mask1.png",dpi=800)
    plt.show()
    return mask1

######################Example######################
mask_ = image_init(image_path = "C:/Gu_lab/space_expr/data/225847_C4/spatial/tissue_hires_image.png",
           barcode_path = "C:/Gu_lab/space_expr/data/225847_C4/spatial/Imginit/barcodes.tsv",
           spot_path = "C:/Gu_lab/space_expr/data/225847_C4/spatial/tissue_positions_list.csv",
           save_path = "C:/Gu_lab/space_expr/data/225847_C4/spatial/Imginit/",
           tissue_hires_scalef = 0.09735203,
           dx = 260*0.09735203,dy = 260*0.09735203,rgb_score = -1)
