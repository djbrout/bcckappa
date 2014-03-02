import matplotlib.pyplot as plt
import matplotlib.axes as axes
import os
import sys
import numpy as np
from scipy.stats.stats import pearsonr
import pyfits as pf
import scipy.signal as signal
import scipy.ndimage.filters as filter
import rdcol2
import subprocess
import dillon_v2 as catalogue

#Helper function to convolve with Gauss and write to new fits file
def convolve_maps(smoothing_grid,some_fit,some_map,filename,isFg):
    for sigma in smoothing_grid:
        print sigma
        if sigma > 0:
            if isFg:
                SmoothKappa = filter.gaussian_filter(some_fit[0].data,sigma)
                SmoothMask = filter.gaussian_filter(np.array(some_fit[1].data).astype("float"),sigma)
            else:
                SmoothKappa = filter.gaussian_filter(some_fit[0].data,sigma)
                SmoothMask = filter.gaussian_filter(np.array(some_fit[2].data).astype("float"),sigma)    

            Mask_ma = trim_mask(SmoothMask)
            SmoothKappa_masked = np.ma.masked_array(SmoothKappa, mask=Mask_ma)
            save_fits_image(SmoothKappa,filename+"_smooth"+str(sigma)+".fits")
            save_fits_image(Mask_ma,"Mask_"+filename+"_smooth"+str(sigma)+".fits") 

        else:
            save_fits_image(some_map,filename+"_smooth"+str(sigma)+".fits")
            save_fits_image(some_fit[-1].data,"Mask_"+filename+"_smooth"+str(sigma)+".fits")

def save_fits_image(image,filename):
    hdu = pf.PrimaryHDU(image)
    if os.path.exists(filename):
        os.remove(filename)
    hdu.writeto(filename)

def trim_mask(SmoothMask):
    Mask = SmoothMask
    Mask[(SmoothMask < .9)] = 0
    Mask[(SmoothMask >= .9)] = 1
    return Mask

def pcc(smoothing_grid,current_magcut,filename,f,config_smoothing,pixel_scale):
    for sigma in smoothing_grid:
        fg_map_p = pf.open("fgmap_fg"+filename+"_smooth"+str(sigma)+".fits")
        k_map_p = pf.open( "kappamap"+filename+"_smooth"+str(sigma)+".fits")
        Masks = pf.open( "Mask_kappamap"+filename+"_smooth"+str(sigma)+".fits")

        ff = np.ma.masked_array(fg_map_p[0].data.ravel(), mask=np.logical_not(Masks[0].data.ravel()))
        kk = np.ma.masked_array(k_map_p[0].data.ravel(), mask=np.logical_not(Masks[0].data.ravel()))

        corr = np.ma.corrcoef(ff,kk)
        print "Cross Corr np.ma:"
        print corr
        f.write(str(corr[0,1])+"\t"+str(current_magcut)+"\t"+str(np.sqrt(sigma**2+config_smoothing**2))+"\t"+str(pixel_scale)+"\n")

def make_smooth_plots():
    cols = rdcol2.read('pcc_table.txt',1,2)
    fig = plt.figure()
    errs = []
    for magc in np.unique(cols["Mag_Cut"]):
        plt.scatter(np.array(cols["Gauss_Smooth"])[np.array(cols["Mag_Cut"]) == magc],
                    np.array(cols["Pearson_Corr"])[np.array(cols["Mag_Cut"]) == magc],
                    label = "Mag Cut: "+str(magc))
        
    plt.xlabel("Smoothing Scale (arcmin)")
    plt.ylabel("$\kappa$ x Fg")
    plt.ylim(0,1)
    #handles, labels = axes.get_legend_handles_labels()
    #import operator
    #hl = sorted(zip(handles, labels),
    #            key=operator.itemgetter(1))
    #handles2, labels2 = zip(*hl)
    #plt.legend(handles2, labels2)
    fig.savefig("cc_vs_smoothing_and_magcut.png")

def make_pix_plots():
    cols = rdcol2.read('pcc_table.txt',1,2)
    fig = plt.figure()
    errs = []
    for magc in np.unique(cols["Mag_Cut"]):
        plt.scatter(np.array(cols["Pixel_Scale"])[np.array(cols["Mag_Cut"]) == magc],
                    np.array(cols["Pearson_Corr"])[np.array(cols["Mag_Cut"]) == magc],
                    label = "Mag Cut: "+str(magc))

    plt.xlabel("Pixel Scale (arcmin)")
    plt.ylabel("$\kappa$ x Fg")
    plt.ylim(0,1)
    #handles, labels = axes.get_legend_handles_labels()
    #import operator
    #hl = sorted(zip(handles, labels),
    #            key=operator.itemgetter(1))
    #handles2, labels2 = zip(*hl)
    #plt.legend(handles2, labels2)
    fig.savefig("cc_vs_pixel_size.png")
                                                                

#get effective smoothing by adding in quadtrature
def get_smoothing(eff_sm,config_sm):
    return_sm = eff_sm
    index = -1
    for eff in eff_sm:
        print eff
        print config_sm
        print np.sqrt(eff**2 - config_sm**2)
        index = index+1
        return_sm[index] = np.sqrt(eff**2 - config_sm**2)
    return return_sm

def edit_config(pixel_scale,smoothing):
    f = open("config.py","r")
    lines = f.readlines()
    lines[0] = "pixel_scale = "+str(pixel_scale)+" #size of the pixel in arc min \n"
    lines[7] = "smooth_size = "+str(smoothing)+" #smoothing scale for background arcmin \n"
    lines[8] = "smooth_size_n = "+str(smoothing)+" #smoothing scale for background arcmin \n"
    g = open("config.py","w")
    for line in lines:
        g.write(line)
    g.close()


def run_ks_mapping():
    a = os.popen("/home/vinu/software/Python2.7/bin/python example.py")
    
def analyze_ks_output(filename,eff_smoothing_grid,mag_cut,config_smoothing,pix,f):
    k_fit = pf.open("kappamap"+filename+".fits")
    k_map = k_fit[0].data
    fg_fit = pf.open("fgmap_fg"+filename+".fits")
    fg_map = fg_fit[0].data
    
    
    smoothing_grid = get_smoothing(eff_smoothing_grid,config_smoothing)
    
    convolve_maps(smoothing_grid,k_fit,k_map,"kappamap"+filename,False)
    convolve_maps(smoothing_grid,fg_fit,fg_map,"fgmap_fg"+filename,True)
    
    pcc(smoothing_grid,mag_cut,filename,f,config_smoothing,pix)






zphot = False
shape_noise = False

if zphot:
    import catalogue_zphot as catalogue
    f = open('pcc_zphot_table.txt', 'w')    
else:
    if shape_noise:
        import catalogue_epsilon as catalogue
        f = open('pcc_epsilon_table.txt','w')
    else:
        import dillon_v2 as catalogue
        f = open('pcc_table.txt', 'w')
    
f.write('Pearson_Corr\t Mag_Cut\t Gauss_Smooth\t Pixel_Scale\n')


#filename = "_im3shape_r_2.0_2.0_g1"
eff_smoothing_grid = [0.5,0.75,1,1.5,2,4,8,15,20,30,40,50,60]
#eff_smoothing_grid = [0.5]

config_smoothing = 0.5
current_mag_cut = 999
pix = 1.0
mag_grid = [23]
edit_config(pix,config_smoothing)
for mag in mag_grid:
    #catalogue.run_catalogue(mag)
    run_ks_mapping()
    filename = "_im3shape_r_"+str(pix)+"_"+str(config_smoothing)+"_g1"
    analyze_ks_output(filename,eff_smoothing_grid,mag,config_smoothing,pix,f)

f.close()
#make_smooth_plots()
#make_pix_plots()
