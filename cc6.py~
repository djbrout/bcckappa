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
import config as c
#import dillon_v2 as catalogue

class cc:
    def __init__(self,config_smoothing = 0.5, pix = 1.0,
                 zphot = False, shape_noise = True,
                 eff_smoothing_grid = [1.0],current_magcut = None,weighted=False):

        self.current_magcut = current_magcut
        self.zphot = zphot
        self.shape_noise = shape_noise
        self.config_smoothing = config_smoothing
        self.pix = pix
        self.eff_smoothing_grid = eff_smoothing_grid
        self.weighted = weighted

        if self.zphot:
            import create_catalogue_zphot as catalogue
            self.catalogue = catalogue
            if self.weighted:
                self.f = open(c.opath+'pcc_fgweighted_zphot_table.txt', 'w')
            else:
                self.f = open('pcc_zphot_table.txt', 'w')
        else:
            if self.shape_noise:
                import create_catalogue_epsilon as catalogue
                self.catalogue = catalogue
                if self.weighted:
                    self.f = open(c.opath+'pcc_fgweighted_epsilon_table.txt','w')
                else:
                    self.f = open('pcc_epsilon_table.txt','w')
            else:
                import create_catalogue as catalogue
                self.catalogue = catalogue
                if self.weighted:
                    self.f = open(c.opath+'pcc_fgweighted_table.txt', 'w')
                else:
                    self.f = open('pcc_table.txt', 'w')

        self.f.write('Pearson_Corr\t Mag_Cut\t Gauss_Smooth\t Pixel_Scale\n')
        
    #Helper function to convolve with Gauss and write to new fits file
    def convolve_maps(self,some_fit,some_map,updated_filename,isFg):
        index = -1
        for sigma in self.return_smoothing_grid:
            index += 1
            print sigma
            if sigma > 0:
                if isFg:
                    SmoothKappa = filter.gaussian_filter(some_map,sigma)
                    #SmoothKappa = filter.gaussian_filter(some_fit[0].data,sigma)
                    SmoothMask = filter.gaussian_filter(np.array(some_fit[1].data).astype("float"),sigma)
                else:
                    SmoothKappa = filter.gaussian_filter(some_map,sigma)
                    #SmoothKappa = filter.gaussian_filter(some_fit[0].data,sigma)
                    SmoothMask = filter.gaussian_filter(np.array(some_fit[2].data).astype("float"),sigma)    

                Mask_ma = self.trim_mask(SmoothMask)
                SmoothKappa_masked = np.ma.masked_array(SmoothKappa, mask=Mask_ma)
                self.save_fits_image(SmoothKappa,c.opath+updated_filename+"_smooth"+str(self.eff_smoothing_grid[index])+".fits")
                self.save_fits_image(Mask_ma,c.opath+"Mask_"+updated_filename+"_smooth"+str(self.eff_smoothing_grid[index])+".fits") 
                    
            else:
                self.save_fits_image(some_map,c.opath+updated_filename+"_smooth"+str(self.eff_smoothing_grid[index])+".fits")
                self.save_fits_image(some_fit[-1].data,c.opath+"Mask_"+updated_filename+"_smooth"+str(self.eff_smoothing_grid[index])+".fits")

    def save_fits_image(self,image,filename):
        hdu = pf.PrimaryHDU(image)
        if os.path.exists(filename):
            os.remove(filename)
        hdu.writeto(filename)

    def trim_mask(self,SmoothMask):
        Mask = SmoothMask
        Mask[(SmoothMask < .9)] = 0
        Mask[(SmoothMask >= .9)] = 1
        return Mask

    def pcc(self):
        for sigma in self.eff_smoothing_grid:
            fg_map_p = pf.open(c.opath+"fgmap_fg"+self.filename+"_smooth"+str(sigma)+".fits")
            k_map_p = pf.open(c.opath+"kappamap"+self.filename+"_smooth"+str(sigma)+".fits")
            Masks = pf.open(c.opath+"Mask_kappamap"+self.filename+"_smooth"+str(sigma)+".fits")

            ff = np.ma.masked_array(fg_map_p[0].data.ravel(), mask=np.logical_not(Masks[0].data.ravel()))
            kk = np.ma.masked_array(k_map_p[0].data.ravel(), mask=np.logical_not(Masks[0].data.ravel()))

            corr = np.ma.corrcoef(ff,kk)
            print "Cross Corr np.ma:"
            print corr
            self.f.write(str(corr[0,1])+"\t"+str(self.current_magcut)+"\t"+str(np.sqrt(sigma**2+self.config_smoothing**2))+"\t"+str(self.pix)+"\n")

    #get effective smoothing by adding in quadtrature
    def get_smoothing(self):
        #self.return_smoothing_grid = self.eff_smoothing_grid
        config_smoothing_grid = np.resize(self.config_smoothing,(len(self.eff_smoothing_grid)))
        self.eff_smoothing_grid = np.array(self.eff_smoothing_grid)
        #index = -1
        #for eff in self.eff_smoothing_grid:
        #    print eff
        #    print self.config_smoothing
        #    print np.sqrt(eff**2 - self.config_smoothing**2)
        #    index = index+1
        #    self.return_smoothing_grid[index] = np.sqrt(eff**2 - self.config_smoothing**2)
        self.return_smoothing_grid = np.sqrt(self.eff_smoothing_grid**2 - config_smoothing_grid**2)
        return
    
    def analyze_ks_output(self):
    
        k_fit = pf.open(c.opath+"kappamap"+self.filename+".fits")
        k_map = k_fit[0].data
        fg_fit = pf.open(c.opath+"fgmap_fg"+self.filename+".fits")
        
        if self.weighted:
            fg_mapl = np.load(c.opath+"kappa_predicted"+self.filename+".npz")
            fg_map = fg_mapl['kappa']
        else:
            fg_map = fg_fit[0].data
            
        self.get_smoothing()
        
        self.convolve_maps(k_fit,k_map,"kappamap"+self.filename,False)
        self.convolve_maps(fg_fit,fg_map,"fgmap_fg"+self.filename,True)
        
        self.pcc()

    def movefiles(self):
        os.system("mv kappamap"+self.filename+".fits "+c.opath+"kappamap"+self.filename+".fits")
        os.system("mv fgmap_fg"+self.filename+".fits "+c.opath+"fgmap_fg"+self.filename+".fits")

    def run_pix_ks_mapping(pixel_scale_grid):
        for pix in pixel_scale_grid:
            filename = "_im3shape_r_"+str(pix)+"_"+str(self.config_smoothing)+"_g1"
            print "Running KS_MAPPING.py for pixel size: "+str(pix)
            edit_config(pix,self.config_smoothing)
            self.filename = "_im3shape_r_"+str(pix)+"_"+str(self.config_smoothing)+"_g1"
            self.run_ks_mapping()
            self.analyze_ks_output()
        return

    def run_smooth_ks_mapping():
        edit_config(self.pix,self.config_smoothing)
        run_ks_mapping()
        self.filename = "_im3shape_r_"+str(self.pix)+"_"+str(self.config_smoothing)+"_g1"
        self.movefiles()##Put files in out dir                                                          
        self.analyze_ks_output()
        os.system("mv ./out/density.npz "+c.bigfilepath+"density.npz")
        return

    def run_ks_mapping():
        if self.weighted:
            a = os.popen("/home/vinu/software/Python2.7/bin/python example_with_weights.py")
        else:
            a = os.popen("/home/vinu/software/Python2.7/bin/python example.py")
            #a = os.popen("/home/vinu/software/Python2.7/bin/python example2.py")    

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
        


if __name__=='__main__':
    
    #THESE ARE THE PARAMETERS YOU CAN TWEAK!############################
    eff_smoothing_grid = [0.5,0.75,1,1.5,2,4,8,15,20,30,40,50,60]
    pixel_scale_grid = [1.0,2.0,3.0,4.0,5.0,10.0,15.0,20.0,40.0,60.0]
    config_smoothing = 0.5
    pix = 1.0
    zphot = False
    shape_noise = True
    weighted = True
    mag_cut_grid = [23.0]
    
    pix_or_smooth = 'pix' #How do you want? Gaussian Smoothing or Pixelization?
    ####################################################################
    
    cc = cc(config_smoothing,pix,zphot,shape_noise,eff_smoothing_grid,weighted,
            pix_or_smooth=pix_or_smooth)        

    for current_mag_cut in mag_cut_grid:
        cc.current_magcut = current_mag_cut
        print "Magnitude Cut: "+str(current_mag_cut)
        cc.catalogue.run_catalogue(current_mag_cut,c.ipath,"/data3/data2/home/dbrout/")
        filename = "_im3shape_r_"+str(pix)+"_"+str(config_smoothing)+"_g1"
        
        if pix_or_smooth == 'pix':
            cc.run_pix_ks_mapping(pixel_scale_grid)
        if pix_or_smooth == 'smooth':
            cc.run_smooth_ks_mapping()

        #edit_config(pix,config_smoothing)
        #run_ks_mapping(weighted)
        #filename = "_im3shape_r_"+str(pix)+"_"+str(config_smoothing)+"_g1"
        #cc.filename = filename
        #cc.current_magcut = current_mag_cut
        #cc.movefiles()##Put files in out dir
        #cc.analyze_ks_output(weighted)
        #os.system("mv ./out/density.npz "+c.bigfilepath+"density.npz")
    cc.f.close()
