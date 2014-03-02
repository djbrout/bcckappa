import rdcol2
import matplotlib.pylab as plt
import numpy as np
from pylab import *
import sys

rcParams['legend.loc'] = 'best'



cols = rdcol2.read('pcc_fgweighted_table.txt',1,2)
x = np.arange(4)
ys = [i+x+(i*x)**2 for i in np.arange(4)]
colors = iter(cm.rainbow(np.linspace(0, 1, len(ys))))
for magc in np.unique(cols["Mag_Cut"]):
    if magc == 23.0:
        figure(1)
        scatter(np.array(cols["Gauss_Smooth"])[np.array(cols["Mag_Cut"]) == magc],
                np.array(cols["Pearson_Corr"])[np.array(cols["Mag_Cut"]) == magc],
                label = "Weighted Fg, Mag Cut: "+str(magc),color=next(colors))
    
plt.xlabel("Smoothing Scale (arcmin)")
plt.ylabel("Pearson_CC  ($\kappa$ x Fg)")
plt.title('Using TrueZ')
plt.ylim(0,1)


cols = rdcol2.read('pcc_table.txt',1,2)
for magc in np.unique(cols["Mag_Cut"]):
    if magc == 23.0:
        scatter(np.array(cols["Gauss_Smooth"])[np.array(cols["Mag_Cut"]) == magc],
                np.array(cols["Pearson_Corr"])[np.array(cols["Mag_Cut"]) == magc],
                label = "Unweighted Fg,  Mag Cut: "+str(magc),color=next(colors),marker='v')

legend(loc=4)
savefig("cc_vs_smoothing_and_magcut_and_weighting.pdf")

#Make plot of smoothing vs w/ & w/o photoz at 23.5 magcut
#make plot of smoothing vs w/ & w/o shape noise at 23.5 magcut
