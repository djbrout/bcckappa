from pylab import *
import config as c
fg_mapl = np.load(c.opath+"/home/dbrout/bcckappa/out/kappa_predicted.npz")
fg_map = fg_mapl['kappa']
figure(1)
print fg_map
imshow(fg_map)
colorbar()
savefig(c.opath+"kappa_predicted_pix_fixed_sz.pdf")
