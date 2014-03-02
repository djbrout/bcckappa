from pylab import *
fg_mapl = np.load("/home/dbrout/kappabias/bias/kappa_predicted.npz")
fg_map = fg_mapl['kappa']
figure(1)
print fg_map
imshow(fg_map)
colorbar()
savefig("kappa_predicted_pix_sourcez_fixzs8.pdf")
