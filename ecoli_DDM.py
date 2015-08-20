#based on Cerbino and Vailati 2009
import sys
sys.path.append('/home/agoldfain/code/ecoli_DDM')
import get_difference as diff
import fit_pxs_for_tau as fpt
reload(diff)
reload(fpt)
import holopy as hp
from os.path import basename, exists
from glob import glob
import divide_background as db
from os import makedirs
from AaronFunctions import A_imsave

from holopy.core.math import fft

#first show normal background divided ecoli

folder='/home/agoldfain/group/agoldfain/data/iSCAT/150422_ecoli/'
#folder = '/tmp/'

name ='*ecoli_01*'
filepaths = glob(folder+name+'.fits')
constscalingsave = True

outbase = '/tmp/results/'
#outbase = None

delta_ts = arange(1,21)
crop=[121,313,138,332,None,None]
abssq = False
t_spacing = .01

for filepath in filepaths:
    subracted = None
    if outbase != None:
        outname = basename(filepath).split('.')[0]
        outfolder = outbase+outname+'/'
        if not exists(outfolder):
            makedirs(outfolder)
    else:
        outfolder = None

    [subtracted,data_header] = diff.get_difference_fits(filepath, delta_ts, outfolder = outfolder, constscalingsave = constscalingsave, crop=crop, return_difference = True, return_header = True, abssq = abssq)

#get fourier transform of difference and take absolute values
subtracted_ft = empty(subtracted.shape)
for i in arange(0,subtracted.shape[2]):
    if i%100 == 0:
        print('fourier transform ' + str(i+1) + ' of ' + str(subtracted.shape[2]))
    subtracted_ft[:,:,i,:] = abs(fft(subtracted[:,:,i,:]))**2
    subtracted[:,:,i,:] = abs(subtracted[:,:,i,:])**2

##################################################################################################
# fit an exponential to each pixel
ts = delta_ts*t_spacing #time spacings

#break video into sections to average over
n_steps = 1 #number of sections to use
n_avg = subtracted.shape[2]/n_steps #number of frames in each average


amps = empty((subtracted.shape[0],subtracted.shape[1],n_steps))
amps_ft = empty((subtracted.shape[0],subtracted.shape[1],n_steps))
taus = empty((subtracted.shape[0],subtracted.shape[1],n_steps))
taus_ft = empty((subtracted.shape[0],subtracted.shape[1],n_steps))

for i in arange(0,n_steps):
    print('avg number ' + str(i+1) + ' of ' + str(n_steps))
    #average video into sections
    i_min = i*n_avg
    i_max = min((i+1)*n_avg, subtracted.shape[2]-delta_ts.max())
    sub_mean = mean(subtracted[:,:,i_min:i_max,:],2)
    sub_mean_ft = mean(subtracted_ft[:,:,i_min:i_max,:],2)

    #fit exponential curve to data
    amps[:,:,i], taus[:,:,i], amps_std, taus_std =  fpt.fit_px_for_tau(sub_mean, ts, norm_amps = True)
    amps_ft[:,:,i], taus_ft[:,:,i], amps_ft_std, taus_ft_std =  fpt.fit_px_for_tau(sub_mean_ft, ts, norm_amps = True)

hp.show(amps)
#savefig(outfolder+'amplitude.tif' )
hp.show(taus)
#savefig(outfolder+'tau.tif' )

'''hp.show(amps_ft)
savefig(outfolder+'amplitude_ft.tif' )
hp.show(taus_ft)
savefig(outfolder+'tau_ft.tif' )

hp.show(amps_std)
savefig(outfolder+'amplitude_std.tif' )
hp.show(taus_std)
savefig(outfolder+'tau_std.tif' )

hp.show(amps_ft_std)
savefig(outfolder+'amplitude_ft_std.tif' )
hp.show(taus_ft_std)
savefig(outfolder+'tau_ft_std.tif' )'''

#Tau analysis from Cerbino and Trappe 2008
px_size = 6.5e-6/100
length_scale = 1./(44./192.*(1/px_size))# speckle size. Found 44/192 factor by looking at fourier transform of subtracted data before taking
                                        # absolute value squared. There is a clear high frequency cutoff at pixel 44 of the 192 pixels.

D_coef = where(taus<0.005, zeros(taus.shape), length_scale**2/(taus))
eta= 1.005e-3
kBT= 295*1.38e-23
diameter = kBT*taus/(length_scale**2*3*pi*eta)*1e9

#########################
#diff coefficeint using q



    
