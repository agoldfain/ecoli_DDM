#based on Cerbino and Trappe 2008
import sys
sys.path.append('/home/agoldfain/code/ecoli_DDM')
import get_difference as diff
reload(diff)
import holopy as hp
from os.path import basename, exists
from glob import glob
import divide_background as db
from os import makedirs
from AaronFunctions import A_imsave
from scipy.optimize import curve_fit

#first show normal background divided ecoli

#folder='/home/agoldfain/group/agoldfain/data/iSCAT/150422_ecoli/'
folder = '/tmp/'

name ='*ecoli_01*'
filepaths = glob(folder+name+'.fits')
constscalingsave = True

outbase = '/tmp/results/'
#outbase = None

delta_ts = arange(1,21)
crop=[121,313,138,332,None,None]
abssq = True
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
    
    
sub_mean = mean(subtracted[:,:,:-delta_ts.max(),:],2)

# fit an exponential to each pixel

def fit_func(t, a, tau):
    return a*(1 - exp(-tau*t) )


ts = delta_ts*t_spacing

amps = zeros(sub_mean.shape[:2]) 
taus = zeros(sub_mean.shape[:2])
amps_var = zeros(sub_mean.shape[:2]) 
taus_var = zeros(sub_mean.shape[:2])
 
for i in arange(sub_mean.shape[0]):
    for j in arange(sub_mean.shape[1]):
        popt, pcov = curve_fit(fit_func, ts, sub_mean[i,j,:], p0 = [sub_mean[i,j,:].max(), 100])
        #plot(ts, sub_mean[i,j,:], ts, fit_func(ts, *popt) )
        
        if isinstance(pcov, numpy.ndarray) and 1/popt[1]<200:
            amps_var[i,j] = pcov[0,0]
            taus_var[i,j] = pcov[1,1]
            amps[i,j] = popt[0]
            taus[i,j] = 1/popt[1]
        
hp.save(outfolder+'amplitude.tif', amps)
hp.save(outfolder+'tau.tif', taus)
    
    
    
    
