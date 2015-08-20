import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def fit_px_for_tau(data, ts, norm_amps = True):

    def fit_func(t, a, tau):
        return  a*(1 - np.exp(-t/tau) )


    amps = np.zeros(data.shape[:2]) 
    taus = np.zeros(data.shape[:2])
    #betas = np.zeros(data.shape[:2])
    amps_std = np.zeros(data.shape[:2]) 
    taus_std = np.zeros(data.shape[:2])
    #betas_std = np.zeros(data.shape[:2])
     
    for i in np.arange(data.shape[0]):
        for j in np.arange(data.shape[1]):
            popt, pcov = curve_fit(fit_func, ts, data[i,j,:], p0 = [data[i,j,:].max(), .01])
            if i ==-1:
                plt.figure()
                plt.plot(ts, data[i,j,:], 'b.', ts, fit_func(ts, *popt), 'b' )
                plt.savefig('/tmp/results/x'+str(i).zfill(3)+'y'+str(j).zfill(3)+'.tif')
                plt.close('all')
            
            if isinstance(pcov, np.ndarray) and popt[1]<1:
                amps_std[i,j] = np.sqrt(pcov[0,0])
                taus_std[i,j] = np.sqrt(pcov[1,1])
                amps[i,j] = popt[0]
                taus[i,j] = popt[1]


    if norm_amps:
        amps = amps/amps.max()
        amps_std = amps_std/amps.max()
    
    return (amps, taus, amps_std, taus_std)
