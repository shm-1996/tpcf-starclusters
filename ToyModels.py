from header import * 

"""
Tools for computing two-point correlation functions.
"""

import warnings
import numpy as np
from sklearn.neighbors import BallTree
from astroML.utils import check_random_state
from TPCF import linear_function,onepowerlaw_function,linear_truncation
from matplotlib.transforms import Bbox
from Fractal_gadi import * 

# Check if scikit-learn's two-point functionality is available.
# This was added in scikit-learn version 0.14
try:
    from sklearn.neighbors import KDTree
    sklearn_has_two_point = True
except ImportError:
    import warnings
    sklearn_has_two_point = False

#TPCF Functions


#############################################################################################
# Toy Model 1 : Spiral Arms
def Create_Spiral(r_thickness,Pitch=10.0,r0=1.0,nsamples=10000):

    #Parameters
    r_0 = r0 # Radius at theta=0.0
    dr = (r_thickness/2.)*np.random.uniform(size=nsamples)
    PA = np.deg2rad(Pitch)  # Pitch angle of spiral arms
    #Generate logarithmic spirals
    # r = r_0exp(theta*tan(PA))
    theta = np.random.uniform(low=-2.0*np.pi,high=2.0*np.pi,size=nsamples)
    #x = rcos(theta)
    x1= (r_0*np.exp(theta*np.tan(PA))+dr)*np.cos(theta)
    y1= (r_0*np.exp(theta*np.tan(PA))+dr)*np.sin(theta)

    #Generate the other side
    x2= (r_0*np.exp(theta*np.tan(PA))+dr)*np.cos(theta+np.pi)
    y2= (r_0*np.exp(theta*np.tan(PA))+dr)*np.sin(theta+np.pi)

    xr,yr = np.append(x1,x2),np.append(y1,y2)
    xr=-xr
    return xr,yr

def Spiral_TPCF(xr,yr,bins):
    data = np.asarray((xr,yr),order='F').T
    corr_lz,dcorr_lz = bootstrap_two_point(data, bins, Nbootstrap=30,
                            method='landy-szalay', return_bootstraps=False,
                            random_state=None)
    return corr_lz,dcorr_lz


def Analyse_Spiral(thickness,Pitch=10,compute_TPCF=False,save=False):
    filename = '../Toy_Models/Spiral_%02d'%(thickness*10)
    xr,yr = Create_Spiral(thickness,Pitch=Pitch)
    bins = np.logspace(np.log10(thickness/100.),np.log10(thickness*2.0),30)
    if(compute_TPCF is True):
        print("Computing TPCF")
        corr_lz,dcorr_lz = Spiral_TPCF(xr,yr,bins)
        saveObj(corr_lz,'{}_01_corr'.format(filename))
        saveObj(dcorr_lz,'{}_01_dcorr'.format(filename))
    else:
        print("Reading")
        corr_lz = loadObj('{}_01_corr'.format(filename))
        dcorr_lz = loadObj('{}_01_dcorr'.format(filename))
    print("Fitting")
    #Fit line to this
    from TPCF import linear_truncation
    separation_bins = (bins[1:]+bins[:-1])/2.
    indices = np.where(corr_lz>0.0)
    dcorr_lz = dcorr_lz[indices]
    separation_bins = separation_bins[indices]
    corr_lz = corr_lz[indices]
    
    bounds = ([-10.0,-3.0,np.min(separation_bins)],[10.0,0.0,np.max(separation_bins)*3.0])
    popt,pcov = curve_fit(linear_truncation,separation_bins,
                        np.log(1+corr_lz),sigma=dcorr_lz/corr_lz,bounds=bounds)
    perr= np.sqrt(np.diag(pcov))
    
    print("Plotting")
    fig,axs = plt.subplots(nrows=2,figsize=(4,8))

    axs[0].scatter(xr,yr,s=0.4,alpha=0.2,color='red')


    #Landy-Szalay version
    
    axs[1].errorbar(separation_bins,1+corr_lz,yerr=dcorr_lz,
                 fmt='.-')
    axs[1].axvline(thickness,ls=':',label=r'$r_\mathrm{thick}$'+r'$= {}$'.format(thickness))

    axs[1].plot(separation_bins,np.exp(linear_truncation(separation_bins,
                                                        popt[0],popt[1],popt[2])),
            ls='--',label=r'$\alpha = {:3.2f} \pm {:3.2f}$'.format(popt[1],perr[1]))

    axs[1].axvline(popt[2],ls='--',color='#F4271C',
                  label=r'$\theta_c = {:3.2f} \pm {:3.2f}$'.format(popt[2],perr[2]))

    axs[1].set_xlabel(r"$\Delta x \, (\mathrm{pixels})$")
    axs[1].set_ylabel(r"$1+\omega_{\mathrm{LS}}\left(\theta \right)$")
    axs[1].legend(loc='best')
    axs[1].set_xscale('log')
    axs[1].set_yscale('log')
    if(save):
        plt.savefig('{}'.format(filename),bbox_inches='tight')
    else:
        plt.show()


#############################################################################################
# Toy Model 2 : Discrete Star Forming Regions 
def Create_SF(r_SF,n_regions=30,n_points=100):
    xr = np.zeros((n_regions,n_points))
    yr = np.zeros((n_regions,n_points))

    #Choose random centres 
    centre_x = np.random.uniform(0.0,1.0,size=n_regions)
    centre_y = np.random.uniform(0.0,1.0,size=n_regions)
    centres = np.vstack((centre_x,centre_y)).T
    
    #Fill circular region randomly around each centre
    for i in range(n_regions):
        centre = centres[i]
        r = np.random.uniform(size=n_points)
        theta = np.random.uniform(low=0, high=2*np.pi, size=n_points)
        x = np.sqrt(r) * r_SF * np.cos(theta)
        y = np.sqrt(r) * r_SF * np.sin(theta)
        centre_x[0],centre_y[0] = 0,0
        x = x+centre_x[i]
        y = y+centre_y[i]
        xr[i] = x
        yr[i] = y



    xr = xr.flatten()
    yr = yr.flatten()
    return xr,yr
    
    

def SF_TPCF(xr,yr,bins):
    data = np.asarray((xr,yr),order='F').T
    corr_lz,dcorr_lz = bootstrap_two_point(data, bins, Nbootstrap=30,
                        method='landy-szalay', return_bootstraps=False,
                        random_state=None)
    return corr_lz,dcorr_lz




def Analyse_SF(r_SF,n_regions=30,n_points=100,compute_TPCF=False,
              save=False):
    filename = '../Toy_Models/SF_Regions/Omega1/SF_{}'.format(r_SF)
    xr,yr = Create_SF(r_SF,n_regions,n_points)
    bins = np.logspace(np.log10(r_SF/100.),np.log10(r_SF*5.0),40)
    
    if(compute_TPCF is True):
        print("Computing TPCF")
        corr_lz,dcorr_lz = SF_TPCF(xr,yr,bins)
        saveObj(corr_lz,'{}_corr'.format(filename))
        saveObj(dcorr_lz,'{}_dcorr'.format(filename))
    else:
        print("Reading")
        corr_lz = loadObj('{}_corr'.format(filename))
        dcorr_lz = loadObj('{}_dcorr'.format(filename))
        
    print("Fitting")
    #Fit line to this
    from TPCF import linear_truncation
    separation_bins = (bins[1:]+bins[:-1])/2.
    indices = np.where(corr_lz>0.0)
    dcorr_lz = dcorr_lz[indices]
    separation_bins = separation_bins[indices]
    corr_lz = corr_lz[indices]
    
    bounds = ([-10.0,-3.0,np.min(separation_bins)],[10.0,0.0,np.max(separation_bins)*3.0])
    popt,pcov = curve_fit(linear_truncation,separation_bins,
                    np.log(1+corr_lz),sigma=dcorr_lz/corr_lz,bounds=bounds)
    perr= np.sqrt(np.diag(pcov))
    
    print("Plotting")
    fig,axs = plt.subplots(nrows=2,figsize=(4,8))

    axs[0].scatter(xr,yr,s=0.4,alpha=0.2,color='red')


    #Landy-Szalay version
    
    axs[1].errorbar(separation_bins,1+corr_lz,yerr=dcorr_lz,
                 fmt='.-',lw=0.2)
    axs[1].axvline(r_SF,ls=':',label=r'$r_\mathrm{SF}$'+r'$= {}$'.format(r_SF))

    axs[1].plot(separation_bins,np.exp(linear_truncation(separation_bins,
                                                        popt[0],popt[1],popt[2])),
            ls='--',label=r'$\alpha = {:3.2f} \pm {:3.2f}$'.format(popt[1],perr[1]))

    axs[1].axvline(popt[2],ls='--',color='#F4271C',
                  label=r'$\theta_c = {:3.2f} \pm {:3.2f}$'.format(popt[2],perr[2]))

    axs[1].set_xlabel(r"$\Delta x \, (\mathrm{pixels})$")
    axs[1].set_ylabel(r"$1+\omega_{\mathrm{LS}}\left(\theta \right)$")
    axs[1].legend(loc='best')
    axs[1].set_xscale('log')
    axs[1].set_yscale('log')
    if(save):
        plt.savefig('{}.pdf'.format(filename),bbox_inches='tight')
        plt.close(fig)
    else:
        plt.show()

#############################################################################################
# Toy Model 3 : Discrete Star Forming Regions with overlying Poisson distribution         

def Create_Poisson(lambda_poisson=10000,xMin=0.0,xMax=1.0,yMin=0.0,yMax=1.0):
    #Simulation window parameters
    xDelta=xMax-xMin;yDelta=yMax-yMin; #rectangle dimensions
    areaTotal=xDelta*yDelta;

    #Point process parameters
    lambda0=lambda_poisson; #intensity (ie mean density) of the Poisson process

    #Simulate Poisson point process
    numbPoints = scipy.stats.poisson( lambda0*areaTotal ).rvs()#Poisson number of points
    xr = xDelta*scipy.stats.uniform.rvs(0,1,((numbPoints,1)))+xMin#x coordinates of Poisson points
    yr = yDelta*scipy.stats.uniform.rvs(0,1,((numbPoints,1)))+yMin#y coordinates of Poisson points
    xr = xr.flatten()
    yr = yr.flatten()
    return xr,yr
    

def Create_SF_Poisson(r_SF,poisson_fraction=1.0,n_regions=100,
                      n_points=100):
    n_sf = n_regions*n_points
    n_poisson = poisson_fraction*n_sf
    xr_sf, yr_sf = Create_SF(r_SF,n_regions,n_points)
    xr_poisson,yr_poisson = Create_Poisson(lambda_poisson=n_poisson)
    xr,yr = np.append(xr_sf,xr_poisson),np.append(yr_sf,yr_poisson)
    return xr,yr


def Analyse_SF_Poisson(r_SF,poisson_fraction=1.0,
                       n_regions=30,n_points=100,compute_TPCF=False,
                      save=False):
    filename = '../Toy_Models/SF_Poisson_Regions/SF_%02d'%(r_SF*100)
    xr,yr = Create_SF_Poisson(r_SF,poisson_fraction,n_regions,n_points)
    bins = np.logspace(np.log10(r_SF/100.),np.log10(r_SF*4.0),30)
    
    if(compute_TPCF is True):
        print("Computing TPCF")
        corr_lz,dcorr_lz = SF_TPCF(xr,yr,bins)
        saveObj(corr_lz,'{}_corr'.format(filename))
        saveObj(dcorr_lz,'{}_dcorr'.format(filename))
    else:
        print("Reading")
        corr_lz = loadObj('{}_corr'.format(filename))
        dcorr_lz = loadObj('{}_dcorr'.format(filename))
        
    print("Fitting")
    #Fit line to this
    from TPCF import linear_truncation
    separation_bins = (bins[1:]+bins[:-1])/2.
    indices = np.where(corr_lz>0.0)
    dcorr_lz = dcorr_lz[indices]
    separation_bins = separation_bins[indices]
    corr_lz = corr_lz[indices]
    
    bounds = ([-10.0,-3.0,np.min(separation_bins)],[10.0,0.0,np.max(separation_bins)*3.0])
    popt,pcov = curve_fit(linear_truncation,separation_bins,
                        np.log(corr_lz),sigma=dcorr_lz/corr_lz,bounds=bounds)
    perr= np.sqrt(np.diag(pcov))
    
    print("Plotting")
    fig,axs = plt.subplots(nrows=2,figsize=(4,8))

    axs[0].scatter(xr,yr,s=0.4,alpha=0.2,color='red')


    #Landy-Szalay version
    
    axs[1].errorbar(separation_bins,corr_lz,yerr=dcorr_lz,
                 fmt='.-')
    axs[1].axvline(r_SF,ls=':',label=r'$r_\mathrm{SF}$'+r'$= {}$'.format(r_SF))

    axs[1].plot(separation_bins,np.exp(linear_truncation(separation_bins,
                                                        popt[0],popt[1],popt[2])),
            ls='--',label=r'$\alpha = {:3.2f} \pm {:3.2f}$'.format(popt[1],perr[1]))

    axs[1].axvline(popt[2],ls='--',color='#F4271C',
                  label=r'$\theta_c = {:3.2f} \pm {:3.2f}$'.format(popt[2],perr[2]))

    axs[1].set_xlabel(r"$\Delta x \, (\mathrm{pixels})$")
    axs[1].set_ylabel(r"$\omega_{\mathrm{LS}}\left(\theta \right)$")
    axs[1].legend(loc='best')
    axs[1].set_xscale('log')
    axs[1].set_yscale('log')
    if(save):
        plt.savefig('{}'.format(filename),bbox_inches='tight')
    else:
        plt.show()

#############################################################################################


def Combined_Spiral_Plots():
    Pitch = 40.0
    Analyse_Spiral(0.1,Pitch = Pitch,compute_TPCF=True,save=True)
    Analyse_Spiral(0.2,Pitch = Pitch,compute_TPCF=True,save=True)
    Analyse_Spiral(0.4,Pitch = Pitch,compute_TPCF=True,save=True)
    Analyse_Spiral(0.8,Pitch = Pitch,compute_TPCF=True,save=True)

    fig,axs = plt.subplots(nrows=2,ncols=4,figsize=(16,8))

    nos = [1,2,4,8]
    for i in range(0,4):
        no = nos[i]
        thickness = np.float(no)/10.0
        file = '../Toy_Models/Spiral_%02d_01_corr'%no
        corr_lz = loadObj(file)
        file = '../Toy_Models/Spiral_%02d_01_dcorr'%no
        dcorr_lz = loadObj(file)
        xr,yr = Create_Spiral(thickness,Pitch=Pitch)
        axs[0,i].scatter(xr,yr,s=0.4,alpha=0.2,color='red')

        
        #Fit
        from TPCF import linear_truncation
        bins = np.logspace(np.log10(thickness/100.),np.log10(thickness*2.0),30)
        separation_bins = (bins[1:]+bins[:-1])/2.
        indices = np.where(corr_lz>0.0)
        dcorr_lz = dcorr_lz[indices]
        separation_bins = separation_bins[indices]
        corr_lz = corr_lz[indices]
        
        bounds = ([-10.0,-3.0,np.min(separation_bins)],[10.0,0.0,np.max(separation_bins)*3.0])
        popt,pcov = curve_fit(linear_truncation,separation_bins,
                            np.log(corr_lz),sigma=dcorr_lz/corr_lz,bounds=bounds)
        perr= np.sqrt(np.diag(pcov))
        
        axs[1,i].errorbar(separation_bins,corr_lz,yerr=dcorr_lz,
                     fmt='.-')
        axs[1,i].axvline(thickness,ls=':',label=r'$r_\mathrm{thick}$'+r'$= {}$'.format(thickness))

        axs[1,i].plot(separation_bins,np.exp(linear_truncation(separation_bins,
                                                            popt[0],popt[1],popt[2])),
                ls='--',label=r'$\alpha = {:3.2f} \pm {:3.2f}$'.format(popt[1],perr[1]))

        axs[1,i].axvline(popt[2],ls='--',color='#F4271C',
                      label=r'$\theta_c = {:3.2f} \pm {:3.2f}$'.format(popt[2],perr[2]))

        axs[1,i].set_xlabel(r"$\Delta x \, (\mathrm{pixels})$")
        if(i == 0):
            axs[1,i].set_ylabel(r"$\omega_{\mathrm{LS}}\left(\theta \right)$")
        axs[1,i].legend(loc='best')
        axs[1,i].set_xscale('log')
        axs[1,i].set_yscale('log')
        
    plt.savefig('../Toy_Models/Pitch40_Spirals',bbox_inches='tight') 



if __name__ == "__main__":
    print("Done")
       



