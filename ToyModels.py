from header import * 

"""
Tools for computing two-point correlation functions.
"""

import warnings
import numpy as np
from sklearn.neighbors import BallTree
from astroML.utils import check_random_state

# Check if scikit-learn's two-point functionality is available.
# This was added in scikit-learn version 0.14
try:
    from sklearn.neighbors import KDTree
    sklearn_has_two_point = True
except ImportError:
    import warnings
    sklearn_has_two_point = False


def uniform_sphere(RAlim, DEClim, size=1):
    """Draw a uniform sample on a sphere

    Parameters
    ----------
    RAlim : tuple
        select Right Ascension between RAlim[0] and RAlim[1]
        units are degrees
    DEClim : tuple
        select Declination between DEClim[0] and DEClim[1]
    size : int (optional)
        the size of the random arrays to return (default = 1)

    Returns
    -------
    RA, DEC : ndarray
        the random sample on the sphere within the given limits.
        arrays have shape equal to size.
    """
    zlim = np.sin(np.pi * np.asarray(DEClim) / 180.)

    z = zlim[0] + (zlim[1] - zlim[0]) * np.random.random(size)
    DEC = (180. / np.pi) * np.arcsin(z)
    RA = RAlim[0] + (RAlim[1] - RAlim[0]) * np.random.random(size)

    return RA, DEC


def ra_dec_to_xyz(ra, dec):
    """Convert ra & dec to Euclidean points

    Parameters
    ----------
    ra, dec : ndarrays

    Returns
    x, y, z : ndarrays
    """
    sin_ra = np.sin(ra * np.pi / 180.)
    cos_ra = np.cos(ra * np.pi / 180.)

    sin_dec = np.sin(np.pi / 2 - dec * np.pi / 180.)
    cos_dec = np.cos(np.pi / 2 - dec * np.pi / 180.)

    return (cos_ra * sin_dec,
            sin_ra * sin_dec,
            cos_dec)


def angular_dist_to_euclidean_dist(D, r=1):
    """convert angular distances to euclidean distances"""
    return 2 * r * np.sin(0.5 * D * np.pi / 180.)


def two_point(data, bins, method='standard',
              data_R=None, random_state=None):
    """Two-point correlation function

    Parameters
    ----------
    data : array_like
        input data, shape = [n_samples, n_features]
    bins : array_like
        bins within which to compute the 2-point correlation.
        shape = Nbins + 1
    method : string
        "standard" or "landy-szalay".
    data_R : array_like (optional)
        if specified, use this as the random comparison sample
    random_state : integer, np.random.RandomState, or None
        specify the random state to use for generating background

    Returns
    -------
    corr : ndarray
        the estimate of the correlation function within each bin
        shape = Nbins
    """
    data = np.asarray(data)
    bins = np.asarray(bins)
    rng = check_random_state(random_state)

    if method not in ['standard', 'landy-szalay']:
        raise ValueError("method must be 'standard' or 'landy-szalay'")

    if bins.ndim != 1:
        raise ValueError("bins must be a 1D array")

    if data.ndim == 1:
        data = data[:, np.newaxis]
    elif data.ndim != 2:
        raise ValueError("data should be 1D or 2D")

    n_samples, n_features = data.shape
    Nbins = len(bins) - 1

    # shuffle all but one axis to get background distribution
    if data_R is None:
        data_R = data.copy()
        for i in range(n_features - 1):
            rng.shuffle(data_R[:, i])
    else:
        data_R = np.asarray(data_R)
        if (data_R.ndim != 2) or (data_R.shape[-1] != n_features):
            raise ValueError('data_R must have same n_features as data')

    factor = len(data_R) * 1. / len(data)

    if sklearn_has_two_point:
        # Fast two-point correlation functions added in scikit-learn v. 0.14
        KDT_D = KDTree(data)
        KDT_R = KDTree(data_R)

        counts_DD = KDT_D.two_point_correlation(data, bins)
        counts_RR = KDT_R.two_point_correlation(data_R, bins)

    else:
        warnings.warn("Version 0.3 of astroML will require scikit-learn "
                      "version 0.14 or higher for correlation function "
                      "calculations. Upgrade to sklearn 0.14+ now for much "
                      "faster correlation function calculations.")

        BT_D = BallTree(data)
        BT_R = BallTree(data_R)

        counts_DD = np.zeros(Nbins + 1)
        counts_RR = np.zeros(Nbins + 1)

        for i in range(Nbins + 1):
            counts_DD[i] = np.sum(BT_D.query_radius(data, bins[i],
                                                    count_only=True))
            counts_RR[i] = np.sum(BT_R.query_radius(data_R, bins[i],
                                                    count_only=True))

    DD = np.diff(counts_DD)
    RR = np.diff(counts_RR)

    # check for zero in the denominator
    RR_zero = (RR == 0)
    RR[RR_zero] = 1

    if method == 'standard':
        corr = factor ** 2 * DD / RR - 1
    elif method == 'landy-szalay':
        if sklearn_has_two_point:
            counts_DR = KDT_R.two_point_correlation(data, bins)
        else:
            counts_DR = np.zeros(Nbins + 1)
            for i in range(Nbins + 1):
                counts_DR[i] = np.sum(BT_R.query_radius(data, bins[i],
                                                        count_only=True))
        DR = np.diff(counts_DR)

        corr = (factor ** 2 * DD - 2 * factor * DR + RR) / RR

    corr[RR_zero] = np.nan

    return corr



def bootstrap_two_point(data, bins, Nbootstrap=10,
                        method='standard', return_bootstraps=False,
                        random_state=None):
    """Bootstrapped two-point correlation function

    Parameters
    ----------
    data : array_like
        input data, shape = [n_samples, n_features]
    bins : array_like
        bins within which to compute the 2-point correlation.
        shape = Nbins + 1
    Nbootstrap : integer
        number of bootstrap resamples to perform (default = 10)
    method : string
        "standard" or "landy-szalay".
    return_bootstraps: bool
        if True, return full bootstrapped samples
    random_state : integer, np.random.RandomState, or None
        specify the random state to use for generating background

    Returns
    -------
    corr, corr_err : ndarrays
        the estimate of the correlation function and the bootstrap
        error within each bin. shape = Nbins
    """
    data = np.asarray(data)
    bins = np.asarray(bins)
    rng = check_random_state(random_state)

    if method not in ['standard', 'landy-szalay']:
        raise ValueError("method must be 'standard' or 'landy-szalay'")

    if bins.ndim != 1:
        raise ValueError("bins must be a 1D array")

    if data.ndim == 1:
        data = data[:, np.newaxis]
    elif data.ndim != 2:
        raise ValueError("data should be 1D or 2D")

    if Nbootstrap < 2:
        raise ValueError("Nbootstrap must be greater than 1")

    n_samples, n_features = data.shape

    # get the baseline estimate
    corr = two_point(data, bins, method=method, random_state=rng)

    bootstraps = np.zeros((Nbootstrap, len(corr)))

    for i in range(Nbootstrap):
        indices = rng.randint(0, n_samples, n_samples)
        bootstraps[i] = two_point(data[indices, :], bins, method=method,
                                  random_state=rng)

    # use masked std dev in case of NaNs
    corr_err = np.asarray(np.ma.masked_invalid(bootstraps).std(0, ddof=1))

    if return_bootstraps:
        return corr, corr_err, bootstraps
    else:
        return corr, corr_err



def two_point_angular(ra, dec, bins, method='standard', random_state=None):
    """Angular two-point correlation function

    A separate function is needed because angular distances are not
    euclidean, and random sampling needs to take into account the
    spherical volume element.

    Parameters
    ----------
    ra : array_like
        input right ascention, shape = (n_samples,)
    dec : array_like
        input declination
    bins : array_like
        bins within which to compute the 2-point correlation.
        shape = Nbins + 1
    method : string
        "standard" or "landy-szalay".
    random_state : integer, np.random.RandomState, or None
        specify the random state to use for generating background

    Returns
    -------
    corr : ndarray
        the estimate of the correlation function within each bin
        shape = Nbins
    """
    ra = np.asarray(ra)
    dec = np.asarray(dec)
    rng = check_random_state(random_state)

    if method not in ['standard', 'landy-szalay']:
        raise ValueError("method must be 'standard' or 'landy-szalay'")

    if bins.ndim != 1:
        raise ValueError("bins must be a 1D array")

    if (ra.ndim != 1) or (dec.ndim != 1) or (ra.shape != dec.shape):
        raise ValueError('ra and dec must be 1-dimensional '
                         'arrays of the same length')

    n_features = len(ra)
    Nbins = len(bins) - 1

    # draw a random sample with N points
    ra_R, dec_R = uniform_sphere((min(ra), max(ra)),
                                 (min(dec), max(dec)),
                                 2 * len(ra))

    data = np.asarray(ra_dec_to_xyz(ra, dec), order='F').T
    data_R = np.asarray(ra_dec_to_xyz(ra_R, dec_R), order='F').T

    # convert spherical bins to cartesian bins
    bins_transform = angular_dist_to_euclidean_dist(bins)

    return two_point(data, bins_transform, method=method,
                     data_R=data_R, random_state=rng)



def bootstrap_two_point_angular(ra, dec, bins, method='standard',
                                Nbootstraps=10, random_state=None):
    """Angular two-point correlation function

    A separate function is needed because angular distances are not
    euclidean, and random sampling needs to take into account the
    spherical volume element.

    Parameters
    ----------
    ra : array_like
        input right ascention, shape = (n_samples,)
    dec : array_like
        input declination
    bins : array_like
        bins within which to compute the 2-point correlation.
        shape = Nbins + 1
    method : string
        "standard" or "landy-szalay".
    Nbootstraps : int
        number of bootstrap resamples
    random_state : integer, np.random.RandomState, or None
        specify the random state to use for generating background

    Returns
    -------
    corr : ndarray
        the estimate of the correlation function within each bin
        shape = Nbins
    dcorr : ndarray
        error estimate on dcorr (sample standard deviation of
        bootstrap resamples)
    bootstraps : ndarray
        The full sample of bootstraps used to compute corr and dcorr
    """
    ra = np.asarray(ra)
    dec = np.asarray(dec)
    rng = check_random_state(random_state)

    if method not in ['standard', 'landy-szalay']:
        raise ValueError("method must be 'standard' or 'landy-szalay'")

    if bins.ndim != 1:
        raise ValueError("bins must be a 1D array")

    if (ra.ndim != 1) or (dec.ndim != 1) or (ra.shape != dec.shape):
        raise ValueError('ra and dec must be 1-dimensional '
                         'arrays of the same length')

    n_features = len(ra)
    Nbins = len(bins) - 1
    data = np.asarray(ra_dec_to_xyz(ra, dec), order='F').T

    # convert spherical bins to cartesian bins
    bins_transform = angular_dist_to_euclidean_dist(bins)

    bootstraps = []

    for i in range(Nbootstraps):
        # draw a random sample with N points
        ra_R, dec_R = uniform_sphere((min(ra), max(ra)),
                                     (min(dec), max(dec)),
                                     2 * len(ra))

        data_R = np.asarray(ra_dec_to_xyz(ra_R, dec_R), order='F').T

        if i > 0:
            # random sample of the data
            ind = np.random.randint(0, data.shape[0], data.shape[0])
            data_b = data[ind]
        else:
            data_b = data

        bootstraps.append(two_point(data_b, bins_transform, method=method,
                                    data_R=data_R, random_state=rng))

    bootstraps = np.asarray(bootstraps)
    corr = np.mean(bootstraps, 0)
    corr_err = np.std(bootstraps, 0, ddof=1)

    return corr, corr_err, bootstraps



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
                        np.log(corr_lz),sigma=dcorr_lz/corr_lz,bounds=bounds)
    perr= np.sqrt(np.diag(pcov))
    
    print("Plotting")
    fig,axs = plt.subplots(nrows=2,figsize=(4,8))

    axs[0].scatter(xr,yr,s=0.4,alpha=0.2,color='red')


    #Landy-Szalay version
    
    axs[1].errorbar(separation_bins,corr_lz,yerr=dcorr_lz,
                 fmt='.-')
    axs[1].axvline(thickness,ls=':',label=r'$r_\mathrm{thick}$'+r'$= {}$'.format(thickness))

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
    filename = '../Toy_Models/SF_Regions/SF_%02d'%(r_SF*100)
    xr,yr = Create_SF(r_SF,n_regions,n_points)
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

if __name__ == "__main__":
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



