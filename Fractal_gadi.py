from TPCF import linear_function,onepowerlaw_function,linear_truncation
from matplotlib.transforms import Bbox
#from ToyModels import *
from header import *
from joblib import Parallel, delayed
from itertools import repeat
import timeit
import time

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



#Some fixed parameters
############################################################################################
max_level = 10
max_level_low = 6
max_level_high = 14
fractal_dims = [0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]
sizes = [0.1,0.2,0.3,0.4,0.5,0.6]
filename = '../Toy_Models/Fractals/Fractal_Analysis.out'
max_points = 50000
MC_min = 2500

############################################################################################

# Toy Model 4 : Fractal Models

def Divide_Region(rand_no,prob,level,current_level_x,current_level_y):
    if(rand_no<prob):
        #Create 4 new centres for this region
        delta = 1./(2**(level+1))/2.
        new_x = [current_level_x-delta,
                current_level_x+delta,
                current_level_x-delta,
                current_level_x+delta]
        new_y = [current_level_y+delta,
                current_level_y+delta,
                current_level_y-delta,
                current_level_y-delta]
        return new_x, new_y
    else:

        return
    



def Create_Fractal(fractal_dim,baselevel=2,max_level=10):

    prob = 2**(fractal_dim-2)
    delta = 1./(2**baselevel)
    delta = 1./(2**(baselevel+1))
    current_level_x, current_level_y = np.meshgrid(np.arange(delta/2.,1.0,
                            delta),np.arange(delta/2.,1.0,delta),indexing='xy')
    current_level_x, current_level_y = current_level_x.flatten(),current_level_y.flatten()
    
    #Inherently serial
    for level in range(baselevel+1,max_level):
        nactive = np.size(current_level_x)
        rand_no = np.random.uniform(size=nactive)
        new_level_x = []
        new_level_y = []
        #Loop over active regions and add regions if condition satisfied
        
        #This part is parallelised
        region = np.arange(0,nactive,1)
        results = Parallel(n_jobs=-1,prefer='processes',verbose=0)(map(delayed(Divide_Region), 
            rand_no,repeat(prob),repeat(level),current_level_x,current_level_y))

        #Remove nones (i.e. non-selected regions) from list
        results = np.array(list(filter(None,results)))

        #Get flattened positions of new level
        new_level_x = results[:,0].flatten()
        new_level_y = results[:,1].flatten()
        
        current_level_x = new_level_x
        current_level_y = new_level_y 
        
    
    return new_level_x,new_level_y    

def Compute_TPCF_Fractal(xfractal,yfractal,no_bins=40,limits=0.0):
    
    indices_1 = np.logical_and(xfractal>limits,yfractal>limits)
    indices_2 = np.logical_and(xfractal<1.-limits,yfractal<1.-limits)
    indices = np.logical_and(indices_1,indices_2)
    size = (1.-limits*2.0)/4.
    bins = np.logspace(np.log10(0.0001),np.log10(size),no_bins+1)
    data = np.asarray((xfractal[indices],yfractal[indices]),order='F').T
    corr_lz,dcorr_lz = bootstrap_two_point(data, bins, Nbootstrap=30,
                            method='landy-szalay', return_bootstraps=False,
                            random_state=None)
    return bins,corr_lz,dcorr_lz

def Analyse_Fractal(fractal_dim,limits=0.0,no_bins=40,max_level=10):

    #Generate Fractal 
    start = time.time()
    xfractal,yfractal = Create_Fractal(fractal_dim,max_level=max_level)
    end = time.time()
    print("Fractal Created. Time taken: {} seconds".format(end-start))

    #Plot Positions
    plt.clf()
    plt.scatter(xfractal,yfractal,s=0.05,c='blue')
    plt.savefig('../Toy_Models/Fractals/Pos_D_{}'.format(int(fractal_dim*10)),
                bbox_inches='tight')
    plt.close()
    pkl_obj = [xfractal,yfractal]
    saveObj(pkl_obj,'../Toy_Models/Fractals/XYPos_{}'.format(int(fractal_dim*10)))

    #Restrict number of points to a max 
    if(np.size(xfractal)>max_points):
        print("Number of points in fractal {} exceeds max points.".format(np.size(xfractal)))
        print("Culling points randomly to bring num_points = {}".format(max_points))
        indices = np.random.randint(0,np.size(xfractal),max_points)
        xfractal = xfractal[indices]
        yfractal = yfractal[indices]

    start = time.time()
    bins,corr_lz,dcorr_lz = Compute_TPCF_Fractal(xfractal,yfractal,no_bins=no_bins)
    end = time.time()
    print("TPCF computed.Time taken: {} seconds".format(end-start))
    save_obj = [xfractal,yfractal,bins,corr_lz,dcorr_lz]
    saveObj(save_obj,'../Toy_Models/Fractals/D_{}'.format(int(fractal_dim*10)))

    #Fit and Plot TPCF
    fig = plot_omega1(bins,corr_lz,dcorr_lz)
    fig.savefig('../Toy_Models/Fractals/TPCF_D_{}'.format(int(fractal_dim*10)),
                bbox_inches='tight')
    plt.clf()
    plt.close(fig)

    #Fit exponential truncation to get truncation radius
    fig = plot_omega(bins,corr_lz,dcorr_lz)
    fig.savefig('../Toy_Models/Fractals/TPCF_cutoff_{}'.format(int(fractal_dim*10)),bbox_inches='tight')        
    plt.clf()
    plt.close(fig)



#############################################################################################
# Toy Model 5 : Fractals with Artificial Boundaries to test boundary effects

def check_contains(x,y,bbox):
    return bbox.contains(x,y)

def MC_fractal(xfractal,yfractal,size,angle=0.0,no_bins=20,N_MC=30,debug=False):

    n = 0
    attempts = 0
    max_attempts = 5
    vfunc = np.vectorize(check_contains,excluded=['bbox'])
    corr_MC = np.zeros((N_MC,no_bins))
    dcorr_MC = np.zeros((N_MC,no_bins))
    pbar = tqdm.tqdm(total=N_MC)
    while n<N_MC:
        #Random corner
        xlow = np.random.uniform(low=size,high=1.-size)
        ylow = np.random.uniform(low=size,high=1.-size)

        left,bottom,width,height = (xlow,ylow,size,size)
        FOV = mpl.patches.Rectangle((left,bottom),width,height,
                                    angle=angle,fill=False,ls='--',
                                   lw=2.0,color='#F9004A')

        bbox = Bbox.from_bounds(left, bottom, width, height)
        
        #Plot
        if(n%10 == 0 and debug == True):
            fig,axs = plt.subplots(ncols=1)
            axs.scatter(xfractal,yfractal,alpha=0.2,s=0.9)
            axs.add_artist(FOV)
            plt.savefig('../Toy_Models/Fractals/Monte_Carlo_Boundaries/Size_{}_{}'.format(int(size*10),np.int(n/10)),
                       bbox_inches='tight')
            plt.clf()
            plt.close(fig)
        
        #Extract out this part of fractal

        indices_true = vfunc(xfractal,yfractal,bbox=bbox)
        x_sec, y_sec = xfractal[indices_true],yfractal[indices_true]
        if(np.size(x_sec)<MC_min):
            #print("Not enough points ({}). Skipping this sample.".format(np.size(x_sec)))
            continue

        #Compute TPCF 
        bins = np.logspace(np.log10(0.0001),np.log10(size/2.),no_bins+1)
        data = np.asarray((xfractal[indices_true],yfractal[indices_true]),order='F').T

        try:
            corr_lz,dcorr_lz = bootstrap_two_point(data, bins, Nbootstrap=50,
                                method='landy-szalay', return_bootstraps=False,
                                random_state=None)
        except ZeroDivisionError: 
                attempts +=1
                if(attempts==5):
                    raise Exception("Maximum attempts for Monte-Carlo Boundaries attempted")
                continue

        corr_MC[n] = corr_lz
        dcorr_MC[n] = dcorr_lz

        n = n+1
        pbar.update(1)
        attempts = 0
    pbar.close()
    return bins,corr_MC,dcorr_MC

def TPCF_MC(xfractal,yfractal,size,angle=0.0,no_bins=20,N_MC=100,MC_directory=None):
    
    if(MC_directory == None):
        MC_directory = '../Toy_Models/Fractals/Monte_Carlo_Boundaries'
    MC_directory = os.path.abspath(MC_directory)

    if(not os.path.exists(MC_directory)):
            os.makedirs(MC_directory)

    bins,corr_MC,dcorr_MC = MC_fractal(xfractal,yfractal,size,angle,no_bins,N_MC)
    corr_avg = np.mean(corr_MC,axis=0)
    dcorr_avg = np.mean(dcorr_MC,axis=0)

    pkl_obj = [bins,corr_avg,dcorr_avg]
    saveObj(pkl_obj,MC_directory+'/Size_{}'.format(size))
    
    #Landy-Szalay version
    fig = plot_omega1(bins,corr_avg,dcorr_avg)
    fig.savefig(MC_directory+'/Size_{}_Fractal.pdf'.format(size),bbox_inches='tight')        
    plt.clf()
    plt.close(fig)


    #Fit exponential truncation to get truncation radius
    fig = plot_omega(bins,corr_avg,dcorr_avg,size=size)
    fig.savefig(MC_directory+'/Size_{}_Cutoff.pdf'.format(size),bbox_inches='tight')        
    plt.clf()
    plt.close(fig)

    return
    

#############################################################################################
# Some miscallaneous functions

def Fractal_Table_OneD(fractal,overwrite=False):
#check if pkl file exists        
    print("Fractal Dimension: {}".format(fractal))
    pkl_file = '../Toy_Models/Fractals/D_{}'.format(int(fractal*10))
    if(os.path.exists(pkl_file+'.pkl') and overwrite == False):
        print("Pickle summary exists reading from it.")
        pkl_obj = loadObj(pkl_file)
        xfractal,yfractal,bins,corr_lz,dcorr_lz = pkl_obj
    #Not computed,
    else:
        if(overwrite):
            print("Analysis being overwritten..")
        else:
            print("Analysis not been done for this fractal dimension. Performing...")
        
        #Set lower number of hierarchial levels for higher fractal dims to save time
        if(fractal>1.7):
            Analyse_Fractal(fractal,max_level=max_level_low)
        elif(fractal<1.0):
            Analyse_Fractal(fractal,max_level=max_level_high)
        else:
            Analyse_Fractal(fractal,max_level=max_level)
        pkl_obj = loadObj(pkl_file)
        xfractal,yfractal,bins,corr_lz,dcorr_lz = pkl_obj

    #Get Fractal dimension fitted to TPCF

    #Remove zeros
    separation_bins = (bins[1:]+bins[:-1])/2.
    indices = np.where(corr_lz>0.0)
    dcorr_lz = dcorr_lz[indices]
    separation_bins = separation_bins[indices]
    corr_lz = corr_lz[indices]

    #Fit fractal dimension to 1+Omega

    #Plot Positions
    plt.clf()
    plt.scatter(xfractal,yfractal,s=0.05,c='blue')
    plt.savefig('../Toy_Models/Fractals/Pos_D_{}'.format(int(fractal*10)),
                bbox_inches='tight')
    plt.close()
    
    plt.clf()
    plt.errorbar(separation_bins,1+corr_lz,yerr=dcorr_lz,
                 fmt='.-',lw=0.2)
    #Fit line to this
    popt,pcov = curve_fit(onepowerlaw_function,separation_bins,
                        np.log(1+corr_lz),sigma=dcorr_lz/(corr_lz))
    perr= np.sqrt(np.diag(pcov))

    plt.plot(separation_bins,np.exp(onepowerlaw_function(separation_bins,
                                                        popt[0],popt[1])),
            ls='-',lw=0.2,c='k',
             label=r'$\alpha = {:2.1f} \pm {:3.2f}$'.format(popt[1],perr[1]))



    plt.xlabel(r"$\Delta x \, (\mathrm{pixels})$")
    plt.ylabel(r"$1+ \omega_{\mathrm{LS}}\left(\theta \right)$")
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.savefig('../Toy_Models/Fractals/TPCF_D_{}'.format(int(fractal*10)),
                bbox_inches='tight')
    plt.clf()
    plt.close()

    D_tpcf = 2+popt[1]
    D_tpcf_err = perr[1]

    #Fit exponential truncation to get truncation radius
    bounds = ([-10.0,-3.0,np.min(separation_bins)],[10.0,0.0,np.max(separation_bins)*3.0])
    popt,pcov = curve_fit(linear_truncation,separation_bins,
                        np.log(corr_lz),sigma=dcorr_lz/corr_lz,bounds=bounds)
    perr= np.sqrt(np.diag(pcov))

    fig,axs = plt.subplots(ncols=1)

    axs.errorbar(separation_bins,corr_lz,yerr=dcorr_lz,
                 fmt='.-',lw=0.2)

    axs.plot(separation_bins,np.exp(linear_truncation(separation_bins,
                                                        popt[0],popt[1],popt[2])),
            ls='--',label=r'$\alpha = {:3.2f} \pm {:3.2f}$'.format(popt[1],perr[1]))

    size = 1.0
    axs.axvline(size/2.,color='#F9004A',lw=1.5,
                label=r'$r_{\mathrm{max}} = $'+r'${}$'.format(size),ls='--')
    axs.axvline(popt[2],ls='--',color='#F4271C',lw=1.5,
                  label=r'$r_c = {:3.2f} \pm {:3.2f}$'.format(popt[2],perr[2]))

    axs.set_xlabel(r"$\Delta x \, (\mathrm{pixels})$")
    axs.set_ylabel(r"$\omega_{\mathrm{LS}}\left(\theta \right)$")
    axs.set_xscale('log')
    axs.set_yscale('log')
    x0,x1 = axs.get_xlim()
    axs.set_xlim(x0,size)
    
    axs.legend(loc='best')
    plt.savefig('../Toy_Models/Fractals/TPCF_cutoff_{}'.format(int(fractal*10)),bbox_inches='tight')        
    plt.clf()
    plt.close(fig)

    #Save zeroth instance of boundary (i.e. full domain)
    L_cut = popt[2]
    L_cut_err = perr[2]
    L_max = 0.5
    D_omega = 2+popt[1]
    D_omega_err = perr[1]

    #Write to file
    file_str = np.column_stack([fractal,D_tpcf,D_tpcf_err,L_max,L_cut,L_cut_err,
        D_omega,D_omega_err])
    file = open(filename,'a')
    np.savetxt(file,file_str,delimiter='\t',fmt='%0.2f')
    file.write("\n")
    file.close()

    #Loop over boundaries
    MC_directory = os.path.abspath('../Toy_Models/Fractals/Monte_Carlo_Boundaries/Fractal_{}'.format(fractal))
            
    print("Looping over Monte Carlo Boundaries")
    for size in sizes:
        print("  Size : {}".format(size))
        
        try:
            TPCF_MC(xfractal,yfractal,size,MC_directory=MC_directory)            
        except Exception:
            print("Max attempts for this size attempted. Moving on..")
            continue

        pkl_obj = MC_directory + '/Size_{}'.format(int(size*10))
        bins,corr_avg,dcorr_avg = loadObj(pkl_obj)

        #Landy-Szalay version
        separation_bins = (bins[1:]+bins[:-1])/2.
        indices = np.where(corr_avg>0.0)
        dcorr_lz = dcorr_avg[indices]
        separation_bins = separation_bins[indices]
        corr_lz = corr_avg[indices]

        #Fit line to this for fractal dimension
        popt,pcov = curve_fit(onepowerlaw_function,separation_bins,
                    np.log(1+corr_lz),sigma=dcorr_lz/(1+corr_lz))
        perr= np.sqrt(np.diag(pcov))
        D_tpcf = 2+popt[1]
        D_tpcf_err = perr[1]
        
        #Fit linear+truncation model to Omega
        bounds = ([-10.0,-3.0,np.min(separation_bins)],[10.0,0.0,np.max(separation_bins)*3.0])
        popt,pcov = curve_fit(linear_truncation,separation_bins,
                            np.log(corr_lz),sigma=dcorr_lz/corr_lz,bounds=bounds)
        perr= np.sqrt(np.diag(pcov))
        
        #Save zeroth instance of boundary (i.e. full domain)
        L_cut = popt[2]
        L_cut_err = perr[2]
        L_max = size/2.
        D_omega = 2+popt[1]
        D_omega_err = perr[1]

        #Write to file
        file_str = np.column_stack([fractal,D_tpcf,D_tpcf_err,L_max,L_cut,L_cut_err,
            D_omega,D_omega_err])
        file = open(filename,'a')
        np.savetxt(file,file_str,delimiter='\t',fmt='%0.2f')
        file.write("\n")
        file.close()
    print("")
    print("")
    print("")
    print("##############################################################################")

def Fractal_Table(overwrite=False): 

    print("Running Fractal Structure Experiments..")
    print("##############################################################################")
    print("")

    #Create output table file    
    file = open(filename,'w')
    header = "#D_in\tD_tpcf\tD_tpcf_error\tL_max\tL_cutoff\tL_cutoff_error\tD_omega\tD_omega_err"
    file.write(header)
    file.write("\n")
    file.close()

    #Get positions and TPCF of fractals
    for fractal in fractal_dims:
        Fractal_Table_OneD(fractal,overwrite=overwrite)

def Fractal_Table_Parallel(overwrite=False): 

    print("Running Fractal Structure Experiments..")
    print("##############################################################################")
    print("")

    #Create output table file    
    file = open(filename,'w')
    header = "#D_in\tD_tpcf\tD_tpcf_error\tL_max\tL_cutoff\tL_cutoff_error\tD_omega\tD_omega_err"
    file.write(header)
    file.write("\n")
    file.close()

    #Get positions and TPCF of fractals
    results = Parallel(n_jobs=-1,prefer='processes',verbose=1)(map(delayed(Fractal_Table_OneD),
            fractal_dims,repeat(overwrite)))    

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

def fill_space(x,y,size=1):
    xlim = [np.min(x),np.max(x)]
    ylim = [np.min(y),np.max(y)]
    xr = xlim[0] + (xlim[1]-xlim[0])*np.random.random(size)
    yr = ylim[0] + (ylim[1]-ylim[0])*np.random.random(size)
    return xr,yr

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
                        random_state=None,random_type='uniform'):
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
    start = time.time()
    corr = two_point(data, bins, method=method, random_state=rng)
    end = time.time()
    #print("One cycle of TPCF done. Time taken: {} seconds".format(end-start))

    bootstraps = np.zeros((Nbootstrap, len(corr)))
    
    data_boot = np.zeros((Nbootstrap,n_samples,n_features))
    data_R = np.zeros((Nbootstrap,n_samples,n_features))
    for i in range(Nbootstrap):
        indices = rng.randint(0, n_samples, n_samples)
        data_boot[i] = data[indices,:]
        if(random_type=='uniform'):
            ra_R, dec_R = uniform_sphere((min(data[:,0]), max(data[:,0])),
                                         (min(data[:,1]), max(data[:,1])),
                                         len(data[:,0]))
            random_arr = np.asarray((ra_R, dec_R), order='F').T
            data_R[i] = random_arr
        elif(random_type == 'random'):
            ra_R,dec_R = fill_space(data[:,0],data[:,1],size=len(data[:,0]))
            random_arr = np.asarray((ra_R, dec_R), order='F').T
            data_R[i] = random_arr
        elif(random_type == 'default'):
            data_Rdefault = data.copy()
            for i in range(n_features - 1):
                rng.shuffle(data_Rdefault[:, i])
            data_R[i] = data_Rdefault
                   
        else:
            raise ValueError("Random type is undefined.")
    results = Parallel(n_jobs=-1,prefer='processes',verbose=0)(map(delayed(two_point),
        data_boot,repeat(bins),repeat(method),data_R,repeat(rng)))
    bootstraps = results
    bootstraps = np.asarray(bootstraps)
    corr = np.mean(bootstraps, 0)
    corr_err = np.std(bootstraps, 0, ddof=1)
    

    if return_bootstraps:
        return corr, corr_err, bootstraps
    else:
        return corr, corr_err

def plot_omega1(bins,corr,dcorr,fit='linear'):
    plt.clf()
    fig,axs = plt.subplots(ncols=1)
    separation_bins = (bins[1:]+bins[:-1])/2.
    indices = np.where(np.logical_and(corr>0.0,corr>dcorr))
    dcorr_lz = dcorr[indices]
    separation_bins = separation_bins[indices]
    corr_lz = corr[indices]
    
    
    axs.errorbar(separation_bins,1+corr_lz,yerr=dcorr_lz,
                 fmt='.-',lw=0.2)
    #Fit line to this
    popt,pcov = curve_fit(onepowerlaw_function,separation_bins,
                        np.log(1+corr_lz),sigma=dcorr_lz/corr_lz)
    perr= np.sqrt(np.diag(pcov))

    axs.plot(separation_bins,np.exp(onepowerlaw_function(separation_bins,
                                                        popt[0],popt[1])),
            ls='-',lw=0.2,c='k',
             label=r'$\alpha = {:2.1f} \pm {:3.2f}$'.format(popt[1],perr[1]))



    axs.set_xlabel(r"$\Delta x \, (\mathrm{pixels})$")
    axs.set_ylabel(r"$1+ \omega_{\mathrm{LS}}\left(\theta \right)$")
    axs.set_xscale('log')
    axs.set_yscale('log')
    axs.legend()
    return fig

def plot_omega(bins,corr,dcorr,size=None):
    plt.clf()
    fig,axs = plt.subplots(ncols=1)
    separation_bins = (bins[1:]+bins[:-1])/2.
    indices = np.where(corr>0.0)
    dcorr_lz = dcorr[indices]
    separation_bins = separation_bins[indices]
    corr_lz = corr[indices]
    bounds = ([-10.0,-3.0,np.min(separation_bins)],[10.0,0.0,np.max(separation_bins)*3.0])

    
    #Fit line to this
    popt,pcov = curve_fit(linear_truncation,separation_bins,
                        np.log(corr_lz),sigma=dcorr_lz/corr_lz,bounds=bounds)
    perr= np.sqrt(np.diag(pcov))

    axs.errorbar(separation_bins,corr_lz,yerr=dcorr_lz,
                 fmt='.-',lw=0.2)
    axs.plot(separation_bins,np.exp(linear_truncation(separation_bins,
                                                        popt[0],popt[1],popt[2])),
            ls='--',label=r'$\alpha = {:3.2f} \pm {:3.2f}$'.format(popt[1],perr[1]))



    axs.axvline(popt[2],ls='--',color='#F4271C',lw=1.5,
                  label=r'$r_c = {:3.2f} \pm {:3.2f}$'.format(popt[2],perr[2]))

    
    axs.set_xlabel(r"$\Delta x \, (\mathrm{pixels})$")
    axs.set_ylabel(r"$\omega_{\mathrm{LS}}\left(\theta \right)$")
    axs.set_xscale('log')
    axs.set_yscale('log')

    if(size):
        axs.axvline(size/2.,color='#F9004A',lw=1.5,
                label=r'$r_{\mathrm{max}}$'+r'$ = {}$'.format(size/2.),ls='--')
        x0,x1 = axs.get_xlim()
        axs.set_xlim(x0,size)

    axs.legend()
    return fig

def filter_stuff(bins,corr,dcorr):
    separation_bins = (bins[1:]+bins[:-1])/2.
    indices = np.where(corr>0.0)
    dcorr_lz = dcorr[indices]
    separation_bins = separation_bins[indices]
    corr_lz = corr[indices]
    return separation_bins,corr_lz,dcorr_lz




if __name__ == "__main__":
    # time the script
    start_time = timeit.default_timer()
    
    fractal_dim = 1.0
    Analyse_Fractal(fractal_dim,max_level=14)

    fractal_dim = 0.9
    Analyse_Fractal(fractal_dim,max_level=14)

    fractal_dim = 0.5
    Analyse_Fractal(fractal_dim,max_level=14)

    fractal_dim = 1.1
    Analyse_Fractal(fractal_dim,max_level=14)
    


    # time the script
    stop_time = timeit.default_timer()
    total_time = stop_time - start_time
    print("***************** time to finish = "+str(total_time)+"s *****************")
    print("Done")
