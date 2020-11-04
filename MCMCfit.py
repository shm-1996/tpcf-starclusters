from header import *
from Galaxy import *
from Plot_Class import *
import emcee
import corner


# beta_limits = [50.0, 300.0]

# distance = galaxy_class.distance*const.Parsec*1.e6
# #beta limits in arcsec
# beta_limits[0]= beta_limits[0]*const.Parsec/distance*u.radian.to(u.arcsec)
# beta_limits[1] = beta_limits[1]*const.Parsec/distance*u.radian.to(u.arcsec)

################################################################################################
#Piecewise power law functions
def model(params,data):
    #Parameters to fit
    A_1,alpha_1,alpha_2,beta = params
    #A_2 in terms of A_1 assuming continuity
    A_2 = A_1
    function = np.piecewise(data,[np.log(data)<beta],[lambda data :  A_1 + alpha_1 * np.log(data), 
        lambda data : A_2 + (alpha_1-alpha_2)*beta + alpha_2*np.log(data)])
    return function

def lnprior(params,distance):
#    A_1,A_2,alpha_1,alpha_2,beta = params
    A_1,alpha_1,alpha_2,beta = params
    beta_limits = [50.0, 300.0]
    beta_limits[0]= beta_limits[0]*const.Parsec/distance*u.radian.to(u.arcsec)
    beta_limits[1] = beta_limits[1]*const.Parsec/distance*u.radian.to(u.arcsec)
    if(-10<A_1<10 and 
       -5<alpha_1<0 and -5<alpha_2<0 and
      np.log(beta_limits[0])<beta<np.log(beta_limits[1])) :
        return 0
    else :
        return -np.inf
def lnprob(params, bins,data,data_error,distance):
    lp = lnprior(params,distance)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(params, bins,data,data_error) 

def lnlike(params,bins,data,data_error):
#    data_model = model_1(params,bins)
    data_model = model(params,bins)
    log_likelihood = -1/2. * np.sum(((data-data_model)/data_error)**2)
    return log_likelihood    

################################################################################################
#Single power law functions
def model_singlepower(params,data):
    A_1,alpha_1 = params
    function = A_1 + np.log(data)*alpha_1
    return function

def lnprior_singlepower(params):
#    A_1,A_2,alpha_1,alpha_2,beta = params
    A_1,alpha_1 = params
    if(-10<A_1<10 and 
       -5<alpha_1<0 ) :
        return 0
    else :
        return -np.inf

def lnprob_singlepower(params, bins,data,data_error):
    lp = lnprior_singlepower(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike_singlepower(params, bins,data,data_error)

def lnlike_singlepower(params,bins,data,data_error):
#    data_model = model_1(params,bins)
    data_model = model_singlepower(params,bins)
    log_likelihood = -1/2. * np.sum(((data-data_model)/data_error)**2)
    return log_likelihood


################################################################################################
#Fitting routines

def fit_MCMC_galaxy(galaxy_name,method='masked_radial',function='piecewise'):
    """
    Fit an MCMC to a galaxy with given name. 
    Parameters: 
        galaxy_name : str 
            Name of galaxy for which to fit
        method : str
            Method of Random distribution
        function: str 
            Piecewise or single PL fits to the TPCF


    """
    galaxy_name = galaxy_name.upper()
    galaxy_class = Galaxy(galaxy_name,verbose=False)
    #Save in subdirectories for below two methods
    if(method == 'masked'):
        galaxy_class.outdir += '/Masked/'
    elif(method == 'uniform'):
        galaxy_class.outdir += '/Uniform/'
    elif(method == 'masked_radial'):
        galaxy_class.outdir += '/Masked_Radial/'
    else:
        raise myError("Method not recognised.")

    galaxy_class = loadObj(galaxy_class.outdir+
            '{}_summary'.format(galaxy_class.name))
    if(function == 'piecewise'):
        fit_MCMC(galaxy_class,save=True)
    else :
        print("Fitting Single PL fit in MCMC.")
        fit_SinglePLMCMC(galaxy_class,save=True)

def fit_MCMC(galaxy_class,save=False):
    """
    Fit an MCMC to the TPCF of a galaxy and create diagnostic plots after. 

    Parameters:
        galaxy_class: Class Galaxy
            Galaxy class object 
        save : Boolean
            Flag to save the figure, else just plot. 
    """


    initial_guess = galaxy_class.fit_values
    beta_limits = [50.0, 300.0]

    distance = galaxy_class.distance*const.Parsec*1.e6
    #beta limits in arcsec
    beta_limits[0]= beta_limits[0]*const.Parsec/distance*u.radian.to(u.arcsec)
    beta_limits[1] = beta_limits[1]*const.Parsec/distance*u.radian.to(u.arcsec)
    initial_guess[3] = np.log((beta_limits[0]+beta_limits[1])/2.)
    nwalkers = 500
    ndim = len(initial_guess)
    step = 1.e-7
    p0 = [np.array(initial_guess) + step * np.random.randn(ndim) for i in range(nwalkers)]


    separation_bins = (galaxy_class.bins[1:]+galaxy_class.bins[:-1])/2
    separation_bins*=(1./arcsec_to_degree)
    corr_fit = galaxy_class.corr[np.where(galaxy_class.corr>0.0)].astype(np.float)
    dcorr_fit = galaxy_class.dcorr[np.where(galaxy_class.corr>0.0)].astype(np.float)
    separation_bins = separation_bins[np.where(galaxy_class.corr>0.0)].astype(np.float)
    distance = galaxy_class.distance*const.Parsec*1.e6
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, 
                                    args=(separation_bins,np.log(corr_fit),
                                         np.log(dcorr_fit),distance))    

    print("Running burn-in...")
    p0, _, _ = sampler.run_mcmc(p0, 100)
    sampler.reset()

    print("Running production...")
    pos, prob, state = sampler.run_mcmc(p0, 10000,progress=True)

    #Get samples
    samples = sampler.flatchain

    #Plot corner plot
    labels = [r'$A_1$',r'$\alpha_1$',r'$\alpha_2$',r'$\beta$']
    fig = corner.corner(samples,show_titles=True,labels=labels,
        plot_datapoints=True,quantiles=[0.16, 0.5, 0.84])
    if(save):
        fig.savefig(galaxy_class.outdir+'MCMC_Posterior_Dist.pdf',bbox_inches='tight')
        plt.close()
    else:
        fig.show()

    #Plot Convergence check
    chain = sampler.get_chain(discard=200,thin=40)
    fig, axes = plt.subplots(4, figsize=(10, 7), sharex=True)
    for i in range(ndim):
        ax = axes[i]
        ax.plot(chain[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(chain))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)

    axes[-1].set_xlabel("step number")
    if(save):
        fig.savefig(galaxy_class.outdir+'MCMC_Convergence.pdf',bbox_inches='tight')
        plt.close()
    else:
        fig.show()

    #Overplot sample fits on data
    fig,axs = plt.subplots(ncols=1)
    ax2 = axs.secondary_xaxis("top",functions=(sep_to_pc,pc_to_sep))

    separation_bins = (galaxy_class.bins[1:]+galaxy_class.bins[:-1])/2
    separation_bins*=(1./arcsec_to_degree)
    chain = sampler.get_chain(discard=200,thin=40,flat=True)
    inds = np.random.randint(len(chain), size=50)
    
    for ind in inds:
        sample = chain[ind]
        axs.plot(separation_bins,np.exp(model(sample,separation_bins)),
             alpha=0.1)
    axs.errorbar(separation_bins,galaxy_class.corr,yerr=galaxy_class.dcorr,
                fmt='.-')

    axs.set_xlabel(r"$\theta \, \left(\mathrm{arcsec} \right)$")
    axs.set_ylabel(r"$\omega_{\mathrm{LS}}\left(\theta \right)$")
    axs.set_xscale('log')
    axs.set_yscale('log')
    axs.callbacks.connect("xlim_changed", axs_to_parsec)
    axs.legend()
    ax2.set_xlabel(r'$\delta x \, \left( \mathrm{pc} \right) $')

    if(save):
        plt.savefig(galaxy_class.outdir+'Sample_fits_MCMC')
        plt.close()
    else:
        plt.show()

    galaxy_class.fit_values = samples[np.argmax(sampler.flatlnprobability)]
    galaxy_class.fit_errors = sampler.flatchain.std(axis=0)
    pl = myPlot(galaxy_class)
    pl.plot_TPCF(save=save,filename=galaxy_class.outdir+'/TPCF_MCMC.pdf')
    saveObj(sampler,galaxy_class.outdir+'/MCMC_sampler')


def fit_SinglePLMCMC(galaxy_class,save=False):
    """
    Fit an MCMC to the TPCF of a galaxy and create diagnostic plots after. 

    Parameters:
        galaxy_class: Class Galaxy
            Galaxy class object 
        save : Boolean
            Flag to save the figure, else just plot. 
    """


    initial_guess = 5.0,-1.0
    
    nwalkers = 500
    ndim = len(initial_guess)
    step = 1.e-7
    p0 = [np.array(initial_guess) + step * np.random.randn(ndim) for i in range(nwalkers)]


    separation_bins = (galaxy_class.bins[1:]+galaxy_class.bins[:-1])/2
    separation_bins*=(1./arcsec_to_degree)
    corr_fit = galaxy_class.corr[np.where(galaxy_class.corr>0.0)].astype(np.float)
    dcorr_fit = galaxy_class.dcorr[np.where(galaxy_class.corr>0.0)].astype(np.float)
    separation_bins = separation_bins[np.where(galaxy_class.corr>0.0)].astype(np.float)
    distance = galaxy_class.distance*const.Parsec*1.e6
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_singlepower, 
                                    args=(separation_bins,np.log(corr_fit),
                                         np.log(dcorr_fit)))    

    if(not os.path.exists(galaxy_class.outdir+'/SinglePL_MCMC')):
        os.makedirs(galaxy_class.outdir+'/SinglePL_MCMC')

    print("Running burn-in...")
    p0, _, _ = sampler.run_mcmc(p0, 100)
    sampler.reset()

    print("Running production...")
    pos, prob, state = sampler.run_mcmc(p0, 10000,progress=True)

    #Get samples
    samples = sampler.flatchain

    #Plot corner plot
    labels = [r'$A_1$',r'$\alpha_1$']
    fig = corner.corner(samples,show_titles=True,labels=labels,
        plot_datapoints=True,quantiles=[0.16, 0.5, 0.84])
    if(save):
        fig.savefig(galaxy_class.outdir+'/SinglePL_MCMC/MCMC_Posterior_Dist.pdf',bbox_inches='tight')
        plt.close()
    else:
        fig.show()

    #Plot Convergence check
    chain = sampler.get_chain(discard=200,thin=40)
    fig, axes = plt.subplots(2, figsize=(8, 7), sharex=True)
    for i in range(ndim):
        ax = axes[i]
        ax.plot(chain[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(chain))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)

    axes[-1].set_xlabel("step number")
    if(save):
        fig.savefig(galaxy_class.outdir+'/SinglePL_MCMC/MCMC_Convergence.pdf',bbox_inches='tight')
        plt.close()
    else:
        fig.show()

    #Overplot sample fits on data
    fig,axs = plt.subplots(ncols=1)
    ax2 = axs.secondary_xaxis("top",functions=(sep_to_pc,pc_to_sep))

    separation_bins = (galaxy_class.bins[1:]+galaxy_class.bins[:-1])/2
    separation_bins*=(1./arcsec_to_degree)
    chain = sampler.get_chain(discard=200,thin=40,flat=True)
    inds = np.random.randint(len(chain), size=50)
    
    for ind in inds:
        sample = chain[ind]
        axs.plot(separation_bins,np.exp(model_singlepower(sample,separation_bins)),
             alpha=0.1)
    axs.errorbar(separation_bins,galaxy_class.corr,yerr=galaxy_class.dcorr,
                fmt='.-')

    axs.set_xlabel(r"$\theta \, \left(\mathrm{arcsec} \right)$")
    axs.set_ylabel(r"$\omega_{\mathrm{LS}}\left(\theta \right)$")
    axs.set_xscale('log')
    axs.set_yscale('log')
    axs.callbacks.connect("xlim_changed", axs_to_parsec)
    axs.legend()
    ax2.set_xlabel(r'$\delta x \, \left( \mathrm{pc} \right) $')

    if(save):
        plt.savefig(galaxy_class.outdir+'/SinglePL_MCMC/Sample_fits_MCMC')
        plt.close()
    else:
        plt.show()

    galaxy_class.fit_values = samples[np.argmax(sampler.flatlnprobability)]
    galaxy_class.fit_errors = sampler.flatchain.std(axis=0)
    pl = myPlot(galaxy_class)
    pl.plot_TPCF(save=save,filename=galaxy_class.outdir+'/SinglePL_MCMC/TPCF_MCMC.pdf',
        function='singlepl')
    saveObj(sampler,galaxy_class.outdir+'/SinglePL_MCMC/MCMC_sampler')


def sep_to_pc(sep) :

    """Angular separation in arcsec to parsec

    Parameters
    ----------
    sep : float
        separation in angular 
    Returns
    -------
    pc : float
         The separation in parsec for given angular separation
    """

    # Distance to NGC 628 = 9.9 Mpc : Olivares et al 2010
    distance_628 = 9.9*const.Parsec*1.e6  
    arcsec_to_pc = u.arcsec.to(u.radian)*distance_628/(const.Parsec)
    return sep*arcsec_to_pc

def pc_to_sep(pc) :
    """Linear separation in parsec to arcsec

    Parameters
    ----------
    pc : float
        separation in linear
    Returns
    -------
    pc : float
         The separation in parsec for given angular separation
    """
    # Distance to NGC 628 = 9.9 Mpc : Olivares et al 2010
    distance_628 = 9.9*const.Parsec*1.e6  
    pc_to_radian = pc*const.Parsec/distance_628
    return pc_to_radian*u.radian.to(u.arcsec)

def axs_to_parsec(axs) :
    """Draw twin axes corresponding to angular limits

    Parameters
    ----------
    axs : Matplotlib axis instance
          Matplotlib axis on which to plot
    Returns
    -------
    None
    """
    xmin,xmax = axs.get_xlim()
    ax2.set_ylim(sep_to_pc(xmin),sep_to_pc(xmax))
    ax2.figure.canvas.draw()


