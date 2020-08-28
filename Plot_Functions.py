from header import *
# from Galaxy_TPCF import linear_function
def Plot_TPCF(galaxy,output=None,save=False) :

    separation_bins = (galaxy.bins[1:]+galaxy.bins[:-1])/2
    separation_bins*=(1./arcsec_to_degree)


    fig,axs = plt.subplots(ncols=1)
    ax2 = axs.secondary_xaxis("top",functions=(sep_to_pc,pc_to_sep))

    axs.errorbar(separation_bins,corr,yerr=dcorr,fmt='.-',label='NGC 628')
    
    #TODO: Fitting
    #Fitting and plot fit
    popt,pcov = fit_power_law(corr,dcorr,bins)
    A_1,A_2,alpha_1,alpha_2,beta = popt
    # perr = np.sqrt(np.diag(pcov))
    #error_A1,error_A2,error_alpha1,error_alpha2,error_beta = perr
    axs.plot(separation_bins,np.exp(linear_function(separation_bins,A_1,A_2,alpha_1,
        alpha_2,beta)),ls=':',label='fit')

    axs.set_xlabel(r"$\theta \, \left(\mathrm{arcsec} \right)$")
    axs.set_ylabel(r"$\omega_{\mathrm{LS}}\left(\theta \right)$")
    axs.set_xscale('log')
    axs.set_yscale('log')
    axs.callbacks.connect("xlim_changed", axs_to_parsec)
    axs.legend()
    ax2.set_xlabel(r'$\delta x \, \left( \mathrm{pc} \right) $')

    

    if(save) :
        filename = output+galaxy
        filename = filename +'_TPCF'
        plt.savefig(filename)
        plt.clf()
        plt.close()
    else :
        plt.show()
    return



### TODO : Get distance to galaxy from data
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