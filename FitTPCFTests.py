from header import *
from Galaxy import *
from Plot_Class import *
from MCMCfit import *


def plot_MCMCfitsall(save=False,outdir='../Results/',indir=None,method='masked_radial',function='piecewise'):
    """
    Read MCMC samplers and plot the TPCF with best fit power law for all galaxies. 

    Parameters: 
        save: 
        Flag to save the plot
        outdir: 
            Output directory in which to store plot. Default is results directory.
        indir :
            Input directory from which to read the results. Default is results directory.
        method : 
            Method for which TPCF's have been computed.
        function:
            Function fitted using MCMC: piecewise or singlepl
    Returns:
        None


    """

    print("Plotting Combined TPCF plot with MCMC fits using a {} function.".format(function))

    #Create figure and axs instance
    fig,axs = plt.subplots(nrows=3,ncols=4,figsize=(16,12))

    if(indir == None):
        indir = os.path.abspath('../Results/Galaxies/')+'/'
        method_dir = ''
    else :
        method_dir = None

    #Directories
    indir = os.path.abspath(indir)+'/'
    outir = os.path.abspath(outdir)+'/'

    if(method_dir is not None):
        if(method.upper() == 'MASKED'):
            method_dir = 'Masked/'
        elif(method.upper() == 'UNIFORM'):
            method_dir = 'Uniform/'
        elif(method.upper() == 'MASKED_RADIAL'):
            method_dir = 'Masked_Radial/'
        else:
            raise myError("Method not recognised.")

    if(function not in ['piecewise','singlepl']):
        raise ValueError("Function should be either piecewise or singlepl.")

    i,j = 0,0            
    #Loop through the galaxies
    for galaxy_name in list_of_galaxies:
        if(method_dir == None):
            galaxy_dir = indir+galaxy_name+'/'            
        else :
            galaxy_dir = indir+galaxy_name+'/'+method_dir

        galaxy_class = loadObj(galaxy_dir+galaxy_name+'_summary')
        if(function == 'piecewise'):
            
            sampler = loadObj(galaxy_dir+'/PiecePL_MCMC/'+'MCMC_sampler')
        else :
            sampler = loadObj(galaxy_dir+'/SinglePL_MCMC/'+'MCMC_sampler')
            
        
        samples = sampler.flatchain
        galaxy_class.fit_values = samples[np.argmax(sampler.flatlnprobability)]
        galaxy_class.fit_errors = samples.std(axis=0)
        plot_class = myPlot(galaxy_class)

        #Secondary axis
        separation_bins = galaxy_class.bin_centres*(1./arcsec_to_degree)
        separation_bins = separation_bins.astype(np.float)
        plot_points = np.linspace(np.min(separation_bins),np.max(separation_bins),1000)

        indices = np.where(galaxy_class.corr>0.0)
        corr_fit = galaxy_class.corr[indices].astype(np.float)
        dcorr_fit = galaxy_class.dcorr[indices].astype(np.float)
        separation_bins = separation_bins[indices].astype(np.float)
        
        try:
            
            axs[i,j].set_yscale('log')
            axs[i,j].errorbar(separation_bins,corr_fit,yerr=dcorr_fit,
                fmt='.-')
            ax2 = axs[i,j].secondary_xaxis("top",functions=(plot_class.sep_to_pc,plot_class.pc_to_sep))
            
            #Fit plot
            if(function == 'piecewise'):
                break_theta = np.exp(galaxy_class.fit_values[3])
                break_theta_error = np.exp(galaxy_class.fit_errors[3])
                axs[i,j].plot(plot_points,np.exp(linear_function(plot_points,galaxy_class.fit_values[0],
                    galaxy_class.fit_values[1],galaxy_class.fit_values[2],galaxy_class.fit_values[3])),
                    ls='--',label='fit')
                axs[i,j].plot(separation_bins,corr_fit,lw=0.0,
                    label=r'$\alpha_1 = {:2.1f} \pm {:3.2f}$'.format(galaxy_class.fit_values[1],galaxy_class.fit_errors[1]))
                axs[i,j].plot(separation_bins,corr_fit,lw=0.0,
                    label=r'$\alpha_2 = {:2.1f} \pm {:3.2f}$'.format(galaxy_class.fit_values[2],galaxy_class.fit_errors[2]))
                axs[i,j].axvline(break_theta,ls=':',label=r'$\beta = {:2.1f} \pm {:2.1f}$'.format(break_theta,
                    break_theta_error))
            else:

                axs[i,j].plot(plot_points,np.exp(onepowerlaw_function(plot_points,galaxy_class.fit_values[0],
                galaxy_class.fit_values[1])),
                ls='--',label='fit')

                axs[i,j].plot(separation_bins,corr_fit,lw=0.0,
                    label=r'$\alpha_1 = {:2.1f} \pm {:3.2f}$'.format(galaxy_class.fit_values[1],
                        galaxy_class.fit_errors[1]))

            #Plot stuff
            #X-labels only on bottom row
            if(i==2):
                axs[i,j].set_xlabel(r"$\theta \, \left(\mathrm{arcsec} \right)$")
            #Y-labels only on left column
            if(j == 0):
                axs[i,j].set_ylabel(r"$\omega_{\mathrm{LS}}\left(\theta \right)$")
            axs[i,j].set_xscale('log')
            axs[i,j].callbacks.connect("xlim_changed", plot_class.axs_to_parsec)
            axs[i,j].legend()

            #Secondary axis label only for top row
            if(i==0):
                ax2.set_xlabel(r'$\delta x \, \left( \mathrm{pc} \right) $')

            
            axs[i,j].text(0.1,0.1,r'$\mathrm{NGC}$'+' '+r'${}$'.format(galaxy_name.split('_')[1]),
                transform=axs[i,j].transAxes)

        except:
            print("Cannot plot TPCF") 

        

        #Get position of subplot
        j +=1
        if(j==4):
            j = 0
            i +=1


    if(save):
        if(function == 'piecewise'):
            filename = outdir+'Combined_TPCF_MCMC_Piecewise.pdf'
        else :
            filename = outdir+'Combined_TPCF_MCMC_SinglePL.pdf'
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else :
        plt.show()



def Test_CombinedFit(save=False,outdir='../Results/',indir=None,method='masked_radial',
    function='piecewise'):
    """
    Plot the TPCF of all galaxies along with fitting as a combined figure.

    Parameters:
        save: 
            Flag to save the plot
        outdir: 
            Output directory in which to store plot. Default is results directory.
        indir :
            Input directory from which to read the results. Default is results directory.
        method : 
            Method for which TPCF's have been computed.
    Returns:
        None


    """
    
    if(function not in ['piecewise','smooth','singlepl']):
        raise ValueError("This functional form does not exist.")


    #Create figure and axs instance
    fig,axs = plt.subplots(nrows=3,ncols=4,figsize=(16,12))

    if(indir == None):
        indir = os.path.abspath('../Results/Galaxies/')+'/'
        method_dir = ''
    else :
        method_dir = None

    #Directories
    indir = os.path.abspath(indir)+'/'
    outir = os.path.abspath(outdir)+'/'

    if(method_dir is not None):
        if(method.upper() == 'MASKED'):
            method_dir = 'Masked/'
        elif(method.upper() == 'UNIFORM'):
            method_dir = 'Uniform/'
        elif(method.upper() == 'MASKED_RADIAL'):
            method_dir = 'Masked_Radial/'
        else:
            raise myError("Method not recognised.")


    i,j = 0,0            
    #Loop through the galaxies
    for galaxy_name in list_of_galaxies:
        if(method_dir == None):
            galaxy_dir = indir+galaxy_name+'/'            
        else :
            galaxy_dir = indir+galaxy_name+'/'+method_dir

        galaxy_class = loadObj(galaxy_dir+galaxy_name+'_summary')
        if(function == 'singlepl'):
            sampler = loadObj(galaxy_dir+'/SinglePL_MCMC/MCMC_sampler')
            samples = sampler.flatchain
            galaxy_class.fit_values = samples[np.argmax(sampler.flatlnprobability)]
            galaxy_class.fit_errors = samples.std(axis=0)
        else:
            galaxy_class.fit_power_law(use_bounds=False,function=function)

        plot_class = myPlot(galaxy_class)

        #Secondary axis
        separation_bins = galaxy_class.bin_centres*(1./arcsec_to_degree)
        plot_points = np.linspace(np.min(separation_bins),np.max(separation_bins),1000)

        indices = np.where(galaxy_class.corr>0.0)
        corr_fit = galaxy_class.corr[indices].astype(np.float)
        dcorr_fit = galaxy_class.dcorr[indices].astype(np.float)
        separation_bins = separation_bins[indices].astype(np.float)
        

        
            
        axs[i,j].set_yscale('log')
        axs[i,j].errorbar(separation_bins,corr_fit,yerr=dcorr_fit,
            fmt='.-')
        ax2 = axs[i,j].secondary_xaxis("top",functions=(plot_class.sep_to_pc,plot_class.pc_to_sep))
        
        #Fit plot
        if(function == 'piecewise'):
            break_theta = np.exp(galaxy_class.fit_values[3])
            break_theta_error = np.exp(galaxy_class.fit_errors[3])
            axs[i,j].plot(plot_points,np.exp(linear_function(plot_points,galaxy_class.fit_values[0],
                galaxy_class.fit_values[1],galaxy_class.fit_values[2],galaxy_class.fit_values[3])),
                ls='--',label='fit')
        elif(function == 'smooth'):
            axs[i,j].plot(plot_points,smooth_function(plot_points,galaxy_class.fit_values[0],
                galaxy_class.fit_values[1],galaxy_class.fit_values[2],galaxy_class.fit_values[3]),
                ls='--',label='fit')
            break_theta = galaxy_class.fit_values[3]
            break_theta_error = galaxy_class.fit_errors[3]

        elif(function == 'singlepl'):
            axs[i,j].plot(plot_points,np.exp(onepowerlaw_function(plot_points,galaxy_class.fit_values[0],
                galaxy_class.fit_values[1])),
                ls='--',label='fit')

        axs[i,j].plot(separation_bins,corr_fit,lw=0.0,
            label=r'$\alpha_1 = {:2.1f} \pm {:2.1f}$'.format(galaxy_class.fit_values[1],galaxy_class.fit_errors[1]))
        
        if(function in ['piecewise','smooth']):
            axs[i,j].plot(separation_bins,corr_fit,lw=0.0,
                label=r'$\alpha_2 = {:2.1f} \pm {:2.1f}$'.format(galaxy_class.fit_values[2],galaxy_class.fit_errors[2]))
            axs[i,j].axvline(break_theta,ls=':',label=r'$\beta = {:2.1f} \pm {:2.1f}$'.format(break_theta,
                break_theta_error))

        #Plot stuff
        #X-labels only on bottom row
        if(i==2):
            axs[i,j].set_xlabel(r"$\theta \, \left(\mathrm{arcsec} \right)$")
        #Y-labels only on left column
        if(j == 0):
            axs[i,j].set_ylabel(r"$\omega_{\mathrm{LS}}\left(\theta \right)$")
        axs[i,j].set_xscale('log')
        axs[i,j].callbacks.connect("xlim_changed", plot_class.axs_to_parsec)
        axs[i,j].legend()

        #Secondary axis label only for top row
        if(i==0):
            ax2.set_xlabel(r'$\delta x \, \left( \mathrm{pc} \right) $')

        
        axs[i,j].text(0.1,0.1,r'$\mathrm{NGC}$'+' '+r'${}$'.format(galaxy_name.split('_')[1]),
            transform=axs[i,j].transAxes)

        

        

        #Get position of subplot
        j +=1
        if(j==4):
            j = 0
            i +=1


    if(save):
        if(function == 'singlepl'):
            filename = outdir+'Combined_TPCF_SinglePL.pdf'
        elif(function == 'smooth'):
            filename = outdir+'Combined_TPCF_Smooth.pdf'    
        elif(function == 'piecewise'):
            filename = outdir+'Combined_TPCF_Piecewise.pdf'    
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else :
        plt.show()


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=
        'Command Line Inputs for tpcf-starclusters. All inputs optional. ')
    ap.add_argument('-method',action='store',type=str,default='masked',
        help='Method to prepare the random catalog: "Uniform","Masked"' +
        '" Masked (default)" ')
    ap.add_argument('-galaxy',action='store',type=str,default=None,
        help = 'Galaxy for which tpcf to be computed. By default done for all.')
    ap.add_argument('-function',action='store',type=str,default='piecewise',
        help='Function to use to fit.')

    args = vars(ap.parse_args())
    if(args['function'] not in ['piecewise','singlepl']):
        raise ValueError("Wrong function type.")

    method = args['method'].lower()
    if(method not in ['uniform','masked','masked_radial']):
        raise ValueError("This method does not exist. Allowed values are "+
            "'Uniform', 'Masked', and 'Masked_Radial'.")

    galaxy_input = args['galaxy'].upper()
    if(galaxy_input == 'ALL'):

        for galaxy_name in list_of_galaxies:

            print("Performing MCMC for galaxy {}".format(galaxy_name))
            fit_MCMC_galaxy(galaxy_name,method=method,function=args['function'])
            plot_MCMCfitsall(save=True,method=method,function=args['function'])

    else :
        fit_MCMC_galaxy(galaxy_input,method=method,function=args['function'])


