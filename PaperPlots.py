"""
Routines to create the plots of the paper : Menon et al 2020. 

AUTHOR:
Shyam Harimohan Menon (2020)
"""
from header import *
from Plot_Class import *

def Combined_TPCF(save=False,outdir='../Results/',indir=None,method='masked_radial'):
    """
    Plot the TPCF of all galaxies as a combined figure.

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
        plot_class = myPlot(galaxy_class)

        #Secondary axis
        separation_bins = galaxy_class.bin_centres*(1./arcsec_to_degree)
        
        try:
            
            axs[i,j].set_yscale('log')
            axs[i,j].errorbar(separation_bins,galaxy_class.corr,yerr=galaxy_class.dcorr,
                fmt='.-')
            ax2 = axs[i,j].secondary_xaxis("top",functions=(plot_class.sep_to_pc,plot_class.pc_to_sep))
            
            #Fit plot
            break_theta = np.exp(galaxy_class.fit_values[3])
            break_theta_error = np.exp(galaxy_class.fit_errors[3])
            axs[i,j].plot(separation_bins,np.exp(linear_function(separation_bins,galaxy_class.fit_values[0],
                galaxy_class.fit_values[1],galaxy_class.fit_values[2],galaxy_class.fit_values[3])),
                ls='--',label='fit')
            axs[i,j].plot(separation_bins,galaxy_class.corr,lw=0.0,
                label=r'$\alpha_1 = {:2.1f} \pm {:2.1f}$'.format(galaxy_class.fit_values[1],galaxy_class.fit_errors[1]))
            axs[i,j].plot(separation_bins,galaxy_class.corr,lw=0.0,
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

        except:
            print("Cannot plot TPCF") 

        

        #Get position of subplot
        j +=1
        if(j==4):
            j = 0
            i +=1


    if(save):
        filename = outdir+'Combined_TPCF_Masked.pdf'
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else :
        plt.show()

def Combined_RGBImages(save=False,outdir='../Results/',indir=None,method='masked_radial'):
    #Create figure and axs instance
    fig = plt.figure(figsize=(16,12),constrained_layout=True)

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
        plot_class = myPlot(galaxy_class)

        #WCS info
        hdu = fits.open(galaxy_class.fits_file)[0]
        wcs = WCS(hdu.header)        

        axs = plt.subplot2grid((4,3),(i,j),projection=wcs)
        with np.errstate(divide='ignore', invalid='ignore'):
            im = axs.imshow(np.log10(hdu.data),vmin=-2.0)

        if(i==3):
            axs.set_xlabel(r"$\mathrm{Right \; Ascension \; (J2000)}$",fontsize=20)
        if(j==0):
            axs.set_ylabel(r"$\mathrm{Declination \; (J2000)}$",labelpad=-0.8,
                       fontsize=20)
        
        axs.text(0.1,0.1,r'$\mathrm{NGC}$'+' '+r'${}$'.format(galaxy_name.split('_')[1]),
            transform=axs.transAxes)

        #Get position of subplot
        #Get position of subplot
        j +=1
        if(j==3):
            j = 0
            i +=1


    if(save):
        filename = outdir+'Combined_RGBImages.pdf'
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else :
        plt.show()


def Combined_Clusters(save=False,outdir='../Results/',indir=None,method='masked_radial'):
    #Create figure and axs instance
    fig = plt.figure(figsize=(16,20),constrained_layout=True)

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
        plot_class = myPlot(galaxy_class)

        #WCS info
        hdu = fits.open(galaxy_class.fits_file)[0]
        wcs = WCS(hdu.header)        

        axs = plt.subplot2grid((4,3),(i,j),projection=wcs,fig=fig)
        im = axs.scatter(galaxy_class.ra,galaxy_class.dec)

        if(i==3):
            axs.set_xlabel(r"$\mathrm{Right \; Ascension \; (J2000)}$",fontsize=20)
        if(j==0):
            axs.set_ylabel(r"$\mathrm{Declination \; (J2000)}$",labelpad=-0.8,
                       fontsize=20)
        
        axs.text(0.1,0.1,r'$\mathrm{NGC}$'+' '+r'${}$'.format(galaxy_name.split('_')[1]),
            transform=axs.transAxes)

        #Get position of subplot
        #Get position of subplot
        j +=1
        if(j==3):
            j = 0
            i +=1


    if(save):
        filename = outdir+'Combined_Clusters.pdf'
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else :
        plt.show()

if __name__ == "__main__":
    print("Preparing plots of paper.")
    Combined_TPCF(save=True,method='masked')

