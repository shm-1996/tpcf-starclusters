"""
Routines to create the plots of the paper : Menon et al 2020. 

AUTHOR:
Shyam Harimohan Menon (2020)
"""
from header import *
from Plot_Class import *
from CreateTables import compare_AIC

#Axes limits in parsec
global_axes_limits = [8,1.e4]

def Combined_TPCF(save=False,outdir='../Results/',indir=None,method='masked'):
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
    print("Plotting TPCF for all galaxies")

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
        separation_bins = separation_bins.astype(np.float)
        plot_points = np.linspace(np.min(separation_bins),np.max(separation_bins),1000)

        
        indices = np.where(galaxy_class.corr>0.0)
        corr_fit = galaxy_class.corr[indices].astype(np.float)
        dcorr_fit = galaxy_class.dcorr[indices].astype(np.float)
        separation_bins = separation_bins[indices].astype(np.float)
            
        axs[i,j].errorbar(separation_bins,corr_fit,yerr=dcorr_fit,
            fmt='.-')

        #Set X-Axis Limits
        distance = galaxy_class.distance*const.Parsec*1.e6
        axs_limits = [0.0,0.0]
        axs_limits[0] =  global_axes_limits[0]*const.Parsec/distance*u.radian.to(u.arcsec)
        axs_limits[1] =  global_axes_limits[1]*const.Parsec/distance*u.radian.to(u.arcsec)
        axs[i,j].set_xlim(axs_limits[0],axs_limits[1])

        ax2 = axs[i,j].secondary_xaxis("top",functions=(plot_class.sep_to_pc,plot_class.pc_to_sep))
        
        #Fit plot
        break_theta = np.exp(galaxy_class.fit_values[3])
        break_theta_error = np.exp(galaxy_class.fit_errors[3])
        axs[i,j].plot(plot_points,np.exp(linear_function(plot_points,galaxy_class.fit_values[0],
            galaxy_class.fit_values[1],galaxy_class.fit_values[2],galaxy_class.fit_values[3])),
            ls='--',label='fit')
        axs[i,j].plot(separation_bins,corr_fit,lw=0.0,
            label=r'$\alpha_1 = {:2.1f} \pm {:2.1f}$'.format(galaxy_class.fit_values[1],galaxy_class.fit_errors[1]))
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
        axs[i,j].set_yscale('log')
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
        filename = outdir+'Combined_TPCF_Masked.pdf'
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else :
        plt.show()


def Combined_TPCF_CompareAges(save=False,outdir='../Results/',indir=None,method='masked'):
    """
    Plot the TPCF of clusters based on age cuts.

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
    print("Plotting TPCFs with Age cuts.")

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

        separation_bins = galaxy_class.bin_centres*(1./arcsec_to_degree)
        separation_bins = separation_bins.astype(np.float)
        plot_points = np.linspace(np.min(separation_bins),np.max(separation_bins),1000)

        #Set axes
        axs[i,j].set_yscale('log')            
        ax2 = axs[i,j].secondary_xaxis("top",functions=(plot_class.sep_to_pc,plot_class.pc_to_sep))
        
        #General TPCF
        indices = np.where(galaxy_class.corr>0.0)
        corr_fit = galaxy_class.corr[indices].astype(np.float)
        dcorr_fit = galaxy_class.dcorr[indices].astype(np.float)
        separation_bins = separation_bins[indices].astype(np.float)
            
        axs[i,j].errorbar(separation_bins,corr_fit,yerr=dcorr_fit,
            fmt='.-',color='#AAF54590',label=r'$\mathrm{All} \, T$')

        ########################################
        ages = galaxy_class.get_cluster_ages()
        galaxy_class.get_ra_dec()
        
        #Clusters < 10Myr
        galaxy_class.ra = galaxy_class.ra[np.where(ages <= 1.e7)]
        galaxy_class.dec = galaxy_class.dec[np.where(ages <= 1.e7)] 
        

        #Young Clusters TPCF            
        corr,dcorr,bootstraps = bootstrap_two_point_angular(galaxy_class,
                        method='landy-szalay',Nbootstraps=100,
                        random_method=method)
        galaxy_class.corr = corr 
        galaxy_class.dcorr = dcorr

        #Isolate non-zero correlation points
        indices = np.where(galaxy_class.corr>0.0)
        corr_fit = galaxy_class.corr[indices].astype(np.float)
        dcorr_fit = galaxy_class.dcorr[indices].astype(np.float)
        separation_bins = galaxy_class.bin_centres*(1./arcsec_to_degree)
        separation_bins = separation_bins.astype(np.float)
        separation_bins = separation_bins[indices].astype(np.float)
        separation_bins = separation_bins.astype(np.float)
        
        
        axs[i,j].errorbar(separation_bins,corr_fit,yerr=dcorr_fit,
            fmt='.-',color='#F56B5C90',label=r'$T < 10 \, \mathrm{Myr}$')
        ########################################
        galaxy_class.get_ra_dec()
        #Old Clusters: T> 10 Myr
        galaxy_class.ra = galaxy_class.ra[np.where(ages > 1.e7)]
        galaxy_class.dec = galaxy_class.dec[np.where(ages > 1.e7)]

        #Compute TPCF
        corr,dcorr,bootstraps = bootstrap_two_point_angular(galaxy_class,
                        method='landy-szalay',Nbootstraps=100,
                        random_method=method)
        galaxy_class.corr = corr 
        galaxy_class.dcorr = dcorr

        #Isolate non-zero correlation points
        indices = np.where(galaxy_class.corr>0.0)
        corr_fit = galaxy_class.corr[indices].astype(np.float)
        dcorr_fit = galaxy_class.dcorr[indices].astype(np.float)
        separation_bins = galaxy_class.bin_centres*(1./arcsec_to_degree)
        separation_bins = separation_bins.astype(np.float)
        separation_bins = separation_bins[indices].astype(np.float)
        separation_bins = separation_bins.astype(np.float)

        axs[i,j].errorbar(separation_bins,corr_fit,yerr=dcorr_fit,
            fmt='.-',color='#4591F590',label=r'$T > 10 \, \mathrm{Myr}$')

        # General Plot stuff
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
        filename = outdir+'Combined_TPCF_AgeCompare.pdf'    
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else :
        plt.show()        

def Combined_TPCF_Ages(save=False,outdir='../Results/',indir=None,method='masked',age_group='young'):
    """
    Plot the TPCF of clusters with age < 10 Myr in all galaxies as a combined figure.

    Parameters:
        save: 
            Flag to save the plot
        outdir: 
            Output directory in which to store plot. Default is results directory.
        indir :
            Input directory from which to read the results. Default is results directory.
        method : 
            Method for which TPCF's have been computed.
        age_group : 
            Whether to do this for young or old clusters where young<10 Myr and old > 10 Myr    
    Returns:
        None


    """
    print("Plotting TPCF for {} clusters for all galaxies".format(age_group))

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

    if(age_group not in ['young','old']):
        raise myError("Please specify either young or old in age group.")


    i,j = 0,0            
    #Loop through the galaxies
    for galaxy_name in list_of_galaxies:
        if(method_dir == None):
            galaxy_dir = indir+galaxy_name+'/'            
        else :
            galaxy_dir = indir+galaxy_name+'/'+method_dir

        galaxy_class = loadObj(galaxy_dir+galaxy_name+'_summary')
        plot_class = myPlot(galaxy_class)

        separation_bins = galaxy_class.bin_centres*(1./arcsec_to_degree)
        separation_bins = separation_bins.astype(np.float)
        plot_points = np.linspace(np.min(separation_bins),np.max(separation_bins),1000)
        
        ages = galaxy_class.get_cluster_ages()
        #Clusters < 10Myr

        if(age_group == 'young'):
            galaxy_class.ra = galaxy_class.ra[np.where(ages <= 1.e7)]
            galaxy_class.dec = galaxy_class.dec[np.where(ages <= 1.e7)] 
        elif(age_group == 'old'):
            galaxy_class.ra = galaxy_class.ra[np.where(ages > 1.e7)]
            galaxy_class.dec = galaxy_class.dec[np.where(ages > 1.e7)]

        corr,dcorr,bootstraps = bootstrap_two_point_angular(galaxy_class,
                        method='landy-szalay',Nbootstraps=100,
                        random_method=method)
        

        galaxy_class.corr = corr 
        galaxy_class.dcorr = dcorr

        galaxy_class.fit_power_law(use_bounds=False)

        #Isolate non-zero correlation points
        indices = np.where(galaxy_class.corr>0.0)
        corr_fit = galaxy_class.corr[indices].astype(np.float)
        dcorr_fit = galaxy_class.dcorr[indices].astype(np.float)
        separation_bins = separation_bins[indices].astype(np.float)
        separation_bins = separation_bins.astype(np.float)
        
        axs[i,j].set_yscale('log')            
        ax2 = axs[i,j].secondary_xaxis("top",functions=(plot_class.sep_to_pc,plot_class.pc_to_sep))
        
        axs[i,j].errorbar(separation_bins,corr_fit,yerr=dcorr_fit,
            fmt='.-')

        #Fit plot
        break_theta = np.exp(galaxy_class.fit_values[3])
        break_theta_error = np.exp(galaxy_class.fit_errors[3])
        axs[i,j].plot(plot_points,np.exp(linear_function(plot_points,galaxy_class.fit_values[0],
            galaxy_class.fit_values[1],galaxy_class.fit_values[2],galaxy_class.fit_values[3])),
            ls='--',label='fit')
        axs[i,j].plot(separation_bins,corr_fit,lw=0.0,
            label=r'$\alpha_1 = {:2.1f} \pm {:2.1f}$'.format(galaxy_class.fit_values[1],galaxy_class.fit_errors[1]))
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
        if(age_group == 'young'):
            filename = outdir+'Combined_TPCF_Young.pdf'
        elif(age_group == 'old'):
            filename = outdir+'Combined_TPCF_Old.pdf'    
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else :
        plt.show()        



def Combined_TPCF_allclass(save=False,outdir='../Results/',indir=None,method='masked'):
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
    print("Comparing TPCF for different classes of clusters for all galaxies.")

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
        separation_bins = separation_bins.astype(np.float)
        
            
        axs[i,j].set_yscale('log')
        ax2 = axs[i,j].secondary_xaxis("top",functions=(plot_class.sep_to_pc,plot_class.pc_to_sep))
        
        color = ['#E13B58','#326EFA','#A5F597']
        for k in range(1,4):         
            separation_bins = galaxy_class.bin_centres*(1./arcsec_to_degree)
            separation_bins = separation_bins.astype(np.float)       
            galaxy_class.Compute_TPCF(cluster_class=k,random_method=method)  
            indices = np.where(galaxy_class.corr>0.0)
            corr_fit = galaxy_class.corr[indices].astype(np.float)
            dcorr_fit = galaxy_class.dcorr[indices].astype(np.float)
            separation_bins = separation_bins[indices].astype(np.float)                  
            axs[i,j].errorbar(separation_bins,corr_fit,yerr=dcorr_fit,
                fmt='.-',label='Class {}'.format(k),color=color[k-1])

        #Plot stuff
        #X-labels only on bottom row
        if(i==2):
            axs[i,j].set_xlabel(r"$\theta \, \left(\mathrm{arcsec} \right)$")
        #Y-labels only on left column
        if(j == 0):
            axs[i,j].set_ylabel(r"$\omega_{\mathrm{LS}}\left(\theta \right)$")
        axs[i,j].set_xscale('log')
        axs[i,j].callbacks.connect("xlim_changed", plot_class.axs_to_parsec)
        if(i+j == 0):
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
        filename = outdir+'Combined_TPCF_Classes.pdf'
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else :
        plt.show()

def Combined_RGBImages(save=False,outdir='../Results/',indir=None,method='masked'):
    
    print("Plotting HST images for all galaxies.")

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


def Combined_Clusters(save=False,outdir='../Results/',indir=None,method='masked'):

    print("Plotting spatial positions of clusters in all galaxies")

    #Create figure and axs instance
    fig = plt.figure(figsize=(20,16),constrained_layout=True)

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
        ages = galaxy_class.get_cluster_ages()
        cmap = cmr.waterlily

        ax1 = plt.subplot2grid((3,4),(i,j),fig=fig)

        x = galaxy_class.ra*3600.0
        y = galaxy_class.dec*3600.0

        #Get central pixel values
        xcen = (np.min(x)+np.max(x))/2.
        ycen = (np.min(y)+np.max(y))/2.

        #Scale offset around bounding box to ~ 5% of axes width
        offset_x = (np.max(x)-np.min(x))*0.05
        offset_y = (np.max(y)-np.min(y))*0.05

        xmin,xmax = np.min(x)-offset_x-xcen, np.max(x)+offset_x-xcen
        ymin,ymax = np.min(y)-offset_y-ycen, np.max(y)+offset_y-ycen
        ax1.set_xlim(xmin,xmax)
        ax1.set_ylim(ymin,ymax)

        im1 = ax1.scatter(x-xcen,y-ycen,s=10,c=np.log10(ages),alpha=0.5,cmap=cmap,lw=0.3)
        cbar = fig.colorbar(im1,ax = ax1,use_gridspec=False,
                            orientation='vertical',pad=0.00,aspect=30)
        

        # #Draw 100 arcsec scale bar
        import matplotlib.lines as lines
        #No of pixels in axes
        total_pixels = np.int(np.floor(ax1.transData.transform((xmax,ymax))[0]) - \
        np.floor(ax1.transData.transform((xmin,ymin))[0]))

        length_per_pixel = (xmax-xmin)/(total_pixels)
        #Convert to parsec 
        length_per_pixel = plot_class.sep_to_pc(length_per_pixel)
        #Scale bar of 50 arcsec
        length = plot_class.sep_to_pc(50)
        no_pixels = np.abs(length/length_per_pixel)
        no_pixels = no_pixels/total_pixels

        scale = lines.Line2D([0.8,0.8+no_pixels],[0.1],
                                         lw=1,color='black',
                                        transform=ax1.transAxes)

        ax1.add_line(scale)
        ax1.annotate(r'$50^{\prime \prime} = %d \, \mathrm{pc}$'%length,(0.65,0.15),xycoords='axes fraction',
                            fontsize=12)


        if(i==2):
            ax1.set_xlabel(r'$\mathrm{separation} \, (\mathrm{arcsec})$',fontsize=16)

        if(j==0):
            ax1.set_ylabel(r'$\mathrm{separation} \, (\mathrm{arcsec})$',fontsize=16)

        if(j==3):    
            cbar.ax.set_ylabel(r"$\log_{10} \, \mathrm{Age} \, (\mathrm{yr})$",
                rotation=90,labelpad=5,fontsize=20)
        
        ax1.text(0.1,0.1,r'$\mathrm{NGC}$'+' '+r'${}$'.format(galaxy_name.split('_')[1]),
            transform=ax1.transAxes)

        #Get position of subplot
        #Get position of subplot
        j +=1
        if(j==4):
            j = 0
            i +=1


    if(save):
        filename = outdir+'Combined_Clusters.pdf'
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else :
        plt.show()

def Scatter_Correlations(save=False,outdir='../Results/',indir=None,method='masked'):

    print("Plotting correlations between galaxy properties and slope of TPCF.")

    #Create figure and axs instance
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
    distance = np.zeros(np.size(list_of_galaxies))
    T_value = np.zeros(np.size(list_of_galaxies))
    SFR = np.zeros(np.size(list_of_galaxies))
    sigma_sfr = np.zeros(np.size(list_of_galaxies))
    R_25 = np.zeros(np.size(list_of_galaxies))
    alpha = np.zeros(np.size(list_of_galaxies))
    alpha_error = np.zeros(np.size(list_of_galaxies))
    index = 0

    for galaxy_name in list_of_galaxies:
        if(method_dir == None):
            galaxy_dir = indir+galaxy_name+'/'            
        else :
            galaxy_dir = indir+galaxy_name+'/'+method_dir



        galaxy_class = loadObj(galaxy_dir+galaxy_name+'_summary')
        galaxy_class.read_galaxyprops()
        distance[index] = galaxy_class.distance
        SFR[index] = galaxy_class.sfr
        sigma_sfr[index] = galaxy_class.sigma_sfr
        R_25[index] = galaxy_class.r25
        T_value[index] = galaxy_class.T_value
        T_value[index] = galaxy_class.T_value

        #Get alpha 
        AIC_single,AIC_piecewise = compare_AIC(galaxy_name)
        if(AIC_single<AIC_piecewise):
            sampler = loadObj(galaxy_dir+'/SinglePL_MCMC/'+'MCMC_sampler')
        else:
            sampler = loadObj(galaxy_dir+'/PiecePL_MCMC/'+'MCMC_sampler')

        samples = sampler.flatchain
        galaxy_class.fit_values = samples[np.argmax(sampler.flatlnprobability)]
        galaxy_class.fit_errors = samples.std(axis=0)


        alpha[index] = galaxy_class.fit_values[1]
        alpha_error[index] = galaxy_class.fit_errors[1]
        index +=1



    fig,axs = plt.subplots(nrows=2,ncols=2,figsize=(12,8))


    axs[0,0].errorbar(R_25,alpha,yerr=alpha_error,elinewidth=3.0,fmt='s',
        ms=8.0,color='#E200E690')
    axs[0,0].set_xlabel(r"$\mathrm{R}_{25} \, \left( \mathrm{kpc} \right)$")
    

    axs[0,1].errorbar(T_value,alpha,yerr=alpha_error,elinewidth=3.0,fmt='s',
        ms=8.0,color='#E200E690')
    axs[0,1].set_xlabel(r"$\mathrm{T_{\mathrm{Hubble}}}$")
    axs[0,1].set_ylabel(r"$\alpha_1$")
    axs[0,1].set_xlim(3,10)

    axs[1,0].errorbar(SFR,alpha,yerr=alpha_error,elinewidth=3.0,fmt='s',
        ms=8.0,color='#E200E690')
    axs[1,0].set_xlabel(r"$SFR_{\mathrm{UV}} \, \left( M_\odot \, \mathrm{yr}^{-1} \right)$")
    axs[1,0].set_xlim(-0.5,12)

    axs[1,1].errorbar(sigma_sfr,alpha,yerr=alpha_error,elinewidth=3.0,fmt='s',
        ms=8.0,color='#E200E690')
    axs[1,1].set_xlabel(r"$\Sigma_{\mathrm{SFR}} \, \left(\mathrm{M}_{\odot} \, \mathrm{pc}^{-2} \, \mathrm{yr}^{-1} \right)$")
    axs[1,1].set_ylabel(r"$\alpha_1$")


    if(save):
        filename = outdir+'Correlations_alpha1.pdf'
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else :
        plt.show()



def Compare_TPCFMethods(save=False,outdir='../Results/',indir=None,method='masked'):

    print("Plotting comparison of methods to compute TPCF for all galaxies.")

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

        print("Computing for {}".format(galaxy_name))

        galaxy_class = loadObj(galaxy_dir+galaxy_name+'_summary')
        plot_class = myPlot(galaxy_class)

        #Secondary axis
        separation_bins = galaxy_class.bin_centres*(1./arcsec_to_degree)
        separation_bins = separation_bins.astype(np.float)

        
        indices = np.where(galaxy_class.corr>0.0)
        corr_fit = galaxy_class.corr[indices].astype(np.float)
        dcorr_fit = galaxy_class.dcorr[indices].astype(np.float)
        separation_bins = separation_bins[indices].astype(np.float)
            
        axs[i,j].errorbar(separation_bins,corr_fit,yerr=dcorr_fit,
            fmt='.-',label='Landy-Szalay')
        ax2 = axs[i,j].secondary_xaxis("top",functions=(plot_class.sep_to_pc,plot_class.pc_to_sep))

        #Compute Standard TPCF
        galaxy_class.Compute_TPCF(random_method='masked',tpcf_method='standard')
        indices = np.where(galaxy_class.corr>0.0)
        corr_fit = galaxy_class.corr[indices].astype(np.float)
        dcorr_fit = galaxy_class.dcorr[indices].astype(np.float)
        separation_bins = galaxy_class.bin_centres*(1./arcsec_to_degree)
        separation_bins = separation_bins.astype(np.float)
        separation_bins = separation_bins[indices].astype(np.float)

        axs[i,j].errorbar(separation_bins,corr_fit,yerr=dcorr_fit,
            fmt='.-',label='Standard')

        #Plot stuff
        #X-labels only on bottom row
        if(i==2):
            axs[i,j].set_xlabel(r"$\theta \, \left(\mathrm{arcsec} \right)$")
        #Y-labels only on left column
        if(j == 0):
            axs[i,j].set_ylabel(r"$\omega_{\mathrm{LS}}\left(\theta \right)$")
        axs[i,j].set_xscale('log')
        axs[i,j].set_yscale('log')
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
        print("Saving...")
        filename = outdir+'Compare_TPCFMethods.pdf'
        plt.savefig(filename,bbox_inches='tight')
        plt.close()
    else :
        plt.show()




if __name__ == "__main__":
    print("Preparing plots of paper.")
    #Combined_TPCF(save=True,method='masked')
    Combined_TPCF_CompareAges(save=True,method='masked')
    #Combined_TPCF_allclass(save=True,method='masked')
    #Combined_TPCF_Ages(save=True,method='masked',age_group='young')
    #Combined_TPCF_Ages(save=True,method='masked',age_group='old')
    #Combined_Clusters(save=True,method='masked')
    #Scatter_Correlations(save=True,method='masked')

