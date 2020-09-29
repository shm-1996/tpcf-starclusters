from header import *
from TPCF import *
from sklearn.neighbors import KDTree

class myPlot():
    """
    Class that provides functionality for different plots for a galaxy
    Parameters
        galaxy: class Galaxy
            Instance of class galaxy that contains its properties
    """
    def __init__(self,galaxy):
        self.galaxy = galaxy
        
        
        
    def plot_TPCF(self,save=False,function='piecewise',filename=None):
        """
        Plot TPCF of a galaxy

        Parameters
        ----------
        None
        
        Returns
        -------
        None
        
        """
        fig,axs = plt.subplots(ncols=1)
        ax2 = axs.secondary_xaxis("top",functions=(self.sep_to_pc,self.pc_to_sep))
        separation_bins = self.galaxy.bin_centres*(1./arcsec_to_degree)

        if(function not in ['piecewise','smooth']):
            raise ValueError("This funtional form does not exist.")
        
        #Try plotting directly if TPCF computed. Else compute.
        try:
            axs.errorbar(separation_bins,self.galaxy.corr,yerr=self.galaxy.dcorr,
            fmt='.-')
        except AttributeError:
            print("TPCF not computed computing first.")
            self.galaxy.Compute_TPCF()
            axs.errorbar(separation_bins,self.galaxy.corr,yerr=self.galaxy.dcorr,
            fmt='.-')

        try:
            
            plot_points = np.linspace(np.min(separation_bins),np.max(separation_bins),1000)
            if(function == 'piecewise'):
                axs.plot(plot_points,np.exp(linear_function(plot_points,self.galaxy.fit_values[0],
                    self.galaxy.fit_values[1],self.galaxy.fit_values[2],self.galaxy.fit_values[3])),
                    ls='--',label='fit')
                break_theta = np.exp(self.galaxy.fit_values[3])
                break_theta_error = np.exp(self.galaxy.fit_errors[3])
            elif(function == 'smooth'):
                axs.plot(plot_points,smooth_function(plot_points,self.galaxy.fit_values[0],
                    self.galaxy.fit_values[1],self.galaxy.fit_values[2],self.galaxy.fit_values[3]),
                    ls='--',label='fit')
                break_theta = self.galaxy.fit_values[3]
                break_theta_error = self.galaxy.fit_errors[3]

                
            axs.plot(separation_bins,self.galaxy.corr,lw=0.0,
                label=r'$\alpha_1 = {:2.1f} \pm {:2.1f}$'.format(self.galaxy.fit_values[1],self.galaxy.fit_errors[1]))
            axs.plot(separation_bins,self.galaxy.corr,lw=0.0,
                label=r'$\alpha_2 = {:2.1f} \pm {:2.1f}$'.format(self.galaxy.fit_values[2],self.galaxy.fit_errors[2]))
            axs.axvline(break_theta,ls=':',label=r'$\beta = {:2.1f} \pm {:2.1f}$'.format(break_theta,
                break_theta_error))
            


        except AttributeError:
            print("Power-law not fitted to TPCF yet. Fitting now.")
            self.galaxy.fit_power_law(function=function)
            plot_points = np.linspace(np.min(separation_bins),np.max(separation_bins),1000)
            break_theta = np.exp(self.galaxy.fit_values[3])
            break_theta_error = np.exp(self.galaxy.fit_errors[3])
            axs.plot(plot_points,np.exp(linear_function(plot_points,self.galaxy.fit_values[0],
                self.galaxy.fit_values[1],self.galaxy.fit_values[2],self.galaxy.fit_values[3])),
                ls='--',label='fit')
            axs.plot(separation_bins,self.galaxy.corr,lw=0.0,
                label=r'$\alpha_1 = {:2.1f} \pm {:2.1f}$'.format(self.galaxy.fit_values[1],self.galaxy.fit_errors[1]))
            axs.plot(separation_bins,self.galaxy.corr,lw=0.0,
                label=r'$\alpha_2 = {:2.1f} \pm {:2.1f}$'.format(self.galaxy.fit_values[2],self.galaxy.fit_errors[2]))

            axs.axvline(break_theta,ls=':',label=r'$\beta = {:2.1f} \pm {:2.1f}$'.format(break_theta,
                break_theta_error))
            
        
        axs.set_xlabel(r"$\theta \, \left(\mathrm{arcsec} \right)$")
        axs.set_ylabel(r"$\omega_{\mathrm{LS}}\left(\theta \right)$")
        axs.set_xscale('log')
        axs.set_yscale('log')
        axs.callbacks.connect("xlim_changed", self.axs_to_parsec)
        axs.legend()
        ax2.set_xlabel(r'$\delta x \, \left( \mathrm{pc} \right) $')
        if(save):
            if(filename == None):
                filename = self.galaxy.outdir+'/{}_TPCF.pdf'.format(self.galaxy.name)
            plt.savefig(filename,bbox_inches='tight')
            plt.close()
        else :
            plt.show()


    def plot_TPCF_allclass(self,random_method = 'masked_radial',
        save=False,filename=None,verbose=False):
        """
        Plot TPCF of a galaxy for all classes comparing between them

       Parameters
        ----------
        save : boolean
            Flag to save the plot, else just show.
        random_method: string
            random method to use
        filename : string
            File to save to, else default filename 
        verbose : string
            Whether to print verbose
        
        Returns
        -------
        None
        
        """
        if(verbose):
            print("Plotting comparison of TPCF for different classes.")
        #Initialise figure
        fig,axs = plt.subplots(ncols=1)
        ax2 = axs.secondary_xaxis("top",functions=(self.sep_to_pc,self.pc_to_sep))
        separation_bins = self.galaxy.bin_centres
        separation_bins*=(1./arcsec_to_degree)

        #Compute TPCF for each class
        for i in range(1,4):
            if(verbose):
                print("Computing TPCF for class {} clusters".format(i))
            self.galaxy.Compute_TPCF(cluster_class=i,random_method=random_method,
                verbose=verbose)
            axs.errorbar(separation_bins,self.galaxy.corr,yerr=self.galaxy.dcorr,
                fmt='.-',label='Class {}'.format(i))
        # Combined
        if(verbose):
            print("Computing TPCF for combined class clusters")

        #TODO: Figure out how to do this. Currently the properties for 
        self.galaxy.Compute_TPCF(cluster_class=-1,random_method=random_method,
            verbose=verbose)
        axs.errorbar(separation_bins,self.galaxy.corr,yerr=self.galaxy.dcorr,
            fmt='.-',label='Class 1+2+3')


        #Rest of plot
        axs.set_xlabel(r"$\theta \, \left(\mathrm{arcsec} \right)$")
        axs.set_ylabel(r"$\omega_{\mathrm{LS}}\left(\theta \right)$")
        axs.set_xscale('log')
        axs.set_yscale('log')

        #Secondary axis
        axs.callbacks.connect("xlim_changed", self.axs_to_parsec)
        
        axs.legend()
        ax2.set_xlabel(r'$\delta x \, \left( \mathrm{pc} \right) $')
        if(save):
            if(filename == None):
                filename = self.galaxy.outdir+'/{}_Classes_TPCF.pdf'.format(self.galaxy.name)
            plt.savefig(filename,bbox_inches='tight')
            plt.close()
        else :
            plt.show()


    def plot_clusters(self,save=False,filename=None):
        """
        Plot spatial distribution of clusters in the galaxy

       Parameters
        ----------
        save : boolean
            Flag to save the plot, else just show.
        filename : string
            File to save to, else default filename 
        
        Returns
        -------
        None
        
        """ 
        hdu = fits.open(self.galaxy.fits_file)[0]
        wcs = WCS(hdu.header)

        fig = plt.figure()
        ax1 = fig.add_subplot(111,projection=wcs)
        im = ax1.scatter(self.galaxy.ra,self.galaxy.dec)
        ax1.set_xlabel(r"$\mathrm{Right \; Ascension \; (J2000)}$",fontsize=20)
        ax1.set_ylabel(r"$\mathrm{Declination \; (J2000)}$",labelpad=-0.8,
                       fontsize=20)
        if(save):
            if(filename == None):
                filename = self.galaxy.outdir+'/{}_clusters'.format(self.galaxy.name)
            plt.savefig(filename,bbox_inches='tight')
            plt.close()
        else :
            plt.show()


    def plot_random(self,save=False,random_method='masked_radial',filename=None):
        """
        Plot spatial distribution of random catalog
        Parameters
        ----------
        save : boolean
            Flag to save the plot, else just show.
        random_method : string
            Random method in use
        filename : string
            File to save to, else default filename 
        
        Returns
        -------
        None
        
        """ 
        #Obtain one instance of the random sample
        if(random_method == 'uniform') :
            #Draw uniform random sample with N points
            ra_R, dec_R = uniform_sphere(self.galaxy,2 * len(self.galaxy.ra))
                                      
        elif(random_method == 'masked') :
            #Draw a masked random sample with N points
            ra_R, dec_R = masked_random_sample(self.galaxy,10*len(self.galaxy.ra))

        elif(random_method == 'masked_radial') :
            #Draw a random sample 
            ra_R, dec_R = masked_radial_random_sample(self.galaxy,1000*len(self.galaxy.ra))

        hdu = fits.open(self.galaxy.fits_file)[0]
        wcs = WCS(hdu.header)
        fig = plt.figure(constrained_layout=True)
        ax1 = fig.add_subplot(111,projection=wcs)
        im = ax1.scatter(ra_R,dec_R)
        ax1.set_xlabel(r"$\mathrm{Right \; Ascension \; (J2000)}$",fontsize=20)
        ax1.set_ylabel(r"$\mathrm{Declination \; (J2000)}$",labelpad=-0.8,
                       fontsize=20)
        if(save):
            if(filename == None) :
                filename = self.galaxy.outdir+'/{}_random'.format(self.galaxy.name) + \
                '_{}'.format(random_method)
            plt.savefig(filename,bbox_inches='tight')
        else :
            plt.show()

    def class_distribution(self,save=False,filename=None):
        #Read file for distribution of classes
        file = np.loadtxt(self.galaxy.catalog_file)
        N0 = np.size(np.where(file[:,33]==0))
        N1 = np.size(np.where(file[:,33]==1))
        N2 = np.size(np.where(file[:,33]==2))
        N3 = np.size(np.where(file[:,33]==3))
        N4 = np.size(np.where(file[:,33]==4))

        #Plot now
        fig,axs = plt.subplots(ncols=1)
        label = ['Class 0', 'Class 1', 'Class 2', 'Class 3','Class 4']
        axs.bar([0,1,2,3,4],[N0,N1,N2,N3,N4],color='#F59005',tick_label=label)
        axs.set_xlabel('Cluster Class')
        axs.set_ylabel('Number')
        if(save):
            if(filename == None) :
                filename = self.galaxy.outdir+'/{}_ClassDist'.format(self.galaxy.name)
            plt.savefig(filename,bbox_inches='tight')
            plt.close()
        else :
            plt.show()


    def mass_histogram(self,save=False,filename=None):
        """
        Plot distribution of masses in each class of clusters.
        Parameters
        ----------
        save : boolean
            Flag to save the plot, else just show.
        filename : string
            File to save to, else default filename 
        
        Returns
        -------
        None
        
        """ 

        #Read file for distribution of masses in each class
        file = np.loadtxt(self.galaxy.catalog_file)
        M1 = file[np.where(file[:,33]==1)][:,19]
        M2 = file[np.where(file[:,33]==2)][:,19]
        M3 = file[np.where(file[:,33]==3)][:,19]
        
        #Plot now
        fig,axs = plt.subplots(ncols=1)
        axs.hist(np.log10(M1),bins='auto',histtype='step',log=True,
            label='Class 1',color='#F51557')
        axs.hist(np.log10(M2),bins='auto',histtype='step',log=True,
            label='Class 2',color='#33A7F4')
        axs.hist(np.log10(M3),bins='auto',histtype='step',log=True,
            label='Class 3',color='#28F56E')

        axs.set_xlabel(r'$\log_{10} \, \mathrm{Mass} \, (M_{\odot})$')
        axs.set_ylabel(r'$\mathrm{Number}$')
        axs.legend()
        if(save):
            if(filename == None) :
                filename = self.galaxy.outdir+'/{}_MassDist'.format(self.galaxy.name)
            plt.savefig(filename,bbox_inches='tight')
            plt.close()
        else :
            plt.show()


    def age_histogram(self,save=False,filename=None):
        """
        Plot distribution of agees in each class of clusters.
        Parameters
        ----------
        save : boolean
            Flag to save the plot, else just show.
        filename : string
            File to save to, else default filename 
        
        Returns
        -------
        None
        
        """ 

        #Read file for distribution of ages in each class
        file = np.loadtxt(self.galaxy.catalog_file)
        A1 = file[np.where(file[:,33]==1)][:,16]
        A2 = file[np.where(file[:,33]==2)][:,16]
        A3 = file[np.where(file[:,33]==3)][:,16]

        # # Ages are in yr. Convert to Myr.
        # A1 /= 1.e6
        # A2 /= 1.e6
        # A3 /= 1.e6

        
        #Plot now
        fig,axs = plt.subplots(ncols=1)
        axs.hist(np.log10(A1),bins='auto',histtype='step',log=True,
            label='Class 1',color='#F51557')
        axs.hist(np.log10(A2),bins='auto',histtype='step',log=True,
            label='Class 2',color='#33A7F4')
        axs.hist(np.log10(A3),bins='auto',histtype='step',log=True,
            label='Class 3',color='#28F56E')

        axs.set_xlabel(r'$\log_{10} \, \mathrm{Age} \, (\mathrm{Myr})$')
        axs.set_ylabel(r'$\mathrm{Number}$')
        axs.legend()
        if(save):
            if(filename == None) :
                filename = self.galaxy.outdir+'/{}_AgeDist'.format(self.galaxy.name)
            plt.savefig(filename,bbox_inches='tight')
            plt.close()
        else :
            plt.show()

    def bin_distribution(self,save=False,filename=None):
        """
        Plot distribution of pairs in each TPCF bin. 
        Parameters
        ----------
        save : boolean
            Flag to save the plot, else just show.
        filename : string
            File to save to, else default filename 
        Returns
        -------
        None
        
        """         

        #Get no of pairs for each class using KDTree two point correlation
        counts_DD = np.zeros((3,self.galaxy.no_bins))
        bins_transform = angular_dist_to_euclidean_dist(self.galaxy.bins)
        for i in range(1,4):
            self.galaxy.get_ra_dec(cluster_class=i)
            xyz_clusters = np.asarray(ra_dec_to_xyz(self.galaxy.ra, 
                                            self.galaxy.dec), order='F').T
            KDT_D = KDTree(xyz_clusters)
            counts = KDT_D.two_point_correlation(xyz_clusters, bins_transform)
            counts_DD[i-1] = np.diff(counts)
        #Reset RA/DEC
        self.galaxy.get_ra_dec(cluster_class=-1)    

        fig,axs = plt.subplots(ncols=1)
        ax2 = axs.secondary_xaxis("top",functions=(self.sep_to_pc,self.pc_to_sep))
        
        axs.bar(self.galaxy.bins_arcsec[:-1],counts_DD[0],
            align='edge',width=np.diff(self.galaxy.bins_arcsec),
            color='#F51557',label='Class 1')
        axs.bar(self.galaxy.bins_arcsec[:-1],counts_DD[1],
                    align='edge',width=np.diff(self.galaxy.bins_arcsec),
                    bottom=counts_DD[0],color='#4983FC',label='Class 2')
        axs.bar(self.galaxy.bins_arcsec[:-1],counts_DD[2],
                    align='edge',width=np.diff(self.galaxy.bins_arcsec),
                    bottom=counts_DD[1],color='#FAC90E',label='Class 3')
            
        axs.set_xlabel(r"$\theta \, \left(\mathrm{arcsec} \right)$")
        axs.set_ylabel(r'$\mathrm{N_{\mathrm{pairs}}}$')
        axs.set_xscale('log')
        axs.set_yscale('log')
        axs.legend(loc='upper left')
        
        axs.callbacks.connect("xlim_changed", self.axs_to_parsec)
        ax2.set_xlabel(r'$\delta x \, \left( \mathrm{pc} \right) $')
        if(save):
            if(filename == None) :
                filename = self.galaxy.outdir+'/{}_BinDist'.format(self.galaxy.name)
            plt.savefig(filename,bbox_inches='tight')
            plt.close()
        else :
            plt.show()


    def galaxy_image(self,save=False,filename=None):

        """
        Plot optical HST image of galaxy. 
        Parameters
        ----------
        save : boolean
            Flag to save the plot, else just show.
        filename : string
            File to save to, else default filename 
        
        Returns
        -------
        None
        
        """         

        hdu = fits.open(self.galaxy.fits_file)[0]
        wcs = WCS(hdu.header)
        fig = plt.figure(constrained_layout=True)
        ax1 = fig.add_subplot(111,projection=wcs)

        with np.errstate(divide='ignore', invalid='ignore'):
            im = ax1.imshow(np.log10(hdu.data),vmin=-2.0)
        cbar = fig.colorbar(im,ax = ax1,orientation='vertical')
        cbar.ax.set_ylabel(r"$\log_{10} \, I$",rotation=90,labelpad=5,fontsize=20)
        ax1.set_xlabel(r"$\mathrm{Right \; Ascension \; (J2000)}$",fontsize=20)
        ax1.set_ylabel(r"$\mathrm{Declination \; (J2000)}$",labelpad=-0.8,
                       fontsize=20)

        if(save):
            if(filename == None) :
                filename = self.galaxy.outdir+'/{}_HSTImage'.format(self.galaxy.name)
            plt.savefig(filename,bbox_inches='tight')
            plt.close()
        else :
            plt.show()




    def sep_to_pc(self,sep):
        """
        Angular separation in arcsec to parsec

        Parameters
        ----------
        sep : float
            angular separation in arcsec
        Returns
        -------
        pc : float
             The separation in parsec for given angular separation
        """

        distance = self.galaxy.distance*const.Parsec*1.e6
        arcsec_to_pc = u.arcsec.to(u.radian)*distance/(const.Parsec)
        return sep*arcsec_to_pc
    def pc_to_sep(self,pc) :
        """
        Linear separation in parsec to arcsec

        Parameters
        ----------
        pc : float
             The separation in parsec 
        Returns
        -------
        sep : float
            angular separation in arcsec
        """
        distance = self.galaxy.distance*const.Parsec*1.e6
        pc_to_radian = pc*const.Parsec/distance
        return pc_to_radian*u.radian.to(u.arcsec)
    def axs_to_parsec(axs) :
        """
        Draw twin axes corresponding to angular limits
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


