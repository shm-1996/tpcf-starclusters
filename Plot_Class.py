from header import *
from TPCF import *

class myPlot():
    """
    Class that provides functionality for different plots for a galaxy
    Parameters
        galaxy: class Galaxy
            Instance of class galaxy that contains its properties
    """
    def __init__(self,galaxy):
        self.galaxy = galaxy
        
        
        
    def plot_TPCF(self,save=False,filename=None):
        """
        Plot TPCF of a galaxy

        Parameters
        ----------
        None
        
        Returns
        -------
        None
        
        """

        if(self.corr is )
        if((self.galaxy.fit_values == [0,0,0,0,0]) or (self.galaxy.fit_errors == [0,0,0,0,0])):
            # Power law not fitted yet. Fit 
            print("Fitting power law. Not fitted yet.")
            self.galaxy.fit_power_law()
        
        fig,axs = plt.subplots(ncols=1)
        ax2 = axs.secondary_xaxis("top",functions=(self.sep_to_pc,self.pc_to_sep))
        separation_bins = (self.galaxy.bins[1:]+self.galaxy.bins[:-1])/2
        separation_bins*=(1./arcsec_to_degree)
        
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
            axs.plot(separation_bins,np.exp(linear_function(separation_bins,self.galaxy.fit_values[0],
                self.galaxy.fit_values[1],self.galaxy.fit_values[2],self.galaxy.fit_values[3],self.galaxy.fit_values[4])),
                ls=':',label='fit')
        except AttributeError:
            print("Power-law not fitted to TPCF yet. Fitting now.")
            self.galaxy.fit_power_law()
            axs.plot(separation_bins,np.exp(linear_function(separation_bins,self.galaxy.fit_values[0],
                self.galaxy.fit_values[1],self.galaxy.fit_values[2],self.galaxy.fit_values[3],self.galaxy.fit_values[4])),
                ls=':',label='fit')
        
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
        else :
            plt.show()

    def plot_TPCF_allclass(self,random_method = 'masked_radial',
        save=False,filename=None,verbose=False):
        """
        Plot TPCF of a galaxy for all classes comparing between them

        Parameters
        ----------
        None
        
        Returns
        -------
        None
        
        """
        if(verbose):
            print("Plotting comparison of TPCF for different classes.")
        #Initialise figure
        fig,axs = plt.subplots(ncols=1)
        ax2 = axs.secondary_xaxis("top",functions=(self.sep_to_pc,self.pc_to_sep))
        separation_bins = (self.galaxy.bins[1:]+self.galaxy.bins[:-1])/2
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
        axs.callbacks.connect("xlim_changed", self.axs_to_parsec)
        axs.legend()
        ax2.set_xlabel(r'$\delta x \, \left( \mathrm{pc} \right) $')
        if(save):
            if(filename == None):
                filename = self.galaxy.outdir+'/{}_Classes_TPCF.pdf'.format(self.galaxy.name)
            plt.savefig(filename,bbox_inches='tight')
        else :
            plt.show()


    def plot_clusters(self,save=False,filename=None):
        """
        Plot spatial distribution of clusters in the galaxy

        Parameters
        ----------
        None
        
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
                filename = self.galaxy.outdir+'/{}_clusters'.format(self.galaxy.galaxy)
            plt.savefig(filename,bbox_inches='tight')
        else :
            plt.show()


    def plot_random(self,save=False,random_method='masked_radial',filename=None):
        """
        Plot spatial distribution of random catalog
        Parameters
        ----------
        None
        
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
                filename = self.galaxy.outdir+'/{}_random'.format(self.galaxy.galaxy) + \
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
        # axs.bar([0,1,2,3,4],[N0,N1,N2,N3,N4],color='#F59005',tick_label=label,
        #     label='{}'.format(self.galaxy.galaxy))
        print('{}'.format(self.galaxy.galaxy))
        axs.bar([0,1,2,3,4],[N0,N1,N2,N3,N4],color='#F59005',tick_label=label)
        axs.set_xlabel('Cluster Class')
        axs.set_ylabel('Number')
        axs.legend()
        if(save):
            if(filename == None) :
                filename = self.galaxy.outdir+'/{}_ClassDist'.format(self.galaxy.galaxy)
            plt.savefig(filename,bbox_inches='tight')
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


