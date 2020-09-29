from header import *
from astropy.coordinates import SkyCoord
from TPCF import *
import copy

# Some defaults 
default_bin_limits = 10,5000
default_no_of_bins = 10

class Galaxy(object):
    """
    A class describing the properties of the galaxy, and
    providing methods to perform calculations using those properties.

    Parameters
        galaxy_name : string
          Name of galaxy 
        filename : string
          name of file from which to read galaxy description
        verbose : Boolean
              print out information about the galaxy as we read it
    Class attributes
       name : string
          Name of the galaxy
       distance : float
          Distance to the galaxy in Mpc 
       inclination : float
          Inclination of the galaxy in the line-of-sight
       inclination_flag : bool
          Flag to perform correction for inclination in the galaxies
       inclination_done : bool
          Flag to inform that inclination correction already done
       centre : tuple
          ra,dec of galaxy centre in fk5 notation
        no_bins : integer
          Number of bins to compute TPCF
       bin_limits : tuple
          limits of the TPCF bins in arcsecs
       catalog_file : string
          filename of the cluster catalog file
       fit_values : tuple
          Obtained parameters (5) for a power-law fit to the TPCF
       fit_errors : tuple
          Obtained errors to the parameters for a power-law fit to the TPCF
          composition of the galaxy
    """


    ####################################################################
    # Method to initialize
    ####################################################################
    def __init__(self, galaxy_name=None, filename=None,verbose=False):
        """
        Parameters
            galaxy_name : string
              Name of galaxy
            filename : string
              name of file from which to read galaxy description
            verbose : Boolean
              print out information about the galaxy as we read it

        Returns
            None

        """
        self.name= None
        self.distance = 0.0
        self.inclination = 0.0
        self.deproject_galaxy = True
        self.inclination_done = False
        self.centre = 0,0
        self.no_bins = 0
        self.bin_limits = [0,0]
        self.catalog_file = None
        self.outdir = None
        self.fits_file = None
        

        # Read from info file if provided
        if(filename != None):
            self.readFile(filename,verbose=verbose)

        #If not provided see if name of galaxy provided
        if(filename is None):
            if(galaxy_name):
                self.name = galaxy_name.strip()
                self.defaultGalaxyInfo(verbose=verbose)
            else:
                raise myError("No provided galaxy or information file.")
        # Obtain RA/DEC of all clusters
        self.get_ra_dec()
        # Set bins based on bin limits & no of bins
        self.set_bins()


    ####################################################################
    # Method to set default galaxy properties
    ####################################################################

    def defaultGalaxyInfo(self,verbose=False):
        """
        Sets default galaxy information by reading its info file, or creating it if
        it does not exist.
        Parameters
            verbose : Boolean
              print out information about the galaxy as we read it

        Returns
            None
        """
        
        if self.name in list_of_galaxies:
            if(verbose):
                print("Starting analysis for "+self.name)

            filename = os.path.abspath('../Data/Galaxy_Information/'+
                self.name+'.info')
            #Check catalog file exists else throw error
            if(not os.path.exists(filename)) :
                #Create file if it does not exist
                if(verbose):
                    print("Creating galaxy info file as it does not exist.")
                self.create_file(verbose=verbose)
            else:
                self.readFile(filename=filename,verbose=verbose)

            self.outdir = os.path.abspath('../Results/Galaxies/{}/'.format(self.name))
            if(verbose):
                print("Setting output directory to {}".format(self.outdir))
                
            if(os.path.exists(self.outdir+'tpcf_bins.pkl')):
                if(verbose):
                   print("Loading TPCF values from file.") 
                self.bins = loadObj(self.outdir+'tpcf_bins')
                self.corr = loadObj(self.outdir+'tpcf_corr')
                self.dcorr = loadObj(self.outdir+'tpcf_dcorr')
                self.bootstraps = loadObj(self.outdir+'tpcf_bootstraps')
        
        else :
            raise myError("The galaxy information+catalog is not available."+
                " Create information file or manually input info.")



    ####################################################################
    # Method to create galaxy info file
    ####################################################################

    def create_file(self,verbose=False):
        """
        Creates info file for galaxy, by reading information from Table I
        of Calzetti et al 2015 (LEGUS). 
        
        Parameters
            verbose : Boolean
              print out information about the galaxy as we read it

        Returns
            None

        """
        info_directory = os.path.abspath('../Data/Galaxy_Information')
        
        #Read galaxy info from Table 1 Calzetti et al 2015
        if(verbose):
            print("Reading galaxy information from Table 1 of Calzetti et al 15.")
        Legus_Table = info_directory+'/Calzetti_Table.txt'

        # Read columns from table
        galaxy_names = np.loadtxt(Legus_Table,usecols=0,delimiter='\t',dtype=str)
        galaxy_inclinations = np.loadtxt(Legus_Table,usecols=4,delimiter='\t',dtype=float)
        galaxy_distances = np.loadtxt(Legus_Table,usecols=5,delimiter='\t',dtype=float)

        #Convert self galaxy name to table name format
        galaxy_formatted = self.name.split("_")[1]
        galaxy_formatted = '%04d'%int(galaxy_formatted)
        galaxy_formatted = self.name.split("_")[0] + ' ' + galaxy_formatted

        index = np.where(galaxy_formatted == galaxy_names)
        distance = galaxy_distances[index][0]
        inclination = galaxy_inclinations[index][0]
        

        #Setting default bin info, defined in header.py
        bin_limits = default_bin_limits
        if(verbose):
            print("Setting default bin limits of {}".format(bin_limits))
        no_bins = default_no_of_bins
        if(verbose):
            print("Setting default no of bins = {}".format(no_bins))
        catalog_file = os.path.abspath('../Data/Cluster_Catalogues/' + 
            self.name+'.tab')
        if(verbose):
            print("Setting Catalog file = {}".format(catalog_file))
        fits_file = os.path.abspath('../Data/HST_Images/' + 
            self.name+'.fits')
        if(verbose):
            print("Setting fits file = {}".format(fits_file))
        region_file = os.path.abspath('../Data/Region_Files/') + '/' +\
                self.name + '.reg'
        if(verbose):
            print("Setting region file = {}".format(region_file))  
        #Write to info file
        info_file = info_directory + '/' + self.name + '.info'
        #Write in required format
        try: 
            fp = open(info_file,'w')
        except IOError:
            raise myError("Cannot write to file "+info_file)


        # Read Position angles from position angle file
        if(verbose):
            print("Reading position angle information from Position_Angles.txt.")
        Position_AngleFile = info_directory+'/Position_Angles.txt'
        galaxy_names = np.loadtxt(Position_AngleFile,usecols=0,delimiter='\t',dtype=str)
        galaxy_pa = np.loadtxt(Position_AngleFile,usecols=1,delimiter='\t',dtype=float)
        #Convert self galaxy name to table name format
        galaxy_formatted = self.name.split("_")[1]
        galaxy_formatted = '%04d'%int(galaxy_formatted)
        galaxy_formatted = self.name.split("_")[0] + ' ' + galaxy_formatted

        index = np.where(galaxy_formatted == galaxy_names)
        pa = galaxy_pa[index][0]

        self.distance = distance 
        self.inclination = inclination
        self.no_bins = no_bins
        self.bin_limits = bin_limits 
        self.catalog_file = catalog_file
        self.fits_file = fits_file
        self.pa = pa
        self.region_file = region_file


        if(verbose):
            print("Writing info file for {}".format(self.name))

        
        #Write header
        fp.write("## Info file for {} \n".format(self.name))
        fp.write('################################################# \n\n')

        #Name of galaxy
        fp.write("# Name of Galaxy \n")
        fp.write("Name = {} \n\n".format(self.name))

        # Distance
        fp.write("# Distance in Mpc \n")
        fp.write("Distance = {} \n\n".format(self.distance))

        # Inclination
        fp.write("# Inclination in degrees \n")
        fp.write("Inclination = {} \n\n".format(self.inclination))

        # Position Angle
        fp.write("# Position Angle in degrees \n")
        fp.write("Position_Angle = {} \n\n".format(self.pa))

        # Number of bins
        fp.write("# Number of bins \n")
        fp.write("No_Bins = {} \n\n".format(self.no_bins))

        #Bin Limits
        fp.write("# Bin Limits in parsec\n")
        fp.write("Bin_Low = {} \n".format(self.bin_limits[0]))
        fp.write("Bin_High = {} \n\n".format(self.bin_limits[1]))

        # Catalog file
        fp.write("# Catalog file \n")
        fp.write("Catalog_File = {} \n\n".format(os.path.abspath(self.catalog_file)))

        # Fits file
        fp.write("# Fits file \n")
        fp.write("Fits_File = {} \n\n".format(os.path.abspath(self.fits_file)))

        # Region file
        fp.write("# Region file \n")
        fp.write("Region_File = {} \n\n".format(os.path.abspath(self.region_file)))

        #Close file
        fp.close()


    ####################################################################
    # Method to read galaxy metadata from a file
    ####################################################################
    def readFile(self,filename,verbose=False): 
        #Try reading file
        try :
            fp = open(filename,'r')
        except IOError :
            raise myError("cannot open file "+filename)
        for line in fp : 
            # Skip empty and comment lines
            if line=='\n':
                continue
            if line.strip()[0] == "#":
                continue
             # Break line up based on equal sign
            linesplit = line.split("=")
            if len(linesplit) < 2:
                raise myError("Error parsing input line: "+line)
            if linesplit[1] == '':
                raise myError("Error parsing input line: "+line)

            # Trim trailing comments from portion after equal sign
            linesplit2 = linesplit[1].split('#')

            # Read stuff based on the token that precedes the equal sign
            if linesplit[0].upper().strip() == 'NAME':
                self.name = str(linesplit2[0]).strip()
                if(verbose) :
                    print("Setting galaxy = "+str(self.name))
            elif linesplit[0].upper().strip() == 'DISTANCE':
                self.distance = float(linesplit2[0])
                if(verbose):
                    print("Setting distance = "+str(self.distance) + " Mpc")
            elif linesplit[0].upper().strip() == 'INCLINATION':
                self.inclination = float(linesplit2[0])
                if(verbose):
                    print("Setting inclination = "+str(self.inclination) + " degrees")
            elif linesplit[0].upper().strip() == 'POSITION_ANGLE':
                self.pa = float(linesplit2[0])
                if(verbose):
                    print("Setting position angle = "+str(self.pa) + " degrees")
            elif linesplit[0].upper().strip() == 'NO_BINS':
                self.no_bins = int(linesplit2[0])
                if(verbose):
                    print("Setting number of bins = "+str(self.no_bins))
            elif linesplit[0].upper().strip() == 'BIN_LOW':
                self.bin_limits[0] = float(linesplit2[0])
                if(verbose):
                    print("Setting lower bin limit for TPCF = "+str(self.bin_limits[0]) +
                         " parsec")
            elif linesplit[0].upper().strip() == 'BIN_HIGH':
                self.bin_limits[1] = float(linesplit2[0])
                if(verbose):
                    print("Setting upper bin limit for TPCF = "+str(self.bin_limits[1]) +
                        " parsec")
            elif linesplit[0].upper().strip() == 'CATALOG_FILE':
                self.catalog_file = os.path.abspath(str(linesplit2[0]).strip())
                if(verbose):
                    print("Setting cluster catalog file = "+str(self.catalog_file))

            elif linesplit[0].upper().strip() == 'FITS_FILE':
                self.fits_file = os.path.abspath(str(linesplit2[0]).strip())

                if(verbose):
                    print("Setting HST fits file = "+str(self.fits_file))

            elif linesplit[0].upper().strip() == 'REGION_FILE':
                self.region_file = os.path.abspath(str(linesplit2[0]).strip())

                if(verbose):
                    print("Setting HST region file = "+str(self.region_file))

            else:
                # Line does not correspond to any known keyword, so
                # throw an error
                raise myError("unrecognized token " +
                    linesplit[0].strip() + " in file " + filename)

        # Close file
        fp.close()

        #Compute ra/dec of centre of galaxy
        ra_dec = SkyCoord.from_name(self.name)
        ra = ra_dec.ra.value
        dec = ra_dec.dec.value
        self.centre = ra,dec

        #Check catalog file exists else throw error
        if(not os.path.exists(self.catalog_file)) :
            raise myError("Catalog file " + self.catalog_file +
                " does not exist.")

        #Check catalog file exists else throw error
        if(not os.path.exists(self.fits_file)) :
            raise myError("Fits file" + self.fits_file +
                " does not exist.")

        #Make sure lower and upper limits of bins are consistent
        if(self.bin_limits[0] >= self.bin_limits[1]) :
            raise myError("Provided bin input for lower limit greater than higher limit.")

    ####################################################################
    # Method to compute TPCF for the galaxy
    ####################################################################

    def Compute_TPCF(self,cluster_class=-1,save=False,
        random_method='masked_radial',verbose=False):
        """
        Parameters
            filename : string
              name of file from which to read galaxy description
            cluster_class : integer
              class of clusters for which TPCF to be computed
              can be 1,2,3 or -1 for 1+2+3 combined
            save : boolean
              flag to save bootstrap realisations of the TPCF
            random_method : string
                method to use to prepare the random catalog
            verbose : boolean
                print out what is being done

        Returns
            None

        """
        if(self.catalog_file is None) :
            raise myError("No catalog file defined for galaxy.")

        if(verbose):
            if(cluster_class==-1):
                print("Computing TPCF for {} for cluster class 1+2+3 using {} random method."
                .format(self.name,random_method))
            else:
                print("Computing TPCF for {} for cluster class {} using {} random method."
                    .format(self.name,cluster_class,random_method))
        
        if(verbose):
            print("Reading cluster positions from cluster catalog in {}"
                .format(self.catalog_file))

        self.get_ra_dec(cluster_class=cluster_class)
        self.set_bins()

        # Safety check for masked_radial method
        if random_method in ['masked_radial','masked']:
            self.region_file = os.path.abspath('../Data/Region_Files/') + '/' +\
                self.name + '.reg'

            #If region file not present, create it
            if(not os.path.exists(self.region_file)):
                self.create_region_file(verbose=verbose)
            else :
                if(verbose):
                    print("Region file exists. Using region file {}".format(
                        self.region_file))

        corr,dcorr,bootstraps = bootstrap_two_point_angular(self,
                            method='landy-szalay',Nbootstraps=100,
                            random_method=random_method)
        if(verbose):
            print("TPCF computation completed.")
        
        self.corr = corr 
        self.dcorr = dcorr
        self.bootstraps = bootstraps


        if(save):
            saveObj(self.bins,self.outdir+'tpcf_bins')
            saveObj(corr,self.outdir+'tpcf_corr')
            saveObj(dcorr,self.outdir+'tpcf_dcorr')
            saveObj(bootstraps,self.outdir+'tpcf_bootstraps')

    ####################################################################
    # Method to obtain ra dec of clusters
    ####################################################################
    def get_ra_dec(self,cluster_class=-1,verbose=False):
        file = np.loadtxt(self.catalog_file)
        Class0_sources = np.where(file[:,33]==0)
        Class1_sources = np.where(file[:,33]==1)
        Class2_sources = np.where(file[:,33]==2)
        Class3_sources = np.where(file[:,33]==3)
        Class4_sources = np.where(file[:,33]==4)
        Cluster_sources = np.append(Class1_sources,Class2_sources)
        Cluster_sources = np.append(Class3_sources,Cluster_sources)
        
        # Compute TPCF for a subset of clusters if required
        if(cluster_class == 1) :
            RA = file[Class1_sources][:,3]
            DEC = file[Class1_sources][:,4]
        elif(cluster_class == 2) :
            RA = file[Class2_sources][:,3]
            DEC = file[Class2_sources][:,4]
        elif(cluster_class == 3) :
            RA = file[Class3_sources][:,3]
            DEC = file[Class3_sources][:,4]
        elif(cluster_class == -1) :
            RA = file[Cluster_sources][:,3]
            DEC = file[Cluster_sources][:,4]
        else :
            raise myError("Invalid cluster class passed.")

        self.ra = RA 
        self.dec = DEC

        self.ra_raw = RA 
        self.dec_raw = DEC
        if(self.deproject_galaxy == True):
            self.correct_inclination(force=True,verbose=verbose)


    ####################################################################
    # Method to obtain set bins.
    ####################################################################
    def set_bins(self,set_no_bins=None,set_bin_limits=None):
        """
        Method to set the number of bins, allowing user to change no of bins
        and bin limits, in which case bins recomputed. 
        Parameters
            set_no_bins : integer
              Number of bins to set
            set_bin_limits : tuple
              Bin limits to set in parsec
        Returns
            None

        """
        
        bin_limits = [0,0]
        if(set_no_bins == None and set_bin_limits == None):
            #Convert bin limits in parsec to arcsec
            distance = self.distance*const.Parsec*1.e6
            bin_limits[0]= self.bin_limits[0]*const.Parsec/distance*u.radian.to(u.arcsec)
            bin_limits[1] = self.bin_limits[1]*const.Parsec/distance*u.radian.to(u.arcsec)

            bin_min,bin_max = np.log10(bin_limits[0]*u.arcsec.to(u.deg)),\
            np.log10(bin_limits[1]*u.arcsec.to(u.deg))
            bins = np.logspace(bin_min,bin_max,self.no_bins+1)
            self.bins = bins
            self.bin_centres = (self.bins[1:]+self.bins[:-1])/2
            self.bins_arcsec = bins*(1./arcsec_to_degree)

        else :
            if(set_no_bins is not None):
                set_no_bins = int(set_no_bins)
                self.no_bins = no_bins
                print("Changing number of bins to {}".format(set_no_bins))
            if(set_bin_limits is not None):
                set_bin_limits = float(set_bin_limits)
                self.bin_limits = set_bin_limits
                print("Changing bin limits to {}".format(set_bin_limits))

            print("Recomputing bins")
            self.no_bins = set_no_bins 
            #Convert bin limits in parsec to arcsec
            distance = self.distance*const.Parsec*1.e6
            bin_limits[0]= self.bin_limits[0]*const.Parsec/distance*u.radian.to(u.arcsec)
            bin_limits[1] = self.bin_limits[1]*const.Parsec/distance*u.radian.to(u.arcsec)

            bin_min,bin_max = np.log10(bin_limits[0]*u.arcsec.to(u.deg)),\
            np.log10(bin_limits[1]*u.arcsec.to(u.deg))
            bins = np.logspace(bin_min,bin_max,self.no_bins+1)
            self.bins = bins
            self.bin_centres = (self.bins[1:]+self.bins[:-1])/2
            self.bins_arcsec = bins*(1./arcsec_to_degree)



    ####################################################################
    # Method to correct for inclination.
    ####################################################################
    def correct_inclination(self,force=False,verbose=True):
        if(force == False):
            if(self.inclination_done == True):
                raise myError("Inclination correction already performed." +
                    " Use flag force =True if you want to force deproject.")
        if(verbose):
            print("Correcting for inclination of galaxy. This galaxy has PA = {}"
                .format(self.pa) + " and inclination = {}".format(self.inclination))

        self.ra_raw = self.ra
        self.dec_raw = self.dec
        # See Eq 1 Grasha et al 2017
        ra_dep = self.ra*np.cos(self.pa) + self.dec*np.sin(self.pa)
        dec_dep = -self.ra*np.sin(self.pa) + self.dec*np.cos(self.pa)
        dec_dep = dec_dep/(np.cos(self.inclination))

        self.ra = ra_dep 
        self.dec = dec_dep

        self.inclination_done = True

    ####################################################################
    # Method to create ds9 region file for FOV footprint
    ####################################################################

    def create_region_file(self,verbose=False):
        """
        Creates region file for galaxy. 
        Parameters
            verbose: Boolean
                Print out information

        Returns
            None

        """
        from footprintfinder import main
        fits_file = self.fits_file
        if(verbose):
            print("Creating region file....")
        #Call footprint finder
        main('-d ' + os.path.abspath(self.fits_file))

        #Move region file to appropriate directory
        region_filebase = self.fits_file.split(os.path.dirname(self.fits_file)+'/')[1]
        region_filebase = region_filebase.split('.fits')[0]
        region_filename = region_filebase + '_footprint_ds9_image.reg'

        # Copy created region file to directory
        self.region_file = os.path.abspath('../Data/Region_Files/' + 
            self.name + '.reg')

        shutil.copy(os.getcwd() + '/' + region_filename,self.region_file)

        if(verbose):
            print("Saved region file in {}".format(self.region_file))

        #Remove other output files 
        if(verbose):
            print("Deleting files created by footprint in local directory.")
        unwanted_file = region_filebase + '_footprint.txt'
        subprocess.run(["rm",os.getcwd()+'/{}'.format(unwanted_file)],
            stdout=subprocess.DEVNULL)
        unwanted_file = region_filebase + '_footprint_ds9_linear.reg'
        subprocess.run(["rm",os.getcwd()+'/{}'.format(unwanted_file)],
            stdout=subprocess.DEVNULL)
        unwanted_file = region_filename
        subprocess.run(["rm",os.getcwd()+'/{}'.format(unwanted_file)],
            stdout=subprocess.DEVNULL)


    ####################################################################
    # Method to fit power law to TPCF
    #TODO: MCMC fit. 
    ####################################################################

    def fit_power_law(self,method='bootstrap',N=1000):
        """
        Parameters
            filename : string
              name of file from which to read galaxy description
            cluster_class : integer
              class of clusters for which TPCF to be computed
              can be 1,2,3 or -1 for 1+2+3 combined
            save : boolean
              flag to save bootstrap realisations of the TPCF

        Returns
            None

        """
        if(method not in ['single','bootstrap','mcmc']) :
            raise ValueError("Method for fitting power law should be one of the"+
            " following: single, bootstrap or mcmc")
        
        bins = self.bins_arcsec
        separation_bins = (bins[1:]+bins[:-1])/2
        
        separation_bins = separation_bins.astype(np.float)
        corr_fit = self.corr[np.where(self.corr>0.0)].astype(np.float)
        dcorr_fit = self.dcorr[np.where(self.corr>0.0)].astype(np.float)
        separation_bins = separation_bins[np.where(self.corr>0.0)].astype(np.float)


        if(method == 'single') :
            popt,pcov = curve_fit(linear_function,separation_bins,
                np.log(corr_fit),sigma=np.log(dcorr_fit))
                  
            self.fit_values = popt
            self.fit_errors = np.sqrt(np.diag(pcov))

        elif(method == 'bootstrap') :
            fit_bootstraps = np.zeros((N,4))
            for i in range(N):
                y_fit = corr_fit + dcorr_fit*np.random.randn(1)

                try:
                #Fit to this
                    popt,pcov = curve_fit(linear_function,separation_bins,
                        np.log(y_fit))
                except :
                    continue

                fit_bootstraps[i] = popt

            self.fit_values = np.median(fit_bootstraps,axis=0)
            self.fit_errors = np.std(fit_bootstraps,axis=0)
            self.fit_distribution = fit_bootstraps

    # def read_summary(self,filename=None,method='masked_radial'):
    #     """
    #     Read pickle file where the class is saved as an object. 
    #     Parameters:
    #         filename: string
    #             Summary file to read
    #         method: string
    #             Random catalog method for correct output directory

        
    #     """
    #     if(filename is None):
    #         if(method.upper() == 'MASKED'):
    #             method_dir = 'Masked/'
    #         elif(method.upper() == 'UNIFORM'):
    #             method_dir = 'Uniform/'
    #         elif(method.upper() == 'MASKED_RADIAL'):
    #             method_dir = 'Masked_Radial/'
    #         else:
    #             raise myError("Method not recognised.")

    #         filename = self.outdir+'/'+ method_dir+'{}_summary.pkl'.format(self.name)


    #     if(os.path.isfile(filename)):

    #         print("Loading class object from {}".format(filename))
    #         file = loadObj(filename.split('.pkl')[0])

            
    #     else :
    #         raise myError("Summary class file {} not present.".format(filename))
            


    

class myError(Exception):
    """
    Class derived from Exception to handle exceptions raised by
    program-specific errors.

    Parameters
       message : string
          the error message
    """

    def __init__(self, message):
        Exception.__init__(self, message)
        self.message = message