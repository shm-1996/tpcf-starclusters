from header import *
from Galaxy import *
from Plot_Class import *


def tpcf_galaxy(galaxy_name,method='masked_radial',outdir=None,overwrite=False):
    """
    Perform all required computations for one galaxy.
    Parameters
        ----------
        galaxy_name : string
            Name of galaxy
        method: string
            Method to prepare random catalog.
        
        Returns
        -------
        galaxy_class : Class Galaxy
            Instance of class galaxy is returned
    """
    print("Computing TPCF for galaxy {}.".format(galaxy_name))
    #Initialise galaxy class
    galaxy_class = Galaxy(galaxy_name,verbose=True)

    #Change output directory if passed by user
    if(outdir is not None):
        galaxy_class.outdir = outdir

    #Save in subdirectories for below two methods
    if(method == 'masked'):
        galaxy_class.outdir += '/Masked/'
    elif(method == 'uniform'):
        galaxy_class.outdir += '/Uniform/'
    elif(method == 'masked_radial'):
        galaxy_class.outdir += '/Masked_Radial/'
    else:
        raise myError("Method not recognised.")
    #Make sure path exists else create
    if(not os.path.exists(galaxy_class.outdir)):
        print("Path for plots {} do not exist. Creating now"
            .format(galaxy_class.outdir))
        os.makedirs(galaxy_class.outdir)

     #First check if pickle file of class already exists here
    if(overwrite is False):
        if(os.path.isfile(galaxy_class.outdir+
            '{}_summary.pkl'.format(galaxy_class.name))):
            
            galaxy_class = loadObj(galaxy_class.outdir+
            '{}_summary'.format(galaxy_class.name))
            return galaxy_class
        else :
            print("TPCF not saved in output directory . Computing TPCF now..")    
    else :
        print("Overwrite flag is provided. Computing and overwriting TPCF.")



    #Compute and fit TPCF for all classes
    #TODO: Implement this

    #Compute TPCF for combined
    galaxy_class.Compute_TPCF(verbose=True,save=True,random_method=method)
    #Fit power law
    galaxy_class.fit_power_law()

    return galaxy_class

def plots_galaxy(pl,method='masked_radial',outdir=None,save=False):
    """
    Plot all the required output plots for each galaxy.
    Parameters
        ----------
        pl : class myPlot
            Instance of class myPlot.
        method: string
            Method to prepare random catalog.
        Returns
        -------
        None
    """

    if(outdir is not None):
        galaxy_class.outdir = outdir

    print("Plotting summary plots for Galaxy {}".format(pl.galaxy.name))
    pl.galaxy_image(save=True)
    pl.plot_clusters(save=save)
    pl.plot_random(save=save,random_method=method)
    pl.class_distribution(save=save)
    pl.mass_histogram(save=save)
    pl.age_histogram(save=save)
    pl.bin_distribution(save=save)
    pl.plot_TPCF(save=save)
    pl.plot_TPCF_allclass(random_method=method,save=save,verbose=True)


def tpcf_allgalaxies(method,overwrite=False,save=False) :
    """
    Plot all the required output plots for all galaxies.
    Parameters
        ----------
        method : string
            Method to prepare random catalog.
        
        Returns
        -------
        None
    """
    print("Computing TPCF for all galaxies in LEGUS Survey.")
    for galaxy_name in list_of_galaxies:
        print("Calculating TPCF for galaxy {}".format(galaxy_name))
        galaxy_class = tpcf_galaxy(galaxy_name,method,overwrite=overwrite)
        plot_class = myPlot(galaxy_class)
        plots_galaxy(plot_class,method,save=save)
        
        #Save galaxy class info as pickle file
        path_exists = os.path.isfile(galaxy_class.outdir+
            '{}_summary.pkl'.format(galaxy_class.name))
        if(overwrite is True or path_exists is False):
            print("Saving class object of {} as pickle file.".format(galaxy_class.name))
            saveObj(galaxy_class,galaxy_class.outdir+'{}_summary'
                .format(galaxy_class.name))
        print("##############################################################")
        print("\n\n")


if __name__ == "__main__":

    #Parsing Arguments
    ############################################################################
    ap = argparse.ArgumentParser(description=
        'Command Line Inputs for tpcf-starclusters. All inputs optional. ')
    ap.add_argument('-method',action='store',type=str,default='masked_radial',
        help='Method to prepare the random catalog: "Uniform","Masked"' +
        '" Masked_radial (default)" ')
    ap.add_argument('-galaxy',action='store',type=str,default=None,
        help = 'Galaxy for which tpcf to be computed. By default done for all.')
    ap.add_argument('-outdir',action='store',type=str,default=None,
        help = 'Alternate output directory for plots and files.')
    ap.add_argument('-overwrite',action='store_true',
        help='overwrite computation of TPCF.Default: False.')
    ap.add_argument('-show_only',action='store_true',
        help='Invoke this flag to only show the plots and not save them.')
    args = vars(ap.parse_args())

    method = args['method'].lower()
    if(method not in ['uniform','masked','masked_radial']):
        raise ValueError("This method does not exist. Allowed values are "+
            "'Uniform', 'Masked', and 'Masked_Radial'.")
    galaxy_name = args['galaxy'].upper()
    if(args['outdir'] is not None):
        output_directory = os.path.abspath(args['outdir']+'/')
        if(not os.path.exists(output_directory)):
            os.makedirs(output_directory)

    else:
        output_directory = None
    #Arguments parsed
    ############################################################################

    #Check if plots to be saved
    if(args['show_only'] is True):
        save=False
    else:
        save=True

    #Only for one galaxy
    if(galaxy_name is not None):
         #for all galaxies
        if(galaxy_name == 'ALL'):
            tpcf_allgalaxies(method,overwrite=args['overwrite'],save=save)
        elif(galaxy_name in list_of_galaxies):

            print("Running tpcf-starclusters for {}.".format(galaxy_name))
            galaxy_class = tpcf_galaxy(galaxy_name,method,output_directory, 
                overwrite=args['overwrite'])
            plot_class = myPlot(galaxy_class)
            plots_galaxy(plot_class,method,output_directory,save=save)
            print("Saving class object of {} as pickle file.".format(galaxy_class.name))
            saveObj(galaxy_class,galaxy_class.outdir+'{}_summary'
            .format(galaxy_class.name))
            #Save galaxy class info as pickle file



        else:
            raise myError("The provided galaxy is not in the list of galaxies"+
                " for which cluster catalogs are available with LEGUS.")

   

    else :
        raise myError("Please provide galaxy name or ALL for all galaxies.")

        #TODO: Prepare galaxy info table
        #galaxy_info(output_directory)

        #TODO: Prepare fit summary table
        #galaxy_fit(output_directory)






    
        
        



