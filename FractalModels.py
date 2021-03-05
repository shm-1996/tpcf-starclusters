from TPCF import linear_function,onepowerlaw_function,linear_truncation
from matplotlib.transforms import Bbox
from ToyModels import *
from joblib import Parallel, delayed
from itertools import repeat
import timeit

#Some fixed parameters
############################################################################################
max_level = 10
max_level_low = 6
max_level_high = 14
fractal_dims = [0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]
sizes = [0.1,0.2,0.3,0.4,0.5,0.6]
filename = '../Toy_Models/Fractals/Fractal_Analysis.out'

############################################################################################

# Toy Model 4 : Fractal Models

def Create_Fractal(fractal_dim,baselevel=2,max_level=10):

    prob = 2**(fractal_dim-2)
    delta = 1./(2**baselevel)
    # xfractal, yfractal = np.meshgrid(np.arange(delta/2.,1.0,
    #                         delta),np.arange(delta/2.,1.0,delta),indexing='xy')
    # xfractal,yfractal = xfractal.flatten(),yfractal.flatten()
    delta = 1./(2**(baselevel+1))
    current_level_x, current_level_y = np.meshgrid(np.arange(delta/2.,1.0,
                            delta),np.arange(delta/2.,1.0,delta),indexing='xy')
    current_level_x, current_level_y = current_level_x.flatten(),current_level_y.flatten()
    for level in range(baselevel+1,max_level):
        nactive = np.size(current_level_x)
        rand_no = np.random.uniform(size=nactive)
        new_level_x = []
        new_level_y = []
        #Loop over active regions and add regions if condition
        #satisfied 
        for region in range(0,nactive):
            if(rand_no[region]<prob): 
                #Add to fractal points
                # xfractal = np.append(current_level_x[region],xfractal)
                # yfractal = np.append(current_level_y[region],yfractal)

                #Add to next level active regions to loop over


                #Create 4 new centres for this region
                delta = 1./(2**(level+1))/2.
                new_x = [current_level_x[region]-delta,
                        current_level_x[region]+delta,
                        current_level_x[region]-delta,
                        current_level_x[region]+delta]
                new_y = [current_level_y[region]+delta,
                        current_level_y[region]+delta,
                        current_level_y[region]-delta,
                        current_level_y[region]-delta]
                new_level_x = np.append(new_x,new_level_x)
                new_level_y = np.append(new_y,new_level_y)

        current_level_x = new_level_x
        current_level_y = new_level_y 
        
    
    return new_level_x,new_level_y    

def Compute_TPCF_Fractal(xfractal,yfractal,no_bins=40,limits=0.0):
    
    indices_1 = np.logical_and(xfractal>limits,yfractal>limits)
    indices_2 = np.logical_and(xfractal<1.-limits,yfractal<1.-limits)
    indices = np.logical_and(indices_1,indices_2)
    size = (1.-limits*2.0)/4.
    bins = np.logspace(np.log10(0.001),np.log10(size),no_bins+1)
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
    start = time.time()
    bins,corr_lz,dcorr_lz = Compute_TPCF_Fractal(xfractal,yfractal,no_bins=no_bins)
    end = time.time()
    print("TPCF computed.Time taken: {} seconds".format(end-start))
    save_obj = [xfractal,yfractal,bins,corr_lz,dcorr_lz]
    saveObj(save_obj,'../Toy_Models/Fractals/D_{}'.format(int(fractal_dim*10)))

    #Plot Positions
    plt.clf()
    plt.scatter(xfractal,yfractal,s=0.05,c='blue')
    plt.savefig('../Toy_Models/Fractals/Pos_D_{}'.format(int(fractal_dim*10)),
                bbox_inches='tight')
    plt.close()

    #Fit and Plot TPCF
    separation_bins = (bins[1:]+bins[:-1])/2.
    indices = np.where(corr_lz>0.0)
    dcorr_lz = dcorr_lz[indices]
    separation_bins = separation_bins[indices]
    corr_lz = corr_lz[indices]
    
    
    plt.errorbar(separation_bins,1+corr_lz,yerr=dcorr_lz,
                 fmt='.-',lw=0.2)
    #Fit line to this
    popt,pcov = curve_fit(onepowerlaw_function,separation_bins,
                        np.log(1+corr_lz),sigma=dcorr_lz/corr_lz)
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
    plt.savefig('../Toy_Models/Fractals/TPCF_D_{}'.format(int(fractal_dim*10)),
                bbox_inches='tight')
    plt.clf()
    plt.close()

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

    size = (1.-limits)/2.
    axs.axvline(size/2.,color='#F9004A',lw=1.5,
                label=r'$r_{\mathrm{max}} =$'+r'${}$'.format(size),ls='--')
    axs.axvline(popt[2],ls='--',color='#F4271C',lw=1.5,
                  label=r'$r_c = {:3.2f} \pm {:3.2f}$'.format(popt[2],perr[2]))

    axs.set_xlabel(r"$\Delta x \, (\mathrm{pixels})$")
    axs.set_ylabel(r"$\omega_{\mathrm{LS}}\left(\theta \right)$")
    axs.set_xscale('log')
    axs.set_yscale('log')
    x0,x1 = axs.get_xlim()
    axs.set_xlim(x0,size)
    
    axs.legend(loc='best')
    plt.savefig('../Toy_Models/Fractals/TPCF_cutoff_{}'.format(int(fractal_dim*10)),bbox_inches='tight')        
    plt.clf()
    plt.close(fig)



#############################################################################################
# Toy Model 5 : Fractals with Artificial Boundaries to test boundary effects

def check_contains(x,y,bbox):
    return bbox.contains(x,y)

def MC_fractal(xfractal,yfractal,size,angle=0.0,no_bins=40,N_MC=30,debug=False):

    n = 0
    attempts = 0
    max_attempts = 5
    vfunc = np.vectorize(check_contains,excluded=['bbox'])
    corr_MC = np.zeros((N_MC,no_bins))
    dcorr_MC = np.zeros((N_MC,no_bins))
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

        #Compute TPCF 
        bins = np.logspace(np.log10(0.001),np.log10(size/4.),no_bins+1)
        data = np.asarray((xfractal[indices_true],yfractal[indices_true]),order='F').T

        try:
            corr_lz,dcorr_lz = bootstrap_two_point(data, bins, Nbootstrap=100,
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
        attempts = 0
    return bins,corr_MC,dcorr_MC

def TPCF_MC(xfractal,yfractal,size,angle=0.0,no_bins=40,N_MC=30,MC_directory=None):
    
    if(MC_directory == None):
        MC_directory = '../Toy_Models/Fractals/Monte_Carlo_Boundaries'
    MC_directory = os.path.abspath(MC_directory)

    if(not os.path.exists(MC_directory)):
            os.makedirs(MC_directory)

    bins,corr_MC,dcorr_MC = MC_fractal(xfractal,yfractal,size,angle,no_bins,N_MC)
    corr_avg = np.mean(corr_MC,axis=0)
    dcorr_avg = np.mean(dcorr_MC,axis=0)

    pkl_obj = [bins,corr_avg,dcorr_avg]
    saveObj(pkl_obj,MC_directory+'/Size_{}'.format(int(size*10)))
    
    #Landy-Szalay version
    separation_bins = (bins[1:]+bins[:-1])/2.
    indices = np.where(corr_avg>0.0)
    dcorr_lz = dcorr_avg[indices]
    separation_bins = separation_bins[indices]
    corr_lz = corr_avg[indices]
    

    #Fit line to this for fractal dimension
    popt,pcov = curve_fit(onepowerlaw_function,separation_bins,
                        np.log(1+corr_lz),sigma=dcorr_lz/(corr_lz))
    perr= np.sqrt(np.diag(pcov))
    
    fig,axs = plt.subplots(ncols=1)

    axs.errorbar(separation_bins,1+corr_lz,yerr=dcorr_lz,
                 fmt='.-',lw=0.2)
    axs.plot(separation_bins,np.exp(onepowerlaw_function(separation_bins,
                                                        popt[0],popt[1])),
            ls='-',lw=0.2,c='k',
             label=r'$\alpha = {:2.1f} \pm {:3.2f}$'.format(popt[1],perr[1]))

    axs.axvline(size/2.,color='#F9004A',lw=1.5,
                label=r'$r_{\mathrm{max}}$',ls='--')

    axs.set_xlabel(r"$\Delta x \, (\mathrm{pixels})$")
    axs.set_ylabel(r"$1+ \omega_{\mathrm{LS}}\left(\theta \right)$")
    axs.set_xscale('log')
    axs.set_yscale('log')
    x0,x1 = axs.get_xlim()
    axs.set_xlim(x0,size)
    
    axs.legend(loc='best')
    plt.savefig(MC_directory+'/Size_{}_Fractal'.format(int(size*10)),bbox_inches='tight')        
    plt.clf()
    plt.close(fig)


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

    axs.axvline(size/2.,color='#F9004A',lw=1.5,
                label=r'$r_{\mathrm{max}}$'+r'$ = {}$'.format(size/2.),ls='--')
    axs.axvline(popt[2],ls='--',color='#F4271C',lw=1.5,
                  label=r'$r_c = {:3.2f} \pm {:3.2f}$'.format(popt[2],perr[2]))

    axs.set_xlabel(r"$\Delta x \, (\mathrm{pixels})$")
    axs.set_ylabel(r"$\omega_{\mathrm{LS}}\left(\theta \right)$")
    axs.set_xscale('log')
    axs.set_yscale('log')
    x0,x1 = axs.get_xlim()
    axs.set_xlim(x0,size)
    
    axs.legend(loc='best')
    plt.savefig(MC_directory+'/Size_{}_Cutoff'.format(int(size*10)),bbox_inches='tight')        
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
if __name__ == "__main__":
    #Fractal_Table_OneD(1.5,overwrite=False)
    
    # time the script
    start_time = timeit.default_timer()
    Analyse_Fractal(1.5,max_level=14)
    # time the script
    stop_time = timeit.default_timer()
    total_time = stop_time - start_time
    print("***************** time to finish = "+str(total_time)+"s *****************")
    print("Done")