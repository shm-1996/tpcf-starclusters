"""
Routines to make summary tables in Latex format.

AUTHOR (modified)
Shyam Harimohan Menon (2020)
"""

from header import *
from Galaxy import *

def GalaxyInfo_Table():

    """
    This method writes a summary latex table of galaxy properties.
    
    Parameters
    ----------
    results_directory : string
        directory onto which to save the summary table. 
    
    Returns
    -------
    None
    """

    print("Making summarised latex table.")
    print("####################################")
    results_directory = '../Results/'
    filename = results_directory+'Summary_Galaxy.tex'

    


    with open(filename,"w") as file :
        #Header and table columns
        print("\\documentclass[a4paper]{article}",file=file)
        print("\\usepackage{booktabs}",file=file)
        print("\\usepackage{fouriernc}",file=file)
        print("\\usepackage[flushleft]{threeparttable}",file=file)
        print("\\renewcommand{\\arraystretch}{1.7}",file=file)
        print("",file=file)
        print("\\begin{document}",file=file)
        print("",file=file)
        print("",file=file)
        print("\\begin{table*}",file=file)
        print("\\centering",file=file)
        print("\\label{tab:galaxyinfo}",file=file)
        print("\\begin{threeparttable}",file=file)
        print("\\caption{Summarised physical quantities of the galaxies in this study.}",file=file) 
        print("\\begin{tabular}{l c c c c c c c c c}",file=file)
        print("\\toprule",file=file)
        #print("Pillar (1)& 2 & 4 & 8 & 20 & 22 & 44\\\\ % Column names row",file=file)
        print("\\multicolumn{1}{l}{Name}& \\multicolumn{1}{c}{Morph.}& \\multicolumn{1}{c}{$T$} &\
 \\multicolumn{1}{c}{Inclin.} & \\multicolumn{1}{c}{P.A} & \\multicolumn{1}{c}{Dist.} & \\multicolumn{1}{c}{$\\mathrm{SFR}_{UV}$} &\
 \\multicolumn{1}{c}{$M_{*}$} & \\multicolumn{1}{c}{$R_{25}$} & \\multicolumn{1}{c}{$\\Sigma_{\\mathrm{SFR}}$} \\\\ % Column names row",file=file)
        print("\\midrule",file=file)
        i = 0

        for galaxy_name in list_of_galaxies:
            galaxy_class = loadObj("../Results/Galaxies/{}/Masked/{}_summary".format(galaxy_name,galaxy_name))
            galaxy_class.read_galaxyprops()

            mstar = galaxy_class.mstar
            mstar = "{:2.1e}".format(mstar)
            mstarpower = int(mstar.split('e')[1]) 
            mstar = float(mstar.split('e')[0])

            sigma_sfr = galaxy_class.sigma_sfr
            sigma_sfr = "{:2.1e}".format(sigma_sfr)
            sigma_sfr_power = int(sigma_sfr.split('e')[1])
            sigma_sfr = float(sigma_sfr.split('e')[0])

            #Write string            
            row_str = ""
            name = galaxy_class.name.split('_')[0] + ' ' + galaxy_class.name.split('_')[1]
            row_str +="{} &".format(name)

            row_str +="{} &".format(galaxy_class.morph_type)
            row_str +="${} $&".format(galaxy_class.T_value)
            row_str +="{:2.1f} &".format(galaxy_class.inclination)
            row_str +="{:2.1f} &".format(galaxy_class.pa)
            row_str +="${:2.1f} $&".format(galaxy_class.distance)
            row_str +="${:3.2f} $&".format(galaxy_class.sfr)
            row_str +="${:2.1f} \\times 10^{{{}}}$&".format(mstar,mstarpower)
            row_str +="${:2.1f} $&".format(galaxy_class.r25)
            row_str +="${:2.1f} \\times 10^{{{}}}$".format(sigma_sfr,sigma_sfr_power)
            row_str +="\\\\"
            print(row_str,file=file)

        print("\\bottomrule",file=file)
        print("\\end{tabular}",file=file)
        print("\\begin{tablenotes}",file=file)
        print("\\small",file=file)
        print("\\item \\textbf{Notes}: Values in each field adopted from Calzetti et al 2015, unless otherwise specified.",
            file=file)
        print("References for position angles: Leroy et al 2019, Kirby et al 2008, Ho et al 2011, Verdes-Montenegro et al 2000,",file=file)
        print(" Vaduvescu et al 2005, SDSS DR5, Swaters et al 2002, Schechtman-Rook et al 2012, Hu et al 2012, ",file=file)
        print("SDSS DR5, Springbob et al 2007, Kirby et al 2008.",file=file)

        print("\\end{tablenotes}",file=file)      
        print("\\end{threeparttable}",file=file) 
        
        print("\\end{table*}",file=file)
        print("\\end{document}",file=file)

    print("Printed table to {}".format(filename))
    print("Creating PDF with pdflatex")
    subprocess.run(["pdflatex","-output-directory={}".format(results_directory),"{}".format(filename)], stdout=subprocess.DEVNULL)
    filename_pdf = os.path.splitext(filename)[0]+'.pdf'
    subprocess.run(["pdfcrop","{}".format(filename_pdf),"{}".format(filename_pdf)], stdout=subprocess.DEVNULL)
    subprocess.run(["rm","{}Summary_Galaxy.aux".format(results_directory)], stdout=subprocess.DEVNULL)
    subprocess.run(["rm","{}Summary_Galaxy.log".format(results_directory)], stdout=subprocess.DEVNULL)
    subprocess.run(["rm","-rf","{}log/".format(results_directory)], stdout=subprocess.DEVNULL)
    return


def AIC_Table():
    print("Making AIC summary table.")
    print("####################################")
    results_directory = '../Results/'
    filename = results_directory+'AIC.tex'

    with open(filename,"w") as file :
        #Header and table columns
        print("\\documentclass[a4paper]{article}",file=file)
        print("\\usepackage{booktabs}",file=file)
        print("\\usepackage{fouriernc}",file=file)
        print("\\usepackage[flushleft]{threeparttable}",file=file)
        print("\\renewcommand{\\arraystretch}{1.7}",file=file)
        print("",file=file)
        print("\\begin{document}",file=file)
        print("",file=file)
        print("",file=file)
        print("\\begin{table*}",file=file)
        print("\\centering",file=file)
        print("\\label{tab:galaxyinfo}",file=file)
        print("\\begin{threeparttable}",file=file)
        print("\\caption{MCMC best fit parameters and associated Akaike Information Criterion (AIC) values.}",file=file) 
        print("\\begin{tabular}{l c c c c}",file=file)
        print("\\toprule",file=file)    
        print("\\multicolumn{1}{l}{Galaxy}& \\multicolumn{1}{c}{$\\mathrm{AIC}_{\\mathrm{PW}}$} &\
\\multicolumn{1}{c}{$\\mathrm{AIC}_{\\mathrm{S}}$} & \\multicolumn{1}{c}{$\\mathrm{AIC}_{\\mathrm{ST}}$}\
 &\\multicolumn{1}{c}{Model} \\\\ % Column names row",file=file)
        #print(" &$\\alpha_1$&$\\alpha_2$&$\\beta$ & & $\\alpha{\\mathrm{S}}$ & & $\\alpha_{\\mathrm{ST}}$& $\\theta_c$ & & \\\\",file=file)
        print("\\midrule",file=file)

        i = 0

        for galaxy_name in list_of_galaxies:
            galaxy_dir = "../Results/Galaxies/{}/Masked".format(galaxy_name)
            galaxy_class = loadObj("{}/{}_summary".format(galaxy_dir,galaxy_name))

            #Read in PW fits
            sampler = loadObj(galaxy_dir+'/PiecePL_MCMC/'+'MCMC_sampler')
            samples = sampler.flatchain

            # PW_median_fit = np.percentile(samples,50,axis=0)
            # PW_error_low = np.percentile(samples,16,axis=0)
            # PW_error_low = PW_median_fit-PW_error_low
            # PW_error_high = np.percentile(samples,84,axis=0)
            # PW_error_high = PW_error_high - PW_median_fit

            # #Account for fact result is log(beta)
            # PW_error_low[3] = np.exp(PW_error_low[3])
            # PW_error_high[3] = np.exp(PW_error_high[3])
            # PW_median_fit[3] = np.exp(PW_median_fit[3])

            #Read in Single PL fit
            sampler = loadObj(galaxy_dir+'/SinglePL_MCMC/'+'MCMC_sampler')
            samples = sampler.flatchain

            # S_median_fit = np.percentile(samples,50,axis=0)
            # S_error_low = np.percentile(samples,16,axis=0)
            # S_error_low = S_median_fit-S_error_low
            # S_error_high = np.percentile(samples,84,axis=0)
            # S_error_high = S_error_high - S_median_fit

            #Read in SinglePL + Truncation fit
            sampler = loadObj(galaxy_dir+'/SingleTrunc_MCMC/'+'MCMC_sampler')
            samples = sampler.flatchain

            # ST_median_fit = np.percentile(samples,50,axis=0)
            # ST_error_low = np.percentile(samples,16,axis=0)
            # ST_error_low = ST_median_fit-ST_error_low
            # ST_error_high = np.percentile(samples,84,axis=0)
            # ST_error_high = ST_error_high - ST_median_fit
            
            #Get AIC values and find best model
            AIC_S, AIC_PW, AIC_ST = compare_AIC(galaxy_name)
            AIC_Values = [AIC_PW, AIC_S, AIC_ST]
            Model = ['PW', 'S', 'ST']
            best_index = np.argmin(AIC_Values)


            #Write to table
            row_str = ""
            name = galaxy_class.name.split('_')[0] + ' ' + galaxy_class.name.split('_')[1]
            row_str +="{} &".format(name)

            # for j in range(1,4):
            #     row_str +="${:2.1f}^{{+{:3.2f}}}_{{-{:3.2f}}} $&".format(PW_median_fit[j],PW_error_high[j],PW_error_low[j])

            row_str +="${:2.1f} $&".format(AIC_PW)

            #row_str +="${:2.1f}^{{+{:3.2f}}}_{{-{:3.2f}}} $&".format(S_median_fit[1],S_error_high[1],S_error_high[1])

            row_str +="${:2.1f} $ &".format(AIC_S)
            # for j in range(1,3):
            #     row_str +="${:2.1f}^{{+{:3.2f}}}_{{-{:3.2f}}} $&".format(ST_median_fit[j],ST_error_high[j],ST_error_high[j])
            row_str +="${:2.1f} $ &".format(AIC_ST)
            row_str +="{}".format(Model[best_index])
            row_str +="\\\\"

            print(row_str,file=file)
        
        print("\\bottomrule",file=file)
        print("\\end{tabular}",file=file)
        print("\\begin{tablenotes}",file=file)
        print("\\small",file=file)
        print("\\item \\textbf{Notes}: Models compared: Piecewise power-law (PW), Single power-law (S) and Single power-law with exponential truncation (ST).",
            file=file)
        print("\\end{tablenotes}",file=file)      
        print("\\end{threeparttable}",file=file)
        print("\\end{table*}",file=file)
        print("\\end{document}",file=file)

        print("Printed table to {}".format(filename))
        print("Creating PDF with pdflatex")
        subprocess.run(["pdflatex","-output-directory={}".format(results_directory),"{}".format(filename)], stdout=subprocess.DEVNULL)
        filename_pdf = os.path.splitext(filename)[0]+'.pdf'
        subprocess.run(["pdfcrop","{}".format(filename_pdf),"{}".format(filename_pdf)], stdout=subprocess.DEVNULL)
        subprocess.run(["rm","{}Summary_Galaxy.aux".format(results_directory)], stdout=subprocess.DEVNULL)
        subprocess.run(["rm","{}Summary_Galaxy.log".format(results_directory)], stdout=subprocess.DEVNULL)
        subprocess.run(["rm","-rf","{}log/".format(results_directory)], stdout=subprocess.DEVNULL)
        return

def GalaxyTPCF_Table():

    return

def compare_AIC(galaxy_name):
    #Read in samplers
    
    sampler_single = loadObj('../Results/Galaxies/{}/Masked/\
SinglePL_MCMC/MCMC_sampler'.format(galaxy_name))

    sampler_piecewise = loadObj('../Results/Galaxies/{}/Masked/\
PiecePL_MCMC/MCMC_sampler'.format(galaxy_name))

    sampler_singletrunc = loadObj('../Results/Galaxies/{}/Masked/\
SingleTrunc_MCMC/MCMC_sampler'.format(galaxy_name))
    
    AIC_single = 2*2.0 - np.max(sampler_single.flatlnprobability)
    AIC_piecewise = 2*4.0 - np.max(sampler_piecewise.flatlnprobability)
    AIC_singletrunc = 2*3.0 - np.max(sampler_singletrunc.flatlnprobability)
    
    return AIC_single,AIC_piecewise, AIC_singletrunc
