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
        print("\\begin{tabular}{l c c c c c}",file=file)
        print("\\toprule",file=file)    
        print("\\multicolumn{1}{l}{Galaxy}& \\multicolumn{1}{c}{$\\mathrm{AIC}_{\\mathrm{PW}}$} &\
\\multicolumn{1}{c}{$\\mathrm{AIC}_{\\mathrm{S}}$} & \\multicolumn{1}{c}{$\\mathrm{AIC}_{\\mathrm{ST}}$}\
 &\\multicolumn{1}{c}{$\\mathrm{AIC}_{\\mathrm{PWT}}$} &\\multicolumn{1}{c}{Model} \\\\ % Column names row",file=file)
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
            AIC_S, AIC_PW, AIC_ST, AIC_PWT = compare_AIC(galaxy_name)
            AIC_Values = [AIC_PW, AIC_S, AIC_ST,AIC_PWT]
            Model = ['PW', 'S', 'ST','PWT']
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
            row_str +="${:2.1f} $ &".format(AIC_PWT)
            row_str +="{}".format(Model[best_index])
            row_str +="\\\\"

            print(row_str,file=file)
        
        print("\\bottomrule",file=file)
        print("\\end{tabular}",file=file)
        print("\\begin{tablenotes}",file=file)
        print("\\small",file=file)
        print("\\item \\textbf{Notes}: Models compared: Piecewise power-law (PW), Single power-law (S),Single power-law with exponential truncation (ST), and Piecewise power-law with exponential truncation (PWT).",
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

def BIC_Table():
    print("Making BIC summary table.")
    print("####################################")
    results_directory = '../Results/'
    filename = results_directory+'BIC.tex'

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
        print("\\caption{MCMC best fit parameters and associated Bayesian Information Criterion (BIC) values.}",file=file) 
        print("\\begin{tabular}{l c c c c c}",file=file)
        print("\\toprule",file=file)    
        print("\\multicolumn{1}{l}{Galaxy}& \\multicolumn{1}{c}{$\\mathrm{BIC}_{\\mathrm{PW}}$} &\
\\multicolumn{1}{c}{$\\mathrm{BIC}_{\\mathrm{S}}$} & \\multicolumn{1}{c}{$\\mathrm{BIC}_{\\mathrm{ST}}$}\
 &\\multicolumn{1}{c}{$\\mathrm{BIC}_{\\mathrm{PWT}}$} &\\multicolumn{1}{c}{Model} \\\\ % Column names row",file=file)
        #print(" &$\\alpha_1$&$\\alpha_2$&$\\beta$ & & $\\alpha{\\mathrm{S}}$ & & $\\alpha_{\\mathrm{ST}}$& $\\theta_c$ & & \\\\",file=file)
        print("\\midrule",file=file)

        i = 0

        for galaxy_name in list_of_galaxies:
            galaxy_dir = "../Results/Galaxies/{}/Masked".format(galaxy_name)
            galaxy_class = loadObj("{}/{}_summary".format(galaxy_dir,galaxy_name))

            #Get number of samples
            corr_fit = galaxy_class.corr[np.where(galaxy_class.corr>0.0)].astype(np.float)
            nsamples = np.size(corr_fit)

            #Read in PW fits
            sampler = loadObj(galaxy_dir+'/PiecePL_MCMC/'+'MCMC_sampler')
            samples = sampler.flatchain

            #Read in Single PL fit
            sampler = loadObj(galaxy_dir+'/SinglePL_MCMC/'+'MCMC_sampler')
            samples = sampler.flatchain

            
            #Read in SinglePL + Truncation fit
            sampler = loadObj(galaxy_dir+'/SingleTrunc_MCMC/'+'MCMC_sampler')
            samples = sampler.flatchain

            #Get BIC values and find best model
            BIC_S, BIC_PW, BIC_ST, BIC_PWT = compare_BIC(galaxy_name,nsamples)
            BIC_Values = [BIC_PW, BIC_S, BIC_ST,BIC_PWT]
            Model = ['PW', 'S', 'ST','PWT']
            best_index = np.argmin(BIC_Values)


            #Write to table
            row_str = ""
            name = galaxy_class.name.split('_')[0] + ' ' + galaxy_class.name.split('_')[1]
            row_str +="{} &".format(name)

            # for j in range(1,4):
            #     row_str +="${:2.1f}^{{+{:3.2f}}}_{{-{:3.2f}}} $&".format(PW_median_fit[j],PW_error_high[j],PW_error_low[j])

            row_str +="${:2.1f} $&".format(BIC_PW)

            #row_str +="${:2.1f}^{{+{:3.2f}}}_{{-{:3.2f}}} $&".format(S_median_fit[1],S_error_high[1],S_error_high[1])

            row_str +="${:2.1f} $ &".format(BIC_S)
            # for j in range(1,3):
            #     row_str +="${:2.1f}^{{+{:3.2f}}}_{{-{:3.2f}}} $&".format(ST_median_fit[j],ST_error_high[j],ST_error_high[j])
            row_str +="${:2.1f} $ &".format(BIC_ST)
            row_str +="${:2.1f} $ &".format(BIC_PWT)
            row_str +="{}".format(Model[best_index])
            row_str +="\\\\"

            print(row_str,file=file)
        
        print("\\bottomrule",file=file)
        print("\\end{tabular}",file=file)
        print("\\begin{tablenotes}",file=file)
        print("\\small",file=file)
        print("\\item \\textbf{Notes}: Models compared: Piecewise power-law (PW), Single power-law (S),Single power-law with exponential truncation (ST), and Piecewise power-law with exponential truncation (PWT).",
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

def AICc_Table():
    print("Making AICc summary table.")
    print("####################################")
    results_directory = '../Results/'
    filename = results_directory+'AICc.tex'

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
        print("\\caption{MCMC best fit parameters and associated Akaike Information Criterion corrected for small sample size (AICc) values.}",file=file) 
        print("\\begin{tabular}{l c c c c c}",file=file)
        print("\\toprule",file=file)    
        print("\\multicolumn{1}{l}{Galaxy}& \\multicolumn{1}{c}{$\\mathrm{AICc}_{\\mathrm{PW}}$} &\
\\multicolumn{1}{c}{$\\mathrm{AICc}_{\\mathrm{S}}$} & \\multicolumn{1}{c}{$\\mathrm{AICc}_{\\mathrm{ST}}$}\
 &\\multicolumn{1}{c}{$\\mathrm{AICc}_{\\mathrm{PWT}}$} &\\multicolumn{1}{c}{Model} \\\\ % Column names row",file=file)
        #print(" &$\\alpha_1$&$\\alpha_2$&$\\beta$ & & $\\alpha{\\mathrm{S}}$ & & $\\alpha_{\\mathrm{ST}}$& $\\theta_c$ & & \\\\",file=file)
        print("\\midrule",file=file)

        i = 0

        for galaxy_name in list_of_galaxies:
            galaxy_dir = "../Results/Galaxies/{}/Masked".format(galaxy_name)
            galaxy_class = loadObj("{}/{}_summary".format(galaxy_dir,galaxy_name))

            #Get number of samples
            corr_fit = galaxy_class.corr[np.where(galaxy_class.corr>0.0)].astype(np.float)
            nsamples = np.size(corr_fit)

            #Read in PW fits
            sampler = loadObj(galaxy_dir+'/PiecePL_MCMC/'+'MCMC_sampler')
            samples = sampler.flatchain

            #Read in Single PL fit
            sampler = loadObj(galaxy_dir+'/SinglePL_MCMC/'+'MCMC_sampler')
            samples = sampler.flatchain

            
            #Read in SinglePL + Truncation fit
            sampler = loadObj(galaxy_dir+'/SingleTrunc_MCMC/'+'MCMC_sampler')
            samples = sampler.flatchain

            #Get AICc values and find best model
            AICc_S, AICc_PW, AICc_ST, AICc_PWT = compare_AICc(galaxy_name,nsamples)
            AICc_Values = [AICc_PW, AICc_S, AICc_ST,AICc_PWT]
            Model = ['PW', 'S', 'ST','PWT']
            best_index = np.argmin(AICc_Values)


            #Write to table
            row_str = ""
            name = galaxy_class.name.split('_')[0] + ' ' + galaxy_class.name.split('_')[1]
            row_str +="{} &".format(name)

            # for j in range(1,4):
            #     row_str +="${:2.1f}^{{+{:3.2f}}}_{{-{:3.2f}}} $&".format(PW_median_fit[j],PW_error_high[j],PW_error_low[j])

            row_str +="${:2.1f} $&".format(AICc_PW)

            #row_str +="${:2.1f}^{{+{:3.2f}}}_{{-{:3.2f}}} $&".format(S_median_fit[1],S_error_high[1],S_error_high[1])

            row_str +="${:2.1f} $ &".format(AICc_S)
            # for j in range(1,3):
            #     row_str +="${:2.1f}^{{+{:3.2f}}}_{{-{:3.2f}}} $&".format(ST_median_fit[j],ST_error_high[j],ST_error_high[j])
            row_str +="${:2.1f} $ &".format(AICc_ST)
            row_str +="${:2.1f} $ &".format(AICc_PWT)
            row_str +="{}".format(Model[best_index])
            row_str +="\\\\"

            print(row_str,file=file)
        
        print("\\bottomrule",file=file)
        print("\\end{tabular}",file=file)
        print("\\begin{tablenotes}",file=file)
        print("\\small",file=file)
        print("\\item \\textbf{Notes}: Models compared: Piecewise power-law (PW), Single power-law (S),Single power-law with exponential truncation (ST), and Piecewise power-law with exponential truncation (PWT).",
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

def compare_AIC(galaxy_name,omega1=False):
    #Read in samplers
    if(omega1):
        sampler_single = loadObj('../Results/Galaxies/{}/Masked/Omega1/\
SinglePL_MCMC/MCMC_sampler'.format(galaxy_name))

        sampler_piecewise = loadObj('../Results/Galaxies/{}/Masked/Omega1/\
PiecePL_MCMC/MCMC_sampler'.format(galaxy_name))

        sampler_singletrunc = loadObj('../Results/Galaxies/{}/Masked/Omega1/\
SingleTrunc_MCMC/MCMC_sampler'.format(galaxy_name))

        sampler_piecewisetrunc = loadObj('../Results/Galaxies/{}/Masked/Omega1/\
PiecewiseTrunc_MCMC/MCMC_sampler'.format(galaxy_name))
    else:
        sampler_single = loadObj('../Results/Galaxies/{}/Masked/\
SinglePL_MCMC/MCMC_sampler'.format(galaxy_name))

        sampler_piecewise = loadObj('../Results/Galaxies/{}/Masked/\
PiecePL_MCMC/MCMC_sampler'.format(galaxy_name))

        sampler_singletrunc = loadObj('../Results/Galaxies/{}/Masked/\
SingleTrunc_MCMC/MCMC_sampler'.format(galaxy_name))

        sampler_piecewisetrunc = loadObj('../Results/Galaxies/{}/Masked/\
PiecewiseTrunc_MCMC/MCMC_sampler'.format(galaxy_name))
    
    AIC_single = 2*2.0 - 2*np.max(sampler_single.flatlnprobability)
    AIC_piecewise = 2*4.0 - 2*np.max(sampler_piecewise.flatlnprobability)
    AIC_singletrunc = 2*3.0 - 2*np.max(sampler_singletrunc.flatlnprobability)
    AIC_piecewisetrunc = 2*5.0 - 2*np.max(sampler_piecewisetrunc.flatlnprobability)
    
    return AIC_single,AIC_piecewise, AIC_singletrunc, AIC_piecewisetrunc

def compare_BIC(galaxy_name,nsamples,omega1=False):
    if(omega1):
        sampler_single = loadObj('../Results/Galaxies/{}/Masked/Omega1/\
SinglePL_MCMC/MCMC_sampler'.format(galaxy_name))

        sampler_piecewise = loadObj('../Results/Galaxies/{}/Masked/Omega1/\
PiecePL_MCMC/MCMC_sampler'.format(galaxy_name))

        sampler_singletrunc = loadObj('../Results/Galaxies/{}/Masked/Omega1/\
SingleTrunc_MCMC/MCMC_sampler'.format(galaxy_name))

        sampler_piecewisetrunc = loadObj('../Results/Galaxies/{}/Masked/Omega1/\
PiecewiseTrunc_MCMC/MCMC_sampler'.format(galaxy_name))
    else:
        sampler_single = loadObj('../Results/Galaxies/{}/Masked/\
SinglePL_MCMC/MCMC_sampler'.format(galaxy_name))

        sampler_piecewise = loadObj('../Results/Galaxies/{}/Masked/\
PiecePL_MCMC/MCMC_sampler'.format(galaxy_name))

        sampler_singletrunc = loadObj('../Results/Galaxies/{}/Masked/\
SingleTrunc_MCMC/MCMC_sampler'.format(galaxy_name))

        sampler_piecewisetrunc = loadObj('../Results/Galaxies/{}/Masked/\
PiecewiseTrunc_MCMC/MCMC_sampler'.format(galaxy_name))
    
    BIC_single = np.log(nsamples)*2.0 - 2*np.max(sampler_single.flatlnprobability)
    BIC_piecewise = np.log(nsamples)*4.0 - 2*np.max(sampler_piecewise.flatlnprobability)
    BIC_singletrunc = np.log(nsamples)*3.0 - 2*np.max(sampler_singletrunc.flatlnprobability)
    BIC_piecewisetrunc = np.log(nsamples)*5.0 - 2*np.max(sampler_piecewisetrunc.flatlnprobability)
    
    return BIC_single,BIC_piecewise, BIC_singletrunc, BIC_piecewisetrunc

def compare_AICc(galaxy_name,nsamples,omega1=False):
    if(omega1):
        sampler_single = loadObj('../Results/Galaxies/{}/Masked/Omega1/\
SinglePL_MCMC/MCMC_sampler'.format(galaxy_name))

        sampler_piecewise = loadObj('../Results/Galaxies/{}/Masked/Omega1/\
PiecePL_MCMC/MCMC_sampler'.format(galaxy_name))

        sampler_singletrunc = loadObj('../Results/Galaxies/{}/Masked/Omega1/\
SingleTrunc_MCMC/MCMC_sampler'.format(galaxy_name))

        sampler_piecewisetrunc = loadObj('../Results/Galaxies/{}/Masked/Omega1/\
PiecewiseTrunc_MCMC/MCMC_sampler'.format(galaxy_name))
    else:
        sampler_single = loadObj('../Results/Galaxies/{}/Masked/\
SinglePL_MCMC/MCMC_sampler'.format(galaxy_name))

        sampler_piecewise = loadObj('../Results/Galaxies/{}/Masked/\
PiecePL_MCMC/MCMC_sampler'.format(galaxy_name))

        sampler_singletrunc = loadObj('../Results/Galaxies/{}/Masked/\
SingleTrunc_MCMC/MCMC_sampler'.format(galaxy_name))

        sampler_piecewisetrunc = loadObj('../Results/Galaxies/{}/Masked/\
PiecewiseTrunc_MCMC/MCMC_sampler'.format(galaxy_name))

    AIC_single = 2*2.0 - 2*np.max(sampler_single.flatlnprobability)
    AIC_piecewise = 2*4.0 - 2*np.max(sampler_piecewise.flatlnprobability)
    AIC_singletrunc = 2*3.0 - 2*np.max(sampler_singletrunc.flatlnprobability)
    AIC_piecewisetrunc = 2*5.0 - 2*np.max(sampler_piecewisetrunc.flatlnprobability)
    
    AICc_single = AIC_single + (2*2*(2+1))/(nsamples - 2 - 1)
    AICc_piecewise = AIC_piecewise + (2*4*(4+1))/(nsamples - 4 - 1)
    AICc_singletrunc = AIC_singletrunc + (2*3*(3+1))/(nsamples - 3 - 1)
    AICc_piecewisetrunc = AIC_piecewisetrunc + (2*5*(5+1))/(nsamples - 5 - 1)
    
    return AICc_single,AICc_piecewise, AICc_singletrunc, AICc_piecewisetrunc

def Prob_AIC(AIC_model,AIC_min):
    P = np.exp((AIC_min - AIC_model)/2.)
    return P

def Get_Cutoff_Scale(galaxy_class,function='best'):
    """
    Returns cutoff scale for a fitted function for a galaxy in arcsec. 
    
    Parameters
    ----------
    results_directory : class instance
        Galaxy class instance
    function           : String
        Function for which to return cutoff scale
    
    Returns
    -------
    None

    """

    indices = np.where(galaxy_class.corr>0.0)
    corr_fit = galaxy_class.corr[indices].astype(np.float)
    nsamples = np.size(corr_fit)
    galaxy_name = galaxy_class.name
    
    if(function == 'best'):
        AIC_single,AIC_piecewise, AIC_single_trunc, AIC_double_trunc = compare_AICc(galaxy_name,nsamples)
        galaxy_functions = ['singlepl','piecewise','singletrunc','doubletrunc']
        galaxy_AIC = [AIC_single,AIC_piecewise,AIC_single_trunc,AIC_double_trunc] 
        galaxy_function = galaxy_functions[np.argmin(galaxy_AIC)]
    else:
        galaxy_function = function

    #Read in best fit depending on best fit
    if(galaxy_function == 'piecewise'):            
        sampler = loadObj('../Results/Galaxies/{}/Masked/\
PiecePL_MCMC/MCMC_sampler'.format(galaxy_name))
    elif(galaxy_function == 'singlepl') :
        sampler = loadObj('../Results/Galaxies/{}/Masked/\
SinglePL_MCMC/MCMC_sampler'.format(galaxy_name))
    elif(galaxy_function == 'singletrunc') :
        sampler = loadObj('../Results/Galaxies/{}/Masked/\
SingleTrunc_MCMC/MCMC_sampler'.format(galaxy_name))
    elif(galaxy_function == 'doubletrunc') :
        sampler = loadObj('../Results/Galaxies/{}/Masked/\
PiecewiseTrunc_MCMC/MCMC_sampler'.format(galaxy_name))
        
    samples = sampler.flatchain
    galaxy_class.fit_values = samples[np.argmax(sampler.flatlnprobability)]
    galaxy_class.fit_errors = samples.std(axis=0)

    #SinglePL: The last non-zero bin
    if(galaxy_function == 'singlepl'):
        indices = np.where(galaxy_class.corr>0.0)
        separation_bins = galaxy_class.bin_centres*(1./arcsec_to_degree)
        separation_bins = separation_bins.astype(np.float)
        separation_bins = separation_bins[indices].astype(np.float)
        cutoff_scale = np.max(separation_bins)
        cutoff_error = np.diff(separation_bins)[-1]/2.
    # Piecewise PL : The break point
    elif(galaxy_function == 'piecewise'):
        cutoff_scale = np.exp(galaxy_class.fit_values[3])
        cutoff_error = np.exp(galaxy_class.fit_errors[3]) 
    # Single Trunc: The truncation scale
    elif(galaxy_function == 'singletrunc'):
        cutoff_scale = galaxy_class.fit_values[2]
        cutoff_error = galaxy_class.fit_errors[2]
    #Double Power Law with Trunc : The truncation scale
    elif(galaxy_function == 'doubletrunc'):
        cutoff_scale = galaxy_class.fit_values[4]
        cutoff_error = galaxy_class.fit_errors[4]

    return cutoff_scale,cutoff_error





if __name__ == "__main__":
    AIC_Table()
    BIC_Table()
    AICc_Table()
