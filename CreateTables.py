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



def GalaxyTPCF_Table():

    return


