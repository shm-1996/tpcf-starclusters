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
        print("\\begin{tabular}{l c c c c c c  }",file=file)
        print("\\toprule",file=file)
        #print("Pillar (1)& 2 & 4 & 8 & 20 & 22 & 44\\\\ % Column names row",file=file)
        print("\\multicolumn{1}{l}{(Name)}& \\multicolumn{1}{c}{Morph.}& \\multicolumn{1}{c}{$T$} &\
 \\multicolumn{1}{c}{Inclin.} & \\multicolumn{1}{c}{Dist.} & \\multicolumn{1}{c}{$\\mathrm{SFR}_{UV}$} &\
 \\multicolumn{1}{c}{$M_{*}$} \\\\ % Column names row",file=file)
        print("\\midrule",file=file)
        i = 0

        for galaxy_name in list_of_galaxies:
            galaxy_class = loadObj("../Results/Galaxies/{}/Masked/{}_summary".format(galaxy_name,galaxy_name))
            #Galaxy info table file
            info_directory = os.path.abspath('../Data/Galaxy_Information')
            Legus_Table = info_directory+'/Calzetti_Table.txt'

            #Convert self galaxy name to table name format
            galaxy_formatted = galaxy_class.name.split("_")[1]
            galaxy_formatted = '%04d'%int(galaxy_formatted)
            galaxy_formatted = galaxy_class.name.split("_")[0] + ' ' + galaxy_formatted

            galaxy_names = np.loadtxt(Legus_Table,usecols=0,delimiter='\t',dtype=str)
            index = np.where(galaxy_formatted == galaxy_names)

            morph_type = np.loadtxt(Legus_Table,usecols=2,delimiter='\t',dtype=str)
            morph_type = morph_type[index][0]

            T_type = np.loadtxt(Legus_Table,usecols=3,delimiter='\t',dtype=str)
            T_type = T_type[index][0]

            if(galaxy_name in ['NGC_1313',]):
                sfr = np.loadtxt(Legus_Table,usecols=9,delimiter='\t',dtype=str)
                mstar = np.loadtxt(Legus_Table,usecols=10,delimiter='\t',dtype=str)
            elif(galaxy_name in ['NGC_3738','NGC_4449','NGC_4656']):
                sfr = np.loadtxt(Legus_Table,usecols=10,delimiter='\t',dtype=str)
                mstar = np.loadtxt(Legus_Table,usecols=11,delimiter='\t',dtype=str)
            else :
                sfr = np.loadtxt(Legus_Table,usecols=11,delimiter='\t',dtype=str)
                mstar = np.loadtxt(Legus_Table,usecols=12,delimiter='\t',dtype=str)
            
            sfr = float(sfr[index][0])
            mstar = float(mstar[index][0])

            #Convert to power 10
            mstar = "{:2.1e}".format(mstar)
            mstarpower = int(mstar.split('e')[1]) 
            mstar = float(mstar.split('e')[0])

            #TODO: Possibly add R_25 to list of properties


            #Write string            
            row_str = ""
            name = galaxy_class.name.split('_')[0] + ' ' + galaxy_class.name.split('_')[1]
            row_str +="{} &".format(name)

            row_str +="{} &".format(morph_type)
            row_str +="${} $&".format(T_type)
            row_str +="{:2.1f} &".format(galaxy_class.inclination)
            row_str +="${:3.2f} $&".format(galaxy_class.distance)
            row_str +="${:3.2f} $&".format(sfr)
            row_str +="${:2.1f} \\times 10^{{{}}}$".format(mstar,mstarpower)
            row_str +="\\\\"
            print(row_str,file=file)

        print("\\bottomrule",file=file)
        print("\\end{tabular}",file=file)
        print("\\begin{tablenotes}",file=file)
        print("\\small",file=file)
        print("\\item \\textbf{Notes}: Info of above here.",file=file)

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


