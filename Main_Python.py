#!/usr/bin/env python 
# -*- coding: utf8 -*-


# Script for manuscript "Detecting Polygenic Adaptation For A Threshold Trait In Experimental Populations Using Time-Series Data"


#Clément BARTHÉLÉMY, La Rochelle University, Station Ifremer de La Tremblade


#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================


#====================== Path ======================#


PATHF = os.path.join('/home/clement/Documents/SimulationPP') #BEWARE : insert your path line 65


#====================== Packages ======================#


import math,os,time,random,sys,subprocess,shutil
import simuPOP as sim
import simuOpt
from simuOpt import setOptions
setOptions(optimized=False, alleleType='short', numThreads=2, quiet=True)
import simuPOP.utils as utils
from simuPOP.utils import saveCSV,export


#====================== Variables ======================#


nb_indiv = int(sys.argv[1]) #Number of individuals in the total population

nb_loci = int(sys.argv[2]) #Number of locii

nb_gen = int(sys.argv[3]) #Number of studied generations

selection_intensity_z1 = int(sys.argv[4]) #Intensity of the selection

repeat = int(sys.argv[5])

zopt_dist = int(sys.argv[6])

micro_env = int(sys.argv[7]) #effet micro env

nb_allele = 2 #Number of allel on each locii, always 2 in our model


#====================== Set date and time for the output file ======================#


global NOM
NOM= time.strftime('%d-%m-%Y-%Hh%Mmin%Ssec',time.localtime())
print NOM


#====================== Creation of the "Results" folder ======================#


chemin = os.path.join('/home/clement/Documents/SimulationPP/Resultats')                #Insert your path here
try :
    os.chdir(chemin)
except :
    os.mkdir(chemin)        # Create "Resultats" folder
    os.chdir(chemin)        # Get in
global CHEMIN
CHEMIN = os.getcwd()


#====================== Initialization (SimuPOP) ======================#

file_to_open_1=os.path.join(PATHF, "Script", "Pop_Evolve_SimuPOP.py")

execfile(file_to_open_1)

#====================== Create the gen column ======================#


if (repeat==1):

    for g in range(nb_gen):

        Lgen=[] #List of the repeat's number

        file_to_open_2=os.path.join(PATHF, "Resultats", "GenCol%d.tmpp"%g)

        ListGEN = open(file_to_open_2, "w")

        line = ListGEN.write("gen"+"\n")

        for l in range(nb_indiv):

            g=str(g)

            line = ListGEN.write(g+"\n")
       
        ListGEN.close()

#====================== Save Informations ======================#


file_to_open_3=os.path.join(PATHF, "savepop0.pop")
file_to_open_4=os.path.join(PATHF, "Resultats", "savepop0.pop")

shutil.copy(file_to_open_3, file_to_open_4)

for g in range(nb_gen):

    file_to_open_5=os.path.join(PATHF, "Resultats", "savepop%d.pop"%g)

    pop = sim.loadPopulation(file_to_open_5)

    file_to_open_6=os.path.join(PATHF, "Resultats", "myPOP_Generation_"+str(g)+".tmp")

    export(pop,format='csv',infoFields=['sum_g1','z1','fitness'],output=file_to_open_6, gui=False) 


#====================== Execute the script wich paste the gen column ======================#


file_to_open_7=os.path.join(PATHF, "Script", "GenCol_Script.sh %s %s")

subprocess.call( [file_to_open_7 % (nb_indiv, nb_gen)], shell=True)


#====================== Time ======================#


start_time = time.time()
interval = time.time()-start_time
print'== > Total time in seconds = ', interval


#====================== Generate the GenePOP file and execute it ======================#


file_to_open_8=os.path.join(PATHF, "Script", "GenePOP_Script.py %s %s %s %s")

subprocess.call( [file_to_open_8 % (nb_indiv, nb_loci, nb_gen, repeat)], shell=True)


#####################################################################################################################





print "+++++++++++++++----------------------------------- | Main_Python.py           | --- | FAIT | ------------------++++++++++++++++++++"





#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================