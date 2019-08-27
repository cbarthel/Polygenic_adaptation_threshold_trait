#!/usr/bin/env python 
# -*- coding: utf8 -*-


# Script for manuscript "Detecting Polygenic Adaptation For A Threshold Trait In Experimental Populations Using Time-Series Data"


#Clément BARTHÉLÉMY, La Rochelle University, Station Ifremer de La Tremblade


#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================


#======================= Packages =======================#


import sys,os,subprocess


#======================= Converting the input argument into integer =======================#


nb_indiv = int(sys.argv[1]) #Number of individuals in the total population

nb_loci = int(sys.argv[2]) #Number of locii

nb_gen = int(sys.argv[3]) #Number of studied generations

repeat = int(sys.argv[4]) #The repetition in which we are


PATHF="/home/clement/Documents/SimulationPP" #Insert your path here


#======================= Fusion of text files (Bash Command) =======================#


os.system("cat ./*.tmp5 | sed '2,${/gen/d;}' > myPOP.tmp4") # Open all the text files, delete all the header except for the first one
os.system("head -1 myPOP.tmp4 > myPOP.txt2") #Save the header in a temporary file
os.system("tail -n+2 myPOP.tmp4 | sort -n >> myPOP.txt2") #Sort the file and add the header
os.system("rm -f ./*.tmp") #Remove the .tmp file
os.system("rm -f ./*.pop") #Remove the .pop file

os.system("sed -e 's/,/, /g' myPOP.txt2 > myPOP.txt")

command1="rm "+os.path.join(PATHF, "Resultats", "myPOP.txt2")

os.system(command1)

os.system("cp myPOP.txt myPOP_GenePOP.tmp") #Copy the file so we can use it for GenePOP

os.system("cp myPOP.txt myPOP_Rtemp.tmp") #Copy the file so we can use it for R


#================ Creation of the sample text file : Format GenePOP ================#


#-----------  Creation of a file with all the name (Population,individual...) (Python Command)  -----------#


Lloci=[] #List of the name of each loci


for i in range(nb_loci):
    Lloci.append(str(" Loc")+str(i+1)) #Add the name of each loci to the list

os.system("touch POPname.tmp") #Bash Command : create a new text file

ListLOCUS = open(os.path.join(PATHF, "Resultats", "POPname.tmp"), "a")

for k in range(nb_gen):
    if (k!=0):
        line = ListLOCUS.write(" Pop"+"\n")
    for j in range(nb_indiv):
        if (k==nb_gen-1 and j==nb_indiv-1):
            line = ListLOCUS.write(str(" Leucothée")+str(j)+'_'+str(k)+"  ,  ")
        else:
            line = ListLOCUS.write(str(" Leucothée")+str(j)+'_'+str(k)+"  ,  \n")

        
ListLOCUS.close()


#-----------  Execute the bash file with the "awk" commands  -----------#


file_to_open_1=os.path.join(PATHF, "Script", "Awk_Script.sh %s %s %s %s")

subprocess.call( [file_to_open_1 % (nb_indiv, nb_loci, nb_gen, repeat)], shell=True)


#-----------  Create a header for the GenePOP file (Python Command)  -----------#


os.system("touch Header.tmp") #Bash Command : create a new text file

ListLOCUS = open(os.path.join(PATHF, "Resultats", "Header.tmp"), "a")

ListLOCUS.write(" Title : Populations under selection pressure \n") #Write a title

for l in Lloci:
    line = ListLOCUS.write(l+"\n") #Write the name of each loci in the text file (\n change the line in the text file)

line = ListLOCUS.write(" Pop"+"\n")


ListLOCUS.close()


#-----------  Merge the file in order to create the final GenePOP file (Bash Command)  -----------#


os.system('paste -d "" POPname.tmp ./myPOP_GenePOP7.tmp 1> WithOutHeader.tmp')

subprocess.call( ["cat ./Header.tmp ./WithOutHeader.tmp > myGenePOP_loci=%s_indd=%s_repeat=%s.gen" % (nb_loci, nb_indiv, repeat)], shell=True)

os.system("rm -f ./*.tmp") #Remove the .tmp file


#####################################################################################################################





print "+++++++++++++++----------------------------------- | GenePOP_Script.py        | --- | FAIT | ------------------++++++++++++++++++++"





#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================