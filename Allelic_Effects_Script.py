#!/usr/bin/env python 
# -*- coding: utf-8 -*-


# Script for manuscript "Detecting Polygenic Adaptation For A Threshold Trait In Experimental Populations Using Time-Series Data"


#Clément BARTHÉLÉMY, La Rochelle University, Station Ifremer de La Tremblade


#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================


#====================== Packages ======================#


import math,random,sys,os
import numpy as np


#====================== Path ======================#


PATHF="/home/clement/Documents/SimulationPP" #Insert your path here


#====================== Variables ======================#


nb_loci = int(sys.argv[1]) #Number of locii

nb_loci_neutral=nb_loci-10

var_allelic_effect = int(sys.argv[2]) #Number of locii

nb_allele = 2 #Number of alleles

weights = [5 for x in xrange(100)] # Weights of each locus (we do it for 100 locii)

allelic_effects = [[0 for x in range(nb_allele)] for x in range(nb_loci)] #Matrix for storage all the allelic effects, initiate at 0

allele_freq = [0 for x in range(nb_allele)]


#====================== Functions ======================#


#-----------  We fill up the matrix of loci weights -----------#


Vector_wei=np.append([1 for x in xrange(10)], [0 for x in xrange(nb_loci_neutral)])




cpt_place = 0

for charac in Vector_wei:
    if (charac!=","):

        number=int(charac)
        weights[cpt_place]=number
        cpt_place=cpt_place+1

        if (cpt_place==100): #BEAWARE : there is a ghost character at the end of the file "Weights_Locii_File.txt"
            break


#-----------  We fill up the matrix of allelic effetcs -----------#


def init_allelic_effects():
    for locus in range(nb_loci):
        for allele in range(nb_allele):
            allelic_effects[locus][allele] = weights[locus] * random.normalvariate(0,var_allelic_effect) #See the associate article
    return (allelic_effects)


a=init_allelic_effects() #Call the function

a=str(a) #Converting to string

Aeffects = open(os.path.join(PATHF, "Resultats", "Allelic_Effects_File.txt"), "w") #Open the file

Aeffects.write(a) #Writing in the text file

Aeffects.close() #Close the file


#####################################################################################################################





print "++++++++++++++++++++------------------ | Allelic_Effect_Script.py    | DONE | ------------------++++++++++++++++++++"





#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================