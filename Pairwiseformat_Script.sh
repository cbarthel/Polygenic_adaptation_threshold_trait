#!/bin/bash
# -*- coding: utf-8 -*-


# Script for manuscript "Detecting Polygenic Adaptation For A Threshold Trait In Experimental Populations Using Time-Series Data"


#Clément BARTHÉLÉMY, La Rochelle University, Station Ifremer de La Tremblade


#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================


#====================== Path ======================#


PATHF="/home/clement/Documents/SimulationPP" #Insert your path here


#====================== Variables ======================#


((nb_gen=$3)) #Argument passing by the GenePOP_Script

((indiv=$1)) #Argument passing by the GenePOP_Script

((locii=$2)) #Argument passing by the GenePOP_Script

nameIN=$4'_'$locii'='$indiv'.tmp'

nameOUT=$4'_'$locii'='$indiv'.csv'


#====================== Process ======================#


sed -e 's/ /,/g' $PATHF/Resultats/$nameIN > $PATHF/Resultats/$nameOUT

sed -i '1s/^/gen,rep,locii,locij,ID,WEI,LD\n/' $PATHF/Resultats/$nameOUT

#####################################################################################################################





echo "++++++++++++++++++++------------------------------ | Pairwiseformat_Script.sh | --- | FAIT | ------------------+++++++++++++++++++++"





#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================