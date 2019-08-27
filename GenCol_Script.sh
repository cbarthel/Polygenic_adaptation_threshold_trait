#!/bin/bash
# -*- coding: utf-8 -*-


# Script for manuscript "Detecting Polygenic Adaptation For A Threshold Trait In Experimental Populations Using Time-Series Data"


#Clément BARTHÉLÉMY, La Rochelle University, Station Ifremer de La Tremblade


#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================


#====================== Global variable ======================#


((nb_gen=$2)) #Argument passing by the GenePOP_Script

((nb_indiv=$1)) #Argument passing by the GenePOP_Script


PATHF="/home/clement/Documents/SimulationPP" #Insert your path here


#====================== Process ======================#


for ((gen=0;gen<=$nb_gen-1;gen+=1))
do

    paste -d", " $PATHF/Resultats/GenCol$gen.tmpp $PATHF/Resultats/myPOP_Generation_$gen.tmp > $PATHF/Resultats/myPOP_Generation2_$gen.tmp5

done


#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================