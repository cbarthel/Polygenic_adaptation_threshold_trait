#!/bin/bash
# -*- coding: utf-8 -*-


# Script for manuscript "Detecting Polygenic Adaptation For A Threshold Trait In Experimental Populations Using Time-Series Data"


#Clément BARTHÉLÉMY, La Rochelle University, Station Ifremer de La Tremblade


#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================


#====================== Path ======================#


PYTHON=$((1)) #Run simulations if 1

R=$((1)) #Run calcul of haplotype frequencies if 1


PATHF="/home/clement/Documents/SimulationPP" #Insert your path here


#====================== Read each scenario ======================#


readarray rows < $PATHF/Script/sc.txt

extension=$((1))

for row in "${rows[@]}";do
    row_array=(${row})

    num_of_repeat=$((row_array[0]))

    nb_gen=$((row_array[1]))

    nb_indiv=$((row_array[2])) #More than 20

    nb_loci=$((row_array[3])) #Less or equal than 100

    zopt_dist=$((50)) #Distance between z and zopt 50

    selection_intensity_z1=$((5)) # 1 or more 5

    nb_cluster=$((row_array[4])) # 1 or more 5

    micro_env=$((row_array[5]))

    var_allelic_effect=$((row_array[6]))


#====================== Create folders ======================#


    mkdir $PATHF/Resultats/

    mkdir $PATHF/Resultats/Courbes

    mkdir $PATHF/Resultats/Matrice

    mkdir $PATHF/Resultats/Donnees


#====================== Save the parameters of the simulation ======================#


    echo "Genetic Data Creation = $PYTHON" > $PATHF/Resultats/Simulation_Parameters.txt
    echo "Geneteic Data Analysis = $R" >> $PATHF/Resultats/Simulation_Parameters.txt
    echo "Number of repeat = $num_of_repeat" >> $PATHF/Resultats/Simulation_Parameters.txt
    echo "Number of generations = $nb_gen" >> $PATHF/Resultats/Simulation_Parameters.txt
    echo "Number of individuals = $nb_indiv" >> $PATHF/Resultats/Simulation_Parameters.txt
    echo "Number of locii = $nb_loci" >> $PATHF/Resultats/Simulation_Parameters.txt
    echo "Distance between z and zopt = $zopt_dist" >> $PATHF/Resultats/Simulation_Parameters.txt
    echo "Selection intensity = $selection_intensity_z1" >> $PATHF/Resultats/Simulation_Parameters.txt
    echo "Number of cluster = $nb_cluster" >> $PATHF/Resultats/Simulation_Parameters.txt
    echo "Micro env = $micro_env" >> $PATHF/Resultats/Simulation_Parameters.txt
    echo "Var allelic effect = $var_allelic_effect" >> $PATHF/Resultats/Simulation_Parameters.txt


#====================== Start simulation ======================#


    bash $PATHF/Script/Main_Bash.sh $num_of_repeat $nb_gen $nb_indiv $nb_loci $selection_intensity_z1 $zopt_dist $nb_cluster $micro_env $var_allelic_effect

    mv $PATHF/Resultats/ $PATHF/Resultats_Scenario_$extension

    extension=$((extension+1))

done


#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================