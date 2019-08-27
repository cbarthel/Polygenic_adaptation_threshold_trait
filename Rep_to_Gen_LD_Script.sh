#!/bin/bash
# -*- coding: utf-8 -*-


# Script for manuscript "Detecting Polygenic Adaptation For A Threshold Trait In Experimental Populations Using Time-Series Data"


#Clément BARTHÉLÉMY, La Rochelle University, Station Ifremer de La Tremblade


#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================


#Explication : In this script we are spliting a file with all the LD for the generations of one repetition in several file for each generation with all the repetition.


#====================== Path ======================#


PATHF="/home/clement/Documents/SimulationPP" #Insert your path here


#====================== Variables ======================#


((repeat=$4)) #Argument passing by Main_Bash_Launch.sh

((nb_gen=$3)) #Argument passing by Main_Bash_Launch.sh

((nb_loci=$2)) #Argument passing by Main_Bash_Launch.sh

((nb_indiv=$1)) #Argument passing by Main_Bash_Launch.sh


#======================  Functions  ======================#


endgen=$(( (nb_loci)*nb_gen ))


for ((rep=1;rep<=$repeat;rep+=1)) #We start at the fisrt generation of the rep
do

    cpt=0

    for ((gen=1;gen<=$endgen;gen+=$nb_loci)) #We start at the fisrt generation of the rep
    do


        cpt=$(( cpt+1 ))

        ((finish=$gen+$nb_loci-1))

        sed -n "$gen,$finish p" $PATHF/Resultats/Matrice/Matrix_LD_File_$nb_loci=$nb_indiv=$rep.tmp2 >> $PATHF/Resultats/Matrice/Matrix_LD_File_$nb_loci=$nb_indiv=$cpt.txt

        sed -n "$gen,$finish p" $PATHF/Resultats/Matrice/Matrix_LambdaStat00_File_$nb_loci=$nb_indiv=$rep.tmp2 >> $PATHF/Resultats/Matrice/Matrix_LambdaStat00_File_$nb_loci=$nb_indiv=$cpt.txt

        sed -n "$gen,$finish p" $PATHF/Resultats/Matrice/Matrix_LambdaStat11_File_$nb_loci=$nb_indiv=$rep.tmp2 >> $PATHF/Resultats/Matrice/Matrix_LambdaStat11_File_$nb_loci=$nb_indiv=$cpt.txt

        sed -n "$gen,$finish p" $PATHF/Resultats/Matrice/Matrix_LambdaStat01_File_$nb_loci=$nb_indiv=$rep.tmp2 >> $PATHF/Resultats/Matrice/Matrix_LambdaStat01_File_$nb_loci=$nb_indiv=$cpt.txt

        sed -n "$gen,$finish p" $PATHF/Resultats/Matrice/Matrix_LambdaStat10_File_$nb_loci=$nb_indiv=$rep.tmp2 >> $PATHF/Resultats/Matrice/Matrix_LambdaStat10_File_$nb_loci=$nb_indiv=$cpt.txt

        sed -n "$gen,$finish p" $PATHF/Resultats/Matrice/Matrix_FHPairwise00_File_$nb_loci=$nb_indiv=$rep.tmp2 >> $PATHF/Resultats/Matrice/Matrix_FHPairwise00_File_$nb_loci=$nb_indiv=$cpt.txt

        sed -n "$gen,$finish p" $PATHF/Resultats/Matrice/Matrix_FHPairwise11_File_$nb_loci=$nb_indiv=$rep.tmp2 >> $PATHF/Resultats/Matrice/Matrix_FHPairwise11_File_$nb_loci=$nb_indiv=$cpt.txt

        sed -n "$gen,$finish p" $PATHF/Resultats/Matrice/Matrix_FHPairwise01_File_$nb_loci=$nb_indiv=$rep.tmp2 >> $PATHF/Resultats/Matrice/Matrix_FHPairwise01_File_$nb_loci=$nb_indiv=$cpt.txt

        sed -n "$gen,$finish p" $PATHF/Resultats/Matrice/Matrix_FHPairwise10_File_$nb_loci=$nb_indiv=$rep.tmp2 >> $PATHF/Resultats/Matrice/Matrix_FHPairwise10_File_$nb_loci=$nb_indiv=$cpt.txt

        sed -n "$gen,$finish p" $PATHF/Resultats/Matrice/Matrix_ProduitFreqAlle00_File_$nb_loci=$nb_indiv=$rep.tmp2 >> $PATHF/Resultats/Matrice/Matrix_ProduitFreqAlle00_File_$nb_loci=$nb_indiv=$cpt.txt

    done
done


#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================