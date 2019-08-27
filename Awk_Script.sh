#!/bin/bash
# -*- coding: utf-8 -*-


# Script for manuscript "Detecting Polygenic Adaptation For A Threshold Trait In Experimental Populations Using Time-Series Data"


#Clément BARTHÉLÉMY, La Rochelle University, Station Ifremer de La Tremblade


#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================


#====================== Global variable ======================#


((repeat=$4)) #Argument passing by the GenePOP_Script

((nb_gen=$3)) #Argument passing by the GenePOP_Script

((nb_loci=$2)) #Argument passing by the GenePOP_Script

((nb_indiv=$1)) #Argument passing by the GenePOP_Script


#====================== Create the GenePOP format ======================#


awk -F$',' -v OFS=',' '{$1=$2=$3=$4=$5=$6="";print}' ./myPOP_GenePOP.tmp | sed '1,${/_1/d;}' 1> ./myPOP_GenePOP2.tmp #Print the column that define the genotype

awk '{gsub(",","",$1);print}' ./myPOP_GenePOP2.tmp 1> ./myPOP_GenePOP3.tmp #Delete all the "," in the first column

awk '{for(i=1;i<=NF;i+=1) gsub(",","",$i)}1' ./myPOP_GenePOP3.tmp 1> ./LD_Genotype_R_$2_$1_$4.txt #Create the text file with the genotype of all the population

awk '{for(i=1;i<=NF;i+=1) {gsub(",","",$i);gsub("1","02",$i);gsub("0","01",$i);gsub("012","02",$i)}}1' ./myPOP_GenePOP3.tmp 1> ./myPOP_GenePOP4.tmp #Operation of substitution for the GenePOP format

awk '{for(i=1;i<NF+1;i+=2) print($i $(i+1))}' ./myPOP_GenePOP4.tmp 1> ./myPOP_GenePOP5.tmp #Stick the column together (2*2)

awk -v nb_loci=$nb_loci 'BEGIN {count=0;}{if(count==nb_loci-1){count=0;printf("%s\n",$0);}else{printf("%s ",$0);count++;}}' ./myPOP_GenePOP5.tmp 1> ./myPOP_GenePOP6.tmp #Reunite the genotype for each individual

sed "0~$nb_indiv G" ./myPOP_GenePOP6.tmp | sed '$d' 1> ./myPOP_GenePOP7.tmp #Skip a line every "nb_indiv" in order to separate the generations and delete the last "skip" line


#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================