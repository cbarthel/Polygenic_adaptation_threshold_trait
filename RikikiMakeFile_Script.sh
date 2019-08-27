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


#====================== Name data files ======================#


nameINLambdaStat00=$PATHF'/Resultats/Pairwise_LambdaStat00_table_'$locii'='$indiv'.csv'
nameINLambdaStat11=$PATHF'/Resultats/Pairwise_LambdaStat11_table_'$locii'='$indiv'.csv'
nameINLambdaStat01=$PATHF'/Resultats/Pairwise_LambdaStat01_table_'$locii'='$indiv'.csv'
nameINLambdaStat10=$PATHF'/Resultats/Pairwise_LambdaStat10_table_'$locii'='$indiv'.csv'


nameINFreqHap00=$PATHF'/Resultats/Pairwise_FHPairwise00_table_'$locii'='$indiv'.csv'
nameINFreqHap11=$PATHF'/Resultats/Pairwise_FHPairwise11_table_'$locii'='$indiv'.csv'
nameINFreqHap01=$PATHF'/Resultats/Pairwise_FHPairwise01_table_'$locii'='$indiv'.csv'
nameINFreqHap10=$PATHF'/Resultats/Pairwise_FHPairwise10_table_'$locii'='$indiv'.csv'


nameINProduitFreqAlle00=$PATHF'/Resultats/Pairwise_ProduitFreqAlle00_table_'$locii'='$indiv'.csv'


#====================== Process ======================#


#----------- Extraction of the data from Lambda00 -----------#


awk -F$',' -v OFS=',' '{print;}' $nameINLambdaStat00 > $PATHF/Resultats/Rikiki_First.tmpR


#----------- Extraction of the data from Lambda11 -----------#


awk -F$',' -v OFS=',' '{print $7}' $nameINLambdaStat11 > $PATHF/Resultats/StorageAWKdata1.tmpR
paste -d',' $PATHF/Resultats/Rikiki_First.tmpR $PATHF/Resultats/StorageAWKdata1.tmpR > $PATHF/Resultats/Rikiki_Second.tmpR


#----------- Extraction of the data from Lambda01 -----------#


awk -F$',' -v OFS=',' '{print $7}' $nameINLambdaStat01 > $PATHF/Resultats/StorageAWKdata2.tmpR
paste -d',' $PATHF/Resultats/Rikiki_Second.tmpR $PATHF/Resultats/StorageAWKdata2.tmpR > $PATHF/Resultats/Rikiki_3.tmpR


#----------- Extraction of the data from Lambda10 -----------#


awk -F$',' -v OFS=',' '{print $7}' $nameINLambdaStat10 > $PATHF/Resultats/StorageAWKdata3.tmpR
paste -d',' $PATHF/Resultats/Rikiki_3.tmpR $PATHF/Resultats/StorageAWKdata3.tmpR > $PATHF/Resultats/Rikiki_4.tmpR


#----------- Extraction of the data from FreqHaplo00 -----------#


awk -F$',' -v OFS=',' '{print $7}' $nameINFreqHap00 > $PATHF/Resultats/StorageAWKdata4.tmpR
paste -d',' $PATHF/Resultats/Rikiki_4.tmpR $PATHF/Resultats/StorageAWKdata4.tmpR > $PATHF/Resultats/Rikiki_5.tmpR


#----------- Extraction of the data from FreqHaplo11 -----------#


awk -F$',' -v OFS=',' '{print $7}' $nameINFreqHap11 > $PATHF/Resultats/StorageAWKdata5.tmpR
paste -d',' $PATHF/Resultats/Rikiki_5.tmpR $PATHF/Resultats/StorageAWKdata5.tmpR > $PATHF/Resultats/Rikiki_6.tmpR


#----------- Extraction of the data from FreqHaplo01 -----------#


awk -F$',' -v OFS=',' '{print $7}' $nameINFreqHap01 > $PATHF/Resultats/StorageAWKdata6.tmpR
paste -d',' $PATHF/Resultats/Rikiki_6.tmpR $PATHF/Resultats/StorageAWKdata6.tmpR > $PATHF/Resultats/Rikiki_7.tmpR


#----------- Extraction of the data from FreqHaplo10 -----------#


awk -F$',' -v OFS=',' '{print $7}' $nameINFreqHap10 > $PATHF/Resultats/StorageAWKdata7.tmpR
paste -d',' $PATHF/Resultats/Rikiki_7.tmpR $PATHF/Resultats/StorageAWKdata7.tmpR > $PATHF/Resultats/Rikiki_8.tmpR


#----------- Extraction of the data from ProduitFreqAlle00 -----------#


awk -F$',' -v OFS=',' '{print $7}' $nameINProduitFreqAlle00 > $PATHF/Resultats/StorageAWKdata8.tmpR
paste -d',' $PATHF/Resultats/Rikiki_8.tmpR $PATHF/Resultats/StorageAWKdata8.tmpR > $PATHF/Resultats/Rikiki_General.txt


#----------- Change the header and remove all the .tmpR files -----------#


sed -i '1d' $PATHF/Resultats/Rikiki_General.txt


sed -i '1s/^/gen,rep,locii,locij,ID,WEI,LDintra00,LDintra11,LDintra01,LDintra10,FreqHaplo00,FreqHaplo11,FreqHaplo01,FreqHaplo10,ProduitFreqAlle00\n/' $PATHF/Resultats/Rikiki_General.txt

rm -f $PATHF/Resultats/*.tmpR


#####################################################################################################################





echo "++++++++++++++++++++------------------------------ | RikikiMakeFile_Script.sh | --- | FAIT | ------------------+++++++++++++++++++++"




#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================