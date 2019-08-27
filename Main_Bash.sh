#!/bin/bash
# -*- coding: utf-8 -*-


# Script for manuscript "Detecting Polygenic Adaptation For A Threshold Trait In Experimental Populations Using Time-Series Data"


#Clément BARTHÉLÉMY, La Rochelle University, Station Ifremer de La Tremblade


#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================


STARTTIME=$(date +%s) #In order to calculate the duration of the Simulation


#====================== Path ======================#


PATHF="/home/clement/Documents/SimulationPP" #Insert your path here


#====================== Global variable ======================#


PYTHON=$((1)) #Argument passing by Main_Launch

R=$((1)) #Argument passing by Main_Launch


((num_of_repeat=$1)) #Argument passing by Main_Launch

((Haplo=1)) #Argument passing by Main_Launch

((nb_gen=$2)) #Argument passing by Main_Launch

((nb_indiv=$3)) #Argument passing by Main_Launch

((nb_loci=$4)) #Argument passing by Main_Launch

((selection_intensity_z1=$5)) #Argument passing by Main_Launch

((zopt_dist=$6)) #Argument passing by Main_Launch

((nb_cluster=$7)) #Argument passing by Main_Launch

((micro_env=$8)) #Argument passing by Main_Launch

((var_allelic_effect=$9)) #Argument passing by Main_Launch


#====================== Filenames (local variable) ======================#


nameIN_Stat00_Statscript="Matrix_LambdaStat00_File_"
nameOUT_Stat00_Statscript="LambdaStat00_Pairwise_"

nameIN_Stat00_table="LambdaStat00_Pairwise_"
nameOUT_Stat00_table="Pairwise_LambdaStat00_table_"

nameIN_Stat00_plot="Pairwise_LambdaStat00_table_"
nameOUT_Stat00_plot="LambdaStat_00_Pairwise_Plot_General"

nameIN_Stat00_format="Pairwise_LambdaStat00_table"

#------------------------------------------------------------------------------------------------------

nameIN_Stat11_Statscript="Matrix_LambdaStat11_File_"
nameOUT_Stat11_Statscript="LambdaStat11_Pairwise_"

nameIN_Stat11_table="LambdaStat11_Pairwise_"
nameOUT_Stat11_table="Pairwise_LambdaStat11_table_"

nameIN_Stat11_plot="Pairwise_LambdaStat11_table_"
nameOUT_Stat11_plot="LambdaStat_11_Pairwise_Plot_General"

nameIN_Stat11_format="Pairwise_LambdaStat11_table"

#------------------------------------------------------------------------------------------------------

nameIN_Stat01_Statscript="Matrix_LambdaStat01_File_"
nameOUT_Stat01_Statscript="LambdaStat01_Pairwise_"

nameIN_Stat01_table="LambdaStat01_Pairwise_"
nameOUT_Stat01_table="Pairwise_LambdaStat01_table_"

nameIN_Stat01_plot="Pairwise_LambdaStat01_table_"
nameOUT_Stat01_plot="LambdaStat_01_Pairwise_Plot_General"

nameIN_Stat01_format="Pairwise_LambdaStat01_table"

#------------------------------------------------------------------------------------------------------

nameIN_Stat10_Statscript="Matrix_LambdaStat10_File_"
nameOUT_Stat10_Statscript="LambdaStat10_Pairwise_"

nameIN_Stat10_table="LambdaStat10_Pairwise_"
nameOUT_Stat10_table="Pairwise_LambdaStat10_table_"

nameIN_Stat10_plot="Pairwise_LambdaStat10_table_"
nameOUT_Stat10_plot="LambdaStat_10_Pairwise_Plot_General"

nameIN_Stat10_format="Pairwise_LambdaStat10_table"

#------------------------------------------------------------------------------------------------------

nameIN_FH00_Statscript="Matrix_FHPairwise00_File_"
nameOUT_FH00_Statscript="FHPairwise00_Pairwise_"

nameIN_FH00_table="FHPairwise00_Pairwise_"
nameOUT_FH00_table="Pairwise_FHPairwise00_table_"

nameIN_FH00_plot="Pairwise_FHPairwise00_table_"
nameOUT_FH00_plot="FrequenceHaplo_00_Pairwise_Plot_General"

nameIN_FH00_format="Pairwise_FHPairwise00_table"
#------------------------------------------------------------------------------------------------------

nameIN_FH11_Statscript="Matrix_FHPairwise11_File_"
nameOUT_FH11_Statscript="FHPairwise11_Pairwise_"

nameIN_FH11_table="FHPairwise11_Pairwise_"
nameOUT_FH11_table="Pairwise_FHPairwise11_table_"

nameIN_FH11_plot="Pairwise_FHPairwise11_table_"
nameOUT_FH11_plot="FrequenceHaplo_11_Pairwise_Plot_General"

nameIN_FH11_format="Pairwise_FHPairwise11_table"

#------------------------------------------------------------------------------------------------------

nameIN_FH01_Statscript="Matrix_FHPairwise01_File_"
nameOUT_FH01_Statscript="FHPairwise01_Pairwise_"

nameIN_FH01_table="FHPairwise01_Pairwise_"
nameOUT_FH01_table="Pairwise_FHPairwise01_table_"

nameIN_FH01_plot="Pairwise_FHPairwise01_table_"
nameOUT_FH01_plot="FrequenceHaplo_01_Pairwise_Plot_General"

nameIN_FH01_format="Pairwise_FHPairwise01_table"

#------------------------------------------------------------------------------------------------------

nameIN_FH10_Statscript="Matrix_FHPairwise10_File_"
nameOUT_FH10_Statscript="FHPairwise10_Pairwise_"

nameIN_FH10_table="FHPairwise10_Pairwise_"
nameOUT_FH10_table="Pairwise_FHPairwise10_table_"

nameIN_FH10_plot="Pairwise_FHPairwise10_table_"
nameOUT_FH10_plot="FrequenceHaplo_10_Pairwise_Plot_General"

nameIN_FH10_format="Pairwise_FHPairwise10_table"

#------------------------------------------------------------------------------------------------------

nameIN_PFA00_Statscript="Matrix_ProduitFreqAlle00_File_"
nameOUT_PFA00_Statscript="ProduitFreqAlle00_Pairwise_"

nameIN_PFA00_table="ProduitFreqAlle00_Pairwise_"
nameOUT_PFA00_table="Pairwise_ProduitFreqAlle00_table_"

nameIN_PFA00_plot="Pairwise_ProduitFreqAlle00_table_"
nameOUT_PFA00_plot="Pairwise_ProduitFreqAlle00_Plot"

nameIN_PFA00_format="Pairwise_ProduitFreqAlle00_table"


#======================  Execute the Python Main  ======================#


if [ $PYTHON -eq 1 ]
then

    for ((locii=$nb_loci;locii<=$nb_loci;locii+=$nb_loci))
    do

        python Allelic_Effects_Script.py "$locii" "$var_allelic_effect"

        for ((indd=$nb_indiv;indd<=$nb_indiv;indd+=$nb_indiv))
        do

            touch $PATHF/Resultats/myRUN.txt #Create a new text file for each combinaison

            for ((repeat=1;repeat<=$num_of_repeat;repeat+=1))
            do

                python Main_Python.py $indd $locii $nb_gen $selection_intensity_z1 $repeat $zopt_dist $micro_env

                paste -d", " $PATHF/Resultats/RepCol.csv $PATHF/Resultats/myPOP.txt >> $PATHF/Resultats/myRUN.txt

            done
            
            sed '2,${/rep/d;}' $PATHF/Resultats/myRUN.txt > $PATHF/Resultats/myRUN_$locii'_'$indd.tmp7

            sed 's/,/, /g' $PATHF/Resultats/myRUN_$locii'_'$indd.tmp7 > $PATHF/Resultats/myRUN_$locii'_'$indd.tmp8

            sed 's/,  /, /g' $PATHF/Resultats/myRUN_$locii'_'$indd.tmp8 > $PATHF/Resultats/myRUN_$locii'_'$indd.csv

            rm $PATHF/Resultats/myRUN.txt

            rm $PATHF/Resultats/myRUN_$locii'_'$indd.tmp7

            rm $PATHF/Resultats/myRUN_$locii'_'$indd.tmp8

            rm $PATHF/Resultats/RepCol.csv

        done    
    done

    rm $PATHF/Resultats/myPOP.txt

    rm -f $PATHF/Resultats/*.tmpp

    rm -f $PATHF/Resultats/*.tmp5

    rm -f $PATHF/Resultats/*.tmp4
fi


#======================  Execute the R analysis scripts  ======================#

date

if [ $R -eq 1 ]
then

    for ((locii=$nb_loci;locii<=$nb_loci;locii+=$nb_loci))
    do
        for ((indd=$nb_indiv;indd<=$nb_indiv;indd+=$nb_indiv))
        do
            for ((repeat=1;repeat<=$num_of_repeat;repeat+=1))
            do

                if [[ $repeat -eq $num_of_repeat ]]
                    then

# -------------------------- Haplotypic frequencies and LD --------------------------

                        date
                        python $PATHF/Script/Matrix_LD_Intra_Script.py $indd $locii $nb_gen $repeat
                        date

                        bash $PATHF/Script/Rep_to_Gen_LD_Script.sh $indd $locii $nb_gen $repeat


                        python Stats_Script.py "$indd" "$locii" "$nb_gen" "$repeat" "$nameIN_Stat00_Statscript" "$nameOUT_Stat00_Statscript"

                        Rscript $PATHF/Script/Pairwisetable_Script.r $indd $locii $nb_gen $repeat $nameIN_Stat00_table $nameOUT_Stat00_table

                        bash $PATHF/Script/Pairwiseformat_Script.sh $indd $locii $nb_gen $nameIN_Stat00_format

                        Rscript $PATHF/Script/Pairwiseplot_Script.r $indd $locii $nb_gen $repeat $nameIN_Stat00_plot $nameOUT_Stat00_plot




                        python Stats_Script.py "$indd" "$locii" "$nb_gen" "$repeat" "$nameIN_FH00_Statscript" "$nameOUT_FH00_Statscript"

                        Rscript $PATHF/Script/Pairwisetable_Script.r $indd $locii $nb_gen $repeat $nameIN_FH00_table $nameOUT_FH00_table

                        bash $PATHF/Script/Pairwiseformat_Script.sh $indd $locii $nb_gen $nameIN_FH00_format

                        Rscript $PATHF/Script/Pairwiseplot_Script.r $indd $locii $nb_gen $repeat $nameIN_FH00_plot $nameOUT_FH00_plot




                        python Stats_Script.py "$indd" "$locii" "$nb_gen" "$repeat" "$nameIN_Stat11_Statscript" "$nameOUT_Stat11_Statscript"

                        Rscript $PATHF/Script/Pairwisetable_Script.r $indd $locii $nb_gen $repeat $nameIN_Stat11_table $nameOUT_Stat11_table

                        bash $PATHF/Script/Pairwiseformat_Script.sh $indd $locii $nb_gen $nameIN_Stat11_format

                        Rscript $PATHF/Script/Pairwiseplot_Script.r $indd $locii $nb_gen $repeat $nameIN_Stat11_plot $nameOUT_Stat11_plot




                        python Stats_Script.py "$indd" "$locii" "$nb_gen" "$repeat" "$nameIN_FH11_Statscript" "$nameOUT_FH11_Statscript"

                        Rscript $PATHF/Script/Pairwisetable_Script.r $indd $locii $nb_gen $repeat $nameIN_FH11_table $nameOUT_FH11_table

                        bash $PATHF/Script/Pairwiseformat_Script.sh $indd $locii $nb_gen $nameIN_FH11_format

                        Rscript $PATHF/Script/Pairwiseplot_Script.r $indd $locii $nb_gen $repeat $nameIN_FH11_plot $nameOUT_FH11_plot




                        python Stats_Script.py "$indd" "$locii" "$nb_gen" "$repeat" "$nameIN_Stat01_Statscript" "$nameOUT_Stat01_Statscript"

                        Rscript $PATHF/Script/Pairwisetable_Script.r $indd $locii $nb_gen $repeat $nameIN_Stat01_table $nameOUT_Stat01_table

                        bash $PATHF/Script/Pairwiseformat_Script.sh $indd $locii $nb_gen $nameIN_Stat01_format

                        Rscript $PATHF/Script/Pairwiseplot_Script.r $indd $locii $nb_gen $repeat $nameIN_Stat01_plot $nameOUT_Stat01_plot




                        python Stats_Script.py "$indd" "$locii" "$nb_gen" "$repeat" "$nameIN_FH01_Statscript" "$nameOUT_FH01_Statscript"

                        Rscript $PATHF/Script/Pairwisetable_Script.r $indd $locii $nb_gen $repeat $nameIN_FH01_table $nameOUT_FH01_table

                        bash $PATHF/Script/Pairwiseformat_Script.sh $indd $locii $nb_gen $nameIN_FH01_format

                        Rscript $PATHF/Script/Pairwiseplot_Script.r $indd $locii $nb_gen $repeat $nameIN_FH01_plot $nameOUT_FH01_plot




                        python Stats_Script.py "$indd" "$locii" "$nb_gen" "$repeat" "$nameIN_Stat10_Statscript" "$nameOUT_Stat10_Statscript"

                        Rscript $PATHF/Script/Pairwisetable_Script.r $indd $locii $nb_gen $repeat $nameIN_Stat10_table $nameOUT_Stat10_table

                        bash $PATHF/Script/Pairwiseformat_Script.sh $indd $locii $nb_gen $nameIN_Stat10_format

                        Rscript $PATHF/Script/Pairwiseplot_Script.r $indd $locii $nb_gen $repeat $nameIN_Stat10_plot $nameOUT_Stat10_plot




                        python Stats_Script.py "$indd" "$locii" "$nb_gen" "$repeat" "$nameIN_FH10_Statscript" "$nameOUT_FH10_Statscript"

                        Rscript $PATHF/Script/Pairwisetable_Script.r $indd $locii $nb_gen $repeat $nameIN_FH10_table $nameOUT_FH10_table

                        bash $PATHF/Script/Pairwiseformat_Script.sh $indd $locii $nb_gen $nameIN_FH10_format

                        Rscript $PATHF/Script/Pairwiseplot_Script.r $indd $locii $nb_gen $repeat $nameIN_FH10_plot $nameOUT_FH10_plot



# -------------------------- Allelic frequencies --------------------------


                        python Stats_Script.py "$indd" "$locii" "$nb_gen" "$repeat" "$nameIN_PFA00_Statscript" "$nameOUT_PFA00_Statscript"

                        Rscript $PATHF/Script/Pairwisetable_Script.r $indd $locii $nb_gen $repeat $nameIN_PFA00_table $nameOUT_PFA00_table

                        bash $PATHF/Script/Pairwiseformat_Script.sh $indd $locii $nb_gen $nameIN_PFA00_format


# -------------------------- Create the RikikiGeneral.txt file --------------------------

                        bash $PATHF/Script/RikikiMakeFile_Script.sh $indd $locii $nb_gen

                fi
            done
        done
    done


    rm -f $PATHF/Resultats/*.tmp2

    rm -f $PATHF/Resultats/Matrice/*.tmp2

    rm -f $PATHF/Resultats/*.tmp


    mv $PATHF/Resultats/Rikiki_General.txt $PATHF/Resultats/Donnees

    mv $PATHF/Resultats/ID_N.csv $PATHF/Resultats/Donnees

    mv $PATHF/Resultats/ID_S.csv $PATHF/Resultats/Donnees

    mv $PATHF/Resultats/ID_SS.csv $PATHF/Resultats/Donnees

    mv $PATHF/Resultats/myRUN* $PATHF/Resultats/Donnees/


    rm -f $PATHF/Resultats/*.gen

    rm -f $PATHF/Resultats/*.txt

    rm -f $PATHF/Resultats/*.csv

    rm -rf $PATHF/Resultats/Matrice

    rm -f $PATHF/Script/Rplots.pdf


fi

mv $PATHF/savepop0.pop $PATHF/Resultats/


ENDTIME=$(date +%s)

echo "Simulation run for  $(($ENDTIME - $STARTTIME)) secondes..."

date

echo "End of simulations !"

#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================