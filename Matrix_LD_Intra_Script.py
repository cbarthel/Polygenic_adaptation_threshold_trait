#!/usr/bin/env python 
# -*- coding: utf-8 -*-


# Script for manuscript "Detecting Polygenic Adaptation For A Threshold Trait In Experimental Populations Using Time-Series Data"


#Clément BARTHÉLÉMY, La Rochelle University, Station Ifremer de La Tremblade


#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  ###################################################
#====================================================================================================================


#====================== Path ======================#


PATHF = os.path.join('/home/clement/Documents/SimulationPP')


#====================== Packages ======================#


import math,os,time,random,sys,subprocess
import numpy


#====================== Global variables ======================#


nb_indiv = int(sys.argv[1]) #Number of individuals in the total population

nb_loci = int(sys.argv[2]) #Number of locii

nb_gen = int(sys.argv[3]) #Number of studied generations

rep = int(sys.argv[4]) #The repetition in which we are

Matrix_Freq0=int(-999)*numpy.ones((nb_gen,nb_loci)) #Matrice des LD

Matrix_Freq1=int(-999)*numpy.ones((nb_gen,nb_loci)) #Matrice des LD

pat_haplo00=numpy.array([0,0]) #Les différents haplotypes
pat_haplo11=numpy.array([1,1])
pat_haplo01=numpy.array([0,1])
pat_haplo10=numpy.array([1,0])


#====================== Functions ======================#


#-----------  Calcul of the allelic frequencies for allele 0 and allele 1  -----------#


def Allfreq(obs1,nb_indiv,frequency0,frequency1,Freq0_chrom,Freq1_chrom,FreqAll0Haplo):

    indicechrom=0
    indicevec=0
    cpt1=0
    cpt0=0
    cpt1loc=0
    cpt0loc=0

    nline=numpy.shape(obs1)[0]
    ncol=numpy.shape(obs1)[1]

    for column in range(0,ncol,1):

        for line in range(0,nline,1):

            if (obs1[line,column]==1):
                cpt1=cpt1+1
                cpt1loc=cpt1loc+1
            else:
                cpt0=cpt0+1
                cpt0loc=cpt0loc+1


        Freq0_chrom[indicechrom]=float(cpt0)/(nb_indiv)

        Freq1_chrom[indicechrom]=float(cpt1)/(nb_indiv)

        cpt1=0
        cpt0=0

        indicechrom=indicechrom+1


        if (indicechrom % 2 == 0):

            frequency0[indicevec]=float(cpt0loc)/(nb_indiv*2)
            frequency1[indicevec]=float(cpt1loc)/(nb_indiv*2)

            cpt0loc=0
            cpt1loc=0

            indicevec=indicevec+1

    indice_FreqAllHaplo=0

    for indy in range(0,nb_loci*2,2):

            FreqAll0Haplo[indice_FreqAllHaplo]=(Freq0_chrom[indy]+Freq0_chrom[indy+1])/2

            indice_FreqAllHaplo=indice_FreqAllHaplo+1

    return (frequency0, frequency1,Freq0_chrom,Freq1_chrom)


#-----------  Creation of the sequence that we use in the "Haplofreq" functions (loop on the chromosome)  -----------#


def Initial_loop_obj1(nb_loci, pos_init_column, loop_obj1):

    for i in range(0,(nb_loci-1)*2,2):
        loop_obj1[i]=pos_init_column
        loop_obj1[i+1]=pos_init_column+1

    return(loop_obj1)


def Initial_Intra_col_obj2(nb_loci, loop_obj2):

    for i in range(0,(nb_loci-1)*2,1):
        t1=i+2
        loop_obj2[i]=t1

    return(loop_Intra_obj2)


#-----------  Calcul of the haplotype frequencies  -----------#


def Haplofreq_intra(obs1, nb_loci, nb_indiv, loop_obj1, loop_Intra_obj2, freq_haplo_pop, pat_haplo00, pat_haplo11, pat_haplo01, pat_haplo10):

    nb_line_data=obs1.shape[0]

    line_matrix_haplo=0

    cpt00=0 #On initie les compteurs des haplotypes (combien d'haplotypes pareils entre deux loci)
    cpt01=0
    cpt10=0
    cpt11=0

    for pos_init_column in range(0,(nb_loci-1)*2,2):

        loop_obj1=Initial_loop_obj1(nb_loci, pos_init_column, loop_obj1)

        cpt_locus=0

        for i in range(0,len(loop_Intra_obj2),1):

            for line_obj in range(0,nb_line_data,1):

                col_obj2=int(loop_Intra_obj2[i])

                val_obj2=obs1[line_obj,col_obj2]

                # print([line_obj,col_obj2])


                col_obj1=int(loop_obj1[i])

                val_obj1=obs1[line_obj,col_obj1]


                # print(numpy.array([line_obj,col_obj1]))

                # print(numpy.array([line_obj,col_obj2]))


                haplo=numpy.array([val_obj1,val_obj2])

                if numpy.array_equal(haplo,pat_haplo00):
                    cpt00=cpt00+1 #Test pour remplir les compteurs

                if numpy.array_equal(haplo,pat_haplo11):
                    cpt11=cpt11+1

                if numpy.array_equal(haplo,pat_haplo01):
                    cpt01=cpt01+1

                if numpy.array_equal(haplo,pat_haplo10):
                    cpt10=cpt10+1


            cpt_locus=cpt_locus+1

            if(cpt_locus%2==0):

                nb_haplo=numpy.array([float(cpt00),float(cpt11),float(cpt01),float(cpt10)]) #Création d'un vecteur avec tous les compteurs

                freq_haplo_pop[line_matrix_haplo,:]=nb_haplo/(nb_indiv*2) #On calcul les fréquences et on rempli une matrice

                line_matrix_haplo=line_matrix_haplo+1

                cpt00=0 #On initie les compteurs des haplotypes (combien d'haplotypes pareils entre deux loci)
                cpt01=0
                cpt10=0
                cpt11=0

        if (len(loop_Intra_obj2)==2):
            break
        else:
            loop_Intra_obj2=loop_Intra_obj2[2:]

    return(freq_haplo_pop)


#-----------  Calcul du LD  -----------#


def LD_Calcul(Matrix_Haplo, nb_loci, FreqAll0ar, FreqAll1ar, Matrix_LD, Matrix_rsquare, Matrix_LambdaStat00, Matrix_LambdaStat11,  Matrix_LambdaStat01, Matrix_LambdaStat10, Matrix_FHPairwise00, Matrix_FHPairwise11, Matrix_FHPairwise01, Matrix_FHPairwise10):

    row_LD_matrix=0 #On commence toujours à la position (0,1)
    col_LD_matrix=1

    LambdaStat00=0 #Pour le LD
    LambdaStat11=0
    LambdaStat01=0
    LambdaStat10=0

    FHPairwise00=0 #Pour la fréquence des haplotypes
    FHPairwise11=0
    FHPairwise01=0
    FHPairwise10=0
    Produit_FreqAllelique=0


    for line in range(0,nb_loci*(nb_loci-1)-(nb_loci-1),1): #On boucle sur les lignes de la matrice des frequences


        LinkageDesi=Matrix_Haplo[line,0]*Matrix_Haplo[line,1]-Matrix_Haplo[line,2]*Matrix_Haplo[line,3] #LD calcul

        FHPairwise00=Matrix_Haplo[line,0] #Haplotype frequencies
        FHPairwise11=Matrix_Haplo[line,1]
        FHPairwise01=Matrix_Haplo[line,2]
        FHPairwise10=Matrix_Haplo[line,3]

        LambdaStat00=Matrix_Haplo[line,0]-(FreqAll0ar[row_LD_matrix]*FreqAll0ar[col_LD_matrix]) #LD calcul
        LambdaStat11=Matrix_Haplo[line,1]-(FreqAll1ar[row_LD_matrix]*FreqAll1ar[col_LD_matrix])
        LambdaStat01=Matrix_Haplo[line,2]-(FreqAll0ar[row_LD_matrix]*FreqAll1ar[col_LD_matrix])
        LambdaStat10=Matrix_Haplo[line,3]-(FreqAll1ar[row_LD_matrix]*FreqAll0ar[col_LD_matrix])

        Produit_FreqAllelique=(FreqAll0ar[row_LD_matrix]*FreqAll0ar[col_LD_matrix])

        Matrix_LD[row_LD_matrix,col_LD_matrix]=LinkageDesi #On met la valeur au bon endroit dans la matrice du LD

        Matrix_LambdaStat00[row_LD_matrix,col_LD_matrix]=LambdaStat00
        Matrix_LambdaStat11[row_LD_matrix,col_LD_matrix]=LambdaStat11
        Matrix_LambdaStat01[row_LD_matrix,col_LD_matrix]=LambdaStat01
        Matrix_LambdaStat10[row_LD_matrix,col_LD_matrix]=LambdaStat10

        Matrix_FHPairwise00[row_LD_matrix,col_LD_matrix]=FHPairwise00
        Matrix_FHPairwise11[row_LD_matrix,col_LD_matrix]=FHPairwise11
        Matrix_FHPairwise01[row_LD_matrix,col_LD_matrix]=FHPairwise01
        Matrix_FHPairwise10[row_LD_matrix,col_LD_matrix]=FHPairwise10

        Matrix_ProduitFreqAlle[row_LD_matrix,col_LD_matrix]=Produit_FreqAllelique

        col_LD_matrix=col_LD_matrix+1 #On incrémente la colonne de la matrice du LD et rsquare

        if (col_LD_matrix==nb_loci): #Si le numéro de colonne dépasse le nombre de loci, alors on doit changer de ligne
            row_LD_matrix=row_LD_matrix+1
            col_LD_matrix=row_LD_matrix+1 #La colonne de départ est toujours la ligne+1 (matrice triangle supérieure)

        if (row_LD_matrix==nb_loci-1): #Lorsque le nombre de ligne à atteint le nombre de loci, on quitte la boucle
            break

    return(Matrix_LD, Matrix_rsquare, Matrix_LambdaStat00, Matrix_LambdaStat11, Matrix_LambdaStat01, Matrix_LambdaStat10, Matrix_FHPairwise00, Matrix_FHPairwise11, Matrix_FHPairwise01, Matrix_FHPairwise10, Matrix_ProduitFreqAlle)


#====================== Import the data ======================#


nfile="myRUN_"+str(nb_loci)+"_"+str(nb_indiv)+".csv"

ndata=os.path.join(PATHF, "Resultats", nfile)


obs = numpy.genfromtxt(
    ndata,                  # file name
    skip_header=1,          # lines to skip at the top
    #skip_footer=0,          # lines to skip at the bottom
    delimiter=',',          # column delimiter
    dtype='float')            # data type
    #filling_values=0)       # fill missing values with 0calcul complexité d'un script


#====================== Initiate the loop ======================#


i=0

init=0 #First value for the second subset (conresponding to the line subset)

end=nb_indiv #Last value

repepepette=1 #Compte les répétitions

colsub= range (7,6+(nb_loci*2)+1,1) #Première subdivision pour les colonnes

fSubset1=obs[:,colsub] #Première subdivision pour les colonnes

cond=nb_gen*rep #On met une condition pour la boucle while, elle permet de parcourir tout le fichier

cpt_gen=0 #Compteur sur les générations


while i < cond:


#----------- Global variables with re-initialization -----------#


    frequency0=numpy.zeros(nb_loci)

    frequency1=numpy.zeros(nb_loci)

    Freq0_chrom=numpy.zeros(nb_loci*2)
    Freq1_chrom=numpy.zeros(nb_loci*2)

    FreqAll0Haplo=numpy.zeros(nb_loci)

    freq_haplo_pop=int(-999)*numpy.ones((nb_loci*(nb_loci-1)/2,4)) #Matrice pour stocker les fréquences des haplotypes

    loop_Intra_obj2=numpy.zeros((nb_loci-1)*2)

    loop_obj1=numpy.zeros((nb_loci-1)*2) #On veut une répétition de "pos-pos+1" pour la boucle de l'obj1 (nb_loci-1)*2 fois

    Matrix_LD=int(-999)*numpy.ones((nb_loci,nb_loci)) #Matrice des LD

    Matrix_rsquare=int(-999)*numpy.ones((nb_loci,nb_loci)) #Matrice des LD

    Matrix_Dprim=int(-999)*numpy.ones((nb_loci,nb_loci)) #Matrice des LD

    Matrix_LambdaStat00=int(-999)*numpy.ones((nb_loci,nb_loci)) #Matrice des Lambda
    Matrix_LambdaStat11=int(-999)*numpy.ones((nb_loci,nb_loci)) #Matrice des Lambda
    Matrix_LambdaStat01=int(-999)*numpy.ones((nb_loci,nb_loci)) #Matrice des Lambda
    Matrix_LambdaStat10=int(-999)*numpy.ones((nb_loci,nb_loci)) #Matrice des Lambda

    Matrix_FHPairwise00=int(-999)*numpy.ones((nb_loci,nb_loci)) #Matrice des fréquences haplotypiques
    Matrix_FHPairwise11=int(-999)*numpy.ones((nb_loci,nb_loci)) #Matrice des fréquences haplotypiques
    Matrix_FHPairwise01=int(-999)*numpy.ones((nb_loci,nb_loci)) #Matrice des fréquences haplotypiques
    Matrix_FHPairwise10=int(-999)*numpy.ones((nb_loci,nb_loci)) #Matrice des fréquences haplotypiques

    Matrix_ProduitFreqAlle=int(-999)*numpy.ones((nb_loci,nb_loci)) #Matrice du produit des fréquences alléliques

    Matrix_Haplo=int(-999)*numpy.ones((nb_loci*(nb_loci-1)-(nb_loci-1),4)) #Matrice des haplotypes


    if (cpt_gen==nb_gen):
        repepepette= repepepette+1
        cpt_gen=0


    loop_obj2=Initial_Intra_col_obj2(nb_loci, loop_Intra_obj2)

    linsub= range (init,end,1) #Lines for subset the data


    fSubset2=fSubset1[linsub,:] #Subseting the data


#----------- Functions calls whitout parellelisation -----------#


    allele=Allfreq(fSubset2,nb_indiv,frequency0,frequency1,Freq0_chrom,Freq1_chrom,FreqAll0Haplo) #Call the "Allfreq" function

    line_matfreq=cpt_gen #Line for the matrix of allele frequencies

    Matrix_Freq0[line_matfreq,]=allele[0] #Fill up the allele frequencies matrix
    Matrix_Freq1[line_matfreq,]=allele[1]


    FreqAll0ar=allele[0]
    FreqAll1ar=allele[1]


    nfile27="Matrix_AlleleFreq_File_"+str(nb_loci)+"="+str(nb_indiv)+"="+str(repepepette)+".txt"

    ndata27=os.path.join(PATHF, "Resultats", "Matrice", nfile27)

    file_mat=open(ndata27, "a")

    numpy.savetxt(
        file_mat,                                 # file name
        allele[0],                                     # array to save
        delimiter=' ',                          # column delimiter
        newline='\n')                           # new line character
        #footer='end of file',                  # file footer
        #comments='# ',                         # character to use for comments
        #header='Data generated by numpy')      # file header


    file_mat.close()


    loop_Intra_obj2=Initial_Intra_col_obj2(nb_loci, loop_Intra_obj2)

    Matrix_Haplo=Haplofreq_intra(fSubset2, nb_loci, nb_indiv, loop_obj1, loop_Intra_obj2, freq_haplo_pop, pat_haplo00, pat_haplo11, pat_haplo01, pat_haplo10)

    LD=LD_Calcul(Matrix_Haplo, nb_loci, FreqAll0ar, FreqAll1ar, Matrix_LD, Matrix_rsquare, Matrix_LambdaStat00, Matrix_LambdaStat11, Matrix_LambdaStat01, Matrix_LambdaStat10, Matrix_FHPairwise00, Matrix_FHPairwise11, Matrix_FHPairwise01, Matrix_FHPairwise10)


#----------- Saving data -----------#


    nfile2="Matrix_LD_File_"+str(nb_loci)+"="+str(nb_indiv)+"="+str(repepepette)+".tmp2"

    ndata2=os.path.join(PATHF, "Resultats", "Matrice", nfile2)

    file_mat=open(ndata2, "a")

    numpy.savetxt(
        file_mat,                                 # file name
        LD[0],                                     # array to save
        fmt='%.18f',                            # formatting, 2 digits in this case
        delimiter=' ',                          # column delimiter
        newline='\n')                           # new line character
        #footer='end of file',                  # file footer
        #comments='# ',                         # character to use for comments
        #header='Data generated by numpy')      # file header


    file_mat.close()




    nfile3="Matrix_LambdaStat00_File_"+str(nb_loci)+"="+str(nb_indiv)+"="+str(repepepette)+".tmp2"

    ndata3=os.path.join(PATHF, "Resultats", "Matrice", nfile3)

    file_mat=open(ndata3, "a")

    numpy.savetxt(
        file_mat,                                 # file name
        LD[2],                                     # array to save
        fmt='%.18f',                            # formatting, 2 digits in this case
        delimiter=' ',                          # column delimiter
        newline='\n')                           # new line character
        #footer='end of file',                  # file footer
        #comments='# ',                         # character to use for comments
        #header='Data generated by numpy')      # file header


    file_mat.close()




    nfile4="Matrix_LambdaStat11_File_"+str(nb_loci)+"="+str(nb_indiv)+"="+str(repepepette)+".tmp2"

    ndata4=os.path.join(PATHF, "Resultats", "Matrice", nfile4)

    file_mat=open(ndata4, "a")

    numpy.savetxt(
        file_mat,                                 # file name
        LD[3],                                     # array to save
        fmt='%.18f',                            # formatting, 2 digits in this case
        delimiter=' ',                          # column delimiter
        newline='\n')                           # new line character
        #footer='end of file',                  # file footer
        #comments='# ',                         # character to use for comments
        #header='Data generated by numpy')      # file header


    file_mat.close()




    nfile5="Matrix_LambdaStat01_File_"+str(nb_loci)+"="+str(nb_indiv)+"="+str(repepepette)+".tmp2"

    ndata5=os.path.join(PATHF, "Resultats", "Matrice", nfile5)

    file_mat=open(ndata5, "a")

    numpy.savetxt(
        file_mat,                                 # file name
        LD[4],                                     # array to save
        fmt='%.18f',                            # formatting, 2 digits in this case
        delimiter=' ',                          # column delimiter
        newline='\n')                           # new line character
        #footer='end of file',                  # file footer
        #comments='# ',                         # character to use for comments
        #header='Data generated by numpy')      # file header


    file_mat.close()




    nfile6="Matrix_LambdaStat10_File_"+str(nb_loci)+"="+str(nb_indiv)+"="+str(repepepette)+".tmp2"

    ndata6=os.path.join(PATHF, "Resultats", "Matrice", nfile6)

    file_mat=open(ndata6, "a")

    numpy.savetxt(
        file_mat,                                 # file name
        LD[5],                                     # array to save
        fmt='%.18f',                            # formatting, 2 digits in this case
        delimiter=' ',                          # column delimiter
        newline='\n')                           # new line character
        #footer='end of file',                  # file footer
        #comments='# ',                         # character to use for comments
        #header='Data generated by numpy')      # file header


    file_mat.close()




    nfile7="Matrix_FHPairwise00_File_"+str(nb_loci)+"="+str(nb_indiv)+"="+str(repepepette)+".tmp2"

    ndata7=os.path.join(PATHF, "Resultats", "Matrice", nfile7)

    file_mat=open(ndata7, "a")

    numpy.savetxt(
        file_mat,                                 # file name
        LD[6],                                     # array to save
        fmt='%.18f',                            # formatting, 2 digits in this case
        delimiter=' ',                          # column delimiter
        newline='\n')                           # new line character
        #footer='end of file',                  # file footer
        #comments='# ',                         # character to use for comments
        #header='Data generated by numpy')      # file header


    file_mat.close()




    nfile8="Matrix_FHPairwise11_File_"+str(nb_loci)+"="+str(nb_indiv)+"="+str(repepepette)+".tmp2"

    ndata8=os.path.join(PATHF, "Resultats", "Matrice", nfile8)

    file_mat=open(ndata8, "a")

    numpy.savetxt(
        file_mat,                                 # file name
        LD[7],                                     # array to save
        fmt='%.18f',                            # formatting, 2 digits in this case
        delimiter=' ',                          # column delimiter
        newline='\n')                           # new line character
        #footer='end of file',                  # file footer
        #comments='# ',                         # character to use for comments
        #header='Data generated by numpy')      # file header


    file_mat.close()




    nfile9="Matrix_FHPairwise01_File_"+str(nb_loci)+"="+str(nb_indiv)+"="+str(repepepette)+".tmp2"

    ndata9=os.path.join(PATHF, "Resultats", "Matrice", nfile9)

    file_mat=open(ndata9, "a")

    numpy.savetxt(
        file_mat,                                 # file name
        LD[8],                                     # array to save
        fmt='%.18f',                            # formatting, 2 digits in this case
        delimiter=' ',                          # column delimiter
        newline='\n')                           # new line character
        #footer='end of file',                  # file footer
        #comments='# ',                         # character to use for comments
        #header='Data generated by numpy')      # file header


    file_mat.close()




    nfile10="Matrix_FHPairwise10_File_"+str(nb_loci)+"="+str(nb_indiv)+"="+str(repepepette)+".tmp2"

    ndata10=os.path.join(PATHF, "Resultats", "Matrice", nfile10)

    file_mat=open(ndata10, "a")

    numpy.savetxt(
        file_mat,                                 # file name
        LD[9],                                     # array to save
        fmt='%.18f',                            # formatting, 2 digits in this case
        delimiter=' ',                          # column delimiter
        newline='\n')                           # new line character
        #footer='end of file',                  # file footer
        #comments='# ',                         # character to use for comments
        #header='Data generated by numpy')      # file header


    file_mat.close()




    nfile11="Matrix_ProduitFreqAlle00_File_"+str(nb_loci)+"="+str(nb_indiv)+"="+str(repepepette)+".tmp2"

    ndata11=os.path.join(PATHF, "Resultats", "Matrice", nfile11)

    file_mat=open(ndata11, "a")

    numpy.savetxt(
        file_mat,                                 # file name
        LD[10],                                     # array to save
        fmt='%.18f',                            # formatting, 2 digits in this case
        delimiter=' ',                          # column delimiter
        newline='\n')                           # new line character
        #footer='end of file',                  # file footer
        #comments='# ',                         # character to use for comments
        #header='Data generated by numpy')      # file header


    file_mat.close()


#----------- Incrementation des variables -----------#


    init=init+nb_indiv

    end=end+nb_indiv

    i=i+1

    cpt_gen=cpt_gen+1


#####################################################################################################################


#----------- Saving data -----------#


nfile12="Matrix_FreqAll_0_File_"+str(nb_loci)+"="+str(nb_indiv)+"="+str(rep)+".txt"

ndata12=os.path.join(PATHF, "Resultats", "Matrice", nfile12)

file_mat=open(ndata12, "a")

numpy.savetxt(
    file_mat,                                 # file name
    Matrix_Freq0,                                     # array to save
    fmt='%.18f',                            # formatting, 2 digits in this case
    delimiter=' ',                          # column delimiter
    newline='\n')                           # new line character
    #footer='end of file',                  # file footer
    #comments='# ',                         # character to use for comments
    #header='Data generated by numpy')      # file header

file_mat.close()




nfile13="Matrix_FreqAll_1_File_"+str(nb_loci)+"="+str(nb_indiv)+"="+str(rep)+".txt"

ndata13=os.path.join(PATHF, "Resultats", "Matrice", nfile13)

file_mat=open(ndata13, "a")

numpy.savetxt(
    file_mat,                                 # file name
    Matrix_Freq1,                                     # array to save
    fmt='%.18f',                            # formatting, 2 digits in this case
    delimiter=' ',                          # column delimiter
    newline='\n')                           # new line character
    #footer='end of file',                  # file footer
    #comments='# ',                         # character to use for comments
    #header='Data generated by numpy')      # file header


file_mat.close()


#####################################################################################################################





print "+++++++++++++++----------------------------------- | Matrix_LD_Intra_Script.py| --- | FAIT | ------------------+++++++++++++++++++++"





#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================