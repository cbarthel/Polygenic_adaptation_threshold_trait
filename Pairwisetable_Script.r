# -*- coding: utf-8 -*-


# Script for manuscript "Detecting Polygenic Adaptation For A Threshold Trait In Experimental Populations Using Time-Series Data"


#Clément BARTHÉLÉMY, La Rochelle University, Station Ifremer de La Tremblade


#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  ###################################################
#====================================================================================================================


rm(list=ls()) #Clear


#====================== Path ======================#


pats=file.path("/home", "clement", "Documents", "SimulationPP")


#====================== Variables ======================#


args <- commandArgs(TRUE) #Import variable from the Main_Bash_Launch.sh script

nb_loci=as.numeric(args[2])
nb_indiv=as.numeric(args[1])
repi=as.numeric(args[4])
nb_gen=as.numeric(args[3])

nameIN=as.character(args[5])
nameOUT=as.character(args[6])


#====================== Import the data ======================#


for (i in seq(1,nb_gen,1)) { #Loop in order to execute this for each generation (in several diffirent file)

    extens=".tmp" #Set the name of the import file

    nfile=paste(nameIN, sep="",nb_loci) #Set the name of the import file

    nfilebis=paste(nfile, sep="=",nb_indiv,i) #Set the name of the import file

    nfilebis2=paste(nfilebis, sep="",extens) #Set the name of the import file

    ndata=file.path(pats, "Resultats", nfilebis2)

    obs1<-read.table(ndata,header=FALSE)

    obs2<-read.table(file.path(pats, "Resultats", "Pos_LambdaStat_Pairwise.txt"),header=FALSE)


#----------- Creation of the generation column -----------#


    gen=array(data = i, dim = nrow(obs1), dimnames = NULL) #Vector for stocking the values

    l1=length(gen)


#----------- Creation of the repietition column -----------#


    repetition=array(data = NA, dim = nrow(obs1), dimnames = NULL) #Vector for stocking the values

    vecrepi=array(data = NA, dim = repi, dimnames = NULL) #Vector for stocking the values

    for (j in seq(1,repi,1)){ #Stock the pattern into a vector
        vecrepi[j]=j
    }

    for (lin in seq(1,(((nb_loci*nb_loci)-nb_loci)/2)*repi,repi)){
      ii=(lin+repi-1)
      repetition[lin:ii]<-vecrepi #Stick the pattern to a vector the exact number of times
    }

    l2=length(repetition)


#----------- Creation of the locii (the first one) column -----------#


    locii=c() #Empty vector

    for (pos in seq(1,nrow(obs2),1)) { #Parse every row of the obs2 file

      hh=rep(obs2[pos,1],repi) #Repeat the value the exact number of times

      locii=append(locii,hh) #Append to the vector the repetition of the value
    }


    l3=length(locii)


#----------- Creation of the locij (the second one) column -----------#


    locij=c()

    for (pos in seq(1,nrow(obs2),1)) {

      kk=rep(obs2[pos,2],repi)

      locij=append(locij,kk)
    }


    l4=length(locij)


#----------- Creation of the LD column -----------#


    vectorLD=obs1


#----------- Creation of the locus ID -----------#


    ID=as.numeric(paste(locii, locij, sep=""))


#----------- Create the color column with the weights of the locii -----------#


    weights<-read.csv(file.path(pats, "Script", "Weights_Locii_File.txt"), header=FALSE)

    weights=as.matrix(weights)

    weights=weights[1:nb_loci]

    weights=as.matrix(weights)


    wei_loci=c() #Create an empty vector

    cpt1=1

    for (i in locii){
        wei_loci[cpt1]=weights[i,] #The value in the vector "locii" is the position on the vector weights
        cpt1=cpt1+1
    }



    wei_locj=c() #Create an empty vector

    cpt2=1

    for (j in locij){
        wei_locj[cpt2]=weights[j,] #The value in the vector "locij" is the position on the vector weights
        cpt2=cpt2+1
    }



    WEI=c() #Create an empty vector

    for (k in seq(1,length(locii),1)){
        WEI[k]=wei_loci[k]+wei_locj[k] #Just add the two column in one 
    }


#----------- Test if the vectors have the same number of row -----------#


    if (l1==l2 & l3==l4 & l1==l4) {

        gen=as.matrix(gen)  #Put the vector into a column
        repetition=as.matrix(repetition)
        locii=as.matrix(locii)
        locij=as.matrix(locij)
        ID=as.matrix(ID)
        WEI=as.matrix(WEI)
        vectorLD=as.matrix(vectorLD)


        obs3=as.matrix(cbind(gen,repetition,locii,locij,ID,WEI,vectorLD)) #assemblate the vector into a matrix


        nfilepp=paste(nameOUT, sep="",nb_loci) #Set the name of the import file

        nfileppbis=paste(nfilepp, sep="=",nb_indiv) #Set the name of the import file

        nfileppbis2=paste(nfileppbis, sep="",extens) #Set the name of the import file

        ndatapp=file.path(pats, "Resultats", nfileppbis2)

        write.table(obs3, file = ndatapp, append=TRUE, row.names = FALSE, col.names = FALSE)

    }


} #End of the loop "for"


#####################################################################################################################





print("+++++++++++++++------------------------------ | Pairwisetable_Script.r   | --- | FAIT | ------------------++++++++++++++++++++")





#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================