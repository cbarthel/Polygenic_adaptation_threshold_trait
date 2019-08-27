# -*- coding: utf-8 -*-


#Clément BARTHÉLÉMY, Janvier 2016, Ifremer Centre Atlantique - Station de La Tremblade, Master Modélisation en Écologie (MODE) - AgroCampus Ouest/OSUR


#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================


rm(list=ls()) #Clear


#----------------------- Functions -----------------------#


Cal_Confusion_Matrix <- function(data,repeti){

    speci=vector()
    sensi=vector()
    Confusion_Matrix=data.frame()

    nb_pair=(nb_loci*(nb_loci-1))/2

    nb_fullcodant=(10*9)/2

    nb_neutre=((nb_loci-10)*(nb_loci-11))/2

    nb_halfcodant=nb_pair-(nb_fullcodant+nb_neutre)

    nb_codant=nb_halfcodant+nb_fullcodant


    for (i in seq(1,repeti,1)){

        data1=subset(data, data$Freq>=i)

        data2=subset(data1, data1$WEI==0) #neutre

        data3=subset(data1, data1$WEI>0) #codant

        FP=as.numeric(nrow(data2))

        TP=as.numeric(nrow(data3))

        TN=as.numeric(nb_neutre-FP)

        FN=as.numeric(nb_codant-TP)

        Matt_Coeff=((TP*TN)-(FP*FN))/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))) #Calcul du coeeficient de matthews

        confu_mat=c(TP,TN,FP,FN,Matt_Coeff,i)

        Confusion_Matrix=rbind(Confusion_Matrix,confu_mat) #Remplir la matrice de confusion

        tpr=TP/(TP+FN) #true positive rate
        tnr=TN/(TN+FP) #true negative rate

        sensi[i]=tpr #sensibility - true positive rate
        speci[i]=1-tnr #specificity - false positive rate
    }

    colnames(Confusion_Matrix) <- c("true_positive","true_negative","false_positive","false_negative","matt_coeff","repetition")

    sensi=t(t(sensi))
    speci=t(t(speci))

    name=seq(1,repeti,1)

    DATA_ROC_Curve=cbind(sensi,speci)
    DATA_ROC_Curve=cbind(DATA_ROC_Curve,name)

    colnames(DATA_ROC_Curve)=(c("SENSI","SPECI","Selection_Rep"))

    DATA_ROC_Curve=as.data.frame(DATA_ROC_Curve)

    return(list(Confusion_Matrix,DATA_ROC_Curve))
}

Define_Thresholds_ROC <- function(data_Confusion_Matrix,Max_Specificity,Optimized_Point){

    data_ROC=as.data.frame(data_Confusion_Matrix[2])
    MC=as.data.frame(data_Confusion_Matrix[1])

    T3=max(MC$matt_coeff,na.rm=TRUE)
    MCsub=subset(MC, MC$matt_coeff==T3)
    T3=MCsub$repetition

    data_ROC=as.data.frame(data_ROC)

    euc.dist <- function(p1, p2) abs((p1[1] - p2[1]) + (p1[2] - p2[2]))

    Distance=vector()

    for (i in seq(1,repi,1)){

        Observed_Point=c(data_ROC[i,1],data_ROC[i,2])
        DEuc=euc.dist(Optimized_Point, Observed_Point)
        Distance[i]=DEuc
    }


    N_Selec_Index=seq(1,repi,1)

    Selected_Limits=vector()

    #Find the most optimized choice

    Data_Distance=cbind(N_Selec_Index,Distance)

    Data_Distance=as.data.frame(Data_Distance)

    Min_Distance=min(Data_Distance$Distance)

    Optimized_Selection_Limit=subset(Data_Distance, Data_Distance$Distance==Min_Distance)

    Optimized_Selection_Limit=as.numeric(Optimized_Selection_Limit)

    Selected_Limits[1]=Optimized_Selection_Limit[1]


    #Find the most pertinent choice

    data_ROC=as.data.frame(data_ROC)

    Powerfull_Selected_Limit=subset(data_ROC, data_ROC$SPECI<=0+Max_Specificity)

    Powerfull_Selected_Limit=min(Powerfull_Selected_Limit$Selection_Rep)

    Selected_Limits[2]=Powerfull_Selected_Limit


    Data_Selected_Limits_Optimized=subset(data_ROC, data_ROC$Selection_Rep==Selected_Limits[1])

    Data_Selected_Limits_Powerfull=subset(data_ROC, data_ROC$Selection_Rep==Selected_Limits[2])

    Data_Selected_Limits=rbind(Data_Selected_Limits_Optimized,Data_Selected_Limits_Powerfull)


    return(Data_Selected_Limits)

}

Select_cluster <- function(data,nb_cluster) {

    tline=max(unique(data$gen)) #Count the number of lines

    rmeds=vector()

    for (i in seq(1,nb_cluster,1)){

        subdatameds=subset(Data_Medoids, Data_Medoids$IDClust==i)
        r=subdatameds$ValueMeds[tline]-subdatameds$ValueMeds[1]

        rmeds[i]=r
    }

    Numclust=seq(1,nb_cluster,1)
    datapos=cbind(Numclust,rmeds)
    as.data.frame(datapos)

    Mxdatapos=max(datapos[,2])

    positiveratiomed <- subset(datapos, datapos[,2]==Mxdatapos) #Create a new matrix with only the positive growth rate
    croissclust=as.integer(positiveratiomed[,1])

    return(croissclust)
}

Extract_DATA_medoids <- function(k,clust) {
    mdo=attributes(clust)$centroids[[k]]
    mdo=t(mdo)
    mdo=t(mdo)
    mdo=as.vector(mdo)

    return(mdo)
}

Create_DATAFRAME_Meds <- function(data,nb_cluster,nbgen) {

    gen=seq(1,nbgen,1)
    gen=t(gen)
    gen=t(gen)

    #Initialization of the dataframe with the data of the first cluster

    IDClust=as.integer(1)
    ValueMeds=Extract_DATA_medoids(1,data)

    clus=cbind(IDClust,ValueMeds)

    DATAFRAME_Meds=cbind(gen,clus)
    colnames(DATAFRAME_Meds)=(c("gen","IDClust","ValueMeds"))

    #Completing the dataframe with the data of the other cluster

    for (i in seq(2,nb_cluster,1)){

        IDClust=as.integer(i)
        ValueMeds=Extract_DATA_medoids(i,data)

        clus=cbind(IDClust,ValueMeds)

        DATAFRAME_Clust=cbind(gen,clus)

        DATAFRAME_Meds=rbind(DATAFRAME_Meds,DATAFRAME_Clust)

    }

    DATAFRAME_Meds=as.data.frame(DATAFRAME_Meds)



    return(DATAFRAME_Meds)
}

Create_ID_haplo <- function(IDpair,haplo){
    ID=vector()

    for (i in seq(1,length(IDpair),1)){
        t=as.character(IDpair[i])
        t=paste(t,haplo,sep="_")
        # t=as.integer(t)
        ID[i]=t
    }
    return(ID)
}

Create_WEI_haplo <- function(IDpair,haplo){

    vecWEI=vector()

    for (i in seq(1,length(IDpair),1)){

        if (IDpair[i] %in% IDSS){
        vecWEI[i]=18
        } else if (IDpair[i] %in% IDS){
        vecWEI[i]=18
        } else {
        vecWEI[i]=0
        }
    }

    WEI=rep(vecWEI,4)
    return(WEI)
}

Create_Index_Selected_Loci <- function(ID, data){
  
    charac1=paste(as.character(ID[1,1]),"00",sep="_")
    charac2=paste(as.character(ID[1,1]),"11",sep="_")
    charac3=paste(as.character(ID[1,1]),"10",sep="_")
    charac4=paste(as.character(ID[1,1]),"01",sep="_")


    Mat_ID=subset(data, data$ID == charac1 | data$ID == charac2 | data$ID == charac3 | data$ID == charac4)


    if (max(Mat_ID$Freq)==0){
        Index_Selected_Loci=Mat_ID[1,]
    } else if (nrow(subset(Mat_ID, Mat_ID$Freq==max(Mat_ID$Freq))) !=1 & max(Mat_ID$Freq) !=0) {
        Mat_ID=subset(Mat_ID, Mat_ID$Freq==max(Mat_ID$Freq))
        Index_Selected_Loci=Mat_ID[1,]
    } else {
        Index_Selected_Loci=subset(Mat_ID, Mat_ID$Freq==max(Mat_ID$Freq))
    }


    for (i in seq(2,nrow(ID),1)){

        charac1=paste(as.character(ID[i,1]),"00",sep="_")
        charac2=paste(as.character(ID[i,1]),"11",sep="_")
        charac3=paste(as.character(ID[i,1]),"10",sep="_")
        charac4=paste(as.character(ID[i,1]),"01",sep="_")


        Mat_ID=subset(data, data$ID == charac1 | data$ID == charac2 | data$ID == charac3 | data$ID == charac4)


        if (max(Mat_ID$Freq)==0){
            SLI=Mat_ID[1,]
        } else if (nrow(subset(Mat_ID, Mat_ID$Freq==max(Mat_ID$Freq))) !=1 & max(Mat_ID$Freq) !=0) {
            SLI=subset(Mat_ID, Mat_ID$Freq==max(Mat_ID$Freq))
            SLI=SLI[1,]
        } else {
            SLI=subset(Mat_ID, Mat_ID$Freq==max(Mat_ID$Freq))
        }

        Index_Selected_Loci=rbind(Index_Selected_Loci,SLI)

    }


    if (nrow(Index_Selected_Loci) - nrow(ID)==0){
        return(Index_Selected_Loci)
    } else {
        return("Index_Selected_Loci et ID n'ont pas la même dimension !")
    }


}

Find_Non_Retained_Loci <- function(ID,data){

    if (length(setdiff(ID,data$ID)) != 0){

        data_out=as.data.frame(setdiff(ID,data$ID))

        wei=rep(0,length(t(data_out)))
        data_out=cbind(data_out,wei,wei)
        colnames(data_out)=c("ID","Freq","WEI")
    } else {

      data_out=data.frame()

    }

    return(data_out)

}

Create_Full_Data <- function(data_coding,data_neutral,ID_neutral,ID_coding){

    N_out=Find_Non_Retained_Loci(ID_neutral,data_neutral)
    C_out=Find_Non_Retained_Loci(ID_coding,data_coding)

    if (nrow(N_out)==0){
        DATA_MF=rbind(data_coding,data_neutral, C_out)
    } else if (nrow(C_out)==0){
        DATA_MF=rbind(data_coding,data_neutral, N_out)
    } else {
        DATA_MF=rbind(data_coding, data_neutral, N_out, C_out)
    }

    DATA_MF=DATA_MF[order(DATA_MF$Freq,DATA_MF$WEI),]

    DATA_MF$WEI[DATA_MF$WEI == "18"] = "Coding"
    DATA_MF$WEI[DATA_MF$WEI == "0"] = "Neutral"

    return(DATA_MF)

}

Transformation_Matrice <- function(data3,Value_Haplo){

    data3=rikiki3

    freqHaploIN=Value_Haplo

    nblineSUB=nrow(data3)

    nblineID=subset(data3, data3$ID==12)
    nblineID=nrow(nblineID)



    if (freqHaploIN==1){
        RikikiFreq=data3$FreqHaplo00
    }

    if (freqHaploIN==2){
        RikikiFreq=data3$FreqHaplo11
    }

    if (freqHaploIN==3){
        RikikiFreq=data3$FreqHaplo01
    }

    if (freqHaploIN==4){
        RikikiFreq=data3$FreqHaplo10
    }


    RikikiFreq=as.vector(RikikiFreq)



    nbgen=max(unique(data3$gen))

    nbpair=length(unique(data3$ID))

    gen=seq(1,nbgen,1) #Sequence of the generations for the plots


    PairFreq=RikikiFreq[1:nblineID] #First pair

    cpt=1

    for (i in seq((nblineID*2)+1,nblineSUB+1,nblineID)){
        start=i-nbgen
        end=i-1
        PairFreqNext=RikikiFreq[start:end]
        PairFreq=rbind(PairFreq,PairFreqNext)
        cpt=cpt+1
    }

    name=seq(1,nbpair,1)
    name=as.character(name)
    rownames(PairFreq)=name
    PairFreq=as.data.frame(PairFreq)

    return(PairFreq)

}


for (i in seq(1,1,1)){

    #====================== Packages ======================#

    library(ggplot2)
    library(dtwclust)
    library(reshape2)
    library(RColorBrewer)
    library(psych)

    #====================== Parameters ======================#
  
    i=37

    sc=read.csv("/home/clement/Documents/clement_barthelemy/SimulationPP/Script/sc.txt",sep=" ",header=FALSE)

    sc=sc[i,]

    repi=as.numeric(sc[1])
    nb_gen=as.numeric(sc[2])
    nb_indiv=as.numeric(sc[3])
    nb_loci=as.numeric(sc[4])
    nb_cluster=as.numeric(sc[5])

    num_sc=37


    #====================== Colors ======================#


    colors5 <- brewer.pal(4, "Set1")
    colors3 <- brewer.pal(9, "Greys")
    colors1 <- brewer.pal(3, "Dark2")
    colors2 <- brewer.pal(10, "Spectral")
    colors <- c(colors1[1],colors1[2])


    #====================== Création du chemin principal ======================#


    namefolder=paste("Resultats_Scenario", sep="_",num_sc)

    print(namefolder)

    pats=file.path("/home", "clement", "Documents", "clement_barthelemy", "SimulationPP", namefolder, "Donnees")


    #====================== Delete data ======================#

    # file.remove(file.path(pats, "Data_ROC_Curve.txt"))
    # file.remove(file.path(pats, "Histogram_data.txt"))
    # file.remove(file.path(pats, "Index_Selected_Loci_Haplo00.txt"))
    # file.remove(file.path(pats, "Matrice_Confusion.txt"))
    # file.remove(file.path(pats, "Selected_loci_index.txt"))
    # file.remove(file.path(pats, "Selection_Pairwise_Terminal_File.txt"))
    

    #----------------------- Data -----------------------#


    rikiki<-read.csv(file.path(pats, "Rikiki_General.txt"),sep=",", header=TRUE)
    
    IDN=subset(rikiki, rikiki$WEI==0)
    IDN=unique(IDN$ID)
    IDN=sort(IDN)
    
    write.table(IDN, file = file.path(pats, "ID_N.csv"), append=FALSE, row.names = FALSE, col.names = FALSE,sep = ",")
    
    
    IDS=subset(rikiki, rikiki$WEI==1)
    IDS=unique(IDS$ID)
    IDS=sort(IDS)
    
    write.table(IDS, file = file.path(pats, "ID_S.csv"), append=FALSE, row.names = FALSE, col.names = FALSE,sep = ",")
    
    IDSS=subset(rikiki, rikiki$WEI==2)
    IDSS=unique(IDSS$ID)
    IDSS=sort(IDSS)
    
    write.table(IDSS, file = file.path(pats, "ID_SS.csv"), append=FALSE, row.names = FALSE, col.names = FALSE,sep = ",")
    
    
    
    
    IDN<-read.csv(file.path(pats, "ID_N.csv"),sep=",", header=FALSE)
    IDS<-read.csv(file.path(pats, "ID_S.csv"),sep=",", header=FALSE)
    IDSS<-read.csv(file.path(pats, "ID_SS.csv"),sep=",", header=FALSE)
    

    
    #----------------------- Creation list match -----------------------#


    for (numrep in seq(5,5,5)){#max(rikiki$rep)

        print(numrep)

        IDN=as.vector(t(IDN))
        IDS=as.vector(t(IDS))
        IDSS=as.vector(t(IDSS))


        rikiki2=subset(rikiki, rikiki$rep==numrep)
        rikiki3=rikiki2[order(rikiki2$ID,rikiki2$gen),]


        nb_gen=max(unique(rikiki3$gen))
        IDpair=unique(rikiki3$ID)


        WEI=Create_WEI_haplo(IDpair)


        ID_00=Create_ID_haplo(IDpair,"00")
        ID_11=Create_ID_haplo(IDpair,"11")
        ID_10=Create_ID_haplo(IDpair,"10")
        ID_01=Create_ID_haplo(IDpair,"01")

        IDpair=c(ID_00,ID_11,ID_01,ID_10)


        IDpair=as.data.frame(IDpair)
        Data=cbind(IDpair,WEI)

        Coding=subset(Data, WEI == 18)
        Neutral=subset(Data, WEI == 0)

        if ((nrow(Coding) + nrow(Neutral)) - nrow(Data)!=0){
            print("ATTENTION --->  codant+neutre!=Data")
        }


        write.table(Coding$IDpair, file = file.path(pats, "ID_Coding.txt"), append=FALSE, row.names = FALSE, col.names = FALSE, sep=",")

        write.table(Neutral$IDpair, file = file.path(pats, "ID_Neutral.txt"), append=FALSE, row.names = FALSE, col.names = FALSE, sep=",")


        res00=Transformation_Matrice(rikiki3,1)
        res11=Transformation_Matrice(rikiki3,2)
        res10=Transformation_Matrice(rikiki3,3)
        res01=Transformation_Matrice(rikiki3,4)

        PairFreq=rbind(res00,res11,res01,res10)

        e=tsclust(series = PairFreq, type = "partitional", k = nb_cluster, distance = "Euclidean", dist.method="L1", centroid = "pam")
        
        plot(e)

        IDcluster=attributes(e)$cluster
        IDcluster=as.data.frame(IDcluster)

        Data=cbind(Data, IDcluster)

        Data_Medoids=Create_DATAFRAME_Meds(e,nb_cluster,nb_gen)
        Clust_Croiss=Select_cluster(Data_Medoids,nb_cluster)
        print(Clust_Croiss)



        if (numrep==1){

            DataMed_GroundZero=subset(Data_Medoids, Data_Medoids$IDClust==Clust_Croiss)

            repit=as.data.frame(rep(numrep,nrow(DataMed_GroundZero)))
            colnames(repit)=c("rep")

            DataMed_GroundZero=cbind(DataMed_GroundZero,repit)

        } else {

            DataMed_choosen=subset(Data_Medoids, Data_Medoids$IDClust==Clust_Croiss)

            repit=as.data.frame(rep(numrep,nrow(DataMed_choosen)))
            colnames(repit)=c("rep")

            DataMed_choosen=cbind(DataMed_choosen,repit)

            DataMed_GroundZero=rbind(DataMed_GroundZero, DataMed_choosen)

        }



        Data_Clust=subset(Data, Data$IDclust==Clust_Croiss)



        List_match=as.character(unique(Data_Clust$IDpair))

        List_match=t(as.vector(List_match))

        write.table(List_match, file = file.path(pats, "Selected_loci_index.txt"), append=TRUE, row.names = FALSE, col.names = FALSE, sep=",")

    }

    # colors4 <- brewer.pal(9, "Set1")

    # ggplot() + geom_path(data=DataMed_GroundZero, aes(gen,ValueMeds, group=rep), size=2) + theme_gray() + xlab("Weight") + ylab("Occurences") + theme_gray() + theme(axis.text=element_text(size=15,face="bold",family="Times",colour="black"), axis.title=element_text(size=18,face="bold",family="Times"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.title=element_blank(),legend.position="none")


    #----------------------- Step 1 -----------------------#


    no_col <- max(count.fields(file.path(pats, "Selected_loci_index.txt"), sep = ","))
    list_SLI=read.csv(file.path(pats, "Selected_loci_index.txt"),sep=",",fill=TRUE,col.names=1:no_col,header=FALSE)

    ID_C <- read.table(file.path(pats, "ID_Coding.txt"), sep = ",")
    ID_N <- read.table(file.path(pats, "ID_Neutral.txt"), sep = ",")



    list_SLI[list_SLI==""] <- NA
    list_SLI[list_SLI=="NA"] <- NA

    list_SLI <- list_SLI[!is.na(list_SLI)]
    tableSLI_Haplo=table(list_SLI)
    tableSLI=as.data.frame(tableSLI_Haplo)
    colnames(tableSLI)=c("ID","Freq")



    ID_C=t(as.vector(ID_C))
    DATA_C <- tableSLI[tableSLI$ID %in% ID_C,]
    wei_C=rep(18, nrow(DATA_C))
    wei_C=t(t(wei_C))
    DATA_C=cbind(DATA_C,wei_C)
    colnames(DATA_C)=(c("ID","Freq", "WEI"))



    ID_N=t(as.vector(ID_N))
    DATA_N <- tableSLI[tableSLI$ID %in% ID_N,]
    wei_N=rep(0, nrow(DATA_N))
    wei_N=t(t(wei_N))
    DATA_N=cbind(DATA_N,wei_N)
    colnames(DATA_N)=(c("ID","Freq", "WEI"))


    #----------------------- Trouver les paires qui ne sont jamais retenues et créer les données avec les quatres haplo pour toutes les paires -----------------------#


    DATA_Full=Create_Full_Data(DATA_C,DATA_N,ID_N,ID_C)


    #----------------------- Sélectionner seulement l'haplotype qui augmente le plus en fréquence sur chaque paire -----------------------#


    ID_S <- read.table(file.path(pats, "ID_S.csv"), sep = ",")
    ID_SS <- read.table(file.path(pats, "ID_SS.csv"), sep = ",")
    ID_N <- read.table(file.path(pats, "ID_N.csv"), sep = ",")

    ID_C=rbind(ID_S,ID_SS)

    Index_Selected_Loci_N=Create_Index_Selected_Loci(ID_N,DATA_Full)
    Index_Selected_Loci_C=Create_Index_Selected_Loci(ID_C,DATA_Full)

    #----------------------- Créer les données finales (Index_Selected_Loci) -----------------------#


    DATA_SLI=rbind(Index_Selected_Loci_C,Index_Selected_Loci_N)
    colnames(DATA_SLI)=c("ID", "Freq", "WEI")

    DATA_SLI=DATA_SLI[order(DATA_SLI$Freq,DATA_SLI$WEI),]


    write.table(DATA_SLI, file = file.path(pats, "Index_Selected_Loci.txt"),sep=",", append=FALSE, row.names = FALSE, col.names = TRUE)

    if (nrow(DATA_SLI) - nrow(ID_N) - nrow(ID_C)!=0){
        print("ATTENTION!!!!")
    }


    #------------------ Calculer les variables de sorties : matrice de confusion -----------------------#


    data1=read.csv(file.path(pats, "Index_Selected_Loci.txt"),sep=",",header=TRUE)
    data1=as.data.frame(data1)
    data1=data1[order(data1$Freq,data1$WEI),]
    data1$WEI=as.character(data1$WEI)

    data1$WEI[data1$WEI=="Coding"] <- "18"
    data1$WEI[data1$WEI=="Neutral"] <- "0"

    data1$WEI=as.numeric(data1$WEI)

    Matrice_Confusion=Cal_Confusion_Matrix(data1,repi)


    write.table(Matrice_Confusion[1], file = file.path(pats, "Matrice_Confusion.txt"), append=FALSE, row.names = FALSE, col.names = TRUE, sep = ",")


    Matrice_confu=as.data.frame(Matrice_Confusion[1])
    MCC=max(Matrice_confu$matt_coeff, na.rm=TRUE)

    MCC_select1=subset(Matrice_confu, Matrice_confu$matt_coeff==MCC)

    MCC_pop=as.numeric(min(MCC_select1$repetition))

    MCC_select2=subset(MCC_select1, MCC_select1$repetition==MCC_pop)

    write.table(MCC_select2, file = "/home/clement/Documents/clement_barthelemy/SimulationPP/MCC_Results_36.txt", sep=",", append=TRUE, row.names = FALSE, col.names = FALSE)


    #----------------------- Plot -----------------------#


    colors4 <- brewer.pal(9, "Set1")

    ggplot() + geom_jitter(data=DATA_SLI, aes(WEI,Freq, group=WEI), size=2) + geom_violin(data=DATA_SLI, aes(WEI,Freq,group=WEI, fill=WEI)) + scale_fill_manual(values=c(colors4[9],colors4[9])) +  geom_hline(yintercept=MCC_pop, color = "red", size=2) + theme_gray() + xlab("Weight") + ylab("Occurences") + theme_gray() + theme(axis.text=element_text(size=15,face="bold",family="Times",colour="black"), axis.title=element_text(size=18,face="bold",family="Times"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.title=element_blank(),legend.position="none")

    ggsave(filename=file.path(pats, "Histo_MCC.pdf"), device="pdf", width=350, height=250, units="mm", dpi=300, limitsize=FALSE)

}



