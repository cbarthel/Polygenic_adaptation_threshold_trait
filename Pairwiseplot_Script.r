# -*- coding: utf-8 -*-


# Script for manuscript "Detecting Polygenic Adaptation For A Threshold Trait In Experimental Populations Using Time-Series Data"


#Clément BARTHÉLÉMY, La Rochelle University, Station Ifremer de La Tremblade


#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================


rm(list=ls()) #Clear


#====================== Path ======================#


pats=file.path("/home", "clement", "Documents", "SimulationPP")


#====================== Packages ======================#

library(ggplot2)
library(reshape2)
library(RColorBrewer)


#====================== Colors ======================#


colors <- brewer.pal(10, "Spectral")
colors2 <- brewer.pal(4, "Paired")
colors3 <- brewer.pal(9, "Greys")
colors4 <- brewer.pal(9, "Set1")


#====================== Variable ======================#


args <- commandArgs(TRUE) #Import variable from the Main_Bash_Launch.sh script

nb_loci=as.numeric(args[2])
nb_indiv=as.numeric(args[1])
repi=as.numeric(args[4])
nb_gen=as.numeric(args[3])

nameIN=as.character(args[5])
nameOUT=as.character(args[6])


#====================== Import the data ======================#


extens=".csv" #Set the name of the import file

nfile=paste(nameIN, sep="",nb_loci) #Set the name of the import file

nfilebis=paste(nfile, sep="=",nb_indiv) #Set the name of the import file

nfilebis2=paste(nfilebis, sep="",extens) #Set the name of the import file

ndata=file.path(pats, "Resultats", nfilebis2) #Set the name of the import file

obs1<-read.csv(ndata,header=TRUE)

obs2=obs1[order(obs1$rep,obs1$gen),] #Order obs1 by rep and gen


#====================== Ploting LD pairwise ======================#


#----------------------- Export the table with all the informations -----------------------#


rikiki=obs2[order(obs2$rep,obs2$locii,obs2$locij,obs2$gen),]


cptt=1 #Count the repetition


for (i in seq(1,nrow(rikiki),nrow(rikiki)/repi)){ #Start a loop for all the repetition (1-2-3-4 ...)

    ii=i+(nrow(rikiki)/repi)-1

    rikiki2<-subset(rikiki[i:ii,]) #Cutting the data (rep, gen, z, sum_g)

    rikiki3=rikiki2[order(rikiki2$WEI,rikiki2$ID,rikiki2$gen),] #Subsetting the data by repetition and order them 


#----------------------- Plot the Stat (Pairwise) for all the generation by repetition -----------------------#


    WEI_N<-rikiki3[which(rikiki3$WEI==0),] ## ATTENTION au poids des locus à modifier le cas échéant !
    WEI_S<-rikiki3[which(rikiki3$WEI==9),]
    WEI_SS<-rikiki3[which(rikiki3$WEI==18),]

    rikiki3$rep <- NULL
    rikiki3$locii <- NULL
    rikiki3$locij <- NULL


    rikiki3.m=melt(rikiki3, id.vars =c("gen", "ID", "WEI"))


    write.table(rikiki3.m, file = file.path(pats, "Resultats", "rikiki3pointm.csv"), append=FALSE, row.names = FALSE, col.names = TRUE,sep = ",")

    write.table(WEI_N, file = file.path(pats, "Resultats", "WEI_N.csv"), append=FALSE, row.names = FALSE, col.names = TRUE,sep = ",")

    write.table(WEI_S, file = file.path(pats, "Resultats", "WEI_S.csv"), append=FALSE, row.names = FALSE, col.names = TRUE,sep = ",")

    write.table(WEI_SS, file = file.path(pats, "Resultats", "WEI_SS.csv"), append=FALSE, row.names = FALSE, col.names = TRUE,sep = ",")


    ID_N<-unique(WEI_N$ID)
    ID_S<-unique(WEI_S$ID)
    ID_SS<-unique(WEI_SS$ID)


    write.table(ID_N, file = file.path(pats, "Resultats", "ID_N.csv"), append=FALSE, row.names = FALSE, col.names = FALSE,sep = ",")

    write.table(ID_S, file = file.path(pats, "Resultats", "ID_S.csv"), append=FALSE, row.names = FALSE, col.names = FALSE,sep = ",")

    write.table(ID_SS, file = file.path(pats, "Resultats", "ID_SS.csv"), append=FALSE, row.names = FALSE, col.names = FALSE,sep = ",")



    headrep=WEI_SS$rep[1]

    titlerep<-as.character(headrep)


    nfile2=paste(nameOUT, sep="", nb_loci) #Set the name of the export file

    nfile2bis=paste(nfile2, sep="=",nb_indiv,cptt) #Set the name of the import file

    nfile2bis2=paste(nfile2bis, sep="","_BoiteAMoustache.pdf") #Set the name of the import file

    nfile2bis3=paste(nfile2bis, sep="","_Courbe.pdf") #Set the name of the import file

    nfile2bis4=paste(nfile2bis, sep="","_Points.pdf") #Set the name of the import file

    ndata_BAM=file.path(pats, "Resultats", "Courbes", nfile2bis2)

    ndata_CL=file.path(pats, "Resultats", "Courbes", nfile2bis3)

    ndata_PT=file.path(pats, "Resultats", "Courbes", nfile2bis4)


    if (grepl("FrequenceHaplo",nameOUT)) {  #Plot with the scale (0,1) for the frequencies

        ggplot(rikiki3.m,aes(x=gen,y=value,group=gen,fill=as.factor(WEI))) + geom_boxplot() + facet_grid(.~WEI) + scale_fill_manual(breaks = c("0", "9", "18"),labels = c("Neutre", "Semi-codant", "Codant"),values = c(colors4[3], colors4[2], colors4[1])) + ylim(-0.1,1.1) + ylab("Frequency") + xlab("Generations") + theme(axis.text=element_text(size=15,face="bold",family="Times"), axis.title=element_text(size=18,face="bold",family="Times",colour=colors3[7]),panel.background = element_rect(fill = colors3[2], colour = colors3[2]),legend.title = element_blank())


        f0=ggplot(rikiki3.m,aes(x=gen,y=value,group=ID,colour=as.factor(WEI))) + geom_path(size=1) + facet_grid(.~WEI,scales="free") + scale_color_manual(breaks = c("0", "9", "18"),labels = c("Neutre", "Semi-codant", "Codant"),values = c(colors4[3], colors4[2], colors4[1])) + ylim(-0.1,1.1) + ylab("Frequency") + xlab("Generations") + theme(axis.text=element_text(size=15,face="bold",family="Times"), axis.title=element_text(size=18,face="bold",family="Times",colour=colors3[7]),panel.background = element_rect(fill = colors3[2], colour = colors3[2]),legend.title = element_blank())

        ggsave(filename=ndata_CL, device="pdf", width=350, height=250, units="mm", dpi=300, limitsize=FALSE)


        f2=ggplot(rikiki3.m,aes(x=gen,y=value,group=ID,colour=as.factor(WEI))) + geom_point(size=2) + facet_grid(.~WEI) + scale_color_manual(breaks = c("0", "9", "18"),labels = c("Neutre", "Semi-codant", "Codant"),values = c(colors4[3], colors4[2], colors4[1])) + ylim(-0.1,1.1) + ylab("Frequency") + xlab("Generations") + theme(axis.text=element_text(size=15,face="bold",family="Times"), axis.title=element_text(size=18,face="bold",family="Times",colour=colors3[7]),panel.background = element_rect(fill = colors3[2], colour = colors3[2]),legend.title = element_blank())

        ggsave(filename=ndata_PT, device="pdf", width=350, height=250, units="mm", dpi=300, limitsize=FALSE)
    }

    cptt=cptt+1
}


#####################################################################################################################





print("+++++++++++++++------------------------------ | Pairwiseplot_Script.r    | --- | FAIT | ------------------++++++++++++++++++++")





#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================