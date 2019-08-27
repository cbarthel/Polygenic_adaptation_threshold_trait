# -*- coding: utf-8 -*-


# Script for manuscript "Detecting Polygenic Adaptation For A Threshold Trait In Experimental Populations Using Time-Series Data"


#Clément BARTHÉLÉMY, La Rochelle University, Station Ifremer de La Tremblade


#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================


#====================== Path ======================#


pats=file.path("/home", "clement", "Documents", "SimulationPP")


#====================== Packages ======================#


library(ggplot2)
library(RColorBrewer)
library(reshape2)


#====================== Import data ======================#


MatC<-read.csv(file.path(pats, "MCC_Results_Intra.txt"),sep=",", header=FALSE)
scenario<-read.csv(file.path(pats, "Script", "sc.txt"),sep=" ", header=FALSE)


#====================== Figure Supplementary S6 ======================#


accuracy=vector()

for (i in seq(1,80,1)){
  
  accuracy[i]=(MatC[i,1]+MatC[i,2])/(MatC[i,1]+(MatC[i,2]+MatC[i,3]+MatC[i,4]))
}

F1score=vector()

for (i in seq(1,80,1)){
  
  F1score[i]=(2*MatC[i,1])/((2*MatC[i,1])+MatC[i,3]+MatC[i,4])
}


DataSenSpe=cbind(accuracy,F1score)

DataMelt=melt(DataSenSpe)

ggplot() + geom_boxplot(data=DataMelt, aes(factor(DataMelt$Var2,levels=c("accuracy","F1score")),value,fill=Var2), size=1.5) + xlab("Classification functions") + ylab("Value") + theme_bw() + theme(axis.text=element_text(size=15,face="bold",family="Times",colour="black"), axis.title=element_text(size=18,face="bold",family="Times"), legend.position="none")


#====================== Figure 8 ======================#


#Output distributions


MCC_Results=MatC[1:80,]
scenario=scenario[1:80,]


x=seq(1,80,1)
y=sort(MCC_Results$V5)
data_Inter=MCC_Results
data_Inter=cbind(data_Inter,x)

data_tot=cbind(scenario,data_Inter)
colnames(data_tot)=c("x1","x2","x3","x4","x5","x6","x7","TP","TN","FP","FN", "MCC","rep", "scenar")
data_tot=as.data.frame(data_tot)


Freq_TP=vector()
Freq_TN=vector()
Freq_FP=vector()
Freq_FN=vector()


for (i in seq(1,80,1)){

    nb_pair=(data_tot[i,4]*(data_tot[i,4]-1))/2

    nb_neutre=((data_tot[i,4]-10)*(data_tot[i,4]-11))/2
    nb_codant=nb_pair-nb_neutre

    Freq_TP[i]=(data_tot[i,]$TP)/nb_codant
    Freq_TN[i]=(data_tot[i,]$TN)/nb_neutre
    Freq_FP[i]=(data_tot[i,]$FP)/nb_neutre
    Freq_FN[i]=(data_tot[i,]$FN)/nb_codant

}


data_tot=cbind(data_tot,Freq_TP,Freq_TN,Freq_FP,Freq_FN)

Freq_MCC=cbind(data_tot$MCC,Freq_TP,Freq_TN,Freq_FP,Freq_FN)
colnames(Freq_MCC)=c("MCC","Sensitivity","Specificity","FPR","FNR")
Freq_MCC=as.data.frame(Freq_MCC)


sd(Freq_MCC$MCC)


Freq_MCC2=as.data.frame(melt(Freq_MCC))
Freq_MCC2$variable <- as.character(Freq_MCC2$variable)

Freq_MCC3=cbind(Freq_MCC,x)
Freq_MCC3=Freq_MCC3[order(Freq_MCC3$MCC),]
Freq_MCC3$x=factor(Freq_MCC3$x, levels = Freq_MCC3$x[order(Freq_MCC3$MCC)])


# Plot

colors4 <- brewer.pal(9, "Set1")

ggplot() + geom_boxplot(data=Freq_MCC2, aes(factor(Freq_MCC2$variable,levels=c("Sensitivity","Specificity","FPR","FNR","MCC")),value,fill=variable), size=1.5) + scale_fill_manual(values=c(colors4[4],colors4[2],colors4[5],colors4[3],colors4[1])) + xlab("Output variables") + ylab("Frequencies") + scale_y_continuous(sec.axis = sec_axis(~.+0, name = "MCC")) + theme_bw() + theme(axis.text=element_text(size=15,face="bold",family="Times",colour="black"), axis.title=element_text(size=18,face="bold",family="Times"), legend.position="none")


#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================