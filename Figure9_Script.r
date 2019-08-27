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
library(RColorBrewer)
library(reshape2)
library(gridExtra)
library(cowplot)


#====================== Function ======================#


SA_output <- function(input,output){

    tell(input, output)


    mu.star <- apply(input$ee, 2, function(input) mean(abs(input)))
    sigma <- apply(input$ee, 2, sd)

    mu.star=as.data.frame(mu.star)
    sigma=as.data.frame(sigma)

    SA_Morris=cbind(mu.star,sigma)

    name=c("x1","x2","x3","x4","x5","x6","x7")
    label=c("nb_rep","nb_gen","nb_ind","nb_locus","nb_cluster","env","alle_eff")


    SA_Morris_MCC=cbind(SA_Morris,name,label)


    return(SA_Morris_MCC)

}


#====================== Import data ======================#


y<-read.csv(file.path(pats, "MCC_Results_Intra.txt"),sep=",",header=FALSE) #Results of the clustering (MCC values)
sc<-read.csv(file.path(pats, "Script", "sc.txt"),sep=" ",header=FALSE) #Scenarios used for the simulations and generated with the morris function


#====================== Create the Morris scenarios ======================#


set.seed(1311)
Morris <- morris(model = NULL, factors=7, r=10,
                 design = list(type = "oat", levels =6 , grid.jump = 3),
                 binf = c(10, 5, 30, 20, 3, 5, 2),
                 bsup = c(110, 20, 230, 90, 18, 40, 12),
                 scale = TRUE)
design <- Morris$X


sc=as.matrix(sc)
Morris$X <- sc


#====================== Outputs (MCC, Se and Sp) ======================#


TP_rate=vector()
TN_rate=vector()

for (i in seq(1,80,1)){
  x=y[i,]
  x_tot=x[,1]+x[,4]
  TP_rate[i]=x[,1]/x_tot
}

for (i in seq(1,80,1)){
  x=y[i,]
  x_tot=x[,2]+x[,3]
  TN_rate[i]=x[,2]/x_tot
}



y1=TP_rate #Se
y2=TN_rate #Sp
y3=y$V5


SA_Morris_MCC=SA_output(Morris,y1)
SA_Morris_Se=SA_output(Morris,y2)
SA_Morris_Sp=SA_output(Morris,y3)


#====================== Plots ======================#


colorsS <- brewer.pal(7, "Dark2")
colors3 <- brewer.pal(9, "Greys")


MCC = ggplot()  + geom_point(data=SA_Morris_MCC, aes(mu.star, sigma),size=11, colour=colorsS) + ylim(0,0.4) + xlim(0,0.7) + geom_text(data=SA_Morris_MCC, aes(SA_Morris_MCC$mu.star,SA_Morris_MCC$sigma, label=SA_Morris_MCC$label, fontface="bold",family="Times"),hjust=0.4, vjust=-1.5) +  ylab(expression(sigma)) + xlab(mu~"*") + theme_bw() + theme(legend.position="none", axis.text=element_text(size=15,face="bold",family="serif",colour="black"), axis.title=element_text(size=18,face="bold",family="serif",colour="black"))

Se = ggplot()  + geom_point(data=SA_Morris_Se, aes(mu.star, sigma),size=11, colour=colorsS) + ylim(0,0.9) + xlim(0,0.7) + geom_text(data=SA_Morris_Se, aes(SA_Morris_Se$mu.star,SA_Morris_Se$sigma, label=SA_Morris_Se$label, fontface="bold",family="Times"),hjust=0.4, vjust=-1.5) +  ylab(expression(sigma)) + xlab(mu~"*") + theme_bw() + theme(legend.position="none", axis.text=element_text(size=15,face="bold",family="serif",colour="black"), axis.title=element_text(size=18,face="bold",family="serif",colour="black"))

Sp = ggplot()  + geom_point(data=SA_Morris_Sp, aes(mu.star, sigma),size=11, colour=colorsS) + ylim(0,0.9) + xlim(0,0.5) + geom_text(data=SA_Morris_Sp, aes(SA_Morris_Sp$mu.star,SA_Morris_Sp$sigma, label=SA_Morris_Sp$label, fontface="bold",family="Times"),hjust=0.4, vjust=-1.5) +  ylab(expression(sigma)) + xlab(mu~"*") + theme_bw() + theme(legend.position="none", axis.text=element_text(size=15,face="bold",family="serif",colour="black"), axis.title=element_text(size=18,face="bold",family="serif",colour="black"))

plot_grid(Se, Sp, MCC, labels=c("A", "B", "C"), ncol = 2, nrow = 2)


#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================