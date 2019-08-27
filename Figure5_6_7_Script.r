# -*- coding: utf-8 -*-


# Script for manuscript "Detecting Polygenic Adaptation For A Threshold Trait In Experimental Populations Using Time-Series Data"


#Clément BARTHÉLÉMY, La Rochelle University, Station Ifremer de La Tremblade


#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================


#====================== Packages ======================#


library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(gridExtra)
library(cowplot)
library(psych)


#====================== Path ======================#


pats=file.path("/home", "clement", "Documents", "SimulationPP")


#====================== Import data ======================#


obs1<-read.csv(file.path(pats, "Resultats_Scenario_XXXX", "myRUN_XX_XX.csv"),header=T,stringsAsFactors = FALSE)

obs2<-read.csv(file.path(pats, "Resultats_Scenario_XXXX", "myRUN_XX_XX.csv"),header=T,stringsAsFactors = FALSE)

obs3<-read.csv(file.path(pats, "Resultats_Scenario_XXXX", "myRUN_XX_XX.csv"),header=T,stringsAsFactors = FALSE)

obs4<-read.csv(file.path(pats, "Resultats_Scenario_XXXX", "myRUN_XX_XX.csv"),header=T,stringsAsFactors = FALSE)

obs5<-read.csv(file.path(pats, "Resultats_Scenario_XXXX", "myRUN_XX_XX.csv"),header=T,stringsAsFactors = FALSE)

obs6<-read.csv(file.path(pats, "Resultats_Scenario_XXXX", "myRUN_XX_XX.csv"),header=T,stringsAsFactors = FALSE)


#====================== Functions ======================#


Graph_mean <- function(data,nb_gen,repi){

    rikiki<-subset(data[,1:4]) #Cutting the data (rep, gen, z, sum_g)

    rikiki<-rikiki[order(rikiki$gen,rikiki$rep),]


    stat6<-describeBy(rikiki,list(rikiki$gen,rikiki$rep),mat=TRUE) #Use the "psych" package to calculate the mean, median, etc.

    from=(nb_gen*repi)*3+1 #Set the first line of the data that interest us (mean)

    to=from+nb_gen*repi-1 #Set the last line of the data that interest us (mean)

    obs6=stat6[from:to,c(2,3,5,6,7,8,11,12)] #Subset the result

    obs6[,1]<-as.integer(levels(obs6[,1]))[obs6[,1]]
    obs6[,2]<-as.integer(levels(obs6[,2]))[obs6[,2]]



    colorsS <- brewer.pal(6, "Dark2")
    colors3 <- brewer.pal(9, "Greys")
    colors4 <- brewer.pal(9, "BuPu")


    f6 = ggplot() + geom_path(data=obs6, aes(obs6$group1,obs6$mean, group=obs6$group2),colour=colors3[8],size=1) + scale_x_continuous(breaks = seq(1,nb_gen,1)) + ylab("Mean of trait z") + xlab("Generations") + labs(colour = "Repetitions") + theme_bw() + theme(legend.position="none", axis.text=element_text(size=15,face="bold",family="serif",colour="black"), axis.title=element_text(size=18,face="bold",family="serif",colour="black")) # Use the loess data to add the 'ribbon' to plot 

    return(f6)

}


Graph_var <- function(data,nb_gen,repi){

    dat2=data[,c(1,2,3,4,5)]

    dat3=dat2[order(dat2$rep,dat2$gen),]

    rep=max(unique(dat3$rep))
    gen=max(unique(dat3$gen))

    varrep=vector()

    varianc=vector()

    cpt=1

    for (j in unique(dat3$gen)){

        dat4<-dat3[which(dat3$rep==1 & dat3$gen==j),]

        vari=var(dat4$z1)

        varrep[cpt]=vari

        cpt=cpt+1
    }

    repp=rep(1,max(unique(dat4$gen))+1)
    repp=as.data.frame(repp)

    genn=seq(1,max(unique(dat4$gen))+1,1)
    genn=as.data.frame(genn)

    dat5=cbind(repp,genn)

    varrep=as.data.frame(varrep)

    dat6=cbind(dat5,varrep)

    for (i in seq(2,max(unique(dat3$rep)),1)){

        cpt=1
        varrep=vector()

        for (j in unique(dat3$gen)){

            dat4<-dat3[which(dat3$rep==i & dat3$gen==j),]

            vari=var(dat4$z1)

            varrep[cpt]=vari

            cpt=cpt+1
        }

        repp=rep(i,max(unique(dat4$gen))+1)
        repp=as.data.frame(repp)

        dat7=cbind(repp,genn)

        varep=as.data.frame(varrep)

        dat8=cbind(dat7,varrep)

        dat6=rbind(dat6,dat8)
    }



    colorsS <- brewer.pal(6, "Dark2")
    colors3 <- brewer.pal(9, "Greys")


    f6 = ggplot(dat6,aes(x=genn,y=varrep,group=repp)) + geom_path(color=colors3[8],size=1) + scale_x_continuous(breaks = seq(1,nb_gen,1)) + ylab("Variance of trait z") + xlab("Generations") + theme_bw() + theme(legend.position="none", axis.text=element_text(size=15,face="bold",family="serif",colour="black"), axis.title=element_text(size=18,face="bold",family="serif",colour="black"))

    return(f6)

}


Graph_threshold <- function(data6,nb_gen,repi){

    fit_survival=vector()
    fit_surv_data=data.frame(matrix(ncol = 0, nrow = nb_gen))
    list_name=vector()

    cpt=1

    for (i in seq(1,max(data6$rep),1)){

        for (j in seq(0,max(data6$gen),1)){

            data=subset(data6, data6$rep==i & data6$gen==j)
            a=length(data$z1)

            data_fit=subset(data, data$fitness > 0)
            b=length(data_fit$fitness)
            rate=b/a
            fit_survival[j+1]=rate
        }

        pop=paste("Pop", sep="_",cpt)
        list_name[cpt]=pop
        cpt=cpt+1

        fit_surv_data=cbind(fit_surv_data,fit_survival)

    }

    colnames(fit_surv_data)=list_name

    genn=rep.int(1:nb_gen,repi)
    genn=as.data.frame(genn)

    fit_surv_data=melt(fit_surv_data, id.vars = NULL)
    fit_surv_data=cbind(fit_surv_data,genn)

    fit_surv_data6<-fit_surv_data[!(fit_surv_data$genn==1 | fit_surv_data$genn==2 | fit_surv_data$genn==3),]

    colorsS <- brewer.pal(3, "Set1")
    colors3 <- brewer.pal(9, "Greys")

    f6 = ggplot(fit_surv_data6, aes(genn,value,group=variable)) + geom_path(size=0.7, colour=colors3[8]) + scale_x_continuous(breaks = seq(1,nb_gen,1)) + ylim(0,1) + xlab("Generations") + ylab("Survival rate") + theme_bw() + theme(axis.text=element_text(size=15,face="bold",family="serif",colour="black"), axis.title=element_text(size=18,face="bold",family="serif"), legend.position="none")

    return(f6)

}


#====================== Function calls ======================#


m1=Graph_mean(obs1,XX,XX)
m2=Graph_mean(obs2,XX,XX)
m3=Graph_mean(obs3,XX,XX)
m4=Graph_mean(obs4,XX,XX)
m5=Graph_mean(obs5,XX,XX)
m6=Graph_mean(obs6,XX,XX)


v1=Graph_var(obs1,XX,XX)
v2=Graph_var(obs2,XX,XX)
v3=Graph_var(obs3,XX,XX)
v4=Graph_var(obs4,XX,XX)
v5=Graph_var(obs5,XX,XX)
v6=Graph_var(obs6,XX,XX)


t1=Graph_threshold(obs1,XX,XX)
t2=Graph_threshold(obs2,XX,XX)
t3=Graph_threshold(obs3,XX,XX)
t4=Graph_threshold(obs4,XX,XX)
t5=Graph_threshold(obs5,XX,XX)
t6=Graph_threshold(obs6,XX,XX)


#====================== Plots and save ======================#


plot_grid(m1, m2, m3, m4, m5, m6, labels=c("A", "B", "C", "D", "E", "F"), ncol = 2, nrow = 3)
ggsave(filename=file.path(pats, "Figure_5.pdf"), device="pdf", width=400, height=400, units="mm", dpi=600, limitsize=FALSE)


plot_grid(v1, v2, v3, v4, v5, v6, labels=c("A", "B", "C", "D", "E", "F"), ncol = 2, nrow = 3)
ggsave(filename=file.path(pats, "Figure_6.pdf"), device="pdf", width=400, height=400, units="mm", dpi=600, limitsize=FALSE)


plot_grid(t1, t2, t3, t4, t5, t6, labels=c("A", "B", "C", "D", "E", "F"), ncol = 2, nrow = 3)
ggsave(filename=file.path(pats, "Figure_7.pdf"), device="pdf", width=400, height=400, units="mm", dpi=600, limitsize=FALSE)


#====================================================================================================================
######################################  Contact : clement.barthelemy@gmail.com  #####################################
#====================================================================================================================