########
#Coded by Cristina Robinson and Kate Synder
#Last Modified 11-8-2017
#Written in R-Studio Version 1.1.383
#R Version 3.4.2
#phytools v0.6-30     ape v4.1    maps v3.2.0
#R.utils_2.6.0  
########

rm(list = objects())
#install.packages("maps") 
#install.packages("ape")
#install.packages("phytools")
#install.packages("R.utils")

#Libraries
library(maps)
library(ape)
library(phytools)
library(R.utils)

#setwd("C:/Users/Kara/Documents/R")##set your directory here
dir <- getwd()
source("Main functions_2B.R")

#load trees and data
ConsensusTree <- LoadPrettyTree("cristinatree.nex") #note that the variable Path is also generated here
birbs <- Loadbirds("AllData 12-7.csv", ConsensusTree)


setwd(file.path(dir, "DataWarehouse"))


#set up thermo variable once to save time
OpenClose <- birbs$O.C
OpenClose <- factor(OpenClose, labels = c("closed","open"))
names(OpenClose) <- birbs$BirdtreeFormat

#loop variables:
call <- names(birbs)[1:(length(birbs))]
call <- call[!call %in% c("O.C", "BirdtreeFormat", "Family", "FormerFamily")]
mod <- rep("log", length(call))
mod[1:2] <- "linear"
FullRate <- RateMatx(ace(x=OpenClose, phy=ConsensusTree, type="discrete", model="ER"))





###Jackknife runs
#Because Acrocephalidae is paraphyletic, we merged the Lucustellidae with Acrocephalidae
birbs$Family[9] <- "Acrocephalidae"

#create a list of indicies belonging to each familiy
Families <- replicate(length(levels(birbs$Family)),NULL)
names(Families) <- levels(birbs$Family)
#get the indicies for species in each family, figure out which families have 4+ species
for(i in 1:length(Families)){Families[[i]]<-which(birbs$Family == names(Families)[i])}
Remove <- names(which(sapply(Families,length)>=4))

#New Loop variables on restricticted subset of significnat variables to jacknife
mod2 <- mod[c(4, 7, 10, 13, 16, 19)]
call2 <- call[c(4, 7, 10, 13, 16, 19)]




#Species Brownie
#Make the empty variables
Fams <- length(Remove)
Families2 <- replicate(length(levels(birbs$Family)),NULL)
names(Families2) <- levels(birbs$Family)
browniedata <- data.frame(Familiy=character(Fams),Pval=numeric(Fams), ERRate=numeric(Fams), ERloglik=numeric(Fams),
                          ERace=numeric(Fams), ARDRate0=numeric(Fams), ARDRate1=numeric(Fams),
                          ARDloglik=numeric(Fams),ARDace=numeric(Fams), convergence=character(Fams),
                          stringsAsFactors = FALSE)
MultiFamSim <- replicate(length(Remove),NULL)
#go where you belong
setwd(file.path(dir, "DataWarehouse/SpecRates"))

for(i in 1:length(call2)){
  #set up place to store brown data and var to test
  Title <- paste0(call2[i], "BrownieFamilyRates")
  
  #setting the variable and dropping NAs
  TestVar <- log(unlist(birbs[call2[i]]))
  names(TestVar) <- birbs$BirdtreeFormat
  noData <- which(is.na(TestVar) == TRUE)
  if(length(noData) > 0){
    TestVar <- TestVar[-noData]
    conse <- drop.tip(ConsensusTree, noData, root.edge = 0)
    dropfam <- birbs$Family[-noData]
  }else{conse <- ConsensusTree
  dropfam <- birbs$Family}
  
  #Create the simmap object for each family
  newPath <- PathFinder2(conse)
  for(j in 1:length(Families2)){Families2[[j]]<-which(dropfam == names(Families2)[j])}
  for(k in 1:length(Remove)){MultiFamSim[[k]] <- SpeciesSimmap(conse, Families2[[Remove[k]]], newPath)}
  names(MultiFamSim) <- Remove
  
  #Run brownie for each family
  for(l in 1:Fams){
    #Do not run if all species in a family have been removed
    numSpec <- length(which(birbs$BirdtreeFormat[Families[[Remove[l]]]] %in% names(noData)))
    if(numSpec==length(Families[[Remove[l]]])){
      break(warning(paste0("No birds for ", call2[i], " in ", Remove[[l]])))
    }
    Results <- brownie.lite(MultiFamSim[[l]], TestVar, maxit=750000)
    browniedata[l,2:9] <- c(Results[[11]], Results[[1]], Results[[4]],  Results[[2]],
                            Results[[6]][1], Results[[6]][2], Results[[9]],
                            Results[[7]])
    browniedata[l,c(1,10)] <-c(Remove[l],Results[[12]])
  }
  write.csv(file = paste(Title,".csv",sep=""), browniedata)
  pdf(paste0(call2[i], "SpecRates.pdf"))
  par(mfrow=c(3,3))
  for(m in 1:Fams){
    plotSimmap(MultiFamSim[[m]], fsize=.6, lwd=1.5, setEnv = FALSE,
               hold=FALSE)
    text(-2,length(TestVar),Remove[m], cex=.7, font=1,adj=0)
    pval <- round(browniedata$Pval[m], digits = 3)
    if(pval == 0){pval <- ">.001"}
    text(-2,length(TestVar)-2,paste0("pVal = ", pval), cex=.7, adj=0)
  }
  dev.off()  
}  