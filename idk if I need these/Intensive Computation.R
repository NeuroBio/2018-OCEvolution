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

setwd("C:/Users/Kara/Documents/R")##set your directory here
dir <- getwd()
source("Main functions_2B.R")

#load trees and data
ConsensusTree <- LoadPrettyTree("cristinatree.nex") #note that the variable Path is also generated here
birbs <- Loadbirds("AllData 12-7.csv", ConsensusTree)

#create output directories and ANOVA printouts
dir.create(file.path(dir, "DataWarehouse"))
dir.create(file.path(dir, "DataWarehouse/Jackknife"))
dir.create(file.path(dir, "DataWarehouse/MimidJackknife"))
dir.create(file.path(dir, "DataWarehouse/SpecRates"))
setwd(file.path(dir, "DataWarehouse"))
file.create("ANOVA.txt")

#make a histogram and boxplot
QuickHistBox(birbs$Syllable.rep.final,birbs$O.C)

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

#Main data (creates RainbowPlots, runs phylANOVA and Brownielite, and plots it)
for(i in 1:length(call)){
  DataExtraction(OpenClose, birbs, ConsensusTree, call[i], mod=mod[i], fullrate=FullRate)
  dataset <- read.csv(paste0(call[i],".csv"))
  datasetFULL <- read.csv(paste0(call[i],"FULL.csv"))
  pdf(paste0(call[i], ".Brownie.pdf"))
  par(mfrow=c(2,2))
  BrowniePlot(dataset,paste0(call[i]))
  BrowniePlot(datasetFULL,paste0(call[i],"FULL"))
  dev.off()
}




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

#first loop (i) pulls out a song variable, makes a folder for it and enters that folder
#second loop (j) cuts out each of the families in turn and runs the dataextraction protocol
#Third loop (k) generate and plot brownie data using partial tree ACE rates
#forth loop (k) generate and plot brownie data using full tree ACE rates

#loop 1
for(i in 1:length(call2)){
  setwd(file.path(dir, "DataWarehouse/Jackknife"))
  dir.create(file.path(dir, "DataWarehouse/Jackknife",call2[i]))
  dir.create(file.path(dir, "DataWarehouse/MimidJackknife"))
  dir.create(file.path(dir, "DataWarehouse/SpeciesBrownie"))
  setwd(file.path(dir, "DataWarehouse/Jackknife",call2[i]))
  file.create("ANOVA.txt")
  
#loop 2: Jackknife using the ACE values from the tree created after species with NAs for
#a song variable were removed and those which were removed by the jacknife procedure itself
  for(j in 1:length(Remove)){
    ConseJack <- drop.tip(ConsensusTree, Families[[Remove[j]]], root.edge = 0)
    Jacks <- birbs[-Families[[Remove[j]]],]
    OC <- OpenClose[-Families[[Remove[j]]]]
    DataExtraction(OC, Jacks, ConseJack, vari=call2[i], 
                   mod=mod2[i], cotitle=paste0("No", Remove[j]), fullrate=FullRate)
  }

  pdf(paste0(call2[i], "Jackknife.pdf"))
  par(mfrow=c(3,3))
  setwd(file.path(dir, "/DataWarehouse"))
  dataset <- read.csv(paste0(call2[i],".csv"))
  BrowniePlot(dataset,paste0(call2[i], " All"))
  setwd(file.path(dir, "DataWarehouse/Jackknife",call2[i]))

#loop 3 Brownie plotting   
  for(k in 1:length(Remove)){
    dataset <- read.csv(paste0(call2[i],"No",Remove[k],".csv"))
    JackLoss <- length(which(is.na(birbs[,call2[i]][Families[[Remove[k]]]])==FALSE))
    BrowniePlot(dataset,paste0(call2[i]," No ",Remove[k]),JackLoss)  
  }
  dev.off()

#Jackknife using the Ace values from the whole tree  (only Brownie)
  pdf(paste0(call2[i], "JackknifeFull.pdf"))
  par(mfrow=c(3,3))
  setwd(file.path(dir, "/DataWarehouse"))
  dataset <- read.csv(paste0(call2[i],"FULL",".csv"))
  BrowniePlot(dataset,paste0(call2[i], " All"))
  setwd(file.path(dir, "DataWarehouse/Jackknife",call2[i]))

#loop 4 Brownie plotting  
  for(l in 1:length(Remove)){
    dataset <- read.csv(paste0(call2[i],"No",Remove[l],"FULL",".csv"))
    JackLoss <- length(which(is.na(birbs[,call2[i]][Families[[Remove[l]]]])==FALSE))
    BrowniePlot(dataset,paste0(call2[i]," No ",Remove[l]),JackLoss)  
  }
  dev.off()
}




###Jacknife with individual Mimids:
#based on the data from above, we decided to repeat the Brownie Analysis with each Mimid removed in turn
setwd(file.path(dir, "DataWarehouse/MimidJackknife"))
speciesIndex <- Families$Mimidae
NAind <- which(is.na(birbs$Syll.song.final)==TRUE)
SylSong <- birbs$Syll.song.final
names(SylSong) <- birbs$BirdtreeFormat
for(i in 1:length(speciesIndex)){
  ConseMime <- drop.tip(ConsensusTree, c(speciesIndex[i],NAind), root.edge = 0)
  Mime <- birbs[-c(speciesIndex[i], NAind),]
  OC <- OpenClose[-c(speciesIndex[i],NAind)]
  sySo <- log(SylSong[-c(speciesIndex[i],NAind)])
  
  BrownieDataGen(ConseMime, OC, sySo, nsim=1300,title=paste(birbs$BirdtreeFormat[speciesIndex[i]]), FullRate)
}
pdf("MimidJackkinfe.pdf")
par(mfrow=c(2,2))
for(i in 1:length(speciesIndex)){
dataset <- read.csv(paste0(birbs$BirdtreeFormat[speciesIndex[i]], ".csv"))
BrowniePlot(dataset, "Syl.Song", paste("No",birbs$BirdtreeFormat[speciesIndex[i]], sep=" "))
}
dev.off()
 



i <- 3 
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
  }else{conse <- ConsensusTree}
  
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
}  
  




  












































  #pdf("Simtest.pdf")
#par(mfrow=c(2,2))
MultiFamSim <- replicate(length(Remove),NULL)
for(i in 1:length(Remove)){MultiFamSim[[i]] <- SpeciesSimmap(ConsensusTree, Families[[Remove[i]]], Path)}
names(MultiFamSim) <- Remove
Fams <- length(Remove)
#dev.off()
#plot(MultiFamSim[[1]])
#vari <- birbs$Syllable.rep.final
#names(vari) <- birbs$BirdtreeFormat
i <- 2
j <- 8
for(i in 1:length(call2)){
#set u place to store brown data and var to test
  Title <- paste0(call2[i], "BrownieFamilyRates")
  browniedata <- data.frame(Familiy=character(Fams),Pval=numeric(Fams), ERRate=numeric(Fams), ERloglik=numeric(Fams),
                            ERace=numeric(Fams), ARDRate0=numeric(Fams), ARDRate1=numeric(Fams),
                            ARDloglik=numeric(Fams),ARDace=numeric(Fams), convergence=character(Fams),
                            stringsAsFactors = FALSE)
  
  TestVar <- log(unlist(birbs[call2[i]]))
  names(TestVar) <- birbs$BirdtreeFormat
  noData <- which(is.na(TestVar) == TRUE)
  if(length(noData) > 0){
    TestVar <- TestVar[-noData]
  }
    
  for(j in 1:Fams){
    print(Remove[j])
#remove any NAs
    if(length(noData) > 0){
      RunBrownie <-  drop.tip(MultiFamSim[[j]], noData, root.edge = 0)
      RunBrownie[[5]] <- RunBrownie[[5]][-]
    }else(RunBrownie <- MultiFamSim[[j]])
    numSpec <- length(which(birbs$BirdtreeFormat[Families[[Remove[j]]]] %in% names(noData)))
    print(numSpec)
    if(numSpec==length(Families[[Remove[j]]])){
      break(warning(paste0("No birds for ", call2[i], " in ", Remove[[j]])))
    }
    Results <- brownie.lite(RunBrownie, TestVar, maxit=750000)
    browniedata[j,2:9] <- c(Results[[11]], Results[[1]], Results[[4]],  Results[[2]],
                        Results[[6]][1], Results[[6]][2], Results[[9]],
                        Results[[7]])
    browniedata[j,c(1,10)] <-c(Remove[j],Results[[12]])
  }
  write.csv(file = paste(Title,".csv",sep=""), browniedata)
}
j <- 6
plot(RunBrownie)
#mapped.edge
puresim[[2]]$mapped.edge
edges <-  paste(ConsensusTree$edge[,1],ConsensusTree$edge[,2], sep=",")
OneIndex <- c(1,2)
MultiFamSim[[1]]$maps
#test for rates in each family
test <- matrix(data=0,nrow=length(birbs$Family),
               ncol=length(levels(birbs$Family)))
rownames(test) <- ConsensusTree$tip.label
for(i in 1:length(levels(birbs$Family))){
  test[birbs$Family == levels(birbs$Family)[i],i] <- 1 
}


MultiFamSim[[1]][[5]]


puresim <- replicate(length(levels(birbs$Family)),NULL)
for(i in 1:ncol(test)){
  puresim[[i]] <- make.simmap(tree=ConsensusTree, x=test[,i], nsim=1, message=FALSE)
}

Conse <- ConsensusTree
#for(i in 1:ncol(test)){
#plot(puresim[[i]])
#}
#rm(SimPaths)


#quick edit graphs!!!
#for(i in 1:length(call)){
#  dataset<- as.data.frame(read.csv(paste0(call[i],".csv")))
#  BrowniePlot(dataset, title=call[i])
#  dataset2<- as.data.frame(read.csv(paste0(call[i],".csv")))
#  BrowniePlot(dataset, title=paste0(call[i], "FULL"))
#}

#lol <-phylANOVA(ConsensusTree,OpenClose,birbs$Syllable.rep.final)
#er <- ace(x=OpenClose, phy=ConsensusTree, type="discrete", model="ER")
#ard <- ace(x=OpenClose, phy=ConsensusTree, type="discrete", model="ARD")
#er

#anova(er, ard)
#qq <- matrix(c(-.0202,.0202, .0202,-.0202),2,2)
#rownames(qq) <- c("closed","open")
#colnames(qq) <- c("closed","open")

#multisimmap <- make.simmap(tree=ConsensusTree,x=OpenClose, nsim=20, message=TRUE,  Q=qq ) 
#plot(multisimmap)

#par(mfrow=c(1,2))
#testvarpretty <- birbs$Syllable.rep.final
#boxplot(testvarpretty)
#quantile(testvarpretty)
#hist(testvarpretty)
#boxplot(testvarpretty)
#quant <- quantile(testvarpretty)
#mean(c(a,b))


#hist(testvarpretty, breaks = seq(0, 10, by=.5))
#range(testvarpretty)
#mean(testvarpretty)
#testvarpretty <- log(birbs$Syllable.rep.final)
#testvarpretty[which(testvarpretty > 4.5)] <- 4.5
#testvarpretty[which(testvarpretty < .5)] <- .5
#rang <- (range(testvarpretty)[2] - range(testvarpretty)[1])
#rang/mean(testvarpretty)
#low <- mean(testvarpretty)-range(testvarpretty)[1]
#high <- range(testvarpretty)[2]-mean(testvarpretty)
#high-low
#shift <- mean(testvarpretty)+(high-low)
#testvarpretty[which(testvarpretty > shift)] <- shift
#names(testvarpretty) <- names(OpenClose)


#QuickHHistBox(birbs$Syllable.rep.final,birbs$O.C)


