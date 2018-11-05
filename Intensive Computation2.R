########
#Coded by Cristina Robinson
#Last Modified 2-2-2018
#Written in R-Studio Version 1.1.419
#R Version 3.4.3
#phytools v0.6-30     ape v4.1    maps v3.2.0
#R.utils_2.6.0  
########


rm(list = objects())
#setwd("C:/Users/Kara/Documents/R/2018-OCEvolution")##set your directory here
#setwd("C:/Users/Kara/Desktop/DataWarehouse")
#setwd("C:/Users/Kara/Desktop")
dir <- getwd()
#Libraries and scripts
source("Main functions_3")

#load trees and data (this takes a while)
ConsensusTree <- LoadPrettyTree("Hack_output.nex") #note that the variable Path is also generated here
birbs <- Loadbirds("AllData 5-24-18.csv", ConsensusTree)


#Set variables and filestructure to output data
#Thermo variables
OpenClose <- birbs$O.C
OpenClose <- factor(OpenClose, labels = c("closed","open"))
names(OpenClose) <- birbs$BirdtreeFormat
#Loop variables for main analysis:
call <- names(birbs)[1:(length(birbs))]
call <- call[!call %in% c("O.C", "BirdtreeFormat", "Family", "FormerFamily")]#remove non-independant vars
mod <- rep("log", length(call))
mod[1] <- "linear"
FullRate <- RateMatx(ace(x=OpenClose, phy=ConsensusTree, type="discrete", model="ER"))

#Loop variables for the Jacknife; based on significnace from the main analysis
mod2 <- mod[c(2, 5, 8, 11, 14, 18)]
call2 <- call[c(2, 5, 8, 11, 14, 18)]
#Create folder Structure
MakeFolderStructure()


 
####FIGURES:




#1) make a histogram and boxplot
QuickScatterBox(vari=cbind(birbs$Syllable.rep.final,birbs$Song.rep.final),
                OC=birbs$O.C,title=c("Syllables","Songs"))

#2) MAIN DATA:
#creates RainbowPlots, runs phylANOVA, outputs text, runs Brownielite, plots it

#Figure for paper
pdf("DoublePlot.pdf")
par(mfrow=c(1,2), mar=c(.1,.1,.1,.1))
DataExtraction(OpenClose, birbs, ConsensusTree, call[2], mod=mod[2], fullrate=FullRate,
               BROWNIE = FALSE, ANOVA = FALSE,DP=TRUE)
DataExtraction(OpenClose, birbs, ConsensusTree, call[8], mod=mod[8], fullrate=FullRate,
               BROWNIE = FALSE, ANOVA = FALSE,DP=TRUE)
dev.off()


ANOVAData <- as.list(1:length(call))
sink(file = "Brownie.txt", append = TRUE, split = FALSE)
for(i in rev(seq_along(call))){
#get data and rainbow plots
  DataExtraction(OpenClose, birbs, ConsensusTree, call[i], mod=mod[i], fullrate=FullRate)
  ANOVAData[[i]] <- ANOVARun
#plotting brownie data
  dataset <- read.csv(paste0(call[i],".csv"))
  datasetFULL <- read.csv(paste0(call[i],"FULL.csv"))
  pdf(paste0(call[i], ".Brownie.pdf"))
  par(mfrow=c(2,2), mgp=c(1.5,.5,0), mar=c(3,3,2,1))
  BrowniePlotRates(dataset, paste0(call[i]), Group=c("Closed-Rate", "Open-Rate"))
  BrowniePlotRates(datasetFULL, paste0(call[i],"FULL"))
  dev.off()
}
sink(file=NULL)

#9 in one brownie :)
pdf("BrownieFullRatesAlltypes.pdf")
par(mfrow=c(3,3), mgp=c(1.5,.5,0), mar=c(3,3,2,1))
for(i in c(1,2,5,8,11,14,17,18)){#FUll rates
  datasetFULL <- read.csv(paste0(call[i],"FULL.csv"))
  title <- unlist(strsplit(call[i],"[.]"))
  title <- paste(toupper(substring(title, 1,1)), substring(title, 2),
        sep="", collapse=" ")
  title <- gsub("Final", "", title)
  BrowniePlotRates(datasetFULL, title)
}
dev.off()
#get Anova results
ANOVAResults(ANOVAData)


#for(i in seq_along(call)){
#  dataset <- read.csv(paste0(call[i],".csv"))
#  datasetFULL <- read.csv(paste0(call[i],"FULL.csv"))
#  BrownieSigTest(dataset, call[i])
#  BrownieSigTest(datasetFULL, call[i])
#}

#3)Jackknife runs
#Because Acrocephalidae is paraphyletic, we merged the Lucustellidae with Acrocephalidae
birbs$Family[which(birbs$BirdtreeFormat == "Locustella_naevia")] <- "Acrocephalidae"
#create a list of indicies belonging to each familiy, get indicied and figure out which have 4+ species
Families <- replicate(length(levels(birbs$Family)),NULL)
names(Families) <- levels(birbs$Family)
for(i in 1:length(Families)){Families[[i]]<-which(birbs$Family == names(Families)[i])}
Remove <- names(which(sapply(Families,length)>=4))
Type <- c("", "FULL")
#first loop (i) enters folder for song variable and sets up ANOVA data
#second loop (j) cuts out each of the families in turn and runs the dataextraction protocol
#Third loop (k) run loop 4 with full and partial rates
#Fourth loop (l) generate and plot brownie data
#loop 1
for(i in seq_along(call2)){
  setwd(file.path(dir, "DataWarehouse/Jackknife",call2[i]))
  ANOVAData <- as.list(1:length(Remove))
#loop 2: Jackknife using the ACE values from the tree created after species with NAs for
#a song variable were removed and those which were removed by the jacknife procedure itself
  for(j in 1:length(Remove)){
    ConseJack <- drop.tip(ConsensusTree, Families[[Remove[j]]], root.edge = 0)
    Jacks <- birbs[-Families[[Remove[j]]],]
    OC <- OpenClose[-Families[[Remove[j]]]]
    DataExtraction(OC, Jacks, ConseJack, vari=call2[i], 
                   mod=mod2[i], cotitle=paste0("No", Remove[j]), fullrate=FullRate)
    ANOVAData[[j]] <- ANOVARun
  }
#loop 3: runs full and partial rates
#for(i in seq_along(call2)){
  #setwd(file.path(dir, "DataWarehouse/Jackknife",call2[i]))
  sink(file = "Brownie.txt", append = TRUE, split = FALSE)
  for(k in 1:2){
    pdf(paste0(call2[i], Type[k],"Jackknife.pdf"))
    par(mfrow=c(3,3), mgp=c(1.5,.5,0), mar=c(3,3,2,1))
#Plot the original brownie run
    setwd(file.path(dir, "DataWarehouse"))
    dataset <- read.csv(paste0(call2[i],Type[k],".csv"))
    BrowniePlotRates(dataset,paste0(call2[i], " All"))
    setwd(file.path(dir, "DataWarehouse/Jackknife",call2[i]))
#loop 4 Brownie plot the jacknife Runs   
    for(l in 1:length(Remove)){
      dataset <- read.csv(paste0(call2[i],"No",Remove[l],Type[k],".csv"))
      JackLoss <- length(which(is.na(birbs[,call2[i]][Families[[Remove[l]]]])==FALSE))
      BrowniePlotRates(dataset,paste0(call2[i]," No ",Remove[l], "(", JackLoss, ")"))  
    }
    dev.off()
  }
  sink(file = NULL)
  ANOVAPrinter(ANOVAData, CritAlpha[which(call == call2[i])-1])
}




#4) Jacknife with individual Mimids:
#based on the data from above, we decided to repeat the Brownie Analysis
#with each Mimid removed in turn for syll.song
setwd(file.path(dir, "DataWarehouse/MimidJackknife"))
speciesIndex <- Families$Mimidae
NAind <- which(is.na(birbs$Syll.song.final)==TRUE)
SylSong <- birbs$Syll.song.final
names(SylSong) <- birbs$BirdtreeFormat
for(i in 1:length(speciesIndex)){
  ConseMime <- drop.tip(ConsensusTree, c(speciesIndex[i],NAind), root.edge = 0)
  Mime <- birbs[-c(speciesIndex[i], NAind),]
  OCmime <- OpenClose[-c(speciesIndex[i],NAind)]
  sySo <- log(SylSong[-c(speciesIndex[i],NAind)])
  BrownieDataGen(ConseMime, OCmime, sySo, nsim=1300,title=paste(birbs$BirdtreeFormat[speciesIndex[i]]), FullRate)
}

pdf("MimidJackkinfe.pdf")
par(mfrow=c(2,2))
sink(file = "Brownie.txt", append = TRUE, split = FALSE)
for(i in 1:length(speciesIndex)){
  dataset <- read.csv(paste0(birbs$BirdtreeFormat[speciesIndex[i]], ".csv"))
  BrowniePlotRates(dataset, paste("Syl.Song","No",birbs$BirdtreeFormat[speciesIndex[i]], sep=" "))
}
sink(file = NULL)
dev.off()

 



#5) Species Brownie
#Make the empty variables
Fams <- length(Remove)
Families2 <- replicate(length(levels(birbs$Family)),NULL)
names(Families2) <- levels(birbs$Family)
MultiFamSim <- replicate(length(Remove),NULL)
TextOutput <- vector(mode="character",5)

#go where you belong
setwd(file.path(dir, "DataWarehouse/SpecRates"))
for(i in seq_along(call2)){

  #set up place to store brown data and var to test
  Title <- paste0(call2[i], "BrownieFamilyRates")
  
  #setting the variables and dropping NAs
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
  for(l in rev(1:Fams)){
#Do not run simmap if all species in a family have been removed
    numSpec <- length(which(birbs$BirdtreeFormat[Families[[Remove[l]]]] %in% names(noData)))
    if(numSpec==length(Families[[Remove[l]]])){
      MultiFamSim[[l]] <- NULL
      warning(paste0("No birds for ", call2[i], " in ", Remove[[l]]))
    }
  }
  BrownieWithWarnings(MultiFamSim, TestVar, Title)
  BrownieData <- read.csv(paste0(Title, ".csv"))
  pdf(paste0(Title, ".pdf"))
  par(mfrow=c(3,3))
  for(m in seq_along(MultiFamSim)){
    plotSimmap(MultiFamSim[[m]], fsize = .6, ftype = "i", lwd=1.5,
               hold=FALSE)
    TextOutput[1] <- paste0("ER=", round(BrownieData$ERRate[m], digits=4))
    TextOutput[2] <- paste0("ARD0=", round(BrownieData$ARDRate0[m], digits=4))
    TextOutput[3] <- paste0("ARD1=", round(BrownieData$ARDRate1[m], digits=4))
    TextOutput[4] <- paste0("pVal=", round(BrownieData$Pval[m], digits=4))
    if(BrownieData$convergence[m] == "Optimization has converged."){
      TextOutput[5] <- ""
    }else{TextOutput[5] <- "May not have converged"}
    text(2,seq(10,4,length=5),TextOutput, cex=.6, adj=0)
  }
  dev.off()
}




#6) Test with Lincolnii Closed
setwd(file.path(dir, "DataWarehouse/ClosedLink"))
OClink <- OpenClose
OClink[which(birbs$BirdtreeFormat == "Melospiza_lincolnii")] <- "closed"
#creates RainbowPlots, runs phylANOVA, outputs text, runs Brownielite, plots it
ANOVAData <- as.list(1:length(call))
sink(file = "Brownie.txt", append = TRUE, split = FALSE)
for(i in rev(seq_along(call))){
  DataExtraction(OClink, birbs, ConsensusTree, call[i], mod=mod[i], fullrate=FullRate)
  ANOVAData[[i]] <- ANOVARun
  #plotting data
  dataset <- read.csv(paste0(call[i],".csv"))
  datasetFULL <- read.csv(paste0(call[i],"FULL.csv"))
  pdf(paste0(call[i], ".Brownie.pdf"))
  par(mfrow=c(2,2))
  BrowniePlotRates(dataset, paste0(call[i]))
  BrowniePlotRates(datasetFULL ,paste0(call[i],"FULL"))
  dev.off()
}
sink(file=NULL)
ANOVAResults(ANOVAData)
  