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
source("Main functions_3.R")

#load trees and data (this takes a while)
ConsensusTree <- LoadPrettyTree("Hack_output.nex") #note that the variable Path is also generated here
birbs <- Loadbirds("AllData 2-2-18.csv", ConsensusTree)


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
#9 in one :)
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
  