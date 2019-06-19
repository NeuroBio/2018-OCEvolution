#get the reverse relationship for EPP Poly
OCswitch <- birbs$O.C
names(OCswitch) <- birbs$BirdtreeFormat

Poly <- birbs$Polygyny.final
Poly <- factor(Poly, labels = c("mono","poly"))
names(Poly) <- birbs$BirdtreeFormat
PolyHold <- which(is.na(Poly))
Poly <- Poly[-PolyHold]
PolyTree <- drop.tip(ConsensusTree, PolyHold)
BigTreeRates <- list(rates=c(0.09313762, 0.01119413), lik.anc=matrix(0, ncol=2, nrow=2, dimnames=list(c(1,2),c("mono", "poly"))))
PolyRate <- RateMatx(BigTreeRates)
#anova(ace(x=Poly, phy=PolyTree, type="discrete", model="ER"),ace(x=Poly, phy=PolyTree, type="discrete", model="ARD"))
#PolyRate <- RateMatx(ace(x=Poly, phy=PolyTree, type="discrete", model="ARD"))
BrownieDataGen(PolyTree, Poly, OCswitch[-PolyHold], title="Polygyny on OC", nsim=1300, rater=PolyRate)

EPP <- birbs$EPP.final
EPP <- factor(EPP, labels = c("absent","present"))
names(EPP) <- birbs$BirdtreeFormat
EPPHold <- which(is.na(EPP))
EPP <- EPP[-EPPHold]
EPPTree <- drop.tip(ConsensusTree, EPPHold)
#anova(ace(x=EPP, phy=EPPTree, type="discrete", model="ER"),ace(x=EPP, phy=EPPTree, type="discrete", model="ARD"))
BigTreeRates <- list(rates=c(0.02155044, 0.04325122), lik.anc=matrix(0, ncol=2, nrow=2, dimnames=list(c(1,2),c("absent", "present"))))
EPPRate <- RateMatx(BigTreeRates)
#EPPRate <- RateMatx(ace(x=EPP, phy=EPPTree, type="discrete", model="ER"))
BrownieDataGen(EPPTree, EPP, OCswitch[-EPPHold], title="EPP on OC", nsim=1300, rater=EPPRate)


pdf("BrownieOCpolyEPPswitched.pdf")
par(mfrow=c(1,2), mgp=c(1.5,.5,0), mar=c(3,3,2,1))
data <- read.csv("Polygyny on OC.csv")
NAind <- which(is.na(data$ARDRate1))
if(length(NAind)>0){
  MAX <- data$ARDRate1[-NAind]
}else{ MAX <- max(data$ARDRate1)}
BrowniePlotRates(data, "Song Stability", Groups=c("Mono", "Poly"),
                 Xlim=c(0, max(.2, MAX)))
data <- read.csv("EPP on OC.csv")
NAind <- which(is.na(data$ARDRate1))
if(length(NAind)>0){
  MAX <- data$ARDRate1[-NAind]
}else{ MAX <- max(data$ARDRate1)}
BrowniePlotRates(data, "Song Stability", Groups=c("Absent", "Present") )
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