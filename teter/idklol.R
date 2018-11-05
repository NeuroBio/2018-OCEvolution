setwd("C:/Users/Kara/Documents/R/2018-OCEvolution/teter")
conse <- ConsensusTree
data <- birbs
OC <- birbs$O.C
names(OC) <- birbs$BirdtreeFormat
vari <- "EPP.final"
vari <- "Polygyny.final"
NAind <- which(is.na(data[,vari])==TRUE)
if(length(NAind) != 0){
  data <- data[-NAind,]
  OC <- OC[-NAind]
  conse <- drop.tip(conse, NAind, root.edge = 0)
}
testvar <- data[,vari]
names(testvar) <- data$BirdtreeFormat

er <- ace(x=testvar, phy=conse, type="discrete", model="ER")
ard <- ace(x=testvar, phy=conse, type="discrete", model="ARD")
results <- anova(er,ard)
ifelse(results$`Pr(>|Chi|)`[2] >= .05, rates <- RateMatx(er), rates <- RateMatx(ard))
BrownieDataGen(conse, testvar, OC, nsim=1300, vari, rates)

par(mfrow=c(1,2))
vari <- "EPP.final"
datasetFULL <- read.csv(paste0(vari,".csv"))
BrowniePlotRates(datasetFULL, vari)
vari <- "Polygyny.final"
datasetFULL <- read.csv(paste0(vari,".csv"))
BrowniePlotRates(datasetFULL, vari)
