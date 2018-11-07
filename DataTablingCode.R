library(xtable)
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}


setwd("D:/Documents/R/2018-OC-EvolutionDataWarehouse")
Data <- read.table("ANOVA.txt", stringsAsFactors = FALSE)
Form <- matrix(Data[,1], ncol=7,byrow=TRUE)
Form <- Form[,-6]
Form <- gsub("Closed=", "", Form)
Form <- gsub("Open=", "", Form)
Form <- gsub("Fval=", "", Form)
Form <- gsub("pVal=", "", Form)
Form <- gsub("Alpha=", "", Form)
Form[,1] <- gsub("[.]", " ", Form[,1])
Form[,1] <- sapply(Form[,1], simpleCap)
Form <- Form[,c(1:4,6,5)]
colnames(Form) <- c("Song Trait", "Closed", "Open", "F-Value", "Corrected alpha","p-Value")
Final <- c(1,4,7,10,13,16,17)
Finals <- Form[Final,]
Finals[,1] <- gsub(" Final", "", Finals[,1])
Finals <- Finals[order(as.numeric(Finals[,6])),]
print(xtable(Finals), include.rownames=FALSE)

MinMax <- Form[-Final,]
MinMax <- MinMax[order(as.numeric(MinMax[,6])),]
print(xtable(MinMax), include.rownames=FALSE)




Data <- read.table(file.path("Brownie.txt"), stringsAsFactors = FALSE, sep="\n")
Form <- matrix(Data[,1], ncol=4,byrow=TRUE)
RowNames <- strsplit(Form[,1], "[.]")
FinalNames <- character(length(RowNames))
for(j in seq_along(FinalNames)){
  FinalNames[j] <- tail(RowNames[[j]], n=1)
}
dump <- which(FinalNames %in% c("final", "finalFULL"))
Form <- Form[-dump,]
Form <- Form[seq(2,nrow(Form),by=2),]
Form <- gsub("FULL", "", Form)
Form <- gsub("OneRate=", "", Form)
Form <- gsub("TwoRates=", "", Form)
Form <- gsub("pVal[=,<]?", "", Form)
Form[,1] <- gsub("[.]", " ", Form[,1])
Form[,1] <- sapply(Form[,1], simpleCap)
colnames(Form) <- c("Song Trait", "One Rate", "Two Rates","p-Value")
Form <- Form[order(as.numeric(Form[,4])),]
print(xtable(Form), include.rownames=FALSE)




setwd("D:/Documents/R/2018-OC-EvolutionDataWarehouse/Jackknife")
for(i in list.files()){
  Data <- read.table(file.path(i,"ANOVA.txt"), stringsAsFactors = FALSE)
  Form <- matrix(Data[,1], ncol=7,byrow=TRUE)
  Form <- Form[,-6]
  Form <- gsub("Closed=", "", Form)
  Form <- gsub("Open=", "", Form)
  Form <- gsub("Fval=", "", Form)
  Form <- gsub("pVal=", "", Form)
  Form <- gsub("Alpha=", "", Form)
  
  Form[,1] <- gsub("[.]", " ", Form[,1])
  RowNames <- strsplit(Form[,1], " ")
  FinalNames <- character(length(RowNames))
  for(j in seq_along(FinalNames)){
    FinalNames[j] <- tail(RowNames[[j]], n=1)
  }
  Form[,1] <- FinalNames
  Form <- gsub("finalNo", "", Form)
  Form <- gsub("rateNo", "", Form)
  Form <- Form[,c(1:4,6,5)]
  colnames(Form) <- c("Removed Family", "Closed", "Open", "F-Value", "Corrected alpha","p-Value")
  Form <- Form[order(as.numeric(Form[,6])),]
  print(i)
  print(xtable(Form), include.rownames=FALSE)
}

for(i in list.files()){
  Data <- read.table(file.path(i,"Brownie.txt"), stringsAsFactors = FALSE, sep="\n")
  Form <- matrix(Data[(nrow(Data)/2+1):nrow(Data),1], ncol=4,byrow=TRUE)
  Form <- gsub("OneRate=", "", Form)
  Form <- gsub("TwoRates=", "", Form)
  Form <- gsub("pVal[=,<]?", "", Form)
  Form <- gsub("\\s*\\([^\\)]+\\)", "", Form)
  Form[,1] <- gsub("[.]final", " ", Form[,1])
  RowNames <- strsplit(Form[,1], " ")
  FinalNames <- character(length(RowNames))
  for(j in seq_along(FinalNames)){
    FinalNames[j] <- tail(RowNames[[j]], n=1)
  }
  Form[,1] <- FinalNames
  colnames(Form) <- c("Removed Family", "One Rate", "Two Rates","p-Value")
  Form <- Form[-1,]
  Form <- Form[order(as.numeric(Form[,4])),]
  print(i)
  print(xtable(Form), include.rownames=FALSE)
}


setwd("D:/Documents/R/2018-OC-EvolutionDataWarehouse/MimidJackknife")
Data <- read.table(file.path("Brownie2.txt"), stringsAsFactors = FALSE, sep="\n")
Form <- matrix(Data[,1], ncol=4,byrow=TRUE)
Form <- gsub("OneRate=", "", Form)
Form <- gsub("TwoRates=", "", Form)
Form <- gsub("pVal[=,<]?", "", Form)
RowNames <- strsplit(Form[,1], " ")
FinalNames <- character(length(RowNames))
for(j in seq_along(FinalNames)){
  FinalNames[j] <- tail(RowNames[[j]], n=1)
}
Form[,1] <- FinalNames
Form <- gsub("_", " ",Form)
colnames(Form) <- c("Removed Mimid", "One Rate", "Two Rates","p-Value")
Form <- Form[order(as.numeric(Form[,4])),]
print(xtable(Form), include.rownames=FALSE)


setwd("D:/Documents/R/2018-OC-EvolutionDataWarehouse/ClosedLink")
Data <- read.table("ANOVA.txt", stringsAsFactors = FALSE)
Form <- matrix(Data[,1], ncol=7,byrow=TRUE)
Form <- Form[,-6]
Form <- gsub("Closed=", "", Form)
Form <- gsub("Open=", "", Form)
Form <- gsub("Fval=", "", Form)
Form <- gsub("pVal=", "", Form)
Form <- gsub("Alpha=", "", Form)
Form[,1] <- gsub("[.]", " ", Form[,1])
Form[,1] <- sapply(Form[,1], simpleCap)
Form <- Form[,c(1:4,6,5)]
colnames(Form) <- c("Song Trait", "Closed", "Open", "F-Value", "Corrected alpha","p-Value")
Final <- c(1,4,7,10,13,16,17)
Finals <- Form[Final,]
Finals[,1] <- gsub(" Final", "", Finals[,1])
Finals <- Finals[order(as.numeric(Finals[,6])),]
print(xtable(Finals), include.rownames=FALSE)



Data <- read.table(file.path("Brownie.txt"), stringsAsFactors = FALSE, sep="\n")
Form <- matrix(Data[,1], ncol=4,byrow=TRUE)
RowNames <- strsplit(Form[,1], "[.]")
FinalNames <- character(length(RowNames))
for(j in seq_along(FinalNames)){
  FinalNames[j] <- tail(RowNames[[j]], n=1)
}
dump <- which(FinalNames %in% c("final", "finalFULL"))
Form <- Form[dump,]
Form <- Form[seq(2,nrow(Form),by=2),]
Form <- gsub("FULL", "", Form)
Form <- gsub("OneRate=", "", Form)
Form <- gsub("TwoRates=", "", Form)
Form <- gsub("pVal[=,<]?", "", Form)
Form[,1] <- gsub("[.]", " ", Form[,1])
Form[,1] <- sapply(Form[,1], simpleCap)
colnames(Form) <- c("Song Trait", "One Rate", "Two Rates","p-Value")
Form <- Form[order(as.numeric(Form[,4])),]
print(xtable(Form), include.rownames=FALSE)
