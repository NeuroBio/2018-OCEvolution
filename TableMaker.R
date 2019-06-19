if(!require("xtable")){
  install.packages("xtable")
}
library("xtable")

setwd("D:/Documents/R/2018-OC-EvolutionDataWarehouse")

ANOVATable <- function(Path, Print=TRUE){
  ANOVA <- read.table(Path, stringsAsFactors = FALSE)[,1]
  ANOVA <- ANOVA[-which(ANOVA=="Corrected")]
  ANOVA <- matrix(ANOVA,ncol=6, byrow = TRUE)
  for(i in 2:6){
    ANOVA[,i] <- sapply(strsplit(ANOVA[,i],"="), function(x) x[2])
  }
  ANOVA[,1] <- gsub("[.]", " ",ANOVA[,1])
  colnames(ANOVA) <- c("Variable","Stable", "Plastic", "F-Value",
                       "p-Value", "Corrected Alpha")
  ANOVA <- ANOVA[order(ANOVA[,5]),]
  Tag <- ifelse(ANOVA[,5] < ANOVA[,6],"*", "")
  ANOVA[,5] <- paste0(ANOVA[,5], Tag)
  ANOVA <- ANOVA[,c(1:4,6,5)]
  if(Print==TRUE){
    print(xtable(ANOVA), include.rownames=FALSE)
  }else{
    return(ANOVA)
  }
}
BrownieTable <- function(Path, Cols=4, Doubs=TRUE, Print=TRUE){
  BROWNIETable <- read.table(Path, header = FALSE, fill=TRUE, row.names = NULL,
                             stringsAsFactors = FALSE)[,1]
  Remove <- which(BROWNIETable == "[1]")
  if(length(Remove)!=0){
    BROWNIETable <- BROWNIETable[-Remove] 
  }
  BROWNIETable <- matrix(BROWNIETable, ncol=Cols,byrow = TRUE)
  if(Doubs == TRUE){
    BROWNIETable <- BROWNIETable[seq(2,nrow(BROWNIETable),2),]
    BROWNIETable[,1] <- gsub("FULL","",BROWNIETable[,1])
  }
  BROWNIETable[,1] <- gsub("[.]"," ",BROWNIETable[,1])
  if(Cols == 4){
    colnames(BROWNIETable) <- c("Group", "One-Rate Model", "Two-Rate Model", "p-Value")  
    for(i in 2:Cols){
      BROWNIETable[,i] <- sapply(strsplit(BROWNIETable[,i],"="), function(x) x[2])
    }
  }else{
    colnames(BROWNIETable) <- c("Group", "Removed Species", "One-Rate Model", "Two-Rate Model", "p-Value")
    for(i in c(3:5)){
      BROWNIETable[,i] <- sapply(strsplit(BROWNIETable[,i],"="), function(x) x[2])
    }
    BROWNIETable[,2] <- gsub("_", " ",BROWNIETable[,2])
  }
  BROWNIETable[,Cols] <- ifelse(BROWNIETable[,Cols]<0.05, paste0(BROWNIETable[,Cols],"*"), BROWNIETable[,Cols])
  BROWNIETable[,Cols] <- ifelse(is.na(BROWNIETable[,Cols]), "<0.001*", BROWNIETable[,Cols])
  if(Print == TRUE){
    print(xtable(BROWNIETable), include.rownames=FALSE) 
  }else{
    return(BROWNIETable)
  }
}

ANOVATable(file.path(getwd(),"ANOVA.txt"))
BrownieTable(file.path(getwd(),"Brownie.txt"))

setwd("C:/Users/Kara/Desktop/DataWarehouse/Jackknife")
Folders <- list.files()
SubGrouping <- rep(c("None","Acrocephalidae","Emberizidae",
      "Fringillidae", "Icteridae", "Mimidae",
      "Muscicapidae", "Parulidae", "Passerellidae"),length(Folders))
ATablelist <- as.list(1:length(Folders))
BTablelist <- ATablelist
for(i in seq_along(Folders)){
  ATablelist[[i]] <- ANOVATable(file.path(getwd(), Folders[i],"ANOVA.txt"),FALSE)
  BTablelist[[i]] <- BrownieTable(file.path(getwd(), Folders[i],"Brownie.txt"),4, TRUE, FALSE)
}
ATableMat <- as.data.frame(do.call(rbind, ATablelist), stringsAsFactors = FALSE)
BTableMat <- as.data.frame(do.call(rbind, BTablelist), stringsAsFactors = FALSE)
ATableMat[,1] <- sapply(strsplit(ATableMat[,1]," "),
                        function(x) paste0(unlist(x)[-length(unlist(x))], collapse = " "))
ATableMat[,"Family Removed"] <- SubGrouping[2:length(FamilyRemove)]
BTableMat[,"Family Removed"] <- SubGrouping
print(xtable(ATableMat[,c(1,7,2:6)]), include.rownames=FALSE)
print(xtable(BTableMat[,c(1,5,2:4)]), include.rownames=FALSE)


setwd("C:/Users/Kara/Desktop/DataWarehouse/MimidJackknife")
BrownieTable(file.path(getwd(),"Brownie2.txt"),5, FALSE)
Path <- file.path(getwd(), Folders[1],"Brownie.txt")

setwd("C:/Users/Kara/Desktop/DataWarehouse/ClosedLink")
ANOVATable(file.path(getwd(),"ANOVA.txt"))
BrownieTable(file.path(getwd(),"Brownie.txt"))
