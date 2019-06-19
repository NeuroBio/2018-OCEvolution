setwd("C:/Users/karar/Desktop/new birbs")
library("xtable")
data <- read.csv("saddleback calc.csv", stringsAsFactors = FALSE)
data[29,1] <- "A-RW"
data[31,] <- NA
data <- data[c(1:2, 31, 3:30),]
data[3,1] <- "individual"

for(i in 1:ncol(data)){
  data[,i] <- as.character(data[,i])
}

print(xtable(data), include.rownames=FALSE)


fortis <- read.csv("fortis.csv", stringsAsFactors = FALSE)
print(xtable(fortis), include.rownames=FALSE)
scadens <- read.csv("scandens.csv", stringsAsFactors = FALSE)
print(xtable(scadens), include.rownames=FALSE)




simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}
SigStar <- function(Tabler, alpha=.05){
  starmie <- which(Tabler[,ncol(Tabler)]<alpha)
  if(length(starmie)>0){
    for(i in starmie){
      Tabler[i,ncol(Tabler)] <- paste0(Tabler[i,ncol(Tabler)],"*")
    }
  }
  return(Tabler)
}
PGLS <- read.csv("2019-05-01 OCpglsLinearLearning.csv", stringsAsFactors = FALSE)
PGLS[,1] <- gsub("[.]final", " ", PGLS[,1])
PGLS[,1] <- gsub("[.]", " ", PGLS[,1])
PGLS[,1] <- sapply(PGLS[,1], simpleCap)
PGLS[,7] <- 0
names(PGLS) <- c("Song Trait", "Slope", "Std Error", "T-Value", "p-Value", "lambda", "Corrected alpha")
PGLS <- PGLS[,c(1:3,6,4,7,5)]
PGLS <- PGLS[order(PGLS$`p-Value`),]
PGLS[,6] <- .05/(7:1)
for(i in 2:7){
  first <- round(PGLS[,i], digits=4)
  first[which(first <.001 & first >-.0001)] <- "<0.001"
  PGLS[,i] <- first
}
#PGLS[1,7] <- "0.0003"
PGLS <- SigStar(PGLS)
print(xtable(PGLS), include.rownames=FALSE)
