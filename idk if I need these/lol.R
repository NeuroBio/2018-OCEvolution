vari <- "lol"; cotitle <- "";testvar <- 1:20; OC <- birbs$O.C[1:20];output <- "lol"
#sink(file = "ANOVA.txt", append = TRUE, split = FALSE)
list(paste0(vari,cotitle),testvar,OC,output)


writeLines(paste0(vari,cotitle))
writeLines(paste0("Closed=", mean(testvar[which(OC=="closed")])))
writeLines(paste0("Open=", mean(testvar[which(OC=="open")])))
writeLines(paste0("Fval=",output$F))
writeLines(paste0("pVal=",output$Pf))
writeLines(paste0("Corrected Alpha =", Alpha))
writeLines(paste("",sep="\n\n"))
writeLines(paste("",sep="\n\n"))
sink(file = NULL);print("ANOVA Complete!")

rm(list("vari", "cotitle", "testvar", "OC", "output"))


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






nsim = 1300






  
  
  




#generate empty dataframe for the data
browniedata <- data.frame(Pval=numeric(nsim), ERRate=numeric(nsim), ERloglik=numeric(nsim),
                          ERace=numeric(nsim), ARDRate0=numeric(nsim), ARDRate1=numeric(nsim),
                          ARDloglik=numeric(nsim),ARDace=numeric(nsim), convergence=character(nsim),
                          stringsAsFactors = FALSE)
write.csv(file = paste(title,".csv",sep=""), browniedata)
#for each simmap, run the brownie code and extract results.  Save every 50 runs and at the end
for (i in 1:nsim) {
  Results <- NULL
  withTimeout(expr={
    Results <- brownie.lite(multisimmap[[i]],testvar, maxit=750000)
    browniedata[i,1:8] <- c(Results[[11]], Results[[1]], Results[[4]],  Results[[2]],
                            Results[[6]][1], Results[[6]][2], Results[[9]],
                            Results[[7]])
    browniedata[i,9] <-Results[[12]]},
    timeout = 5, substitute = FALSE, onTimeout = "warning")
  #if a run takes too long, quit the brownie and throw NAs
  if(is.null(Results)==TRUE){print(paste("Timeout warning; data not generated for simmap", i, sep=" "));browniedata[i,1:8] <- NA}
  if (i %in% seq(50,2000,by=50)){
    write.csv(file = paste(title,".csv",sep=""), browniedata);print("Saved data")
    print(paste("End brownie loop iteration",i,Sys.time()))
  }
  
  
  
  
Conse <- conse
Fam <- Families2[[Remove[k]]]
Path <- newPath  







conse <- ConsensusTree2
PathFinder <- function(conse){
  #there is a minor issue with the consensus tree code in that it does not
  #ensure that the full length of each path is the same (so that all of tips end at the same place).
  
  #This code first recreates the node paths from root to tip for each species as a vector of nodes (the loops).
  #Once this exist, it calculates the length of each path and finds the max length.
  #It then subtracts each path length from the max length to get a value of how much
  #length needs to be added to have the tips line up.
  #That value is added to the last edge in the path (the one connecting to the tip).
  
  #Based on my testing, this does not affect future calculations and makes the trees look better.
  #If you are unsure, comment out the "PrettyConsensus" line of LoadPrettyTree;
  #the rest of the code will run just fine. :)
  
  i <- 0
  
  numspec <- length(conse$tip.label)
  speciespath <- rep(list(NA),numspec)
  #find new nodes counting forward from i
  for (j in 1:numspec){
    i <- i+1
    startedge <- i
    while(conse$edge[i,2] > numspec){
      #conse[[1]][i,2] > numspec){
      #Counts up sequentially counting nodes
      while(conse$edge[i,2] == conse$edge[i+1,1] && i+1 <= numspec){
        i <- i+1
      }
      #Counts up nonsequentually counting internodes
      if(conse$edge[i,2] > numspec){
        i <- i+1
        
      }
    }
    endedge <- i
    nodepath <- startedge:endedge
    
    #if i did not get to the root (i=1) count backwards from i/k to find path
    if(conse$edge[startedge,1] != numspec + 1){
      k <- startedge
      while(k != 1  && conse$edge[k,1] != numspec + 1){
        #finds where to jump to in order to go back
        linker <- which(conse$edge[,2] == conse$edge[k,1])
        k <- linker
        #count down until reaches a nonsequential number greater than or equal to 1(root)
        
        while(conse$edge[k,1] == conse$edge[k-1,2] && k > 1){
          k <- k-1
        }
        nodepath <- c(nodepath, k:linker)
      }
    }
    
    #get the full path length based on indicies above, get tip index, load into database
    nodepath <- sort(nodepath)
    speciespath[j]<- list(nodepath)
    names(speciespath[j]) <- conse$tip.label[j]
  }
}
PathFinder(conse)
plot(conse)
  