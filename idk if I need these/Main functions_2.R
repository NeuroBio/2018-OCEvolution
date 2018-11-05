########
#Coded by Cristina Robinson and Kate Synder
#Last Modified 11-8-2017
#Written in R-Studio Version 1.1.383
#R Version 3.4.2
#phytools v0.6-30     ape v4.1    maps v3.2.0
#R.utils_2.6.0  
########

################
#Functions for the 
################
LoadPrettyTree <- function(filename){
  #what it says: runs the fuctions necessary to create a nice consensus tree
  #RawTree=read.nexus("6001.tre")
  #ConsensusTree <- consensus.edges(RawTree)
  ConsensusTree <- read.nexus(filename)
  ConsensusTree <- PrepareConsensus(ConsensusTree)
  ConsensusTree <- drop.tip(ConsensusTree, "Emberiza_elegans", root.edge = 0)
  path <- PathFinder(ConsensusTree)
  ConsensusTree <- PrettyConsensus(ConsensusTree, path)
  assign("Path", path, pos=1)#this will be a global variable
  return(ConsensusTree)
}
Loadbirds <- function(filename, conse){
  #reading in the info  
  birdorder <- c(conse$tip.label)
  allinfo <- read.table(filename, sep=",", header=TRUE)
  #reorder to match the concensus tree
  allinfo[,1] <- gsub(" ", "_", allinfo[,"BirdtreeFormat"])
  reorder <- 1:length(birdorder)
  for(i in 1:length(birdorder)){
    reorder[i] <- which(allinfo[,"BirdtreeFormat"] == birdorder[i])
  }
  allinforeorder <- allinfo[reorder,]
  
  #clean up the data
  #  variable <- c("Syll.song", "Syllable.rep", "Song.rep", "Duration", "Interval")
  #  names(allinforeorder)[c(14, 17)] <- c("Duration.final", "Interval.final")
  #  for(i in 1:length(variable)){
  #    allinforeorder <- meanminmaxNA(variable[i], allinforeorder)
  #  }
  
  return(allinforeorder)
}
PrettyConsensus <- function(conse, path){
  if(missing(path)==TRUE){
    path <- PathFinder(conse)
  }
#get the full path length based on indicies above, get tip index, load into database
  speciespath <- matrix(NA, ncol=3, nrow=length(path))
  for(i in 1:length(path)){
    endedge <- max(path[[i]])
    pathleng <- sum(conse$edge.length[path[[i]]])
    finaledge <- conse$edge.length[max(path[[i]])]
    speciespath[i,]<- cbind(endedge, pathleng, finaledge)
  }
  tipstop <- max(speciespath[,2])
  increments <- (speciespath[,2]-tipstop)*-1
  conse$edge.length[speciespath[,1]] <- conse$edge.length[speciespath[,1]] + increments
  return(conse)
}

#Goes through the tree and creates a list for each species.
#This list containes the indicies for $edge that a species uses
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
  return(speciespath)
}
conse <- ConsensusTree

PathFinder2 <- function(conse){
  edges <- conse$edge
  specindex <- which(edges[,2] < edges[1,1])
  path <- rep(list(NA),length(specindex))
  path[[1]] <- 1:specindex[1]
  for(i in 2:length(specindex)){
#find new nodes to next tip and merge with last path
    new <- (specindex[i-1]+1):specindex[i]
    twopath <- c(edges[path[[i-1]]],edges[new])
#find the last common node and use it to make full node path
    common <- which(table(twopath) == 2)
#case 1 deals with the fork at the root as no path is shared
    if(names(common) == edges[1,1]){
      path[[i]] <- c(new)
    }else{
      path[[i]] <- c(path[[i-1]][1:(common-1)],new)
  }}
return(path)
}

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

#Pull the edges from the consensus and pull the first path
PrepareConsensus <-function(conse){
  #phytools becomes irate about edge lengths of 0
  conse$edge.length[conse$edge.length==0] = 0.0000001
  #the tree needs to be unrooted for phytools.
  #We used sayornis phoebe (suboscine) as the root and remove it
  #to create our unrooted tree
  conse=drop.tip(conse, 1, root.edge = 0)
  return(conse)
}

QuickHistBox <- function(vari, OC){
  vari <- log10(vari)
  open <- vari[which(OC == 1)]
  closed <- vari[which(OC == 0)]
  pdf("QuickHistBox.pdf")
  par(mfrow = c(1,2))
  plot(OC, vari, pch=20, xaxt = 'n', yaxt = 'n',xlim = c(-.5, 1.5),
       xlab = "", ylab = "Repertoire Size")
  axis(1, at=0:1, labels=c("Closed", "Open"))
  axis(2, at=seq(0,3, by=.5), labels=c(0,5,10,50,100,500,1000), las=1)
  boxplot(closed, open, ylab = "Repertoire Size",
          names=c("Closed", "Open"), yaxt = 'n')
  axis(2, at=seq(0,3, by=.5), labels=c(0,5,10,50,100,500,1000), las=1)
  dev.off()
}

DataExtraction <- function(OC, data, conse, vari, mod = "linear", cotitle="",
                           fullrate, RAIN=TRUE, ANOVA=TRUE, BROWNIE=TRUE){
#cut NAs
  NAind <- which(is.na(data[,vari])==TRUE)
  if(length(NAind) != 0){
    data <- data[-NAind,]
    OC <- OC[-NAind]
    conse <- drop.tip(conse, NAind, root.edge = 0)
  }

#get Vari
  testvar <- data[,vari]
  names(testvar) <- data$BirdtreeFormat
  if(mod == "log"){
    testvar <- log(testvar)
    testvarpretty <- Squish(testvar)
  }else(testvarpretty <- testvar)
#plot it
  if(RAIN==TRUE){
    rates <- Rainbowplot(conse, testvarpretty, OC, title=paste0(vari,cotitle))
    print("RainbowPlot Finished!")
  }else(rates <- RateMatx(ace(x=OC, phy=conse, type="discrete", model="ER")))
#test for diff
  if(ANOVA==TRUE){
    if(length(unique(testvar))>2){
      set.seed(49) #lucky number :)
      sink(file = "ANOVA.txt", append = TRUE, split = FALSE)
      print(paste0(vari,cotitle));print(paste0("Closed=", mean(testvar[which(OC=="closed")])))
      print(paste0("Open=", mean(testvar[which(OC=="open")])))
      output<-phylANOVA(conse, OC, testvar, nsim = 1000);print(paste0("Fval=",output$F))
      print(paste0("pVal=",output$Pf));print(paste("",sep="\n\n"));print(paste("",sep="\n\n"))
      sink(file = NULL);print("ANOVA Complete!")
    }else{warning("Dependant variable is binary; phylANOVA will not be generated")}
  }
#Brownie 
  if(BROWNIE == TRUE){
    Browniedata <- BrownieDataGen(conse, OC, testvar, title=paste0(vari,cotitle), nsim=1300, rater=rates)
    #Browniedata <- as.data.frame(read.csv(paste0(vari,cotitle, ".csv")))
    BrowniePlot (Browniedata, title = paste0(vari,cotitle))
    if(missing(fullrate) == FALSE){
      print("Full Rate Sims started!")
      Browniedata <- BrownieDataGen(conse, OC, testvar, title=paste0(vari,cotitle,"FULL"), nsim=1300, rater=fullrate)
      #Browniedata<- as.data.frame(read.csv(paste0(vari,cotitle, "FULL", ".csv")))
      #BrowniePlot (Browniedata, title = paste0(vari,cotitle, "FULL"))
    }
  }
}

Rainbowplot <- function(conse, testvar, OCdata, title="RainbowPlot"){
  #cheating way to offset species labels
  conse$tip.label <- paste("     ", conse$tip.label, sep = "")
  names(testvar) <- paste("     ", names(testvar), sep = "")
  names(OCdata) <- conse$tip.label
  #make pretty map
  co <- c("black","white")
  if(length(unique(testvar))>2){
    Bin <- FALSE
    #continuous data
    AB<-contMap(conse,x=testvar,type="fan",fsize=0.05,lwd=2,plot=FALSE)
    names(AB$cols)=rev(names(AB$cols))
    #AB<-setMap(AB,colors=rev(c("#000000", "#6a3a2d", "#b84a29", "#e66731", "#ed9445", "#f6c15b", "#fdeb73", "#ffffd1", "#FFFFFF")))
    AB<-setMap(AB,colors=c("#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026"))
  }else{
#Binary Data
    Bin <- TRUE
#test which model is better and pull the rates for the simmaps
    er <- ace(x=testvar, phy=conse, type="discrete", model="ER")
    ard <- ace(x=testvar, phy=conse, type="discrete", model="ARD")
    results <- anova(er,ard)
    ifelse(results$`Pr(>|Chi|)`[2] >= .05, rates <- RateMatx(er), rates <- RateMatx(ard))
    AB <- make.simmap(conse, testvar, Q=rates, nsim=1000, message=FALSE)
  }
#get the OC data
  testopenclosed = ace(x=OCdata, phy=conse, type="discrete", model="ER")
  returnrates <- RateMatx(testopenclosed)
  
  pdf(paste(title, ".pdf",sep=""))
  if(Bin == TRUE){
    AB <- densityMap(AB,fsize=0.7,legend=FALSE)
    add.color.bar(leg=5,cols=rev(AB$cols),lims=c(0,1),
                  digits=3,prompt=FALSE,x=0,
                  y=(length(conse$tip.label)-3),lwd=4,fsize=.7,title = "", states = c(1,0), subtitle=title)
  }else{
    plot(AB,lwd=2,fsize=0.7, label.offset=.5, legend=FALSE)
    add.color.bar(5,rev(AB$cols), subtitle=title, title="",
                  lims=AB$lims,digits=3,prompt=FALSE,x=0,
                  y=(length(conse$tip.label)-3),lwd=4,fsize=.7)
  }
  nodelabels(thermo=testopenclosed$lik.anc, piecol=co, cex=0.2)
  tiplabels(pch = 21, bg = co[as.numeric(OCdata)], cex = 0.7, adj = 1.5)
  dev.off()
  return(returnrates)
}
BrownieDataGen <- function(conse, OC, testvar, nsim =10, title, rater){
  print("Making simmaps for brownie.") 
  #make simmaps
  multisimmap <- make.simmap(tree=conse,x=OC, nsim=nsim, message=FALSE, Q=rater) 
  print(paste0("Simmaps generated", "  Starting a for loop with ", nsim, " loops."))
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
  }
  print("End brownie loop"); write.csv(file = paste(title,".csv",sep=""), browniedata)
  return(browniedata)
}


BrowniePlot <- function(dataset, title, Jackloss) {
  
#pull in data, remove bad runs
  #dataset<- as.data.frame(read.csv(paste0("Final.polygyny", ".csv")))
  ifelse(missing(Jackloss)==TRUE,Jackloss <- "", Jackloss <- paste0(" (", Jackloss, ")"))
  dataset <- data.frame(dataset)
  dataset<- dataset[!is.na(dataset$Pval),]
  dataset<- dataset[dataset$convergence == "Optimization has converged.",]
#get a few plotting variables
  ARDratio0to1 <- dataset$ARDRate0/dataset$ARDRate1
  equal <- which(ARDratio0to1 == 1)
  if(length(equal) != 0){
    print(equal)
    ARDratio0to1 <- ARDratio0to1[-equal]
    dataset <- dataset[-equal]
  }
  sig <- length(which(dataset$Pval < 0.05))
  tot <- length(dataset$Pval)
  
  #title <- "lol"
  #pdf(file=paste0(title, " Brownie pvals.pdf"))
  #par(mfrow=c(2,2))
#Standard histogram
  #dataspot <- hist(dataset$Pval,xlab="Brownie pValues", ylab = "Frequency", main=NULL,
  #                 breaks = seq(from=0,to=1,by=0.025))
  #segments(.05,0,.05,max(dataspot$counts), col="cyan3", lwd=1.5)
  #ARD ratio plot
  percents <- numeric(4)
  percents[1] <- round(length(which(dataset$Pval[which(ARDratio0to1 > 1)]<.05))/length(ARDratio0to1)*100, digits =1)
  percents[2] <- round(length(which(dataset$Pval[which(ARDratio0to1 < 1)]<.05))/length(ARDratio0to1)*100, digits =1)
  percents[3] <- round(length(which(dataset$Pval[which(ARDratio0to1 > 1)]>=.05))/length(ARDratio0to1)*100, digits =1)
  percents[4] <- round(length(which(dataset$Pval[which(ARDratio0to1 < 1)]>=.05))/length(ARDratio0to1)*100, digits =1)

  plot(dataset$Pval, ARDratio0to1, xlab = "Brownie pValues",ylab="ARD Ratio: Rate 0/1",log="y",
       main=paste(title, Jackloss,"\nSigni Runs/Total Runs:", sig, "/", tot), cex = .7, cex.main=1,
       col = ifelse(ARDratio0to1 > 1, "indianred2","cornflowerblue"),
       xlim = c(0,1), ylim = c(.01,100), yaxt = "n",#col = c("red", "blue", "black"),
       panel.first = {abline(v=.05, col="grey68", lwd=1);abline(h=1,col="black", lwd=1.5)})
  text(.6,c(10, .1), paste0(c("Closed","Open"), " rate faster"), col=c("indianred2","cornflowerblue"), font=2)
  axis(2, c(.01,.1,1,10,100),c(.01,.1,1,10,100), las=2)
  #axis(4, c(75, .025), paste0(c("Closed","Open"), " rate faster"))
  text(c(-.035,-.035,.15,.15), c(65, .02,65, .02), paste0(percents, "%"), cex = .8, adj=0, font=2)
}

meanminmaxNA <- function(variable, Birdtemp){
  NAind <- which(is.na(Birdtemp[,paste0(variable, ".min")])==TRUE)
  Birdtemp[NAind,paste0(variable, ".min")] <- Birdtemp[NAind,paste0(variable, ".final")]
  NAind <- which(is.na(Birdtemp[,paste0(variable, ".max")])==TRUE)
  Birdtemp[NAind,paste0(variable, ".max")] <- Birdtemp[NAind,paste0(variable, ".final")]
  return(Birdtemp)
}
Squish <- function(variable){
  quant <- quantile(variable)
  low<-quant[2] - (quant[3]-quant[2])
  high<-quant[4] + (quant[4]-quant[3])
  variable[which(variable > high)] <- high
  variable[which(variable < low)] <- low
  return(variable)
}
RateMatx <- function(rawrate){
  matsize <- length(rawrate$rates)
  #handles ER cases
  if(matsize == 1){matsize<- matsize+1}
  #create matrix of correct size; make diag negative; multiple by rates
  qrates <- matrix(rep(1,matsize^2),matsize,matsize)
  diag(qrates) <- -1
  qrates <- qrates*rawrate$rates  
  #pulls the state names and sets as col and row names
  rownames(qrates) <- dimnames(rawrate$lik.anc)[[2]]
  colnames(qrates) <- dimnames(rawrate$lik.anc)[[2]]
  return(qrates)
}

DistributionPlots <- function(dataset){
  #pull data, combine, kill dups, sort, separate out closed-ended learners
  #stuff for graphing
  #familiytyp <- unique(as.character(data$Family))
  colovec <- integer(length(dataset$Family))
  for(i in 1:length(colovec)){
    if(i == 1){
      counter <- 1
      colovec[i] <- counter
    }else{
      if(dataset$Family[i]==dataset$Family[i-1]){
        colovec[i] <- counter
      }else{
        counter <- counter+1
        colovec[i] <- counter
      }
    }
  }
  dataset <- cbind(dataset, colovec)
  dataset <- dataset[order(data$Syllable.rep.final),]
  
  Reper <- data$Syllable.rep.final
  familiy <- unique(as.character(data$Family))
  
  #Species <- data$BirdtreeFormat
  ticks <- length(data$BirdtreeFormat)
  par(mfrow=c(1,2), mar = c(5,3.5,.5,1), mgp =c(2.2,.6,0))
  legtop <- max(Reper) - .025*(max(Reper))
  col=rainbow(c(1,1,2,1))
  #plot for linear graph
  plot(0, type="n", xlim = c(1,ticks), ylim = c(min(Reper), max(Reper)),
       las = 2, xaxt = "n", xlab = "", ylab ="Repertoire Size",  cex = .6, cex.axis = .75)
  segments(2:ticks-1, Reper[2:ticks-1], 2:ticks, Reper[2:ticks],
           lty = 1, lwd = 1.5, col = "black")
  points(1:ticks, Reper, pch=19, col=rainbow(length(familiy))[dataset$colovec])
  points(Cind, Close$repsize, pch=19, col = "blue")
  axis(1, 1:ticks, labels = nodups$species.name, cex.axis = .6, las=2, font = 3)
  legend(1, legtop, legend=c("Closed Learners", "Open Learners"),
         col=c("blue", "red"), pch = 19, cex=.8)
  
  #plot for ln() graph
  par(mar=c(0,0,2.1,0), mar = c(5,3.5,.5,1), mgp =c(2.2,.6,0), fig = c(.5, 1, .5, 1), new = T)
  nodups$repsize <- log(nodups$repsize)
  Close$repsize <- log(Close$repsize)
  plot(0, type="n", xlim = c(1,ticks), ylim = c(min(nodups$repsize), max(nodups$repsize)),
       las = 2, xaxt = "n", xlab = "", ylab ="ln(Repertoire Size)", pch=21, col = "red",  cex = .6,
       cex.axis = .75,
       panel.first = {segments(loc, 0, loc, max(nodups$repsize), lty = 2, lwd = 1.5, col = "grey60")})
  segments(2:ticks-1, nodups$repsize[2:ticks-1], 2:ticks, nodups$repsize[2:ticks],
           lty = 1, lwd = 1.5, col = "black")
  points(1:ticks, nodups$repsize, pch=19, col = "red")
  points(Cind, Close$repsize, pch=19, col = "blue")
  axis(1, 1:ticks, labels = nodups$species.name, cex.axis = .6, las=2, font  = 3)
}

SpeciesSimmap <- function(Conse, Fam, Path){
#make a simmap object
  SimPaths <- Conse
  numrow <- length(SimPaths$edge.length)
  SimPaths[[5]] <- as.list(rep(NA, numrow))
  SimPaths[6] <- list(cbind(Conse$edge.length, rep(0,numrow)))
  edgenames <-  paste(Conse$edge[,1],Conse$edge[,2], sep=",")
  rownames(SimPaths[[6]]) <- edgenames 
  colnames(SimPaths[[6]]) <- c(0,1)
  SimPaths[c(7:8)] <- NA
  names(SimPaths)[5:8] <- c("maps", "mapped.edge", "Q", "LogL")
  class(SimPaths) <- c("simmap", "phylo")
#make a holder for the map data
  edgelen <- SimPaths$edge.length
  names(edgelen) <- rep(0, numrow)

#turn all of nodes for Fam into the 1 state and load into simmap object
  OneIndex <- CladeHighlight(Fam, Path)
  names(edgelen)[OneIndex] <- "1"
  SimPaths[[6]][OneIndex,] <- cbind(rep(0,length(OneIndex)),SimPaths[[6]][OneIndex,1])
  for(i in 1:numrow){SimPaths[[5]][[i]] <- edgelen[i]}
return(SimPaths)
}
CladeHighlight <- function(Fam, Path){
#set the nodes assocaited with a family to 1
#---get the index for species of interest. pull the path values and count
#---the largest value present num of species times is the first common node
  Speci <- Fam
  numSpec <- length(Fam)
  Counts <- table(unlist(Path[Speci]))
  CommonPath <- which(Counts==numSpec)
  ComnNode <- max(as.numeric(names(CommonPath)))
#pull all nodes after the common node
  PostCNodes <- rep(list(NA),numSpec)
  for(i in 1:numSpec){
    WorkingData <- unlist(Path[Speci[i]])
    NewNodes <- WorkingData[which(WorkingData > ComnNode)]
    PostCNodes[[i]] <- NewNodes
  }
  AllOneNodes <- c(unique(unlist(PostCNodes)), ComnNode)
  return(AllOneNodes)
}


