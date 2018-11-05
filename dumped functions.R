#Not currently using:
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

BrownieSigTestGuts <- function(pvals, title){
  #no data
  if(length(pvals)<50){
    writeLines(paste(title, "not enough data"))
    return()
  }
  #check if uniform distribution; if not...
  #if(ks.test(pvals, "punif")$p.value < .05){
  #run a chisquare
  sig <- length(which(pvals < .05))
  notsig <- length(pvals)-sig
  results <- chisq.test(c(sig,notsig), p=c(.05,.95))
  if(results$p.value < .05 && sig > length(pvals)*.05){
    writeLines(paste(title, "skewed significant", "Chi2:", results$statistic,
                     "DF:", results$parameter, "pval:", results$p.value))
  }else{writeLines(paste(title, "not skewed significant","Chi2:", results$statistic,
                         "DF:", results$parameter, "pval:", results$p.value))}
}



BrownieSigTest <- function(dataset, title){
  sink(file = "Brownie.txt", append = TRUE, split = FALSE)
  writeLines(paste(title))
  negind <- which(dataset$ARDRate0/dataset$ARDRate1 < 1)
  #open faster data
  BrownieSigTestGuts(dataset$Pval[negind], "Open Faster:")
  #closed faster data
  BrownieSigTestGuts(dataset$Pval[-negind],"Closed Faster:")
  writeLines(paste("\r\n"))
  sink(file = NULL)
}

if(file.exists("Brownie.txt") == FALSE){
  file.create("Brownie.txt")
}

#clean up the data
#  variable <- c("Syll.song", "Syllable.rep", "Song.rep", "Duration", "Interval")
#  names(allinforeorder)[c(14, 17)] <- c("Duration.final", "Interval.final")
#  for(i in 1:length(variable)){
#    allinforeorder <- meanminmaxNA(variable[i], allinforeorder)
#  }