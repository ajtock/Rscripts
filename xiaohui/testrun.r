######This function is used to get the coverage for given sliding window, eg. 10kb, 100kb.
#### dat.file is the coverage file, window is the 10kb..., chr.end is the length of chrs
### cov.col is the column of coverage, normally is 1,
### sub.read is the total number of alignment reads.
### out.file is the output file for given window size.
window.merge1 <- function(dat.file, window, chr.end, cov.col, sub.read, out.file){
  dat <- read.table(dat.file, header = F)
  dim.dat <- dim(dat)[1]
  start.loc <- seq(1, chr.end, window)
  start.len <- length(start.loc)
  coords <- seq(1,chr.end, 1)
  total.loc <- c(start.loc, chr.end)
  totalm <- NULL
  for(k in 1:start.len){
    index <- which(coords>=total.loc[k]&coords <= total.loc[k+1])
    index.len <- length(index)
    #print(index.len);
    if(index.len!=0){
      subd <- mean(dat[index[1]:index[index.len],cov.col])
      # the normalised here is using the length within given window 
      # and also the total number of reads
      meanc <- subd/sub.read
    } else{meanc <- 0}
    totalm <- c(totalm, meanc)
  }
  
  write.table(totalm, file = out.file, row.names = F, col.names = F)
  gc()
}
## This funcition is a little bit differ from the window.merge1 by adding a line to trim the outliers by a given trim num.
window.merge <- function(dat.file, window, chr.end, cov.col, sub.read, out.file, trim.num){
	dat <- read.table(dat.file, header = F)
	dim.dat <- dim(dat)[1]
	## extra line to trim the outliers for coverage higher than 100
	## refined on 23/03/2015.
	dat[dat[,cov.col]>=trim.num,cov.col] <- trim.num
	start.loc <- seq(1, chr.end, window)
	start.len <- length(start.loc)
	coords <- seq(1,chr.end, 1)
	total.loc <- c(start.loc, chr.end)
	totalm <- NULL
	for(k in 1:start.len){
		index <- which(coords>=total.loc[k]&coords <= total.loc[k+1])
		index.len <- length(index)
		#print(index.len);
		if(index.len!=0){
			subd <- mean(dat[index[1]:index[index.len],cov.col])
			meanc <- subd/sub.read
			} else{meanc <- 0}
		totalm <- c(totalm, meanc)
		}
   
   write.table(totalm, file = out.file, row.names = F, col.names = F)
   gc()
}
### The difference for window.merge2 is changing the normalised method
## Here is just using the total align reads as the window is fixed.
window.merge2 <- function(dat.file, window, chr.end, cov.col, sub.read, out.file, cutoff){
  dat <- read.table(dat.file, header = F)
  cindex <- which(dat[,1]>=cutoff)
  clen <- length(cindex)
  if(clen!=0){
  dat[cindex,1] <- 0}else{
    dat <- dat
  }
  dim.dat <- dim(dat)[1]
  start.loc <- seq(1, chr.end, window)
  start.len <- length(start.loc)
  coords <- seq(1,chr.end, 1)
  total.loc <- c(start.loc, chr.end)
  totalm <- NULL
  for(k in 1:start.len){
    index <- which(coords>=total.loc[k]&coords < total.loc[k+1])
    index.len <- length(index)
    #print(index.len);
    if(index.len!=0){
      subd <- mean(dat[index[1]:index[index.len],cov.col])
    
      # and also the total number of reads
      meanc <- subd/sub.read
    } else{meanc <- 0}
    totalm <- c(totalm, meanc)
  }
  
  write.table(totalm, file = out.file, row.names = F, col.names = F)
  gc()
}
######This is the funciton used to extract the coverage from bwa output bam file
#### The coverage chrs name is given as 1,2,3,4,...	
chrdat_function <- function(index, libfiles, covfiles){
  #un1 <- paste(rep.name, "/", rep.name, "_trim_hc11_c30.srt.bam", sep = "")
  un1 <- libfiles[index]
  reads1 <- readGAlignments(un1)
  cov1 <- coverage(reads1)
  chrs <- list()

  chrs[[1]] <- as.integer(cov1$'1'); chrs[[2]] <- as.integer(cov1$'2')
  chrs[[3]] <- as.integer(cov1$'3'); chrs[[4]] <- as.integer(cov1$'4')
  chrs[[5]] <- as.integer(cov1$'5')
  for(a in 1:5){
    out.file <- covfiles[[index]][a]
    write.table(chrs[[a]], file = out.file, row.names = F, col.names = F)
    gc()
  }	
}	
##
chrdat_pair_function <- function(index, libfiles, covfiles){
  #un1 <- paste(rep.name, "/", rep.name, "_trim_hc11_c30.srt.bam", sep = "")
  un1 <- libfiles[index]
  reads1 <- readGAlignmentPairs(un1)
  cov1 <- coverage(reads1)
  chrs <- list()
  
  chrs[[1]] <- as.integer(cov1$'1'); chrs[[2]] <- as.integer(cov1$'2')
  chrs[[3]] <- as.integer(cov1$'3'); chrs[[4]] <- as.integer(cov1$'4')
  chrs[[5]] <- as.integer(cov1$'5')
  for(a in 1:5){
    out.file <- covfiles[[index]][a]
    write.table(chrs[[a]], file = out.file, row.names = F, col.names = F)
    gc()
  }	
}	
######This is the funciton used to extract the coverage from bowtie output bam file
#### The coverage chrs name is given as Chr1, Chr2....
chrdat_bt_function <- function(index, libfiles, covfiles){
  #un1 <- paste(rep.name, "/", rep.name, "_trim_hc11_c30.srt.bam", sep = "")
  un1 <- libfiles[index]
  reads1 <- readGAlignments(un1)
  cov1 <- coverage(reads1)
  chrs <- list()
  
  chrs[[1]] <- as.integer(cov1$'Chr1'); chrs[[2]] <- as.integer(cov1$'Chr2')
  chrs[[3]] <- as.integer(cov1$'Chr3'); chrs[[4]] <- as.integer(cov1$'Chr4')
  chrs[[5]] <- as.integer(cov1$'Chr5')
  for(a in 1:5){
    out.file <- covfiles[[index]][a]
    write.table(chrs[[a]], file = out.file, row.names = F, col.names = F)
    gc()
  }	
}	
##
chrdat_pair_bt_function <- function(index, libfiles, covfiles){
  #un1 <- paste(rep.name, "/", rep.name, "_trim_hc11_c30.srt.bam", sep = "")
  un1 <- libfiles[index]
  reads1 <- readGAlignmentPairs(un1)
  cov1 <- coverage(reads1)
  chrs <- list()

  chrs[[1]] <- as.integer(cov1$'Chr1'); chrs[[2]] <- as.integer(cov1$'Chr2')
  chrs[[3]] <- as.integer(cov1$'Chr3'); chrs[[4]] <- as.integer(cov1$'Chr4')
  chrs[[5]] <- as.integer(cov1$'Chr5')
  for(a in 1:5){
    out.file <- covfiles[[index]][a]
    write.table(chrs[[a]], file = out.file, row.names = F, col.names = F)
    gc()
  }
}

#####This is used to check the total number of reads function to use for later normalised
## the coverage by the total align reads of five chrs.
## index is the number of the files you interested,
## libfiles is the full directory path of the bam files.
align_num_reads <- function(index, libfiles){
  un1 <- libfiles[index]
  reads1 <- readGAlignments(un1)
  subr <- table(seqnames(reads1))
  sumr <- sum(subr[1:5])
  sumr
}
###################################################
## plot function
chr.plot <- function(chr.ind, dat.files, cen.start, cen.end, cens, filter.win,swin){
  dlen <- length(dat.files)
  ma <- rep(1, filter.win)/filter.win
  datas <- list();
  filterd <- list();
  newfilterd <- NULL
  for(a in 1:dlen){
    datas[a]<- read.table(dat.files[a], header = F) 
    dim.dat <- length(datas[[a]])
    filterd[[a]] <- filter(datas[[a]], ma)
    submax <- max(filterd[[a]][-which(is.na(filterd[[a]]))])
    newfilterd <- c(newfilterd, submax)
    }
  ymax <- max(newfilterd)
  
  xplot <- c(1:dim.dat)
  filter.col <- colorRampPalette(c("red", "blue"))(length(dat.files))
  plot(xplot, filterd[[1]], main =paste("Chr", chr.ind,sep=""), xlab = "Coords", ylab = "coverage", type="l", col=filter.col[1], ylim=c(0,ymax))
  for(b in 2:dlen){
    par(new=T)
    lines(xplot, filterd[[b]], col=filter.col[b])
  }
  rect(cen.start[chr.ind]/swin, 0, cen.end[chr.ind]/swin, ymax, col=rgb(0.5,0.5,0.5,1/4), border = NA)
  abline(v = centrameres[chr.ind]/swin, lty = 2)
  }
#######################################
##coverage files using the 5'end instead of full coverage
##revised on 5/01/2016 to add the step dropping the spo11oligo width shorter than 20bp.

chrdat_5end_function <- function(index, libfiles, strand.col, chr.col, start.col, end.col, width.col, chr.lens, covfiles){
  un1 <- libfiles[index]
  reads1 <- readGAlignments(un1)
  newd <- cbind(as.character(seqnames(reads1)),as.character(strand(reads1)),start(reads1),end(reads1), width(reads1))
  wind <- which(newd[,width.col]<=20)
  
  newdat <- newd[-wind,]
  nstart <- ifelse(newdat[,strand.col]=="+", newdat[,start.col], newdat[,end.col])
  cov1 <- cbind(newdat[,c(chr.col,strand.col)],nstart, newdat[,width.col])
  for(i in 1:5){
    subcov1 <- cov1[cov1[,chr.col]==as.character(i),]
    tab1 <- table(subcov1[,start.col])
    ntab1 <- data.frame(tab1)
    loc <- as.numeric(as.character(ntab1[,1]))
    covs <- ntab1[,2]
    tab2 <- cbind(loc,covs)
    chr.len <- c(1:chr.lens[i])
    mind <- chr.len%in%loc
    subloc <- chr.len[mind==F]
    tab3 <- cbind(subloc, rep(0,length(subloc)))
    colnames(tab3) <- c("loc", "covs")
    finald <- rbind(tab2, tab3)
    finald <- finald[order(finald[,1]),]
    write.table(finald, file=covfiles[[index]][i], row.names=F, col.names=F)
  }
  gc()
}
####add the function for wheat specifically.20/01/2017
#### due to Tom's readBam does not allow the one replicate. 
chrdat_pair_wheat_function <- function(lib.file, out.file){
  reads1 <- readGAlignmentPairs(lib.file)
  cov1 <- coverage(reads1)
  chrs <- as.integer(cov1[[1]])
  write.table(chrs, file = out.file, row.names = F, col.names = F)
  gc()
  }	

  
