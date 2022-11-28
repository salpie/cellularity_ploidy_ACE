#calculate ploidy and cellularity from lowpassWGS
### Load Libraries
x<-c("QDNAseq", "methods", "dplyr", "stringr", "CGHcall", "ACE", "copynumber", "zoo", "ggplot2", "tidyverse", "MALDIquant")
lapply(x, require, character.only = TRUE)

# functions

peaky <- function(x, y, w=1, ...) {
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}


calculateCellPloid <- function(copynumber_df,w,span,lowest) {

	copynumber_df <- as.data.frame(copynumber_df)

	# get names of samples
	samples <- names(copynumber_df)
	#make column for ploidy and cellularity
	samples <- as.data.frame(samples)
	samples$ploidy <- 0
	samples$cellularity <- 0


	# calcualting cellularity and ploidy
	for (samp in 1:ncol(copynumber_df)) {

		print(paste("getting cellularity and ploidy from sample", names(copynumber_df)[samp]))
	
		dh <- hist(copynumber_df[,samp],plot=FALSE, breaks=500)

		y = dh$counts
		x = dh$mids

		peaks <- peaky(x, y, w=w, span=span) #get peaks from noisy data

		logRs_of_peaks = x[peaks$i[peaks$y.hat[peaks$i]>5]] #get distance between peaks to use as cellularity
		if (length(logRs_of_peaks)==1) {
			samples$cellularity[samp] = 0.02 #set cellularity to be extremely low if there is only one peak
			samples$ploidy[samp] = 2
		} else {
			mean_distance = mean(diff(logRs_of_peaks)) # cellularity
			samples$cellularity[samp] = mean_distance
			samples$ploidy[samp] = which.max(peaks$y.hat[peaks$i[peaks$y.hat[peaks$i]>5]])	# find location of the highest peak to guage ploidy
			# ploidy based on how many peaks before diploid peak
		}

		pdf(paste("Cellularity_Ploidy_", samples$samples[samp], ".pdf", sep=""), width=10, height=8)
		plot_this <- plot(x, y, cex=0.75, col="Gray", main=paste("w = ", w, ", span = ", span, sep=""))
		print(plot_this)
		lines(x, peaks$y.hat,  lwd=2) #$
		y.min <- min(y)
		sapply(peaks$i, function(i) lines(c(x[i],x[i]), c(y.min, peaks$y.hat[i]), col="Red", lty=2))
		points(x[peaks$i[peaks$y.hat[peaks$i]>lowest]], peaks$y.hat[peaks$i[peaks$y.hat[peaks$i]>lowest]], col="Red", pch=19, cex=1.25)
		dev.off()
	}

	# scale cellularity 
	samples$cellularity <- sapply(samples$cellularity, function(i) (((i - min(samples$cellularity)) * (0.8 - 0.3)) / (max(samples$cellularity) - min(samples$cellularity))) + 0.3)
	samples$cellularity[samples$cellularity==0.3] <- 0
	return(samples)
}

runACEnow <- function(readCounts=1, copyNumbersSegmented,copynumber_df,segmented_df,w=1,span=0.125,lowest=3) {
	segmented_df <- as.data.frame(segmented_df)
	copynumber_df <- as.data.frame(copynumber_df)
	# Split out row names by chr, start, end
	locations_split = lapply(strsplit(rownames(copynumber_df), "-"), function(i) unlist(strsplit(i, split = ":")))

	# Add location columns
	chr   = as.numeric(unlist(lapply(locations_split, function(i) i[1])))
	start = as.numeric(unlist(lapply(locations_split, function(i) i[2])))
	stop   = as.numeric(unlist(lapply(locations_split, function(i) i[3])))
	test_chr = na.omit(chr)
	if (sum(test_chr>23) > 1) {
		info = readCounts@featureData@data
		chr   = as.data.frame(info$chr)
		start = as.data.frame(info$start)
		stop = as.data.frame(info$end)
	}

	bin = 1:nrow(copynumber_df)
	joint_calls <- cbind(bin, chr, start, stop, copynumber_df)
	samples = calculateCellPloid(copynumber_df,w,span,lowest)

	if (class(copyNumbersSegmented)=="QDNAseqCopyNumbers"){
		for (sample in 1:nrow(samples)){
			pdf(paste("ACE_", samples$samples[sample], ".pdf", sep=""), width=10, height=8)
			printplot = ACEcall(copyNumbersSegmented,sample,cellularity=samples$cellularity[sample],ploidy = samples$ploidy[sample],title = paste(samples$samples[sample], ",     Ploidy = ", samples$ploidy[sample],sep=""), standard = 1, cap=6)$calledplot
			print(printplot)
			dev.off()
		}
	}

	#change ploidy
	samples$ploidy[samples$ploidy==1] <- 2

	#save the ploidy and cellularity calls
	write.table(samples, file="Cellularity_Ploidy_calls.txt", sep="\t", row.names=F, quote=F)


	for (i in 5:ncol(joint_calls)) {

		print(paste("getting ACE calls for sample", names(joint_calls)[i]))

		sample_ID <- i-5
		sample1 = cbind(joint_calls[,c(1:4,i)], segmented_df[,i-4])
		sample1 <- sapply(sample1, as.numeric)
		sample1 = as.data.frame(sample1)
		names(sample1) = c("bins", "chr", "start", "end", "copynumbers", "segments")
	
		#get info for segments and copynumber etc
		segmentdf <- getadjustedsegments(sample1, cellularity = samples$cellularity[i-4], ploidy=samples$ploidy[i-4])
		new_segments_ACE <- as.data.frame(lapply(segmentdf, rep, segmentdf$Num_Bins))
		write.table(segmentdf, file=paste(names(joint_calls)[i], "_single_pcf_segments.txt", sep=""), row.names=F, sep="\t", quote=F)

		# plot
		pdf(paste("ACE_copynumber_calls_", names(joint_calls)[i], ".pdf", sep=""), width=10, height=8)
		plotted <- singleplot(sample1, cellularity = samples$cellularity[i-4], ploidy = samples$ploidy[i-4])
		print(plotted)
		dev.off()
	}
}


runACEnewploidy <- function(samples, copyNumbersSegmented,copynumber_df,segmented_df, squaremodel_t) {
	segmented_df <- as.data.frame(segmented_df)
	copynumber_df <- as.data.frame(copynumber_df)
	# Split out row names by chr, start, end
	locations_split = lapply(strsplit(rownames(copynumber_df), "-"), function(i) unlist(strsplit(i, split = ":")))

	# Add location columns
	#chr   = as.numeric(unlist(lapply(locations_split, function(i) i[1])))
	#start = as.numeric(unlist(lapply(locations_split, function(i) i[2])))
	#stop   = as.numeric(unlist(lapply(locations_split, function(i) i[3])))
	#if (sum(chr>23) > 1) {
		info = readCounts@featureData@data
		chr   = as.data.frame(info$chr)
		start = as.data.frame(info$start)
		stop = as.data.frame(info$end)
	#}

	bin = 1:nrow(copynumber_df)
	joint_calls <- cbind(bin, chr, start, stop, copynumber_df)

	for (i in 5:ncol(joint_calls)) {

		print(paste("getting ACE calls for sample", names(joint_calls)[i]))

		sample_ID <- i-5
		sample1 = cbind(joint_calls[,c(1:4,i)], segmented_df[,i-4])
		sample1 <- sapply(sample1, as.numeric)
		sample1 = as.data.frame(sample1)
		names(sample1) = c("bins", "chr", "start", "end", "copynumbers", "segments")
		
		if (squaremodel_t==T) {
			sqmodel <- squaremodel(sample1)
		
			samples$cellularity[i-4] <- sqmodel$minimadf$cellularity[which.min(sqmodel$minimadf$error)]
			samples$ploidy[i-4] <- sqmodel$minimadf$ploidy[which.min(sqmodel$minimadf$error)]
		
			samples$cell[i-4] <- sqmodel$minimadf$cellularity[which.min(abs((samples$cellularity[i-4]-samples$cellularity[i-4])+sqmodel$minimadf$error)+abs(samples$ploidy[i-4]-sqmodel$minimadf$ploidy))]
		
			samples$ploidy[i-4] <- sqmodel$minimadf$ploidy[which.min(abs((samples$ploidy[i-4]-samples$ploidy[i-4])+sqmodel$minimadf$error)+abs(samples$ploidy[i-4]-sqmodel$minimadf$ploidy))]

			pdf(paste("ACE_best_calls_", names(joint_calls)[i], ".pdf", sep=""), width=10, height=8)
			printplot=ACEcall(copyNumbersSegmented, QDNAseqobjectsample = i-4,cellularity=samples$cell[i-4],ploidy = sqmodel$minimadf$ploidy[which.min(abs((samples$cellularity[i-4]-samples$cellularity[i-4])+sqmodel$minimadf$error)+abs(samples$ploidy[i-4]-sqmodel$minimadf$ploidy))],title = paste(names(joint_calls)[i], sep="_"), standard = 1, cap=6)$calledplot
			print(printplot)
			dev.off()
		} else {
		actual_info=samples[names(joint_calls)[i]==samples$samples,]

		pdf(paste("ACE_calls_yours_", names(joint_calls)[i], ".pdf", sep=""), width=10, height=8)
		printplot=ACEcall(copyNumbersSegmented, QDNAseqobjectsample = which(names(joint_calls)[i]==copyNumbersSegmented$name),cellularity= actual_info$cellularity,ploidy = actual_info$ploidy,title = paste(names(joint_calls)[i], sep="_"), standard = 1, cap=6)$calledplot
		print(printplot)
		dev.off()

	
		#get info for segments and copynumber etc
		segmentdf <- getadjustedsegments(sample1, cellularity = actual_info$cellularity, ploidy=actual_info$ploidy)
		new_segments_ACE <- as.data.frame(lapply(segmentdf, rep, segmentdf$Num_Bins))
		write.table(segmentdf, file=paste(names(joint_calls)[i], "_single_pcf_segments.txt", sep=""), row.names=F, sep="\t", quote=F)

		# plot
		pdf(paste("ACE_copynumber_calls_", names(joint_calls)[i], ".pdf", sep=""), width=10, height=8)
		plotted <- singleplot(sample1, cellularity = actual_info$cellularity, ploidy = actual_info$ploidy)
		print(plotted)
		dev.off()
		}
	}
}

runACEnewsampleploidy <- function(samples, copyNumbersSegmented,copynumber_df,segmented_df,samplenam) {
	segmented_df <- as.data.frame(segmented_df)
	copynumber_df <- as.data.frame(copynumber_df)
	# Split out row names by chr, start, end
	locations_split = lapply(strsplit(rownames(copynumber_df), "-"), function(i) unlist(strsplit(i, split = ":")))

	bin = 1:nrow(copynumber_df)
	joint_calls <- cbind(bin, chr, start, stop, copynumber_df)


		print(paste("getting ACE calls for sample", samplenam))

		sample_ID <- which(names(joint_calls)== samplenam)
		i=which(names(joint_calls)== samplenam)
		sample1 = cbind(joint_calls[,c(1:4,i)], segmented_df[,i-4])
		sample1 <- sapply(sample1, as.numeric)
		sample1 = as.data.frame(sample1)
		names(sample1) = c("bins", "chr", "start", "end", "copynumbers", "segments")
		
		actual_info=samples[names(joint_calls)[i]==samples$samples,]

		pdf(paste("ACE_calls_yours_", names(joint_calls)[i], ".pdf", sep=""), width=10, height=8)
		printplot=ACEcall(copyNumbersSegmented, QDNAseqobjectsample = which(names(joint_calls)[i]==copyNumbersSegmented$name),cellularity= actual_info$cellularity,ploidy = actual_info$ploidy,title = paste(names(joint_calls)[i], sep="_"), standard = 1, cap=6)$calledplot
		print(printplot)
		dev.off()

	
		#get info for segments and copynumber etc
		segmentdf <- getadjustedsegments(sample1, cellularity = actual_info$cellularity, ploidy=actual_info$ploidy)
		new_segments_ACE <- as.data.frame(lapply(segmentdf, rep, segmentdf$Num_Bins))
		write.table(segmentdf, file=paste(names(joint_calls)[i], "_single_pcf_segments.txt", sep=""), row.names=F, sep="\t", quote=F)

		# plot
		pdf(paste("ACE_copynumber_calls_", names(joint_calls)[i], ".pdf", sep=""), width=10, height=8)
		plotted <- singleplot(sample1, cellularity = actual_info$cellularity, ploidy = actual_info$ploidy)
		print(plotted)
		dev.off()

}


