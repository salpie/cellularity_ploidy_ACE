##how to run get ploidy and cellularity 


# load functions
source("ploidy_cellularity.R")

#use the output from QDNAseq - copynumbersegmented
#the output from this line copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
copynumber_df <- copyNumbersSegmented@assayData$copynumber
segmented_df <- copyNumbersSegmented@assayData$segmented

#run "runACEnow" gives you cellularity and ploidy estimates and plots these out and runs ACE for you to get absolute 
#copy number calls 
runACEnow(copynumber_df, segmented_df)

