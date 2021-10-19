# Calculate cellularity and ploidy for ACE

This R wrapper automatically calculates the cellularity and ploidy of partly processed QDNAseq data to get segmented absolute copy number shallowWGS (sWGS) data using ACE. 

## Getting Started

First process your lowpassWGS (lpWGS) data with [QDNAseq] [https://github.com/ccagc/QDNAseq]
to obtain a [CopyNumbersSegmented object] [https://rdrr.io/bioc/ACE/man/copyNumbersSegmented.html] for use with ACE [https://github.com/tgac-vumc/ACE]

### Prerequisites

You need to have the following packages
    
    x<-c("QDNAseq", "methods", "dplyr", "stringr", "CGHcall", "ACE", "copynumber", "zoo", "ggplot2", "tidyverse", "MALDIquant")
    lapply(x, require, character.only = TRUE)


## Usage

    ##how to run get ploidy and cellularity 


    # load functions
    source("ploidy_cellularity.R")

    #use the output from QDNAseq - copynumbersegmented
    #the output from this line copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

    copynumber_df <- copyNumbersSegmented@assayData$copynumber
    segmented_df <- copyNumbersSegmented@assayData$segmented

    #run "runACEnow" gives you cellularity and ploidy estimates and plots these out and runs ACE for you to get absolute 
    #copy number calls 
    #runACEnow(copynumber_df, segmented_df)

    #default settings 
    runACEnow(readCounts, copyNumbersSegmented,copynumber_df,segmented_df,2,0.125,2)

    # read in your altered cellularity ploidytext file as samples then run this
    samples <- read.table("Cellularity_Ploidy_calls.txt",header=T)
    runACEnewploidy(samples, copyNumbersSegmented,copynumber_df,segmented_df)



    
