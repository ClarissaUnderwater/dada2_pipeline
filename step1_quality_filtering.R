library(dada2)
library(DECIPHER) # check yes box in package list
packageVersion("DECIPHER")
packageVersion("dada2") 

setwd("/vortexfs1/scratch/ckarthauser/pellets")
list.files() # make sure what we think is here is actually here
samples <- scan("samples", what="character") # list samples in directory

###################################
# 16S
###################################

# variables with filenames for files we will generate
forward_reads <- paste0(samples, "_16S_R1_trimmed.fq.gz")
reverse_reads <- paste0(samples, "_16S_R2_trimmed.fq.gz")
filtered_forward_reads <- paste0(samples, "_16S_R1_filtered.fq.gz")
filtered_reverse_reads <- paste0(samples, "_16S_R2_filtered.fq.gz")

pdf("/vortexfs1/home/ckarthauser/dada2_output/plotQualityProfile_forward_reads_16S.pdf")
plotQualityProfile(forward_reads)
dev.off()

pdf("/vortexfs1/home/ckarthauser/dada2_output/plotQualityProfile_reverse_reads_16S.pdf")
plotQualityProfile(reverse_reads)
dev.off()

filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              reverse_reads, filtered_reverse_reads, maxEE=c(2,2),
                              rm.phix=TRUE, minLen=175, truncLen=c(280,280))

write.table(filtered_out, "/vortexfs1/home/ckarthauser/dada2_output/filtered_out_16S.tsv", sep = "\t", quote=F, col.names=NA)

pdf("/vortexfs1/home/ckarthauser/dada2_output/plotQualityProfile_filtered_fwd_reads_16S.pdf")
plotQualityProfile(filtered_forward_reads)
dev.off()

pdf("/vortexfs1/home/ckarthauser/dada2_output/plotQualityProfile_filtered_rvs_reads_16S.pdf")
plotQualityProfile(filtered_reverse_reads)
dev.off()


#############################################################################
########### 18S
#############################################################################

# variables with filenames for files we will generate
forward_reads_18S <- paste0(samples, "_18S_R1_trimmed.fq.gz")
reverse_reads_18S <- paste0(samples, "_18S_R2_trimmed.fq.gz")
filtered_forward_reads_18S <- paste0(samples, "_18S_R1_filtered.fq.gz")
filtered_reverse_reads_18S <- paste0(samples, "_18S_R2_filtered.fq.gz")

pdf("/vortexfs1/home/ckarthauser/dada2_output/plotQualityProfile_forward_reads_18S.pdf")
plotQualityProfile(forward_reads_18S)
dev.off()

pdf("/vortexfs1/home/ckarthauser/dada2_output/plotQualityProfile_reverse_reads_18S.pdf")
plotQualityProfile(reverse_reads_18S)
dev.off()

filtered_out_18S <- filterAndTrim(forward_reads_18S, filtered_forward_reads_18S,
                                  reverse_reads_18S, filtered_reverse_reads_18S, maxEE=c(2,2),
                                  rm.phix=TRUE, minLen=175, truncLen=c(280,280))

# checking filtering matrix for how much we cut off
filtered_out_18S

write.table(filtered_out_18S, "/vortexfs1/home/ckarthauser/dada2_output/filtered_out_18S.tsv", sep = "\t", quote=F, col.names=NA)

pdf("/vortexfs1/home/ckarthauser/dada2_output/plotQualityProfile_filtered_fwd_reads_28S.pdf")
plotQualityProfile(filtered_forward_reads_18S)
dev.off()

pdf("/vortexfs1/home/ckarthauser/dada2_output/plotQualityProfile_filtered_rvs_reads_18S.pdf")
plotQualityProfile(filtered_reverse_reads_18S)
dev.off()



