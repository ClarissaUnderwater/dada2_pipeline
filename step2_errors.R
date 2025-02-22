library(dada2)
library(DECIPHER) # check yes box in package list
packageVersion("DECIPHER")
packageVersion("dada2") 

setwd("/vortexfs1/scratch/ckarthauser/pellets")
list.files() # make sure what we think is here is actually here
samples <- scan("samples", what="character") # list samples in directory

print("starting error models")

# always do in between
load(file = "/vortexfs1/home/ckarthauser/dada2_output/dada2.RData")

# 17:13 started
err_forward_reads <- learnErrors(filtered_forward_reads)  
save.image("/vortexfs1/home/ckarthauser/dada2_output/dada2.RData")
print("16S fwd done")

err_reverse_reads <- learnErrors(filtered_reverse_reads)
save.image("/vortexfs1/home/ckarthauser/dada2_output/dada2.RData")
print("16S rvs done")

err_forward_reads_18S <- learnErrors(filtered_forward_reads_18S)  
save.image("/vortexfs1/home/ckarthauser/dada2_output/dada2.RData")
print("18S fwd done")

err_reverse_reads_18S <- learnErrors(filtered_reverse_reads_18S) # running
save.image("/vortexfs1/home/ckarthauser/dada2_output/dada2.RData")
print("18S rsv done")

pdf("/vortexfs1/home/ckarthauser/dada2_output/err_forward_reads_16S.pdf")
plotErrors(err_forward_reads, nominalQ=TRUE)
dev.off()

pdf("/vortexfs1/home/ckarthauser/dada2_output/err_reverse_reads_16S.pdf")
plotErrors(err_reverse_reads, nominalQ=TRUE)
dev.off()

pdf("/vortexfs1/home/ckarthauser/dada2_output/err_forward_reads_18S.pdf")
plotErrors(err_forward_reads_18S, nominalQ=TRUE)
dev.off()

pdf("/vortexfs1/home/ckarthauser/dada2_output/err_reverse_reads_18S.pdf")
plotErrors(err_reverse_reads_18S, nominalQ=TRUE)
dev.off()
