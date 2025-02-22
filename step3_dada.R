library(dada2)
library(DECIPHER) # check yes box in package list
packageVersion("DECIPHER")
packageVersion("dada2") 

load(file = "/vortexfs1/home/ckarthauser/dada2_output/dada2_step2.RData")

setwd("/vortexfs1/scratch/ckarthauser/pellets")
list.files() # make sure what we think is here is actually here
samples <- scan("samples", what="character") # list samples in directory


print("starting derep")

## 16S 
derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
names(derep_forward) <- samples # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo")

save.image("/vortexfs1/home/ckarthauser/dada2_output/dada2.RData")

derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
names(derep_reverse) <- samples
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo")

save.image("/vortexfs1/home/ckarthauser/dada2_output/dada2.RData")

merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
                               derep_reverse, trimOverhang=TRUE, minOverlap=170)

save.image("/vortexfs1/home/ckarthauser/dada2_output/dada2.RData")

seqtab <- makeSequenceTable(merged_amplicons)
class(seqtab) # matrix
dim(seqtab) # 20 2521
write.table(seqtab, "seqtab_16S.tsv", sep = "\t", quote=F, col.names=NA)

seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T) 
sum(seqtab.nochim)/sum(seqtab) 
write.table(seqtab.nochim, "seqtab_16S_nochim.tsv", sep = "\t", quote=F, col.names=NA)

# write a summary of the filtering steps:
getN <- function(x) sum(getUniques(x))

summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
                          dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))

write.table(summary_tab, "summary_16S.tsv", sep = "\t", quote=F, col.names=NA)

save.image("/vortexfs1/home/ckarthauser/dada2_output/dada2.RData")

print("starting derep 18S")

## 18S 
derep_forward_18S <- derepFastq(filtered_forward_reads_18S, verbose=TRUE)
names(derep_forward_18S) <- samples # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
dada_forward_18S <- dada(derep_forward_18S, err=err_forward_reads_18S, pool="pseudo")

save.image("/vortexfs1/home/ckarthauser/dada2_output/dada2.RData")

derep_reverse_18S <- derepFastq(filtered_reverse_reads_18S, verbose=TRUE)
names(derep_reverse_18S) <- samples
dada_reverse_18S <- dada(derep_reverse_18S, err=err_reverse_reads_18S, pool="pseudo")

save.image("/vortexfs1/home/ckarthauser/dada2_output/dada2.RData")


merged_amplicons_18S <- mergePairs(dada_forward_18S, derep_forward_18S, dada_reverse_18S,
                                   derep_reverse_18S, trimOverhang=TRUE, minOverlap=170)

save.image("/vortexfs1/home/ckarthauser/dada2_output/dada2.RData")

seqtab_18S <- makeSequenceTable(merged_amplicons_18S)
class(seqtab_18S) # matrix
dim(seqtab_18S) # 16 7245
write.table(seqtab_18S, "seqtab_18S.tsv", sep = "\t", quote=F, col.names=NA)

seqtab_18S.nochim <- removeBimeraDenovo(seqtab_18S, verbose=T) 
sum(seqtab_18S.nochim)/sum(seqtab_18S) # 0.9654984 # in this case we barely lost any in terms of abundance
write.table(seqtab_18S.nochim, "seqtab_18S_nochim.tsv", sep = "\t", quote=F, col.names=NA)

# write a summary of the filtering steps:
getN <- function(x) sum(getUniques(x))

summary_tab_18S <- data.frame(row.names=samples, dada2_input=filtered_out_18S[,1],
                              filtered=filtered_out_18S[,2], dada_f=sapply(dada_forward_18S, getN),
                              dada_r=sapply(dada_reverse_18S, getN), merged=sapply(merged_amplicons_18S, getN),
                              nonchim=rowSums(seqtab_18S.nochim),
                              final_perc_reads_retained=round(rowSums(seqtab_18S.nochim)/filtered_out_18S[,1]*100, 1))

write.table(summary_tab_18S, "summary_18S.tsv", sep = "\t", quote=F, col.names=NA)

save.image("/vortexfs1/home/ckarthauser/dada2_output/dada2_step3.RData")

print("done step 3")
