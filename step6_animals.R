library(dada2)
library(DECIPHER) # check yes box in package list
packageVersion("DECIPHER")
packageVersion("dada2") 

load(file = "/vortexfs1/home/ckarthauser/dada2_output/dada2_step5.RData")

setwd("/vortexfs1/scratch/ckarthauser/pellets")
list.files() # make sure what we think is here is actually here
samples <- scan("samples", what="character") # list samples in directory

print("running MZG")

############################################################
print("loading animals")

load("/vortexfs1/home/ckarthauser/dada2_output/MZG.RData")

#### super slow!
ids_18S_MZG <- IdTaxa(asv_fasta_18S,trainingSet_MZG, strand="top", processors=NULL, verbose=FALSE) # use all processors
####
write.table(ids_18S_MZG, "/vortexfs1/home/ckarthauser/dada2_output/ids_18S_MZG.tsv", sep = "\t", quote=F, col.names=NA)

print("done MZG IDs")

############################################################
# OR

taxa <- assignTaxonomy(asv_fasta_18S, trainingSet_MZG, multithread=FALSE)

write.table(taxa, "/vortexfs1/home/ckarthauser/dada2_output/taxa.tsv", sep = "\t", quote=F, col.names=NA)

#ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
############################################################


# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks_MZG <- c("rootrank", "species")
ranks_MZG <- c("rootrank", "genus", "species")

asv_tax_18S_MZG <- t(sapply(ids_18S_MZG, function(x) {
  m <- match(ranks_MZG, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

colnames(asv_tax_18S_MZG) <- ranks_MZG
rownames(asv_tax_18S_MZG) <- gsub(pattern=">", replacement="", x=asv_headers_18S)

write.table(asv_tax_18S_MZG, "/vortexfs1/home/ckarthauser/dada2_output/ASVs_taxonomy_18S_MZG.tsv", sep = "\t", quote=F, col.names=NA)

print("done MZG")

save.image("/vortexfs1/home/ckarthauser/dada2_output/pellets_dada2_step456.RData")

print("done step 4 5 6")
