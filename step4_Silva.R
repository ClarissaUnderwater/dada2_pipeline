library(dada2)
library(DECIPHER) # check yes box in package list
packageVersion("DECIPHER")
packageVersion("dada2") 

load(file = "/vortexfs1/home/ckarthauser/dada2_output/dada2_step3.RData")

setwd("/vortexfs1/scratch/ckarthauser/pellets")
list.files() # make sure what we think is here is actually here
samples <- scan("samples", what="character") # list samples in directory


print("loading data base")

load("/vortexfs1/home/ckarthauser/SILVA_SSU_r138_2019.RData") 

########################################################

dna_16S <- DNAStringSet(getSequences(seqtab.nochim))
dna_18S <- DNAStringSet(getSequences(seqtab_18S.nochim))

########################################################
print("starting 16S taxonomy")

#### super slow!
ids_16S <- IdTaxa(dna_16S, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
####

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "/vortexfs1/home/ckarthauser/dada2_output/ASVs_16S.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "/vortexfs1/home/ckarthauser/dada2_output/ASVs_counts_16S.tsv", sep="\t", quote=F, col.names=NA)
asv_16S <- read.table("/vortexfs1/home/ckarthauser/dada2_output/ASVs_counts_16S.tsv", sep="\t", header=TRUE, quote="\"", stringsAsFactors=FALSE)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax <- t(sapply(ids_16S, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(asv_tax, "/vortexfs1/home/ckarthauser/dada2_output/ASVs_taxonomy_16S.tsv", sep = "\t", quote=F, col.names=NA)
asv_tax_16S <- read.table("/vortexfs1/home/ckarthauser/dada2_output/ASVs_taxonomy_16S.tsv", sep="\t", header=TRUE, quote="\"", stringsAsFactors=FALSE)

########################################################
print("starting 18S taxonomy")

ids_18S <- IdTaxa(dna_18S, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
# super slow!

asv_seqs_18S <- colnames(seqtab_18S.nochim)
asv_headers_18S <- vector(dim(seqtab_18S.nochim)[2], mode="character")

for (i in 1:dim(seqtab_18S.nochim)[2]) {
  asv_headers_18S[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta_18S <- c(rbind(asv_headers_18S, asv_seqs_18S))
write(asv_fasta_18S, "/vortexfs1/home/ckarthauser/dada2_output/ASVs_18S.fa")

# count table:
asv_18S <- t(seqtab_18S.nochim)
row.names(asv_18S) <- sub(">", "", asv_headers_18S)
write.table(asv_18S, "/vortexfs1/home/ckarthauser/dada2_output/ASVs_counts_18S.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax_18S <- t(sapply(ids_18S, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(asv_tax_18S) <- ranks
rownames(asv_tax_18S) <- gsub(pattern=">", replacement="", x=asv_headers_18S)

write.table(asv_tax_18S, "/vortexfs1/home/ckarthauser/dada2_output/ASVs_taxonomy_18S.tsv", sep = "\t", quote=F, col.names=NA)

print("done with Silva step 4")

save.image("/vortexfs1/home/ckarthauser/dada2_output/dada2_step4.RData")
