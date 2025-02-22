library(dada2)
library(DECIPHER) # check yes box in package list
packageVersion("DECIPHER")
packageVersion("dada2") 

load(file = "/vortexfs1/home/ckarthauser/dada2_output/dada2_step5.RData")

print("starting PR step b")

setwd("/vortexfs1/scratch/ckarthauser/pellets")
list.files() # make sure what we think is here is actually here
samples <- scan("samples", what="character") # list samples in directory


### version b of PR step
#Note: PR2 has different taxLevels than the dada2 default. 
#When assigning taxonomy against PR2, use the following: 
#assignTaxonomy(..., taxLevels = c("Domain","Supergroup","Division","Subdivision", "Class","Order","Family","Genus","Species"))

ranks_b <- c("Domain","Supergroup","Division","Subdivision", "Class","Order","Family","Genus","Species")
asv_tax_18S_PR2_b <- t(sapply(ids_18S_PR2, function(x) {
  m <- match(ranks_b, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

colnames(asv_tax_18S_PR2_b) <- ranks_b
rownames(asv_tax_18S_PR2_b) <- gsub(pattern=">", replacement="", x=asv_headers_18S)

write.table(asv_tax_18S_PR2_b, "/vortexfs1/home/ckarthauser/dada2_output/ASVs_taxonomy_18S_PR2_b.tsv", sep = "\t", quote=F, col.names=NA)
asv_tax_18S_PR2_b <- read.table("/vortexfs1/home/ckarthauser/dada2_output/ASVs_taxonomy_18S_PR2_b.tsv", sep="\t", header=TRUE, quote="\"", stringsAsFactors=FALSE)

save.image("/vortexfs1/home/ckarthauser/dada2_output/dada2_step5b.RData")


print("done PR step b")