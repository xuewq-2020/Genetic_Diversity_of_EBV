# Set the working directory
setwd("/Users/data/EBV_genoms/")

# Load the Biostrings library
library(Biostrings)

#get file list of the EBV genomes 
fastafile_name = list.files(,pattern=".fasta")

#get both ends of the cds for target gene (eg.BALF2)
BALF2 <- readDNAStringSet("NC007605_BALF2.txt")
tail_cds_balf2<-reverseComplement(subseq(BALF2,1,70))
head_cds_balf2<-reverseComplement(subseq(BALF2,-70,-1))


#extract the coding sequences of target gene (eg.BALF2) from EBV genomes
  BALF2.ext<-function(ebv){
  a1<-vmatchPattern(head_cds_balf2$NC007605_BALF2,ebv,fixed=TRUE,max.mismatch = 10)
  a2<-vmatchPattern(tail_cds_balf2$NC007605_BALF2,ebv,fixed=TRUE,max.mismatch = 10)
  BALF2.re<-subseq(ebv,start=a1@ends[[1]]-a1@width0+1,end=a2@ends[[1]]) 
  BALF2<- reverseComplement(BALF2.re)
  writeXStringSet(BALF2,"BALF2_cds_seq_2.fasta",format = "fasta",append = TRUE)}

#extract the cds of target gene (eg.BALF2) from ALL EBV genomes
for(i in 1:length(fastafile_name)) {
BALF2.ext(readDNAStringSet(fastafile_name[i]))} 

#translate of nucleotide sequences to protein sequences
BALF2 <- readDNAStringSet("BALF2_cds_seq.fasta")
BALF2.protein<-translate(BALF2,if.fuzzy.codon = "X")
writeXStringSet(BALF2.protein,"BALF2_protein_seq.fasta",format = "fasta",append = TRUE) 