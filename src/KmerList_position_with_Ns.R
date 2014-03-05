KmerList_position = function(kmer_length,pwd,fasta){
#################
# This function creates a list of names kmers and elements kmer position in a given genome, and saves it in the folder of interest
# kmer_length is the length of kmer we want to look for in the genome
# pwd is the folder in which we want to save the output of the function
# fasta is the location and name of fasta file of the genome of interest
# Positive strand is negative values of positions, negative strand is positive values of positions
#################
library(seqinr)
library(gtools)
library(stringr)
options(stringsAsFactors = F)

setwd(pwd)
ref = read.fasta(fasta)
refpos = unlist(getElement(ref,1))
l = length(refpos)
seqmat = matrix(data = NA,nrow = l,ncol = 1)
seqmat[which(refpos == "c")] = 0
seqmat[which(refpos == "t")] = 1
seqmat[which(refpos == "a")] = 2
seqmat[which(refpos == "g")] = 3
seqmat[which(refpos != "c" & refpos != "a" & refpos != "t" & refpos != "g")] = 4
refpos2 = seqmat
l2 = length(refpos2) 
refpos3 = matrix(data = NA,nrow = length(refpos2),ncol = 1)
inda = which(refpos2 == 1)
indc = which(refpos2 == 3)
indg = which(refpos2 == 0)
indt = which(refpos2 == 2)
indn = which(refpos2 == 4)
refpos3[inda] = 2
refpos3[indc] = 0
refpos3[indt] = 1
refpos3[indg] = 3
refpos3[indn] = 4

kmerList = sapply(seq(1,(l-kmer_length+1),1), function(x){
  motifpos = paste(refpos2[x:(x+kmer_length-1)], collapse = '')
  motifneg = paste(refpos3[(x+kmer_length-1):x], collapse = '')
  c(motifpos,-x,motifneg,(x+kmer_length-1))
})

kmerList = t(kmerList)
kmer_with_no_N = gregexpr("4",kmerList[,1])
kmer_with_no_N=lapply(kmer_with_no_N,function(x){
	if (x==-1){TRUE}else({FALSE})
	})
kmerList=kmerList[unlist(kmer_with_no_N),]
kmerMat = data.frame(rbind(cbind(kmerList[,2],kmerList[,1]),cbind(kmerList[,4],kmerList[,3])))
position_per_kmer = unstack(kmerMat)

save(position_per_kmer,file = paste("KmerList_position_",as.character(kmer_length),".RData",sep = ""))
}
