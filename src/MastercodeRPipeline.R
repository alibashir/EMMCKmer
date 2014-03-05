rm(list = ls())
##########
# run entire R part from bash script, see help in modified kmer detection for description of input variables
# for test:
# pwd_bacteria_reads_nat = "/projects/KmerIPD/test/nativedata.txt"
# pwd_bacteria_reads_wga = "/projects/KmerIPD/test/wgadata.txt"
# bacteria_name = "Geobacter_metallireducens"
# kmer_length = 4:7
# fasta = "/projects/KmerIPD/test/First_50K.fasta"
# pwd = "/projects/KmerIPD/test/"
###########
args = commandArgs(trailingOnly = T)
pwd_bacteria_reads_nat = args[[1]]
pwd_bacteria_reads_wga = args[[2]]
bacteria_name = args[[3]]
kmer_length_start = as.numeric(args[[4]])
kmer_length_end = as.numeric(args[[5]])
fasta = args[[6]]
pwd = args[[7]]

source("src/modified_kmer_detection.R")
modified_kmer_detection(pwd_bacteria_reads_nat,pwd_bacteria_reads_wga,bacteria_name,kmer_length_start,kmer_length_end,fasta,pwd)
