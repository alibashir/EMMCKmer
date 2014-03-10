modified_kmer_detection = function(pwd_bacteria_reads_nat,pwd_bacteria_reads_wga,bacteria_name,kmer_length_start,kmer_length_end,fasta,pwd,save_intermediary_files=FALSE){
########################
# This function goes through the  pipeline of kmer modification detection from the output of the 
# python code reading and normalizing the cmph5 file to the output of the llr values per position 
# of each kmer present in the genome of interest.
#
# pwd_bacteria_reads_nat is the path to the output of csv readout of previous python code for native data
# pwd_bacteria_reads_wga is the path to the output of csv readout of previous python code for WGA data
# bacteria_name is the name of the bacteria of interest we are analyzing
# kmer_length is a vector of length of kmer we want to consider for this analyses
# fasta is the path and file name to the fasta file of the bacterial DNA sequence
# pwd is the path in which the output of every step of this function will be saved
# save_intermediary_files is a boolean that specify whether you want to save intermediary output files or not. 
# By default it is set to FALSE
# Intermediary outputs that would be saved along the way: 
# - the csv files transformed for R usage in "BacteriaName"_native_reads.RData and "BacteriaName"_wga_reads.RData
# - the kmer positions per kmer size in KmerList_position_"KmerLength".RData
# - the unstacked IPD data per position for both native and WGA in unstacked_ipds.RData
# - the median IPD per position for both native and WGA in median_of_reads_per_position.RData
# - the llr values per kmer position in llr_values_"KmerLength".RData
# - a table of sorted max llr per kmer in llr_median_"KmerLength".txt
# - the median IPD value per kmer position in median_of_medians_"KmerLength".RData
# - a table of the llr value per kmer per letter in llr_median_"KmerLength"_each_letter.txt
# - the catted txt file of the previous output, also input for the next python code in catkmer_llr_"BacteriaName".txt
#########################

#####
# for test:
#####
# run python code on: "/projects/KmerIPD/test/out_native_final.cmp.h5" and "/projects/KmerIPD/test/out_wga_final.cmp.h5"
# pwd_bacteria_reads_nat = "/projects/KmerIPD/test/nativedata.txt"
# pwd_bacteria_reads_wga = "/projects/KmerIPD/test/wgadata.txt"
# bacteria_name = "Geobacter_metallireducens"
# kmer_length = 5
# fasta = "/projects/KmerIPD/test/First_50K.fasta"
# pwd = "/projects/KmerIPD/test/"
#####
  kmer_length=kmer_length_start:kmer_length_end
  source("src/data_formatting.R")
  source("src/KmerList_position_with_Ns.R")
  source("src/unstacking.R")
  source("src/median_per_pos.R")
  source("src/stat_test.R")
  source("src/create_table_each_letter.R")

  #######
  #run only if new csv file
  #######
  pwd_csv_files = pwd_bacteria_reads_nat
  bacteria_name = bacteria_name
  data_formatting(pwd_csv_files,bacteria_name,"nat",pwd,save_intermediary_files)

  pwd_csv_files = pwd_bacteria_reads_wga
  bacteria_name = bacteria_name
  data_formatting(pwd_csv_files,bacteria_name,"wga",pwd,save_intermediary_files)

  #######
  # run only if kmers have not yet been computed for that specific genome
  #######
  for (j in 1:length(kmer_length)){
    KmerList_position(kmer_length[j],pwd,fasta)
  }

  #######
  # run only if unstacking have not yet been computed for that specific bacteria
  #######
  pwd_bacteria_reads_nat = paste(pwd,"/",bacteria_name,"_native_reads.RData",sep = "")
  pwd_bacteria_reads_wga = paste(pwd,"/",bacteria_name,"_wga_reads.RData",sep = "")

  unstacking(pwd,pwd_bacteria_reads_nat,pwd_bacteria_reads_wga,fasta,save_intermediary_files)

  #######
  # run only if medianing have not yet been computed for that specific bacteria
  #######

  median_per_pos(pwd,save_intermediary_files)

  ########
  # run llr test on data
  ########
  pwd_unstacking = paste(pwd,"/unstacked_ipds.RData",sep = "")
  pwd_kmer_files = pwd

  stat_test(kmer_length,pwd,pwd_unstacking,pwd_kmer_files)

  #######
  # run to create table input for Python code
  #######

  create_table_each_letter(kmer_length,pwd)
     
  setwd(pwd)
  # 
  system(paste("cat ",paste('llr_median_',kmer_length,'_each_letter.txt',sep = "", collapse = " "),
    " > catkmer_llr_",bacteria_name,".txt",sep = ""))
}




