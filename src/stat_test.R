stat_test = function(length_kmer_vector,pwd,pwd_bacteria_reads,pwd_kmer_files,save_intermediary_files){
##################
# This function calculates the llr on the distribution of IPDs per kmer
# length_kmer_vector is the length of the all the kmers to consider
# pwd is the location in which all the output files will be saved
# pwd_bacteria_reads is the location and file name of the bacteria reads after formatting from the csv file
# pwd_kmer_files is the location of the files of the kmer position per bacteria
# save_intermediary_files is a boolean that specifies whether to save the intermediary files or not
##################

  library(permute)
  setwd(pwd)

  llr_test = function(md_ipd_nat,md_ipd_wga){
    ###################
    # This function calculates the llr test value for 2 distributions iof IPDs
    # md_ipd_nat are the median IPD values of the native sample
    # md_ipd_wga are the median IPD values of the WGA sample
    # If the median IPD of native is larger than the median IPD of WGA, the llr will be negative
    ####################

    md_log_nats = log(md_ipd_nat)
    md_log_wgas = log(md_ipd_wga)
    md_combine = c(md_log_nats,md_log_wgas)
    md_lwgas = length(md_log_wgas)
    md_lnats = length(md_log_nats)
    md_mean_wgas=mean(md_log_wgas)
    md_mean_nats=mean(md_log_nats)
    md_A = sum((md_combine-mean(md_combine))^2)
    md_B = sum((md_log_wgas - md_mean_wgas)^2)
    md_C = sum((md_log_nats - md_mean_nats)^2)
    md_W = ((md_lwgas-1)*var(md_log_wgas)+(md_lnats-1)*var(md_log_nats))/(md_lwgas+md_lnats-2)
    llr_value = (md_A-md_B-md_C)/md_W
    if (md_mean_nats > md_mean_wgas){llr_value =- llr_value}
    llr_value
  }

  median_values_per_kmer = function(median_tmp){
    ############
    # This function extracts the median IPD values at each occurence of one kmer letter, 
    # removes the NA (absent values at kmer letter position), 
    # adds 0.01 at the IPD values equal to 0 
    # and remove the Inf values from the IPD data
    # median_tmp is a vector with all the IPD values for a specific kmer letter accross the entire genome
    ############

    mediantmp = median_tmp[which(is.na(median_tmp) == FALSE)]
    md_ipd = mediantmp
    md_ipd[md_ipd == 0] = md_ipd[md_ipd == 0]+0.01
    md_ipd = md_ipd[which(abs(md_ipd) != Inf)]
  }

  if (save_intermediary_files){
    load(paste("median_of_reads_per_position.RData",sep = ""))
  }
  
  for (count_kmer in 1:length(length_kmer_vector)){
    length_kmer = length_kmer_vector[count_kmer]
    load(paste(pwd_kmer_files,"/KmerList_position_",as.character(length_kmer),".RData",sep = ""))

    kmer_names = names(position_per_kmer) 
    kmer_names_length = length(kmer_names) 
    llr_median_test = list()
    median_mat = matrix(data = list(),nrow = kmer_names_length,ncol = 3)

    for (i in 1:kmer_names_length) {
      lkmer = nchar(kmer_names[i])
      md_tmp_for_second_loop_llr = matrix(data = NA,nrow = lkmer,ncol = 1)
      letters = unlist(strsplit(kmer_names[i],''))
      positionsFirstLetter = as.numeric(position_per_kmer[[i]])
      positionsFirstLetterPos = abs(positionsFirstLetter[positionsFirstLetter<0])
      positionsFirstLetterNeg = positionsFirstLetter[positionsFirstLetter>0]
      for (j in 1:lkmer){
        median_ipd_nat_for_graph = c()
        median_ipd_wga_for_graph = c()
        positionsPos = positionsFirstLetterPos+j-1
        positionsNeg = positionsFirstLetterNeg-j+1
        if (letters[j] == '2'){
          median_tmp_pos = median_nat[positionsPos,1]
          median_tmp_neg = median_nat[positionsNeg,2]
          median_tmp = c(median_tmp_pos,median_tmp_neg)
          md_ipd_nat_A = median_values_per_kmer(median_tmp)
          median_tmp_pos = median_wga[positionsPos,1]
          median_tmp_neg = median_wga[positionsNeg,2]
          median_tmp = c(median_tmp_pos,median_tmp_neg)
          md_ipd_wga_A = median_values_per_kmer(median_tmp)
          if (length(md_ipd_nat_A)>3 & length(md_ipd_wga_A)>3){
            md_tmp_for_second_loop_llr[j] = llr_test(md_ipd_nat_A,md_ipd_wga_A)
            median_ipd_nat_for_graph = c(median_ipd_nat_for_graph,mean(md_ipd_nat_A))
            median_ipd_wga_for_graph = c(median_ipd_wga_for_graph,mean(md_ipd_wga_A))
          }else{
            md_tmp_for_second_loop_llr[j] = NA
            median_ipd_nat_for_graph = c(median_ipd_nat_for_graph,NA)
            median_ipd_wga_for_graph = c(median_ipd_wga_for_graph,NA)
          }
        }
        if (letters[j] == '0'){
          median_tmp_pos = median_nat[positionsPos,1]
          median_tmp_neg = median_nat[positionsNeg,2]
          median_tmp = c(median_tmp_pos,median_tmp_neg)
          md_ipd_nat_C = median_values_per_kmer(median_tmp)
          median_tmp_pos = median_wga[positionsPos,1]
          median_tmp_neg = median_wga[positionsNeg,2]
          median_tmp = c(median_tmp_pos,median_tmp_neg)
          md_ipd_wga_C = median_values_per_kmer(median_tmp)
          if (length(md_ipd_nat_C)>3 & length(md_ipd_wga_C)>3){
            md_tmp_for_second_loop_llr[j] = llr_test(md_ipd_nat_C,md_ipd_wga_C)
            median_ipd_nat_for_graph = c(median_ipd_nat_for_graph,mean(md_ipd_nat_C))
            median_ipd_wga_for_graph = c(median_ipd_wga_for_graph,mean(md_ipd_wga_C))
          }else{
            md_tmp_for_second_loop_llr[j] = NA
            median_ipd_nat_for_graph = c(median_ipd_nat_for_graph,NA)
            median_ipd_wga_for_graph = c(median_ipd_wga_for_graph,NA)
          }
        }
        if (letters[j] == '1'){
          median_tmp_pos = median_nat[positionsPos,1]
          median_tmp_neg = median_nat[positionsNeg,2]
          median_tmp = c(median_tmp_pos,median_tmp_neg)
          md_ipd_nat_T = median_values_per_kmer(median_tmp)
          median_tmp_pos = median_wga[positionsPos,1]
          median_tmp_neg = median_wga[positionsNeg,2]
          median_tmp = c(median_tmp_pos,median_tmp_neg)
          md_ipd_wga_T = median_values_per_kmer(median_tmp)
          if (length(md_ipd_nat_T)>3 & length(md_ipd_wga_T)>3){
            md_tmp_for_second_loop_llr[j] = llr_test(md_ipd_nat_T,md_ipd_wga_T)
            median_ipd_nat_for_graph = c(median_ipd_nat_for_graph,mean(md_ipd_nat_T))
            median_ipd_wga_for_graph = c(median_ipd_wga_for_graph,mean(md_ipd_wga_T))
          }else{
            md_tmp_for_second_loop_llr[j] = NA
            median_ipd_nat_for_graph = c(median_ipd_nat_for_graph,NA)
            median_ipd_wga_for_graph = c(median_ipd_wga_for_graph,NA)
          }
        }
        if (letters[j] == '3'){
          median_tmp_pos = median_nat[positionsPos,1]
          median_tmp_neg = median_nat[positionsNeg,2]
          median_tmp = c(median_tmp_pos,median_tmp_neg)
          md_ipd_nat_G = median_values_per_kmer(median_tmp)
          median_tmp_pos = median_wga[positionsPos,1]
          median_tmp_neg = median_wga[positionsNeg,2]
          median_tmp = c(median_tmp_pos,median_tmp_neg)
          md_ipd_wga_G = median_values_per_kmer(median_tmp)
          if (length(md_ipd_nat_G)>3 & length(md_ipd_wga_G)>3){
            md_tmp_for_second_loop_llr[j] = llr_test(md_ipd_nat_G,md_ipd_wga_G)
            median_ipd_nat_for_graph = c(median_ipd_nat_for_graph,mean(md_ipd_nat_G))
            median_ipd_wga_for_graph = c(median_ipd_wga_for_graph,mean(md_ipd_wga_G))
          }else{
            md_tmp_for_second_loop_llr[j] = NA
            median_ipd_nat_for_graph = c(median_ipd_nat_for_graph,NA)
            median_ipd_wga_for_graph = c(median_ipd_wga_for_graph,NA)
          }
        }
      }
      z = unlist(strsplit(kmer_names[i],""))
      z[which(z == "3")] = "G"
      z[which(z == "2")] = "A"
      z[which(z == "1")] = "T"
      z[which(z == "0")] = "C"
      z = paste(z,collapse = "")
      median_mat[i,1] = z
      median_mat[i,2] = median_ipd_wga_for_graph
      median_mat[i,3] = median_ipd_nat_for_graph
      llr_median_test[[z]] = md_tmp_for_second_loop_llr
    }

    save(llr_median_test,file = paste("llr_values_",length_kmer,".RData",sep = ""))
    max_median_llr = llr_median_test

    for (i in 1:length(llr_median_test)){
      max_median_llr[i] = max(abs(na.omit(unlist(llr_median_test[i]))))
    }

    sort_max_median_llr = unlist(max_median_llr)
    ind = is.finite(sort_max_median_llr)
    sort_max_median_llr = sort(abs(sort_max_median_llr[ind]),decreasing = FALSE)

    if (save_intermediary_files){
      write.table(sort_max_median_llr, file = paste0("llr_median_",length_kmer,".txt"))
      
      save(median_mat,file = paste0("median_of_medians_",length_kmer,".RData"))
    }
  }
}
