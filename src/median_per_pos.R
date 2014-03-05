median_per_pos = function(pwd){
  #################
  # This function computes the median of the sequenced sample for the coverage wanted
  # median_nat is the matrix with column 1 being positions and column 2 being IPDs
  # median_wga is the matrix with column 1 being positions and column 2 being IPDs
  # genome size is the genome size of the genome of interest
  # pwd is the pathway to the unstacked IPD files of interest
  # Positive strand is in column 1, negative strand is in column 2
  ##################
  load(paste(pwd,"/unstacked_ipds.RData",sep=""))

  median_nat_tmp = matrix(data = NA,nrow = genome_size,ncol = 2)
  index = as.numeric(names(nat_reads_tmp))
  index_pos = abs(index[index<0])
  index_neg = index[index>0]
  medians_pos = sapply(seq(1,length(index_pos)),function(x){
    median(nat_reads_tmp[[x]])
  })
  medians_neg = sapply(seq(1+length(index_pos),length(index)),function(x){
    median(nat_reads_tmp[[x]])
  })
  median_nat_tmp[index_pos,1] = medians_pos
  median_nat_tmp[index_neg,2] = medians_neg
  median_nat = median_nat_tmp
  colnames(median_nat) = c("positive_strand","negative_strand")

  median_wga_tmp = matrix(data = NA,nrow = genome_size,ncol = 2)
  index = as.numeric(names(wga_reads_tmp))
  index_pos = abs(index[index<0])
  index_neg = index[index>0]
  medians_pos = sapply(seq(1,length(index_pos)),function(x){
    median(wga_reads_tmp[[x]])
  })
  medians_neg = sapply(seq(1+length(index_pos),length(index)),function(x){
    median(wga_reads_tmp[[x]])
  })
  median_wga_tmp[index_pos,1] = medians_pos
  median_wga_tmp[index_neg,2] = medians_neg
  median_wga = median_wga_tmp
  colnames(median_wga) = c("positive_strand","negative_strand")

  save(median_nat,median_wga,genome_size,file = paste("median_of_reads_per_position.RData",sep = ""))
}
