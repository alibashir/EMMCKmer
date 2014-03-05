create_table_each_letter = function(length_kmers,pwd){
###############
# This function creates a table with all the scores per kmer per letter for each of the 3 tests done in stat_test
# length_kmer is the length of the kmer of interest
# pwd is the location in which the files will be saved
###############


   for (j in 1:length(length_kmers)){
    	length_kmer = length_kmers[j]
    	setwd(pwd)
    	load(paste("llr_values_",length_kmer,".RData",sep = ""))

		llr_test = unlist(llr_median_test)
		ind = is.finite(llr_test)
		llr_test = llr_test[ind]

		write.table(llr_test, file = paste0('llr_median_',length_kmer,'_each_letter.txt'),col.names=F)
	} 
}