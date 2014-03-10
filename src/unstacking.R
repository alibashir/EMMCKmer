unstacking = function(pwd,pwd_bacteria_reads_nat,pwd_bacteria_reads_wga,fasta,save_intermediary_files){
#################
# This function is unstacking the IPD file
# pwd is the location in which we want to save our output
# pwd_bacteria_reads_nat is the location of the native IPDs RData file
# pwd_bacteria_reads_wga is the location of the WGA IPDs RData file
# fasta is the location of the fasta file
#################

	setwd(pwd)
	if (save_intermediary_files){
		load(pwd_bacteria_reads_nat)
		load(pwd_bacteria_reads_wga)
	}
	
	ref = read.fasta(fasta)
	genome_size = getLength(ref)
	wga_reads_tmp = unstack(wga_reads[2:1])
	nat_reads_tmp = unstack(nat_reads[2:1])
	if (save_intermediary_files){
		save(nat_reads_tmp,wga_reads_tmp,genome_size,file = paste("unstacked_ipds.RData",sep = ""))
	}
}
