data_formatting = function(pwd_csv_files,bacteria_name,nat_wga,pwd){
#############
# This function takes as input the ipd per position csv file after Python code
# The output is 2 RData files containing each a matrix of the positions and corresponding IPDs 
# pwd_csv_files is the location of the folder and IPD file
# Bacteria_name is the bacteria name of the bacteria of interest
# nat_wga is the description of native versus wga data (needs to be equal to "nat" or "wga")
# Positive strand is negative values of positions, negative strand is positive values of positions
# pwd is the location in which to save the output files
#############
	setwd(pwd)
	t = read.table(pwd_csv_files)
	ind = which(t[,3] == 1)
	t[,1] = t[,1]+1
	t[ind,1] = -t[ind,1]

	if (nat_wga == "nat"){
		nat_reads = t[,1:2]
	save(nat_reads,file = paste(bacteria_name,"_native_reads.RData",sep = ""))
	}
	if (nat_wga == "wga"){
		wga_reads = t[,1:2]
	save(wga_reads,file = paste(bacteria_name,"_wga_reads.RData",sep = ""))
	}
}