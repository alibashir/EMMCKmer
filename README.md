EMMCKmer
========

**E**pigenetic **M**etagenomic **M**otif **C**haracterization by **Kmer** analysis

## Overview
This tool compares the distribution of IPDs between NATIVE and WGA data for every existing kmer in the genome (using log-likelihood rations) and then uses an iterative greedy procedure to identify significant motifs.  The output of this tool is a list of the significantly modified kmers.

##Running
To run the code simply type:

<code>sh EMMCkmer.sh [wga.cmp.h5] [nat.cmp.h5] [sample id] [fasta file] [output directory] </code>

### INPUTS

This wrapper script takes as input the cmp.h5 files for both NATIVE and WGA data, the SAMPLE name, the FASTA file of the reference genome and the OUTPUT directory.  

### INSTALLATION/REQUIREMENTS
You should not need to install anything if you have satisifed the following requirements:

<ul>
<li>SMRTANALYSIS package distributed by PacBio </li>
<li>PYTHON (include in smrtanalysis) </li>
<li>R (included in smrtanalysis) </li>
</ul>
 
## Citation
To cite the code, please refer to ...
