#!/bin/bash
############
# BASH script mastercode:
############
WGA=$1
NAT=$2
SAMPLE=$3
FASTA=$4
OUTDIR=$5
FASTAMERGE=$SAMPLE.merged.fa

source  /etc/profile.d/modules.bash


echo "python src/cmph5SubreadNormMultiRef.py $FASTA $OUTDIR/$FASTAMERGE $WGA 7 2 > $OUTDIR/wgadata.txt"
python cmph5SubreadNormMultiRef.py $FASTA $OUTDIR/$FASTAMERGE $WGA 7 2 > $OUTDIR/wgadata.txt

echo "python src/cmph5SubreadNormMultiRef.py $FASTA $OUTDIR/$FASTAMERGE $NAT 7 2 > $OUTDIR/nativedata.txt"
python cmph5SubreadNormMultiRef.py $FASTA $OUTDIR/$FASTAMERGE $NAT 7 2 > $OUTDIR/nativedata.txt

echo "Rscript src/MastercodeRPipeline.R nativedata.txt wgadata.txt $SAMPLE 4 7 $FASTA $OUTDIR"
Rscript MastercodeRPipeline.R $OUTDIR/nativedata.txt $OUTDIR/wgadata.txt $SAMPLE 4 7 $FASTAMERGE $OUTDIR/


echo "python src/kmerSignficance.py --scoreThresh=0  --kvpos --simple --nozscore --usePosAsCut=0 --removeNeighbors  --dontRescoreShiftChildren   --useGammaProbAsCut=0.000001 --neighborCut  7.6198530241604696e-24 --out $OUTDIR/significant_kmers_$SAMPLE.txt  $OUTDIR/catkmer_llr_$SAMPLE.txt"
python kmerSignficance.py --scoreThresh=0  --kvpos --simple --nozscore --usePosAsCut=0 --removeNeighbors  --dontRescoreShiftChildren   --useGammaProbAsCut=0.000001 --neighborCut  7.6198530241604696e-24 --out $OUTDIR/significant_kmers_$SAMPLE.txt  $OUTDIR/catkmer_llr_$SAMPLE.txt


