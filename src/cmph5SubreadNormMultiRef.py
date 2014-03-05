import os,sys
from pbcore.io.CmpH5Reader import *
from pbcore.io.FastaIO import FastaReader 
import numpy as np

def writeLinesFromCmph5 (cmph5, leftAnchor, rightAnchor, offsetDict):
     reader           = CmpH5Reader(cmph5)
     alignments_list  = [r for r in reader]
     #refInfoTable = reader.referenceInfoTable
     #refDict = {}
     #for i in range (len(refInfoTable)):
     #     rid = refInfoTable[i][0]
     #     rn = refInfoTable[i][2]
     #     rname = refInfoTable[i][3]
     #     refDict[rn] = rname
          #print refInfoTable[i]
          
     for i, alignment in enumerate(alignments_list):

          #movieID       = str(alignment.movieInfo[0])
          alignedLength = alignment.alignedLength
          fps           = alignment.movieInfo[2]
          #refName       = alignment.referenceName
          #refName = refDict[refName]
          refName       = str(alignment.referenceInfo[3])
          #refGroupID    = alignment.refGroupID
          #refName = refDict[refGroupID]
          #zmw           = str(alignment.HoleNumber)
          #mol           = str(alignment.MoleculeID)
          if alignment.isForwardStrand:
              strand = str(0)
          else:
              strand = str(1)
          ref_bases  = alignment.reference()
          read_calls = alignment.transcript()
          ref_pos    = list(alignment.referencePositions())
          IPD        = list(alignment.IPD())

          delim           = " "

          error_mk = []
          for read_call in read_calls:
              # Go through all entries and flag which positions are MM/indels
              if read_call != "M":
                  # Mismatch or indel at this position!
                  error_mk.append(1)
              else:
                  error_mk.append(0)

          # Get the indices of all the non-matches
          error_idx = [i for (i,val) in enumerate(error_mk) if val == 1]
          for error_id in error_idx:
              try:
                  for j in range(leftAnchor):
                      error_mk[error_id - (j+1)] = 1
                  for j in range(rightAnchor):
                      error_mk[error_id + (j+1)] = 1
              except IndexError:
                  pass
          error_mk = np.array(error_mk)

          ipds       = np.array(IPD) / fps
          strands    = np.array([strand]     * alignedLength)

          ref_bases  = np.array(list(ref_bases))
          ref_pos    = np.array(ref_pos)
          read_calls = np.array(list(read_calls))

          ref_bases  =  ref_bases[error_mk==0]
          ref_pos    =    ref_pos[error_mk==0]
          read_calls = np.array(read_calls)[error_mk==0]
          ipds       =       ipds[error_mk==0]
          strands    =    strands[error_mk==0]
          ipds = ipds/np.median(ipds)
          for i in range (ipds.size):
              newpos = ref_pos[i] + offsetDict[refName]
              print newpos, ipds[i], strand


"python cmph5SubreadNormMultiRef.py comb_50k.fa comb_50k_merged.fa out_csal_wga.cmp.h5,out_wga_final.cmp.h5  7 2 > test_input.txt"

fastaFile            = sys.argv[1]
fastaOut             = open(sys.argv[2], 'w')
cmph5s               = sys.argv[3].split(",")
if len(sys.argv) > 2:
     leftAnchor      = int(sys.argv[4])
     rightAnchor     = int(sys.argv[5])
else:
     leftAnchor      = 1
     rightAnchor     = 1


fastaOut.write(">merged\n")

count = 0
offset = 0
offsetDict = {}
for entry in FastaReader(fastaFile):
     if count != 0:
          fastaOut.write("NNNNNNNNNN")
          offset += 10
     fastaOut.write(entry.sequence)
     offsetDict[entry.name] = offset
     offset += len(entry.sequence)
     count += 1
fastaOut.write("\n")
fastaOut.close()

for cmph5 in cmph5s:
     writeLinesFromCmph5 (cmph5, leftAnchor, rightAnchor, offsetDict)
