#!/usr/bin/env python
""" Get the most signficant IPDs """

import sys, os
import optparse, logging
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import scipy.stats as ss

class KmerIPDGraph:
    def __init__( self ):
        self.__parseArgs( )
        self.__initLog( )


    def __parseArgs( self ):
        """Handle command line argument parsing"""
        
        usage = "%prog [--help] [options] <catted_LLR_kmers_for_each_letter>"
        parser = optparse.OptionParser( usage=usage, description=__doc__ )

        parser.add_option( "-l", "--logFile", help="Specify a file to log to. Defaults to stderr." )
        parser.add_option( "-d", "--debug", action="store_true", help="Increases verbosity of logging" )
        parser.add_option( "-i", "--info", action="store_true", help="Display informative log entries" )
        parser.add_option( "-p", "--profile", action="store_true", help="Profile this script, dumping to <scriptname>.profile" )
        parser.add_option( "-s", "--scoreThresh", help="score threshold to filter data (DEFAULT = 100)" )
        parser.add_option( "-t", "--top", help="take the top X hits at each kmer size")
        parser.add_option( "--topPercent", help="take the top % hits at each kmer size")
        parser.add_option( "--pval", action="store_true", help="handle values as pvals")        
        parser.add_option( "--kvpos", action="store_true", help="includes the position in kmer (REQUIRED IF YOU PASS THE KMER FILE WITH POS)")
        parser.add_option( "--simple", action="store_true", help="simple motif set method instead of graph")        
        parser.add_option( "--hitsToReturn", help="The total umber of motifs to return in simple mode")        
        parser.add_option( "--writePngs", action="store_true", help="write out local kmer distribion to png")        
        parser.add_option( "--nozscore", action="store_true", help="just do llr subtraction")
        parser.add_option( "--removeNeighbors", action="store_true", help="dont rescore parents and childrens just eliminate")
        parser.add_option( "--usePosAsCut", help="float to set as multiple for cutoff")
        parser.add_option( "--useMadAsCut", help="float to set as multiple for MAD cutoff")
        parser.add_option( "--useGammaProbAsCut", help="float to set as gamma 1-cdf cutoff")

        parser.add_option("--dontRescoreShiftChildren", action="store_true", help="don't rescore children of prefix/suffix shifts of a true motif")
        parser.add_option("--neighborCut", help="don't rescore children of prefix/suffix shifts of a true motif")
        parser.add_option("--out", help="name of outfile")

        parser.set_defaults( logFile=None, debug=False, info=False, profile=False, scoreThresh=100, top=False, 
                            topPercent=False, kvpos=False, pval=False, simple=False, hitsToReturn=1000, writePngs=False,
                            removeNeighbors=False, usePosAsCut=None, useMadAsCut=None, dontRescoreShiftChildren=False,
                            neighborCut=0, useGammaProbAsCut=None, out=None)
        
        self.opts, args = parser.parse_args( )

        if len(args) != 1:
            parser.error( "Expected a single argument." )

        self.inFN = args[0]
        self.scoreThresh = float(self.opts.scoreThresh)
        self.topN = int(self.opts.top)
        self.topP = float(self.opts.topPercent)
        self.kvpos = self.opts.kvpos
        self.maxKmerLen = 7
        self.zcutoff = 20
        # the raw score = llr
        self.sindex = 1 
        # value to reset when we observe LLRs in a neighborhood
        self.zindex = 2
        # check to see if this value has already been rescored as a neighbor
        # always assume that the first (more significant) neighborhood should
        # be used
        self.findex = 3
        # the neighbor rescore
        self.dindex = 4
        # not used by default - but looks at Median Absolute Deviation
        self.mindex = 5
        # the gamma probablity score
        self.gindex = 6
        
        self.minKVal = {}
        for i in range (10):
            self.minKVal[i] = -1.0

        self.hitsToReturn = int(self.opts.hitsToReturn)
        if self.opts.writePngs:
            try:
                os.mkdir("pngs")
            except:
                pass
        self.kmerSet = set()
        self.usePosAsCut = 0
        if self.opts.usePosAsCut != None:
            self.usePosAsCut = float(self.opts.usePosAsCut)
        self.useMadAsCut = 0
        if self.opts.useMadAsCut != None:
            self.useMadAsCut = float(self.opts.useMadAsCut)
        if self.opts.useGammaProbAsCut != None:
            self.useGammaProbAsCut = float(self.opts.useGammaProbAsCut)
                           
        self.neighborCut = float(self.opts.neighborCut)
        
        if self.opts.out:
            self.out = open(self.opts.out, 'w')
        else:
            self.out = sys.stdout

    def __initLog( self ):
        """Sets up logging based on command line arguments. Allows for three levels of logging:
        logging.error( ): always emitted
        logging.info( ): emitted with --info or --debug
        logging.debug( ): only with --debug"""

        logLevel = logging.DEBUG if self.opts.debug else logging.INFO if self.opts.info else logging.ERROR
        logFormat = "%(asctime)s [%(levelname)s] %(message)s"
        if self.opts.logFile != None:
            logging.basicConfig( filename=self.opts.logFile, level=logLevel, format=logFormat )
        else:
            logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
                                                                 
    def run( self ):
        """Executes the body of the script."""
    
        logging.info("Log level set to INFO")
        logging.debug("Log Level set to DEBUG")
        if self.opts.simple:
            logging.info("reading kmer file")
            # motifDict - keys = k, vals = list of motif full infos
            # motifSet - a dict wih keys = motif, vals = motif full info
            motifDict, motifSet = self.readKmerFileToSimpleKmerDict (self.inFN)
            logging.info("findng top hits")
            self.out.write("kmer_tuple llr gamma_prob neighbor_prob\n")
            self.iterativeFindTopHits (motifDict, motifSet)
        else:
            print >>sys.stderr, "Currently simple is the only supported mode"
        return 0


    def readKmerFileToSimpleKmerDict (self, inputFile):
        '''Read the input file to a dictionary where keys 
        are kmer tuples and values correspond to lists with llrs 
        and other scores'''
        motifDict = {}
        motifSet = {}
        maxKmerLen = 0
        for line in open(inputFile).xreadlines():
            ll = line.strip().split(" ")
            if self.kvpos:
                kmer = ll[0][1:-2]
                pos = int(ll[0][-2:-1])
            else:
                kmer = ll[0][1:-1]
                pos = 0
            llr_raw = -1*(float(ll[1]))
            kmerTup = [(kmer,pos), llr_raw, 0, False, 1, 0, 0]
            k = len(kmer)
            if k > self.maxKmerLen:
                continue
            if self.minKVal[k] > llr_raw:
                self.minKVal[k] = llr_raw
                
            motifDict.setdefault(k, []).append(kmerTup)
            motifSet[(kmer,pos)] = kmerTup

        self.simpleKmerDistOverLengths (motifDict)
        logging.info(str(self.minKVal))
        return motifDict, motifSet
            

    def simpleKmerDist ( self, motifs):
        ''' Try to get a simple distribution of significance for 
        kmers within a given set'''
        allvals = map(lambda x: x[1], motifs)
        allvals.sort()
        vlen = len(allvals)

        s, e = 5, len(allvals)-5
        vals = allvals[s:e]
        v_mean = np.mean(vals)
        v_std = np.std(vals)
        logging.info(" mean = %f, stdev = %f" %(v_mean, v_std))
        
        d = np.abs(allvals - np.median(allvals))
        mdev = np.median(d)
        med_allvals = np.median(allvals)
        
        for motif in motifs:
            if self.opts.nozscore:
                motif[2] = motif[1]
            else:
                motif[2] = abs(motif[1]-v_mean)/v_std
            motif[self.mindex] = abs(motif[1] - med_allvals)/ mdev

    

    def simplePlotter (self, motifs, titlePrefix, gammaPlot = False):
        ''' For plotting distribution of kmers within a set to assess
        significance'''
        if len(motifs) == 0:
            return
        vals = map(lambda x: x[1], motifs)
                
        plt.hist(map(lambda x: x[1], motifs), bins=1000, color='b', alpha=.3, log=True)
        ymin, ymax = plt.ylim()
        plt.ylim(0.5, ymax)
        plt.title("%s llr\n(%i total values)" %(titlePrefix, len(motifs)))
        plt.savefig("pngs/%s_llr.png" %(titlePrefix))
        plt.clf()
                
        plt.hist(map(lambda x: x[2], motifs), color='r', alpha=.3, log=True)
        plt.title("%s z-score\n(%i total values)" %(titlePrefix, len(motifs)))
        ymin, ymax = plt.ylim()
        plt.ylim(0.5, ymax)
        plt.savefig("pngs/%s_zscore.png" %(titlePrefix))
        plt.clf()
        
    


    def simpleKmerDistOverLengths (self, motifDict):
        for k, motifs in motifDict.items():
            logging.info("For kmer of size %i:" %(k))
            self.simpleKmerDist(motifs)
            # perform gamma corretion
            vals = map(lambda x: x[1], motifs)
            min_val = min(vals)
            pos_vals = filter(lambda x: x>0, vals)
            pos_vals.sort()
            neg_vals = filter(lambda x: x<0, vals)
            neg_vals = map(lambda x: -1*x, neg_vals)
            pos_vals.sort()

            cutoff = max(-1*min_val, pos_vals[int(.99*len(pos_vals))-1])
            pos_vals_sub = filter(lambda x: x<(cutoff), pos_vals)
            fit_alpha,fit_loc,fit_beta=ss.gamma.fit(pos_vals_sub)
            rv = ss.gamma(fit_alpha, fit_loc, fit_beta)

            for motif in motifs:
                motif[self.gindex] = rv.sf(motif[1])

            if self.opts.writePngs:
                self.simplePlotter (motifs, "All motifs of size %i" %(k), gammaPlot=True)
 
    def rescoreSet (self, motifs):
        if len(motifs) == 0:
            return
        vals = map(lambda x: x[self.sindex], motifs)        
        if len(motifs) > 2:            
            vals.sort()

            mid = np.mean(vals)
            allowedDiff = mid-min(vals)
            subsetvals = filter(lambda x: abs(x-mid) <= allowedDiff, vals)
        
            correction = np.mean(subsetvals)
            stdev = np.std(subsetvals)
        else:

            correction = np.mean(vals)
            stdev = np.std(vals)
        logging.info("motif len: %i, motif count: %i, mean: %f, stdev: %f" %(len(motifs[0][0][0]), len(motifs), correction, stdev))

        d = np.abs(vals - np.median(vals))
        mdev = np.median(d)
        med_vals = np.median(vals)

        for motif in motifs:
            if self.opts.removeNeighbors:
                motif[self.mindex] = 0
                motif[self.zindex] = 0
                if motif[self.dindex] == 1:
                    std =  abs(motif[self.sindex]-correction)/stdev
                    motif[self.dindex] = (len(motifs)*ss.norm.sf(std))

                    
            elif self.opts.nozscore:

                if not motif[self.findex]:
                    motif[self.zindex] = motif[self.sindex]-correction
                    motif[self.findex] = True
            else:
                motif[self.zindex] = abs(motif[self.sindex]-correction)/stdev


    def rescoreKmerNeighbors ( self, kmerTup, motifSet ):
        kmer, modpos = kmerTup        
        logging.info("rescoring (%s, %i)" %(kmer,modpos))
        kmerPre = kmer[0:-1]
        kmerSuf  = kmer[1:]        
        logging.info("kmerPre: %s; kmerSuf: %s" %(kmerPre, kmerSuf))
        motifSubset = []
        if (kmerPre, modpos) in motifSet:
            motifSubset.append(motifSet[(kmerPre, modpos)])
        if (kmerSuf, modpos-1) in motifSet:
            motifSubset.append(motifSet[(kmerSuf, modpos-1)])
        self.rescoreSet(motifSubset)        

        motifSubset = []
        for nuc in "ACGT":
            kmer1 = nuc+kmerPre
            kmer2 = kmerSuf+nuc
            if kmer1 != kmer:
                if (kmer1, modpos+1) in motifSet:
                    motifSubset.append(motifSet[(kmer1, modpos+1)])
            if kmer2 != kmer:
                if (kmer2, modpos-1) in motifSet:
                    motifSubset.append(motifSet[(kmer2, modpos-1)])
        self.rescoreSet(motifSubset)        
        
        motifKDictSet = {}
        self.recursiveRescoreByK(motifSet, motifKDictSet, kmerTup)
        self.rescoreKDictSet(motifKDictSet, kmerTup)
        motifKDictSet = {}
        preTup = (kmerPre, modpos)
        sufTup = (kmerSuf, modpos-1)
        if not self.opts.dontRescoreShiftChildren:
            logging.info(" for rescoring, kmerpre: %s, kmersuf: %s" %(preTup, sufTup))
            self.recursiveRescoreByK(motifSet, motifKDictSet, preTup)
            self.rescoreKDictSet(motifKDictSet, preTup)
            motifKDictSet = {}
            self.recursiveRescoreByK(motifSet, motifKDictSet, sufTup)
            self.rescoreKDictSet(motifKDictSet, sufTup)



    def removeKmerNeighbors ( self, kmerTup, motifSet ):
        kmer, modpos = kmerTup        
        logging.info("rescoring (%s, %i)" %(kmer,modpos))
        kmerPre = kmer[0:-1]
        kmerSuf  = kmer[1:]        
        motifSubset = []
        if (kmerPre, modpos) in motifSet:
            motifSubset.append(motifSet[(kmerPre, modpos)])
        if (kmerSuf, modpos-1) in motifSet:
            motifSubset.append(motifSet[(kmerSuf, modpos-1)])
        self.rescoreSet(motifSubset)        

        motifSubset = []
        for nuc in "ACGT":
            kmer1 = nuc+kmerSuf
            kmer2 = kmerPre+nuc
            if kmer1 != kmer:
                if (kmer1, modpos+1) in motifSet:
                    motifSubset.append(motifSet[(kmer1, modpos+1)])
            if kmer2 != kmer:
                if (kmer2, modpos) in motifSet:
                    motifSubset.append(motifSet[(kmer2, modpos)])
        self.rescoreSet(motifSubset)        
        
        motifKDictSet = {}
        self.recursiveRescoreByK(motifSet, motifKDictSet, kmerTup)
        self.rescoreKDictSet(motifKDictSet, kmerTup)
        
    def rescoreKDictSet (self, motifKDictSet, kmerTup):
        for k in motifKDictSet:
            motifSubset = motifKDictSet[k].values()
            if self.opts.writePngs:
                self.simplePlotter(motifSubset, "(%s, %i) containing motifs of size %i" %(kmerTup[0], kmerTup[1], k))            
            self.rescoreSet(motifSubset)

    def recursiveRescoreByK (self, motifSet, motifKDictSet, kmerTup):
        kmer, modpos = kmerTup
        k = len(kmer)
        motifSubset = motifKDictSet.setdefault(k+1, {})
        
        for nuc in "ACGT":
            kmer1 = nuc+kmer
            kmer2 = kmer+nuc
            kmerTup1 = (kmer1, modpos+1) 
            if kmerTup1 in motifSet:
                motifSubset[kmerTup1] = motifSet[kmerTup1]
            kmerTup2 = (kmer2, modpos)                 
            if kmerTup2 in motifSet:
                motifSubset[kmerTup2] = motifSet[kmerTup2]
                    
        if k < self.maxKmerLen:
            for nuc in "ACGT":
                self.recursiveRescoreByK(motifSet, motifKDictSet, (nuc+kmer, modpos+1))
                self.recursiveRescoreByK(motifSet, motifKDictSet, (kmer+nuc, modpos))                
            
        return
        
    def recursiveRescore (self, motifSet, kmerTup):
        motifSubset = []
        kmer, modpos = kmerTup
        for nuc in "ACGT":
            kmer1 = nuc+kmer
            kmer2 = kmer+nuc
            if kmer1 != kmer:
                if (kmer1, modpos+1) in motifSet:
                    motifSubset.append(motifSet[(kmer1, modpos+1)])
            if kmer2 != kmer:
                if (kmer2, modpos) in motifSet:
                    motifSubset.append(motifSet[(kmer2, modpos)])
        self.rescoreSet(motifSubset)        
        if len(kmer) < self.maxKmerLen:
            for nuc in "ACGT":
                self.recursiveRescore(motifSet, (nuc+kmer, modpos+1))
                self.recursiveRescore(motifSet, (kmer+nuc, modpos))                
        return
        

    def iterativeFindTopHits (self, motifDict, motifSet, count = 0, count2 = 0, sigMotifSet=set()):
        if count >= self.hitsToReturn:
            if count2 >= self.hitsToReturn:
                return
            #if count2 == 0:
                #self.out.write("neighbor hits\n")
                #self.out.write("kmer_tuple llr gamma_prob neighbor_prob\n")
            maxz = 1
            maxmotif = None
            for k in motifDict:
                if k < 4 or k > self.maxKmerLen:
                    continue
                motifs = motifDict[k]
                for i in range (len(motifs)):
                    motif = motifs[i]
                    z = motif[self.sindex]
                    
                    if self.opts.useMadAsCut != None:
                        if self.useMadAsCut < motif[self.mindex]:
                            pass
                        else:
                            continue               
                    if self.opts.usePosAsCut != None:
                        if -self.usePosAsCut*self.minKVal[k] < z:
                            pass
                        else:
                            continue
                    if self.opts.useGammaProbAsCut != None:
                        if self.useGammaProbAsCut/(float(4.0*k**4)) > motif[self.gindex]:
                            pass
                        else:
                            continue
                    neighbor_val = motif[self.dindex]
                    if  neighbor_val < maxz:
                        maxmotif = motifs[i]
                        maxz = neighbor_val
            if self.neighborCut > maxz:
                logging.info("found neighbor new motif, with zscore: %f" %(maxz))
                maxmotif[self.dindex] = 1
                if self.checkIfScoreHigherThanParents (maxmotif[0], maxmotif[self.sindex], motifSet, sigMotifSet, printFlag=False):
                    if maxmotif[0] in sigMotifSet:
                        pass
                    else:
                        self.out.write("%s %f %s %s\n" %(maxmotif[0], maxmotif[self.sindex], str(maxmotif[self.gindex]), str( maxz)))  
                    sigMotifSet.add(maxmotif[0])
                    return self.iterativeFindTopHits(motifDict, motifSet, count, count2+1, sigMotifSet)

                else:
                    return self.iterativeFindTopHits(motifDict, motifSet, count, count2, sigMotifSet)
            return
        maxz = 0
        maxmotif = None
        for k in motifDict:
            if k < 4 or k > self.maxKmerLen:
                continue
            motifs = motifDict[k]

            for i in range (len(motifs)):
                motif = motifs[i]
                z = motif[self.zindex]
                
                if self.opts.useMadAsCut != None:
                    if self.useMadAsCut < motif[self.mindex]:
                        pass
                    else:
                        continue               
                if self.opts.usePosAsCut != None:
                    if -self.usePosAsCut*self.minKVal[k] < z:
                        pass
                    else:
                        continue
                if self.opts.useGammaProbAsCut != None:
                    if self.useGammaProbAsCut/(float(4.0*k**4)) > motif[self.gindex]:
                        pass
                    else:
                        continue

                if z > maxz:
                    maxmotif = motifs[i]
                    maxz = z
                    
        if maxz > self.scoreThresh:
            logging.info("found new motif, with zscore: %f" %(maxz))
            
            if self.checkIfScoreHigherThanParents (maxmotif[0], maxz, motifSet):
                self.out.write("%s %f %s N\A\n" %(maxmotif[0], maxmotif[self.sindex],str(maxmotif[self.gindex]))  )
                self.rescoreKmerNeighbors ( maxmotif[0], motifSet )
                newcount = count + 1
                sigMotifSet.add(maxmotif[0])
            else:
                newcount = count
            maxmotif[self.zindex] = 0

            return self.iterativeFindTopHits(motifDict, motifSet, newcount, 0, sigMotifSet)
                
        #rescore subet
        maxz = 0
        maxmotif = None
        return self.iterativeFindTopHits(motifDict, motifSet, self.hitsToReturn, 0, sigMotifSet)

    def checkIfScoreHigherThanParents (self, motif, mz, motifSet, sigMotifSet=set(), printFlag=False):
        for i in range ((len(motif[0])-4)+1):
            # make sure the modified position in our motif
            if i >= motif[1]:
                continue
            # select a k from allowed range (takinginto account remaining motif length)
            for k in range(4, min(self.maxKmerLen, len(motif[0])-i)+1):
                # create a new valid parent motif
                newmotif= (motif[0][i:i+k], motif[1]-i)

                # make sure the modified position in our motif
                if motif[1]-i > k:
                    continue

                # make sure the new motif isn't in the already signifcant motif set
                if newmotif in sigMotifSet:
                    continue

                # check if current parent motif score is greatr than observed motif
                if newmotif in motifSet and motifSet[newmotif][self.sindex] > mz:
                    if printFlag:
                        print "found a new motif %s with higher score than the original %s" %(str(newmotif), str(motif))
                    return False
        return True

if __name__ == "__main__":
    app = KmerIPDGraph()
    if app.opts.profile:
        import cProfile
        cProfile.run( 'app.run()', '%s.profile' % sys.argv[0] )
    sys.exit( app.run() )
