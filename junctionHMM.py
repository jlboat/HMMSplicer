#! /usr/bin/env python
#
#    Copyright 2009,2010 Michelle Dimon
#
#    This file is part of HMMSplicer.
#
#    HMMSplicer is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    HMMSplicer is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with HMMSplicer.  If not, see <http://www.gnu.org/licenses/>.
#
# junctionHMM.py
#
# This module contains higher level functions used by HMMSplicer.  These functions accomplish the
# major pieces of functionality involved in HMMSplicer, for example seeding the reads in the genome
# based on the bowtie results, training the HMM, running the HMM, and finding the location of the
# second part of the read.
#
# Author: Michelle Dimon, May 2009
#

from numpy import array, zeros, argmax
import hmmWithQuality
import hmmUtils
import configVals
import hmmErrors
import os

OBSERVE_LIST = ["M", "X"]
STATE_LIST = ["a", "j"]

def seedReads(bowtieSeed, genomeDict, outputSeedFile, solexaFormat=False, colorspaceFormat=False):
    """Given bowtie seed results (from half the read) figures out the right seed location for each
    read and output the results in a tab delimited file, one line per read, with these fields:
    -- name
    -- seq
    -- quality string
    -- genome seq (aligns with the read seq)
    -- chr
    -- tStart
    -- tStop
    -- which half matched (1 or 2)
    -- score of the seed alignment
    -- score of the non-aligning seed
    
    The input GenomeDict is a dictionary of the genome so the correct sequence can be looked up.
    """
    
    if solexaFormat and colorspaceFormat:
        raise hmmErrors.InvalidInputException("Reads cannot be both in Solexa format and ColorSpace format")
    
    out = open(outputSeedFile, "w")
    count = 0
    countGood = 0
    
    for line in open(bowtieSeed):
        count += 1
        #if count > 10:
        #    break

        try:
            [name, strand, chr, chrStart, seq, qual, x, mm] = line.split("\t")
            chrStart = int(chrStart)
        except ValueError:
            print line
            print line.split()
            print len(line.split())
            raise hmmErrors.InvalidFileFormatException("Invalid line in junctionHMM.seedReads")
        
        numWritten = _processOneBowtieRead(seq, qual, chr, chrStart, strand, genomeDict, name, out, colorspaceFormat)
        countGood += numWritten
        
    #print "Done!  Out of %s reads, there were %s good lines" % (count, countGood)
    
    return countGood


def _processOneBowtieRead(seq, qual, chr, chrStart, strand, genomeDict, name, out, colorspaceFormat):
    """Writes one read-half alignment to the output in dbf (DeRisi Bed file) format.
    """
    
    #genomeSeq = genomeDict[chr][chrStart:chrStart+len(seq)]
    
    # this happens when the read is at the very beginning or very end of the chromosome.  We can ignore these.
    #if len(seq) != len(genomeSeq):
    #    print "Read sequence and genome sequence are different lengths!"
    #    print "Strand: %s, Seq: %s, name: %s, chr: %s, chrStart: %s" % (strand, seq, name, chr, chrStart)
    #    print "----------"
    #    return 0
    
    # from the name, parse the half, the rest of the seq and the rest of the quals
    [title, readHalf, seqHalf, qualHalf] = name.split("|")
    #print 
    #print title
    if strand=="+" and readHalf=="First":
        #print "+, First"
        #print seq, seqHalf
        if colorspaceFormat:
            seqHalf = hmmUtils.translateColorspace(seq[-1], seqHalf)
        #print seq[-1], seqHalf
        fullSeq = seq + seqHalf
        fullQual = qual + qualHalf
        fullGen = genomeDict[chr][chrStart:chrStart+len(fullSeq)]
        startPos = chrStart
        splitPos = chrStart+len(seq)-1
        seedHalf = "1"
    elif strand=="-" and readHalf=="First":
        #print "-, First"
        seq = hmmUtils.reverseComplement(seq)
        if colorspaceFormat:
            seqHalf = hmmUtils.translateColorspace(seq[-1], seqHalf)
        fullSeq = seq + seqHalf
        fullQual = qual[::-1] + qualHalf
        fullGen = hmmUtils.reverseComplement(genomeDict[chr][chrStart-len(seqHalf):chrStart+len(seq)])
        startPos = chrStart + len(fullSeq) - len(seqHalf)
        splitPos = chrStart+1
        seedHalf = "2"
    elif strand=="+" and readHalf=="Second":
        #print "+, Second"
        if colorspaceFormat:
            seqHalf = hmmUtils.translateColorspace(seq[0], seqHalf, False)
        fullSeq = hmmUtils.reverseComplement(seqHalf + seq)
        fullQual = qual[::-1] + qualHalf[::-1]
        fullGen = hmmUtils.reverseComplement(genomeDict[chr][chrStart-len(seqHalf):chrStart+len(seq)])
        startPos = chrStart + len(fullSeq) - len(seqHalf) 
        splitPos = chrStart+1
        seedHalf = "2"
    elif strand=="-" and readHalf=="Second":
        #print "-, Second"
        seq = hmmUtils.reverseComplement(seq)
        if colorspaceFormat:
            seqHalf = hmmUtils.translateColorspace(seq[0], seqHalf, False)
        fullSeq = hmmUtils.reverseComplement(seqHalf + seq)
        fullQual = qual + qualHalf[::-1]
        fullGen = genomeDict[chr][chrStart:chrStart+len(fullSeq)]
        startPos = chrStart
        splitPos = chrStart+len(seq)
        seedHalf = "1"
    else:
        raise hmmErrors.UnexpectedException("Didn't find combo for strand '%s' and readHalf '%s'" % (strand, readHalf))
    
    if colorspaceFormat:
        fullQual = "#" * len(fullSeq)

    vals = [title, strand, fullSeq, fullQual, fullGen, chr, startPos, splitPos, seedHalf, 0, 0]
    out.write("\t".join(str(x) for x in vals))
    out.write("\n")

        
    return 1



def _getObservationList(seedFile, numToRead=-1):
    """Converts the seed file data into a list of observations and qualities.
    Each observation is the list of matches/mismatches for one seed result
    and the qualities list is a parallel list with the quality score for each
    position.
    """
    
    obs = []
    quals = []
    
    count = 0
    for line in open(seedFile):
        count += 1
        if numToRead > 0 and count > numToRead:
            break
        
        try:
            [name, strand, readSeq, qual, genomeSeq, thisChr, start, stop, half, x, y] = line.split("\t")
        except ValueError:
            print "Illegal line:", line
            print line
            raise hmmErrors.InvalidInputException("Invalid line in seedFile")
        
        obsStr, qual, start =_convertOneToObs(readSeq, qual, genomeSeq)
        if obsStr == None:
            count -= 1
            continue
        obs.append(obsStr)
        quals.append(qual)
        
    return obs, quals

def _convertOneToObs(readSeq, qual, genomeSeq):
    """Takes in a read sequence, a genome sequence, and the associated quality string
    and returns an observation string and a quality list.
        
    For example:
    ACATGAGTACTA
    ||||  || |||
    ACATCCGTTCTA
    hhhhhhhhhhhh
    
    Would become "MMMMXXMMXMMM" and [40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40]
    """
    matchList = []
    qualList = hmmUtils.convertQualityStr(qual)
    qualList = [max(0, round(x, -1)/10) for x in qualList]
    isBegin = True
    startSite = len(readSeq)/3 #0 #len(readSeq)/2
    for i in range(startSite, len(readSeq)):
        try:
            if readSeq[i] == genomeSeq[i]:
                matchList.append("M")
            else:
                matchList.append("X")
        except IndexError:
            # this will happen if the read is at the end of the chromosome.  We can't use those reads
            return None, None, None

    return "".join(matchList), qualList[startSite:], startSite

def trainHMM(inputFile, initialHMM, trainedHMM, numObs=100000):
    """ Trains the input HMM with numObs from the input file.  Saves the trained HMM in the
    trainedHMM file.
    """
    
    obsList, quals = _getObservationList(inputFile, numObs)

    for i in range(len(obsList)):
        if len(obsList[i]) != len(quals[i]):
            raise hmmErrors.InvalidInputException("ACK!  the obs list and qual list are not the same for %s (%s vs %s): %s and %s" % (i, len(obsList[i]), len(quals[i]), obsList[i], quals[i]))
            
        
    hmm = hmmWithQuality.HMM(STATE_LIST, OBSERVE_LIST, 5)
    hmm.loadHMM(open(initialHMM))

    hmm.multiple_learn(obsList, quals, 100, True)

    print "Trained Values:"  
    hmm.dump() 
    
    hmm.saveHMM(open(trainedHMM, "w"))
    
def runHMM(inputFile, inputHMM, outputNoJunction, outputSmallJunction, outputJunction, anchorSize, chrName=None):
    """Runs the HMM on the input file using the inputHMM (which should be the trained HMM generated
    by trainHMM).  
    inputFile : the file of seeds from seedReads
    inputHMM : pickled trained HMM from trainHMM
    outputNoJunction : the output file for all the seeds which never transfer to State 2, i.e. are not junctions
    outputSmallJunction : the output file for the seeds which transfer to State 2 very late so the second part
        of the read is < configVals.ANCHOR_SIZE big.  These are used in the later rescue step.
    outputJunction : the junction reads -- the data we move forward with.
    chrName : optional.  If chrName is specified then only seeds on the given chromosome will be analyzed.
    
    This is also where we do our line conversions.  Read name becomes readName.split()[0]. Etc.
    """
    hmm = hmmWithQuality.HMM(STATE_LIST, OBSERVE_LIST, 5)
    inputHMM = inputHMM.replace("/", os.sep)
    hmm.loadHMM(open(inputHMM))

    outNo = open(outputNoJunction, "w")
    outSmall = open(outputSmallJunction, "w")
    outJunct = open(outputJunction, "w")
    
    countNo = 0
    countSmall = 0
    countJunct = 0
    
    originalPosition = 0
    for line in open(inputFile):
        originalPosition += 1
        try:
            [name, strand, readSeq, qual, genomeSeq, thisChr, start, stop, half, x, y] = line.split("\t")
        except ValueError:
            print "Illegal line:", line
            print line
            raise hmmErrors.InvalidFileFormatException("Invalid line in seed file.")
        
        if chrName:
            if thisChr != chrName:
                continue
            # save our line position.  Used later to re-sort the file into the original order
            name = "%s|%s" % (name.split()[0], originalPosition)
        else:
            name = name.split()[0]
            
        # re-create the line but with the fixed name
        line = "\t".join([name, strand, readSeq, qual, genomeSeq, thisChr, start, stop, half, x, y.strip()])

        
        obs, quals, start = _convertOneToObs(readSeq, qual, genomeSeq)

        if obs == None:
            continue
        
        try:
            traj = hmm.analyze(obs, quals)
        except:
            print line
            print obs
            print quals
            raise
            
        trajStr = "."*start + "".join(traj)
        switchPosition = trajStr.find("j")
        #print line.strip(), switchPosition, trajStr
        
        if switchPosition < 0:
            countNo += 1
            outNo.write(line)
        elif (len(trajStr) - switchPosition) < anchorSize:
            countSmall += 1
            outSmall.write("%s\t%s\t%s\n" % (line, switchPosition, trajStr))
        else:
            countJunct += 1
            outJunct.write("%s\t%s\t%s\n" % (line, switchPosition, trajStr))
            
    #print "Done! Found %s reads with no junctions, %s reads with junctions shorter than %s bases, and %s junction reads" % (countNo, countSmall, configVals.ANCHOR_SIZE, countJunct)
    
    return countNo, countSmall, countJunct
    
    
def matchSecondHalf(inputJunction, genomeDict, outputNotFound, outputMultiple, outputJunctionBed, outputShortIntronBed, minIntron,
                    maxIntron, maxWiggle, readLength, anchorSize):
    """After the HMM has determined where the splice junction begins in a read, this function is called to determine where the
    'second half' (the remainder of the read) is placed in the genome."
    inputJunction : the junction data from runHMM above
    genomeDict : a dictionary of chrName:chrSequence
    outputNotFound : the output file for reads where the second part is not able to be found
    outputMultiple : the output file for reads where the second part is found multiple times.  These are rescued later.
    outputJunctionBed : a bed file with the junctions found when the second part has a unique match
    outputShortIntronBed : a bed file with the junctions found when the second part has a unique match but it is less than 
            the minIntron size
    minIntron, maxIntron : the min and max intron sizes
    maxWiggle : the maximum number of base pairs that the intron will be slide left or right to try to find canonical splice junctions
    readLength : the length of the read.
    """
    
    outNo = open(outputNotFound, "w")
    outMult = open(outputMultiple, "w")
    outJunct = open(outputJunctionBed, "w")
    outShort = open(outputShortIntronBed, "w")
    
    if readLength % 2 == 0:
        theoryMax = readLength/2 * readLength/2 
    else:
        theoryMax = readLength/2 * (readLength/2 + 1)
    multiplier = 1200.0 / theoryMax
    
    juncTypes = {"no":0, "mult":0, "junc":0, "short":0}
    for line in open(inputJunction):
        matchType = _matchOneLine(genomeDict, line, outNo, outMult, outJunct, outShort, minIntron, maxIntron, maxWiggle, 
                                  multiplier, anchorSize)
        if not juncTypes.has_key(matchType):
            juncTypes[matchType] = 0
        juncTypes[matchType] += 1
        
    return juncTypes["no"], juncTypes["mult"], juncTypes["junc"], juncTypes["short"]
        
        
def _matchOneLine(genomeDict, line, outNo, outMult, outJunct, outShort, minIntron, maxIntron, maxWiggle, multiplier, anchorSize):
    """Called from matchSecondHalf, this method finds the second half for one read."""
    try:
        [name, strand, readSeq, qual, genomeSeq, chr, start, stop, half, score1, score2, splitPosition, trajStr] = line.split("\t")
    except ValueError:
        print line
        raise
    
    name, splitPosition, secondHalf, strand, isForward, start = _getParams(name, splitPosition, readSeq, strand, half, start)
    
    # get the next 10k of sequence
    genomeSeq = _getGenomeSeq(genomeDict, chr, start, maxIntron, isForward, strand)
    
    # get all the location where the seed matches and count the matches at that location
    allLocs = _getSeeds(splitPosition, secondHalf, genomeSeq, anchorSize)
    
    # if no matches, write to no output
    if len(allLocs) == 0:
        #name = "%s|noSeed" % (name)
        # Don't print out the genome seq anymore -- too big!  7/12/09
        outNo.write("\t".join([str(x) for x in [name, strand, readSeq, qual, "", 
                               chr, start, stop, half, score1, score2, splitPosition, trajStr]]))
        outNo.write("\n")
        return "no"
    
    # find all the indexes that have the maximum score
    bestKey = _findBestKey(allLocs)
    
    # if no best key, then there are multiple matches
    if bestKey == -1:
        outMult.write(line)
        return "mult"

    if maxWiggle > 0:
        splitPosition, bestKey = _doWiggle(readSeq, genomeSeq, qual, splitPosition, bestKey, maxWiggle, configVals.WIGGLE_SPLICE_SITES)
        secondHalf = readSeq[splitPosition:]
        
    start, stop, blockSizes, blockStarts, intronLength = _getBedVals(isForward, start, stop, len(readSeq), splitPosition, 
                                                                    strand, bestKey, len(secondHalf))
    
    score = _calculateScore(readSeq, genomeSeq, qual, splitPosition, bestKey, intronLength, minIntron, maxIntron,
                           multiplier)
    
    bedVals = [chr, start, stop, name, score, strand, start, stop, "0", "2", blockSizes, blockStarts]
    if intronLength < minIntron:
        outShort.write("\t".join(str(x) for x in bedVals))
        outShort.write("\n")
        return "short"
    else:
        #print "\t".join(str(x) for x in bedVals)
        outJunct.write("\t".join(str(x) for x in bedVals))
        outJunct.write("\n")
        return "junc"
        
def _getParams(name, splitPosition, readSeq, strand, half, start):
    """Helper function for _matchOneLine converts input parameters and returns usable values."""
    name = name.split()[0]
    
    splitPosition = int(splitPosition) 
    secondHalf = readSeq[splitPosition:]

    if strand == "F":
        strand = "+"
    elif strand == "R":
        strand = "-"

    # are we going 'forward' or 'reverse on the genome
    if half == "1":
        isForward = True
    else:
        isForward = False
        
    start = int(start)
    
    return name, splitPosition, secondHalf, strand, isForward, start

def _getGenomeSeq(genomeDict, chr, start, maxIntron, isForward, strand):
    """Returns the genome sequence from the beginning of the intron to the max intron size, in
    the correct direction."""
    
    if isForward:
        genomeSeq = genomeDict[chr][start:start+maxIntron]
    else:
        if strand == "+":
            genomeSeq = hmmUtils.reverseComplement(genomeDict[chr][start-maxIntron:start])
        else:
            genomeSeq = hmmUtils.reverseComplement(genomeDict[chr][start-maxIntron:start])
            
    return genomeSeq

def _getSeeds (splitPosition, secondHalf, genomeSeq, anchorSize):
    """Finds all the matches for the secondHalf within the genomeSeq."""
    allLocs = {}
    
    seed = secondHalf[:anchorSize]
    findStart = splitPosition
    while (True):
        nextLoc = genomeSeq.find(seed, findStart)
        if nextLoc < 0:
            break
        numMatch = 0
        numMM = 0

        try:
            for i in range(len(secondHalf)):
                if secondHalf[i] == genomeSeq[nextLoc+i]:
                    numMatch += 1
                else:
                    numMM += 1
            if numMM <= configVals.SECOND_HALF_MM_ALLOWED:
                allLocs[nextLoc] = numMatch
        except IndexError:
            # this happens when it continues past the end of the genome fragment,
            # we're done looking for matches for this guy!
            break
        findStart = nextLoc+1
        
        
    # TEST WITH A SECOND SEED
    if len(secondHalf) >= 2*anchorSize:
        seed2 = secondHalf[anchorSize:2*anchorSize]
        findStart = splitPosition + anchorSize
        while(True):
            nextLoc = genomeSeq.find(seed2, findStart)
            if nextLoc < 0:
                break
            # already found this with the other seed
            if allLocs.has_key(nextLoc - anchorSize):
                findStart = nextLoc + 1
                continue
            
            numMatch = 0
            numMM = 0
            try:
                for i in range(len(secondHalf)):
                    if secondHalf[i] == genomeSeq[nextLoc-anchorSize+i]:
                        numMatch += 1
                    else:
                        numMM += 1
                if numMM <= configVals.SECOND_HALF_MM_ALLOWED:
                    allLocs[nextLoc-anchorSize] = numMatch
            except IndexError:
                break
            findStart = nextLoc + 1
            
    return allLocs

def _findBestKey(allLocs):
    """Finds the location of the best match within the set of all locations."""
    maxMatches = max(allLocs.values())
    bestIndexList = []
    for k,v in allLocs.iteritems():
        if v == maxMatches:
            bestIndexList.append(k)

    # if there is one index with the best score, it's obvious.  Otherwise
    # if there are multiple with the best score but only one less than SECOND_HALF_ALT_INTRON, go with 
    # that one.  Otherwise, throw our hands up and say there are multiple best
    if len(bestIndexList) == 1:
        return(bestIndexList[0])
    else:
        # sorting can take a while -- some of these have thousands of matches.  Do it the simple way!
        # make sure the min val is less than SECOND_HALF_ALT_INTRON
        minVal = min(bestIndexList)
        if minVal > configVals.SECOND_HALF_ALT_INTRON:
            return -1
        # then make sure it's the ONLY one less than SECOND_HALF_ALT_INTRON
        for item in bestIndexList:
            if item > minVal and item < configVals.SECOND_HALF_ALT_INTRON:
                return -1
            
        # if it makes it past these checks, return the minVal
        return minVal

 
def _getBedVals(isForward, start, stop, readSeqLen, splitPosition, strand, bestKey, secondHalfLen):
    """Converts from our usable values into values ready to be written to a bed file."""
    if isForward:
        stop = bestKey + start + secondHalfLen
        blockSizes = "%s,%s," % (splitPosition, readSeqLen-splitPosition)
        blockStarts = "0,%s," % (bestKey)
        intronLength = stop - start - readSeqLen
    else:
        if strand == "+":
            #start = start + 1
            stop = start - bestKey - secondHalfLen
            blockSizes = "%s,%s," % (readSeqLen-splitPosition, splitPosition)
            blockStarts = "0,%s," % (start-splitPosition-stop)

            temp = stop
            stop = start
            start = temp
            intronLength = stop-start-readSeqLen
        if strand == "-":
            #start = start - 1    
            stop = start - bestKey - secondHalfLen
            blockSizes = "%s,%s," % (readSeqLen-splitPosition, splitPosition)
            blockStarts = "0,%s," % (start-splitPosition-stop)

            temp = stop
            stop = start
            start = temp
            intronLength = stop-start-readSeqLen
        
    return start, stop, blockSizes, blockStarts, intronLength
            
            
def _doWiggle(readSeq, genomeSeq, qual, readSplit, genomeSplit, maxWiggle, spliceList):
    """"wiggles" the edge of the intron (intron size remains the same but the start/stop is shifted
    either right or left up to maxWiggle bases).  The wiggle checks the sequence to make sure the
    wiggle is valid.  Wiggling is done to all 'canonical' splice sites found in the configValues.
    
    This version is the slow but prioritized version.  The wiggling is done in this order:
    1)  can we wiggle this junction to the first item in the spliceList using sequence as the guide?
    2)  can we wiggle this junction to the second item in the spliceList using sequence as the guide?
    4)  etc. though all the items in the spliceList...
    
    This function returns the read split, genome split, and splice site:
    read split:  where to divide the read, i.e. where the intron begins.
    genome splice: where to divide the genome, i.e. where the intron ends and the next exon begins.
    splice site: which splice site was matched.  For example "GT-AG"
    
    if none of the items in spliceList can be found, then None, None, site is returned.
    
    """

    leftSeq = genomeSeq[readSplit:readSplit+2]
    rightSeq = genomeSeq[genomeSplit-2:genomeSplit]
    seq = "%s-%s" % (leftSeq, rightSeq)
    
    for site in spliceList:
        siteRev = hmmUtils.reverseComplement(site) 
        if seq == site:
            return readSplit, genomeSplit
        elif seq == siteRev:
            return readSplit, genomeSplit
        
        x, y, s = _doWiggleWork(site, siteRev, readSeq, genomeSeq, readSplit, genomeSplit, maxWiggle, True)
        if x and y:
            return x, y
        # always check sequence when doing the wiggle
        #x, y, s = _doWiggleWork(site, siteRev, readSeq, genomeSeq, readSplit, genomeSplit, maxWiggle, False)
        #if x and y:
        #    return x,y
        
    return readSplit, genomeSplit

def _doWiggleWork(site1, site2, readSeq, genomeSeq, readSplit, genomeSplit, maxWiggle, useSeq=True):
    wiggle = 1
    foundSpliceSite = False
    foundLeft = False
    foundRight = False
    while (not foundSpliceSite and wiggle <= maxWiggle):
        # first wiggle left
        wigLeft = "%s-%s" % (genomeSeq[readSplit-wiggle:readSplit+2-wiggle], genomeSeq[genomeSplit-2-wiggle:genomeSplit-wiggle])
        #print "       wiggling seq1 %s left: %s" % (wiggle, wigLeft)
        if useSeq:
            foundSpliceSite = (wigLeft == site1 or wigLeft==site2) and (readSeq[readSplit-wiggle:readSplit] == genomeSeq[genomeSplit-wiggle:genomeSplit]) 
        else:
            foundSpliceSite = (wigLeft == site1 or wigLeft==site2)
        if foundSpliceSite:
            if wigLeft == site1:
                return readSplit-wiggle, genomeSplit-wiggle, site1
            else:
                return readSplit-wiggle, genomeSplit-wiggle, site2
        
        # then wiggle right
        wigRight = "%s-%s" % (genomeSeq[readSplit+wiggle:readSplit+2+wiggle], genomeSeq[genomeSplit-2+wiggle:genomeSplit+wiggle])
        #print "       wiggling seq1 %s right: %s" % (wiggle, wigRight)
        if useSeq:
            foundSpliceSite = (wigRight == site1 or wigRight==site2) and (readSeq[readSplit:readSplit+wiggle] == genomeSeq[genomeSplit:genomeSplit+wiggle])
        else:
            foundSpliceSite = (wigRight == site1 or wigRight==site2)
        if foundSpliceSite:
            if wigRight == site1:
                return readSplit+wiggle, genomeSplit+wiggle, site1
            else:
                return readSplit+wiggle, genomeSplit+wiggle, site2
        
        wiggle += 1
        
    return None, None, None

            
def _calculateScore(readSeq, genomeSeq, qual, splitPosition, genomeSecondPosition, intronLength, minIntron, maxIntron,
    multiplier):
    """Determines the score for a junction based on the length of each side along the junction and
    the quality scores/mismatches in the junction.  The score also takes into account the difference
    between the junction and a non-junction 'full-length-match'."""
    
    readFirst = readSeq[:splitPosition]
    readSecond = readSeq[splitPosition:]
    genomeFirst = genomeSeq[:splitPosition]
    genomeSecond = genomeSeq[genomeSecondPosition:genomeSecondPosition+len(readSeq)-splitPosition]
    genomeLeftMovedRight = genomeSeq[genomeSecondPosition-splitPosition:genomeSecondPosition]
    genomeRightMovedLeft = genomeSeq[splitPosition:len(readSeq)]
    qualFirst = qual[:splitPosition]
    qualSecond = qual[splitPosition:]
    
    #print "ReadFirst then genomeFirst then genomeLeftMovedRight"
    #print readFirst
    #print genomeFirst
    #print genomeLeftMovedRight
    #print "-"*60
    #print "ReadSecond then genomeSecond then genomeRightMovedLeft"
    #print readSecond
    #print genomeSecond
    #print genomeRightMovedLeft
    #print "=*"*60
    #print
    
    
    scoreFirst = _scoreHalf(readFirst, genomeFirst, qualFirst)
    scoreSecond = _scoreHalf(readSecond, genomeSecond, qualSecond)
    
    scoreLeftMovedRight = _scoreHalf(readFirst, genomeLeftMovedRight, qualFirst)
    scoreRightMovedLeft = _scoreHalf(readSecond, genomeRightMovedLeft, qualSecond)
    
    maxUngapped = max( scoreFirst*scoreRightMovedLeft, scoreSecond*scoreLeftMovedRight)

    if intronLength < minIntron:
        return 10
    elif intronLength > maxIntron:
        return 10
    
    #print"=" * 50
    #print readFirst
    #print genomeFirst
    #print qualFirst
    #print "-" * 30
    #print readSecond
    #print genomeSecond
    #print qualSecond
    #print "-" * 30
    #print scoreFirst * scoreSecond
    #print "=" * 50
    
    #return scoreFirst * scoreSecond
    #print intronPenalty, intronLength
    return ((scoreFirst * scoreSecond) - (0.5*maxUngapped) ) * multiplier
    
def _scoreHalf(readSeq, genomeSeq, qual, useQuality=True):
    """Scores one alignment.  At this point the strands have been assumed to already be taken care of."""
    
    #if genomeSeq != genomeSeq.upper():
    #    print "genomeSeq %s contains lower case letters!" % (genomeSeq)
    #    raise hmmErrors.UnexpectedException("Invalid input into _scoreHalf")
    
    # if things are masked out, they shouldn't count as matched (they don't count elsewhere)
    #genomeSeq = genomeSeq.upper()
    
    if len(readSeq) != len(genomeSeq):
        print "Invalid input to _scoreHalf"
        print readSeq
        print genomeSeq
        print qual
        raise hmmErrors.UnexpectedException("Invalid input into _scoreHalf")
    
    total = 0
    qualInts = hmmUtils.convertQualityStr(qual)
    for i in range(len(readSeq)):
        if readSeq[i] == genomeSeq[i]:
            # for now, no quality in the seeding
            if useQuality:
                total +=  (1-hmmUtils.qualityToProbability(qualInts[i], isSolexa=True))
            else:
                total += 1
        
    #print isFirst, isForward
    #print readSeq
    #print genomeSeq
    #print total
    #print "="*60
    
    return total

def removeDupReads(junctionBedIn, junctionBedOut, dupJuncBedOut=None):
    """The seeds can include multiple seeds for the same read.  These should all be right next to each other in the
    output file.  Here we reduce down to one junction for each read name:
    -- if there is a single junction to begin with, great
    -- if there are multiple junction and they are the same: great  (this will happen when both halves of the read match 
            in the same location)
    -- if there are multiple junctions and they are different: write to duplicate junction file, if any
    """
    
    out = open(junctionBedOut, "w")
    dupOut = None
    if dupJuncBedOut != None:
        dupOut = open(dupJuncBedOut, "w")
    currentRead = ""
    currentLines = []
    
    classifications = {}
    count = 0
    for line in open(junctionBedIn):
        count += 1
        try:
            readName = line.split("\t")[3]
        except:
            print "ERROR interpreting line %s" % count
            print line
            raise
        
        if readName == currentRead:
            currentLines.append(line)
        else:
            retType = _processOneReadSet(currentLines, out, dupOut)
            if not classifications.has_key(retType):
                classifications[retType] = 0
            classifications[retType] += 1
            currentRead = readName
            currentLines = [line]
            
    retType = _processOneReadSet(currentLines, out, dupOut)
    if not classifications.has_key(retType):
        classifications[retType] = 0
    classifications[retType] += 1
    
    #for k,v in classifications.iteritems():
    #    print k, v
    
            
def _processOneReadSet(currentLines, out, dupOut=None):
    """Handles all the junctions found for one read sequence, writing the read
    out to the single or duplicate file as described in removeDupReads."""
    if len(currentLines) == 0:
        return "0"
            
    if len(currentLines) == 1:
        #print "ONLY 1", currentLines[0]
        out.write(currentLines[0])
        return "1"
            
    (chr, start, stop, name, score1) = currentLines[0].split("\t")[:5]
    start = int(start)
    stop = int(stop)
    score1 = float(score1)
    wasDifferent = False
    maxScore = score1
    maxLine = currentLines[0]
    scores = [score1]
    for i in range(1, len(currentLines)):
        #print "Working on", currentLines[i]
        (x,y,z, name, score2) = currentLines[i].split("\t")[:5]
        score2 = float(score2)
        
        if currentLines[i].split("\t")[3] != currentLines[0].split("\t")[3]:
            raise hmmErrors.InvalidInput("Mismatched read names in _processOneReadSet: %s vs. %s" % (currentLines[i].split("\t")[3],currentLines[0].split("\t")[3]))
            
        if score2 > maxScore:
            maxScore = score2
            maxLine = currentLines[i]
            
        # need to make a decision -- if scores are equal, dump, but if one 'much better' than 
        # other than we should keep that one
        if (chr!=x or abs(start-int(y))>3 or abs(stop-int(z))>3):
            scores.append(score2)
     
    # now we have the max scoring line and all the scores for different lines
    # if there is only one score, all the lines were in the same place
    # if there are multiple scores, then if the top two scores are very similar, then we have different
    #       places, same scores, so we don't use it
    # but if the second score is much less than the first score, then just use the best scoring line
    if len(scores) == 1:
        # they are all in the same place
        out.write(maxLine)
        return "2-samePlace"
    else:   
        scores.sort()
        score1 = scores[-1]
        score2 =  scores[-2]   
        if (abs(score1 - score2) < configVals.MIN_SCORE_DIFF) and (chr!=x or start !=y or stop != z):
            for line in currentLines:
                score = float(line.split("\t")[4])
                if dupOut and abs(score-score1) < configVals.MIN_SCORE_DIFF:
                    dupOut.write(line)
            return "2-difPlace-eqScore"
        else:
            out.write(maxLine)
            return "2-difPlace-difScore"

        
def rescueMultipleSeeds(multipleFile, initialJunctions, outputFound, genomeDict, minIntron,
                    maxIntron, readLength, appendToFound):
    """Implements the rescue portion of HMMSplicer.  In this part of the algorithm, reads which 
    were unable to have a definite junction placed earlier are revisited.  The other junctions
    found in matchSecondHalf are used as a guide.  For any of the ambiguous reads, if one edge
    of the read matches to a known good read, then the second part of the read is checked to see
    if it supports the same junction.
    multipleFile : the set of reads which was set aside before.  These are both the outputSmallJunction
        from runHMM and the outputMultiple from matchSecondHalf
    initialJunctions : the set of junctions found uniquely in matchSecondHalf with the duplicates removed
        in removeDupReads.  These are taken as the 'known good' set of junctions
    outputFound : the set of junctions that we rescued
    genomeDict : the dictionary of chrName:chrSequence
    minIntron, maxIntron : the min and max intron values allowed
    readLength : the length of the reads, used to calculate the score
    appendToFound : whether the rescued reads are being written to a new file or appended to an existing file.
    """
    
    if appendToFound:
        # do append instead of "w" so we can easily combine small second side with multiple matches
        outFound = open(outputFound, "a")
    else:
        outFound = open(outputFound, "w")
        
    if readLength % 2 == 0:
        theoryMax = readLength/2 * readLength/2 
    else:
        theoryMax = readLength/2 * (readLength/2 + 1)
    multiplier = 1200.0 / theoryMax
    
    # read known junctions into a dictionary of chr|one edge|other edge = # (do collapsing at the same time)
    knownJunctionsLeft, knownJunctionsRight = _createInitialJunctionDicts(initialJunctions)
        
    count = 0
    for line in open(multipleFile):
        if line.startswith("track"):
            continue

        if len(line) < 2:
            continue
    
        try:
            [name, strand, readSeq, qual, genomeSeq, chr, start, stop, half, score1, score2, splitPosition, trajStr] = line.split("\t")
        except ValueError:
            print line
            raise
    
        if name.find("|") >= 0:
            name, type = name.split("|")
            if type == "noSeed":
                continue

        # get info required to write output, including first edge
        name, splitPosition, secondHalf, strand, isForward, start = _getParams(name, splitPosition, readSeq, strand, half, start)
        
        secondLength = len(secondHalf)
        allowedMM = secondLength / 6
        # is this the same as another edge?  (remember to look for wiggles, too)
        secondPosition = -1
        if isForward:
            firstEdge = start + splitPosition
            if knownJunctionsLeft.has_key(chr) and knownJunctionsLeft[chr].has_key(firstEdge):
                for possibleRight in knownJunctionsLeft[chr][firstEdge]:
                    genomeSeq = genomeDict[chr][start:start+maxIntron+readLength]
                    if len(genomeSeq) == 0:
                        raise hmmErrors.InvalidInputException("possibleRight in knownJunctionsLeft. %s, %s, %s, %s, %s, %s" % (chr,
                                                                len(genomeDict[chr]), start, maxIntron, readLength, splitPosition))
                    #print chr, start, start+maxIntron+readLength, len(genomeSeq)
                    # do stuff: look at all actual seeds for match to second half
                    #print start, firstEdge, knownJunctionsLeft[chr][firstEdge][0]
                    possibleSecondPosition = knownJunctionsLeft[chr][firstEdge][0] - start
                    if possibleSecondPosition < splitPosition:
                        continue
                    if possibleSecondPosition > len(genomeSeq):
                        continue
                    testGenome = genomeSeq[possibleSecondPosition:(possibleSecondPosition+secondLength)] 
                    numMM = 0
                    for i in range(len(testGenome)):
                        if testGenome[i] != secondHalf[i]:
                            numMM += 1
                    if numMM <= allowedMM:
                        secondPosition = possibleSecondPosition
                        break
        
        else:
            if strand == "+":
                start = start + 1
                firstEdge = start - splitPosition
            elif strand == "-":
                start = start - 1
                firstEdge = start - splitPosition
                
            if knownJunctionsRight.has_key(chr) and knownJunctionsRight[chr].has_key(firstEdge):
                for possibleLeft in knownJunctionsRight[chr][firstEdge]:
                    if start > (maxIntron+readLength):
                        genomeSeq = hmmUtils.reverseComplement(genomeDict[chr][start-maxIntron-readLength:start])
                    else:
                        genomeSeq = hmmUtils.reverseComplement(genomeDict[chr][:start])
                    #print origGenomeSeq
                    #print genomeSeq[:50]
                    # do stuff: look at all actual seeds for match to second half
                    possibleSecondPosition = start - knownJunctionsRight[chr][firstEdge][0]
                    if possibleSecondPosition < splitPosition:
                        continue
                    if possibleSecondPosition > len(genomeSeq):
                        continue
                    testGenome = genomeSeq[possibleSecondPosition:(possibleSecondPosition+secondLength)] 
                    numMM = 0
                    for i in range(len(testGenome)):
                        if testGenome[i] != secondHalf[i]:
                            numMM += 1
                    if numMM <=allowedMM:
                        secondPosition = possibleSecondPosition
                        break


        # if not found, write to not found and continue
        if secondPosition < 0:
            continue

        #print isForward
        #print line
        #print secondPosition
        #print origGenomeSeq
        #print genomeSeq[:50]

        # back to the broken start for now
        if not isForward:
            if strand == "+":
                start = start - 1
            else:
                start = start + 1

        start, stop, blockSizes, blockStarts, intronLength = _getBedVals(isForward, start, stop, len(readSeq), splitPosition, 
                                                                    strand, secondPosition, len(secondHalf))
        score = _calculateScore(readSeq, genomeSeq, qual, splitPosition, secondPosition, intronLength, minIntron, maxIntron,
                           multiplier)
    
        bedVals = [chr, start, stop, name, score, strand, start, stop, "0", "2", blockSizes, blockStarts]
        if intronLength >= minIntron:
            outFound.write("\t".join(str(x) for x in bedVals))
            outFound.write("\n")
        

def _createInitialJunctionDicts(initialJunctions):
    """Creates the set of 'known good' junctions from the bed file of initial junctions."""
    knownJunctionsLeft = {}
    knownJunctionsRight = {}
    
    count = 0
    for line in open(initialJunctions):
        if line.startswith("track"):
            continue
        #count +=1 
        #if count > 10:
        #    break
        pieces = line.split()
        chr = pieces[0]
        leftEdge, rightEdge = hmmUtils.getEdges(int(pieces[1]), pieces[10], pieces[11])
        
        if not knownJunctionsLeft.has_key(chr):
            knownJunctionsLeft[chr] = {}
            knownJunctionsRight[chr] = {}
            
        if leftEdge < rightEdge:
            if not knownJunctionsLeft[chr].has_key(leftEdge):
                knownJunctionsLeft[chr][leftEdge] = []
            knownJunctionsLeft[chr][leftEdge].append(rightEdge)
            if not knownJunctionsRight[chr].has_key(rightEdge):
                knownJunctionsRight[chr][rightEdge] = []
            knownJunctionsRight[chr][rightEdge].append(leftEdge)
        else:
             raise hmmErrors.InvalidInputException("leftEdge > rightEdge.  %s > %s in %s" % (leftEdge, rightEdge, line.strip()))

    return knownJunctionsLeft, knownJunctionsRight
        



             
