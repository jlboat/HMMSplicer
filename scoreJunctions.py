#!/usr/bin/env python
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
# processJunctions.py
#
# Author: Michelle Dimon
# Created: September, 2009
#
# Functions not used by the HMMSplicer program that might be useful for interpreting the results.

import processJunctions
import hmmUtils
import commands
import random


##################################################################################################
#  Matching Results to Known Junctions (Annotations / ESTs)
##################################################################################################

def measureSpecificity(foundJunctionsBed, estJunctionsBed, wiggle=0, goodFile=None, badFile=None):
    """Measures the number and percent of junctions in the foundJunctionsBed file that overlap
    the known junctions in the estJunctionsBed.  goodFile and badFile are optional, but allow
    the user to divide the foundJunctionsBed into those that do and do not match known."""
    # read in all the est junctions
    ests = processJunctions.readJunctionsFromBed(estJunctionsBed)
    
    return measureSpecificityUsingDict(foundJunctionsBed, ests, wiggle, goodFile, badFile)
    
def measureSpecificityUsingDict(foundJunctionsBed, ests, wiggle=0, goodFile=None, badFile=None):
    """Prints the percent of found junctions that overlap with ests.
    This version takes in a dictionary of ests created by processJunctions.readJunctionsFromBed.
    This dictionary can take a while to create, especially if you have a lot of EST data to load,
    so if you are calling this function over and over (for example, a series of hmmUtils.filterBedFile
    then scoreJunctions.measureSpecificityUsingDict is useful for seeing how score thresholds affect
    junction accuracy), then this version of is recommended over measureSpecificity.
    """
    
    overlaps = 0
    noOverlaps = 0
    savedLines = None
    if goodFile or badFile:
        savedLines = {}
    
    # do it slow for now
    # read in all the junctions we found
    found = {}  # found[chr][(left, right)] = count
    for line in open(foundJunctionsBed):
        if line.startswith("track"):
            continue

        if len(line) < 3:
            continue
    
        pieces = line.split("\t")
        #print pieces
        if len(pieces) == 0:
            continue
        
        if pieces[0].startswith("Plasmodium_falciparum"):
            pieces[0] = pieces[0].split("|")[1].replace("MAL", "chr")

        if pieces[0].startswith("psu|Pf"):
            pieces[0] = "chr" + str(int(pieces[0].split()[0].split("|")[1].split("_")[1]))

        
        if not found.has_key(pieces[0]):
            found[pieces[0]] = {}
        
        leftEdge, rightEdge = hmmUtils.getEdges(int(pieces[1]), pieces[10], pieces[11])
        if not found[pieces[0]].has_key( (leftEdge,rightEdge) ):
            found[pieces[0]][(leftEdge,rightEdge)] = 0
        found[pieces[0]][(leftEdge,rightEdge)] += 1
        
        if goodFile or badFile:
            if not savedLines.has_key(pieces[0]):
                savedLines[ pieces[0]] = {}
            savedLines[pieces[0]][(leftEdge,rightEdge)] = line
    
    # for every one of our junction, do they overlap with an est?
    if goodFile != None and badFile != None:
        goodOut = open(goodFile, "w")
        badOut = open(badFile, "w")
        
    for chr, edgeDict in found.iteritems():
        for (leftEdge, rightEdge) in edgeDict.keys():
            foundOne = False
            for x in range(leftEdge-wiggle, leftEdge+wiggle+1):
                for y in range(rightEdge-wiggle, rightEdge+wiggle+1):
                    if ests.has_key(chr):
                        if ests[chr].has_key( (x, y) ):
                            foundOne = True
            if foundOne:
                overlaps += found[chr][(leftEdge,rightEdge)]
                if goodFile != None:
                    #goodOut.write("\t".join([chr, str(leftEdge), str(rightEdge), str(overlaps)]))
                    #goodOut.write("\n")
                    goodOut.write(savedLines[chr][(leftEdge, rightEdge)])
            else:
                noOverlaps += found[chr][(leftEdge,rightEdge)]
                if badFile != None:
                    #badOut.write("\t".join([chr, str(leftEdge), str(rightEdge), str(noOverlaps)]))
                    #badOut.write("\n")
                    badOut.write(savedLines[chr][(leftEdge,rightEdge)])
            
    if (noOverlaps + overlaps) > 0:
        print "%s overlapped but %s did not.  %.2d%% overlapped" % (overlaps, noOverlaps, (overlaps*100.0)/(noOverlaps+overlaps))

        return (overlaps*100.0) / (noOverlaps + overlaps)
    else:
        print "No junctions found!"
        return 0

def measureSensitivity(foundJunctionsBed, estJunctionsBed, wiggle=3):
    """Counts the number of junctions in the estJunctionsBed that are covered by
    junctions in the foundJunctionsBed.  Designed for use with the uncollapsed results
    file, creates a histogram of the number of junctions and how many reads cover them
    (i.e. there are 230 junctions with 0 predicted junctions covering them, 4530 junctions
    with 1 predicted junction covering them, etc.)
    """
    # read in all the est junctions
    ests = processJunctions.readJunctionsFromBed(estJunctionsBed)
    
    return measureSensitivityUsingDict(foundJunctionsBed, ests, wiggle)
    
    
def measureSensitivityUsingDict(foundJunctionsBed, ests, wiggle=3):
    """Prints a histogram of coverage level per junction in the EST.  For example, the 
    ESTs with no coverage will be listed as zero, but the ESTs where 10 junctions found the
    EST junction will be listed as 10.
    This version takes in a dictionary of ests created by processJunctions.readJunctionsFromBed.
    This dictionary can take a while to create, especially if you have a lot of EST data to load,
    so if you are calling this function over and over (for example, a series of hmmUtils.filterBedFile
    then scoreJunctions.measureSensitivityUsingDict is useful for seeing how score thresholds affect
    junction coverage), then this version of is recommended over measureSensitivity.
    """
    
    counts = {}
    
    # do it slow for now
    # read in all the junctions we found
    found = {}  # found[chr][(left, right)] = count
    for line in open(foundJunctionsBed):
        if line.startswith("track"):
            continue
        if len(line) < 3:
            continue
        pieces = line.split("\t")
        #print pieces
        if len(pieces) == 0:
            continue
        
        if pieces[0].startswith("Plasmodium_falciparum"):
            pieces[0] = pieces[0].split("|")[1].replace("MAL", "chr")

        if pieces[0].startswith("psu|Pf"):
            pieces[0] = "chr" + str(int(pieces[0].split()[0].split("|")[1].split("_")[1]))
        
        if not found.has_key(pieces[0]):
            found[pieces[0]] = {}

        #print pieces[1], pieces[10], pieces[11]
        leftEdge, rightEdge = hmmUtils.getEdges(int(pieces[1]), pieces[10], pieces[11])
        if not found[pieces[0]].has_key( (leftEdge,rightEdge) ):
            found[pieces[0]][(leftEdge,rightEdge)] = 0
        found[pieces[0]][(leftEdge,rightEdge)] += 1
    
    # for every est junction, how many of ours cover it
    numJunct = 0
    for chr, edgeDict in ests.iteritems():
        numJunct += len(edgeDict)
        for (leftEdge, rightEdge) in edgeDict.keys():
            numFound = 0
            for x in range(leftEdge-wiggle, leftEdge+wiggle+1):
                for y in range(rightEdge-wiggle, rightEdge+wiggle+1):
                    if found.has_key(chr):
                        if  found[chr].has_key( (x, y) ):
                            numFound += found[chr][(x, y)]
            if not counts.has_key(numFound):
                counts[numFound] = 0
            counts[numFound] += 1
            
    allkeys = counts.keys()
    allkeys.sort()
    
    print "Total number of junctions in the dict is: %s" % (numJunct)
    moreThanZero = 0
    for k in allkeys:
        if k > 0:
            moreThanZero += counts[k]
        print "%s\t%s" % (k,counts[k])
    print "total > zero is %s" % moreThanZero

    return moreThanZero


def generateROCcurve(junctionBed, ests, stepSize=100, wiggle=0, ungappedPairDbf=None):
    """Generates an ROC curve to measure good/bad junctions in the junction bed file based on the location
    of the paired end in the ungapped Dbf (DeRisi Bed file).
    
    The ROC curve is calculated from the 'good' and 'bad' matches, with the criteria being the score cut-off.
    
    The 'good' or 'bad' can take paired end data into account (reads are only judged if their pair mapped and
    'good' reads must map facing their pair on the same chromosome and bad must map on a different chromosome.
    """
    
    if ungappedPairDbf:
        goodNames, badNames = findGoodBadReadsWithPairedEnd(junctionBed, ungappedPairDbf, estBed)
    else:
        goodNames, badNames = findGoodBadReads(junctionBed, ests, wiggle)
 
    scores = []
    scoredDict = {}
    countGood = 0
    countBad = 0
    countStacked = 0
    countUnclear = 0

    for line in open(junctionBed):
        if line.startswith("track"):
            continue
        try:
            [chr, start, stop, name, score] = line.split()[:5]
        except ValueError:
            print "2"
            continue
        score = float(score)
        #score = random.randint(0, 10)
        if goodNames.has_key(name):
            isGood = True
        elif badNames.has_key(name):
            isGood = False
        else:
            continue
    
        scores.append(score)
        
        if isGood:
            countGood += 1
            
        else:
            countBad += 1
            
        scoredDict[name] = (score, isGood)
        
    print "Counted %s good and %s bad" % (countGood, countBad)
    
    print scores[:3]
    scores.sort(reverse=True)
    print scores[:5]
    # use every 100th score
    for x in scores[::stepSize]:
        numGoodAbove = 0
        numBadAbove = 0
        for k, (score, real) in scoredDict.iteritems():
            if score >= x:
                if real:
                    numGoodAbove += 1
                else:
                    numBadAbove += 1
                    
        print x, numGoodAbove, numBadAbove, numGoodAbove / float(countGood) , numBadAbove / float(countBad)
        
def findGoodBadReads(junctionBed, ests, wiggle):
    """Finds 'good' and 'bad' reads, only using the estBed as a guide."""
    
    #ests = processJunctions.readJunctionsFromBed(estBed)
    
    goodNames = {}
    badNames = {}
    for line in open(junctionBed):
        if line.startswith("track"):
            continue
        
        try:
            [chr, start, stop, name, score, strand, thStart, thStop, rgb, blockCount, blockSizes, blockStarts] = line.split()
        except ValueError:
            print "1"
            break
        
        hasEst = testIfEst(ests, chr, int(start), int(blockCount), blockSizes, blockStarts, wiggle)
        
        if hasEst:
            goodNames[name] = 1
        else:
            badNames[name] = 1
        
    return goodNames, badNames
    
def findGoodBadReadsWithPairedEnd(junctionBed, ungappedPairDbf, estBed):
    """Goes through a set of junctions and determines which are 'good' and which are 'bad' based
    on the paired end info and ESTs.  Most reads are thrown out as unsure, but the clearly good and
    clearly bad ones are written out and can be used in further tests.

    For each junction, there are several possibilities:
    -- the other end of the paired read was not mapped ungapped
    -- the other end of the paired read matches to a 'good' spot.  Good is defined as same chr, correct direction,
          with 0-300 bp between the read ends (=100-400 insert size)
    -- other end of the paired read matches to a 'bad' spot.  Bad is defined as a different chr, wrong direction, or
          more than 10kb downstream from the junction read
    -- other end of the paired read matches to an 'unclear' spot.  This means same chr, correct direction, but >300 bp
          between the read ends (essentially these are probably 'bad' but we're being generous)

    """
    
    ests = processJunctions.readJunctionsFromBed(estBed)
    
    # the junctions will be much smaller than all ungapped, so read in all junctions, root name : (chr, end, strand, score)
    # (Note: end is the stop for strand F but the start for strand R.  It's essentially where the arrow would be drawn in a picture)
    junctionDict = {}
    for line in open(junctionBed):
        if line.startswith("track"):
            continue

        [chr, start, stop, name, score, strand, thStart, thStop, rgb, blockCount, blockSizes, blockStarts] = line.split()
        root = name.split("#")[0]
        score = float(score)
        
        hasEst = testIfEst(ests, chr, int(start), int(blockCount), blockSizes, blockStarts)
        
        if strand == "+":
            junctionDict[root] = (chr, int(stop), strand, score, hasEst)
        else:
            junctionDict[root] = (chr, int(start), strand, score, hasEst)
            
    print "Read in %s junctions" % (len(junctionDict))
    #print junctionDict.keys()[:5]
            
    # go through the ungapped file and, if the root name is found in the junction dict, add to new dict (pairs) of root name: (chrJ, 
    # endJ, strandJ, scoreJ, chrP, endP, strandP) where J is junction read and P is paired end.  Remove junction from junction dict.
    pairsDict = {}
    for line in open(ungappedPairDbf):
        if line.startswith("track"):
            continue
        
        [name, seq, qual, status, qStart, qEnd, chr, chrStart, chrEnd, strand, 
                    blockCount, blockSizes, blockStarts] = line.split()
        root = name.split("#")[0]
        if junctionDict.has_key(root):
            if strand == "F":
                pairsDict[root] = (junctionDict[root][0], junctionDict[root][1], junctionDict[root][2], junctionDict[root][3],
                              chr, int(chrEnd), "+", junctionDict[root][4])
            else:
                pairsDict[root] = (junctionDict[root][0], junctionDict[root][1], junctionDict[root][2], junctionDict[root][3],
                              chr, int(chrStart), "-", junctionDict[root][4])
            del junctionDict[root]
            
    # print out the number of reads still in the junction dict.
    print "After filling the pairs dict, there were %s junctions still remaining in the junction dict (unmapped paired end)" % (len(junctionDict))
    
    # for each read in the pairs dict, judge if good or not.  build new dictionary (scored) of name:(score, isGood) 
    # count unclear as you go, but toss out, don't include is dict.
    # keep track of max and min scores
    # count number of actual good and actual bad
    # Note: at this point the dict is only 'good' and 'bad'.  all 'unmapped' and 'unclear' are tossed.
    scores = []
    scoredDict = {}
    countEstGood = 0
    countEstBad = 0
    countGood = 0
    countBad = 0
    countStacked = 0
    countUnclear = 0
    countExcludedForEst = 0 
    for k, (chrJ, arrowJ, strandJ, scoreJ, chrP, arrowP, strandP, estResult) in pairsDict.iteritems():
        test = testIfGood(chrJ, arrowJ, strandJ, scoreJ, chrP, arrowP, strandP)
        
        if estResult:
            countEstGood += 1
        else:
            countEstBad += 1
        
        if test == "G":
            if estResult:
                countGood += 1
                scoredDict[k] = (scoreJ, True)
                #maxScore = max(maxScore, scoreJ)
                #minScore = min(minScore, scoreJ)
                scores.append(scoreJ)
            else:
                countExcludedForEst += 1
        elif test == "U":
            countUnclear += 1
            continue
        elif test == "S":
            countStacked += 1
            continue
        else:
            countBad += 1
            scoredDict[k] = (scoreJ, False)
            
        
        
    print "After scoring, there were %s good, %s unclear, %s stacked, %s excluded for no matching EST and %s bad junctions" % (countGood, countUnclear, countStacked, countExcludedForEst, countBad)
    print "And there were %s good EST junctions and %s bad ESTjunctions" % (countEstGood, countEstBad)
    
    goodNames = {}
    badNames = {}
    countFoundInScored = 0
    for line in open(junctionBed):
        if line.startswith("track"):
            continue
        
        [chr, start, stop, name, score, strand, thStart, thStop, rgb, blockCount, blockSizes, blockStarts] = line.split()
        if not chr.startswith("chr"):
            continue
        
        root = name.split("#")[0]
        if scoredDict.has_key(root):
            countFoundInScored += 1
            if scoredDict[root][1]:
                #name = "%s|GOOD" % (name)
                goodNames[name] = 1
                #goodOut.write("\t".join([chr, start, stop, name, score, strand, thStart, thStop, rgb, blockCount, blockSizes, blockStarts]))
                #goodOut.write("\n")
            else:
                #name = "%s|bad" % (name)
                badNames[name] = 1
                #badOut.write("\t".join([chr, start, stop, name, score, strand, thStart, thStop, rgb, blockCount, blockSizes, blockStarts]))
                #badOut.write("\n")
        
    print "Found %s lines in the scored dict" % countFoundInScored
    
    return goodNames, badNames

def testIfEst(ests, chr, start, blockCount, blockSizes, blockStarts, wiggle=0):
    """Determines whether the ests dictionary contains the junction described by the
    chr, start, blockCount, blockSizes."""
    if int(blockCount) != 2:
        print "ERROR!  the block count isn't 2!  %s, %s, %s, %s, %s" % (chr, start, blockCount, blockSizes, blockStarts)
        return False
    
    (leftEdge, rightEdge) = hmmUtils.getEdges(start, blockSizes, blockStarts)
    
    return processJunctions.hasJunction(ests, chr, leftEdge, rightEdge, wiggle)


def testIfGood(chrJ, arrowJ, strandJ, scoreJ, chrP, arrowP, strandP):
    """Tests if the junction is 'good' -- i.e. same chr, facing each other, 0-300 bp between arrows.
    Returns "G", "B", or "U" for good, bad or unclear."""
    
    # same chromosome
    if chrJ != chrP:
        #print "diff chr", chrJ, arrowJ, strandJ, scoreJ, chrP, arrowP, strandP
        return "B"
    
    # facing each other (this and next set)
    if strandJ == strandP:
        #print "same strand", chrJ, arrowJ, strandJ, scoreJ, chrP, arrowP, strandP
        return "S"
    
    if arrowJ <= arrowP:
        if strandJ == "-":
            if abs(arrowJ - arrowP) < 100:
                return "S"
            else:
                #print "Not facing 1", chrJ, arrowJ, strandJ, scoreJ, chrP, arrowP, strandP
                return "S"
    else:
        if strandJ == "+":
            if abs(arrowJ - arrowP) < 100:
                return "S"
            else:
                #print "Not facing 2", chrJ, arrowJ, strandJ, scoreJ, chrP, arrowP, strandP
                return "S"
        
    # check distance between arrows
    arrowDistance = abs(arrowJ - arrowP)
    if arrowDistance < 300:
        return "G"
    elif arrowDistance < 10000:
        return "U"
    else:
        #print "Distance > 10000", chrJ, arrowJ, strandJ, scoreJ, chrP, arrowP, strandP
        return "B"
    
def getAllEdges(start, blockSizes, blockStarts):
    """Returns a list of all the edges in the line.  Safe if there are any
    number of junctions in the line."""
    
    if blockSizes.endswith(","):
        sizes = [int(x) for x in blockSizes.split(",")[:-1]]
        starts = [int(x) for x in blockStarts.split(",")[:-1]]
    else:
        sizes = [int(x) for x in blockSizes.split(",")]
        starts = [int(x) for x in blockStarts.split(",")]
        
    retval = []
    for i in range(len(starts)-1):
        leftEdge = start + starts[i] + sizes[i]
        rightEdge = start + starts[i+1]
        retval.append( (leftEdge, rightEdge) )
        
    return retval


##################################################################################################
#  Comparing Result files
##################################################################################################

def vennDiagram(bed1File, bed2File, only1Output=None, only2Output=None, bothOutput=None):
    """Counts the number of splice junctions unique to bed1, unique to bed2, and shared in both.
    The optional output files allow the splice junctions to be divided so they can be further processed.
    For example, measureSpecificity can be run on each output file to determine if the accuracy is 
    different for junctions found in both vs. those found in a single result."""
    
    bed1 = processJunctions.readJunctionsFromBed(bed1File, True)
    bed2 = processJunctions.readJunctionsFromBed(bed2File, True)
    
    count1 = 0
    count2 = 0
    countBoth = 0
    
    out1 = None
    if only1Output:
        out1 = open(only1Output, "w")
    out2 = None
    if only2Output:
        out2 = open(only2Output, "w")
    both = None
    if bothOutput:
        both = open(bothOutput, "w")

    for chr, chrJunct in bed1.iteritems():
        for (start,stop) in chrJunct:
            if bed2.has_key(chr):
                if bed2[chr].has_key( (start, stop) ):
                    if both:
                        for line in bed1[chr][(start,stop)]:
                            both.write(line)
                            both.write("\n")
                    del bed2[chr][(start,stop)]
                    countBoth += 1
                else:
                    count1 += 1
                    if out1:
                        line = bed1[chr][(start,stop)][0]
                        pieces = line.split()
                        bedVals = [chr, start-10, stop+10, pieces[3], pieces[4], pieces[5], start-10, stop+10, pieces[8], pieces[9],
                                   "10,10", "0,%s"%(stop-start+10)]
                        out1.write("\t".join(str(x) for x in bedVals))
                        out1.write("\n")
            else:
                count1 += 1
                if out1:
                    line = bed1[chr][(start,stop)][0]
                    pieces = line.split()
                    bedVals = [chr, start-10, stop+10, pieces[3], pieces[4], pieces[5], start-10, stop+10, pieces[8], "2",
                               "10,10", "0,%s"%(stop-start+10)]
                    out1.write("\t".join(str(x) for x in bedVals))
                    out1.write("\n")
            
                
    count2 = sum( len(chrJunct) for chrJunct in bed2.values())
    if out2:
        for chr, chrJunct in bed2.iteritems():
            for (start,stop) in chrJunct:
                line = bed2[chr][(start,stop)][0]
                pieces = line.split()
                bedVals = [chr, start-10, stop+10, pieces[3], pieces[4], pieces[5], start-10, stop+10, pieces[8], "2",
                           "10,10", "0,%s"%(stop-start+10)]
                out2.write("\t".join(str(x) for x in bedVals))
                out2.write("\n")
    
    print "There were %s in both, %s in the first one and %s in the second one" % (countBoth, count1, count2)


##################################################################################################
#  Handling GFF and other file formats
##################################################################################################

def makeGeneDict(gffInput):
    """Reads a gffInput file into a dictionary of d[id] = [chromosome, strand, start, stop, [exons]]
    """
    d = {}

    for line in open(gffInput):
        if line.startswith("##FASTA"):
            break

        # skip the line if it starts with "##" (it's a comment)
        if line.startswith("##") or len(line) < 1:
            continue

        # figure out if it's an exon
        mainPieces = line.split("\t")

        if len(mainPieces) != 9:
            print "Illegal line found!  Ignoring this line:"
            print line
            continue
        
        if not mainPieces[2] == "exon" and not mainPieces[2] == "gene":
            continue

        id, chromosome, strand, start, stop = parseGFFLine(mainPieces, isParent=mainPieces[2]=="gene")

        if not d.has_key(id):
            # this is a new id in the dictionary
            d[id] = [chromosome, strand, 0, 0, []]
            
        if mainPieces[2] == "gene":
            d[id][2] = start
            d[id][3] = stop
        else:
            # if the dictionary already has this id, just add this exon.  If the strand is + then this exon
            # comes AFTER the previous line(s) read in.  If the strand is - then this exon is before the 
            # line(s) we've already read it.  
            if strand == "+":
                d[id][4].append((start, stop))  
            else:
                d[id][4].insert(0, (start, stop))


    return d

def findGene(geneDict, chr, start, stop, blockSizes, blockStarts):
    """Given a gene dictionary (such as one created with makeGeneDict above, and a 
    junction defined by the chr, start, stop, blockSizes, and blockStarts, returns
    the gene record for the gene with that junction.  Returns None if none can be found.
    Will return the first one found, so the geneDict should not have overlapping genes."""
    for id,(geneChr, strand, geneStart, geneStop, exons) in geneDict.iteritems():
        if geneChr != chr:
            continue
        
        if start > geneStart and stop < geneStop:
            # here, check the left and right egdes match
            left, right = hmmUtils.getEdges(start, blockSizes, blockStarts)
            right += 1
            #print "Looking for %s-%s from (%s      %s      %s)" % (left, right, start, blockSizes, blockStarts)
            #print "    in %s: %s" % (id, exons)
            foundLeft = False
            foundRight = False
            for (eStart, eStop) in exons:
                if left == eStop:
                    foundLeft = True
                if right == eStart:
                    foundRight = True
                    
            if foundLeft and foundRight:
                #print "GOT IT!!!!!!!!!!!!!!!!!!!!!"
                return id
            else:
                #if we found the gene but not the exon edges, then stop processing
                break

    return None


def parseGFFLine(pieces, isParent):
    # the chromosome is the part after the "|" symbol from the first piece
    chromosome = pieces[0]
    start = int(pieces[3])
    stop = int(pieces[4])
    strand = pieces[6]
    
    # now the hardest part.  The ID is parsed out of the text piece at the end.
    if isParent:
        # the first item is the ID
        id = pieces[8].split(";")[0].split("=")[1].strip()
    else:
        # take the "PARENT" which is second
        id = pieces[8].split(";")[2].split("=")[1].strip()

    return id, chromosome, strand, start, stop


##################################################################################################
#  RPKM Related
##################################################################################################

def measureRPKMperJunction(junctionBed, tophatRPKMresults, gffFile):
    """Given the TopHat RPKM results and a junction bed and a gffFile with annotations, returns
    the RPKM value for each gene in the gff file."""
    coveragePerJunction = readTophatRPKM(tophatRPKMresults)
    
    histAll = {}
    minRPKM = 0
    maxRPKM = 0
    print "reading in coverage per junction"
    for gene, rpkm in coveragePerJunction.iteritems():
        if not histAll.has_key(rpkm):
            histAll[rpkm] = 0
        histAll[rpkm] += 1
        minRPKM = min(minRPKM, rpkm)
        maxRPKM = max(maxRPKM, rpkm)
        
    print "making gene dict"
    geneDict = makeGeneDict(gffFile)
    
    print "getting coverage per gene"
    histFound = {}
    for line in open(junctionBed):
        if line.startswith("track"):
            continue
        
        pieces = line.split()
        id = findGene(geneDict, pieces[0], int(pieces[1]), int(pieces[2]), pieces[10], pieces[11])
        if id:
            if coveragePerJunction.has_key(id):
                rpkm = coveragePerJunction[id]
                if not histFound.has_key(rpkm):
                    histFound[rpkm] = 0
                histFound[rpkm] += 1
                #print rpkm

    
    print "printing results"
    #for x in xrange(minRPKM, maxRPKM):
    for x in xrange(0, 20):
        allVal = 0
        if histAll.has_key(x):
            allVal = histAll[x]
        foundVal = 0
        if histFound.has_key(x):
            foundVal = histFound[x]
        if allVal > 0:
            percent = foundVal * 100.0 / allVal
        else:
            percent = "N/A"
        print x, allVal, foundVal, percent
        
    sumAll = 0
    sumFound = 0
    for x in xrange(21, maxRPKM):
        if histAll.has_key(x):
            sumAll += histAll[x]
        if histFound.has_key(x):
            sumFound += histFound[x]
    print ">200", sumAll, sumFound, sumFound * 100.0 / sumAll
    
    
def readTophatRPKM(tophatRPKMresults):
    
    results = {}
    minval = 100
    maxval = 0

    for line in open(tophatRPKMresults):
        [gene, rpkm] = line.split()
        rpkm = float(rpkm)
        minval = min(minval, rpkm)
        maxval = max(maxval, rpkm)
        rpkm = int( round(rpkm, 0)) 
        bin = rpkm / 10
        results[gene] = bin


    print "The min rpkm was %s and the max was %s" % (minval, maxval)
    return results 


    



        
    
            
            
    
    
    
    

    
        

    
