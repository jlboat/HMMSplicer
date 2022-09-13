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
# Created: April, 2009
#
# Functions related to finding and outputting junctions
#   (This does not include the HMM code, which is in junctionHMM)

import hmmErrors
import configVals
import tempfile


            
def divideByGTAG(inputBed, outputGood, outputBad, genome):
    """Divides the input bed into junctions that match configVals.SPLICE_SITES
    and junctions that don't.  The genome dictionary is required
    to check the intron edges for each junction."""
    
    outGood = open(outputGood, "w")
    outBad = open(outputBad, "w")
    
    # dictionary is much faster than a list for matching
    spliceDict = {}
    for item in configVals.SPLICE_SITES:
        spliceDict[item] = 0
    
    for line in open(inputBed):
        if line.startswith("track") or len(line) < 2:
            continue
        
        try:
            [chr, start, stop, name, score, strand, thStart, thStop, rgb, blockCount, blockSizes, blockStarts] = line.split("\t")
        except:
            print(line)
            raise
        
        if int(blockCount) != 2:
            raise hmmErrors.InvalidFileFormatException("Illegal line does not have 2 blocks.\n%s" % line)

        start = int(start)
        stop = int(stop)
        [size1, size2] = [int(x) for x in blockSizes.split(",")[:2]]
        [start1, start2] = [int(x) for x in blockStarts.split(",")[:2]]
        leftEdge = start + size1
        rightEdge = start + start2  # start2 is relative to chr start
        
        leftSeq = genome[chr][leftEdge:leftEdge+2]
        rightSeq = genome[chr][rightEdge-2:rightEdge]
        seq = "%s-%s" % (leftSeq, rightSeq)
        
        if seq in spliceDict:
            outGood.write(line)
        else:
            outBad.write(line)


def collapseCloseJunctionsMemory(inputBed, outputBed, withinBp):
    """collapseCloseJunctions must read the whole set of junctions into memory to determine
    if any are duplicated.  This can cause memory problems if the junction set is too large.
    So, this function breaks the input bed into a set of temporary bed files, each at most
    1,000,000 lines long.  These are collapsed separately then combined together
    and collapsed again."""

    inputBedList = []

    tempSubset = tempfile.NamedTemporaryFile(delete=False)
    inputBedList.append(tempSubset.name)
    count = 0
    for line in open(inputBed):
        count += 1
        if count % 50000 == 0:
            #print "gotcha", count
            tempSubset.close()
            tempSubset = tempfile.NamedTemporaryFile(delete=False)
            inputBedList.append(tempSubset.name)
        tempSubset.write(line)

    tempSubset.close() 
    #print inputBedList
    if len(inputBedList) == 1:
        collapseCloseJunction(inputBedList[0], outputBed, withinBp)
    else:
        collapsedList = []
        for f in inputBedList:
            collapsedName = f+".collapsed"
            #print collapsedName
            collapseCloseJunctions(f, collapsedName, withinBp)
            collapsedList.append(collapsedName)

        collapseCloseJunctionList(collapsedList, outputBed, withinBp)


def collapseCloseJunctionList(inputList, outputBed, withinBp):
    """Collapses together many bed files, saving the overall results to outputBed."""

    # concatenate all the files
    tempConcat = tempfile.NamedTemporaryFile(delete=False)
    for infile in inputList:
        for line in open(infile):
            if infile != inputList[0] and line.startswith("track"):
                continue
        
            tempConcat.write(line)
        
    tempConcatName = tempConcat.name
    tempConcat.close()
    
    # then collapse
    collapseCloseJunctions(tempConcatName, outputBed, withinBp)
        


def collapseCloseJunctions(inputBed, outputBed, withinBp):
    """Collapses junctions that are within 'withinBp' bp of each other.  The most probable junction
    of the close junctions are used.  Probabilities are converted appropriately (see collapseJunctions
    for description).
    Also collapses identical junctions in the same step.
    Use withinBp=0 to only collapse identical junctions
    """
    
    junct = _readAndCombine(inputBed, withinBp)
        
    # then go through each dictionary item
    #   -- if the list is length 1, just write the line
    #   -- if the list is longer than 1, combine probabilities, adjust the name to include the number, and write 1
    out = open(outputBed, "w")
    out.write("track name=collapsedJunctions description='Junctions' useScore=1\n")
    for chr, junctions in junct.items():
        for (leftEdge, x,rightEdge, y, intronLength), junctionList in junctions.items():
            if len(junctionList) == 1:
                out.write(junctionList[0][1].strip())
                #pieces = junctionList[0][0].split()
                #pieces.pop(4)
                #pieces.insert(4, "100")
                #out.write("\t".join(pieces))
                out.write("\n")
            else:
                out.write(_combineLines(junctionList, leftEdge, rightEdge))
                out.write("\n")
                
def _readAndCombine(inputBed, withinBp):
    """Helper for collapseCloseJunctions reads in the inputBed into a dictionary."""
    junct = {}
    
    # collapse a
    count = 0
    for line in open(inputBed):
        count += 1
        #if count % 100000==0:
        #    print count
        if line.startswith("track"):
            #out.write(line.strip())
            #out.write(" useScore=1\n")
            continue
        
        [chr, start, stop, name, score, strand, thStart, thStop, rgb, blockCount, blockSizes, blockStarts] = line.split("\t")
        score = float(score)
        if chr not in junct:
            junct[chr] = {}
        
        if int(blockCount) != 2:
            #print "Illegal line does not have 2 blocks"
            #print line
            continue
        
        start = int(start)
        stop = int(stop)
        [size1, size2] = [int(x) for x in blockSizes.split(",")[:2]]
        [start1, start2] = [int(x) for x in blockStarts.split(",")[:2]]
        leftEdge = start + size1
        rightEdge = start + start2  # start2 is relative to chr start
        intronLength = rightEdge - leftEdge
        
        toCombine = []
        for (other) in list(junct[chr].keys()):
            (otherMinLeft, otherMaxLeft, otherMinRight, otherMaxRight, otherLength) = other
            if otherLength != intronLength:
                continue
            
            if otherMaxLeft < (leftEdge-withinBp) or otherMinLeft > (leftEdge+withinBp):
                continue
            
            if otherMaxRight < (rightEdge-withinBp) or otherMinRight > (rightEdge+withinBp):
                continue
            
            toCombine.append(other)
            
        allLines = [ (score, line, leftEdge, rightEdge) ]
        minLeft = maxLeft = leftEdge
        minRight = maxRight = rightEdge
        for (other) in toCombine:
            (otherMinLeft, otherMaxLeft, otherMinRight, otherMaxRight, intronLength) = other
            minLeft = min(minLeft, otherMinLeft)
            maxLeft = max(maxLeft, otherMaxLeft)
            minRight = min(minRight, otherMinRight)
            maxRight = max(maxRight, otherMaxRight)
            
            allLines.extend(junct[chr][other])
            del junct[chr][other]
            
        junct[chr][ (minLeft, maxLeft, minRight, maxRight, intronLength) ] = allLines
        
    return junct

def _combineLines(junctionList, leftEdge, rightEdge):
    """Helper method for collapseCloseJunctions.
    Takes several junction lines and combines them into a single line, then returns that line.
    The score for the combined junction is increased accordingly.
    """
    
    # the only things we have to adjust are the score and the name
    scoreSoFar = 0
    previousNumBasesCovered = 0
    minStart = -1
    maxStop = -1
    # keep a count of all the left/right edges so we can pick the most common
    leftEdgeCount = {}
    rightEdgeCount = {}
    
    # first we need to sort the junctionList by score, in reverse (highest score first)
    junctionList.sort(reverse=True)
    name = junctionList[0][1].split("\t")[3]
    countPlus = 0
    countMinus = 0
    numJunctions = 0
    for (score, line, leftEdge, rightEdge) in junctionList:
        pieces = line.split("\t")
        #print "previous start: %s, this start: %s" % (minStart, pieces[1])
        start = int(pieces[1])
        stop = int(pieces[2])
        strand = pieces[5]
        if strand == "+":
            countPlus += 1
        elif strand == "-":
            countMinus += 1
        else:
            print(line)
            raise hmmErrors.InvalidInputException("ERROR!  strand value %s not valid. " % strand)
        
        name = pieces[3]
        if pieces[3].find("|junc=") > 0:
            numCollapsed = int(pieces[3].split("junc=")[-1])
            numJunctions += numCollapsed
        else:
            numJunctions += 1
        
        if previousNumBasesCovered == 0:
            previousNumBasesCovered = (leftEdge - start) + (stop - rightEdge)
            scoreSoFar = float(pieces[4])
        else:
            # we only want to count bases grown on the outer sides (start and stop) and ignore bases on the inner edges
            newBases = max(0, minStart-start) + max(0, stop-maxStop)
            scoreSoFar = scoreSoFar + (newBases / float(newBases+previousNumBasesCovered) ) * float(pieces[4])
            previousNumBasesCovered = previousNumBasesCovered + newBases
            

        if minStart > 0:
            minStart = min(minStart, start)
        else:
            minStart = start
            
        maxStop = max(maxStop, stop)
        
        if leftEdge not in leftEdgeCount:
            leftEdgeCount[leftEdge] = 0
        leftEdgeCount[leftEdge] += 1
        if rightEdge not in rightEdgeCount:
            rightEdgeCount[rightEdge] = 0
        rightEdgeCount[rightEdge] += 1
        
    maxLeft = max(leftEdgeCount.values())
    for k, v in leftEdgeCount.items():
        if v == maxLeft:
            useLeft = k
            break
    maxRight = max(rightEdgeCount.values())
    for k, v in rightEdgeCount.items():
        if v == maxRight:
            useRight = k
            break
    if countPlus >= countMinus:
        strand = "+"
    else:
        strand = "-"
        
    pieces = junctionList[0][1].split("\t")
    namePieces = pieces[3].split("|")
    rootName = ""
    for piece in namePieces:
        if not piece.startswith("junc="):
            rootName += piece + "|"
    finalName = rootName + ("junc=%s" % numJunctions)
    blockStarts = "0,%s," % (useRight - minStart)
    blockSizes = "%s,%s," % ( (useLeft-minStart), (maxStop-useRight) )

    return "\t".join(str(x) for x in [pieces[0], minStart, maxStop, finalName, scoreSoFar, 
                      strand, minStart, maxStop, "0", "2", blockSizes, blockStarts
                      ])
    

def convertScoreToProb(score):
    """Converts a score in a bed file (between 0 and 1000) into a probability (between 0 and 1)."""
    return float(score) / 1000.0

def convertProbToScore(prob):
    """Converts a probability (between 0 and 1) into a score for a bed file (between 0 and 1000)."""
    return str(prob * 1000)

    
def readJunctionsFromBed(inputBed, saveWholeLine=False, wiggle=3):
    """Reads the junctions from a bed file into a dictionary.  
    For example, used by several methods in scoreJunctions to read in EST/gene annotation files."""
    junct = {}
    for line in open(inputBed):
        if line.startswith("track"):
            continue
        
        if len(line) < 2:
            continue
        
        [chr, start, stop, name, score, strand, thStart, thStop, rgb, blockCount, blockSizes, blockStarts] = line.split()
        if chr not in junct:
            junct[chr] = {}
            
        start = int(start)
        stop = int(stop)
        blockCount = int(blockCount)
        if blockSizes.endswith(","):
            sizes = [int(x) for x in blockSizes.split(",")[:-1]]
            starts = [int(x) for x in blockStarts.split(",")[:-1]]
        else:
            sizes = [int(x) for x in blockSizes.split(",")]
            starts = [int(x) for x in blockStarts.split(",")]
        for i in range(1, blockCount):
            leftEdge = start + starts[i-1] + sizes[i-1]
            rightEdge = start + starts[i]  
            if saveWholeLine:
                if (leftEdge, rightEdge) not in junct[chr]:
                    junct[chr][(leftEdge, rightEdge)] = []
                junct[chr][(leftEdge, rightEdge)].append(line.strip())
            else:
                if rightEdge - leftEdge > 15:
                    found = False
                    for wL in range(-wiggle,wiggle):
                        for wR in range(-wiggle,wiggle):
                            if (leftEdge+wL, rightEdge+wR) in junct[chr]:
                                found = True
                                break
                    if not found:
                        junct[chr][ (leftEdge, rightEdge) ] = (name, score, strand)
    
    return junct

def hasJunction(junc, chr, leftEdge, rightEdge, wiggle):
    """Determines whether the given junction dictionary has a junction.
    wiggle is the amount plus/minus to count as the 'same place'."""
    
    for i in range(leftEdge-wiggle, leftEdge+wiggle+1):
        for j in range(rightEdge-wiggle, rightEdge+wiggle+1):
            try:
                if (i, j) in junc[chr]:
                    return True
            except KeyError:
                return False
            
    return False

def readGeneModelsBed(inputBed):
    """Reads in a EST or PlasmoDB bed file and returns a junction dictionary in the format junction[chr][(start,stop)] = (name, score)
    where name is the plasmodb/est name + "_" + intron number and score is always 1000.  
    
    This function can handle cases where the 'read' (bed line) has more than one junction, and the case where multiple ESTs cover the
    same junction (in this case the later ones are ignored for now)."""
    
    d = {}
    
    for line in open(inputBed):
        [chr, start, stop, name, score, strand, x, y, z, blockCount, blockSizes, blockStarts] = line.split()
        
        if chr not in d:
            d[chr] = []
            
        blockCount = int(blockCount)
        blockSizeList = [int(x) for x in blockSizes.split(",")[:-1]]
        blockStartList = [int(x) for x in blockStarts.split(",")[:-1]]
        
        # one block = no introns = we aren't interested
        if blockCount == 1:
            continue
            
        # if this exact junction has been found before (will heppen for ESTs) then just skip for now

        for i in range(1, blockCount):
            rightEdge = blockStartList[i] + start
            leftEdge = blockStartList[i-1] + blockSizeList[i-1]
            #print "        ", rightEdge, leftEdge
            
            if (leftEdge, rightEdge) in d[chr]:
                continue
            
            d[chr][(leftEdge, rightEdge)] = (name+"_"+str(i-1), 1000)

             
    return d


