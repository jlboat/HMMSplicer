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
##! /usr/bin/env python
#
# hmmUtils.py
#
# This script contains base-level utilities used in HMMSplicer.  For example, this script is
# where the utilities to parse a fasta file into a dictionary are kept. 
#
# Author: Michelle Dimon, May 2009
#

from types import StringTypes
import random
import string
import hmmErrors
import hmmWithQuality
from numpy import array
from math import log

##################################################################################################
#  Utilities to work with input files (fasta, fastq, etc)
##################################################################################################






def convertQualityStr(strVersion, offset=33):
    """Converts a quality string to a list of integers.
    an offset of 33 is Phred style and an offet of 64 is Solexa style"""
    return map( lambda x: ord(x)-offset, strVersion )

def convertToQualityStr(intList, offset=33):
    """Converts a list of integers to a quality string."""
    return "".join( map(lambda x: chr(x+offset), intList) )

def convertSolexaToPhred(solexaQual):

    qualInts = convertQualityStr(solexaQual, 64)
    adjusted = [ int(round(solexaQualToPhredQual(x), 0)) for x in qualInts]
    return convertToQualityStr(adjusted, 33)
    
def solexaQualToPhredQual(solexaQual):
    return 10 * log(1+10 ** (solexaQual / 10.0)) / log(10)


def qualityToProbability(qual, isSolexa=True):
    """Given a score, returns the probability associated with that score.
    This is the probability that the base call is an error, so a very low number
    (like 10 ^ -5) is good quality.
    Calculations are done differently for solexa and sanger, default is solexa"""
    
    if isSolexa:
        return 1 / (1 + 10 ** (qual/10.0))
    else:
        return 10 ** (-qual / 10.0)
    
def convertFastqOldSolexaToPhred(inputSolexaFastq, outputPhredFastq):
    """For Solexa quality values from Solexa Pipeline before 1.3.
    The integer values for Solexa go from -5 to 40 and are 
    'asymptotically similar' to Phred quality values but derived
    differently.  Plus the offset is different.
    Converts an entire file with old solexa quality score into
    a phred quality score file."""
    
    out = open(outputPhredFastq, "w")
    for (title, seq, qual) in FastqIterator(inputSolexaFastq):
        qual = convertSolexaToPhred(qual)
        out.write("@%s\n%s\n+%s\n%s\n" % (title, seq, title, qual))

def convertFastqNewSolexaToPhred(inputSolexaFastq, outputPhredFastq):
    """For Solexa quality values from Solexa Pipeline 1.3 and above.
    The integer values mean the same thing but the offset is different.
    Converts an entire file with new solexa quality scores into a 
    phred qualityscore file.
    """
    
    out = open(outputPhredFastq, "w")
    for (title, seq, qual) in FastqIterator(inputSolexaFastq):
        intVals = convertQualityStr(qual, 64)
        qual = convertToQualityStr(intVals, 33)
        out.write("@%s\n%s\n+%s\n%s\n" % (title, seq, title, qual))


def colorspaceToPhred(inputColorspaceFastq, outputPhredFastq):
    """Translates a colorspace fastq file into a phred-based 
    regular sequence file."""

    out = open(outputPhredFastq, "w")
    for (title, csSeq, quality) in FastqIterator(inputColorspaceFastq):
        seq = translateColorspace(csSeq[0], csSeq[1:])
        out.write("@%s\n%s\n+%s\n%s\n" % (title, seq, title, quality[1:]))

def colorspaceFastaToPhred(inputColorspaceFasta, outputPhredFastq):

    out = open(outputPhredFastq, "w")
    for (title, csSeq) in FastaIterator(inputColorspaceFasta):
        seq = translateColorspace(csSeq[0], csSeq[1:])
        out.write("@%s\n%s\n+%s\n%s\n" % (title, seq, title, "#"*len(seq)))


def fastaToFastq(inputFasta, outputFastq):

    out= open(outputFastq, "w")
    for(title, seq) in FastaIterator(inputFasta):
        out.write("@%s\n%s\n+%s\n%s\n" % (title, seq, title, "#"*len(seq)))

def translateColorspace(initialLetter, csSeq, goForward=True):
    """Given the initial nucleotide and a colorspace string, returns the                            
     corresponding nucleotide string, NOT INCLUDING the initial letter.                                           
    If goForward is False then the "initialLetter" is actual taken to be                                             
    the final letter.  To handle these cases, the reverse complement of                                              
    the initial letter is appended to the REVERSED csSeq.  The result                                      
    is translated and then the reverseComplement is returned.                                                
    """

    if not goForward:
        #print "%s + %s =" % (csSeq, initialLetter),                                                          
        initialLetter = reverseComplement(initialLetter)
        csSeq = csSeq[::-1]

    seq = ""
    previousLetter = initialLetter
    for s in csSeq:
        thisLetter = translateOneFromColorspace(previousLetter, s)
        seq += thisLetter
        previousLetter = thisLetter

    if not goForward:
        seq = reverseComplement(seq)

    return seq

def translateOneFromColorspace(previous, color):
    """Given the previous nucleotide and a color, returns this nucleotide."""

    if previous == "A":
        if color == "0" or color == "A":
            return "A"
        elif color == "1" or color == "C":
            return "C"
        elif color == "2" or color == "G":
            return "G"
        elif color == "3" or color == "T":
            return "T"
        else:
            raise hmmErrors.InvalidInputException("Given previousNt of '%s' and thisColor of '%s', unable to assign nucleotide" % (previous, color))

    elif previous== "C":
        if color == "0"or color == "A":
            return "C"
        elif color == "1" or color == "C":
            return "A"
        elif color == "2" or color == "G":
            return "T"
        elif color == "3" or color == "T":
            return "G"
        else:
            raise hmmErrors.InvalidInputException("Given previousNt of '%s' and thisColor of '%s', unable to assign nucleotide"% (previous, color))


    elif previous== "G":
        if color == "0"or color == "A":
            return "G"
        elif color == "1" or color == "C":
            return "T"
        elif color == "2" or color == "G":
            return "A"
        elif color == "3" or color == "T":
            return "C"
        else:
            raise hmmErrors.InvalidInputException("Given previousNt of '%s' and thisColor of '%s', unable to assign nucleotide"% (previous, color))


    elif previous== "T":
        if color == "0"or color == "A":
            return "T"
        elif color == "1" or color == "C":
            return "G"
        elif color == "2" or color == "G":
            return "C"
        elif color == "3" or color == "T":
            return "A"
        else:
            raise hmmErrors.InvalidInputException("Given previousNt of '%s' and thisColor of '%s', unable to assign nucleotide"% (previous, color))

    else:
        raise hmmErrors.InvalidInputException("Given previousNt of '%s' and thisColor of '%s', unable to assign nucleotide"% (previous, color))





def reverseComplement(seq):
    """return the reverse complement of a sequence.
    Original code from Kael and Dale."""
    seq=seq.upper()
    # complement
    compl = complement(seq)
    # reverse
    return compl[::-1]

def complement(seq,transl=None):
    """Return complement of seq.
    Original code from Kael and Dale.
    """
    transl = string.maketrans('aAcCgGtTnNxX-\t\n ','tTgGcCaAnNxX-\t\n ')
    compl = seq.translate(transl)
    return compl

def randomlySelectFromFile(inputFile, outputFile, numToSelect, numLines):
    """Randomly selects 'numToSelect' lines from inputFile and puts them
    in outputFile.
    
    numLines is the total number of lines in the file.
    """
    
    # pick 'numToSelect' numbers from range(0, numLines)
    linesToUse = random.sample(xrange(numLines), numToSelect)
    
    linesDict = {}
    for lineNum in linesToUse:
        linesDict[lineNum] = 0
    
    # read through inputFile and if line number is in random set, write to output file
    count = 0
    out = open(outputFile, "w")
    for line in open(inputFile):     
        if linesDict.has_key(count):
            out.write(line)
        count += 1
        
    out.close()
    
def randomlySelectFromFastqFile(inputFastq, outputFile, numToSelect, numRecords, startPos=None, stopPos=None, useFasta=False):
    """Randomly selects 'numToSelect' fastq records from inputFile and puts them
    in outputFile.
    
    numRecords is the total number of records in the fastq file.
    """
    
    # pick 'numToSelect' numbers from range(0, numLines)
    linesToUse = random.sample(xrange(numRecords), numToSelect)
    
    linesDict = {}
    for lineNum in linesToUse:
        linesDict[lineNum] = 0
    
    # read through inputFile and if line number is in random set, write to output file
    count = 0
    out = open(outputFile, "w")
    if useFasta:
        for (title, seq) in FastaIterator(inputFastq):
            if linesDict.has_key(count):
                if startPos != None:
                    seq = seq[startPos:]
                if stopPos != None:
                    seq = seq[:stopPos]
                out.write(">%s\n%s\n" % (title, seq))
            count += 1
    else:
        for (title, seq, qual) in FastqIterator(inputFastq):
            if linesDict.has_key(count):
                if startPos!= None:
                        seq = seq[startPos:]
                if stopPos != None:
                    seq = seq[:stopPos]
                out.write("@%s\n%s\n+%s\n%s\n" % (title, seq, title, qual))
            count += 1

        
    out.close()
    
def convertFastqToFasta(inputFastq, outputFasta):
    """Converts a fastq file into a fasta file."""
    out = open(outputFasta, "w")
    for (titleStr, seqStr, qualityStr) in FastqIterator(inputFastq):
        out.write(">%s\n%s\n" % (titleStr, seqStr))
        
        
def getEdges(start, blockSizes, blockStarts):
    """Returns the left and right edge of a junction.  There must be exactly
    one junction in the line or an exception is raised.
    
    The start, blockSizes, and blockStarts are from BED file format."""
    if blockSizes.endswith(","):
        sizes = [int(x) for x in blockSizes.split(",")[:-1]]
        starts = [int(x) for x in blockStarts.split(",")[:-1]]
    else:
        sizes = [int(x) for x in blockSizes.split(",")]
        starts = [int(x) for x in blockStarts.split(",")]
        
    if len(starts) > 2:
        raise hmmErrors.InvalidInputException("getEdges: ERROR! one junction per line.  Input was: %s, %s, %s" % (start, blockSizes, blockStarts))
    
    leftEdge = start + starts[0] + sizes[0]
    rightEdge = start + starts[1]
    
    return leftEdge, rightEdge  
        
def divideReads(inputFastq, outputFastq, solexaFormat=False):
    """Divides the reads in half and writes both halves out to the outputFastq file.
    This is the funciton used to create read-halves for the seeding step.
    The information about the other half of the read (seq/quality) is saved in the
    title so it doesn't have to be looked up later."""
    
    out = open(outputFastq, "w")
    
    # put the other half of the sequence/quality in the fasta title so we don't have to look it up later
    # (with big fastq files, we may not even be able to put it all in memory to look it up!)
    for (titleStr, seqStr, qualityStr) in FastqIterator(inputFastq):
        # we are going to skip reads with "N" later, may as well not process them through bowtie, etc.
        if seqStr.find("N") >= 0:
            continue
        if solexaFormat:
            qualityStr = convertToQualityStr( convertQualityStr(qualityStr, 64) )
        titleStr = titleStr.replace("|", "_")
        if len(seqStr) % 2 == 0:
            le = len(seqStr)/2
            out.write("@%s|First|%s|%s\n%s\n+\n%s\n" % (titleStr, seqStr[le:], qualityStr[le:], seqStr[:le], qualityStr[:le])) 
            out.write("@%s|Second|%s|%s\n%s\n+\n%s\n" % (titleStr, seqStr[:le], qualityStr[:le], seqStr[le:], qualityStr[le:]))
        else:
            # for reads of uneven length, use a smaller half but include the extra base in the fastq title.  i.e. if
            # the read is 45 bp long, each half is only 22 bases, but have the other 23 bases in the title 
            le = len(seqStr)/2
            out.write("@%s|First|%s|%s\n%s\n+\n%s\n" % (titleStr, seqStr[le:], qualityStr[le:], seqStr[:le], qualityStr[:le])) 
            out.write("@%s|Second|%s|%s\n%s\n+\n%s\n" % (titleStr, seqStr[:le+1], qualityStr[:le+1], seqStr[le+1:], qualityStr[le+1:]))
            

def divideReadsColorspace(inputFastq, outputFastq):
    """Divides the reads in half and writes both halves out to the outputFastq file.
    This is the funciton used to create read-halves for the seeding step.
    The information about the other half of the read (seq/quality) is saved in the
    title so it doesn't have to be looked up later."""
    
    out = open(outputFastq, "w")
    
    # put the other half of the sequence/quality in the fasta title so we don't have to look it up later
    # (with big fastq files, we may not even be able to put it all in memory to look it up!)
    for (titleStr, seqStr, qualityStr) in FastqIterator(inputFastq):
        # we are going to skip reads with "N" later, may as well not process them through bowtie, etc.
        if seqStr.find("N") >= 0:
            continue
        titleStr = titleStr.replace("|", "_")
        le = (len(seqStr)-1)/2
        out.write("@%s|First|%s|%s\n%s\n+\n%s\n" % (titleStr, seqStr[le+1:], qualityStr[le+1:], seqStr[:le+1], qualityStr[:le+1])) 
        out.write("@%s|Second|%s|%s\n%s\n+\n%s\n" % (titleStr, seqStr[:le+1], qualityStr[:le+1], seqStr[le+1:], qualityStr[le+1:]))
        
        
        
def fastaDictionary(inFile, chrName=None):
    """return a dictionary of fasta title to sequence strings.
    This loads the entire fasta file into memory.  This can be quite large.  For the human genome,
    this requires almost 3 GB of memory.
    """

    d = {}
    for (title, seq) in FastaIterator(inFile):
        title = title.split()[0]
        if not chrName:
            d[title] = seq
        elif chrName == title:
            d[title] = seq
            return d

    if chrName:
        raise hmmErrors.InvalidInputException("fastaDictionary: Not able to find input chr %s in the fasta file." % chrName)
    
    return d

def trimFastq(infile, outfile, trimFront, trimRear):
    """Trims trimLen number of bases from the end of the fastq seq and quality."""
    out = open(outfile, "w")
    
    for (title, seq, qual) in FastqIterator(infile):
        #print seq
        if trimFront > 0:
            seq = seq[trimFront:]
            qual = qual[trimFront:]
        #print seq
        if trimRear > 0:
            seq = seq[:-trimRear]
            qual = qual[:-trimRear]
        #print seq
        #break
        out.write("@%s\n%s\n+%s\n%s\n" % (title, seq, title, qual))
        
    out.close()
    
    
def filterFasta(inputFasta, outputFasta, titleStr):
    """Writes lines with the titleStr string in them into the output file.
    Does a simple titleline.find(titleStr) > 0 to determine if titleStr is 
    in the title."""
    out = open(outputFasta, "w")
    
    writing = False
    for line in open(inputFasta):
        if line.startswith(">"):
            if line.find(titleStr) > 0:
                writing = True
                out.write(line)
            else:
                writing = False
        else:
            if writing:
                out.write(line)

#
# Iterator
#
def FastaIterator(fh):
    """return an iterator of Records found in file handle, fh.
        Original code from Kael and Dale.
    """
    def readTotitle(fh):
        """returns a tuple ([lines before the next title line], next tile line)
        """
        preLines = []
        while True:
            l = fh.readline().strip()
            if l.startswith('>'):
                return (preLines,l)
            elif l == '':
                return preLines,None
            else:
                preLines.append(l)


    if type(fh) in StringTypes:
        fh = file(fh)
    
    preLines,nextTitleLine =readTotitle(fh)

    while nextTitleLine != None:
        title = nextTitleLine[1:].rstrip()
        preLines,nextTitleLine=readTotitle(fh)
        yield (title,''.join(map(lambda x: x.rstrip(),preLines)))

        
def FastqIterator(fh):
    """return an iterator of Records found in file handle, fh.
    Original code from Kael and Dale.
    """
    def readTotitle(fh, titleChar):
        """returns a tuple ([lines before the next title line], next tile line)
        """
        preLines = []
        while True:
            l = fh.readline().strip()
            if l.startswith(titleChar):
                return (preLines,l)
            elif l == '':
                return preLines,None
            else:
                preLines.append(l)

    if type(fh) in StringTypes:
        fh = file(fh)
    
    preLines,nextTitleLine =readTotitle(fh,'@')

    while nextTitleLine != None:
        seqTitle = nextTitleLine[1:].rstrip()
        preLines,nextTitleLine=readTotitle(fh,'+')
        qualTitle = nextTitleLine[1:].rstrip()
        if len(qualTitle.strip()) > 0 and seqTitle != qualTitle:
            print seqTitle
            print preLines
            print qualTitle
            raise hmmErrors.InvalidFastq("Error in parsing: @title sequence entry must be immediately followed by corresponding +title quality entry.")
        seqLines = preLines
        qualLines = []
        for i in range(len(seqLines)): # Quality characters should be the same length as the sequence
            qualLines.append( fh.readline().strip() )
        
        preLines,nextTitleLine=readTotitle(fh,'@')

        yield (seqTitle, ''.join(seqLines), ''.join(qualLines))
        
        
##################################################################################################
#  Utilities to work with output files
##################################################################################################


def filterBedFile(inputBed, outputBed, scoreFilterSingle, scoreFilterMultiple, newName=""):
    """Refilters a set of junction results based on the score.
    inputBed: Should be the full (non-score-filtered) junction results
    outputBed: The filter results
    multipleScore : the score for junctions with support from two or more reads
    singletonScore : the score for junctions with only a single read support
    """
    
    out = open(outputBed, "w")
    
    count = 0
    for line in open(inputBed):
        count += 1
        if line.startswith("track"):
            if count > 1:
                continue
            
            if newName != "":
                pieces = line.split()
                pieces.pop(1)
                pieces.insert(1, "name='%s'" % newName)
                pieces.pop(2)
                pieces.insert(2, "description='%s'" % newName) 
            newTrack = " ".join(pieces)
            out.write(newTrack)
            out.write("\n")
            continue
        
        pieces = line.split("\t")
        
        numCollapsed = 0
        if pieces[3].find("|junc=") > 0:
            numCollapsed = int(pieces[3].split("junc=")[-1])
            
        score = float(pieces[4])
        if numCollapsed < 2 and score <= scoreFilterSingle:
                continue
        elif score <= scoreFilterMultiple:
                continue
        
        out.write("\t".join(pieces))
        # if split on "\t" then "\n" still there.  otherwise need this.
        #out.write("\n")
        
def getCoveragePerBpFixed(wigFileName):
    """Given a wiggle file name, creates a dictionary with the data in the form
    dictionary[contig][position] = coverage at that position.
    Assumes the wiggle file is in fixedStep format (one column per data line)
    """
    
    currentContig = ""
    currentPosition = 0
    d = {}

    for line in open(wigFileName):
        if line.startswith("track") or line.startswith("browser"):
            continue
        
        if line.startswith("fixedStep"):
            currentContig = line.split()[1].split("=")[1]
            currentPosition = 1
            d[currentContig] = []
            
        else:
            d[currentContig].append(int(line))
            currentPosition = currentPosition + 1
            
    return d

def getCoveragePerBpVariable(wigFileName):
    """Given a wiggle file name, creates a dictionary with the data in the form
    dictionary[contig][position] = coverage at that position.
    Assumes the wiggle file is in variableStep format (two columns per data line).
    """
    
    d = {}

    for line in open(wigFileName):
        if line.startswith("track") or line.startswith("browser"):
            continue

        [chr, start, stop, level] = line.split()
        chr = chr.split("|")[1].replace("MAL", "chr")
        level = int(level)
        start = int(start)
        stop = int(stop)
        
        if not d.has_key(chr):
            d[chr] = {}
        
        for i in range(start, stop+1):
            d[chr][i] = level
            
    return d

def getCoverageLevel(inputBed, inputWig, outputTxt, wigTypeIsVariable=True):
    """For each gene (entry) in the bed file, the coverage level in the wig file
    is determined and written to the output."""
    
    if wigTypeIsVariable:
        cpbp = getCoveragePerBp2(inputWig)
    else:
        cpbp = getCoveragePerBpFixed(inputWig)
    
    out = open(outputTxt, "w")
    for line in open(inputBed):
        [chr, startStr, stopStr] = line.split()[:3]
        start = int(startStr)
        stop = int(stopStr)
        
        initialCov = cpbp[chr][start]
        endCov = cpbp[chr][stop]
        covList = []
        for i in range(start, stop):
            covList.append(cpbp[chr][i])
        avgCov = (initialCov + endCov) / 2.0
        
        out.write("\t".join([chr, startStr, stopStr, str(avgCov), str(max(covList))]))
        out.write("\n")
        
        
def removeMitoPlastid(inputBed, outputBed):
    """Removes all the results that do not have "chr" in their location titles."""
    out = open(outputBed, "w")
    
    for line in open(inputBed):
        pieces = line.split()
        
        # as long as we are here, remove mito/plastid
        if not pieces[0].find("chr") >= 0:
            continue
        
        out.write(line)


##################################################################################################
#  Other Utilities
##################################################################################################

def createHMMfile(outputHMMname):
    """Creates the 'initialHMM.txt' file used by HMMSplicer in the training step."""
    hmm = hmmWithQuality.HMM(["a", "j"], ["M", "X"], 5,
                             array( [ [ 0.5, 0.5], [ 0.0,  1.0 ] ]),
                             array( [[[ 0.5, 0.5],
                                      [ 0.5,  0.5]],
                                     [[ 0.5, 0.5],
                                      [ 0.5, 0.5]],
                                     [[ 0.5,  0.5],
                                      [ 0.5,  0.5]],
                                     [[ 0.5, 0.5],
                                      [ 0.5,  0.5]],
                                     [[ 0.5,  0.5],
                                      [ 0.5,  0.5]]]))
    hmm.dump()
    hmm.saveHMM(open(outputHMMname, "w"))
    
def writeHMMToFile(logfile, hmmFile):
    
    hmm = hmmWithQuality.HMM([], [], 0)
    hmm.loadHMM( open(hmmFile) )
    
    logfile.write(hmm.getDumpStr())
    logfile.flush()
    
def checkFastqReadLengths(inputFastq):
    
    sizes = {}
    
    #count = 0
    for (title, seq, qual) in FastqIterator(inputFastq):
        #count +=1 
        #if count > 50:
        #    break
        thisSize = len(seq)
        if not sizes.has_key(thisSize):
            sizes[thisSize] = 0
        sizes[thisSize] += 1
        
    for k, v in sizes.iteritems():
        print k, "=", v

def getMinMaxIntron(inputBed, percentBetweenLow, percentBetweenHigh):
    """Analyzes a genome annotation to determine basic intron information. 
    Prints the size of the largest and smallest intron, as well as the total number 
    of introns and the number of introns between the 'percentBetweenLow' and 'percentBetweenHigh'
    paramters.  Also prints the average intron size, average exon size, and the 
    average number of introns per gene.
    """
    
    minIntron = 1000000
    maxIntron = 0
    countTotal = 0
    countBetween = 0
    sizeList = []
    exonSizeList = []
    numGenes = 0
    for line in open(inputBed):
        if line.startswith("track"):
            continue
        
        numGenes += 1
        pieces = line.split()
        sizes = [int(x) for x in pieces[10].split(",")[:-1]]
        exonSizeList.extend(sizes)
        starts = [int(x) for x in pieces[11].split(",")[:-1]]
        for i in range(1, len(starts)):
            countTotal += 1
            intronSize = starts[i] - starts[i-1] - sizes[i-1]
            sizeList.append(intronSize)
            if intronSize >= percentBetweenLow and intronSize <= percentBetweenHigh:
                countBetween += 1
            minIntron = min(minIntron, intronSize)
            maxIntron = max(maxIntron, intronSize)
            
    print "The largest intron was %s and the smallest intron was %s" % (maxIntron, minIntron)
    print "There were %s introns total in %s genes, and %s (%s%%) were between %s and %s (inclusive)" % (countTotal, numGenes, countBetween, 
                                        (countBetween*100.0/countTotal), percentBetweenLow, percentBetweenHigh)
    print "The average intron size is %.2f" % ( float(sum(sizeList))/len(sizeList))
    print "Average number of introns per gene is %.2f" % ( float(countTotal) / numGenes) 

    print "Average exon size is %.2f" % ( float(sum(exonSizeList)) / len(exonSizeList))
    

def getFastaFromBed(inputBed, inputGenomeFasta, outputGeneFasta, chrToUse=None):
    """Given a genome annotation and the associated genome fasta file, this function
    creates as fasta file with the sequences of the genes in the input bed file.
    If chrToUse is specified then only the genes on chrToUse are saved."""
    gd = {}
    for (title, seq) in FastaIterator(inputGenomeFasta):
        #rd[title] = seq
        #print title
        if chrToUse and title == chrToUse:
            gd[chrToUse] = seq
            #print "     used"
            break
    
    out = open(outputGeneFasta, "w")
    for line in open(inputBed):
        if line.startswith("track"):
            continue
        
        pieces = line.split()
        chr = pieces[0]
        start = int(pieces[1])
        sizes = [int(x) for x in pieces[10].split(",")[:-1]]
        starts = [int(x) for x in pieces[11].split(",")[:-1]]
        seq = ""
        for i in range(len(starts)):
            seq = seq + gd[chr][start+starts[i]:start+starts[i]+sizes[i]]
        out.write(">%s\n" % (pieces[3]))
        for i in range(0, len(seq), 60):
            out.write(seq[i:i+60])
            out.write("\n")

def graphQualityPerPosition(inputFastq):
    """Creates a histogram counting the quality value per position in the input fastq file."""
    
    histD = {}
   
    count = 0
    for (titleStr, seqStr, qualityStr) in FastqIterator(inputFastq):
        count += 1
        #if count < 200000:
        #    continue
        #if count > 1200000:
        #    break
    
        qInts = convertQualityStr(qualityStr) 
        for i in range(len(qInts)):
            q = qInts[i]
            if q < 0 or q > 40:
                raise hmmErrors.InvalidQualityException("Invalid quality value %s at position %s of %s in %s" % (q, i, qualityStr, titleStr))
            
            if not histD.has_key(i):
                histD[i] = [0]*41
                
            histD[i][q] += 1
    
    print "Found %s records" % (count)
    print "Histogram of quality score per position"
    allk = histD.keys()
    allk.sort()
    for k in allk:
        print "%s|" % k, "|".join(str(x) for x in histD[k])
        
        
        
def bowtieToWig(bowtieFile, trackName, trackDescription, outputName):
    """Takes in a bowtie results file and creates a fixedStep wig file 
    (continuous data for coverage per bp in the genome)."""
    #print "Reading file in..."
    allCounts = _readBowtieInput(bowtieFile)
    #print "Chromosomes found:"
    #print allCounts.keys()
    _writeWiggle(trackName, trackDescription, allCounts, outputName)        
    
    
def _readBowtieInput(inputFile):
    """Helper method for bowtieToWig and reads the bowtie output into
    memory -- a dictionary of chromosomes, and for each chromosome, a 
    dictionary of bp indexes : coverage level.  
    i.e. allCounts["chr7"][12345] = 15 means chr7 at position 12345 is
    covered by 15 reads."""
    allCounts = {}
    count = 0
    countIgnored = 0
    for line in open(inputFile):
        count += 1
        #if count % 500000 == 0:
        #    print count

        # save up the number                                                                                    
        [read, strand, chrom, start, seq] = line.split("\t")[:5]         
        start = int(start)
        end = start + len(seq)  
        
        if not allCounts.has_key(chrom):
            allCounts[chrom] = {}
                                                                  
        
        for i in range(start, end):
            if allCounts[chrom].has_key(i):
                allCounts[chrom][i] += 1
            else:
                allCounts[chrom][i] = 1

    return allCounts
        
def _writeWiggle(trackName, trackDescription, allCounts, wigOut):
    """Writes the allCounts dictionary out to a wiggle file."""
    wigFile = open(wigOut, "w")
    wigFile.write("track type=wiggle_0 name='%s' description='%s' visibility=2\n" % (trackName,
                                                              trackDescription))
    for name in allCounts.keys():
        start = 0
        end = max(allCounts[name].keys())

        wigFile.write("fixedStep chrom=%s start=%s step=1 span=1\n" % (name, start+1))

        for i in range(start, end):
            if allCounts[name].has_key(i):
                curValue = allCounts[name][i]
            else:
                curValue = 0
            wigFile.write("%s\n" % (curValue))

    wigFile.close()
    
    








    
  



                

        
