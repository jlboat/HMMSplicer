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
# runHMM.py
#
# This script runs the HMM Splice Finder.
#
# Author: Michelle Dimon, May 2009
#

HELP_STRING = """
    HMM Splice Finder finds splice junctions in a set of sequence reads.  
    
    Options:
        -o    output directory
        -i    input fastq-format reads
        -g    input genome in fasta format
        
        -a    anchor size.  Default=8
        -x    delete the tmp directory.  Default=True
        -c    collapse similar junctions.  Default=True
        -t    temporary directory.  Default=<outputDirectory>/tmp
        -b    bowtie index of the genome.  Default=None (will be built from genome in the same folder as genome)
        -j    min intron size.  Default=5
        -k    max intron size.  Default=80000
        -p    number of processors.  Default=1
        -f    put junctions without canonical edges (GT-AG, etc) to a separate file.  Default=True
        -w    number of bases to move when trying to 'wiggle' to GT-AG.  Zero will turn off wiggling. Default=5.
        -m    score threshold for multiple (clustered) splice junctions.  Default=400
        -n    score threshold for single splice junctions.  Default=600
        -q    quality format.  p = Phred format.  s = old solexa format (pre1.3).  t = new solexa format (pipeline 1.3+).
        -d    include duplicates (put in a separate final bed file).  Default=False
        -e    num mismatches in the first 28 bp of the second half matching.  Default=2 (bowtie default)
        -r    large genome.  Required for large genomes (e.g. human) when running multiple processors.  Default=False
 """
 

import sys
import hmmErrors

# major version number checking
if sys.version_info[0] == (3):
    raise hmmErrors.InvalidPythonVersion("\n\nHMMSplicer does not support Python version 3.x yet.")

if sys.version_info[:2] < (2, 5):
    raise hmmErrors.InvalidPythonVersion("\n\nHMMSplicer does not support Python versions before 2.5.")

try:
    import numpy
except:
    raise hmmErrors.InvalidPythonVersion("\n\nUnable to import numpy.  Please download and install numpy 1.3.0 (http://sourceforge.net/projects/numpy/files/)")

import time
import commands
from getopt import getopt
import os
import glob

from threading import Thread

if sys.version_info[:2] >= (2, 6):
    import multiprocessing
    #from multiprocessing import Queue, Process, Pool
    #from Queue import Empty

import subprocess
import junctionHMM
import processJunctions
import hmmUtils
import configVals



# Defaults for the values which can be adjusted by input parameters
DEFAULT_ANCHOR_SIZE = 8
DEFAULT_MIN_INTRON_SIZE = 5
DEFAULT_MAX_INTRON_SIZE = 80000
DEFAULT_NUM_PROCESSORS = 1
DEFAULT_WIGGLE_AMOUNT = 5
DEFAULT_SCORE_FILTER_MULTIPLE = 400
DEFAULT_SCORE_FILTER_SINGLE = 600
PHRED_QUALITY = 0
SOLEXA_OLD_QUALITY = 1
SOLEXA_NEW_QUALITY = 2
COLORSPACE = 3

def main(argv=None):
    """The starting point for the code.  Parses the arguments, verifies file paths are valid, then
    calls the function to start the actual work."""
    
    if argv is None:
        argv = sys.argv

        
    outputDir = None
    inputFastq = None
    inputGenome = None
    
    anchorSize = DEFAULT_ANCHOR_SIZE
    doCollapse = True
    tmpDir = None
    inputBowtie = None
    inputBlat = None
    minIntron = DEFAULT_MIN_INTRON_SIZE
    maxIntron = DEFAULT_MAX_INTRON_SIZE
    numProcessors = DEFAULT_NUM_PROCESSORS
    separateNonCanonical = False
    wiggleAmount = DEFAULT_WIGGLE_AMOUNT
    filterToGtAg = True
    scoreFilterMultiple = DEFAULT_SCORE_FILTER_MULTIPLE
    scoreFilterSingle = DEFAULT_SCORE_FILTER_SINGLE
    deleteTmp = True
    includeDuplicates=False
    numMMinSeed = configVals.MISMATCHES_PER_HALF_BT_ANCHOR
    isLargeGenome = False
    qualityFormat = PHRED_QUALITY
         
    # get all the parameters from the user
    try:
        optlist, args = getopt(argv[1:], "ho:i:g:a:c:t:b:a:j:k:p:f:w:m:n:q:d:e:r:x:")
    except:
        print ""
        print HELP_STRING
        sys.exit(1)
        
    if len(optlist) == 0:
        print ""
        print HELP_STRING
        sys.exit(1)
        
    for (opt, opt_arg) in optlist:
        if opt == '-h':
            print ""
            print HELP_STRING
            sys.exit(1)
        elif opt == "-o":
            outputDir = opt_arg
        elif opt == "-i":
            inputFastq = opt_arg
        elif opt == "-g":
            inputGenome = opt_arg
            # make sure it ends in ".fa" or ".fasta" or ".fna"
        elif opt == "-c":
            doCollapse = opt_arg.startswith("T") or opt_arg.startswith("t")
        elif opt == "-t":
            tmpDir = opt_arg
        elif opt == "-b":
            inputBowtie = opt_arg
        elif opt == "-a":
            try:
                anchorSize = int(opt_arg)
            except ValueError:
                print "The Anchor Size (-a) value must be an integer."
                print
                print HELP_STRING
                sys.exit(1)
        elif opt == "-c":
            inputBlat = opt_arg
        elif opt == "-j":
            try:
                minIntron = int(opt_arg)
            except ValueError:
                print "The Minimum Intron (-i) value must be an integer."
                print
                print HELP_STRING
                sys.exit(1)
            if minIntron < 0:
                print "The Minimum Intron (-i) value cannot be negative."
                print
                print HELP_STRING
                sys.exit(1)               
        elif opt == "-k":
            try:
                maxIntron = int(opt_arg)
            except ValueError:
                print "The Maximum Intron (-j) value must be an integer."
                print
                print HELP_STRING
                sys.exit(1)
            if maxIntron < 0:
                print "The Maximum Intron (-j) value cannot be negative."
                print
                print HELP_STRING
                sys.exit(1)  
        elif opt == "-p":
            try:
                numProcessors = int(opt_arg)
            except ValueError:
                print "The number of Processors (-p) value must be an integer."
                print
                print HELP_STRING
                sys.exit(1)
            if numProcessors < 0:
                print "The number of Processors (-p) value cannot be negative."
                print
                print HELP_STRING
                sys.exit(1)  
        elif opt == "-f":
            filterToGtAg = opt_arg.startswith("T") or opt_arg.startswith("t")
        elif opt == "-w":
            try:
                wiggleAmount = int(opt_arg)
            except ValueError:
                print "The Wiggle Amount (-w) value must be an integer"
                print
                print HELP_STRING
                sys.exit(1)
            if wiggleAmount < 0:
                print "The Wiggle Amount (-w) value cannot be negative."
                print
                print HELP_STRING
                sys.exit(1)  
        elif opt == "-m":
            try:
                scoreFilterMultiple = int(opt_arg)
            except ValueError:
                print "The Score Filter for Multiple Splice Junctions (-m) value must be an integer"
                print
                print HELP_STRING
                sys.exit(1)
        elif opt == "-n":
            try:
                scoreFilterSingle = int(opt_arg)
            except ValueError:
                print "The Score Filter for Single Splice Junctions (-m) value must be an integer"
                print
                print HELP_STRING
                sys.exit(1)
        elif opt == "-e":
            try:
                numMMinSeed = int(opt_arg)
            except ValueError:
                print "The number of mismatches in the first 28 bp of the Second Half Match (-e) must be an integar"
                print
                print HELP_STRING
                sys.exit(1)
        elif opt == "-x":
            deleteTmp = opt_arg.startswith("T") or opt_arg.startswith("t")

        elif opt == "-q":
            if opt_arg.startswith("P") or opt_arg.startswith("p"):
                qualityFormat = PHRED_QUALITY
            elif opt_arg.startswith("S") or opt_arg.startswith("s"):
                qualityFormat = SOLEXA_OLD_QUALITY
            elif opt_arg.startswith("T") or opt_arg.startswith("t"):
                qualityFormat = SOLEXA_NEW_QUALITY
            elif opt_arg.startswith("C") or opt_arg.startswith("c"):
                qualityFormat = COLORSPACE
        elif opt == "-d":
            includeDuplicates = opt_arg.startswith("T") or opt_arg.startswith("t")
        elif opt == "-r":
            isLargeGenome = opt_arg.startswith("T") or opt_arg.startswith("t")
            
            
    if numProcessors > 1 and sys.version_info[:2] < (2, 6):
        raise hmmErrors.InvalidPythonVersion("\n\nTo run with multiple processors, HMMSplicer requires Python version 2.6.x or 2.7.x.  Please re-run HMMSplicer with either a single processor or an upgraded Python version.")
        
    # verify the required parameters were set    
    if outputDir == None:
        print "Output Directory (-o) is a required parameter"
        print
        print HELP_STRING
        sys.exit(1)
    elif inputFastq == None:
        print "Input Fastq (-i) is a required parameter"
        print
        print HELP_STRING
        sys.exit(1)
    elif inputGenome == None:
        print "Input Genome (-g) is a required parameter"
        print
        print HELP_STRING
        sys.exit(1)
    
    # put in defaults for path/genome variables
    if not outputDir.endswith("/"):
        outputDir += "/"
    if inputBowtie == None:
        inputBowtie = ".".join(inputGenome.split(".")[:-1])
    if inputBlat == None:
        inputBlat = ".".join(inputGenome.split(".")[:-1]) + ".2bit"
    if tmpDir == None:
        tmpDir = outputDir + "tmp/"
        
#    blatOOCname = ".".join(inputGenome.split(".")[:-1]) + (".%s.ooc" % (configVals.BLAT_TILE_SIZE))
    logfile = doSetup(outputDir, inputFastq, inputGenome, tmpDir, inputBowtie, qualityFormat == COLORSPACE)
    
    readLength = verifyQualityType(inputFastq, qualityFormat)
    
    originalFastq = inputFastq
    if qualityFormat == SOLEXA_OLD_QUALITY or qualityFormat == SOLEXA_NEW_QUALITY:
        printStatus("Converting to Phred quality format")
        log(logfile, "Converting to Phred quality format")
        newInputFastq = tmpDir + "convertedInput.fastq"
        convertFastqFile(inputFastq, newInputFastq, qualityFormat)
        inputFastq = newInputFastq
        
    logParams(logfile, outputDir, originalFastq, qualityFormat, inputGenome, doCollapse, tmpDir, inputBowtie,
           minIntron, maxIntron, numProcessors, wiggleAmount, filterToGtAg, 
           scoreFilterMultiple, scoreFilterSingle, includeDuplicates,
           readLength, numMMinSeed, isLargeGenome, anchorSize)
            
    runHMM(outputDir, inputFastq, inputGenome, doCollapse, tmpDir, inputBowtie,
           anchorSize, minIntron, maxIntron, numProcessors, wiggleAmount, filterToGtAg, 
           scoreFilterMultiple, scoreFilterSingle, includeDuplicates,
           readLength, numMMinSeed, isLargeGenome, logfile, qualityFormat, deleteTmp)
    
    # always remove the "part" files
    for f in glob.glob("%s*part*" % tmpDir):
        os.remove(f)
        
    # now delete the files as they are done being used
    # but leave this here just to take care of any that we missed
    if deleteTmp:
        printStatus("Deleting temporary files")
        log(logfile, "Deleting temporary files")
        # remove all the "dbf" files
        for f in glob.glob("%s*.dbf" % tmpDir):
            os.remove(f) 
        # remove all the "bowtie" files
        for f in glob.glob("%s*owtie*" % tmpDir):
            os.remove(f) 
        # remove all the "seed" files
        for f in glob.glob("%s*seed*" % tmpDir):
            os.remove(f) 
        if os.path.isfile("%sconvertedInput.fastq" % tmpDir):
            os.remove("%sconvertedInput.fastq" % tmpDir)
    

##################################################################################################
#  Set Up Functions (create output directory, build Bowtie index, etc)
##################################################################################################

def doSetup(outputDir, inputFastq, inputGenome, tmpDir, inputBowtie, isColorspace):
    """Verifies the files/directories exist, creates required genome indexes, and creates output
    folders."""

    verifyOneFile(inputFastq, "Input Fastq", "i")
    verifyOneFile(inputGenome, "Input Genome", "g")
    
    hadToCreate = False
    if not os.path.exists(outputDir):
        hadToCreate = True
        os.mkdir(outputDir)
        
    logfile = open(outputDir+"log.txt", "w")
    if hadToCreate:
        log(logfile, "Creating output directory %s" % (outputDir))
    if not os.path.exists(tmpDir):
        log(logfile, "Creating temporary directory %s" % (tmpDir))
        os.mkdir(tmpDir)
        
    if not os.path.isfile(inputBowtie+".1.ebwt"):
        buildBowtie(logfile, inputGenome, inputBowtie, isColorspace)
        
    return logfile
    
def verifyOneFile(fileToCheck, variableName, paramLetter):
    """Verifies that one file exists and errors out with nice message if the file doesn't exist."""
    if not os.path.isfile(fileToCheck):
        print "The file '%s' was not found." % (fileToCheck)
        print "Given as parameter %s (-%s)" % (variableName, paramLetter)
        print
        print HELP_STRING
        sys.exit(1)
        
def verifyQualityType(inputFastq, qualityFormat):
    """Verifies that the quality offset is correct."""
    # go through the first 50 quality values and check that the integers are in the expected range
    count = 0
    readLength = 0
    for (titleStr, seqStr, qualityStr) in hmmUtils.FastqIterator(inputFastq):
        #seqStr = seqStr.strip()
        #qualityStr = qualityStr.strip()
        count += 1
        if count > 50:
            break
        
        if readLength == 0:
            readLength = len(seqStr.strip())
        else:
            if readLength != len(seqStr.strip()):
                raise hmmErrors.InvalidFastq("The read sizes must all be the same")
        
        if qualityFormat != COLORSPACE:
            if qualityFormat == SOLEXA_OLD_QUALITY or qualityFormat == SOLEXA_NEW_QUALITY:
                qualInts = hmmUtils.convertQualityStr(qualityStr, 64)
            else:
                qualInts = hmmUtils.convertQualityStr(qualityStr, 33)
            
            for i in range(0, len(qualInts)):
                x = qualInts[i]
                if x > 40 or (x < 0 and qualityFormat != SOLEXA_OLD_QUALITY) or (x<-5 and qualityFormat == SOLEXA_OLD_QUALITY):
                    raise hmmErrors.InvalidQuality("Found an illegal quality value of %s in %s (quality string: %s, position: %s).\nAre you sure the quality flag is set right?" % (x, titleStr, qualityStr, i))

    return readLength
    
def buildBowtie(logfile, inputGenome, bowtiePrefix, isColorspace):
    """Builds the bowtie index from the inputGenome."""
    
    printStatus("Building bowtie index")
    log(logfile, "Building bowtie index")
    params = ""
    if isColorspace:
        params = "-C"
    cmd = "%s-build %s %s %s" % (configVals.PATH_TO_BOWTIE, params, inputGenome, bowtiePrefix)
    log(logfile, cmd)
    
    try:
        p = subprocess.Popen([cmd], shell=True)
        p.wait()
    except:
        raise   
    
    # Mac/Linux/Unix only:
    #(status, output) = commands.getstatusoutput(cmd)
    #if status != 0:
    #    print output
    #    raise hmmErrors.CommandLineException("Error building bowtie index with cmd %s" % (cmd))
    
def convertFastqFile(inputFastq, outputFastq, qualityFormat):

    if qualityFormat == SOLEXA_NEW_QUALITY:
        hmmUtils.convertFastqNewSolexaToPhred(inputFastq, outputFastq)
    elif qualityFormat == SOLEXA_OLD_QUALITY:
        hmmUtils.convertFastqOldSolexaToPhred(inputFastq, outputFastq)
    else:
        raise hmmErrors.UnexpectedException("Illegal quality format type %s in convertFastqFile" % qualityFormat)
    
    
##################################################################################################
#  Go!  (Actually run the HMMSplicer algorithm)
##################################################################################################
    
def runHMM(outputDir, inputFastq, inputGenome, doCollapse, tmpDir, inputBowtie, anchorSize,
           minIntron, maxIntron, numProcessors, wiggleAmount, filterToGtAg, 
           scoreFilterMultiple, scoreFilterSingle, includeDuplicates,
           readLength, numMMinSeed, isLargeGenome, logfile, qualityFormat, doDeleteTmp,
           genomeDict=None):
    """Runs the actual Splice-Finding pipeline.
        -- runs the bowtie
        -- runs blat on the reads that didn't match with bowtie
        -- makes the seeds
        -- trains the HMM on a subset of the data
        -- runs the HMM on all the data
        -- finds the second half of the match
        -- final processing (collapsing if required, etc)
    """
    
    testStr = ""

    # create all the file names 
    bowtieOutput = tmpDir + "bowtie.txt"
    bowtieNotFound = tmpDir + "notFoundBowtie.fastq"
    bowtieNotFoundFasta = tmpDir + "notFoundBowtie.fasta"
    bowtieHalfSeeds = tmpDir + "bowtieHalf.fastq"
    bowtieSeedOutput = tmpDir + "bowtieHalf.bowtie"
    blatOutput = tmpDir + "blat.txt"
    seedFile = tmpDir + "hmm.seed"
    seedSelectionFile = tmpDir + "hmm.training.seed"
    trainedHMM = tmpDir + "trained.hmm"
    noJunction = tmpDir + "hmm.noJunction.dbf"
    smallJunction = tmpDir + "hmm.smallJunction.dbf"
    initialJunction = tmpDir + "hmm.junction.dbf"
    multMatch = tmpDir + "hmm.multMatch.dbf"
    noMatch = tmpDir + "hmm.noMatch.dbf"
    junctionBed = tmpDir + "hmm.junction.initial.bed"
    junctionIndel = tmpDir + "hmm.shortIntron.bed"
    rescued = tmpDir + "hmm.rescued.bed"
    
    junctionBedNoDups = tmpDir + "junction"
    junctionBedNoDups2 = tmpDir + "junctionAgain"
    junctionBedWithDups = tmpDir + "junction_withDups"
    junctionBedWithDups2 = tmpDir + "junction_withDups2"
    junctionFiltered = outputDir + "junction.final.bed"
    junctionNonCannonical = outputDir + "junction.nonCanonical.bed"
    junctionFilteredWithDups = outputDir + "junction_withDups.final.bed"

    
    # do bowtie
    printStatus("Running bowtie")
    params = "-p %s --un %s" % (numProcessors, bowtieNotFound)
    if qualityFormat == COLORSPACE:
        params += " -C"
    #if solexaFormat:
    #    params += " --solexa-quals"
    cmd = "%s %s %s %s %s" % (configVals.PATH_TO_BOWTIE, params, inputBowtie, inputFastq, bowtieOutput)
    log(logfile, "Running bowtie with command: %s" % cmd)
    #print "ALERT!!! FIRST BOWTIME COMMENTED OUT!"
    try:
        p = subprocess.Popen([cmd], shell=True)
        p.wait()
    except:
        raise   
    #(status, output) = commands.getstatusoutput(cmd)
    #if status != 0:
    #    print output
    #    raise hmmErrors.CommandLineException("Error running bowtie with cmd %s" % (cmd))
    
    # do bowtie-half
    printStatus("Running bowtie-half")
    #print "DIVIDE READS COMMENTED OUT"
    if qualityFormat == COLORSPACE:
        hmmUtils.divideReadsColorspace(bowtieNotFound, bowtieHalfSeeds)
    else:
        hmmUtils.divideReads(bowtieNotFound, bowtieHalfSeeds, False)
    
    # delete as you go
    if doDeleteTmp:
        os.remove(bowtieNotFound)
        os.remove(bowtieOutput)
        if os.path.isfile("%sconvertedInput.fastq" % tmpDir):
            os.remove("%sconvertedInput.fastq" % tmpDir)
        
    params = "-n %s -m %s -k %s -p %s" % (numMMinSeed, configVals.MAX_SEED_REPEATS, 
                                          configVals.MAX_SEED_REPEATS, numProcessors)
    if qualityFormat == COLORSPACE:
        params += " -C --col-cqual --col-keepends"
    # no solexaFormat here -- we converted it in divideReads
    #if solexaFormat:
    #    params += " --solexa-quals"
    cmd = "%s %s %s %s %s" % (configVals.PATH_TO_BOWTIE, params, inputBowtie, bowtieHalfSeeds, bowtieSeedOutput)
    log(logfile, "Running bowtie for read-halves with command: %s" % cmd)
    #print "ALERT!!! SECOND BOWTIME COMMENTED OUT!"
    try:
        p = subprocess.Popen([cmd], shell=True)
        p.wait()
    except:
        raise   
    #(status, output) = commands.getstatusoutput(cmd)
    #if status != 0:
    #    print output
    #    raise hmmErrors.CommandLineException("Error running bowtie with cmd %s" % (cmd))
    
    # delete as you go
    if doDeleteTmp:
        os.remove(bowtieHalfSeeds)
    
    printStatus("Making genome dictionary")
    if not genomeDict:
        genomeDict = hmmUtils.fastaDictionary(inputGenome) 
    
    # make seeds
    printStatus("Making seeds")
    numSeeds = junctionHMM.seedReads(bowtieSeedOutput, genomeDict, seedFile, False, qualityFormat==COLORSPACE)
    log(logfile, "Found %s seeds" % numSeeds)
    
    # delete as you go
    if doDeleteTmp:
        os.remove(bowtieSeedOutput)
    
    if numSeeds < configVals.HMM_TRAINING_SIZE:
        #print "There were %s seeds, but %s are required for HMM training" % (numSeeds, HMM_TRAINING_SIZE)
        #sys.exit(1)
        printStatus("Training HMM with ALL reads")
        junctionHMM.trainHMM(seedFile, configVals.HMM_INITIAL_NAME, trainedHMM, numObs=numSeeds)
    else:
        # train HMM
        printStatus("Training HMM with subset")
        hmmUtils.randomlySelectFromFile(seedFile, seedSelectionFile, configVals.HMM_TRAINING_SIZE, numSeeds)
        junctionHMM.trainHMM(seedSelectionFile, configVals.HMM_INITIAL_NAME, trainedHMM, numObs=configVals.HMM_TRAINING_SIZE)
        
    log(logfile, "Trained HMM")
    hmmUtils.writeHMMToFile(logfile, trainedHMM)

    
    #run HMM and match second half
    if numProcessors == 1 or isLargeGenome:
        #don't bother with extra files if they aren't required
        printStatus("Running HMM")
        countNo, countSmall, countJunc = junctionHMM.runHMM(seedFile, trainedHMM, noJunction, smallJunction, initialJunction, anchorSize)
        log(logfile, "Ran HMM.  Found %s junctions, %s reads had small second pieces, and %s reads had no junction" % (countJunc, countSmall, countNo))
        printStatus("Matching second half")
        countNo, countMult, countJunc, countShort = junctionHMM.matchSecondHalf(initialJunction, genomeDict, noMatch, multMatch, 
                                    junctionBed, junctionIndel, minIntron, 
                                    maxIntron, wiggleAmount, readLength, anchorSize)
        log(logfile, "Matched Second half.  Found %s valid junctions, but %s did not match, %s matched multiple places, and %s had introns too short" % (countJunc,
                                                                            countNo, countMult, countShort))
        
    else:
        # break seed file into # processors sections.
        printStatus("Running HMM and matching second half")
        numSeedsPerFile = numSeeds / numProcessors
        seedFileName = tmpDir + "hmm.part%s.seed"
        noJunctionName = tmpDir + "hmm.noJunction.part%s.dbf"
        multJunctionName = tmpDir + "hmm.multMatch.part%s.dbf"
        smallJunctionName = tmpDir + "hmm.smallJunction.part%s.dbf"
        initialJunctionName = tmpDir + "hmm.junction.part%s.dbf"
        noMatchName = tmpDir + "hmm.noMatch.part%s.dbf"
        junctionBedName = tmpDir + "hmm.junction.initial.part%s.bed"
        junctionIndelName = tmpDir + "hmm.shortIntron.part%s.bed"
        
        f = open(seedFile)
        for i in range(numProcessors-1):
            out = open(seedFileName % i, "w")
            for j in range(numSeedsPerFile):
                out.write(f.readline())
            out.close()
            
        # last one -- read until the file is empty
        line = f.readline()
        out = open(seedFileName % (numProcessors-1), "w")
        while(len(line)>2):
            out.write(line)
            line = f.readline()
        out.close()
        
        # for each section, runHMM and matchSecondHalf in a new process
        #currtime = time.time()
        #po = Pool(processes=numProcessors)
        #for i in :
        #    print i
        #    po.apply_async(runAndMatch,((seedFileName%i, trainedHMM, noJunctionName%i, multJunctionName%i, 
        #                                 smallJunctionName%i, initialJunctionName%i, genomeDict, 
        #                                 noMatchName%i, junctionBedName%i, junctionIndelName%i, 
        #                                 minIntron, maxIntron, wiggleAmount, readLength)),callback=cb)
        #print "done with loop"
        #po.close()
        #po.join()
        #print "after close(), join()"
        count = 0
        countFinished = 0
        runningJobs = []
        outOfJobs = False
        keepGoing = True
        i = 0
            
        while keepGoing:
            # first deal with any completed jobs
            for job in runningJobs:
                if job.isDone():
                    countFinished += 1
                    jobNum = job.results
                    if countFinished % 1 == 0:
                        #listener.reportImportantInfo("At %s, %s genomes have completed.  Last one was %s long." % (
                        #                datetime.datetime.now().strftime("%m/%d/%Y at %H:%M:%S"), countFinished, len(job.viralGenome)))
                        log(logfile, "Job %s completed." % (jobNum))
                    runningJobs.remove(job)
            
            # then add enough jobs to bring up to numProcessors
            while len(runningJobs) < numProcessors and not outOfJobs:
                job = hmmSplicerJob(i, seedFileName%i, trainedHMM, noJunctionName%i, multJunctionName%i, 
                                         smallJunctionName%i, initialJunctionName%i, genomeDict, 
                                         noMatchName%i, junctionBedName%i, junctionIndelName%i, 
                                         minIntron, maxIntron, wiggleAmount, readLength, anchorSize)
                job.start()
                runningJobs.append(job)
                if i < (numProcessors-1):
                    i += 1
                else:
                    outOfJobs = True
                
            # don't stop until we are out of jobs AND all the running jobs have completed and been processed
            if len(runningJobs) == 0 and outOfJobs:
                keepGoing = False
            else:
                #if len(runningJobs) < numProcessors:
                #    print len(runningJobs)
                #    for i in len(runningJobs):
                #        print i, runningJobs[i].viralName, len(runningJobs[i].viralGenome)
                    
                # sleep a bit before repeating this loop
                time.sleep(10)
            
        # concatenate all junctionBed files
        out = open(junctionBed, "w")
        outSmall = open(smallJunction, "w")
        outMult = open(multMatch, "w")
        countJunc = 0
        countSmall = 0
        countMult = 0
        for i in range(numProcessors):
            for line in open(junctionBedName%i):
                if len(line) > 1:
                    countJunc += 1
                    out.write(line)
            for line in open(smallJunctionName%i):
                if len(line) > 1:
                    countSmall += 1
                    outSmall.write(line)
            for line in open(multJunctionName%i):
                if len(line) > 1:
                    countMult += 1
                    outMult.write(line)
        print "files combined"
                    
        log(logfile, "Ran HMM and matched second half.  There were %s junctions, %s reads with introns less than the minimum size and %s with multiple second matches" % (countJunc, 
                                                                                                    countSmall, countMult))
                    
        out.close()
        outSmall.close()
        outMult.close()
        
    # delete as you go
    if doDeleteTmp:
        for f in glob.glob("%s*part*" % tmpDir):
            os.remove(f)
        # remove all the "seed" files
        for f in glob.glob("%s*seed*" % tmpDir):
            os.remove(f) 

    if includeDuplicates:
        junctionHMM.removeDupReads(junctionBed, junctionBedNoDups+".bed", junctionBedWithDups+".bed")
    else:
        junctionHMM.removeDupReads(junctionBed, junctionBedNoDups+".bed")

    # rescue: first the small second side then the multiple matches                                      
    junctionHMM.rescueMultipleSeeds(smallJunction, junctionBedNoDups+".bed", rescued, 
                                    genomeDict, minIntron,     
                                    maxIntron, readLength, False)
    
    junctionHMM.rescueMultipleSeeds(multMatch, junctionBedNoDups+".bed", rescued, 
                                    genomeDict, minIntron,         
                                    maxIntron, readLength, True)                                                              

    # delete as you go
    if doDeleteTmp:
        # remove all the "dbf" files
        for f in glob.glob("%s*.dbf" % tmpDir):
            os.remove(f) 

    initial = open(junctionBedNoDups+".bed", "a")                                                                
    for line in open(rescued):                                                                          
        initial.write(line)                                                                            
    initial.close()  

    if includeDuplicates:
        junctionHMM.removeDupReads(junctionBedNoDups+".bed", junctionBedNoDups2+".bed", junctionBedWithDups2+".bed")
    else:
        junctionHMM.removeDupReads(junctionBedNoDups+".bed", junctionBedNoDups2+".bed")
    
    currentNameRoot = junctionBedNoDups2
    currentDupNameRoot = junctionBedWithDups2
    currentNonCanRoot = None
    description = "Splice Junctions"
    dupDesc = "Duplicate Splice Junctions"
        
    if filterToGtAg:
        printStatus("Filtering for GT-AG")
        nextName = currentNameRoot + ".can"
        nextNonCanName = currentNameRoot + ".nonCan"
        processJunctions.divideByGTAG(currentNameRoot+".bed", nextName+".bed", nextNonCanName+".bed", genomeDict)
        description += ", filtered to GT-AG"
        currentNameRoot = nextName
        currentNonCanRoot = nextNonCanName
        if includeDuplicates:
            nextDupName = currentDupNameRoot + ".gtag"
            processJunctions.divideByGTAG(currentDupNameRoot+".bed", nextDupName+".bed", 
                                          currentDupNameRoot+".notGtag.bed", genomeDict)
            dupDesc += ", filtered to GT-AG"
            currentDupNameRoot = nextDupName
       
    if doCollapse:
        printStatus("Collapsing junctions")
        nextName = currentNameRoot + ".collapsed"
        processJunctions.collapseCloseJunctions(currentNameRoot+".bed", nextName+".bed", 0)
        description += ", collapsed"
        currentNameRoot = nextName
        if currentNonCanRoot:
            nextNonCanName = currentNonCanRoot + ".collapsed"
            processJunctions.collapseCloseJunctions(currentNonCanRoot+".bed", nextNonCanName+".bed", 0)
            currentNonCanRoot = nextNonCanName
        if includeDuplicates:
            nextDupName = currentDupNameRoot + ".collapsed"
            processJunctions.collapseCloseJunctions(currentDupNameRoot+".bed", nextDupName+".bed", 0)
            dupDesc += ", collapsed"
            currentDupNameRoot = nextDupName
        
    hmmUtils.filterBedFile(currentNameRoot+".bed", junctionFiltered, scoreFilterSingle, scoreFilterMultiple, 
                           newName="HMMSplicer Junctions")
    hmmUtils.filterBedFile(currentNonCanRoot+".bed", junctionNonCannonical, scoreFilterSingle, 
                           scoreFilterMultiple, newName="HMMSplicer Non-cannonical Junctions")
    
    if includeDuplicates:
        hmmUtils.filterBedFile(currentDupNameRoot+".bed", junctionFilteredWithDups, scoreFilterSingle, scoreFilterMultiple,
                                newName=dupDesc)
    
    printStatus("Done.")


def printStatus(msg):
    print "%s: %s" % (time.strftime("%X %x"), msg)
    
def log(f, msg):
    f.write("%s: %s\n" % (time.strftime("%X %x"), msg))
    f.flush()

    
def logParams(logfile, outputDir, inputFastq, qualityFormat, inputGenome, doCollapse, tmpDir, inputBowtie,
           minIntron, maxIntron, numProcessors, wiggleAmount, filterToGtAg, 
           scoreFilterMultiple, scoreFilterSingle, includeDuplicates,
           readLength, numMMinSeed, isLargeGenome, anchorSize):
    """Write all the input parameters and config vals out to the log so we can have them later
    for debugging."""
    
    log(logfile, "Output directory: %s" % (outputDir))
    log(logfile, "Input fastq: %s" % (inputFastq))
    log(logfile, "Input quality format: %s" % (qualityFormat))
    log(logfile, "Input genome: %s" % (inputGenome))
    log(logfile, "Do collapse: %s" % (doCollapse))
    log(logfile, "Temporary directory: %s" % (tmpDir))
    log(logfile, "Input bowtie: %s" % (inputBowtie))
    log(logfile, "Minimum intron size: %s" % (minIntron))
    log(logfile, "Maximum intron size: %s" % (maxIntron))
    log(logfile, "Number of processors: %s" % (numProcessors))
    log(logfile, "Wiggle amount: %s" % (wiggleAmount))
    log(logfile, "Filter to canonical: %s" % (filterToGtAg))
    log(logfile, "Score threshold for multiple junctions: %s" % (scoreFilterMultiple))
    log(logfile, "Score threshold for single junctions: %s" % (scoreFilterSingle))
    log(logfile, "Include a duplicates file: %s" % (includeDuplicates))
    log(logfile, "Read length: %s" % (readLength))
    log(logfile, "Number of mismatches in the bowtie seed: %s" % (numMMinSeed))
    log(logfile, "Is large genome: %s" % (isLargeGenome))
    log(logfile, "Anchor size: %s" % (anchorSize))
    
    log(logfile, "Max seed repeats: %s" % (configVals.MAX_SEED_REPEATS))
    log(logfile, "Mismatches allowed in second piece: %s" % (configVals.SECOND_HALF_MM_ALLOWED))
    log(logfile, "Second piece alt intron: %s" % (configVals.SECOND_HALF_ALT_INTRON))
    log(logfile, "Minimum score difference: %s" % (configVals.MIN_SCORE_DIFF))
    log(logfile, "HMM Training Size: %s" % (configVals.HMM_TRAINING_SIZE))
    log(logfile, "Wiggle splice sites: %s" % (configVals.WIGGLE_SPLICE_SITES))
    log(logfile, "Filter splice sites: %s" % (configVals.SPLICE_SITES))
    
    
class hmmSplicerJob:
    """Modeled after Graham's AssemblyJob code for multiprocessing."""
    
    def __init__(self, jobNum, seedFile, trainedHMM, noJunction, multJunction, smallJunction, initialJunction, genomeDict, 
                noMatch, junctionBed, junctionIndel, minIntron, maxIntron, wiggleAmount, readLength, anchorSize): 
        self.jobNum = jobNum
        self.seedFile = seedFile
        self.trainedHMM = trainedHMM
        self.noJunction = noJunction
        self.multJunction = multJunction
        self.smallJunction = smallJunction
        self.initialJunction = initialJunction
        self.genomeDict = genomeDict
        self.noMatch = noMatch
        self.junctionBed = junctionBed
        self.junctionIndel = junctionIndel
        self.minIntron = minIntron
        self.maxIntron = maxIntron
        self.wiggleAmount = wiggleAmount
        self.readLength = readLength
        self.anchorSize = anchorSize
        
        self._started = False
        self._done = False
        
    def start(self):
        if not self._started:
            self._started = True
            self._donePipeIn, self._donePipeOut = multiprocessing.Pipe()
            self._resultsPipeIn, self._resultsPipeOut = multiprocessing.Pipe()
            self._job = multiprocessing.Process(target=self._run, args=(self._donePipeIn, self._resultsPipeIn,))
            self._job.start() 
            # for testing
            #self._run(self._donePipeIn, self._resultsPipeIn)
    
    def isDone(self):
        if not self._started:
            return False
        
        if not self._done:
            if not self._donePipeOut.poll():
                return False
            else:
                self._donePipeIn.close()
                self._donePipeOut.close()
                self.results = self._resultsPipeOut.recv()
                self._resultsPipeIn.close()
                self._resultsPipeOut.close()
                self._done = True
        
        return True
    
    def _run(self, donePipeIn, resultsPipeIn):
        #printStatus("Running HMM %s" % (self.jobNum))
        junctionHMM.runHMM(self.seedFile, self.trainedHMM, self.noJunction, self.smallJunction, self.initialJunction, self.anchorSize)
        #printStatus("Matching second half")
        junctionHMM.matchSecondHalf(self.initialJunction, self.genomeDict, self.noMatch, self.multJunction, 
                                    self.junctionBed, self.junctionIndel, self.minIntron, self. maxIntron, 
                                    self.wiggleAmount, self.readLength, self.anchorSize)
        #print "Done with runAndMatch"
                
        donePipeIn.send(True)
        resultsPipeIn.send(self.jobNum)

    
    
#def runAndMatch(seedFile, trainedHMM, noJunction, multJunction, smallJunction, initialJunction, genomeDict, 
#                noMatch, junctionBed, junctionIndel, minIntron, maxIntron, wiggleAmount, readLength):
#    printStatus("Running HMM")
#    junctionHMM.runHMM(seedFile, trainedHMM, noJunction, smallJunction, initialJunction)
#    printStatus("Matching second half")
#    junctionHMM.matchSecondHalf(initialJunction, genomeDict, noMatch, multJunction, junctionBed, junctionIndel, minIntron, 
#                               maxIntron, wiggleAmount, readLength)
#   print "Done with runAndMatch"
    
        
#def cb(r): #optional: callback function
#    pass
#    #print r





##############################################
if __name__ == "__main__":
    sys.exit(main(None)) 
