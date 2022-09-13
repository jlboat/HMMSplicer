# HMMSplicer Documentation

The code and documentation represented here have been updated by Lucas Boatwright for execution using Python 3. The original work was done by Dimon et al. 2010 as described below.

Table of contents:
1.  Requirements
2.  Quick setup
3.  Running test set
4.  Typical set up for an organism
5.  Custom configurations and parameters
6.  Usage Examples
7.  Result files
8.  Analyzing results
9.  Changes per version

The manuscript describing HMMSplicer can be found at http://www.ncbi.nlm.nih.gov/pubmed/21079731

Dimon MT, Sorber K, Derisi JL.  HMMSplicer: A Tool for Efficient and Sensitive Discovery of 
Known and Novel Splice Junctions in RNA-Seq Data.  PLoS One. 2010 Nov 8;5(11):e13875.

---

1.  Requirements
* Python 3
* numpy
* Bowtie

HMMSplicer was developed on Mac OSX.  Limited tests have been performed successfully on Linux, Unix, and Windows.

---

2.  Quick setup

* Install Python 3 and Numpy

* Install Bowtie

* Download HMMSplicer

* Open a terminal window and change into the directory you want to be the parent for the code:
      ```cd ~/software/```

* Clone HMMSplicer
      ```git clone https://github.com/jlboat/HMMSplicer```

* Change directory into the hmmSplicer directory 
      ```cd HMMSplicer```

* Set up hmmSplicer to be able to run bowtie.  Do one of the following:
    - option 1:  Edit HMMSplicer/configVals.py, variable "PATH_TO_BOWTIE" to the 
                     full path to bowtie (e.g. '/Users/mdimon/software/bowtie-0.12.3/bowtie'
                     for Mac or 'C:\software\bowtie-0.12.3\bowtie' for Windows)
    - option 2:  Add bowtie to your path.  For Mac, either copy (or sym-link) the
                     executable to /usr/local/bin or by add the full path to bowtie 
                     to your PATH (e.g. 
                     export PATH=$PATH:/Users/mdimon/software/bowtie-0.12.3).  For Windows,
                     you can add the full path to bowtie to your PATH environment variable.

* To run HMMSplicer, type "python runHMM.py".  Running without any options will print a help message with the required and optional parameters
      ```python runHMM.py```

* see "Typical Set up for an organism" (below) for instructions on preparing a genome for HMMSplicer.

---

3.  Running test set

HMMSplicer provides a small test dataset along with the P. falciparum genome (a relatively small 23 MB genome).  To test your installation, you can run HMMSplicer on this dataset.  The dataset is the first 25,000 reads from the SRX001454 dataset, which are randomly scattered across the P. falciparum transcriptome.

* cd into the HMMSplicer directory:
      ```cd ~/software/HMMSplicer/```

* run the test set:
      ```python runHMM.py -o pfTest/results/ -i pfTest/sraSRX001454_partial.fastq -g pfTest/pfGenome.fa -j 10 -k 1000```

* Check the results:
    HMMSplicer took 15 minutes to run with a single processor on the test machine (a Mac Pro).  The output is:
```
      15:24:39 03/01/10: Building bowtie index
 [... bowtie output related to building index ]
      15:25:37 03/01/10: Running bowtie
 [... bowtie output related to full-length alignment ]
      15:25:40 03/01/10: Running bowtie-half
 [... bowtie output related to read-half alignment ]
      15:25:45 03/01/10: Making genome dictionary
      15:25:46 03/01/10: Making seeds
      15:25:47 03/01/10: Training HMM with subset
      Converged in 29 iterations
      Trained Values:
      ================================================================================
      <hmmWithQuality.HMM instance at 0x1953e40>
      States:  2
      Observations:  2
      ------------------------------------------------------------------------
      State transition probabilities:
      [[ 0.94664406  0.05335594]
       [ 0.          1.        ]]
      ------------------------------------------------------------------------
      Observation probabilities:
      [[[ 0.65430797  0.27883193]
        [ 0.34569203  0.72116807]]

       [[ 0.71954824  0.34963291]
        [ 0.28045176  0.65036709]]

       [[ 0.70842223  0.33450779]
        [ 0.29157777  0.66549221]]

       [[ 0.70793786  0.34840234]
        [ 0.29206214  0.65159766]]

       [[ 0.92408317  0.3697206 ]
        [ 0.07591683  0.6302794 ]]]
      ------------------------------------------------------------------------
      15:39:10 03/01/10: Running HMM
      15:39:36 03/01/10: Matching second half
      15:39:37 03/01/10: Filtering for GT-AG
      15:39:37 03/01/10: Collapsing junctions
      15:39:37 03/01/10: Done.
```

The important files are junction.final.bed and junction.nonCanonical.bed.
* junction.final.bed:  The found junctions that match the canonical splice sites of GT-AG and GC-AG.  There are 36 junctions returned.
* junction.nonCanonical.bed: The found junctions that do not match GT-AG or GC-AG edges.  There are approximately 73 junctions returned.  

To examine the results in the UCSC Genome Browser:
* Go to the Malaria UCSC Genome Browser at http://areslab.ucsc.edu/cgi-bin/hgGateway
* Upload junction.final.bed (click the "manage custom tracks" button, then click the "add custom tracks" button on the next page.  Click the "browse" button and browse to junction.final.bed.  Pick it then click OK.)
* To view any particular junction, go to the genome browser (click the "go to genome browser" button).  Open the junction.final.bed file (in TextEdit or Notepad) and select the first 3 columns (i.e. "chr7	917534	917725") and paste this into the "position/search" box and click "jump".  Zoom out and adjust tracks as necessary.
-  Example:  The second junction "SRR005491.1542" at chr7	:396981-397136 shows up as novel in the genome browser, i.e. it does not align to any known introns.  This is because the UCSC Genome Browser is using older annotations, not the new version 6.3 of the P. falciparum annotations.  Going to PlasmoDB and looking up the closest gene, MAL7P1.208, reveals that this gene now is annotated with a new initial exon and three new trailing exons.  The junction discovered by HMMSplicer aligns with the middle trailing intron (the third intron in the gene).  PlasmoDB link: http://plasmodb.org/plasmo/showRecord.do?name=GeneRecordClasses.GeneRecordClass&source_id=MAL7P1.208&project_id=PlasmoDB

---

4.  Typical set up for an organism

* Download and prepare genome
HMMSplicer requires the genome in fasta format.  The initial part of the fasta title (up to the first whitespace character) will be used in the results output as the contig/chr location of the junction (the first column of the BED file).  Make sure the fasta title (up to the first whitespace) is unique.  Also, if you plan to upload to UCSC Genome Browser, make sure the fasta title is the name used by the UCSC Genome Browser for your organism and genome (usually 'chrX' where X is 1, 2, 3, X, etc).  For other tools, make sure the contig/chr in the fasta file is the one you wish to use for your BED analysis.  Also, if you don't download your genome from UCSC but you want to visualize the results against UCSC, make sure the genome version is exactly the same, otherwise your position values will be off.

The pre-built indexes on the Bowtie website do not contain the fasta sequence so they are not sufficient for running HMMSplicer. 

For organisms with UCSC Genome Browser setups, download the chromFa.tar.gz file.  For example, for humans, click the "ftp" link at the end of the organism overview, then click the bigZips folder, then the "chromFa.tar.gz" link.  Unzip and unarchive the results, and then concatenate all the chromosome fasta files together to create one file for the genome.

For large genomes such as human, the bowtie-build step can take several hours.  This step only has to be done once and HMMSplicer will initiate it.

* Run HMMSplicer
See examples below.

* Viewing results with UCSC genome browser
See description above for viewing the results of the test dataset.  Several other viewing options are available, depending on the organism.  The Broad Institute's Integrated Genome Viewer (IGV; http://www.broadinstitute.org/igv/) can view the BED files produced by HMMSplicer.

* Viewing results with Integrative Genome Viewer (IGV)
The Broad Institute's IGV is a great alternative to the UCSC Genome Browser, especially if your organism of interest is not hosted on the UCSC Genome Browser.  IGV can import BED files, so it can use the output of HMMSplicer without any modification.  IGV data is especially nice for viewing high throughput datasets at a 'zoomed in' level, where it shows the base for each read, colored by quality, making SNP confirmation very easy (this view requires SAM input, not BED format which does not include sequence or quality information).

---

5.  Custom configurations and parameters
Command Line Parameters:
```
        -o    output directory
        -i    input fastq-format reads
        -g    input genome in fasta format
        -a    anchor size.  Default=8
        -x    delete the tmp directory.  Default=True
        -c    collapse similar junctions.  Default=True
        -t    temporary directory.  Default=<outputDirectory>/tmp
        -b    bowtie index of the genome.  Default=None (will be built from genome in the 
                  same folder as genome)
        -j    min intron size.  Default=5
        -k    max intron size.  Default=80000
        -p    number of processors.  Default=1
        -f    put junctions without canonical edges (GT-AG, etc) to a separate file.  Default=True
        -w    number of bases to move when trying to 'wiggle' to GT-AG.  Zero will turn off 
                   wiggling. Default=5.
        -m    score threshold for multiple (clustered) splice junctions.  Default=400
        -n    score threshold for single splice junctions.  Default=600
        -q    quality format.  p = Phred format.  s = old solexa format (pre1.3).  t = new solexa format (pipeline 1.3+).
        -d    include duplicates (put in a separate final bed file).  Default=False
        -e    num mismatches in the first 28 bp of the second half matching.  Default=2 
                 (bowtie default)
        -r    large genome.  Required for large genomes (e.g. human) when running 
                  multiple processors.  Default=False
```

ConfigVals:
* ConfigVals.py contains configuration variables that are more rarely changed than those listed as command line parameters.  These variables are:
```
# the number mismatches allowed in bowtie when matching each read-half
MISMATCHES_PER_HALF_BT_ANCHOR = 2

# the maximum number of hits allowed for a read-half before all matches are dropped (the read-half is considered a repeat match)
MAX_SEED_REPEATS = 50

# the number of mismatches allowed in the second half for it to be considered a 'match'
SECOND_HALF_MM_ALLOWED = 3

# If there are two best matches for an intron and one is less than this amount long
# and the second is more than this many bp long then the first one is considered a 
# best match.  If both are less than this many bp long than the half 'cannot be matched'
# and no junction is selected, unless the read is able to be rescued by another splice
# junction at the same location.
#ECOND_HALF_ALT_INTRON = 1000
#
# the minimum difference between two scores for the junctions to be considered as having the
# 'same score'.  Two junctions in different places with the same score are considered duplicates
# and only included if the duplicate flag is on.  Two junctions in different places with different
# scores are not considered as duplicates and only the top scoring junction is retained.
MIN_SCORE_DIFF = 20

# the number of reads used in training the HMM
HMM_TRAINING_SIZE = 10000

# the splice sites to wiggle to.  The reverse complement will automatically be included also.
# sites will be searched in the order given.  For example, if the default wiggle is 5 bp and
# the junction can either wiggle to GT-AG that is 4 bp away or to GC-AG that is 1 bp away,
# the GT-AG site will be selected if GT-AG is listed before GC-AG in the list.
WIGGLE_SPLICE_SITES = ["GT-AG", "GC-AG", "AT-AC"]

# the splice sites considered 'good'.  Add both the splice junction and it's reverse complement ("GT-AG" and "CT-AC")
SPLICE_SITES = ["GT-AG", "CT-AC", "GC-AG", "CT-GC"]

# The program name you would type on the command line to run bowtie
PATH_TO_BOWTIE = "bowtie"

# The initial HMM values, before training
# If HMMSplicer is NOT run from the install folder, then this file name will have to be
# adjusted to the full path of the install folder.
HMM_INITIAL_NAME = "initialHMM.txt"
```

---

6.  Usage Examples

* Command line used for A. thaliana analysis in the manuscript:
```
python runHMM.py -o /Volumes/scratch/hmmSplicerTests/cress/al_Jan5 -g /Users/Shared/data/arabidopsis/allchr.fa -i /Volumes/mirror/SRX002554/allsix.fastq -j 5 -k 6000 -p 4 -d True
```

* Command line used for H. sapiens analysis in the manuscript:
```
python runHMM.py -o /Volumes/scratch/hmmSplicerTests/human/al_Jan27 -i /Volumes/mirror/SRX011550/both.fastq -g /Users/Shared/data/human/chromFa/allchr.unmasked.fa -j 5 -k 80000 -p 4 -d True -r True
```

* Command line used for P. falciparum analysis in the manuscript:
```
python runHMM.py -o /Volumes/scratch/hmmSplicerTests/pf/al_Jan27/ -j 5 -k 1000 -g /Users/Shared/data/plasmodium_downloads/falciparum/PfalciparumGenomic_PlasmoDB-6.0.chrNames.fa -i /Volumes/mirror/SRX001454/allfour.bcRemoved.fastq -p 4 -d True
```

* Running with shorter reads than recommended
The HMMSplicer algorithm requires that read-halves be able to align relatively uniquely within the genome.  The exact size requirement can vary depending on the size and complexity of the genome, but, as an example, reads less than 40 bp long do not perform as well when aligning to the human genome.  If your reads are relatively short (less than 45 bp for human, less than 40 bp for smaller genomes such as P. falciparum), I recommend trying HMMSplicer with the following adjustments to configVals.py:
```
MISMATCHES_PER_HALF_BT_ANCHOR = 1 
MAX_SEED_REPEATS = 25
```

---

7.  Result files

HMMSplicer produces the following output files:

* junction.final.bed:
    These are the final canonical junctions found by HMMSplicer.  The file is in BED format (http://genome.ucsc.edu/FAQ/FAQformat.html#format1).  When run with the default options, these junctions are collapsed (reads covering the same intron are collapsed into a single predicted junction and the score is increased accordingly) and filtered for canonical (GT-AG and GC-AG) splice sites.

* junction.nonCanonical.bed:
This file is identical to junction.final.bed except these are the junctions without canonical (GT-AG and GC-AG) splice sites.  

* tmp/junctionAgain.bed:
Initial results before dividing by canonical vs. non-canonical splice sites

* tmp/junctionAgain.can.bed and tmp/junctionAgain.nonCan.bed:
The result of dividing junctionAgain.bed according to canonical (default = GT-AG and GC-AG) vs. non-canonical splice sites.

* tmp/junctionAgain.can.collapsed.bed and tmp/junctionAgain.nonCan.collapsed.bed:
The result of collapsing junctionAgain.can.bed and junctionAgain.nonCan.bed, respectively.  This is the file to use for filterBedFile and generateROCcurve (see below) if you choose to do this further analysis.  The only difference between these files and junction.final.bed and junction.nonCanonical.bed is the filter by score.


As examples, the result sets from the manuscript are included in the HMMSplicer package.
For each organism:
* [organism]_results.zip : the compressed result files.  Each contains:
    - canonical.bed : The canonical (GT-AG and GC-AG) junctions above the score threshold (600 for singles, 400 for junctions with multiple instances)
    - noncanonical.bed : The junction edges without GT-AG or GC-AG that are above the score threshold
    - canonical.unfiltered.bed : The canonical junctions at all score thresholds (this is used to create the ROC plots)
    - noncanonical.unfiltered.bed : The noncanonical junctions at all score thresholds

If you want to run HMMSplicer and re-create one of these result sets, you will need to download the genome and raw FASTQ read file as follows:
A. thaliana:
    - The whole chromosomes can be found at: ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/
       These individual fasta files must be combined into a single fasta file to run HMMSplicer
    - The dataset can be found at:  http://www.ncbi.nlm.nih.gov/sra/SRX002554

P. falciparum:
    - The genome can be downloaded from PlasmoDB at:  http://plasmodb.org/common/downloads/release-6.0/Pfalciparum/
       (The version 6.0 was used for this analysis.  There are no changes between 6.0 and 6.3.  I haven't tested the latest 6.4 release yet,
but I don't anticipate any changes.  The mito and plastid contigs were not used for this analysis.)
       To be able to upload the data and view it against the Ares Lab Malaria Genome Browser (at http://areslab.ucsc.edu/index.html?org=P.+falciparum&db=pf5&hgsid=15702), the chromosome names must be edited to "chr1", "chr2", etc before HMMSplicer is run.
    - The dataset can be found at: http://www.ncbi.nlm.nih.gov/sra/SRP000384
       These sequences have an initial 2 nucleotide barcode which must be removed before the data can be analyzed with HMMSplicer.

H. sapiens:
     - The genome can be downloaded from UCSC Genome Browser at: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/
        The chromFa.tar.gz file was downloaded, then the individual chromosome files were combined into a single file and unmasked (lower case 
        repeats were converted back to upper case).
     - The dataset can be found at:  http://www.ncbi.nlm.nih.gov/sra/SRX011550

---

8.  Analyzing results
HMMSplicer comes with several scripts to help you analyze your data.  Some of the most helpful functions include:

* hmmUtils.fastaDictionary -- takes in a fasta file and returns a dictionary where the key is the title and the value is the sequence.  Takes a while (and a good amount of memory) for a large genome like human.

* hmmUtils.filterBedFile -- filters the bed file based on a score threshold.
The default results are filtered to a score of 600 for single junctions and a score of 400 for multiple junctions.  Where you want to set the score threshold depends on your tolerance for false positives, based on your experiment goals, as well as on the parameters of your run (scores are based partially on quality so runs with very low quality may want to adjust the score threshold downwards.  Also, if you have extremely high coverage, you may want to adjust your scores upwards).  'filterBedFile' in conjunction with 'measureSpecificity' is probably the best way to get a feel for the accuracy of your data, with the caveat that the specificity returned is only as good as the annotation used 

* processJunctions.readJunctionsFromBed -- takes in a bed file of junctions (i.e. the results or a set of EST annotations) and reads them into a dictionary

* scoreJunctions.measureSpecificity -- takes in a results bed file and a bed file of ESTs (or RefSeq or whatever genome annotations you have) and returns the number of junctions in the results file that align with known junctions.  
I've done testing on human, arabidopsis and Plasmodium falciparum, and generally about 80-95% of the junctions match known annotations, depending on the score threshold,  accuracy of the annotations used, and the extent of the annotations used (human % is low if I only include RefSeq but goes up dramatically if I include all EST evidence also).

* scoreJunctions.measureSpecificityUsingDict -- same as above but takes in the dictionary (from readJunctionsFromBed) instead of re-creating the EST dictionary every time.  
Use this if you are going to be calling this function over and over to save time.  For example, I'll often do repeats of filterBedFile / measureSpecificityUsingDict to see how changing the score threshold affects the percent of junctions matching known introns.

* scoreJunctions.vennDiagram -- takes two bed files (for example, results from two experiments) and counts the number in only the first file, only the second file or in both.  Optionally can output each set to a bed file for further analysis.

* scoreJunctions.generateROCcurve -- generates a ROC curve using the junction scores.  
You want to use this on the full result set.  Gives an idea of accuracy of the scoring algorithm (with the caveat that for real results, as opposed to simulated results, the ability to accurately decipher a 'false' junction is limited so this ROC curve will be an underestimate of the true predictive power of the score.)

* scoreHistogram.scoreHistogram : In subsequent analysis, we found that creating a score histogram for known vs. novel splice junctions is a great way to visualize the correct score threshold.  The scoreHistogram.py module contains the scoreHistogram function which prints out the data you need to generate this type of graph.  For smaller datasets, adjusting the score threshold isn't as necessary, but when you start pooling a lot of RNA-Seq data together you will probably want to adjust the score thresholds higher.


Details on how to use these scripts:
Python comes with a built-in interpreter that I tend to use for my analysis so the scripts are set up to take advantage of this functionality.  (If you've used R to analyze data, it is comparable to the R interface.)  You can either use the interpreter to run these scripts, or you can write your own python wrapper that can run a set of analysis.  Which you prefer depends on how you like to work.  I like to use the python interpreter to play around with the data and get a feel for what analysis we want to do. On the other hand, writing python wrappers enables you to keep very meticulous records of exactly what analysis you performed, especially if you have good file names (with dates) on the wrapper scripts.  I often do this for the final analysis after I'm done "playing around" when I want to generate the publication-ready datasets.

Interpreter Specifics:
---
The nice thing about the interpreter, especially for larger genomes such as mouse and human, is that you can load in the data structures once and then play around with different settings.  For example, the first part of the example for the interpreter demonstrates how you can use different score thresholds to affect specificity and sensitivity in your HMMSplicer results.  It can take a long time to load in the library of known splice variants (the human one took me about 10-20 minutes, as I recall), but then each filtering / measuring command is very fast.

Here's a transcript of an example interpreter session:
```python
# Python 2.6.2 (r262:71600, Apr 16 2009, 09:17:39) 
# [GCC 4.0.1 (Apple Computer, Inc. build 5250)] on darwin
# Type "help", "copyright", "credits" or "license" for more information.
import hmmUtils
import processJunctions
import scoreJunctions

pfKnown = processJunctions.readJunctionsFromBed("/Users/Shared/data/plasmodium_downloads/falciparum/pdbAndEsts.bed")

scoreJunctions.measureSpecificityUsingDict("/Volumes/backup2/hsExample/junction.final.bed", pfKnown, wiggle=3)
# 3018 overlapped but 444 did not.  87% overlapped
# 87.175043327556324

hmmUtils.filterBedFile("/Volumes/backup2/hsExample/junction.final.bed", "/Volumes/backup2/hsExample/junction.filter700.bed", 800, 700, "filter 700")
scoreJunctions.measureSpecificityUsingDict("/Volumes/backup2/hsExample/junction.filter700.bed", pfKnown, wiggle=3)
# 2640 overlapped but 275 did not.  90% overlapped
# 90.566037735849051
hmmUtils.filterBedFile("/Volumes/backup2/hsExample/junction.final.bed", "/Volumes/backup2/hsExample/junction.filter800.bed", 800, 800, "filter 800")
scoreJunctions.measureSpecificityUsingDict("/Volumes/backup2/hsExample/junction.filter800.bed", pfKnown, wiggle=3)
# 2564 overlapped but 256 did not.  90% overlapped
# 90.921985815602838
hmmUtils.filterBedFile("/Volumes/backup2/hsExample/junction.final.bed", "/Volumes/backup2/hsExample/junction.filter900.bed", 900, 900, "filter 900")
scoreJunctions.measureSpecificityUsingDict("/Volumes/backup2/hsExample/junction.filter900.bed", pfKnown, wiggle=3)
# 2064 overlapped but 157 did not.  92% overlapped
# 92.93111211166142
hmmUtils.filterBedFile("/Volumes/backup2/hsExample/junction.final.bed", "/Volumes/backup2/hsExample/junction.filter1000.bed", 1000, 1000, "filter 1000") 
scoreJunctions.measureSpecificityUsingDict("/Volumes/backup2/hsExample/junction.filter1000.bed", pfKnown, wiggle=3)
# 1416 overlapped but 65 did not.  95% overlapped
# 95.611073598919646

scoreJunctions.measureSpecificityUsingDict("/Volumes/backup2/hsExample/junction.filter800.bed", pfKnown, wiggle=3, goodFile="/Volumes/backup2/hsExample/junction.filter800.matchesKnown.bed", badFile="/Volumes/backup2/hsExample/junction.filter800.novel.bed")
# 2564 overlapped but 256 did not.  90% overlapped
# 90.921985815602838
```

Wrapper Specifics:
---
For the wrapper example, I'm going to assume we've decided on a score threshold and the goal is to print out some statistics and generate a set of files to upload to the UCSC genome browser.  

Set up for wrapperScript.py:
```
# BEFORE doing anything -- python can't find hmmUtils
echo $PYTHONPATH
# /Users/mdimon/src:/Users/mdimon/src/lab
python wrapperScript.py 
# Traceback (most recent call last):
#   File "wrapperScript.py", line 1, in <module>
#     import hmmUtils
# ImportError: No module named hmmUtils
# ADD hmmSplicer to your python path
export PYTHONPATH=$PYTHONPATH:/Users/mdimon/src/hmmSplicer
```

Here's what the wrapper.py script looks like:

```python
import hmmUtils
import processJunctions
import scoreJunctions
import scoreHistogram
# read in the known EST set
pfKnown = processJunctions.readJunctionsFromBed("/Users/Shared/data/plasmodium_downloads/falciparum/pdbAndEsts.bed")
# Assume we've already decided we want to use 800 as the score threshold
# generate the new filtered file
hmmUtils.filterBedFile("/Volumes/backup2/hsExample/junction.final.bed", "/Volumes/backup2/hsExample/junction.filter800.bed", 800, 800, "filter 800")
# measure the specificity at this level and generate the known and novel files
print("Specificity for score filter 800")
scoreJunctions.measureSpecificityUsingDict("/Volumes/backup2/hsExample/junction.filter800.bed", pfKnown, wiggle=3, goodFile="/Volumes/backup2/hsExample/junction.filter800.matchesKnown.bed", badFile="/Volumes/backup2/hsExample/junction.filter800.novel.bed")
# generate the score histograms to view the data
print()
print("Score Histogram for filter800 for KNOWN junctions")
scoreHistogram.scoreHistogram("/Volumes/backup2/hsExample/junction.filter800.matchesKnown.bed")
print()
print("Score Histogram for filter800 for NOVEL junctions")
scoreHistogram.scoreHistogram("/Volumes/backup2/hsExample/junction.filter800.novel.bed")
```

The files are generated as expected and the analysis information (print statements) get printed to the command line
 (the output skips the scoreHistogram results because they are long):

```
python wrapperScript.py
# Specificity for score filter 800
# 2564 overlapped but 256 did not.  90% overlapped
```

---

9.  Changes per version

9/12/2022: Version 1.0.0
* Conversion of scripts for Python 3

11/23/2010: Version 0.9.5
* If the user chooses to delete temporary files then they will be deleted as soon as they are no longer needed, instead of waiting until the end of the run.
* Improved description of how to use the extra script tools (i.e. section 8 above)
* Included the scoreHistograms.py module
* Bug fix in scoreJunctions.measureSpecificity to correctly create BED files for junctions that match and don't match known junctions.
* Improved the handling of multiple processor jobs (without large genomes)

