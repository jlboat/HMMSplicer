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
#
# This file defines constants used in the HMMSplicer.

# the number mismatches allowed in bowtie when matching each read-half
MISMATCHES_PER_HALF_BT_ANCHOR = 2

# the maximum number of hits allowed for a read-half before all matches are dropped (the read-half is considered a repeat match)
MAX_SEED_REPEATS = 50

# the minimum exact size required when finding the location of the second part of the read
# This parameter changed frequently enough that it was moved to the command line input parameters.
#ANCHOR_SIZE = 6
#ANCHOR_SIZE = 8

# the number of mismatches allowed in the second half for it to be considered a 'match'
SECOND_HALF_MM_ALLOWED = 3

# If there are two best matches for an intron and one is less than this amount long
# and the second is more than this many bp long then the first one is considered a 
# best match.  If both are less than this many bp long than the half 'cannot be matched'
# and no junction is selected, unless the read is able to be rescued by another splice
# junction at the same location.
SECOND_HALF_ALT_INTRON = 1000

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
#SPLICE_SITES = ["GT-AG", "CT-AC", "GC-AG", "CT-GC", "AT-AC", "GT-AT"]
SPLICE_SITES = ["GT-AG", "CT-AC", "GC-AG", "CT-GC"]
#SPLICE_SITES = ["GT-AG", "CT-AC"]

# The program name you would type on the command line to run bowtie
PATH_TO_BOWTIE = "bowtie"

# The initial HMM values, before training
# If HMMSplicer is NOT run from the install folder, then this file name will have to be
# adjusted to the full path of the install folder.
HMM_INITIAL_NAME = "initialHMM.txt"

