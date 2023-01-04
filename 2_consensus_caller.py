#!/usr/bin/env python3

"""UnifiedConsensusMaker.py
Version 3.1
May 29, 2019
Programs by Scott Kennedy(1)
V3.1 Update by Brendan Kohrn(1)
(1) Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195
  
Changelog:
V3.1: Add loclen option to accommodate dumis, fix last-read skip issue, 
      and update to python3.  Also, add docstring to program

This program is intended to convert raw sequencing reads from Duplex 
Sequencing libraries into PCRcs and cDNAcs sequences for later analysis.  
It has three parts: first, move barcode sequences to the read name, 
second, sort based on read names, and third, make consensuses based on 
those barcode sequences.  

Settings considerations:
--input is a paired-end, unaligned bam file, such as from PicardTools 
    FastqToSam
--taglen, --spacerlen, and --loclen should be set based on the library 
    preperation protocol used.  
        If dumis are being used, --loclen must 
    be >0.  
        For the method in Nature Protocols, --taglen is 12, 
    --spacerlen is 5, and --loclen can be 0.  

After this, further steps might include:
    Aligning (with BWA, etc)
    Read trimming (fixed-length trimming and/or overlap trimming)
    Local realignment (such as with IndelRealigner from GATK 3.8)
    Variant calling

When doing variant calling, remember that, for Duplex Sequencing, 
    samples should not be considered diploid
"""

import datetime
import gzip
import multiprocessing as mp
import os
import random
import sys
from argparse import ArgumentParser
from collections import defaultdict

import pysam

class iteratorWrapper:
    def __init__(self, inIterator, finalValue):
        self.it = inIterator
        self.finalValue = finalValue
        self.endIter = False
    def __iter__(self):
        return(self)
    def __next__(self):
        try:
            temp = next(self.it)
        except StopIteration:
            if self.endIter == False:
                temp = self.finalValue
                self.endIter = True
            else:
                raise(StopIteration)
        return(temp)
    next = __next__


def consensus_caller(input_reads, cutoff, tag, length_check):

    nuc_identity_list = [0, 0, 0, 0, 0, 0]  
    # In the order of T, C, G, A, N, Total
    nuc_key_dict = {0: 'T', 1: 'C', 2: 'G', 3: 'A', 4: 'N'}
    consensus_seq = ''

    if length_check is True:

        for read in input_reads[1:]:
            if len(read) != len(input_reads[0]):
                raise Exception((f"Read lengths for tag {tag} used for "
                                 f"calculating the PCRcs are not uniform!!!"
                                 ))

    for i in range(len(input_reads[0])):  
        # Count the types of nucleotides at a position in a read.
        # i is the nucleotide index within a read in groupedReadsList
        for j in range(len(input_reads)):  
        # Do this for every read that comprises a tag family.
        # j is the read index within groupedReadsList
            try:
                if input_reads[j][i] == 'T':
                    nuc_identity_list[0] += 1
                elif input_reads[j][i] == 'C':
                    nuc_identity_list[1] += 1
                elif input_reads[j][i] == 'G':
                    nuc_identity_list[2] += 1
                elif input_reads[j][i] == 'A':
                    nuc_identity_list[3] += 1
                elif input_reads[j][i] == 'N':
                    nuc_identity_list[4] += 1
                else:
                    nuc_identity_list[4] += 1
                nuc_identity_list[5] += 1
            except Exception:
                break
        try:
            for j in [0, 1, 2, 3, 4]:
                if (float(nuc_identity_list[j])
                        /float(nuc_identity_list[5])
                        ) >= cutoff:
                    consensus_seq += nuc_key_dict[j]
                    break
                elif j == 4:
                    consensus_seq += 'N'
        except Exception:
            consensus_seq += 'N'
        nuc_identity_list = [0, 0, 0, 0, 0, 0]  
        # Reset for the next nucleotide position

    return consensus_seq


def qual_calc(qual_list):
    return [sum(qual_score) for qual_score in zip(*qual_list)]

def main():
    startTime = datetime.datetime.now()
    parser = ArgumentParser()
    parser.add_argument(
        '--input', 
        dest = 'in_bam', 
        required = True,
        help = 'Path to unaligned, paired-end, bam file.'
        )
    parser.add_argument(
        '--taglen', 
        dest = 'tag_len', 
        type = int, 
        default = 16+11,
        help = 'Length in bases of the duplex tag sequence.[16+11]'
        )
    parser.add_argument(
        '--spacerlen', 
        dest = 'spcr_len', 
        type = int, 
        default = 0,
        help = (f'Length in bases of the spacer sequence between'
                f'duplex tag and the start of target DNA. [0]'
                )
        )
    parser.add_argument(
        "--loclen", 
        dest = 'loc_len', 
        type = int, 
        default = 0, 
        action = "store",
        help = (f"Number of base pairs to add to barcode for location "
                f"specificity.  Bases are not removed from read.  [0]"
                )
        )
    parser.add_argument(
        "--tagstats", 
        dest = 'tagstats', 
        action = "store_true",
        help = "Output tagstats file"
        )
    parser.add_argument(
        '--minmem1', 
        dest = 'minmem1', 
        type = int, 
        default = 1,
        help = "Minimum number of reads allowed to comprise a PCR consensus. [1]"
        )
    parser.add_argument(
        '--maxmem1', 
        dest = 'maxmem1', 
        type = int, 
        default = 1000,
        help = "Maximum number of reads allowed to comprise a PCR consensus. [1000]"
                        )
    parser.add_argument(
        '--cutoff1', 
        dest = 'cutoff1', 
        type = float, 
        default = .99,
        help = (f"Percentage of nucleotides at a given position "
                f"in a read that must be identical in order "
                f"for a PCR consensus to be called at that position. "
                f"[0.99]"
                )
        )
    parser.add_argument(
        '--minmem2', 
        dest = 'minmem2', 
        type = int, 
        default = 3,
        help = "Minimum number of reads allowed to comprise a cDNA consensus. [3]"
        )
    parser.add_argument(
        '--maxmem2', 
        dest = 'maxmem2', 
        type = int, 
        default = 1000,
        help = "Maximum number of reads allowed to comprise a cDNA consensus. [1000]"
                        )
    parser.add_argument(
        '--cutoff2', 
        dest = 'cutoff2', 
        type = float, 
        default = .7,
        help = (f"Percentage of nucleotides at a given position "
                f"in a read that must be identical in order "
                f"for a cDNA consensus to be called at that position. "
                f"[0.7]"
                )
        )
    parser.add_argument(
        '--Ncutoff', 
        dest = 'Ncutoff', 
        type = float, 
        default = .9,
        help = (f"With --filt 'n', maximum fraction of Ns allowed in a "
                f"consensus [0.9]"
                )
        )
    parser.add_argument(
        '--write-PCRcs', 
        dest = 'write_PCRcs', 
        action = "store_true",
        help = "Print the PCRcs reads to file in FASTQ format"
        )
    parser.add_argument(
        '--without-cDNAcs', 
        dest = 'without_cDNAcs', 
        action = "store_true",
        help = "Don't print final cDNAcs reads"
        )
    parser.add_argument(
        "--rep_filt", 
        action = "store",  
        type = int, 
        dest = 'rep_filt',
        default = 9,
        help = (f"Remove tags with homomeric runs of nucleotides of length "
                f"x. [9]"
                )
        )
    parser.add_argument(
        '--prefix', 
        action = "store",
        dest = 'prefix', 
        type = str, 
        required = True,
        help = "Sample name to uniquely identify samples"
        )
    parser.add_argument(
        '--numCores', 
        action = "store",
        dest = "cores", 
        type = int, 
        default = 1, 
        help = "Number of cores to use for sorting UMI-processed reads."
        )
    parser.add_argument(
        '--numAlignmentReads', 
        dest = "numAlignReads", 
        action = "store",
        type = int, 
        default = 500000, 
        help = "Number of read pairs to output as fastq to align for determining raw reads on target.  Set to 0 to skip this output.  Will stop when it reaches the end of the file or this number."
        )
    o = parser.parse_args()
    
    # adjust number of cores
    if o.cores >= mp.cpu_count() and o.cores > 1:
        o.cores = mp.cpu_count() - 1
    
    dummy_header = {'HD': {'VN': '1.0'}, 
                    'SQ': [{'LN': 1575, 'SN': 'chr1'}, 
                           {'LN': 1584, 'SN': 'chr2'}
                           ]
                    }
    in_bam_file = pysam.AlignmentFile(o.in_bam, "rb", check_sq=False)
    temp_bam = pysam.AlignmentFile(f"{o.prefix}.temp.bam", 
                                   'wb', 
                                   header=dummy_header
                                   )
    # Initialize Counters:
    # Counter for raw reads that were UMI-processed
    read_count = 0
    #Counter for number of read pairs with matching UMIs
    read_pair_count = 0	
    # Counter for number of families
    familyCtr = 0
    # Counter for cDNAcs UMIs with bad UMIs
    nUMIs = 0
    monoNtUMIs = 0
    # Counter for reads processed
    readsCtr = 0
    # Counter for low familiy size families
    smallFamilySize = 0
    # Counter for unrepresented families
    zeroFamilySize = 0
    # Counter for high-N PCRcs filtered
    highN_PCRcs = 0
    # Counter for number of PCRcs made
    numPCRcs = 0
    # counter for number of families that fail to find their partner
    failedcDNAcs = 0
    # Counter for number of cDNAcs made
    numcDNAcs = 0
    # Counter for number of high-N cDNAcs filtered
    highN_cDNAcs = 0
    # Counter for final number of cDNAcs reads passing all filters
    num_filt_cDNAcs = 0

    # Open Files
    if o.write_PCRcs is True:
        read1_PCRcs_fq_file = gzip.open(f"{o.prefix}_read1_PCRcs.fq.gz", 'wt')
        read2_PCRcs_fq_file = gzip.open(f"{o.prefix}_read2_PCRcs.fq.gz", 'wt')

    if o.without_cDNAcs is False:
        read1_cDNAcs_fq_file = gzip.open(f"{o.prefix}_read1_cDNAcs.fq.gz", 'wt')
        read2_cDNAcs_fq_file = gzip.open(f"{o.prefix}_read2_cDNAcs.fq.gz", 'wt')

    # This block of code takes an unaligned bam file, extracts the tag 
    # sequences from the reads, and converts them to to "ab/ba" format 
    # where 'a' and 'b' are the tag sequences from Read 1 and Read 2, 
    # respectively. Conversion occurs by putting the tag with the "lesser" 
    # value in front of the tag with the "higher" value. The original 
    # tag orientation is denoted by appending #ab or #ba to the end of 
    # the tag. After conversion, the resulting temporary bam file is then
    # sorted by read name.
    
    print("Parsing UMIs...")
    alignedReadCount = 0
    if o.numAlignReads != 0:
        fAlign1 = gzip.open(f"{o.prefix}_aln_seq1.fq.gz", 'wt')
        fAlign2 = gzip.open(f"{o.prefix}_aln_seq2.fq.gz", 'wt')
    
    for line in in_bam_file.fetch(until_eof=True):

        if read_count % 2 == 0:

            temp_read1_entry = pysam.AlignedSegment()
            temp_read1_entry.query_name = line.query_name
            temp_read1_entry.query_sequence = line.query_alignment_sequence
            temp_read1_entry.query_qualities = line.query_alignment_qualities

            if line.has_tag("RX"):
                temp_read1_entry.set_tag("RX", line.get_tag("RX"), 'Z')
            else:
                raise IOError("RX tag not present in read %s" % read_count)

        if read_count % 2 == 1:

            temp_bam_entry = pysam.AlignedSegment()
            if line.has_tag("RX"):
                if temp_read1_entry.get_tag("RX") != line.get_tag("RX"):
                    raise IOError("RX tags for read 1 and read 2 different at line %s" % read_count)
#New:
                if temp_read1_entry.get_tag("RX") == line.get_tag("RX"):
                    read_pair_count += 1
            else:
                raise IOError("RX tag not present in read %s" % read_count)
            
            
            #~ print(temp_read1_entry.get_tag("RX"))
            #~ print(line.get_tag("RX"))
            #~ exit()
            temp_bam_entry.query_name = temp_read1_entry.get_tag("RX").split('-')[0] + "-" + temp_read1_entry.get_tag("RX").split('-')[1] + '#'
            
            # Write entries for Read 1
            temp_bam_entry.query_name += ":1"
            temp_bam_entry.query_sequence = temp_read1_entry.query_sequence
            temp_bam_entry.query_qualities = temp_read1_entry.query_qualities
            temp_bam_entry.set_tag('X?', temp_read1_entry.query_name, 'Z')
            temp_bam.write(temp_bam_entry)

            # Write entries for Read 2
            temp_bam_entry.query_name = temp_bam_entry.query_name.replace('1', '2')
            temp_bam_entry.query_sequence = line.query_sequence
            temp_bam_entry.query_qualities = line.query_qualities
            temp_bam_entry.set_tag('X?', line.query_name, 'Z')
            temp_bam.write(temp_bam_entry)

        read_count += 1
        if read_count % 200000 == 0:
            print(f"{int(read_count/2)} read pairs processed...")

    in_bam_file.close()
    temp_bam.close()

    if o.numAlignReads != 0:
        fAlign1.close()
        fAlign2.close()

    print("Sorting reads on tag sequence...")

    pysam.sort("-n", "-@", f"{o.cores}", "-o", f"{o.prefix}.temp.sort.bam", f"{o.prefix}.temp.bam")
    # Sort by read name, which will be the tag sequence in this case.
    os.remove(f"{o.prefix}.temp.bam")

    #Extracting tags and sorting based on tag sequence is complete. 
    #This block of code now performs the consensus calling on the tag 
    #families in the temporary name sorted bam file.
    
    seq_dict = {}
    qual_dict = {}
    fam_size_x_axis = []
    fam_size_y_axis = []

    read1_cDNAcs_len = 0
    read2_cDNAcs_len = 0
    in_bam_file = pysam.AlignmentFile(
        f"{o.prefix}.temp.sort.bam", "rb", check_sq=False
        )
    first_line = next(in_bam_file)
    readsCtr += 1
    FinalValue = pysam.AlignedSegment()
    FinalValue.query_name = "FinalValue-0#0"

    seq_dict[first_line.query_name.split('#')[0].split('-')[1]] = {':1': [], ':2': []}
    qual_dict[first_line.query_name.split('#')[0].split('-')[1]] = {':1': [], ':2': []}
    seq_dict[first_line.query_name.split('#')[0].split('-')[1]][first_line.query_name.split('#')[1]].append(first_line.query_sequence)
    qual_dict[first_line.query_name.split('#')[0].split('-')[1]][first_line.query_name.split('#')[1]].append(list(first_line.query_qualities))
    tag_count_dict1 = defaultdict(lambda: 0)
    tag_count_dict2 = defaultdict(lambda: 0)

    print("Creating consensus reads...")

    for line in iteratorWrapper(in_bam_file.fetch(until_eof=True), FinalValue): 
        tag, subtag_order = first_line.query_name.split('#')[0].split('-')[0], first_line.query_name.split('#')[1]
        
        if line.query_name.split('#')[0].split('-')[0] == tag:
            readsCtr += 1
            if line.query_name.split('#')[0].split('-')[1] in seq_dict:
                seq_dict[line.query_name.split('#')[0].split('-')[1]][line.query_name.split('#')[1]].append(line.query_sequence)
                qual_dict[line.query_name.split('#')[0].split('-')[1]][line.query_name.split('#')[1]].append(list(line.query_qualities))
            else:
                seq_dict[line.query_name.split('#')[0].split('-')[1]] = {':1': [], ':2': []}
                qual_dict[line.query_name.split('#')[0].split('-')[1]] = {':1': [], ':2': []}
                seq_dict[line.query_name.split('#')[0].split('-')[1]][line.query_name.split('#')[1]].append(line.query_sequence)
                qual_dict[line.query_name.split('#')[0].split('-')[1]][line.query_name.split('#')[1]].append(list(line.query_qualities))

        else:
            for subtag in seq_dict:
                if len(seq_dict[subtag][':1']) != len(seq_dict[subtag][':2']):
                    raise Exception('ERROR: Read counts for Read1 and Read 2 do not match for tag %s' % tag)
                
                for tag_subtype in seq_dict[subtag]:

                    if len(seq_dict[subtag][tag_subtype]) > 0:
                        tag_count_dict1[len(seq_dict[subtag][tag_subtype])] += 1
                        familyCtr += 1
                    if len(seq_dict[subtag][tag_subtype]) < o.minmem1:
                        seq_dict[subtag][tag_subtype] = []
                        qual_dict[subtag][tag_subtype] = []
                        smallFamilySize += 1
                    elif o.minmem1 <= len(seq_dict[subtag][tag_subtype]) <= o.maxmem1:  # Tag types w/o reads should not be submitted
                        #  as long as minmem is > 0
                        seq_dict[subtag][tag_subtype] = [consensus_caller(seq_dict[subtag][tag_subtype], o.cutoff1, tag, True),
                                                str(len(seq_dict[subtag][tag_subtype]))]
                        qual_dict[subtag][tag_subtype] = qual_calc(qual_dict[subtag][tag_subtype])
                        numPCRcs += 1
                    elif len(seq_dict[subtag][tag_subtype]) > o.maxmem1:
                        seq_dict[subtag][tag_subtype] = [consensus_caller(seq_dict[subtag][tag_subtype][:o.maxmem1], o.cutoff1, tag, True),
                                                str(len(seq_dict[subtag][tag_subtype]))]
                        qual_dict[subtag][tag_subtype] = qual_calc(qual_dict[subtag][tag_subtype])
                        numPCRcs += 1

                if o.write_PCRcs is True:

                    if len(seq_dict[subtag][':1']) != 0 and len(seq_dict[subtag][':2']) != 0:
                        corrected_qual_score = map(lambda x: x if x < 41 else 41, qual_dict[subtag][':1'])
                        read1_PCRcs_fq_file.write('@%s-%s#ab/1\n%s\n+%s\n%s\n' %
                                                (tag, subtag, seq_dict[subtag][':1'][0], seq_dict[subtag][':1'][1], "".join(chr(x + 33)
                                                                                        for x in corrected_qual_score)))

                        corrected_qual_score = map(lambda x: x if x < 41 else 41, qual_dict[subtag][':2'])
                        read2_PCRcs_fq_file.write('@%s-%s#ab/2\n%s\n+%s\n%s\n' %
                                                (tag, subtag, seq_dict[subtag][':2'][0], seq_dict[subtag][':2'][1], "".join(chr(x + 33)
                                                                                        for x in corrected_qual_score)))
                    else:
                        sys.stderr.write("%s or %s not greater than 0\n" % (len(seq_dict[subtag][':1']), len(seq_dict[subtag][':2'])))

            if o.without_cDNAcs is False:
                consensusList = [seq_dict[x][':1'][0] for x in seq_dict if seq_dict[x][':1'] != []]
                if len(consensusList) > 0:
                    numcDNAcs += 1
                    tag_count_dict2[len(consensusList)] += 1
                if len(consensusList) < o.minmem2:
                    cDNAcs_read_1=[]
                    read1_cDNAcs_len = len(cDNAcs_read_1)
                    smallFamilySize += 1
                elif o.minmem2 <= len(consensusList) <= o.maxmem2:
                    cDNAcs_read_1 = [consensus_caller(consensusList, o.cutoff2, tag, False),
                                str(sum([int(seq_dict[x][':1'][1]) for x in seq_dict if seq_dict[x][':1'] != []]))]
                    cDNAcs_read_1_qual = map(lambda x: x if x < 41 else 41, qual_calc([qual_dict[x][':1'] for x in qual_dict if qual_dict[x][':1'] != []]))
                    read1_cDNAcs_len = len(cDNAcs_read_1)
                    Z = cDNAcs_read_1[0]
                    lenZ = len(Z)
#                    print (len(z))
#                    print (z.count('N'))
###New lines below. If N count exceeds cutoff, will change every cs nt to an 'N' and update qual to '!'.
                    if Z.count('N')/float(lenZ) > o.Ncutoff:
                        highN_cDNAcs +=1
                        cDNAcs_read_1 = ['N' * lenZ,
                                str(sum([int(seq_dict[x][':1'][1]) for x in seq_dict if seq_dict[x][':1'] != []]))] 
                        cDNAcs_read_1_qual = map(lambda x: x if x < 41 else 41, qual_calc([qual_dict[x][':1'] for x in qual_dict if qual_dict[x][':1'] != []]))
#By changing the read qual scores to reads that fail the Ncutoff filter, there are alignment problems downstream (incomplete read pairs). To preserve
#	function of the Ncutoff filter I've opted to not edit qual score for failing reads, although technically the updated score should be 0/!
#                        cDNAcs_read_1_qual = map(lambda x: 0 if x > 0 else 0, qual_calc([qual_dict[x][':1'] for x in qual_dict if qual_dict[x][':1'] != []]))
###Original lines below did not count the read length properly. Updated to allow for functional Ncutoff for cDNAcs reads.
#                    if cDNAcs_read_1.count('N')/float(read1_cDNAcs_len) > o.Ncutoff:
#                        highN_cDNAcs += 1
#                        cDNAcs_read_1 = 'N' * read1_cDNAcs_len
#                        cDNAcs_read_1_qual = '!' * read1_cDNAcs_len
                elif len(consensusList) > o.maxmem2:
                    cDNAcs_read_1 = [consensus_caller(consensusList[:o.maxmem2], o.cutoff2, tag, False),
                                str(sum([int(seq_dict[x][':1'][1]) for x in seq_dict if seq_dict[x][':1'] != []]))]
                    cDNAcs_read_1_qual = map(lambda x: x if x < 41 else 41, qual_calc([qual_dict[x][':1'] for x in qual_dict if qual_dict[x][':1'] != []]))
                    read1_cDNAcs_len = len(cDNAcs_read_1)
                    Q = cDNAcs_read_1[0]
                    lenQ = len(Q)
###New lines again.
                    if Q.count('N')/float(lenQ) > o.Ncutoff:
                        highN_cDNAcs +=1
                        cDNAcs_read_1 = ['N' * lenQ,
                                str(sum([int(seq_dict[x][':1'][1]) for x in seq_dict if seq_dict[x][':1'] != []]))]
                        cDNAcs_read_1_qual = map(lambda x: x if x < 41 else 41, qual_calc([qual_dict[x][':1'] for x in qual_dict if qual_dict[x][':1'] != []]))
#                        cDNAcs_read_1_qual = map(lambda x: 0 if x > 0 else 0, qual_calc([qual_dict[x][':1'] for x in qual_dict if qual_dict[x][':1'] != []]))

###Original lines again.
#                    if cDNAcs_read_1.count('N')/float(read1_cDNAcs_len) > o.Ncutoff:
#                        highN_cDNAcs += 1
#                        cDNAcs_read_1 = 'N' * read1_cDNAcs_len
#                        cDNAcs_read_1_qual = '!' * read1_cDNAcs_len
                
                consensusList = [seq_dict[x][':2'][0] for x in seq_dict if seq_dict[x][':2'] != []]
                if len(consensusList) > 0:
                    numcDNAcs += 1
                    tag_count_dict2[len(consensusList)] += 1
                if len(consensusList) < o.minmem2:
                    cDNAcs_read_2=[]
                    read2_cDNAcs_len = len(cDNAcs_read_2)
                elif o.minmem2 <= len(consensusList) <= o.maxmem2:
                    cDNAcs_read_2 = [consensus_caller(consensusList, o.cutoff2, tag, False),
                                str(sum([int(seq_dict[x][':2'][1]) for x in seq_dict if seq_dict[x][':2'] != []]))]
                    cDNAcs_read_2_qual = map(lambda x: x if x < 41 else 41, qual_calc([qual_dict[x][':2'] for x in qual_dict if qual_dict[x][':2'] != []]))
                    read2_cDNAcs_len = len(cDNAcs_read_2)
                    P = cDNAcs_read_2[0]
                    lenP = len(P)
###New lines again.
                    if P.count('N')/float(lenP) > o.Ncutoff:
                        highN_cDNAcs +=1
                        cDNAcs_read_2 = ['N' * lenP,
                                str(sum([int(seq_dict[x][':2'][1]) for x in seq_dict if seq_dict[x][':2'] != []]))]
                        cDNAcs_read_1_qual = map(lambda x: x if x < 41 else 41, qual_calc([qual_dict[x][':1'] for x in qual_dict if qual_dict[x][':1'] != []]))
#                        cDNAcs_read_2_qual = map(lambda x: 0 if x > 0 else 0, qual_calc([qual_dict[x][':2'] for x in qual_dict if qual_dict[x][':2'] != []]))
###Original lines again.
#                    if cDNAcs_read_2.count('N')/float(read2_cDNAcs_len) > o.Ncutoff:
#                        highN_cDNAcs += 1
#                        cDNAcs_read_2 = 'N' * read2_cDNAcs_len
#                        cDNAcs_read_2_qual = '!' * read2_cDNAcs_len
                elif len(consensusList) > o.maxmem2:
                    cDNAcs_read_2 = [consensus_caller(consensusList[:o.maxmem2], o.cutoff2, tag, False),
                                str(sum([int(seq_dict[x][':2'][1]) for x in seq_dict if seq_dict[x][':2'] != []]))]
                    cDNAcs_read_2_qual = map(lambda x: x if x < 41 else 41, qual_calc([qual_dict[x][':2'] for x in qual_dict if qual_dict[x][':2'] != []]))
                    read2_cDNAcs_len = len(cDNAcs_read_2)
                    S = cDNAcs_read_2[0]
                    lenS = len(S)
###New lines again.
                    if S.count('N')/float(lenS) > o.Ncutoff:
                        highN_cDNAcs +=1
                        cDNAcs_read_2 = ['N' * lenS,
                                str(sum([int(seq_dict[x][':2'][1]) for x in seq_dict if seq_dict[x][':2'] != []]))]
                        cDNAcs_read_1_qual = map(lambda x: x if x < 41 else 41, qual_calc([qual_dict[x][':1'] for x in qual_dict if qual_dict[x][':1'] != []]))
#                        cDNAcs_read_1_qual = map(lambda x: 0 if x > 0 else 0, qual_calc([qual_dict[x][':2'] for x in qual_dict if qual_dict[x][':2'] != []]))
###Original lines again.
#                    if cDNAcs_read_2.count('N')/float(read2_cDNAcs_len) > o.Ncutoff:
#                        highN_cDNAcs += 1
#                        cDNAcs_read_2 = 'N' * read2_cDNAcs_len
#                        cDNAcs_read_2_qual = '!' * read2_cDNAcs_len
                if tag.count('N') >0:
                    nUMIs += 1
                if ('A' * o.rep_filt in tag or 'C' * o.rep_filt in tag or 'G' * o.rep_filt in tag or 'T' * o.rep_filt in tag):
                    monoNtUMIs += 1
                if read1_cDNAcs_len != read2_cDNAcs_len:
                    failedcDNAcs += 1
                if read1_cDNAcs_len != 0 and read2_cDNAcs_len != 0 and tag.count('N') == 0 and \
                                        'A' * o.rep_filt not in tag and 'C' * o.rep_filt not in tag and \
                                        'G' * o.rep_filt not in tag and 'T' * o.rep_filt not in tag:
                    num_filt_cDNAcs += 2
                    read1_cDNAcs_fq_file.write('@%s/1\n%s\n+%s\n%s\n' % (tag, cDNAcs_read_1[0], cDNAcs_read_1[1],
                                                                        "".join(chr(x + 33) for x in cDNAcs_read_1_qual)))
                    read2_cDNAcs_fq_file.write('@%s/2\n%s\n+%s\n%s\n' % (tag, cDNAcs_read_2[0], cDNAcs_read_2[1],
                                                                        "".join(chr(x + 33) for x in cDNAcs_read_2_qual)))

            # reset conditions for next tag family
            if line != FinalValue:
                readsCtr += 1
                first_line = line
                seq_dict = {}
                qual_dict = {}
                read1_cDNAcs_len = 0
                read2_cDNAcs_len = 0
                cDNAcs_read_1 = ''
                cDNAcs_read_2 = ''
                Z = ''
                lenZ = 0
                Q = ''
                lenQ = 0
                P = ''
                lenP = 0
                S = ''
                lenS = 0
                
                seq_dict[line.query_name.split('#')[0].split('-')[1]] = {':1': [], ':2': []}
                qual_dict[line.query_name.split('#')[0].split('-')[1]] = {':1': [], ':2': []}
                seq_dict[line.query_name.split('#')[0].split('-')[1]][line.query_name.split('#')[1]].append(line.query_sequence)
                qual_dict[line.query_name.split('#')[0].split('-')[1]][line.query_name.split('#')[1]].append(list(line.query_qualities))

    if o.write_PCRcs is True:
        read1_PCRcs_fq_file.close()
        read2_PCRcs_fq_file.close()

    if o.without_cDNAcs is False:
        read1_cDNAcs_fq_file.close()
        read2_cDNAcs_fq_file.close()
    
# Try to plot the tag family sizes
    if o.tagstats is True:
        tag_stats_file1 = open(o.prefix + ".tagstats1.txt", 'w')

        x_value1 = []
        y_value1 = []
        total_reads = sum([tag_count_dict1[tag_family_size] * tag_family_size for tag_family_size
                            in tag_count_dict1])

        for tag_family_size in sorted(tag_count_dict1):
            fraction = (tag_count_dict1[tag_family_size] * tag_family_size) / float(total_reads)
            tag_stats_file1.write('%d\t%d\t%f\n' % (tag_family_size, tag_count_dict1[tag_family_size], fraction))
            x_value1.append(tag_family_size)
            y_value1.append(fraction)
        
        tag_stats_file2 = open(o.prefix + ".tagstats2.txt", 'w')

        x_value2 = []
        y_value2 = []
        total_reads = sum([tag_count_dict2[tag_family_size] * tag_family_size for tag_family_size
                            in tag_count_dict2])
        
        for tag_family_size in sorted(tag_count_dict2):
            fraction = (tag_count_dict2[tag_family_size] * tag_family_size) / float(total_reads)
            tag_stats_file2.write('%d\t%d\t%f\n' % (tag_family_size, tag_count_dict2[tag_family_size], fraction))
            x_value2.append(tag_family_size)
            y_value2.append(fraction)

        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt

            plt.figure(1)
            plt.bar(x_value1, y_value1)
            plt.xlabel('Family Size')
            plt.ylabel('Proportion of Total Reads')
            plt.savefig(o.prefix + 'family_size1.png', bbox_inches='tight')
            
            plt.figure(2)
            plt.bar(x_value2, y_value2)
            plt.xlabel('Family Size')
            plt.ylabel('Proportion of Total Reads')
            plt.savefig(o.prefix + 'family_size2.png', bbox_inches='tight')
			

        except ImportError:
            sys.stderr.write(
                'matplotlib not present. Only tagstats file will be generated.'
                )

        tag_stats_file1.close()
        tag_stats_file2.close()
		
    endTime = datetime.datetime.now()
    # Print consensus making statistics
    startTimeStr = startTime.strftime("%A, %d. %B %Y %I:%M%p")
    endTimeStr = endTime.strftime("%A, %d. %B %Y %I:%M%p")
    cmStatsFile = open(f"{o.prefix}_cmStats.txt", 'w')
    cmStatsFile.write(
        f"Consensus Making Statistics:\n"
        f"Command: {' '.join(sys.argv)}\n"
        f"Started at {startTimeStr}\n"
        f"Finished at {endTimeStr}\n"
        f"{read_count} raw reads UMI processed\n"
        f"{read_pair_count} total read pairs\n"
        f"{readsCtr} reads processed for consensus calling\n"
        f"{familyCtr/2} families processed (consisting of reads 1 & 2 with shared 16UMI-11UMI barcodes)\n"
#        f"\t{zeroFamilySize} unrepresented families\n"
#        f"\t{badUMIs} families (cDNAcs pairs) filtered for UMIs with mononucleotide repeats\n"
        f"{numPCRcs} PCRcs made\n"
        f"\t{numPCRcs/2} unique 11 nt UMIs discovered\n"
#        f"\t{highN_PCRcs} PCRcs filtered for excessive Ns\n"
        f"{numcDNAcs} cDNAcs made\n"
        f"\t{numcDNAcs/2} unique 16 nt UMIs discovered\n"
#        f"\t{failedcDNAcs} cDNAcs failed due to missing PCRcs\n"
        f"\t{smallFamilySize} cDNA families failed with family size < {o.minmem2}, but > 0\n"
        f"\t{nUMIs + monoNtUMIs} cDNAcs failed due to bad UMIs\n"
        f"\t{highN_cDNAcs} cDNAcs filtered for excessive Ns\n"
        f"{num_filt_cDNAcs} final cDNAcs reads passing filters\n"
        )
    cmStatsFile.close()
    sys.stderr.write(
        f"Consensus Making Statistics:\n"
        f"Command: {' '.join(sys.argv)}\n"
        f"Started at {startTimeStr}\n"
        f"Finished at {endTimeStr}\n"
        f"{read_count} raw reads UMI processed\n"
        f"{read_pair_count} total read pairs\n"
        f"{readsCtr} reads sorted and processed for consensus calling\n"
        f"{familyCtr/2} families processed (consisting of reads 1 & 2 with shared 16UMI-11UMI barcodes)\n"
#        f"\t{zeroFamilySize} unrepresented families\n"
#        f"\t{badUMIs} families (cDNAcs pairs) filtered for UMIs with mononucleotide repeats\n"
        f"{numPCRcs} PCRcs made\n"
        f"\t{numPCRcs/2} unique 11 nt UMIs discovered\n"
#        f"\t{highN_PCRcs} PCRcs filtered for excessive Ns\n"
        f"{numcDNAcs} cDNAcs made\n"
        f"\t{numcDNAcs/2} unique 16 nt UMIs discovered\n"
#        f"\t{failedcDNAcs} cDNAcs failed due to missing PCRcs\n"
        f"\t{smallFamilySize} cDNA families failed with family size < {o.minmem2}, but > 0\n"
        f"\t{nUMIs + monoNtUMIs} cDNAcs failed due to bad UMIs\n"
        f"\t{highN_cDNAcs} cDNAcs filtered for excessive Ns\n"
        f"{num_filt_cDNAcs} final cDNAcs reads passing filters\n"
        )

if __name__ == "__main__":
    main()
