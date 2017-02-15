#!/usr/bin/env python
__author__ = 'Shlomo Shenzis'
from Bio import Seq
from Bio.Alphabet import generic_rna
import re
import os
import argparse
import bisect
import subprocess
import copy
from itertools import islice
from collections import Counter


#outputBedToolsFilename = "getFastaOut_randCon.fa"
#annotationsFilename = "exons_rand_ctrl_EXACT.bed"
POSITIVE = '+'
NEGATIVE = '-'
# forcing the middle to be divided by 3 and possible of endless
POS_REG = "(?=((ATG)([A-Z][A-Z][A-Z])*?(TGA|TAA|TAG|$)))"
INF_REG = "((M)([^*])*?($))"

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'


# ------------------------------------
# --------FindExons-------------------
# ------------------------------------
def find(args):
    output = open('exons', 'w')
    circsH = open(args.circ, 'r')
    anotH = open(args.annotation, 'r')
    anot = anotH.read().strip().split('\n')
    out_circ = open('_CIRC_TEMP_.bed', 'w')
    filtered_out = open('filtered_out_circles.bed', 'w')
    filtered_out_ex = open('filtered_out_circles_EXACT.bed', 'w')
    # counter for filtered out circles
    numFiltered = 0
    # counter for filtered out circles in the EXACT case.
    numFiltered_ex = 0

    # Create temporary modified circles .bed file:
    circles = list()
    # We want each row in the new file have the <chr>\t<locA>\t<locB>\t<ID>
    # so we create a counter.
    counter = 0
    first = True
    for circle in circsH:
        # First row is the header, throw it away.
        if first:
            first = False
            continue
        else:
            circle = circle.strip('\n').strip('\t').split('\t')
            # Filter out circles with zero reads.
            circles.append(circle[0:6])
            # write the rows!
            out_circ.write(circle[0]+'\t'+circle[1]+'\t'+circle[2]+'\t'+str(counter)+'\n')
            counter += 1
    out_circ.close()

    # Intersect circ and annotation.
    # In the intersection file each circle has the annotation lines containing it (-f 1.0).
    os.system('bedtools intersect -a _CIRC_TEMP_.bed -b '+args.annotation+' -wa -wb -f 1.0 > _INTERSECTION_.bed')
    # The following row of code is garbage deletion but leaving it is good for debugging.
    #os.remove('_CIRC_TEMP_.bed')

    # Open the intersection file we created.
    inerxH = open('_INTERSECTION_.bed', 'r')
    inerx = inerxH.read().strip('\n').split('\n')
    inerxH.close()
    # No circles scanned yet (@ = 'none'):
    current_circle = '@'
    winner_exons = list()
    winner_introns = list()
    winner_score = -1
    winner_name = ''
    for row in inerx:
        row = row.strip('\n').strip('\t').split('\t')
        # We scan the rows, each time when the circle changes (ID)
        # We make the new circle the current circle and update
        # the winning transcript as the transcript of the prev circle exons.
        if current_circle != row[3]:
            # Just a small check that this is not the first circle
            # when we dont yet have any previous circles.
            if current_circle != '@':
                # add data to circles:
                # calculate blockStarts (offsets, one of the columns in the BED specification)
                blockStarts = list()
                blockStarts.append(0)
                for i in range(1, len(winner_exons)):
                    blockStarts.append(blockStarts[i-1]+int(winner_exons[i-1])+int(winner_intorns[i-1]))
                # Save the number of exons in the winner transcript.
                ex_len = len(winner_exons)
                # Convert all data to comma separated strings.
                winner_exons = ",".join(winner_exons)
                blockStarts = ",".join(map(str, blockStarts))
                # Add all data to the circles(circle per row) array
                # It is - [0,           1,                 2,           3,...] -> circles list of lists.
                #          |             |                  |           |,...
                #          v             v                  v           v,...
                #  [circ. columns]  [circ. columns]  [circ. columns] ...
                circles[int(current_circle)].append(str(ex_len))
                circles[int(current_circle)].append(winner_exons)
                circles[int(current_circle)].append(blockStarts)
                circles[int(current_circle)].append(str(winner_score))
                circles[int(current_circle)].append(winner_name)
                #Set up next circle:
                winner_exons = list()
                winner_introns = list()
                winner_score = -1
                winner_name = ''
            current_circle = row[3]
            locA = row[1]
            locB = row[2]

        #scan exons and introns
        exons = list()
        introns = list()
        score = 0
        name = row[8]
        starts = row[12].strip(',').split(',')
        ends = row[13].strip(',').split(',')
        # Find the exon in the current transcript that starts
        # as close as possible from the "left"(lower) to the circle
        # start location.
        placeStart = bisect.bisect_left(map(int,starts), int(locA))
        # StartLoop -  keeps track if the exons scanning loop
        # should start at all, as sometimes we finish the whole
        # circle exons scan in the special case below.
        startLoop = True
        # The case where the circle starts after all the exons in this transcript
        # may be either IN the last exon (Internal) or not at all.
        if placeStart == len(starts):
            placeStart -= 1
        # The awesome bull's eye case (starts exactly at exon)
        # In this case the score is at least 1.
        if starts[placeStart] == locA:
            score += 1
            low = placeStart
        else: # does not start at exon
            low = placeStart-1
            if int(ends[low]) < int(locA):  # we fell in intron-before
                #print "fell in intron before..."
                if placeStart+1 < len(starts):
                    # Next exon is out of our scope (finished here)
                    if int(starts[placeStart+1]) >= int(locB):
                        # maybe it ENDS at end of exon though...
                        if int(ends[placeStart]) == int(locB):
                            score += 1
                        exons.append(str(int(locB) - int(locA)))
                        introns.append('0')
                        # no need to loop over other exons, already done.
                        startLoop = False
                if startLoop:
                    # the case where the end is at the intron after the first exon.
                    if int(ends[placeStart]) >= int(locB):
                        # or might be just a bull's eye at the end.
                        if int(ends[placeStart]) == int(locB):
                            score += 1
                        exons.append(str(int(locB) - int(locA)))
                        introns.append('0')
                        startLoop = False
                    else:
                        # intron and exon as first circ exon.
                        exons.append(str(int(ends[placeStart]) - int(locA)))
                        introns.append(str(int(starts[placeStart+1]) - int(ends[placeStart])))
                        low = placeStart + 1
            else:  # circle start in middle of exon
                if low+1 < len(starts):
                    if int(starts[low+1]) >= int(locB):
                        if int(ends[low]) == int(locB):
                            score += 1
                        exons.append(str(int(locB) - int(locA)))
                        introns.append('0')
                        startLoop = False
                if startLoop:
                    if int(ends[low]) >= int(locB):
                        if int(ends[low]) == int(locB):
                            score += 1
                        #and the end is in the same exon!
                        exons.append(str(int(locB) - int(locA)))
                        introns.append('0')
                        startLoop = False
                    else:
                        # Sanity check:
                        if (int(ends[low]) - int(locA)) < 0:
                            print "Negative!!! (start mid ex) at " + locA +":"+locB
                        exons.append(str(int(ends[low]) - int(locA)))
                        introns.append(str(int(starts[placeStart]) - int(ends[low])))
                        low = placeStart
        # Now the actual loop. After all the start special cases,
        # we scan exon by exon until we finish the circle.
        # However, we also have to check here the END special cases.
        while startLoop:
            # No more exons left
            if low >= len(starts):
                break
            # Just some array out-of-bound problems fix.
            if low + 1 < len(starts):
                nextStart = starts[low+1]
            else:
                nextStart = ends[low]
            # ended in middle of exon
            if int(locB) < int(ends[low]):
                # Sanity check:
                if (int(locB) - int(starts[low])) < 0:
                    print "Negative!!! (end mid ex) at " + locA +":"+locB
                exons.append(str(int(locB) - int(starts[low])))
                introns.append('0')
                break
            # last exon
            if int(locB) < int(nextStart) or int(locB) == int(ends[low]):
                # Maybe bull's eye at end?
                if int(locB) == int(ends[low]):
                    score += 1
                # Another sanity check:
                if (int(locB) - int(starts[low])) < 0:
                    print "Negative!!! (last ex) at " + locA +":"+locB
                # exon and intron as last exon
                exons.append(str(int(locB) - int(starts[low])))
                introns.append('0')
                break
            #Just another exon:
            # sanity check:
            if (int(ends[low]) - int(starts[low])) < 0:
                    print "Negative!!! (reg ex) at " + locA +":"+locB
            exons.append(str(int(ends[low]) - int(starts[low])))
            introns.append(str(int(nextStart) - int(ends[low])))
            low += 1
        # Keep track of the winner transcript, if this one is better
        # make it the winner.
        if score > winner_score:
            winner_score = score
            winner_name = name
            winner_exons = copy.deepcopy(exons)
            winner_intorns = copy.deepcopy(introns)

    # Now when all the calculations are done, we can write output
    # and find statistics:
    numOK = 0
    for circle in circles:
        # If a circle has less than 7 columns, it means it should be filtered
        # As no transcripts where found for it:
        if len(circle) > 7:
            # Make the EXACT filtration.
            if circle[9] != '2':
                filtered_out_ex.write(circle[0]+'\t'+circle[1]+'\t'+circle[2]+'\t'+circle[3]+'\t'+circle[4]+'\n')
                numFiltered_ex += 1
            # reorganize all circle data into BED format columns:
            # (4 and 10 are merged - those are the CG name and gene-name)
            circle_line = circle[0]+'\t'+circle[1]+'\t'+circle[2]+'\t'+circle[4]+':'+circle[10]+'\t'+circle[5]+'\t'+\
                circle[3]+'\t'+circle[1]+'\t'+circle[2]+'\t0,0,255\t'+circle[6]+'\t'+circle[7]+'\t'+\
                circle[8]+'\t'+circle[9]
            output.write(circle_line + '\n')
            numOK += 1
        else: # filter it
            circle_line = '\t'.join(circle)
            filtered_out.write(circle_line + '\n')
            filtered_out_ex.write(circle_line + '\n')
            numFiltered += 1
            numFiltered_ex += 1
    filtered_out.close()
    filtered_out_ex.close()
    output.close()
    # Write all possible outputs one may need:
    print "Creating EXACT (only score 2) output..."
    os.system("cat exons | awk '$13 == 2 {print $0}' | cut -f 1-12 > exons_EXACT.bed")
    print "Creating BED12 file..."
    os.system("cat exons | cut -f 1-12 > exons_bed12.bed")
    print "Creating BED6 file..."
    os.system("cat exons | cut -f 1-12 | bed12ToBed6 > exons_bed6.bed")
    print "Creating modified circles file..."
    os.system("grep -F -v -f filtered_out_circles.bed "+args.circ+" > filtered_"+args.circ)
    print "Creating modified EXACT circles file..."
    os.system("grep -F -v -f filtered_out_circles_EXACT.bed "+args.circ+" > filtered_EXACT_"+args.circ)
    print "Finished."
    print str(numOK) + " circles where annotated and splitted to exons."
    print "Of them, "+str(numFiltered_ex-numFiltered) + " where INTERNAL (non-EXACT)."
    print "Output can be found at -> exons.bed"
    print str(numFiltered) + " circles where filtered in the process, they can be found at -> filtered_out_circles.bed"
    print str(numFiltered_ex) + " circles where filtered in the process for EXACT, they can be found at -> filtered_out_circles_EXACT.bed"


def get_fasta(args):
    os.system('bedtools getfasta -fi ' + args.system + ' -bed '+args.exons+' -s -split -name -fo getFastaOutput.fa')
# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------

def add_introns(circ_start, circ_end, orf_start_ref, orf_end_ref, STRAND_SIGN, gene_name, circ_name, chr):
    #print "Counting introns for: ", gene_name, circ_name
    lines = anot_dict[gene_name]
    line = lines[-1]
    for elem in lines:
        if str(circ_start) in elem.split()[8].strip(',').split(','):
            if str(circ_end) in elem.split()[9].strip(',').split(','):
                line = elem
                break
    exon_count = int(line.split()[7])
    exons_starts = map(int, line.split()[8].strip(',').split(','))  # [int(x) for x in line.split()[9].strip(',').split(',')]
    exons_ends = map(int, line.split()[9].strip(',').split(','))  # [int(x) for x in line.split()[10].strip(',').split(',')]
    # sanity assertions
    if STRAND_SIGN == POSITIVE and not int(circ_start) < int(circ_end): (circ_start, circ_end) = (circ_end, circ_start)
    if STRAND_SIGN == NEGATIVE and not int(circ_end) < int(circ_start): (circ_start, circ_end) = (circ_end, circ_start)
    sign = 1
    if STRAND_SIGN == NEGATIVE:
        sign = -1
    #if STRAND_SIGN == POSITIVE and not int(orf_start_ref) > int(orf_end_ref): (orf_start_ref, orf_end_ref) = (orf_end_ref, orf_start_ref)

    # compute start and end position for the ORF.
    orf_start_pos = int(circ_start)
    orf_end_pos = int(circ_start)

    counter_start = orf_start_ref
    counter_end = orf_end_ref

    for i in (range(bisect.bisect_left(map(int,exons_starts), int(circ_start)), len(exons_starts)) if STRAND_SIGN == POSITIVE else
              range(bisect.bisect_left(map(int,exons_ends), int(circ_start)), -1, -1)):
        if counter_start <= 0:
            break
        len_exon = int(exons_ends[i])-int(exons_starts[i])
        if len_exon > counter_start:
            orf_start_pos += sign*counter_start
            break
        else:
            counter_start -= len_exon
            orf_start_pos += sign*len_exon
            if STRAND_SIGN == POSITIVE and i+1 < len(exons_starts):
                orf_start_pos += sign*(int(exons_starts[i+1])-int(exons_ends[i]))
            elif STRAND_SIGN == NEGATIVE and i > 0:
                orf_start_pos += sign*(int(exons_starts[i])-int(exons_ends[i-1]))
            else:
                break

    for i in (range(bisect.bisect_left(map(int,exons_starts), int(circ_start)), len(exons_starts)) if STRAND_SIGN == POSITIVE else
              range(bisect.bisect_left(map(int,exons_ends), int(circ_start)), -1, -1)):
        if counter_end <= 0:
            break
        len_exon = int(exons_ends[i])-int(exons_starts[i])
        if len_exon > counter_end:
            orf_end_pos += sign*counter_end
            break
        else:
            counter_end -= len_exon
            orf_end_pos += sign*len_exon
            if STRAND_SIGN == POSITIVE and i+1 < len(exons_starts):
                orf_end_pos += sign*(int(exons_starts[i+1])-int(exons_ends[i]))
            elif STRAND_SIGN == NEGATIVE and i > 0:
                orf_end_pos += sign*(int(exons_starts[i])-int(exons_ends[i-1]))
            else:
                break

    return orf_start_pos, orf_end_pos


# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------

# Perform checking over each of the ORFs that were found. It should hold that:
# 1. The ORF should start before the junction.
# 2. The ORF should end after the junction
# 3.
def perform_checks(pro, circ_start, circ_end, circ_data):
    len_data = len(circ_data)

    if len(pro.group(1)) % 3 != 0:  # The ORF is devisable by 3 as full amino acids.
        #print "Error"
        return False
    if len(pro.group(1)) < 63:  # 21 Amino acids. as the stop is truncated
        return False
    if (len_data - pro.start(1) <= 0):  # ORF starts at first repetition of sequence
        #print bcolors.WARNING + "1'st condition" + bcolors.ENDC
        #print bcolors.WARNING + "len_data", str(len_data) + bcolors.ENDC
        return False
    lenFromStart = pro.start(1) + len(pro.group(1))
    if (lenFromStart - len_data <= 3):  # The junction is not part of the ORF
        #print bcolors.WARNING + "2'nd condition" + bcolors.ENDC
        return False
    #if (pro.start(1) <= (lenFromStart - len_data)):  # stop not passing start - not needed apparently...
    #    #print bcolors.WARNING + "3'rd condition" + bcolors.ENDC
    #    return False


    #print bcolors.OKGREEN + "An ORF were found!" + bcolors.ENDC
    return True


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# This function find ORF, and write them into outputFile.
# It computes the value of the circle start and end positions in both the nucleotide position and amino-acid position
def findORF(circle, circleData, sign, circ_start, circ_end, chr, firstOut):
    #if sign == '-':
    #    circleData = circleData[::-1]
    doubled_sequence = circleData + circleData + circleData + circleData
    #if sign == POSITIVE:
    p = re.compile(POS_REG, re.I)
    #num =0
    find = False
    bigger = dict()
    bigger[0] = (0,0)
    bigger[1] = (0,0)
    bigger[2] = (0,0)
    for pro in p.finditer(doubled_sequence):
        #if circle.split(":")[0] == 'mbl' and circle.split(":")[1].split('__')[0]=='CG33197-RC':
        #        print pro.group(1)
        if perform_checks(pro, circ_start, circ_end, circleData):
            #if circle.split(":")[0] == 'mbl' and circle.split(":")[1].split('__')[0]=='CG33197-RC':
            #    print 'OK',pro.group(1)
            #    print 'OK',int(pro.start(1))%3
            find = True
            tempStart = int(pro.start(1))
            frame = int(tempStart % 3)
            if bigger[frame] == (0,0):
                bigger[frame] = (pro.group(1),pro.start(1))
                #num = num+1
            else:
                if len(pro.group(1)) > len(bigger[frame][0]):
                    bigger[frame] = (pro.group(1),pro.start(1))
    # inf:
    pinf = re.compile(INF_REG, re.I)
    for pro in pinf.finditer(Seq.translate(doubled_sequence)):
        if len(pro.group())>21 and 1<pro.start()<len(circleData)/3:
            find = True
            if bigger[0] == (0,0):
                bigger[0] = (doubled_sequence[pro.start()*3:],pro.start()*3)
            else:
                if len(pro.group()) > len(bigger[0][0]):
                    bigger[0] = (doubled_sequence[pro.start()*3:],pro.start()*3)
    for pro in pinf.finditer(Seq.translate(doubled_sequence[1:])):
        if len(pro.group())>21 and 1<pro.start()<len(circleData)/3:
            find = True
            if bigger[1] == (0,0):
                bigger[1] = (doubled_sequence[pro.start()*3+1:],pro.start()*3+1)
            else:
                if len(pro.group()) > len(bigger[1][0]):
                    bigger[1] = (doubled_sequence[pro.start()*3+1:],pro.start()*3+1)
    for pro in pinf.finditer(Seq.translate(doubled_sequence[2:])):
        if len(pro.group())>21 and 1<pro.start()<len(circleData)/3:
            find = True
            if bigger[2] == (0,0):
                bigger[2] = (doubled_sequence[pro.start()*3+2:],pro.start()*3+2)
            else:
                if len(pro.group()) > len(bigger[2][0]):
                    bigger[2] = (doubled_sequence[pro.start()*3+2:],pro.start()*3+2)
    if find:
        if args.collapse_frames:
            tuple_largest = max(bigger.values(), key=lambda r: len(r[0]) if r != (0, 0) else -1)
            bigger = {0: tuple_largest, 1: (0, 0), 2: (0, 0)}
        for frame in bigger.keys():
            if not bigger[frame] == (0, 0):
                protein = Seq.translate(bigger[frame][0])
                #print protein
                orf_start_nuc = bigger[frame][1]
                lenFromStart = (bigger[frame][1]-1 + len(bigger[frame][0])-1)
                orf_end_nuc = lenFromStart - len(circleData)*(lenFromStart/len(circleData)) - 1  # - (len(circleData)+1)
                if orf_end_nuc > len(circleData):
                    print 'ERROR'
                    exit(-1)
                lenORF = len(bigger[frame][0])
                #if circle.split(":")[0] == 'esn':
                #    print orf_start_nuc, lenFromStart, orf_end_nuc, lenORF

                firstOut.write(circle.split(":")[0]+";"+circle.split(":")[1].split('__')[0]+";"+chr+":"+circ_start+"-"+circ_end+";"+sign+";")
                    #remove non-translation matches:

                # add the introns so far to the start and end positions we found
                circ_name = circle.split(':')[1]
                gene_name = (circle.split(':')[1]).split('__')[0]
                (orf_start_nuc, orf_end_nuc) = add_introns(int(circ_start), int(circ_end), orf_start_nuc, orf_end_nuc, sign, gene_name, circ_name, chr)
                #if circle.split(":")[0] == 'esn':
                #    print orf_start_nuc, lenFromStart, orf_end_nuc, lenORF
                #    exit(0)
                #print "orf_start_nuc", bcolors.OKBLUE + str(orf_start_nuc) + bcolors.ENDC, "orf_end_nuc", bcolors.OKBLUE + str(orf_end_nuc) + bcolors.ENDC
                firstOut.write(str(orf_start_nuc)+":"+str(orf_end_nuc)+';'+protein+';'+str(len(circleData))+'\n')

        #exit()

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


def extract_ORF(args):

    #open files
    try:
        firstOut = open('output_ORFs.txt', 'w')
        output_bed_file = open(args.fasta, 'r')
    except:
        print bcolors.FAIL + "Cannot open file " + outputBedToolsFilename + ". Exiting..." + bcolors.ENDC

    # create genes dictionary
    genes_dict = create_genes_dict(args.exons)
    #------------------------
    circleName = ''
    circleData = ''
    #circleName_counter = 0
    namesDict = dict()
    for line in output_bed_file.readlines():
        if line[0] == '>':
            circleName = line.strip()[1:]
            if circleName.strip() not in namesDict:
                namesDict[circleName.strip()] = 0
            else:
                namesDict[circleName.strip()] += 1
            circleName = circleName + "__" + str(namesDict[circleName.strip()])
        else:
            circleData = line.strip()
           # if genes_dict[circleName][0] == NEGATIVE: continue
            findORF(circleName, circleData, genes_dict[circleName][0], genes_dict[circleName][1], genes_dict[circleName][2], genes_dict[circleName][3],firstOut)
    firstOut.close()


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------


def create_genes_dict(annotations_filename):

    # try to open the file
    try:
        anot_file = open(annotations_filename, "r")
    except:
        print bcolors.FAIL + "Can\'t open the file " + annotations_filename + bcolors.ENDC
        exit(2)

    # create dict
    dict_genes = {}

    # for each line, add the (CDS_start, CDS_end, STRAND_SIGN) tuple to the "name" place in the dictionary
    lastName_counter = 0
    lastName = ''
    namesDict = dict()
    for line in anot_file.readlines():
        line_split = line.split()
        name = line_split[3].strip()
        if name not in namesDict:
            namesDict[name] = 0
        else:
            namesDict[name] += 1
        name = name + "__" + str(namesDict[name])
        sign = line_split[5]
        circ_start = line_split[1]
        circ_end = line_split[2]
        chr = line_split[0]

        dict_genes[name] = (sign, circ_start, circ_end, chr)

    anot_file.close()
    #print dict_genes
    return dict_genes

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

def create_genes_dict2(annotations_filename):

    # try to open the file
    try:
        anot_file = open(annotations_filename, "r")
    except:
        print "Can\'t open the file " + annotations_filename + ". Exiting..."
        exit(2)

    # create dict
    dict_genes = {}

    for line in anot_file.readlines():
        line_split = line.split()

        gene_and_transc_name = line_split[4]

        sign = line_split[3]
        chrom = line_split[0]
        CDS_start = line_split[5]
        CDS_end = line_split[6]

        dict_genes[gene_and_transc_name] = (sign, chrom, CDS_start, CDS_end)

    anot_file.close()
    return dict_genes


def main2(args):
    try:
        f = open("output_ORFs.txt", 'r')
    except:
        print "---ERROR--- Cannot open 'output_revised.txt'. Exiting..."
        exit(2)

    #create output file
    try:
        o = open("output_CDS_info.txt", 'w')
    except:
        print "---ERROR--- Cannot open 'output_CDS_info.txt'. Exiting..."
        exit(2)

    #-----

    gene_dict = create_genes_dict2(args.annotation)
    counter = 0
    counter_ORF_AT_CDS_START = 0
    counter_BB = 0
    counter_BI = 0
    counter_BA = 0
    counter_IB = 0
    counter_II = 0
    counter_IA = 0
    counter_AB = 0
    counter_AI = 0
    counter_AA = 0
    for line in f.readlines():

        # acquire info of the line. Change it according to new format.
        line_s = line.split(';')
        gene_and_transc_name = line_s[1].strip()
        circ_name = line_s[0].strip()
        chrom_name = line_s[2].split(':')[0].strip()
        circ_start = line_s[2].split(':')[1].split('-')[0].strip()
        circ_end = line_s[2].split(':')[1].split('-')[1].strip()
        orf_start = line_s[4].split(':')[0].strip()
        orf_end = int(line_s[4].split(':')[1].strip())
        sign = line_s[3].split()[0]
        protein = line_s[5].strip()
        circleLen = line_s[6].strip()

        (sign1, chrom, CDS_start, CDS_end) = gene_dict[gene_and_transc_name]
        #print (sign, chrom, CDS_start, CDS_end)
        #-----
        counter += 1
        orf_start = int(orf_start) + 1
        o.write(gene_and_transc_name + ' ; '+circ_name + ' ; '+chrom_name+':'+circ_start+'-'+circ_end+'\n')
        if (sign == '-' and int(orf_start) == int(CDS_end)+1) or (sign == '+' and int(orf_start) == int(CDS_start)):
            atstart = 'YES'
            counter_ORF_AT_CDS_START += 1
        else:
            atstart = 'NO'
        o.write(protein+'\n')
        o.write(sign +' IS_AT_START: ' + atstart + '\n')
        o.write(sign +' ORF absolute start position: ' + str(orf_start)+ '\n')
        o.write(sign +' ORF absolute end position: ' + str(orf_end) + '\n')

        orf_start_CDS_relation = ''
        orf_end_CDS_relation = ''
        #start position relation against the CDS
        if (sign == '+' and (int(orf_start) < int(CDS_start))) or (sign == '-' and (int(orf_start) > int(CDS_end)+1)):
            orf_start_CDS_relation = "BEFORE"
        elif int(CDS_start) <= int(orf_start) <= (int(CDS_end) if sign == '+' else int(CDS_end)+1):
            orf_start_CDS_relation = "INSIDE"
        elif (sign == '+' and int(CDS_end) < int(orf_start)) or (sign == '-' and int(CDS_start) > int(orf_start)):
            orf_start_CDS_relation = "AFTER"
        else:
            print "1)--ERROR-- an unknown case!", sign, CDS_start, CDS_end, orf_start
            exit(2)


        #end position relation against the CDS
        if (sign == '+' and (int(orf_end) < int(CDS_start))) or (sign == '-' and (int(orf_end) > int(CDS_end)+1)):
            orf_end_CDS_relation = "BEFORE"
        elif int(CDS_start) <= int(orf_end) <= (int(CDS_end) if sign == '+' else int(CDS_end)+1):
            orf_end_CDS_relation = "INSIDE"
        elif (sign == '+' and int(CDS_end) < int(orf_end)) or (sign == '-' and int(CDS_start) > int(orf_end)):
            orf_end_CDS_relation = "AFTER"
        else:
            print "2)--ERROR-- an unknown case!", sign, CDS_start, CDS_end, orf_end
            exit(2)



        #statistics computation
        if orf_start_CDS_relation == "BEFORE" and orf_end_CDS_relation == "BEFORE":
            counter_BB += 1
        elif orf_start_CDS_relation == "BEFORE" and orf_end_CDS_relation == "INSIDE":
            counter_BI += 1
        elif orf_start_CDS_relation == "BEFORE" and orf_end_CDS_relation == "AFTER":
            counter_BA += 1
        elif orf_start_CDS_relation == "INSIDE" and orf_end_CDS_relation == "BEFORE":
            counter_IB += 1
        elif orf_start_CDS_relation == "INSIDE" and orf_end_CDS_relation == "INSIDE":
            counter_II += 1
        elif orf_start_CDS_relation == "INSIDE" and orf_end_CDS_relation == "AFTER":
            counter_IA += 1
        elif orf_start_CDS_relation == "AFTER" and orf_end_CDS_relation == "BEFORE":
            counter_AB += 1
        elif orf_start_CDS_relation == "AFTER" and orf_end_CDS_relation == "INSIDE":
            counter_AI += 1
        elif orf_start_CDS_relation == "AFTER" and orf_end_CDS_relation == "AFTER":
            counter_AA += 1
        else:
            print "--ERROR-- an unknown case!"
            exit(2)

        o.write(sign+" ORF starts " + orf_start_CDS_relation + ' the CDS\n')
        o.write(sign+" ORF ends " + orf_end_CDS_relation + " the CDS\n")
        o.write(sign+" Circle exonic length is: " + circleLen + "\n\n")



    print "Number of ORFs starts exactly at the CDS start: ", str(counter_ORF_AT_CDS_START), "out of ", str(counter) +". Percentage:", 100.0*float(counter_ORF_AT_CDS_START)/float(counter)
    print "counter_BB", counter_BB
    print "counter_BI", counter_BI
    print "counter_BA", counter_BA
    print "counter_IB", counter_IB
    print "counter_II", counter_II
    print "counter_IA", counter_IA
    print "counter_AB", counter_AB
    print "counter_AI", counter_AI
    print "counter_AA", counter_AA

def main3(args):
    if not args.output.lower().endswith('.gary'):
        args.output += '.GARY'
    output = open(args.output, 'w')
    orfs = open('output_CDS_info.txt' ,'r').read().strip('\n\n').split('\n\n')
    anot = open(args.annotation, 'r').read().strip('\n').split('\n')
    anot_dict = dict()
    yes=0
    no=0
    for row in anot:
        row = row.strip('\t').split('\t')
        anot_dict[row[4]] = row
    counter = -1
    for orf_line in orfs:
        counter += 1
        orf = copy.deepcopy(orf_line.strip('\n').split('\n'))
        if len(orf) > 3:
            orf_gene = orf[0].split(';')[0].strip(" ")
            orf_start = int(orf[3].split(':')[1].strip())
            # orf_end = int(orf[4].split(':')[1].strip())
            orf_strand = orf[3][0]
            if "ORF starts INSIDE the CDS" in orf[5] and ("ORF ends BEFORE the CDS" in orf[6] or "ORF ends INSIDE the CDS" in orf[6]):
                if orf[5].startswith('+'):
                    sign = '+'
                    #if orf[1][-1] != '*':
                        #print "booya"
                    exonic_length = 0
                    cds_start = int(anot_dict[orf_gene][5])+1  # One nucleotide forward in UCSC
                    exon_starts = anot_dict[orf_gene][8].strip(',').split(',')
                    exon_ends = anot_dict[orf_gene][9].strip(',').split(',')
                    placeEnd = bisect.bisect_left(map(int,exon_ends), orf_start)
                    placeStart = bisect.bisect_left(map(int,exon_starts), cds_start)
                    stop = False
                    if placeStart == len(exon_starts):
                        placeStart -= 1
                    if int(exon_starts[placeStart]) != cds_start:
                        placeStart -= 1
                        if int(exon_ends[placeStart]) >= cds_start:
                            exonic_length += min(orf_start, int(exon_ends[placeStart])) - cds_start
                            if min(orf_start, int(exon_ends[placeStart])) == orf_start:
                                stop = True
                        placeStart += 1
                    #print exonic_length
                    if not stop:
                        if placeEnd == len(exon_ends):
                            placeEnd -= 1
                        if exon_ends[placeEnd] != orf_start:
                            placeEnd -= 1
                            #print "next to last "+str(int(exon_starts[placeEnd+1]) <= orf_start)
                            #print "next to last ",int(exon_starts[placeEnd+1]), orf_start
                            if int(exon_starts[placeEnd+1]) <= orf_start:
                                exonic_length += orf_start-max(cds_start, int(exon_starts[placeEnd+1]))
                        #print exonic_length
                        for i in range(placeStart, placeEnd+1):
                            exonic_length += int(exon_ends[i])-int(exon_starts[i])
                else:
                    sign = '-'
                    #print "We've got a  minus!"
                    #if orf[1][-1] != '*':
                        #print "booya"
                    exonic_length = 0
                    cds_end = int(anot_dict[orf_gene][6]) - 1  # One nucleotide forward in UCSC
                    exon_starts = anot_dict[orf_gene][8].strip(',').split(',')
                    exon_ends = anot_dict[orf_gene][9].strip(',').split(',')
                    placeEnd = bisect.bisect_left(map(int,exon_ends), cds_end)
                    placeStart = bisect.bisect_left(map(int,exon_starts), orf_start)
                    stop = False
                    if placeStart == len(exon_starts):
                        placeStart -= 1
                    if int(exon_starts[placeStart]) != orf_start:
                        placeStart -= 1
                        if int(exon_ends[placeStart]) >= orf_start:
                            exonic_length += min(cds_end, int(exon_ends[placeStart])) - orf_start - 1
                            if min(cds_end, int(exon_ends[placeStart])) == cds_end:
                                stop = True
                        placeStart += 1
                    #print exonic_length
                    if not stop:
                        if placeEnd == len(exon_ends):
                            placeEnd -= 1
                        if exon_ends[placeEnd] != cds_end:
                            placeEnd -= 1
                            #print "next to last "+str(int(exon_starts[placeEnd+1]) <= cds_end)
                            if int(exon_starts[placeEnd+1]) <= cds_end:
                                exonic_length += cds_end-max(orf_start+1, int(exon_starts[placeEnd+1]))
                        #print exonic_length
                        for i in range(placeStart, placeEnd+1):
                            exonic_length += int(exon_ends[i])-int(exon_starts[i])
                # Check what strand it is for direction:
                exonic_length = max(0, exonic_length)
                #print exonic_length
                if exonic_length % 3 == 0:
                    #print orf_start
                    #print orf_gene
                    #print "In Open Frame of CDS"
                    yes += 1
                    orfs[counter] += "\n"+sign+" In CDS ORF, exonic length is: "+str(exonic_length)#+" start:"+str(placeStart)
                else:
                    #print orf_start
                    #print orf_gene
                    #print "NOT in Open Frame of CDS"
                    orfs[counter] += "\n"+sign+" Not in CDS ORF, exonic length is: "+str(exonic_length)#+" start:"+str(placeStart)
                    no += 1
                orflen = str(len(orf[1])) if orf[1][-1] == '*' else "INF"
                orfs[counter] += "\n"+sign+" The ORF is of length: "+orflen
    output.write("\n\n".join(orfs))
    output.close()
    print str(yes)+" yes and "+str(no)+" are no!"
    print "Meaning "+ str(float(yes)/float(yes+no)) + "% are in CDS ORF."
    print "FINISHED!"

def main(args):
    extract_ORF(args)


def check_orfs(args):
    orfs = open(args.output, 'r').read().strip('\n\n').split('\n\n')
    fasta = open(args.system, 'r').read().split('>')
    index = [chrom.strip('\n').split('\n')[0] for chrom in fasta]
    index = index[1:]
    fasta = fasta[1:]
    numOKStart = 0
    numOKEnd = 0
    numNoStart = 0
    numNoEnd = 0
    edStart = 0
    fasta = [fasta[i][len(index[i])+1:].replace('\n', '') for i in range(len(fasta))]
    for orf in orfs:
        orf = orf.split('\n')
        if len(orf) < 3:
            continue
        gstart = int(orf[3].split(':')[1].strip())
        gend = int(orf[4].split(':')[1].strip())
        gchrm = orf[0].split(';')[2].split(':')[0].strip()
        sign = orf[2][0]
        inf = (orf[1][-1] != '*')
        if sign == '+':
            ok = True
            if fasta[index.index(gchrm)][gstart-1:gstart+2].upper() == 'ATG':
                numOKStart += 1
                ok = ok and True
            else:
                numNoStart += 1
                ok = ok and False
            ending = fasta[index.index(gchrm)][gend:gend+3]
            if ending.upper() == 'TAA' or ending.upper() == 'TGA' or ending.upper() == 'TAG' or inf:
                numOKEnd += 1
                ok = ok and True
            else:
                numNoEnd += 1
                ok = ok and False
        else:
            ok = True
            starting = fasta[index.index(gchrm)][gend-3:gend]
            if starting.upper() == 'TTA' or starting.upper() == 'TCA' or starting.upper() == 'CTA' or inf:
                numOKStart += 1
                ok = ok and True
            else:
                numNoStart += 1
                ok = ok and False
            ending = fasta[index.index(gchrm)][gstart-4:gstart-1]
            if ending.upper() == 'CAT':
                numOKEnd += 1
                ok = ok and True
            else:
                numNoEnd += 1
                ok = ok and False
    percentStart = 100.0*float(numOKStart)/float(numOKStart+numNoStart)
    percentEnd = 100.0*float(numOKEnd)/float(numOKEnd+numNoEnd)
    print "           Starts"
    print "Correct | Incorrect | % Correct"
    print numOKStart, "     ", numNoStart, "        ", percentStart
    print "            --"
    print "           Ends"
    print "Correct | Incorrect | % Correct"
    print numOKEnd, "     ", numNoEnd, "        ", percentEnd
    print
    if percentStart < 90:
        print "WARNING! Many of the ORF start locations found are incorrect."
    if percentEnd < 90:
        print "WARNING! Many of the ORF stop locations found are incorrect."
    print "          ~~~~~~~"
    if percentStart > 90 and percentEnd > 90:
        print "CHECK PASSED SUCCESSFULLY!"
    else:
        print "CHECK FAILED!"
    print "          ~~~~~~~"
    print
    print "Some ORF starts/ends fall on the junction, hence, this check may consider them falsely as incorrect."
    print "On average more than 99% means no actual errors where observed. (test threshold = 90%)"


def only_exact(args):
    exons = open(args.exons, 'r').read().strip('\n').split('\n')
    output = open(args.exons.replace('.bed','_autoFiltEx.bed'), 'w')
    for exon in exons:
        exon_s = exon.strip('\t').split('\t')
        if exon_s[3].split(':')[1] not in anot_dict:
            continue
        all_from_annot = anot_dict[exon_s[3].split(':')[1]]
        line = all_from_annot[-1]
        for elem in all_from_annot:
            if exon_s[1] in elem.split()[8].strip(',').split(','):
                if exon_s[2] in elem.split()[9].strip(',').split(','):
                    output.write(exon+'\n')
                    break
    args.exons = args.exons.replace('.bed', '_autoFiltEx.bed')
    output.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-circ","-ci", "-c", help="A file with the circles locations.", default=None)
    group.add_argument("-exons", "-ex", "-e", help="Exact circles exons.", default=None)
    parser.add_argument("--force_exact", help="force the circles to be exact (If you suspect the circle exons"
                                             " might contain non-EXACT circles)", action='store_true')
    parser.add_argument("--collapse_frames", help="force the circles to be exact (If you suspect the circle exons"
                                             " might contain non-EXACT circles)", action='store_true')
    parser.add_argument("-system","-sys","-Sys", help="Full path to genome fasta file to use.", required=True)      
    parser.add_argument("-output", "-o", "-O", help="Output file name", default='cons_stats', required=False)
    parser.add_argument("-annotation", "-a", "-A", help="Annotation file.", required=True)
    args = parser.parse_args()
    #open file
    anot_dict = dict()
    try:
        anot_file = open(args.annotation, 'r')
    except:
        print bcolors.FAIL + "Can\'t open the file \"annotation\"..." + bcolors.ENDC
        exit(2)

    # acquire the number of exons and the starts and ends of them, for the specified gene.
    for line in anot_file.readlines():
        if line.split()[4] not in anot_dict:
            anot_dict[line.split()[4]] = [line]
        else:
            anot_dict[line.split()[4]].append(line)
    # start main:
    os.system("clear")
    if not args.exons:
        print "-----------------------------"
        print "No exons given. Starting from cRNA locations"
        print "-----------------------------"
        print "Finding circular exons..."
        find(args)
        args.exons = 'exons_EXACT.bed'
        print "-----------------------------"
        print
    else:
        print "-----------------------------"
        print "Exons given."
        print "-----------------------------"
        print
    if args.force_exact:
        print "Filtering only exact circles... (force_exact: ON)"
        only_exact(args)
        print "-----------------------------"
        print
    print "Extracting cRNA sequences from genome..."
    print "-----------------------------"
    get_fasta(args)
    print
    args.fasta = 'getFastaOutput.fa'
    print "-----------------------------"
    print "Searching for ORFs...", ("Reporting largest ORF among all frames. (collapse_frames: ON)" if args.collapse_frames
                                    else "Reporting ORFs in all frames. (collapse_frames: OFF)")
    main(args)
    print "-----------------------------"
    print
    print "Analyzing CDS data..."
    main2(args)
    print "-----------------------------"
    print
    print "Annotating in/out frame ORFs..."
    print "Creating .GARY file..."
    main3(args)
    print "-----------------------------"
    print
    print "Performing sanity checks on ORFs starts/stops..."
    check_orfs(args)
    print "-----------------------------"
    print
    print "ALL DONE! Thank you, come again :)"
    print "-----------------------------"
