#!/usr/bin/env python3
###########################################################
#####IMPORT-STATEMENTS#####
import argparse
###########################################################
#####FUNCTION_DEFINITIONS#####
def avg_phred(qual_scores):
    """
    The function takes a character sequence of Illumina Quality scores
    and converts it into an average Quality Score for the sequence.
    """
    qual_scores = qual_scores.strip()
    N = len(qual_scores)
    Total = 0
    for char in qual_scores:
        phred_score = ord(char)-33
        Total += phred_score
    avg = Total/N
    return avg

def reverse_complement(sequence):
    """
    Takes a sequence of nucleotides (A,C,G,T) and returns 
    the reverse complement sequence.
    """
    reverse_sequence = sequence[::-1]
    comp_dict = {"A":"T","C":"G","G":"C","T":"A"}
    rev_comp = ""
    for nuc in reverse_sequence:
        rev_comp += comp_dict[nuc]
    return rev_comp

def hammax(a,b):
    """
    description: hammax takes 2 strings of equal length and returns
    their hamming distance.
    INPUT: 2 STRINGS
    OUTPUT: Hamming Distance <int>
    """
    if isinstance(a,str) != True:
        what = str(type(a))
        raise Exception("hammax takes STRINGS as input. Your 1st input was a {}.".format(what))
    if isinstance(b,str) != True:
        why = str(type(b))
        raise Exception("hammax takes STRINGS as input. Your 2nd input was a {}.".format(why))
    if len(a) != len(b):
        raise Exception("{} and {} must be of equal length".format(a,b))
    mismatch = 0
    x = len(a)
    for i in range(x):
        if a[i] == b[i]:
            continue
        else:
            mismatch+= 1
            continue
    return mismatch

def keepORnot(seq, a_set):
    """
    INPUT:  seq=A STRING
            a_set=AN ITERABLE OBJECT (set, list, dict, tuple) CONTAINING
                  A SET OF DISTINCT STRINGS.
    OUTPUT: RETURNS 
    """
    if isinstance(seq,str) != True:
        what = str(type(a))
        raise Exception("hammax takes a STRING as the first input. Your input was a {}.".format(what))
    if isinstance(a_set,set) != True:
        why = str(type(a_set))
        raise Exception("hammax takes a SET as the second input. Your input was a {}.".format(why))
    for s in a_set:
        if hammax(seq,s) < 2:
            return s
        else:
            continue
    return ""
###########################################################
#####ARGUMENT_PARSER#####
# Create the argument parser.
parser = argparse.ArgumentParser(
        description=
        """
        INPUT:  Forward & Reverse Biological FASTQs (required) <read1>, <read2>\n
                Forward & Reverse Dual Matched Index FASTQs (required) <index1>, <index2>\n
                TSV of Known Indexes (required) <known>\n
                An INTEGER: Filter reads whose average quality score is below <qc_cutoff> (optional, default=20)\n
        OUTPUT: For N known indexes, this script outputs\n
                (2N + 2) files:\n
                1 forward read FASTQ file for each index, containing reads w/ MATCHED indexes.\n
                1 reverse read FASTQ file for each index, containing reads w/ MATCHED indexes.\n
                1 forward read FASTQ file for UNKNOWN/MISMATCHED/LOW-QUALITY indexes.\n
                1 reverse read FASTQ file for UNKNOWN/MISMATCHED/LOW-QUALITY indexes.\n
                1 txt file reporting number of matched indexes per index, number of all indexes.\n
                that contain an "N", and number of observed index-hopping events per index.\n
        """)
# Inform the parser of command line options.
parser.add_argument('-r1', '--read1', help='Your forward read FASTQ file', 
                    type=str, required=True)
parser.add_argument('-r2', '--read2', help='Your reverse read FASTQ file',
                    type=str, required=True)
parser.add_argument('-i1', '--index1', help='Your forward index FASTQ file', 
                    type=str, required=True)
parser.add_argument('-i2', '--index2', help='Your reverse index FASTQ file',
                    type=str, required=True)
parser.add_argument('-k', '--known', help='A tsv file w/header, containing known index IDs in column 4\
                    and known index sequences in column 5.', type=str, required=True)
parser.add_argument('-c', '--qc_cutoff', help='An integer: Filter reads whose index quality-score is below \
                    --qc_cutoff <int>, default=20', type=int, default=20)
# Call the parser.
args = parser.parse_args()
# Save the qc_cutoff value in a variable cutoff.
cutoff = args.qc_cutoff
###########################################################
#####BEFORE OPENING ANY FILES#####
#####COUNTERS#####
# Initiate a counter for indexes with "N".
N_cnt = 0
# Initiate a counter for unmatched reads.
unmatched_cnt = 0
# Initiate a dictionary for low quality matches per index, key=INDEX, value=<count>.
low_qual = {}
# Initiate a dictionary for matches per index, key=INDEX, value=<count>
match_cnt = {}
# Initiate a dictionary for index-hopping-events (IHEs) per index, key=INDEX, value=<# of IHEs>.
I_hop = {}
# Initiate a dictionary for total reads per index, key=INDEX, value=<readcount>.
total_rds ={}
# Initiate a dictionary for true total reads, all indexes.
true_total = 0
# Initiate a counter for total low-quality-score reads.
QC = 0
# Initiate a counter for total index-hopping among all reads.
InHop = 0
###########################################################
#####DICTIONARIES#####
# Initiate an empty dictionary for MATCHING: Key=INDEX, Value=Reverse_Complement.
match_dict = {}
# Initiate a dictionary of FILENAMES to be created: Key=INDEX Value=[ <Read 1 FILENAME>, <Read 2 FILENAME> ].
output_files = {}
###########################################################
#####SETS#####
# Initiate an empty set for forward read indexes.
indexes = set()
# Initate an empty set for index reverse complements.
indexes_rev = set()
###########################################################
#####LISTS#####
# Initiate lists for input and output filenames, used with list comprehension, below.
output_fns = []
input_fns = [args.read1,args.index1,args.index2,args.read2]
unmatched = ['Unmatched_1.fq','Unmatched_2.fq']
###########################################################
#####OPEN-KNOWN-INDEX-FILE#####
# Open the text file containing the index sequences.
with open(args.known) as ID:
    # Iterate through each line in the text file.           
    for line in ID:
        # If the line is not the header, do the following:
        if line.startswith('sample'):
            continue
        # Strip whitespace and hidden characters from the line, then split by column
        # into a list, 'line'.
        line = line.strip().split()
        # Set the index to the variable "barcode".
        barcode = line[4]
        # Set the indexes reverse complement to "barcodeR".
        barcodeR = reverse_complement(barcode)
        # Add key=INDEX, value=# of Matches.
        match_cnt[barcode] = 0
        low_qual[barcode] = 0
        I_hop[barcode] = 0
        total_rds[barcode] = 0
        # Add key=INDEX, value=REVERSE_COMPLEMENT to match_dict.
        match_dict[barcode] = barcodeR
        indexes.add(barcode)
        indexes_rev.add(barcodeR)
        # Add key=INDEX, value=FILENAME to output_files.
        output_files[barcode] = [barcode+"_1.matched.fq",barcode+"_2.matched.fq"]
###########################################################
#####CREATE OUTPUT FILENAMES#####
# Create the list of output filenames to be used with list comprehension, below.
for i in output_files:
    for j in output_files[i]:
        output_fns.append(j)
###########################################################
#####LIST COMPREHENSION OPEN FILES#####
# Open input files for reading using list comprehension from input_fns list.
input_data = {fn_in: open(fn_in,'r') for fn_in in input_fns}
# Create and open output files for appending using list comprehension from output_fns list.
output_data = {fn_out: open(fn_out,'a') for fn_out in output_fns}
# Create and open output files for forward and reverse reads with unmatched indices.
unmatched_data = {um_out: open(um_out,'a') for um_out in unmatched}
###########################################################
#####PROCESS READS#####
while True:
    ##########
    # Read 1
    R1L1 = input_data[args.read1].readline().strip()
    # Break 'while loop' if line 1 is empty.
    if R1L1 == "":
        break
    R1L2 = input_data[args.read1].readline().strip()
    R1L3 = input_data[args.read1].readline().strip()
    R1L4 = input_data[args.read1].readline().strip()
    ##########
    # Index 1
    I1L1 = input_data[args.index1].readline().strip()
    # Break while loop if line 1 is empty.
    if I1L1 == "":
        break
    I1L2 = input_data[args.index1].readline().strip()
    I1L3 = input_data[args.index1].readline().strip()
    I1L4 = input_data[args.index1].readline().strip()
    ##########
    # Index 2
    I2L1 = input_data[args.index2].readline().strip()
    # Break while loop if line 1 is empty.
    if I2L1 == "":
        break
    I2L2 = input_data[args.index2].readline().strip()
    I2L3 = input_data[args.index2].readline().strip()
    I2L4 = input_data[args.index2].readline().strip()
    ##########
    # Read 2
    R2L1 = input_data[args.read2].readline().strip()
    # Break while loop if line 1 is empty.
    if R2L1 == "":
        break
    R2L2 = input_data[args.read2].readline().strip()
    R2L3 = input_data[args.read2].readline().strip()
    R2L4 = input_data[args.read2].readline().strip()
###########################################################
    # Increment true_total.
    true_total +=1
    # If "N" in either index, Increment N-count.
    if I1L2.count("N") > 0 or I2L2.count("N") > 0:
        N_cnt += 1
#######################################################
#####BEGIN MAIN CONTROL FLOW#####
    # Check quality-scores.
    if avg_phred(I1L4) < cutoff or avg_phred(I2L4) < cutoff:
        # No assumptions can be made if the quality-score
        # is below the quality-score cutoff. Indexes may 
        # be incorrect due to sequencing error.
        
        # If Indexes are still in sets, despite low QS:
        if I1L2 in indexes and I2L2 in indexes_rev:
            # If they match, record the index pair on header line,
            # And information about low QS on "+" line.
            if I2L2 == match_dict[I1L2]:
                # Increment QC.
                QC += 1
                # Increment unmatched.
                unmatched_cnt +=1
                # Increment total_rds.
                total_rds[I1L2] += 1
                # Increment low_qual.
                low_qual[I1L2] += 1
                # Add Index1:Index2 to Read Headers for R1 and R2.
                R1L1 += ":{}:{}".format(I1L2,I2L2)
                R2L1 += ":{}:{}".format(I1L2,I2L2)
                # Input NOTES on R1L3 and R2L3
                R1L3+=" Matched Index: AVG INDEX q-score<{:d}: Index1 Q-Score: {:.1f} Index2 Q-Score: {:.1f}".format(cutoff,
                                                                                                                    avg_phred(I1L4),
                                                                                                                    avg_phred(I2L4))
                R2L3+=" Matched Index: AVG INDEX q-score<{:d}: Index1 Q-Score: {:.1f} Index2 Q-Score: {:.1f}".format(cutoff,
                                                                                                                    avg_phred(I1L4),
                                                                                                                    avg_phred(I2L4))                
                # Create a list of the 4 lines to be joined, plus empty string to space.
                Read1temp = [R1L1,R1L2,R1L3,R1L4,""]
                Read2temp = [R2L1,R2L2,R2L3,R2L4,""]
                # Join the reads with newlines.
                Read1 = "\n".join(Read1temp)
                Read2 = "\n".join(Read2temp)
                # Write the reads to their unknown files.
                unmatched_data['Unmatched_1.fq'].write(Read1)
                unmatched_data['Unmatched_2.fq'].write(Read2)        
            # If both indexes are in the set, but not a match,
            # report index-hopping.
            else:
                # Increment QC.
                QC += 1
                # Increment InHop.
                InHop += 1
                # Increment total reads.
                total_rds[I1L2] += 1
                # Increment low_qual.
                low_qual[I1L2] += 1
                # Increment Unmatched.
                unmatched_cnt += 1
                # Increment the index-hop counter.
                I_hop[I1L2] += 1
                # Add Index to Read Headers for R1 and R2.
                R1L1 += ":{}:{}".format(I1L2,reverse_complement(I2L2))
                R2L1 += ":{}:{}".format(I1L2,reverse_complement(I2L2))
                # Input NOTES on R1L3 and R2L3
                R1L3 += " Mismatched Index: Observed Index-Hopping Event: {}:{}".format(I1L2,reverse_complement(I2L2))
                R2L3 += " Mismatched Index: Observed Index-Hopping Event: {}:{}".format(I1L2,reverse_complement(I2L2))
                Read1IHtemp = [R1L1,R1L2,R1L3,R1L4,""]
                Read2IHtemp = [R2L1,R2L2,R2L3,R2L4,""]
                Read1IH = "\n".join(Read1IHtemp)
                Read2IH = "\n".join(Read2IHtemp)
                # Write the reads to their unknown files.
                unmatched_data['Unmatched_1.fq'].write(Read1IH)
                unmatched_data['Unmatched_2.fq'].write(Read2IH)        
        # If 1 or both of the index-reads are NOT in the set(s):
        else:        
            # Check if some index reads are off by a hamming distance of 1.
            realIndex1 = keepORnot(I1L2,indexes)
            realIndex2 = keepORnot(I2L2, indexes_rev)
            # If hamming distance = 1:
            if realIndex1 != "" and realIndex2 != "":
                # If a match is retrieved:
                if realIndex2 == match_dict[realIndex1]:    
                    # Increment QC.
                    QC += 1
                    # Increment unmatched.
                    unmatched_cnt +=1
                    # Increment total_rds.
                    total_rds[realIndex1] += 1
                    # Increment low_qual.
                    low_qual[realIndex1] += 1
                    # Add INDEX1:INDEX2 to R1 and R2 headers.
                    # Also add message to line 3 of R1 and R2, informing user read has low quality-score.
                    R1L1 += ":{}:{}".format(realIndex1,realIndex2)
                    R1L3 += " Index read(s) contained 1 mismatch: Filtered for Low Quality-Score"
                    R2L1 += ":{}:{}".format(realIndex1,realIndex2)
                    R2L3 += " Index read(s) contained 1 mismatch: Filtered for Low-Quality-Score"
                    # Create a list of the 4 lines to be joined, plus empty string to space.
                    Read1temp = [R1L1,R1L2,R1L3,R1L4,""]
                    Read2temp = [R2L1,R2L2,R2L3,R2L4,""]
                    # Join the reads with newlines.
                    Read1 = "\n".join(Read1temp)
                    Read2 = "\n".join(Read2temp)
                    # Write the reads to their unknown files.
                    unmatched_data['Unmatched_1.fq'].write(Read1)
                    unmatched_data['Unmatched_2.fq'].write(Read2)

                else:
                    # Increment QC.
                    QC += 1
                    # Increment InHop.
                    InHop += 1
                    # Increment total reads.
                    total_rds[realIndex1] += 1
                    # Increment low_qual.
                    low_qual[realIndex1] += 1
                    # Increment Unmatched.
                    unmatched_cnt += 1
                    # Increment the index-hop counter.
                    I_hop[realIndex1] += 1
                    # Add Index to Read Headers for R1 and R2.
                    R1L1 += ":{}:{}".format(realIndex1,reverse_complement(realIndex2))
                    R2L1 += ":{}:{}".format(realIndex1,reverse_complement(realIndex2))
                    # Input NOTES on R1L3 and R2L3
                    R1L3 += " Mismatched Index: Observed Index-Hopping Event: {}:{}".format(I1L2,reverse_complement(realIndex2))
                    R2L3 += " Mismatched Index: Observed Index-Hopping Event: {}:{}".format(I1L2,reverse_complement(realIndex2))
                    Read1IHtemp = [R1L1,R1L2,R1L3,R1L4,""]
                    Read2IHtemp = [R2L1,R2L2,R2L3,R2L4,""]
                    Read1IH = "\n".join(Read1IHtemp)
                    Read2IH = "\n".join(Read2IHtemp)
                    # Write the reads to their unknown files.
                    unmatched_data['Unmatched_1.fq'].write(Read1IH)
                    unmatched_data['Unmatched_2.fq'].write(Read2IH)
            # Unknown file: If the hamming distance is greater than 1:
            else:
                # Increment QC.
                QC += 1
                # Increment unmatched.
                unmatched_cnt +=1
                # Input NOTES on R1L3 and R2L3
                R1L3 += " Unmatched Index: AVG INDEX q-score<{:d}: Index1 Q-Score: {:.1f} Index2 Q-Score: {:.1f}".format(cutoff,
                                                                                                                        avg_phred(I1L4),
                                                                                                                        avg_phred(I2L4))
                R2L3 += " Unmatched Index: AVG INDEX q-score<{:d}: Index1 Q-Score: {:.1f} Index2 Q-Score: {:.1f}".format(cutoff,
                                                                                                                        avg_phred(I1L4),
                                                                                                                        avg_phred(I2L4))
                # Create a list of the 4 lines to be joined, plus empty string to space.
                Read1QCtemp = [R1L1,R1L2,R1L3,R1L4,""]
                # Create a list of the 4 lines to be joined, plus empty string to space.
                Read2QCtemp = [R2L1,R2L2,R2L3,R2L4,""]
                # Join the lines into a read.
                Read1QC = "\n".join(Read1QCtemp)
                # Join the lines into a read.
                Read2QC = "\n".join(Read2QCtemp)
                # Write the reads to their unknown files.
                unmatched_data['Unmatched_1.fq'].write(Read1QC)
                unmatched_data['Unmatched_2.fq'].write(Read2QC)        
    # If both indexes are above cutoff and match known indexes.
    elif I1L2 in indexes and I2L2 in indexes_rev:
        # If it's a correct match, do the following:
        if I2L2 == match_dict[I1L2]:
            # Add Index and Index Reverse Complement to Read Headers for R1 and R2.
            R1L1 += ":{}:{}".format(I1L2,I2L2)
            R2L1 += ":{}:{}".format(I1L2,I2L2)
            # Create a list of the 4 lines to be joined, plus empty string to space.
            Read1temp = [R1L1,R1L2,R1L3,R1L4,""]
            Read2temp = [R2L1,R2L2,R2L3,R2L4,""]
            # Join the reads with newlines.
            Read1 = "\n".join(Read1temp)
            Read2 = "\n".join(Read2temp)
            # Increment total reads.
            total_rds[I1L2] += 1
            # Increment match_cnt for the Index.
            match_cnt[I1L2] += 1
            # Access output_files for filenames, save to variables.
            newread_1 = output_files[I1L2][0]
            newread_2 = output_files[I1L2][1]
            # Write the reads to the appropriate output files.
            output_data[newread_1].write(Read1)
            output_data[newread_2].write(Read2)
        # If index-hopping is observed, do the following:
        else:
            # Increment InHop.
            InHop += 1
            # Increment total reads.
            total_rds[I1L2] += 1
            # Increment Unmatched.
            unmatched_cnt += 1
            # Increment the index-hop counter.
            I_hop[I1L2] += 1
            # Add Index to Read Headers for R1 and R2.
            R1L1 += ":{}:{}".format(I1L2,reverse_complement(I2L2))
            R2L1 += ":{}:{}".format(I1L2,reverse_complement(I2L2))
            # Input NOTES on R1L3 and R2L3
            R1L3 += " Mismatched Index: Observed Index-Hopping Event: {}:{}".format(I1L2,reverse_complement(I2L2))
            R2L3 += " Mismatched Index: Observed Index-Hopping Event: {}:{}".format(I1L2,reverse_complement(I2L2))
            Read1IHtemp = [R1L1,R1L2,R1L3,R1L4,""]
            Read2IHtemp = [R2L1,R2L2,R2L3,R2L4,""]
            Read1IH = "\n".join(Read1IHtemp)
            Read2IH = "\n".join(Read2IHtemp)
            # Write the reads to their unknown files.
            unmatched_data['Unmatched_1.fq'].write(Read1IH)
            unmatched_data['Unmatched_2.fq'].write(Read2IH)
    # If both index-reads fail to match known indexes, check hamming distance:
    else:
        # Check if some index reads are off by a hamming distance of 1.
        realIndex1 = keepORnot(I1L2,indexes)
        realIndex2 = keepORnot(I2L2, indexes_rev)
        # If hamming distance = 1:
        if realIndex1 != "" and realIndex2 != "":
            # If a match is discovered:
            if realIndex2 == match_dict[realIndex1]:    
                # Add Index and Index Reverse Complement to Read Headers for R1 and R2.
                # Also add message to line 3 of R1 and R2, informing user read was salvaged.
                R1L1 += ":{}:{}".format(realIndex1,realIndex2)
                R1L3 += " Index read(s) contained 1 mismatch: Salvaged Record"
                R2L1 += ":{}:{}".format(realIndex1,realIndex2)
                R2L3 += " Index read(s) contained 1 mismatch: Salvaged Record"
                # Create a list of the 4 lines to be joined, plus empty string to space.
                Read1temp = [R1L1,R1L2,R1L3,R1L4,""]
                Read2temp = [R2L1,R2L2,R2L3,R2L4,""]
                # Join the reads with newlines.
                Read1 = "\n".join(Read1temp)
                Read2 = "\n".join(Read2temp)
                # Increment total reads.
                total_rds[realIndex1] += 1
                # Increment match_cnt for the Index.
                match_cnt[realIndex1] += 1
                # Access output_files for filenames, save to variables.
                newread_1 = output_files[realIndex1][0]
                newread_2 = output_files[realIndex1][1]
                # Write the reads to the appropriate output files.
                output_data[newread_1].write(Read1)
                output_data[newread_2].write(Read2)
            else:
                # Increment InHop.
                InHop += 1
                # Increment total reads.
                total_rds[realIndex1] += 1
                # Increment Unmatched.
                unmatched_cnt += 1
                # Increment the index-hop counter.
                I_hop[realIndex1] += 1
                # Add Index to Read Headers for R1 and R2.
                R1L1 += ":{}:{}".format(realIndex1,reverse_complement(realIndex2))
                R2L1 += ":{}:{}".format(realIndex1,reverse_complement(realIndex2))
                # Input NOTES on R1L3 and R2L3
                R1L3 += " Mismatched Index: Observed Index-Hopping Event: {}:{}".format(I1L2,reverse_complement(realIndex2))
                R2L3 += " Mismatched Index: Observed Index-Hopping Event: {}:{}".format(I1L2,reverse_complement(realIndex2))
                Read1IHtemp = [R1L1,R1L2,R1L3,R1L4,""]
                Read2IHtemp = [R2L1,R2L2,R2L3,R2L4,""]
                Read1IH = "\n".join(Read1IHtemp)
                Read2IH = "\n".join(Read2IHtemp)
                # Write the reads to their unknown files.
                unmatched_data['Unmatched_1.fq'].write(Read1IH)
                unmatched_data['Unmatched_2.fq'].write(Read2IH)
        else:
            # Increment unmatched.
            unmatched_cnt += 1
            # Input NOTES on R1L3 and R2L3
            R1L3 += " Unmatched Index: Ambigious Index(es) {}:{}".format(I1L2,I2L2)
            R2L3 += " Unmatched Index: Ambiguous Index(es) {}:{}".format(I1L2,I2L2)
            Read1UKtemp = [R1L1,R1L2,R1L3,R1L4,""]
            Read2UKtemp = [R2L1,R2L2,R2L3,R2L4,""]
            Read1UK = "\n".join(Read1UKtemp)
            Read2UK = "\n".join(Read2UKtemp)
            # Write the reads to their unknown files.
            unmatched_data['Unmatched_1.fq'].write(Read1UK)
            unmatched_data['Unmatched_2.fq'].write(Read2UK)

###########################################################
#####CLOSE ALL FILES#####
# Close input files.
for f in input_data:
    input_data[f].close()
# Close output files.
for f in output_data:
    output_data[f].close()
# Close unmatched files.
for f in unmatched_data:
    unmatched_data[f].close()
###########################################################
#####PRINT INFORMATION OUTPUT FILE#####
# Print table of informative output.
header = "Index Seq\tTotal\t%ofAll\tMatch#\tMatch%\tQC#\tQC%\tInHop#\tInHop%\n"
with open("deMULTI_output.tsv", mode='a') as t:
    # Write the header to the file.
    t.write(header)
    template = "{}\t{:d}\t{:.2%}\t{:d}\t{:.2%}\t{:d}\t{:.2%}\t{:d}\t{:.2%}\n"
    # Iterate through indexes, writing each indexes stats to file.
    for key in match_dict:
        this_index = key
        total_reads = total_rds[key]
        per_tot = total_rds[key]/true_total
        num_matches = match_cnt[key]
        if total_reads == 0:
            continue
        match_percent = num_matches/total_reads
        num_low = low_qual[key]
        low_percent = num_low/total_reads
        num_index_hop = I_hop[key]
        index_hop_percent = num_index_hop/total_reads
        t.write(template.format(this_index,total_reads,per_tot,num_matches,match_percent,num_low,low_percent,num_index_hop,index_hop_percent))
    all_unmatched = unmatched_cnt/true_total*100
    all_Ns = N_cnt/true_total*100
    all_low = QC/true_total*100
    all_hopping = InHop/true_total*100
    t.write(f"Total Number of Reads: {true_total}\n")
    t.write(f"Number of Unmatched Reads: {unmatched_cnt}\n")
    t.write(f"Percent of All Reads that were Unmatched: {all_unmatched:.2f}%\n")
    t.write(f"Number of Reads with an index containing an N: {N_cnt}\n")
    t.write(f"Percent of All Indexes Containing at least 1 N: {all_Ns:.2f}%\n")
    t.write(f"Number of Reads Removed for Low Quality-Scores: {QC}\n")
    t.write(f"Percent of All Reads with a Low Quality-Score: {all_low:.2f}%\n")
    t.write(f"Number of Reads Removed for Index-Hopping: {InHop}\n")
    t.write(f"Percent of All Reads with Index-Hopping: {all_hopping:.2f}%\n")

