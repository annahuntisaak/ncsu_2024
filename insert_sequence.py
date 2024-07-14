import math
from locating_insert_locations import *

#6/21/24
#given the upstream and downstream sequences bounding our insert, want to replace the current intermediate region with our insert (GFP)

#filepath variables
locus_filepath = "/Users/annahuntisaak/Desktop/PreMiEr/E10.txt"
original_genome_filepath = "/Users/annahuntisaak/Desktop/PreMiEr/E10.txt"
insert_filepath = "/Users/annahuntisaak/Desktop/PreMiEr/insert.txt"
upstream_downstream_filepath = "/Users/annahuntisaak/Desktop/PreMiEr/ampE_up_down.txt"
genome_with_insert_filepath = "/Users/annahuntisaak/Desktop/PreMiEr/ampE_inserted.txt"


#function that takes string containing locus information and returns locus name (contig_#)
def locus_name(locus_string):
    locus_name = ''
    #using function written in the locating_insert_locations file
    no_letters = remove_letters(locus_string)
    #getting the locus number and the string containing all other numerical info
    locus_number, rest = retrieve_locus_value(no_letters)
    locus_name = 'contig_' + str(locus_number)

    return locus_name


#function that takes a line containing sequence we want and returns a version of that line
#in all caps, with any numbers, spaces removed, and newlines removed
def format_sequence_line(line):
    if line.find("//") != -1:
        return ""
    else:
        sequence_only = ""
        desired_characters = ['a', 'c', 'g', 't']
        for character in line:
            for i in desired_characters:
                if character == i:
                    sequence_only += character
        return sequence_only.upper()


#function that takes as input filepath to file with complete genome information (including loci info and genome sequence)
#returns the genome sequence (all caps), with lines that delineate different loci (> contig_#)
def genome_sequence(original_genome_filepath):
    genome_information = open(original_genome_filepath, 'r')
    original_genome = ""
    origin_line_found = False
    sequence_started = False

    for line in genome_information:
        if line.find('LOCUS') != -1:
            locus = locus_name(line.strip("\n"))
            original_genome = original_genome + "\n" + "> " + locus + "\n"
        #finding where the lines labeled ORIGIN are, after this line there will be the sequence we want
        elif line.find('ORIGIN') != -1:
            origin_line_found = True
        #finding line the immediately follows an ORIGIN line
        elif (origin_line_found == True) & (sequence_started == False):
            original_genome = original_genome + format_sequence_line(line)
            sequence_started = True
            #changing origin_line_found variable to False to reset it for next time
            origin_line_found = False
        elif (origin_line_found == False) & (sequence_started == True):
            if format_sequence_line(line) != "":
                original_genome = original_genome + format_sequence_line(line)
            else:
                sequence_started = False
    
    genome_information.close()
    return original_genome


#function that takes 1 input: filepath to file with contents of the complete sequence to be inserted
#returns string of the sequence in all caps
#note: requires that the inputted file only contains the sequence and nothing else
def format_insert(insert_file):
    unformatted_insert = open(insert_file, 'r')
    formatted = ""
    for line in unformatted_insert:
        formatted += (line.strip("\n")).upper()
    unformatted_insert.close()
    return formatted


#function that takes 2 inputs: filepath to file with upstream and downstream sequences, keyword string (all caps) that determines where desired content starts
#returns string of the entire desired sequence (i.e. upstream or downstream segment) in caps
#note: there should be at least one completely blank line that separates the upstream and downstream sections
def find_up_or_down(file, keyword):
    up_down_file = open(file, 'r')
    sequence = ""
    start_found = False
    end_found = False

    for line in up_down_file:
        #finding line that contains the keyword, "UPSTREAM" or "DOWNSTREAM"
        if (line.upper()).find(keyword) != -1:
            start_found = True
        #finding line that ends our desired segment (i.e. a completely blank line)
        elif (start_found == True) & (len(line.strip()) == 0):
            end_found = True
        elif (start_found == True) & (end_found == False):
            sequence += (line.strip("\n")).upper()

    up_down_file.close()
    return sequence


#function that takes as input 4 filepaths to files containing: 1) original genome, 2) insert sequence, 3) upstream and downstream sequences
#and 4) file where you want to write down the sequence of the genome with the insert
#doesn't return anything, just writes to the appropriate file
def write_insert_genome(original_genome_filepath, insert_filepath, upstream_downstream_filepath, genome_with_insert_filepath):
    
    #getting the entire original genome sequence
    original_genome = genome_sequence(original_genome_filepath)

    #getting the desired insert sequence in the correct format
    insert_sequence = format_insert(insert_filepath)

    #obtaining the correctly formatted upstream sequence
    upstream_sequence = find_up_or_down(upstream_downstream_filepath, "UPSTREAM")
    upstream_length = len(upstream_sequence)

    #getting the correctly formatted downstream sequence
    downstream_sequence = find_up_or_down(upstream_downstream_filepath, "DOWNSTREAM")
    downstream_length = len(downstream_sequence)

    #finding where the upstream and downstream sequences are located in the original genome
    upstream_start_index = original_genome.find(upstream_sequence)
    upstream_end_index = upstream_start_index + upstream_length - 1

    #don't need to know where end of downstream sequence is located since we only care about inserting between upstream and downstream
    downstream_start_index = original_genome.find(downstream_sequence)


    #can now piece together the final sequence we want
    #everything up to and including upstream (from original_genome) - sequence to insert - everything including and beyond downstream (from original_genome)
    genome_with_insert = ""
    genome_with_insert = original_genome[0:upstream_end_index+1]
    genome_with_insert = genome_with_insert + insert_sequence
    genome_with_insert = genome_with_insert + original_genome[downstream_start_index:]


    #writing this new string to a separate txt file
    altered_genome = open(genome_with_insert_filepath, 'w')
    altered_genome.write(genome_with_insert.strip("\n"))
    altered_genome.close