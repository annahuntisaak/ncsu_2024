#goal is to compare Prokka annotations of original genome and genome with desired insertion made
#want to identify any discrepencies in CDS (named genes or hypothesized genes)
#were any new CDS identified following insertion?  were any lost?

import math
from locating_insert_locations import *
from insert_sequence import *

original_genome_filepath = "/Users/annahuntisaak/Desktop/PreMiEr/E10.txt"
insert_genome_filepath = "/Users/annahuntisaak/Desktop/PreMiEr/ampE_insert_genome.txt"


#want to look at possible differences across the genomes in terms of hypothesized genes (CDS without explicit gene listed)

#function that takes as input an integer (index of a CDS's start or end) and a labeled sequence line, returns boolean indicating 
#if this is the correct line to begin extracting sequence information from
def check_to_start(desired_start, sequence_line):
    if (desired_start / 60).is_integer():
        desired_start -= 1
    line_start = (math.floor(desired_start / 60) * 60) + 1
    line_number = int(''.join(character for character in sequence_line if character.isdigit()))
    return (line_start == line_number)


#function that takes as input a CDS info list, which index of the list to use, desired sequence length, and sequence line
#returns the current state of desired sequence (string), number of bases left to be added (int), and whether sequence is complete (boolean)
def first_line_sequence_extraction(cds_list, index, desired_length, sequence_line):
    sequence_to_return = ""
    remaining = desired_length
    sequence_complete = False
    index_start = 0

    #getting the sequence alone (i.e. line number and spaces removed)
    sequence_from_current_line = ''.join(character for character in sequence_line if character.isalpha())

    #determining where in the sequence's string to go, depends on whether we're forming the start or end sequence
    line_number = int(''.join(character for character in sequence_line if character.isdigit()))
    #start sequence
    if index == 2:
        index_start = cds_list[2] - line_number
    #end sequence
    elif index == 3:
        index_start = (cds_list[3] - desired_length + 1) - line_number

    #checking that the entire desired sequence can be taken from this single line
    if (index_start + (desired_length - 1)) < (len(sequence_from_current_line)):
        sequence_to_return = sequence_from_current_line[index_start : (index_start + desired_length)]
        sequence_complete = True
        remaining = 0
    #case where our desired sequence will spill over into the next line
    else:
        sequence_to_return = sequence_from_current_line[index_start : len(sequence_from_current_line)]
        #recording the remainder left to be added to the sequence at the next line
        remaining -= (len(sequence_to_return))
    
    return sequence_to_return, remaining, sequence_complete


#function that takes as input the sequence to which we want to add (string), number of bases to be added,
#and the sequence line from which we will extract these remaining bases
#returns the completed sequence string, remainder to be added (should be 0), and indicator of the sequence being complete
def second_line_sequence_extraction(incomplete_sequence, to_add, sequence_line):
    sequence_to_return = incomplete_sequence
    remaining = 0
    sequence_complete = True

    sequence_from_current_line = ''.join(character for character in sequence_line if character.isalpha())
    sequence_to_return += sequence_from_current_line[0 : to_add]

    return sequence_to_return, remaining, sequence_complete


#function that takes as inputs a list containing key info about a CDS (of correct format) and filepath to the relevant annotated genome file
#outputs the first and last 12 bases of its sequence
def CDS_sequence_start_end(key_cds_info, filepath):
    sequence_start = ""
    sequence_end = ""

    locus_found = False
    origin_found = False
    sequence_start_found = False
    sequence_start_complete = False
    sequence_end_found = False
    sequence_end_complete = False

    desired_start_end_length = 12
    start_left = desired_start_end_length
    end_left = desired_start_end_length

    all_genome_data = open(filepath, 'r')

    while (sequence_start_complete == False) & (sequence_end_complete == False) :
        for line in all_genome_data:
            #encountering a locus line before we have identified the correct one to be in
            if (locus_found == False) & (line.find('LOCUS') != -1):
                filtered_locus_line = remove_letters(line)
                locus_number, rest_of_line = retrieve_locus_value(filtered_locus_line)
                #comparing this line's locus number to the one we want
                if locus_number == key_cds_info[0]:
                    locus_found = True
            #encountering the desired origin line within the correct locus
            elif (locus_found == True) & (origin_found == False) & (line.find('ORIGIN') != -1):
                origin_found = True

            elif (origin_found == True) & (sequence_start_found == False):
                #determining if this is the line to start recording the beginning of the CDS sequence
                if check_to_start(key_cds_info[2], line):
                    sequence_start_found = True
                    sequence_start, start_left, sequence_start_complete = first_line_sequence_extraction(key_cds_info, 2, desired_start_end_length, line)
            #filling in the rest of the start sequence, if needed
            elif (sequence_start_found == True) & (sequence_start_complete == False) & (start_left > 0):
                sequence_start, start_left, sequence_start_complete = second_line_sequence_extraction(sequence_start, start_left, line)

            #now creating the end sequence through a very similar process
            elif (sequence_start_complete == True) & (sequence_end_found == False):
                #figuring out where to start based on the end index - desired length this time
                if check_to_start(key_cds_info[3] - (desired_start_end_length - 1), line):
                    sequence_end_found = True
                    sequence_end, end_left, sequence_end_complete = first_line_sequence_extraction(key_cds_info, 3, desired_start_end_length, line)
            elif (sequence_end_found == True) & (sequence_end_complete == False) & (end_left > 0):
                sequence_end, end_left, sequence_end_complete = second_line_sequence_extraction(sequence_end, end_left, line)

    all_genome_data.close()
    return sequence_start, sequence_end


#function that takes as input 2 CDS lists and the appropriate filepaths needed to extract their corresponding sequences
#uses the CDS_sequence_start_end function to obtain the start and end of both CDS sequences --> compare them to determine if they're the same
#returns a boolean indicating whether the 2 CDS are equivalent (which we assume if their starts and ends match)
def CDS_same_or_not(cds_1, cds_2, filepath_1, filepath_2):
    #getting start and end of sequence for both CDS
    start_1, end_1 = CDS_sequence_start_end(cds_1, filepath_1)
    start_2, end_2 = CDS_sequence_start_end(cds_2, filepath_2)
    return (start_1 == start_2) & (end_1 == end_2)


#function that takes as input a complete list of a genome's CDSs and returns a reduced list of just the hypothesized genes
def hypothesized_gene_list(all_key_data_list):
    all_hypothesized_gene_key_data = []
    hypothesized_gene_containing = 0
    for entry in all_key_data_list:
        if entry[5] == 'N/A':
            all_hypothesized_gene_key_data.append(entry)
            hypothesized_gene_containing += 1

    return all_hypothesized_gene_key_data


#function that takes as input filepaths to the 2 files with our genomes of interest (original and insert-containing)
#systematically compares hypothesized genes and returns 2 lists – hypothesized genes that were lost and those that are new
def compare_hypothesized_genes(original_genome_filepath, insert_genome_filepath):
    #obtaining lists of just the hypothesized genes
    complete_original_CDS = create_CDS_list(original_genome_filepath)
    hypothesized_original = hypothesized_gene_list(complete_original_CDS)

    complete_insert_CDS = create_CDS_list(insert_genome_filepath)
    hypothesized_insert = hypothesized_gene_list(complete_insert_CDS)

    lost_hypothesized = []
    new_hypothesized = []

    located = False
    located_count = 0
    lost_count = 0

    #conducting pairwise comparison of these hypothesized gene lists, using that from the original genome as the "template"
    for hg_1 in hypothesized_original:
        located = False
        start_1, end_1 = CDS_sequence_start_end(hg_1, original_genome_filepath)
        hg_1_length = hg_1[3] - hg_1[2] + 1

        for hg_2 in hypothesized_insert:
            hg_2_length = hg_2[3] - hg_2[2] + 1

            #if they're in the same locus AND have the exact same length
            if (hg_1[0] == hg_2[0]) & (hg_1_length == hg_2_length):

                start_2, end_2 = CDS_sequence_start_end(hg_2, insert_genome_filepath)

                #determining if CDS are equal
                if (start_1 == start_2) & (end_1 == end_2):
                    #if we have found a matching hypothesized gene, remove it from the list we're checking against
                    hypothesized_insert.remove(hg_2)
                    located = True
                    located_count += 1

                if located == True:
                    break

        #seeing if we ever found a matching gene
        if located == False:
            #if we haven't found a match for this hypothesized gene in the new genome, we assume it was "lost"
            lost_hypothesized.append(hg_1)
            lost_count += 1

    new_hypothesized = hypothesized_insert
    return lost_hypothesized, new_hypothesized


#function that takes as input 2 filepaths to text files with genomes we want to compare
#first should be original genome, second should be the altered one
#returns 2 lists: genes that were "lost" (present in original and not in altered),
#genes that are "new" (weren't present in original, showed up in altered)
def compare_known_genes(original_genome_filepath, insert_genome_filepath):
    #getting lists of all CDS in both genomes
    original_CDS = create_CDS_list(original_genome_filepath)
    insert_CDS = create_CDS_list(insert_genome_filepath)

    #now, getting lists that only contain the known genes
    original_genes = gene_list(original_CDS)
    insert_genes = gene_list(insert_CDS)

    #lists to be returned at end
    lost_genes = []
    new_genes = []

    #going through all the genes in the original genome and seeing if we can find it in the altered genome
    found = False
    for gene in original_genes:
        found = False
        for possible_match in insert_genes:
            if gene[5] == possible_match[5]:
                found = True

                original_genes_counter += 1

                #removing this element from the altered genome's gene list
                insert_genes.remove(possible_match)

        #determining whether a match was found or not
        if found == False:
            lost_genes.append(gene)
    
    #all genes that are left in the insert_genes list are ones that weren't in the orignal genome
    #they weren't matched to any of the original genes, so they weren't removed
    new_genes = insert_genes

    return lost_genes, new_genes


#function that takes as input 2 lists: lost (known) genes, new (known) genes
#and 2 filepaths: original genome, insert genome
#checks whether any of the genes that have the same base name are actually the same
#returns 2 updated lists for the lost and new genes
def checking_lost_new_similarities(lost_genes, new_genes, original_genome_filepath, insert_genome_filepath):
    for possible_lost in lost_genes:
        lost_base_name = core_name(possible_lost[5])
        lost_start, lost_end = CDS_sequence_start_end(possible_lost, original_genome_filepath)
        for possible_new in new_genes:
            new_base_name = core_name(possible_lost[5])
            if lost_base_name == new_base_name:
                new_start, new_end = CDS_sequence_start_end(possible_new, insert_genome_filepath)
                if (lost_start == new_start) & (lost_end == new_end):
                    lost_genes.remove(possible_lost)
                    new_genes.remove(possible_new)
    
    return lost_genes, new_genes
