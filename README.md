# Finding Possible Insertion Positions Within Genome and Insertion Analysis
*Developed in Summer 2024, during the PreMiEr REU*\
Anna Rose Hunt-Isaak, McGill University\
Personal Email: annahuntisaak@gmail.com\
Student Email: anna.hunt-isaak@mail.mcgill.ca\
\
**Important Note: All three files should be placed in the same folder on your computer to allow various import statements (and thus the use of functions across the files) to work properly.**

## File 1: locating_insert_locations
### Overview
This program is intended to help you find possible locations for insertion within a genome (which has been previously annotated via Prokka).  It allows you to perform three central tasks:
1) Identify any inter-CDS intervals that meet your desired minimum length threshold
2) Locate specific, individual genes and determine if the genome contains multiple copies of a given gene
3) Record the sequence of an appropriately-sized region containing your identified position as well as the start and end points of the position within said region

### Use of Important Functions
#### create_CDS_list
##### Input: file path to your Prokka annotated genome (should be a .txt file!)
##### Output: returns a list with information about each CDS in the genome
The format of this output is important for a variety of downstream functions.  Each element of the list that is returned is itself a list and contains the following information:\
[contig number *(int)*, total length of the contig *(int)*, CDS start *(int)*, CDS end *(int)*, if this CDS is on the complementary strand or not *(bool)*, gene name (if CDS is associated with a known gene) or N/A (if this is a hypothesized gene) *(str)*]\
\
This function should generally be called prior to any of the others contained within this file.

#### find_possible_intervals 
##### Inputs: list with information about each CDS in the genome, minimum desired length for an inter-CDS interval
##### Output: returns a list with information about all relevant intervals
Before calling this function, you should first call create_CDS_list → use its output as the first input for this function.  In addition, there is a constant named minimum_interval (located above this function’s definition) that can be changed to reflect your desired, minimum interval size and used as the second input.\
\
The inter-CDS intervals that will be selected and returned are those that fulfill the following requirements:
1) Located between two CDSs of the same contig that are on opposite strands
2) Located between two CDSs that are **not** associated with known genes
3) Is of a sufficient length
The output is a list in which each element summarizes the key information of an interval.  The format of these elements is: [ [information about the upstream CDS], [information about the downstream CDS], total length of the interval ]\
*Note: Lists containing information about the upstream and downstream CDSs are in the same format as those described previously in the create_CDS_list section.*

#### find_specific genes
##### Inputs: list with information about all CDSs (or just known genes) in the genome, list of gene names you want to search for
##### Output: returns a list with information about the genes that were successfully found
Before calling this function, you should first call create_CDS_list (and, optionally, can call gene_list to produce a list containing only CDSs associated with *known* genes).\
*Note: Function detects all genes whose names **contain** the name of interest (not just those that are an exact match).\
Multiple copies example: if you search for “geneX” and there are 2 copies of it within the genome → geneX_1 and geneX_2 will both be returned\
Similar names example: if you search for “abc” and “abcA” exists → all copies of “abc” **and** “abcA” will be returned*

#### find_multiple_copies
##### Input: list with information about all *known* genes
##### Output: returns a dictionary of the genes which appear to have multiple copies, each key is a gene name and each corresponding value is a list of all instances of that gene
Before calling this function, you should first call create_CDS_list → call gene_list.

#### bounding_each_position
##### Inputs: list of your positions of interest (either known genes or inter-CDS intervals), desired length for the sequences upstream and downstream of position
##### Output: returns a list where each element tells you if the corresponding position from the input list can be embedded in upstream and downstream sequences of the desired length (i.e. each element is a boolean)
The list you pass as an input to this function can be obtained by calling create_CDS_list → choosing one of the two following routes
1) If you want to focus on certain known genes: call gene_list → find_specific_genes
2) If you want to focus on inter-CDS intervals: call find_possible_intervals
There is a constant called desired_length_on_ends (located above all of the bounding-related function definitions, currently set to 1500 base pairs) that can be passed as the second input.

#### writing_sequence 
##### Inputs: single list with information about one position of interest, desired length for the upstream and downstream sequences, file path to your Prokka annotated genome (.txt file advised), file path indicating where you want to store this sequence (.txt file advised)
##### Output: doesn’t return anything, simply writes the desired sequence to a file
Before calling this function, you should call a bounding-related function to check that it is indeed possible to extract a continuous sequence containing your position and upstream/downstream segments of the desired length.\
*Note: This function has been written to identify lines containing “contig_#” to help determine where the desired sequence is located.  Any annotated genome file passed as input should thus use “contig_#” to label the different continuous sections.*

#### writing_start_end
##### Inputs: single list with information about one position of interest, file path indicating where you want to store its start and end points
##### Output: doesn’t return anything, simply writes the position’s start and end points to a file
Can be used alongside the writing_sequence function, allowing you to keep track of where the position is within the larger sequence being written down.


## File 2: insert_sequence
### Overview
The shortest of the three files, this program’s central purpose is the production of an updated genome sequence – i.e. contains your insert at the desired position.  This modified genome can subsequently be re-annotated and compared to the original genome.

### File Paths and Use Important Functions
#### File Paths
There are several file paths located at the top of this file that should be changed to suit your own.  They should be sufficient to represent all of the different files you may need access to when calling functions in this program.

#### write_insert_genome
##### Inputs: file path to your Prokka annotated file, file path to your desired insert’s sequence, file path to sequences that are upstream and downstream of your chosen position, file path indicating where you want to record the altered genome sequence
##### Output: doesn’t return anything, simply writes the down the new genome sequence (i.e. with insert sequence in the desired location) to a file
The file containing both the upstream and downstream sequences must have the following format:
1) Each sequence labeled accordingly, with the label in in all caps on the line directly above the sequence
2) Upstream and downstream sequences separated from each other by at least 1 blank line


## File 3: insertion_analysis
### Overview
This program allows you to compare two annotated genomes (original genome and post-insertion genome), identifying any losses or gains of both known and hypothesized genes.

### File Paths and Use of Important Functions
#### File Paths
There are two file paths located at the top of this file that should be changed to suit your own.

#### CDS_sequence_start_end
##### Inputs: single list with information about a single CDS, file path to Prokka annotated genome that contains the CDS’s sequence
##### Outputs: returns two strings, the beginning of the CDS’s sequence and the end
In general, this function doesn’t need to be called (beyond calls within other functions) or altered by the user.\
The only part of the function designed to be changed is its desired_start_end_length constant; this integer determines the number of base pairs that will be retrieved from the beginning and end of the CDS’s complete sequence. 

#### compare_hypothesized_genes
##### Inputs: file path to your Prokka annotated original genome, file path to your Prokka annotated altered genome (genome containing your insert at the correct location)
##### Outputs: returns two lists, one containing information about hypothesized genes that were “lost” (appear in original genome but *not* in the altered genome) and another containing information about hypothesized genes that are “new” (appear in altered genome but *not* in the original genome)
This function compares the beginning and end of two hypothesized genes’ sequences to determine if they are the same.  It only compares sequences if two hypothesized genes (one in the original genome, one in the altered genome) meet the following requirements:
1) Are located in the same contig
2) Have the exact same total length

#### compare_known_genes 
##### Inputs: file path to your Prokka annotated original genome, file path to your Prokka annotated altered genome (genome containing your insert at the correct location)
##### Outputs: returns two lists, one containing information about known genes that were “lost” (appear in original genome but *not* in the altered genome) and another containing information about known genes that are “new” (appear in altered genome but *not* in the original genome)
This function determines “lost” and “new” genes based solely on their names.  You should call checking_lost_new_similarities (see below) after this one to ensure that a gene which has simply been re-indexed (ex: label of the same gene slightly changed across the genomes, such as geneX to geneX_1) isn’t accidentally categorized as “lost” or “new”.

#### checking_lost_new_similarities
##### Inputs: list with information about “lost” known genes, list with information about “new” known genes, file path to your Prokka annotated original genome, file path to Prokka annotated altered genome
##### Outputs: returns updated versions of both the “lost” known genes list and the “new” known genes list
Should ideally be used directly after compare_known_genes.
