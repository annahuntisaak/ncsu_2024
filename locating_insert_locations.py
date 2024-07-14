import math

#takes string and returns string that doesn't have letters (or parentheses)
def remove_letters(string):
    filtered_line = ''
    #list of all characters that I want to keep around for now (all numbers, periods separating CDS start and end, spaces)
    letters = ['0','1','2','3','4','5','6','7','8','9','.',' ']
    for character in string:
        for i in letters:
            if character == i:
                filtered_line += character
    #should return string that doesn't have additional words
    return filtered_line.strip()


#idea: add characters to string until reaching a space, then stop, convert string to integer, and return that value
#along the way, update string to not include that character
def retrieve_locus_value(string):
    value = ''
    while string[0] != ' ':
        value += string[0]
        #string is now everything after this character, chopping off this character from it
        string = string[1: ]
    return int(value), string


#similar to retrieve_locus_value function, but now we care about hitting a period to separate first and second numbers
def retrieve_CDS_values(string):
    value_one = ''
    value_two = ''
    period_reached = False
    for character in string:
        #determining when we hit a period
        if character == '.':
            period_reached = True
        #prior to reaching period, add all characters to value 1
        if period_reached == False:
            value_one += character
        #after hitting a period, as long as we're not on a period character itself, add to value 2
        if (period_reached == True) & (character != '.'):
            value_two += character

    return int(value_one), int(value_two)


#function to get the name for the gene out of the line
def retrieve_gene_name(string):
    name = ''
    equal_sign_reached = False
    #only adding characters to our gene_name string once they follow the = symbol
    for character in string:
            if character == '=':
                equal_sign_reached = True
            elif equal_sign_reached == True:
                name += character
    #getting rid of the "" to clean up string
    return name.replace('"', '')


#function that reads through annotated Prokka file to retrieve list containing key details
#each element of the list is a list itself of the form: [locus #, locus length, CDS start, CDS end, complement strand or not, gene name]
def create_CDS_list(prokka_filepath):
    data = open(prokka_filepath, 'r')
    all_lines = data.readlines()
    all_key_data = []
    
    locus_line = ''
    CDS_line = ''
    gene_line = ''
    
    #important variables to keep track of, each will be an element of the specific CDS's list
    locus_number = 0
    locus_length = 0
    CDS_start = 0
    CDS_end = 0
    complement = False
    gene_name = ''

    #checking that the find() is working correctly
    counter_locus = 0
    counter_CDS = 0
    counter_gene = 0
    counter_nongene = 0

    #iterating through all lines
    for line in all_lines:
        #want to isolate first (locus #) and second (locus length) numbers
        if line.find('LOCUS') != -1:
            counter_locus += 1
            locus_line = line

            #getting locus line that only has numbers and spaces
            filtered_locus_line = remove_letters(locus_line)
            
            #getting first number in string, locus number
            locus_number, cut_locus_line = retrieve_locus_value(filtered_locus_line)

            #taking away the leading spaces from the string now that the characters of locus number have been removed
            cut_locus_line = cut_locus_line.strip()

            #getting second number in string, locus length
            locus_length, date_remaining = retrieve_locus_value(cut_locus_line)


        #determines if there's a CDS, and if so, creates list (element) to be added to overall list
        #adds both the locus and CDS lines to this new list
        #note: also had to account for the single case where a gene was identified under "tmRNA" instead of "CDS"
        if (line.find(' CDS ') != -1) or (line.find(' tmRNA ') != -1):
            counter_CDS += 1
            CDS_line = line

            #want to determine if it's complement or not before stripping away all the letters
            #if line includes complement
            if line.find('complement') != -1:
                complement = True
            else:
                complement = False
            
            #can now remove all letters from this line
            filtered_CDS_line = remove_letters(CDS_line)

            CDS_start, CDS_end = retrieve_CDS_values(filtered_CDS_line)
            
            #adding all the relevant elements to this list
            key_data = [locus_number]
            key_data.append(locus_length)
            key_data.append(CDS_start)
            key_data.append(CDS_end)
            key_data.append(complement)

            all_key_data.append(key_data)
    
        #updates gene each time a new one is encountered
        #we always want to add the gene name in, so add it to the last element of overall list
        if line.find(' /gene=') != -1:
            counter_gene += 1
            gene_line = line.strip()
            gene_name = retrieve_gene_name(gene_line)

            #adding gene name as the last entry in the last element of the list (i.e. most recently created list, should be corresponding CDS)
            all_key_data[-1].append(gene_name)
        
    # print(counter_locus)
    # print(counter_CDS)
    # print(counter_gene)

    data.close()

    #going through our total list now, looking at each element (which are also lists)
    #any list that doesn't have 6 elements (locus #, locus length, CDS start, CDS end, complement, gene name)
    #must be missing gene name, so assign to it N/A at that last index
    for list in all_key_data:
        if len(list) == 5:
            counter_nongene += 1
            list.append('N/A')

    return all_key_data




#now want to find all relevant positions using the all_key_data list
#IMPORTANT: format of all lists is [locus #, locus length, CDS start, CDS end, comp, gene]

#function that takes as input a complete list of all the CDS (both known and hypothesized genes) in the genome and a minimum desired interval size
#returns a list of positions (and their sizes) that can be found between 2 CDS (on opposite strands where neither is a KNOWN gene)
#each position is defined by its upstream and downstream CDS (so 2 lists returned for each found position)

#constant to set the minimum distance we want between potential CDS
minimum_interval = 1300

def find_possible_intervals(all_key_data, minimum_interval):
    previous_locus = all_key_data[0][0]
    previous_locus_length = all_key_data[0][1]
    previous_CDS_end = all_key_data[0][3]
    previous_strand = all_key_data[0][4]

    #elements of this list will be lists themselves
    #format of those lists is: [index of first CDS, index of second CDS], where indices are in terms of the all_key_data list
    #area betwen first and second is what we care about
    position_indices = []

    index = 0
    for list in all_key_data:
        #checking if we're in the same locus as previous, if this one's start is ahead of last one's end by minimum amount,
        #and if they're on opposite strands
        if (list[0] == previous_locus) & (list[2] - previous_CDS_end > minimum_interval) & (list[4] != previous_strand):
            position_indices.append([index-1, index])

        #updating variables to match this list before moving to next one
        previous_locus = list[0]
        previous_locus_length = list[1]
        previous_CDS_end = list[3]
        previous_strand = list[4]
        
        index += 1

    #further filtering, removing any entries in the position_indices list that have known, labeled genes involved
    no_gene_positions = []
    for pair in position_indices:
        #only keeping entire where both of the CDS bounding the region don't have known genes associated with them
        if (all_key_data[pair[0]][5] == 'N/A') & (all_key_data[pair[1]][5] == 'N/A'):
            no_gene_positions.append(pair)


    #last bit of clean-up, including size of interval for each interval
    positions_and_sizes = []
    size = 0
    for position in no_gene_positions:
        #if the CDS are in same locus
        if all_key_data[position[0]][0] == all_key_data[position[1]][0]:
            #interval between them is just difference between second one's start - first one's end
            size = all_key_data[position[1]][2] - all_key_data[position[0]][3] - 1
        
        positions_and_sizes.append([all_key_data[position[0]], all_key_data[position[1]], size])
    
    return positions_and_sizes




#want to look for genes that we don't want (antibiotic resistance, not important for cell survival, multiple copies, etc.)

#function that takes as input a list of lists with our key data and outputs a list with only 
#note: input list must have elements that are all themselves lists with 6 elements, where the last is the gene name (or 'NA')
def gene_list(all_key_data_list):
    all_gene_key_data = []
    gene_containing = 0
    for entry in all_key_data_list:
        if entry[5] != 'N/A':
            all_gene_key_data.append(entry)
            gene_containing += 1

    return all_gene_key_data


#function that helps us find instances of specific genes which we know by name
#inputs: all_genes is a list containing info about all genes in genome (i.e. list of lists that are of the form [locus #, locus length, CDS start, CDS end, comp, gene]),
#and genes_to_find is a list of strings where each entry is a specific gene name or core gene name
#returns the list of genes that were successfully found
def find_specific_genes(all_genes, genes_to_find):
    genes_found = []
    for gene in all_genes:
        for to_find in genes_to_find:
            if gene[5].find(to_find) != -1:
                genes_found.append(gene)
    return genes_found



#function that returns only the core part of a gene's name (i.e. removes everything after and including underscore)
def core_name(gene_name):
    #splitting gene name into sections before and after underscore
    gene_name_segments = gene_name.split("_")
    core_gene_name = gene_name_segments[0]
    return core_gene_name


#function that finds genes that have multiple copies (i.e. same name appears more than once)
#takes as input a list of all of the genes in the genome (and their associated key data)
#returns dictionary with gene names as the keys and list of lists of all of instances of that gene as the values
def find_multiple_copies(all_gene_key_data):
    multiple_copies = {}

    fully_checked = 0
    set_to_check = []
    current_gene_name = ''
    comparison_gene_name = ''

    for current_gene in all_gene_key_data:
        current_gene_name = current_gene[5]
        
        #if a gene name already exists as key we know previous iteration captured instances of that gene, can skip past it
        if current_gene_name in multiple_copies:
            fully_checked += 1
            continue
        #if gene name hasn't shown up before, check all other genes to see if there are matches with it
        else:
            set_to_check = all_gene_key_data[fully_checked + 1: len(all_gene_key_data)]

            for comparison_gene in set_to_check:
                comparison_gene_name = comparison_gene[5]

                #retreiving core names of the current and comparison genes
                core_comparison_name = core_name(comparison_gene_name)
                core_current_name = core_name(current_gene_name)

                #checking whether current and comparison genes have same core names
                if core_comparison_name == core_current_name:

                    #checking if we've already found at least 1 match with our current gene name
                    #if we have, means gene name should be a key in dictionary
                    if core_current_name in multiple_copies:
                        #updating relevant list to now include this additional instance of the gene
                        multiple_copies[core_current_name].append(comparison_gene)

                    #means this is the first match we've found to our current gene name
                    else:
                        multiple_copies[core_current_name] = [current_gene, comparison_gene]
                
            fully_checked += 1

    return multiple_copies
    



#want to get sequence of larger intervals around/including our identified positions

#recall: we previously created a list of non-gene-bounded positions that we care about (positions_and_sizes list)
#in this list, each element was a list itself that included the indices of the upper and lower CDS bounds (in terms of all_key_data list)
#and the total length of the interval
#also have specific genes we care about in list obtained from find_specific_genes

#we want at least 1.5 kb (1500) on either side of our identified position's start and end, setting this as a variable 
desired_length_on_ends = 1500

#function that simply determines whether or not it's possible to bound our chosen position with sequences of the desired length on both ends
#inputs: list containing all the key data for the gene of interest, length we want on either end
#output: returns True if there's sufficient room on both ends, returns False if not
def bounding_possible_gene(gene_data, length):
    #if the start of the gene CDS is too close to the beginning
    if gene_data[2] < length:
        return False
    #if there isn't enough space at the end (difference between total locus length and where gene ends)
    elif (gene_data[1] - gene_data[3]) < length:
        return False
    
    return True


#function that determines if intervals are boundable with sequence of desired length
#inputs: list containing the upstream and downstream info for an interval, length we want on either end
#output: returns True if there's sufficient room on both ends, returns False if not
def bounding_possible_nongene(interval_data, length):
    #if the upstream CDS's end is too close to the beginning
    if interval_data[0][3] < length:
        return False
    #if difference between total locus length and downstream CDS's start is too short
    elif (interval_data[1][1] - interval_data[1][2]) < length:
        return False
    
    return True


#function that returns a list where each entry records whether its possible to get a sequence of desired length using the corresponding
#entry from the original, input list
#input: list of identified positions (i.e. list of lists), length we want on either end
#output: list of same length as input, elements reflect if corresponding position from input list can be bounded
def bounding_each_position(position_list, length):
    bounding_positions = []
    possible = False
    for i in range(len(position_list)):
        #checking if the given position refers to a gene or an inter-CDS interval
        #if it's an interval, the first 2 elements of the position's list should be lists themselves (for the upstream and downstream CDS)
        if type(position_list[i][0]) == list:
            possible = bounding_possible_nongene(position_list[i], length)
        else: 
            possible = bounding_possible_gene(position_list[i], length)
        
        bounding_positions.append(possible)

    return bounding_positions


#function that writes the desired sequence to a file
#inputs: list with position data, length of bounding sequence, filepath containing sequence to be read, filepath to store desired sequence
#output: doesn't return anything, just writes to a text file
def writing_sequence(position, length, read_filepath, store_filepath):
    data = open(read_filepath, 'r')
    write_sequence = open(store_filepath, 'w')
    
    locus_found = False
    origin_found = False
    start_found = False
    end_found = False

    #calculating the correct sequence lines to write down
    #based on the fact that the sequence lines are labeled numerically, starting from 1 and increasing by 60 each line
    #also takes into account whether the given position is a gene or an inter-CDS interval

    #if our position is an inter-CDS interval
    if type(position[0]) == list:
        start = str((math.floor((position[0][3] - length) / 60) * 60) + 1)
        end = str((math.ceil((position[1][2] + length) / 60) * 60) + 1)
    #otherwise, we're looking at a specific, known gene
    else:
        start = str((math.floor((position[2] - length) / 60) * 60) + 1)
        end = str((math.ceil((position[3] + length) / 60) * 60) + 1)

    for line in data:
        if type(position[0]) == list:
            locus_number = 'contig_' + str(position[0][0])
        else:
            locus_number = 'contig_' + str(position[0])
        #looking for line that contains correct contig
        if (line.find('LOCUS') != -1) & (line.find(locus_number) != -1):
            locus_found = True
        #once we've gotten into the correct locus, looking for the origin
        elif (locus_found == True) & (line.find('ORIGIN') != -1):
            origin_found = True
        #if we're in correct locus and have found the origin, can start looking for correct sequence lines to record
        elif (locus_found == True) & (origin_found == True):
            if line.find(start) != -1:
                write_sequence.write(line)
                start_found = True
            elif (start_found == True) & (end_found == False):
                write_sequence.write(line)
                if (line.find(end) != -1):
                    end_found = True
            elif (start_found == True) & (end_found == True):
                break

    data.close()
    write_sequence.close()


#function that just writes the start and end positions of a given gene to a text file
#inputs: gene data, file path to text file to be written to
#output: doesn't return anything, just writes to file
def writing_start_end(position, filepath):
    start_end = open(filepath, 'w')
    if type(position[0]) == list:
        start_end.write(str(position[0][3] + 1) + '\n')
        start_end.write(str(position[1][2] - 1))
    else:
        start_end.write(str(position[2]) + '\n')
        start_end.write(str(position[3]))
    start_end.close()
