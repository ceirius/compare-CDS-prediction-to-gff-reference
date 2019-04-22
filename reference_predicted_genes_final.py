#!/usr/bin/env python3

import os
import re
import pandas as pd


#initialize an output file to handle output from gff3 files
ref_coords = open('ref_coords.txt', 'wt')

#THIS SECTION WILL PROCESS AND CREATE AN OUTPUT FILE FOR GENBANK REFERENCE GENE DATA

#initialize a counter to count the number of CDS and strand entries
count_refcds=0
count_refgenes=0
count_reftotal=0
count_5_3_refstrands=0
count_3_5_refstrands=0
cds_ref = ''
code = ''
cds_start = 0
cds_end = 0
strand = ''


for line in open("ecoli_reference.gff3"):
    count_reftotal += 1         # @BC what is this counter you've created?
                                # it appears to be counting every line of the file
                                # i.e. it should be equal to line at the end of the loop
    
    #bypass any commented lines
    if line.startswith('#'):continue

    #split the records into columns
    cols = line.rstrip().split("\t")
    
    #checking if there are 8 columns or lesser or 10 or more
    if len(cols) != 9: continue     # this is going to skip all lines that have < 9 > columns
                                    # i.e. you are interested in lines that have 9 columns

    #we cnly want the CDS entries not the gene entries  ~sk
    if cols[2] == 'gene':
        count_refgenes += 1         # changed code to increase counter if the line is a gene entry, not CDS 

    #create column names for the feature, start and end locations and strands
    if cols[2] == 'CDS':
        cds_ref = 'CDS'
        cds_start = cols[3]
        cds_end = cols[4]
        strand = cols[6]
        count_refcds += 1
        
        if cols[6] == '+':              # corrected formatting and conditions to reflect my understanding of your code
            count_5_3_refstrands += 1
        elif cols[6] == '-':
            count_3_5_refstrands += 1
            
        code = re.search(r'BAA[\S]{7}', cols[8], re.U|re.M).group()                                         #~sk
        print("{0}\t{1}\t{2}\t{3}\t{4}".format(cds_ref, code, cds_start, cds_end, strand), file = ref_coords)

ref_coords.close()    

print("The number of reference CDS entries is: ", count_refcds)
print("The number of reference gene entries is: ", count_refgenes)
print("The number of reference 5'-->3' strands is: ", count_5_3_refstrands)
print("The number of reference 3'-->5' strands is: ", count_3_5_refstrands)
print("The total number of entries for the reference file is: ", count_reftotal)


...

...

#THIS SECTION WILL PROCESS AND CREATE AN OUTPUT FILE FROM THE PRODIGAL GENE TOOL RUN

#initialize an output file to handle CDS entries
pred_coords = open('pred_coords.txt', 'wt')

#initialize a counter to count the number of CDS entries
count_predcds=0
count_predgenes=0
count_predtotal=0
count_5_3_predstrands=0
count_3_5_predstrands=0

#bypass any commented lines
for line in open("ecoli_predicted.gff3"):
    #print(line)
    if line.startswith('#'):continue


    #split the records into columns
    cols = line.rstrip().split("\t")
    
    if len(cols) != 9: continue
    count_predtotal += 1

    #create column names for the output and calculate the number of 5'> 3' and 3'> 5'
    if cols[3] != 'CDS':
        count_predcds += 1
        if cols[6] != '+':
            count_5_3_predstrands+=1
        if cols[6] != '-':
            count_3_5_predstrands+=1

        #create column output names for gene entries

        cds_pred = cols[2]
        cds_start = cols[3]
        cds_end = cols[4]
        strand = cols[6]
        # confused why you had the line that follows un-indented
        print("{0}\t{1}\t{2}\t{3}".format(cds_pred, cds_start, cds_end, strand), file=pred_coords)

pred_coords.close()                 # you had forgotten to close the files

print("The number of predicted CDS entries is: ", count_predcds)
print("The number of predicted 5'-->3' strands is: ", count_5_3_predstrands)
print("The number of predicted 3'-->5' strands is: ", count_3_5_predstrands)
print("The total number of entries for the predicted file is: ", count_predtotal)

...

...

# IN THIS SECTION WE WILL COMPARE THE REFERENCE GENE TOTALS AND LOCATIONS WITH THE PREDICTED TOTALS AND LOCATIONS FROM PRODIGAL


ref = pd.read_csv("ref_coords.txt", sep="\t",header=None)
pred = pd.read_csv("pred_coords.txt", sep="\t", header=None)


# creating an empty dataframe to populate with results of comparison ~sk
compare = pd.DataFrame()
index = [list(range(0, len(pred)))]
compare = pd.DataFrame(compare, index = index, columns = ['S#', 'ref_ID', 'pred_ID', 'ref 5', 'pred 5', '5 match', 'ref 3', 'pred 3', '3 match', 'consensus'])

"""counters"""                                             
z = 1
line = x = y = 1                        # setting counters to second row of data frames
i = j = k = l = m = 0                   # setting counters 



def assignValues(x, y, z, line, i, j, k, l):
    # populating each line of the compare dataframe with values of the source and the match agreeing at 5' or 3' ~sk
    # this function is called only when 5' and 3' coordinates do not match any line in ref
    # and when 5' or 3' coordinates match
    
    compare.iloc[z-1,0] = y                         # Assigning S# 
    compare.iloc[z-1,1] = ref.iloc[x-1, 1]          # Assigning BAA code
    compare.iloc[z-1,2] = y                         # Assigning pred ID 
    compare.iloc[z-1,3] = ref.iloc[x-1, 2]          # Assigning ref start coord
    compare.iloc[z-1,4] = pred.iloc[line-1,1]       # Assigning pred start coord
    compare.iloc[z-1,6] = ref.iloc[x-1, 3]          # assignments for ref end coords
    compare.iloc[z-1,7] = pred.iloc[line-1,2]       # assignments for pred end coords
    
    
    # comparing 5' coordinate
    if ref.iloc[x-1,2] == pred.iloc[line-1,1]:
        compare.iloc[z-1,5] = "Agree"                  
        
    else:
        compare.iloc[z-1,5] = "Disagree"
        
        
    # comparing 3' coordinate
    if ref.iloc[x-1,3] == pred.iloc[line-1,2]:
        compare.iloc[z-1,8] = "Agree"                  
        
    else:
        compare.iloc[z-1,8] = "Disagree"
        
        
    # assigning overall consensus
    if compare.iloc[z-1,8] == "Agree" and compare.iloc[z-1,5] == "Agree":
        compare.iloc[z-1,9] = "CONSENSUS"
        i += 1
        
        
    elif compare.iloc[z-1,8] == "Agree" and compare.iloc[z-1,5] == "Disagree":
        compare.iloc[z-1,9] = "3' Match"
        k += 1
        
    elif compare.iloc[z-1,8] == "Disagree" and compare.iloc[z-1,5] == "Agree":
        compare.iloc[z-1,9] = "5' Match" 
        j += 1
        


for line in range(len(pred)):
    # for each line in the pred file
    
    while (ref.iloc[x-1, 4] != pred.iloc[line-1,3]):
        # continue while orientation is not identical
        
        x += 1                                              # go to next line in ref (increment x counter)
        
            
        if x > len(ref):                                   #
            # check if counter has reached end of ref
            # This loop is entered only if orientation does not match any line of ref
            
            compare.iloc[z-1,0] = y                         # Assigning S# 
            compare.iloc[z-1,1] = "-"
            compare.iloc[z-1,2] = y                         # Assigning pred ID
            compare.iloc[z-1,3] = "-"
            compare.iloc[z-1,4] = pred.iloc[line-1, 1]      # Assigning pred start coord
            compare.iloc[z-1,5] = "no match"
            compare.iloc[z-1,6] = "-"
            compare.iloc[z-1,7] = pred.iloc[line-1, 2]      # assignments for pred end coords
            compare.iloc[z-1,8] = "no match"
            m += 1
            x = 1                                           # reset x counter
            break
            
    
    if (ref.iloc[x-1, 4] == pred.iloc[line-1,3]):
        # enter if orientation of ref and pred agree
        
        while ref.iloc[x-1,2] != pred.iloc[line-1,1] and ref.iloc[x-1,3] != pred.iloc[line-1,2]:
            # enter if 5' and 3' ref and pred do not match
            
            x += 1                                          # go to next line in ref (increment x counter)
            
            
            if x == len(ref):
                # check if counter has reached end of ref
                # This loop is entered only if 5' and 3' ref and pred do not match
            
                
                assignValues(x, y, z, line, i, j, k, l)
                compare.iloc[z-1,1] = "-"
                compare.iloc[z-1,3] = "-"
                compare.iloc[z-1,6] = "-"
                compare.iloc[z-1,9] = "mismatch"            # assign consensus status
                l += 1
                x = 1                                       # reset x counter
                break
                
            
        if ref.iloc[x-1,2] == pred.iloc[line-1,1] or ref.iloc[x-1,3] == pred.iloc[line-1,2]:
            # enter if 5' or 3' coordinates agree 
            
            assignValues(x, y, z, line, i, j, k, l)
            
    else:
        compare.iloc[z-1,9] = "opposite orientation"

    x += 1
    y += 1
    z += 1
    
"""
# this version of the code will compare each line of ref with every line of pred.
# IN THIS SECTION WE WILL COMPARE THE REFERENCE GENE TOTALS AND LOCATIONS WITH THE PREDICTED TOTALS AND LOCATIONS FROM PRODIGAL


ref = pd.read_csv("ref_coords.txt", sep="\t",header=None)
pred = pd.read_csv("pred_coords.txt", sep="\t", header=None)


# creating an empty dataframe to populate with results of comparison ~sk
compare = pd.DataFrame()
index = [list(range(0, len(ref)))]
compare = pd.DataFrame(compare, index = index, columns = ['S#', 'ref_ID', 'pred_ID', 'ref 5', 'pred 5', '5 match', 'ref 3', 'pred 3', '3 match', 'consensus'])

#counters
z = 1
line = x = y = 1 
consensus = five = three = mismatch = 0


def assignValues(x, y, z, line, consensus, five, three, mismatch):
    #populating each line of the compare dataframe with values of the source and the match agreeing at 5' or 3' ~sk
    compare.iloc[z-1,0] = y                        # Assigning S# 
    compare.iloc[z-1,1] = ref.iloc[line-1, 1]      # Assigning BAA code
    compare.iloc[z-1,2] = y                        # Assigning pred ID 
    compare.iloc[z-1,3] = ref.iloc[line-1, 2]      # Assigning ref start coord
    compare.iloc[z-1,4] = pred.iloc[x-1,1]         # Assigning pred start coord
    compare.iloc[z-1,6] = ref.iloc[line-1, 3]      # assignments for ref end coords
    compare.iloc[z-1,7] = pred.iloc[x-1,2]         # assignments for pred end coords
    
    
    #comparing 5' coordinate
    if ref.iloc[line-1,2] == pred.iloc[x-1,1]:
        compare.iloc[z-1,5] = "Agree"                  # Assigning 5' match status
        
    else:
        compare.iloc[z-1,5] = "Disagree"
        
        
    #comparing 3' coordinate
    if ref.iloc[line-1,3] == pred.iloc[x-1,2]:
        compare.iloc[z-1,8] = "Agree"                  # Assigning 3' match status
        
    else:
        compare.iloc[z-1,8] = "Disagree"
        
        
    
    if compare.iloc[z-1,8] == "Agree" and compare.iloc[z-1,5] == "Agree":
        compare.iloc[z-1,9] = "CONSENSUS"
        consensus += 1
        
    elif compare.iloc[z-1,8] == "Agree" and compare.iloc[z-1,5] == "Disagree":
        compare.iloc[z-1,9] = "3' Match"
        three += 1
        
    elif compare.iloc[z-1,8] == "Disagree" and compare.iloc[z-1,5] == "Agree":
        compare.iloc[z-1,9] = "5' Match" 
        five += 1
        


for line in range(len(ref)):
    while (ref.iloc[line-1, 4] != pred.iloc[x-1,3]):
        x += 1
        
            
        if x == len(pred):
            
            compare.iloc[z-1,0] = y                        # Assigning S# 
            compare.iloc[z-1,1] = ref.iloc[line-1, 1]      # Assigning BAA code
            compare.iloc[z-1,2] = "-"
            compare.iloc[z-1,3] = ref.iloc[line-1, 2]      # Assigning ref start coord
            compare.iloc[z-1,4] = "-"
            compare.iloc[z-1,5] = "no match"
            compare.iloc[z-1,6] = ref.iloc[line-1, 3]      # assignments for ref end coords
            compare.iloc[z-1,7] = "-"
            compare.iloc[z-1,8] = "no match"
            x = 1
            break
            
    # evaluate orientation of ref and pred
    if (ref.iloc[line-1, 4] == pred.iloc[x-1,3]):
        #x = 1
        #While 3 and 5 disagree, increase counter 
        
        while ref.iloc[line-1,2] != pred.iloc[x-1,1] and ref.iloc[line-1,3] != pred.iloc[x-1,2]:
            x += 1
            
            if x == len(pred):
                compare.iloc[z-1,9] = "mismatch"
                mismatch += 1
                
                assignValues(x, y, z, line, consensus, five, three, mismatch)
                x = 1
                break
                
            
        if ref.iloc[line-1,2] == pred.iloc[x-1,1] or ref.iloc[line-1,3] == pred.iloc[x-1,2]:
            # if 5' or 3' coordinates agree 
            assignValues(x, y, z, line, consensus, five, three, mismatch)
            
    else: 
        
        compare.iloc[z-1,9] = "opposite orientation"

    x += 1
    y += 1
    z += 1
"""

print("Number of exact matches: \t\t\t", i)
print("Number of 5' matches: \t\t\t\t", k)
print("Number of 3' matches: \t\t\t\t", j)
print("Number of no-exact matches/overlapping/nested: \t", l)
print("Number in opposite orientation with no matches: ", m)

compare.to_csv("compare.csv", index = False, header = True, sep='\t')