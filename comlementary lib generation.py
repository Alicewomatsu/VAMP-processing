# -*- coding: utf-8 -*-
import csv
'''
Author: @ZeyiHuang
Finish date: 2023-6-1

Function: Create a library. In this library, you will have 2 sublibraries, 
and (part of) the sequences of these two sublibraries are complementary to each other 
these complementary bases are continuos
e.g. anticodon/codon
'''
'''
Degenrate bases table
N = [A, C, G, T]
B = [C, G, T]
D = [A, G, T]
H = [A, C, T]
V = [A, C, G]
K = [G, T]
M = [A, C]
R = [A, G]
S = [C, G]
W = [A, T]
Y = [C, T]
'''
def complementary_library_generation (sub1, sub2, c1_start,  c2_start, comp_len, head, body, tail):
    #sub1: the first sublibrary
    #sub2: the second sublibrary
    #c1_start: the first complementary base position in sub1
    #c2_start: the first complementary base position in sub2
    #comp_len: the length of complementary region
    #head: the 5' constant overhang before sub1
    #body: the insert region between sub1 and sub2. If no, please provide as ''.
    #tail: the 5' constant overhang after sub2
    #pay attention to the sequence direction. In this code, direction is 5'-sub1-sub2-3',
    # and both sub1 and sub2 are in 5'-3' direction of the genes.
    # This is a anticodon/codon library. c1_start complement to c2_end.

    Bases_table = {'N' : ['A', 'C', 'G', 'T'],
    'Y' : ['C', 'T'],
    'R': ['A', 'G'],
    'H' : ['A', 'C', 'T'],
    'A':['A'],'T':['T'],'C':['C'],'G':['G']}
    Complement_table = {'A':'T','T':'A','C':'G','G':'C'}
    #Here is the bases table you want to use in the lab
    lib_Sequence = sub1+sub2

    holder_1 = []
    holder_2 = []
    holder_3 = []
    holder_4 = []
    for base in Bases_table[lib_Sequence[0]]:
        holder_1.append(base)
    #create the list with all the randomized  bases based on the corresponding possibility at that position.
        # Loop through and add the next bases one at a time
    for l in range(1, len(lib_Sequence)):
        for seq in holder_1:
            for base in Bases_table[lib_Sequence[l]]:
                holder_2.append(seq + base)
        holder_1.clear()
        holder_1 = holder_2.copy()
        holder_2.clear()

    #revised the generated sequences so that they have the complementary region
    for seq in holder_1:
        comp_seq = [seq[e] for e in list(range(len(seq)))]
        # get the sequence into a list of each randomized sequence. Ready for revision.
        for i in range (0, comp_len):
            comp_seq[len(sub1)+c2_start+comp_len-2-i] = Complement_table[comp_seq[c1_start-1+i]]
        #Go through the randomized sequence till the first complementary pair. Then correct the bases in sub2 to complement for sub1
        #Until the whole complementary region length is over.
        comp_seq_strip = ''.join(comp_seq).strip()
        #Strip the list of this revised complement sequence back to a string.
        if comp_seq_strip not in holder_3:
            holder_3.append(comp_seq_strip)
        # Only append the not repeated sequence.
    for seq in holder_3:
        oligo = [seq[e] for e in list(range(len(seq)))]
        # get the sequence into a list. Ready for converting the library into the oligo for ordering.
        oligo.insert(len(sub1),body)
        #Insert the sequence between two sub library
        oligo.insert(0,head)
        #Insert the 5' sequence of sub1
        oligo.append(tail)
        #Insert the 3' sequence of sub2
        oligo_seq = ''.join(oligo).strip()
        holder_4.append(oligo_seq)

    holder_4 = sorted(holder_4)
    # Sort alphabetically
    return holder_4

head = 'CTTTATATATCTTGTGGAAAGGACGAAACACCGGGGGACGGTCcggcGACCaGCGGGT' #This is 3'U6-tRNA before HTS25-anticodon loop
body = 'ACCTaGCcagCGGGGttcgactCCCCGGTCTCTCgTTTTTTgctagCGGAAGAGCCGACGAGGCTCTTCcttgagcagaacaaacactccaagtgga'
#This is the sequence of 3'tRNA-SapI GGA cassette-VP3 sequence before lib
tail = 'accacgcagtcaaggcttcagttttctcaggcc'#VP3 sequence after lib
Quadruplet = complementary_library_generation('YTNNNNRH','NNNN',3,1,4,head,body,tail)

print(len(Quadruplet))
output_file = r'D:\PhD PROJECT\Experiment\AAV\Quadruplet codon\2023_06_01_Quadruplet_MaPyltR_G55T_AAV454_sequence.csv'
with open(output_file, 'w+', encoding='utf-8') as output_file:
    for i in range(len(Quadruplet)):
        output_file.write(f"{Quadruplet[i]}\n")
    output_file.close()