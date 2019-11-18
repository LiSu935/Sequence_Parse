#!/home/mico935/anaconda3/bin/python3.7
# --------- 2019-10-26 ----------#
# This file was originallt coppied from Exam.py, which was made to extract informations, \
# such as how many sequences are in a fasta file, what is the length of the longest sequence in the file.
# To make it easier, I am going to import SeqIO from Bio and to parse the sequences in the file.

import sys
import getopt
def usage():
    print("""
    Generate_Dic.py : reads a FASATA file and builds a dictionary with all sequences bigger \ 
    than a given length
    Generate_Dir.py [-h] [-l <length>] <filename>
        -h              print this message
        -l <length>     filter all seqs with a length small than <length> (default <length>=0)
        <filename>      the file has to be in FASATA format
        """)

o, a = getopt.getopt(sys.argv[1:], 'l:h') # o = list of optional arguments; a = list of required arguments
# when ':' is added just after, this means that the option expects a value
opts = {}
seqlen = 0;

for k,v in o:
    opts[k] = v

if '-h' in opts.keys():
    usage(); sys.exit()

if len(a) < 1:
    usage()
    sys.exit("input fasta file is missing")

if '-l' in opts.keys():
    if int(opts['-l']) < 0:
        print("Length of sequence should be positive!")
        sys.exit(0)
    seqlen=int(opts['-l'])

filename = sys.argv[-1]
print(filename)
try:
    f=open(filename)
except IOError:
	print("File does not exist, please check the path")
#seqs = {}
#for line in f:
#    # let's discard the newline at the end (if any)
#    line = line.rstrip()
#    #distinguish header from sequence
#    print(line)
#    if line.startswith('>'):
#        words=line.split()
#        name=words[0][1:]
#        seqs[name]=''
#    else: # sequence, not header
#        seqs[name] = seqs[name] + line
#print("The number of total sequences is %int" % len(seqs))
#seq_filted = {}
#for name,seq in seqs.items(): #items() method is used to retrieve the key and corresponding values from dic
#    if len(seq) >= seqlen:
#        seq_filted[name]=seq
#print(seq_filted)	
#print(len(seq_filted))
#seq_len = {}
#for name,seq in seqs.items():
#    seq_len[name]=len(seq)
#print(seq_len)
#print(sorted(seq_len.values()))
 
#f.close()

from Bio import SeqIO
seq_records = [seq_record for seq_record in SeqIO.parse(filename, "fasta")]
identifiers = [seq_record.id for seq_record in SeqIO.parse(filename, "fasta")]
#seqs       =  [seq_record.seq for seq_record in SeqIO.parse(filename, "fasta")]
seq_len     = [len(seq_record.seq) for seq_record in SeqIO.parse(filename, "fasta")]

print(identifiers)
print("Found %i records" % len(seq_records))
print(seq_len)
print(max(seq_len))

longest_length = max(seq_len)
for seq_record in seq_records:
    if len(seq_record.seq) == longest_length:
        longest_record = seq_record
        print(' longest_length =',longest_record.id)
        print(' longest_length', longest_length)
        print("the longest record in the file is with ID {}, with length {}" \
                .format(longest_record.id, longest_length))
        print("the longest record in the file is with ID %s, with length %i" \
                % (longest_record.id, longest_length))

seq_longest_orf_len = {}
for seq_record in seq_records:
    orf_len = {}
    for frame in range(3):
        length = 3 * ((len(seq_record.seq) - frame) // 3) 
        
        a = 0
    
        for pro in seq_record.seq[frame:frame+length].translate(1).split("*"):
            if 'M' in pro:
                pro = pro[pro.find('M'):]
               # print(pro)
                pos = seq_record.seq[frame:frame+length].translate(1).find('*',a)
                a = pos + 1
                # Because pos-codonposition+1 = len(pro) +1, codonposition = pos - len(pro), noticing that \
                        # pro does not contain the stop codon. so protein length should + 1.
                aa_start = (pos -len(pro)) * 3 + frame
                orf_len[aa_start] = (len(pro)+1) * 3
                print("%s has ORF in frame %s position %i with length %i" %\
                        (seq_record.id, frame, aa_start + 1, orf_len[aa_start]))
            else:
                pos = seq_record.seq[frame:frame+length].translate(1).find('*',a)
                a = pos + 1
    
    seq_longest_orf_length = max(orf_len.values())
    for name, value in orf_len.items():
        if value == seq_longest_orf_length:
            print(name)

    print("In position need to +1, the length of the longest ORF in %s is %i long"\
            % (seq_record.id, seq_longest_orf_length)) 
    seq_longest_orf_len[seq_record.id] = seq_longest_orf_length
file_longest_orf_len = max(seq_longest_orf_len.values())
for name, value in seq_longest_orf_len.items():
    if value == file_longest_orf_len:
        print(name)
print("has the longest ORF with length %s" % file_longest_orf_len) 

