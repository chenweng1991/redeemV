import sys
import gzip


I2_fqfile=sys.argv[1]

def Mismatch(a,b):
    return sum(a[i]!=b[i] for i in range(len(a)))


def reverse(seq):
    """Returns a reversed string"""
    return seq[::-1]


def complement(seq):
    """Returns a complement DNA sequence"""
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
    seq_list = list(seq)
    # I can think of 3 ways to do this but first is faster I think ???
    # first list comprehension
    seq_list = [complement_dict[base] for base in seq_list]
    # second complicated lambda
    # seq_list = list(map(lambda base: complement_dict[base], seq_list))
    # third easy for loop
    # for i in range(len(seq_list)):
    #    seq_list[i] = complement_dict[seq_list[i]]
    return ''.join(seq_list)

def reverse_complement(seq):
    """"Returns a reverse complement DNA sequence"""
    seq = reverse(seq)
    seq = complement(seq)
    return seq


#I2_fqfile='Test_I2.fastq'
i=0
with gzip.open(I2_fqfile) as f:
    while True:
        L1=f.readline().strip()
        L2=f.readline().strip()
        L3=f.readline().strip()
        L4=f.readline().strip()
        if L1==b"":
            break
        ReadName=L1.split()[0].decode()
        Sequence=L2
        CellBC=reverse_complement(Sequence[0:16].decode())
        sys.stdout.write(ReadName+"\t"+CellBC+"\n")
