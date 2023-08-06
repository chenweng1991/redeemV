"""
Script Name: FastQ to SQLite Processor
Author: Chen Weng
Date: 2023-07-29
Description: This script processes an indexed gzipped fastq file (i5) and stores the data in an SQLite database.
             It includes functions to reverse, complement, and reverse-complement DNA sequences, 
             and the main function reads the fastq file and writes the data into a SQLite database.

Usage: python MakeBC.SQLite.Nextseq_Nova1.5.py [path where the SQLite db will be write] [Path_To_I2_fqfile]
"""

import sys
import gzip
import sqlite3
import argparse

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
def main(dbpath,I2_fqfile):
    """ Create a SQLite database on disk to the path and database is named BC_database.db 
        And create a table named kv
    """
    conn = sqlite3.connect(dbpath+'/BC_database.db')
    c = conn.cursor()
    c.execute('CREATE TABLE kv (key text, value text)')
    with gzip.open(I2_fqfile) as f:
        while True:
            L1=f.readline().strip()
            L2=f.readline().strip()
            L4=f.readline().strip()
            if L1==b"":
                break
            ReadName=L1.split()[0].decode()
            Sequence=L2
            CellBC=reverse_complement(Sequence[8:24].decode())
            c.execute('INSERT INTO kv VALUES (?, ?)', (ReadName, CellBC))
            #sys.stdout.write(ReadName+"\t"+CellBC+"\n")
    c.execute('CREATE INDEX idx_key ON kv (key)') 
    conn.commit()
    conn.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process the indexed gzipped fastq file (i5) and store data in an SQLite database.')
    parser.add_argument('dbpath', help='Path to the directory where the database file will be stored.')
    parser.add_argument('I2_fqfile', help='Path to the gzipped fastq file.')

    args = parser.parse_args()

    main(args.dbpath,args.I2_fqfile)

