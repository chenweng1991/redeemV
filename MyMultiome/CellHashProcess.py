import sys
import gzip

##
fqR1=sys.argv[1]
fqR2=sys.argv[2]

##CellHash Dictionary
HashDic={
"TGTCTTTCCTGCCAG":"A0257",
"CTCCTCTGCAATTAC":"A0258",
"CAGTAGTCACGGTCA":"A0259",
"ATTGACCCGCGTTAG":"A0260",
"TAACGACCAGCCATA":"A0262",
"AAATCTCTCAGGCTC":"A0263",
"CTGTATGTCCGATTG":"A0264",
"TAAGATTCAGAGCGA":"A0265"}

##Define hamming_distance
def hd(s1, s2):
    return sum(ch1 != ch2 for ch1,ch2 in zip(s1,s2))

R1=gzip.open(fqR1,"rb")
R2=gzip.open(fqR2,"rb")

while True:
    R1Name=R1.readline()
    if R1Name==b"":
        break
    R1Seq=R1.readline()
    R1Opt=R1.readline()
    R1Qua=R1.readline()
    R2Name=R2.readline()
    R2Seq=R2.readline()
    R2Opt=R2.readline()
    R2Qua=R2.readline()
    CB=R1Seq[0:16]
    UMI=R1Seq[16:28]
    R2Hash=R2Seq[0:15]
    for Barcode in HashDic.keys():
        if(hd(R2Hash.decode(),Barcode)<=1):
            print(CB.decode()+"\t"+HashDic[Barcode]+"\t"+UMI.decode())
            break


R1.close()
R2.close()
