import pyfastx
import gzip
import csv
import sys
from progress.bar import Bar
##Get the fq file name
fqFile=sys.argv[1]
RN=int(sys.argv[2]) ## Estimated reads number
barcodecsv=sys.argv[3]
##Define hamming_distance
def hamming_distance(s1, s2):
    return sum(ch1 != ch2 for ch1,ch2 in zip(s1,s2))

##Make indexies dictionary from indexbarcode to a name
IndexDic={}
with open(barcodecsv) as csvf:
    c=csv.reader(csvf,delimiter='\t')
    for row in c:
        print(row)
        IndexDic[row[1].encode()]=row[0]



# fq = pyfastx.Fastx('210706_Multiome_Donor01_BMMC_ATAC_Nova_R1.fastq.gz')

## Create the demultiplexed files gz handle
Prefix=set(IndexDic.values())
WriteFileIn={}
for p in Prefix:
    WriteFileIn[p]=gzip.open(p+"_"+fqFile,'wb')

## Start to read the fq file
Read=0
bar = Bar('ReadingFq', max=RN)
with gzip.open(fqFile) as f:
    while True:
        Name=f.readline()
        if Name==b"":
            break
        Seq=f.readline()
        Opt=f.readline()
        Qua=f.readline()
        idx=Name.split()[1].split(b":")[3]
        for i in IndexDic.keys():
            if(hamming_distance(i,idx)<=1):
                WriteFileIn[IndexDic[i]].write(Name)
                WriteFileIn[IndexDic[i]].write(Seq)
                WriteFileIn[IndexDic[i]].write(Opt)
                WriteFileIn[IndexDic[i]].write(Qua)
                break
        # Read+=1
        print(Read)
        bar.next()

bar.finish()

##Close the gz files
for file in WriteFileIn.values():
    file.close()
