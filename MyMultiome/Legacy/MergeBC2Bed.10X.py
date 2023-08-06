import sys


##Make Barcode dictionary
BCfilename=sys.argv[1]
Bedfilename=sys.argv[2]
BCcelldic={}
# i=0
with open(BCfilename) as f:
    while True:
        # print(i)
        # i=i+1
        Line=f.readline().strip()
        if Line=="":
            break
        else:
            LineArray=Line.split()
            BCcelldic[LineArray[0][1:]]=LineArray[1]

##
## Read in the Bed fikle and search the reads name for Cell Barcode
with open(Bedfilename) as f:
    while True:
        Line=f.readline().strip()
        if Line=="":
            break
        else:
            Record=Line.split()
            if(Record[3] in BCcelldic.keys()):
                #print(BCcelldic[Record[3]]+'\t'+Line)
                sys.stdout.write(BCcelldic[Record[3]]+'\t'+Line+'\n')
