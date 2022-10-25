import os
import sys

path=sys.argv[1]
prefix=sys.argv[2]

files=os.listdir(path)
MergedNW=""
for fname in files:
  print(path+fname)
  with open(path+fname) as f:
    NW=f.readline().strip()
    MergedNW=MergedNW+NW+"\n"
with open(path+prefix+".nw","w") as f:
    f.writelines(MergedNW)
