import sys

out_dir = sys.argv[1]

StrandBiaseBlackList_file = out_dir + "/final/StrandBiaseBlackList"
Total = out_dir + "/final/RawGenotypes.Total"
VerySensitive = out_dir + "/final/RawGenotypes.VerySensitive"
Sensitive = out_dir + "/final/RawGenotypes.Sensitive"
Specific = out_dir + "/final/RawGenotypes.Specific"

Total_o_file = Total + ".StrandBalance"
VerySensitive_o_file = VerySensitive + ".StrandBalance"
Sensitive_o_file = Sensitive + ".StrandBalance"
Specific_o_file = Specific + ".StrandBalance"

# Read strand bias blacklist into a set (faster lookup than a list/tuple)
StrandBiasDic = set()

with open(StrandBiaseBlackList_file) as f:
    for line in f:
        StrandBiasDic.add(line.strip())

# Open output files
Total_o = open(Total_o_file, "w")
VerySensitive_o = open(VerySensitive_o_file, "w")
Sensitive_o = open(Sensitive_o_file, "w")
Specific_o = open(Specific_o_file, "w")

def process_file(input_file, output_file):
    with open(input_file) as f:
        header = f.readline()  # Skip header line
        for line in f:
            content = line.strip().split()
            if len(content) <= 9:  # Ensure there are enough columns
                continue
            V = content[3]
            try:
                DB = float(content[9])
            except ValueError:
                continue  # Skip malformed lines
            if V in StrandBiasDic and DB == 0:
                continue
            output_file.write(line)

# Process each input file
process_file(Total, Total_o)
process_file(VerySensitive, VerySensitive_o)
process_file(Sensitive, Sensitive_o)
process_file(Specific, Specific_o)

# Close output files
Total_o.close()
VerySensitive_o.close()
Sensitive_o.close()
Specific_o.close()
