import sys

out_dir = sys.argv[1]

# Define file paths for regular files
StrandBiaseBlackList_file = out_dir + "/final/StrandBiaseBlackList"
Total = out_dir + "/final/RawGenotypes.Total"
VerySensitive = out_dir + "/final/RawGenotypes.VerySensitive"
Sensitive = out_dir + "/final/RawGenotypes.Sensitive"
Specific = out_dir + "/final/RawGenotypes.Specific"

Total_o_file = Total + ".StrandBalance"
VerySensitive_o_file = VerySensitive + ".StrandBalance"
Sensitive_o_file = Sensitive + ".StrandBalance"
Specific_o_file = Specific + ".StrandBalance"

# Define file paths for combined files
StrandBiaseBlackList_combined_file = out_dir + "/final/StrandBiaseBlackList.combined"
Total_combined = out_dir + "/final/RawGenotypes.Total.combined"
VerySensitive_combined = out_dir + "/final/RawGenotypes.VerySensitive.combined"
Sensitive_combined = out_dir + "/final/RawGenotypes.Sensitive.combined"
Specific_combined = out_dir + "/final/RawGenotypes.Specific.combined"

Total_combined_o_file = Total_combined + ".StrandBalance"
VerySensitive_combined_o_file = VerySensitive_combined + ".StrandBalance"
Sensitive_combined_o_file = Sensitive_combined + ".StrandBalance"
Specific_combined_o_file = Specific_combined + ".StrandBalance"

def load_blacklist(blacklist_file):
    """Load strand bias blacklist into a set for fast lookup"""
    blacklist = set()
    try:
        with open(blacklist_file) as f:
            for line in f:
                blacklist.add(line.strip())
    except FileNotFoundError:
        print(f"Warning: {blacklist_file} not found, proceeding without blacklist filtering...")
    return blacklist

def process_file(input_file, output_file, blacklist):
    """Process a single file with strand bias filtering"""
    try:
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
                if V in blacklist and DB == 0:
                    continue
                output_file.write(line)
    except FileNotFoundError:
        print(f"Warning: {input_file} not found, skipping...")

# Load blacklists
print("Loading regular blacklist...")
StrandBiasDic = load_blacklist(StrandBiaseBlackList_file)

print("Loading combined blacklist...")
StrandBiasDic_combined = load_blacklist(StrandBiaseBlackList_combined_file)

# Open output files for regular processing
Total_o = open(Total_o_file, "w")
VerySensitive_o = open(VerySensitive_o_file, "w")
Sensitive_o = open(Sensitive_o_file, "w")
Specific_o = open(Specific_o_file, "w")

# Open output files for combined processing
Total_combined_o = open(Total_combined_o_file, "w")
VerySensitive_combined_o = open(VerySensitive_combined_o_file, "w")
Sensitive_combined_o = open(Sensitive_combined_o_file, "w")
Specific_combined_o = open(Specific_combined_o_file, "w")

print("Processing regular files...")
# Process regular files
process_file(Total, Total_o, StrandBiasDic)
process_file(VerySensitive, VerySensitive_o, StrandBiasDic)
process_file(Sensitive, Sensitive_o, StrandBiasDic)
process_file(Specific, Specific_o, StrandBiasDic)

print("Processing combined files...")
# Process combined files
process_file(Total_combined, Total_combined_o, StrandBiasDic_combined)
process_file(VerySensitive_combined, VerySensitive_combined_o, StrandBiasDic_combined)
process_file(Sensitive_combined, Sensitive_combined_o, StrandBiasDic_combined)
process_file(Specific_combined, Specific_combined_o, StrandBiasDic_combined)

# Close all output files
print("Closing output files...")
Total_o.close()
VerySensitive_o.close()
Sensitive_o.close()
Specific_o.close()

Total_combined_o.close()
VerySensitive_combined_o.close()
Sensitive_combined_o.close()
Specific_combined_o.close()

print("Processing complete!")