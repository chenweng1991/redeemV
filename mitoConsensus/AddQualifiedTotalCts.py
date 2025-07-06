import sys
out_dir=sys.argv[1]

# Initialize dictionaries for both regular and combined files
Total={}
VerySensitive={}
Sensitive={}
Specific={}

Total_combined={}
VerySensitive_combined={}
Sensitive_combined={}
Specific_combined={}

def process_file(filename, target_dict):
    """Helper function to process a file and populate the target dictionary"""
    try:
        with open(filename) as f:
            for line in f:
                line = line.strip() 
                if not line: 
                    continue
                content = line.split()
                CellPos = content[1] + content[2]
                if CellPos in target_dict.keys():
                    target_dict[CellPos].append(line.strip())
                else:
                    target_dict[CellPos] = [line.strip()]
    except FileNotFoundError:
        print(f"Warning: {filename} not found, skipping...")

## Read in the genotypes
print("Start...")

# Process regular files
process_file(out_dir+"/final/RawGenotypes.Total", Total)
print("Dic1 In")

process_file(out_dir+"/final/RawGenotypes.VerySensitive", VerySensitive)
print("Dic2 In")

process_file(out_dir+"/final/RawGenotypes.Sensitive", Sensitive)
print("Dic3 In")

process_file(out_dir+"/final/RawGenotypes.Specific", Specific)
print("Dic4 In")

# Process combined files
process_file(out_dir+"/final/RawGenotypes.Total.combined", Total_combined)
print("Dic1 Combined In")

process_file(out_dir+"/final/RawGenotypes.VerySensitive.combined", VerySensitive_combined)
print("Dic2 Combined In")

process_file(out_dir+"/final/RawGenotypes.Sensitive.combined", Sensitive_combined)
print("Dic3 Combined In")

process_file(out_dir+"/final/RawGenotypes.Specific.combined", Specific_combined)
print("Dic4 Combined In")

# Open output files for both regular and combined (overwriting input files like original)
o_Total = open(out_dir+"/final/RawGenotypes.Total","w")
o_VerySensitive = open(out_dir+"/final/RawGenotypes.VerySensitive","w")
o_Sensitive = open(out_dir+"/final/RawGenotypes.Sensitive","w")
o_Specific = open(out_dir+"/final/RawGenotypes.Specific","w")

o_Total_combined = open(out_dir+"/final/RawGenotypes.Total.combined","w")
o_VerySensitive_combined = open(out_dir+"/final/RawGenotypes.VerySensitive.combined","w")
o_Sensitive_combined = open(out_dir+"/final/RawGenotypes.Sensitive.combined","w")
o_Specific_combined = open(out_dir+"/final/RawGenotypes.Specific.combined","w")


def write_genotypes(cell_pos, content, dicts, output_files):
    """Helper function to write genotypes for a given cell position"""
    total_dict, vs_dict, s_dict, sp_dict = dicts
    o_total, o_vs, o_s, o_sp = output_files
    
    if cell_pos in total_dict.keys():
        for molecule in total_dict[cell_pos]:
            o_total.write(molecule+"\t"+content[2]+"\n")
        if cell_pos in vs_dict.keys():
            for molecule in vs_dict[cell_pos]:
                o_vs.write(molecule+"\t"+content[3]+"\n")
            if cell_pos in s_dict.keys():
                for molecule in s_dict[cell_pos]:
                    o_s.write(molecule+"\t"+content[4]+"\n")
                if cell_pos in sp_dict.keys():
                    for molecule in sp_dict[cell_pos]:
                        o_sp.write(molecule+"\t"+content[5]+"\n")

with open(out_dir+"/final/QualifiedTotalCts.combined") as f:
    for line in f:
        line = line.strip()
        if not line: 
            continue
        content = line.split()
        CellPos = content[0] + content[1]
        
        # Process regular dictionaries
        write_genotypes(CellPos, content, 
                       (Total, VerySensitive, Sensitive, Specific),
                       (o_Total, o_VerySensitive, o_Sensitive, o_Specific))
        
        # Process combined dictionaries
        write_genotypes(CellPos, content,
                       (Total_combined, VerySensitive_combined, Sensitive_combined, Specific_combined),
                       (o_Total_combined, o_VerySensitive_combined, o_Sensitive_combined, o_Specific_combined))

# Close all output files
o_Total.close()
o_VerySensitive.close()
o_Sensitive.close()
o_Specific.close()

o_Total_combined.close()
o_VerySensitive_combined.close()
o_Sensitive_combined.close()
o_Specific_combined.close()