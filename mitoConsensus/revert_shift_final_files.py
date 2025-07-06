#!/usr/bin/env python3

import argparse
import csv
import glob
import os
import sys
import pandas as pd
from collections import defaultdict


# ─── CONFIGURATION: define exactly what columns you expect ────────────────────

# For merging RawGenotypes:
RAW_PREFIX_COLS = [
    "MoleculeID", "CellBC", "Position", "Variant", "Call", "Ref"
]
RAW_SUM_COLS = [
    "FamSize", "GT_Cts", "CSS", "DB_Cts", "SG_Cts",
    "ForwardStrand", "ReverseStrand"
]

# For merging QualifiedTotalCts:
QUAL_KEY_COLS = ["CellBC", "Position"]
# automatically sum every other column present in the header,
# but drop any stray columns 

# ──────────────────────────────────────────────────────────────────────────────


### set base dir to find files in 
BASE_DIR = "./final"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Shift ALT files and merge RawGenotypes & QualifiedTotalCts")
    parser.add_argument(
        "-f", "--ref-fasta", required=True,
        help="Reference FASTA – use length of its first record as ref length")
    parser.add_argument(
        "-b", "--barcodes", required=True,
        help="Path to your barcodes file (one barcode per line)")
    parser.add_argument(
        "-o", "--offset", type=int, default=8000,
        help="Coordinate shift offset (default: 8000)")
    
    return parser.parse_args()

def read_reference_length(fasta_path):
    sequence_lengths = []
    current_name = None
    current_length = 0

    try:
        with open(fasta_path) as fasta:
            for line in fasta:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_name is not None:
                        sequence_lengths.append((current_name, current_length))
                    current_name = line[1:].split()[0]
                    current_length = 0
                else:
                    current_length += len(line)
            if current_name is not None:
                sequence_lengths.append((current_name, current_length))
    except FileNotFoundError:
        sys.exit(f"ERROR: Cannot open FASTA '{fasta_path}'")

    if not sequence_lengths:
        sys.exit(f"ERROR: No sequences found in FASTA '{fasta_path}'")

    if len(sequence_lengths) > 1:
        name, length = sequence_lengths[0]
        sys.stderr.write(
            f"WARNING: {len(sequence_lengths)} sequences in FASTA; "
            f"using first: {name} (length={length})\n"
        )
        return length

    return sequence_lengths[0][1]

def shift_coordinate(position, offset, ref_length):
    return ((position - 1 + offset) % ref_length) + 1

def shift_raw_genotypes(input_file, offset, ref_length):
    output_file = input_file + ".shifted"
    with open(input_file, newline="") as in_handle, \
         open(output_file, "w", newline="") as out_handle:

        reader = csv.DictReader(in_handle, delimiter="\t")
        # just copy all columns through – merging will pick only the whitelist
        writer = csv.DictWriter(
            out_handle,
            fieldnames=reader.fieldnames,
            delimiter="\t",
            lineterminator="\n"
        )
        writer.writeheader()

        for record in reader:
            # MoleculeID format: "<sequence>_<start>_<end>"
            try:
                seq_name, start_str, end_str = record["MoleculeID"].rsplit("_", 2)
                start_pos = int(start_str)
                end_pos   = int(end_str)
            except Exception:
                sys.stderr.write(f"SKIP invalid MoleculeID: {record.get('MoleculeID')}\n")
                continue

            record["MoleculeID"] = (
                f"{seq_name}_"
                f"{shift_coordinate(start_pos, offset, ref_length)}_"
                f"{shift_coordinate(end_pos,   offset, ref_length)}"
            )

            # Position
            orig_pos = int(record["Position"])
            record["Position"] = str(shift_coordinate(orig_pos, offset, ref_length))

            # Variant: "<pos>_<ref>_<alt>"
            parts = record["Variant"].split("_")
            if len(parts) == 3 and parts[0].isdigit():
                new_pos = shift_coordinate(int(parts[0]), offset, ref_length)
                record["Variant"] = f"{new_pos}_{parts[1]}_{parts[2]}"
            else:
                sys.stderr.write(
                    f"WARNING: Could not parse Variant='{record['Variant']}' in {input_file}\n"
                )

            writer.writerow(record)

    print(f"Shifted RAW    : {input_file} → {output_file}")

def shift_qualified_counts(input_file, offset, ref_length):
    output_file = input_file + ".shifted"
    with open(input_file) as in_handle, \
         open(output_file, "w") as out_handle:

        header_line = in_handle.readline().rstrip("\n")
        out_handle.write(header_line + "\n")
        columns = header_line.split("\t")

        if "CellBC" not in columns or "Position" not in columns:
            sys.stderr.write(f"ERROR: Missing 'CellBC' or 'Position' in {input_file}\n")
            return

        pos_index = columns.index("Position")

        for line in in_handle:
            line = line.rstrip("\n")
            if not line:
                out_handle.write("\n")
                continue
            fields = line.split("\t")
            orig_pos = int(fields[pos_index])
            fields[pos_index] = str(shift_coordinate(orig_pos, offset, ref_length))
            out_handle.write("\t".join(fields) + "\n")

    print(f"Shifted QUAL   : {input_file} → {output_file}")



def merge_raw_genotypes(original_file, shifted_file):
    """
    Merge RawGenotypes:
      - Load each file (with or without header)
      - Concatenate both
      - Group by MoleculeID
      - Sum numeric fields for duplicates
      - Write merged output to <original_file>.combined
    """
    fields = RAW_PREFIX_COLS + RAW_SUM_COLS

    def load_raw_file(path, expected_fields):
        try:
            df = pd.read_csv(path, sep="\t", usecols=expected_fields, dtype=str)
        except ValueError:
            # Fallback: assume no header
            df = pd.read_csv(path, sep="\t", header=None, dtype=str)
            if df.shape[1] < len(expected_fields):
                sys.exit(f"ERROR: {path} has only {df.shape[1]} columns, expected at least {len(expected_fields)}")
            df = df.iloc[:, :len(expected_fields)]
            df.columns = expected_fields
        return df
    
   

   

    df1 = load_raw_file(original_file, fields)
    df2 = load_raw_file(shifted_file, fields)

    print(f"[DEBUG] {original_file}: {len(df1)} rows, {df1['MoleculeID'].nunique()} unique MoleculeIDs")
    print(f"[DEBUG] {shifted_file} : {len(df2)} rows, {df2['MoleculeID'].nunique()} unique MoleculeIDs")

    df = pd.concat([df1, df2], ignore_index=True)

    dup_mask = df.duplicated(subset=["MoleculeID", "Variant"], keep=False)
    dups = df[dup_mask].copy()
    unique = df[~dup_mask].copy()

    print(f"[DEBUG] Combined total rows: {len(df)}")
    print(f"[DEBUG] Unique MoleculeIDs in combined: {len(unique)}")
    print(f"[DEBUG] Duplicated MoleculeIDs in combined: {dups['MoleculeID'].nunique()}")

    for col in RAW_SUM_COLS:
        dups[col] = dups[col].astype(float)

    summed = dups.groupby(["MoleculeID", "Variant"], as_index=False).agg({
    **{col: "first" for col in RAW_PREFIX_COLS if col not in ["MoleculeID", "Variant"]},
    **{col: "sum" for col in RAW_SUM_COLS}
    })


    result = pd.concat([unique, summed], ignore_index=True)
    result = result[fields]
    result[RAW_SUM_COLS] = result[RAW_SUM_COLS].astype(float).astype(str)
    result.sort_values(by=["CellBC", "Position"], inplace=True)

    output_file = original_file + ".combined"
    result.to_csv(output_file, sep="\t", index=False)
    print(f"[DEBUG] Merge complete: {output_file}")



def merge_qualified_counts(original_file, shifted_file):
    """
    Simple merge: sum counts for matching (CellBC, Position) pairs
    """
    # Load both files
    df1 = pd.read_csv(original_file, sep='\t')
    df2 = pd.read_csv(shifted_file, sep='\t')
    
    # Merge on CellBC and Position, summing the count columns
    merged = pd.merge(df1, df2, on=['CellBC', 'Position'], how='outer', suffixes=('', '_y'))
    
    # Sum the count columns
    for col in ['Total', 'VerySensitive', 'Sensitive', 'Specific']:
        merged[col] = merged[col].fillna(0) + merged[f'{col}_y'].fillna(0)
        merged.drop(f'{col}_y', axis=1, inplace=True)
    
    # Save result
    combined_file = original_file + ".combined"
    merged.to_csv(combined_file, sep='\t', index=False)
    
    print(f"Merged: {original_file} + {shifted_file} → {combined_file}")


def filter_genotypes_from_combined(combined_file, ref_length,barcode_list):
    """
    Read the combined genotype file and create filtered versions based on
    consensus sequence support (CSS) and family size thresholds.
    
    Parameters:
    - combined_file: Path to the .RawGenotypes.Total.combined file
    """
    
    if not os.path.exists(combined_file):
        sys.exit(f"Error: Combined file {combined_file} does not exist.")
    
    # Generate output filenames - keep .combined and replace .Total with filter type
    out_very_sensitive = combined_file.replace('.Total.combined', '.VerySensitive.combined')
    out_sensitive = combined_file.replace('.Total.combined', '.Sensitive.combined')
    out_specific = combined_file.replace('.Total.combined', '.Specific.combined')
    
    print(f"Creating filtered files from: {combined_file}")
    print(f"  - Very Sensitive: {out_very_sensitive}")
    print(f"  - Sensitive: {out_sensitive}")
    print(f"  - Specific: {out_specific}")
    
  

    # Read combined file and apply filters
    with open(combined_file, 'r') as infile, \
         open(out_very_sensitive, 'w') as out_vs, \
         open(out_sensitive, 'w') as out_s, \
         open(out_specific, 'w') as out_sp:
        
        # Write hardcoded headers to all output files
        header = "MoleculeID\tCellBC\tPosition\tVariant\tCall\tRef\tFamSize\tGT_Cts\tCSS\tDB_Cts\tSG_Cts\tForwardStrand\tReverseStrand\n"
        out_vs.write(header)
        out_s.write(header)
        out_sp.write(header)
        
        # Skip the header line in input file
        infile.readline()
        
        # Process each line
        line_count = 0
        vs_count = 0
        s_count = 0
        sp_count = 0
        
        for line in infile:
            line_count += 1
            
            # Parse the line
            fields = line.strip().split('\t')
            if len(fields) < 13:
                print(f"Warning: Line {line_count} has fewer than expected columns, skipping.")
                continue
            
            try:
                # Extract relevant fields (0-indexed)
                fam_size = int(float(fields[6]))  # Convert to int, handling potential float strings
                css = float(fields[8])
                db_cts = int(float(fields[9]))    # Convert to int, handling potential float strings
                



                # Apply filtering criteria based on original logic
                if db_cts == 0:  # Single strand
                    # Very Sensitive: CSS > 0.75 and FamSize >= 2
                    if css > 0.75 and fam_size >= 2:
                        out_vs.write(line)
                        vs_count += 1
                    
                    # Sensitive: CSS > 0.75 and FamSize >= 3
                    if css > 0.75 and fam_size >= 3:
                        out_s.write(line)
                        s_count += 1
                    
                    # Specific: CSS > 0.9 and FamSize >= 4
                    if css > 0.9 and fam_size >= 4:
                        out_sp.write(line)
                        sp_count += 1
                        
                else:  # Double strand
                    # Very Sensitive: CSS > 0.75 and FamSize >= 1
                    if css > 0.75 and fam_size >= 1:
                        out_vs.write(line)
                        vs_count += 1
                    
                    # Sensitive: CSS > 0.75 and FamSize >= 2
                    if css > 0.75 and fam_size >= 2:
                        out_s.write(line)
                        s_count += 1
                    
                    # Specific: CSS > 0.9 and FamSize >= 3
                    if css > 0.9 and fam_size >= 3:
                        out_sp.write(line)
                        sp_count += 1
                        
            except (ValueError, IndexError) as e:
                print(f"Warning: Error parsing line {line_count}: {e}")
                continue

    print(f"\nFiltering complete:")
    print(f"  - Total variants processed: {line_count}")
    print(f"  - Very Sensitive variants: {vs_count}")
    print(f"  - Sensitive variants: {s_count}")
    print(f"  - Specific variants: {sp_count}")


def main():
    args       = parse_args()
    ref_length = read_reference_length(args.ref_fasta)
    barcodes_path = args.barcodes
    # Read barcodes file into a list
    with open(barcodes_path, 'r') as f:
        barcode_list = [line.strip() for line in f if line.strip()]
    print(f"Reference length = {ref_length}")

    # Gather ALT files

    
    alt_files = glob.glob(os.path.join(BASE_DIR,"ALT.RawGenotypes.Total*")) \
               + glob.glob(os.path.join(BASE_DIR,"ALT.QualifiedTotalCts*"))
    # Exclude already‐shifted
    alt_files = [f for f in alt_files if not f.endswith(".shifted")]
    if not alt_files:
        sys.exit("ERROR: No ALT files found to shift.")

    # 1) SHIFT step
    for fname in sorted(alt_files):
        if "QualifiedTotalCts" in fname:
            print("shifting qualified counts")
            shift_qualified_counts(fname, args.offset, ref_length)
        else:
            print(f"shifting: {fname}")
            shift_raw_genotypes(fname, args.offset, ref_length)

    # 2) MERGE RawGenotypes
    for shifted in glob.glob(os.path.join(BASE_DIR,"ALT.RawGenotypes.Total.shifted")):
        base     = os.path.basename(shifted).split("ALT.", 1)[1]
        original = os.path.join(BASE_DIR, base.replace(".shifted",""))
        if os.path.exists(original):
            merge_raw_genotypes(original, shifted)
        else:
            sys.stderr.write(f"WARNING: RawGenotypes original '{original}' not found\n")

    # 3) MERGE QualifiedTotalCts
    for shifted in glob.glob(os.path.join(BASE_DIR,"ALT.QualifiedTotalCts*.shifted")):
        base = os.path.basename(shifted).split("ALT.", 1)[1]
        original = os.path.join(BASE_DIR,base.replace(".shifted", ""))
        if os.path.exists(original):
            merge_qualified_counts(original, shifted)
        else:
            sys.stderr.write(f"WARNING: QualifiedTotalCts original '{original}' not found\n")
    # 4) CREATE FILTERED VERSIONS from combined Total files and recreate qualifiedtotalcts combined file
    print("\n" + "="*50)
    print("CREATING FILTERED GENOTYPE FILES")
    print("="*50)
    combined_file = os.path.join(BASE_DIR, "RawGenotypes.Total.combined")
    if os.path.exists(combined_file):
        filter_genotypes_from_combined(combined_file,ref_length,barcode_list)
    else:
        print(f"Warning: {combined_file} not found")

if __name__ == "__main__":
    main()