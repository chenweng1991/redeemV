'''
File: AddBC2Fastq.py
Author: Chen Weng
Date: 2023-07-30
Updated: 2023-07-30

Description:
    This script is designed to append barcodes from an i5 barcode fastq file to the names of sequences in paired-end FASTQ reads (R1 and R2).
    It accepts two input FASTQ files (or gzipped FASTQ files), an i5 barcode FASTQ file, and two output FASTQ files as command-line arguments.
    By extracting and transforming the barcodes, it creates two output files with the barcodes appended to the sequence identifiers.

Usage:
    python3 AddBC2Fastq.py <input1> <input2> <barcodes> <output1> <output2>
'''


import gzip
import argparse
import sys

CANDIDATE_RULES = (
    ("rc(8:24)", 8, 24, True),
    ("rc(10:26)", 10, 26, True),
    ("rc(6:22)", 6, 22, True),
)

def is_gzip_file(filepath):
    try:
        with gzip.open(filepath, 'rb') as f:
            f.read(1)
        return True
    except OSError:
        return False
    
def open_file(filepath, mode, buffering=-1):
    if is_gzip_file(filepath):
        return gzip.open(filepath, mode)
    else:
        return open(filepath, mode,buffering=buffering)

# def reverse(seq):
#     """Returns a reversed string"""
#     return seq[::-1]


# def complement(seq):
#     """Returns a complement DNA sequence"""
#     complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
#     seq_list = list(seq)
#     seq_list = [complement_dict[base] for base in seq_list]
#     return ''.join(seq_list)

# def reverse_complement(seq):
#     """"Returns a reverse complement DNA sequence"""
#     seq = reverse(seq)
#     seq = complement(seq)
#     return seq
TRANSLATION_TABLE = str.maketrans('ATCGN', 'TAGCN')
def reverse_complement(seq):
    return seq[::-1].translate(TRANSLATION_TABLE)

def extract_barcode(sequence, start, end, revcomp):
    barcode = sequence[start:end]
    if revcomp:
        barcode = reverse_complement(barcode)
    return barcode

def load_whitelist(path):
    whitelist = set()
    with open(path, 'rt') as handle:
        for line in handle:
            barcode = line.strip()
            if barcode:
                whitelist.add(barcode)
    return whitelist

def choose_rule(barcodes_path, whitelist_path, detect_reads, min_default_match_rate, min_selected_match_rate):
    whitelist = load_whitelist(whitelist_path)
    if not whitelist:
        raise ValueError(f"Whitelist is empty: {whitelist_path}")

    counts = {name: 0 for name, _, _, _ in CANDIDATE_RULES}
    matches = {name: 0 for name, _, _, _ in CANDIDATE_RULES}
    sampled = 0

    with open_file(barcodes_path, 'rt') as bc_file:
        while sampled < detect_reads:
            l1 = bc_file.readline()
            if l1 == "":
                break
            sequence = bc_file.readline().strip()
            bc_file.readline()
            bc_file.readline()

            sampled += 1
            for name, start, end, revcomp in CANDIDATE_RULES:
                barcode = extract_barcode(sequence, start, end, revcomp)
                if len(barcode) != (end - start) or "N" in barcode:
                    continue
                counts[name] += 1
                if barcode in whitelist:
                    matches[name] += 1

    default_name = CANDIDATE_RULES[0][0]
    default_rate = matches[default_name] / counts[default_name] if counts[default_name] else 0.0
    best_name = max(
        (name for name, _, _, _ in CANDIDATE_RULES),
        key=lambda name: (matches[name], counts[name], name == default_name),
    )
    rule_lookup = {name: (start, end, revcomp) for name, start, end, revcomp in CANDIDATE_RULES}
    if default_rate >= min_default_match_rate or best_name == default_name:
        return default_name, rule_lookup[default_name], sampled, counts, matches

    best_rate = matches[best_name] / counts[best_name] if counts[best_name] else 0.0
    if best_rate < min_selected_match_rate:
        raise ValueError(
            "No supported barcode extraction rule matched the whitelist well enough. "
            f"Default rc(8:24) matched {matches[default_name]}/{counts[default_name]} "
            f"({default_rate:.4f}); best candidate {best_name} matched "
            f"{matches[best_name]}/{counts[best_name]} ({best_rate:.4f}), below required "
            f"selected-rule threshold {min_selected_match_rate:.4f}."
        )
    print(
        "WARNING: default rc(8:24) whitelist match rate "
        f"{matches[default_name]}/{counts[default_name]} ({default_rate:.4f}) is below "
        f"threshold {min_default_match_rate:.4f}; switching to {best_name} with "
        f"{matches[best_name]}/{counts[best_name]} ({best_rate:.4f}) based on {sampled} barcode reads.",
        file=sys.stderr,
    )
    return best_name, rule_lookup[best_name], sampled, counts, matches


def write_fastq_file(in_file, out_file, barcode):
    read_line = next(in_file).strip().split()[0]
    out_file.write(read_line + '|' + barcode + '\n')
    out_file.write(next(in_file))  # sequence
    out_file.write(next(in_file))  # separator
    out_file.write(next(in_file))  # quality

def append_barcodes(input1, input2, barcodes, output1, output2, buffering, start=8, end=24, revcomp=True):
    with open_file(input1, 'rt',buffering = buffering) as in1, \
         open_file(input2, 'rt',buffering = buffering) as in2, \
         open_file(barcodes, 'rt',buffering = buffering) as bc_file, \
         open_file(output1, 'wt',buffering = buffering) as out1, \
         open_file(output2, 'wt',buffering = buffering) as out2:

        while True:
            # Read 4 lines from the barcode file
            L1 = bc_file.readline().strip()
            L2 = bc_file.readline().strip()
            L3 = bc_file.readline().strip()
            L4 = bc_file.readline().strip()

            if L1 == "":
                break

            Sequence = L2
            CellBC = extract_barcode(Sequence, start, end, revcomp)

            # Process input1 and input2
            write_fastq_file(in1, out1, CellBC)
            write_fastq_file(in2, out2, CellBC)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get Cell barcode from i5 and Append barcodes to FASTQ file names.')
    parser.add_argument('input1', help='Path to the R1 input FASTQ (or gzipped FASTQ) file.')
    parser.add_argument('input2', help='Path to the R2 input FASTQ (or gzipped FASTQ) file.')
    parser.add_argument('barcodes', help='Path to the i5 barcode FASTQ (or gzipped FASTQ) file.')
    parser.add_argument('output1', help='Path to the R1 output FASTQ file.')
    parser.add_argument('output2', help='Path to the R2 output FASTQ file.')
    parser.add_argument('--buffering', type=int, default=-1, help='Buffer size for reading/writing files (in bytes). Default is -1.')
    parser.add_argument('--whitelist', help='Optional whitelist of expected cell barcodes used to validate barcode slicing.')
    parser.add_argument('--detect-reads', type=int, default=50000, help='Number of barcode reads to sample when validating extraction against a whitelist.')
    parser.add_argument('--min-default-match-rate', type=float, default=0.05, help='Minimum whitelist match rate required to keep the default rc(8:24) rule.')
    parser.add_argument('--min-selected-match-rate', type=float, default=0.50, help='Minimum whitelist match rate required to switch to a non-default barcode extraction rule.')

    args = parser.parse_args()

    start = 8
    end = 24
    revcomp = True
    if args.whitelist:
        rule_name, rule, sampled, counts, matches = choose_rule(
            args.barcodes,
            args.whitelist,
            args.detect_reads,
            args.min_default_match_rate,
            args.min_selected_match_rate,
        )
        start, end, revcomp = rule
        print(
            f"Barcode extraction rule selected: {rule_name}; sampled={sampled}; "
            f"default_matches={matches['rc(8:24)']}/{counts['rc(8:24)']}",
            file=sys.stderr,
        )

    append_barcodes(
        args.input1,
        args.input2,
        args.barcodes,
        args.output1,
        args.output2,
        buffering=args.buffering,
        start=start,
        end=end,
        revcomp=revcomp,
    )
