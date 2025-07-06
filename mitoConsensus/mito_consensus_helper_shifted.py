import os
import sys
import pysam
import numpy as np
import pandas as pd
import argparse
import logging
from collections import defaultdict
from typing import Dict, Any, List
from mito_consensus_nf_updated_shift import generate_genotype_matrices

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)


def build_molecule_dict(bam_file: str, barcode_tag: str):
    """
    Build read_pair_dict: {readName → [read0, read1]}
    and molecule_dict: {moleculeID → [readName, ...]} based on ALT BAM.
    """
    local_read_pair_dict = defaultdict(list)
    local_molecule_dict = defaultdict(list)

    with pysam.AlignmentFile(bam_file, "rb") as bam_in:
        for read in bam_in:
            local_read_pair_dict[read.query_name].append(read)

    for read_name, reads in local_read_pair_dict.items():
        if len(reads) != 2:
            continue
        r0, r1 = reads
        if r0.is_reverse and not r1.is_reverse:
            fwd, rev = r1, r0
        elif not r0.is_reverse and r1.is_reverse:
            fwd, rev = r0, r1
        else:
            continue

        if fwd.has_tag(barcode_tag):
            cell_bc = fwd.get_tag(barcode_tag)
            molecule_id = f"{cell_bc}_{fwd.pos}_{fwd.pos + abs(fwd.tlen)}"
            local_molecule_dict[molecule_id].append(read_name)

    return local_read_pair_dict, local_molecule_dict


def generate_alt_matrices(
    molecule_dict: Dict[str, List[str]],
    read_pair_dict: Dict[str, List[Any]],
    barcodes: List[str],
    alt_ref_df: pd.DataFrame,
    baseq_thld_hi: int
):
    """
    For each molecule in ALT BAM:
    - allocate SG_Genotypes (max_bp × 5),
    - DB_Genotypes (max_bp × 5),
    - Strand_mtx (max_bp × 2),
    and fill them by iterating read pairs.
    """
    dna_letters = ["A", "C", "G", "T", "N"]
    max_bp = alt_ref_df.shape[0]
    stats = {
        "molecules_total": len(molecule_dict),
        "molecules_processed": 0
    }

    SG_matrices = {}
    DB_matrices = {}
    Strand_matrices = {}

    for idx, mol_id in enumerate(molecule_dict):
        if idx % 100 == 0:
            logger.info(f"Alt molecule {idx+1}/{stats['molecules_total']}")
        sg = np.zeros((max_bp, 5), dtype=int)
        db = np.zeros((max_bp, 5), dtype=int)
        strand = np.zeros((max_bp, 2), dtype=int)

        for read_name in molecule_dict[mol_id]:
            try:
                read0, read1 = read_pair_dict[read_name]
                seq0, seq1 = read0.seq, read1.seq
                qual0, qual1 = read0.query_qualities, read1.query_qualities

                pos0 = np.asarray(read0.get_aligned_pairs(matches_only=True))
                pos1 = np.asarray(read1.get_aligned_pairs(matches_only=True))
                if pos0.size == 0 or pos1.size == 0:
                    continue

                overlap_positions = np.intersect1d(pos0[:, 1], pos1[:, 1])
                spec0 = pos0[~np.isin(pos0[:, 1], overlap_positions)]
                spec1 = pos1[~np.isin(pos1[:, 1], overlap_positions)]
                ov0 = pos0[np.isin(pos0[:, 1], overlap_positions)]
                ov1 = pos1[np.isin(pos1[:, 1], overlap_positions)]

                # non-overlap left (read0)
                for base0 in spec0:
                    i0, refpos = base0
                    if qual0[i0] > baseq_thld_hi:
                        nt = seq0[i0]
                        sg[refpos, dna_letters.index(nt)] += 1
                        strand[refpos, int(read0.is_reverse)] += 1
                    else:
                        sg[refpos, 4] += 1

                # overlap
                for b0, b1 in zip(ov0, ov1):
                    i0, p = b0
                    i1, _ = b1
                    if seq0[i0] == seq1[i1]:
                        if qual0[i0] > baseq_thld_hi or qual1[i1] > baseq_thld_hi:
                            nt = seq0[i0]
                            db[p, dna_letters.index(nt)] += 1
                            strand[p, 0] += 1
                            strand[p, 1] += 1
                        else:
                            db[p, 4] += 1
                    else:
                        db[p, 4] += 1

                # non-overlap right (read1)
                for base1 in spec1:
                    i1, refpos = base1
                    if qual1[i1] > baseq_thld_hi:
                        nt = seq1[i1]
                        sg[refpos, dna_letters.index(nt)] += 1
                        strand[refpos, int(read1.is_reverse)] += 1
                    else:
                        sg[refpos, 4] += 1

            except Exception as e:
                logger.error(f"Error in molecule {mol_id}, read {read_name}: {e}")
                continue

        SG_matrices[mol_id] = sg
        DB_matrices[mol_id] = db
        Strand_matrices[mol_id] = strand
        stats["molecules_processed"] += 1

    logger.info(
        f"Alt matrices built: {stats['molecules_processed']}/{stats['molecules_total']} molecules"
    )
    return SG_matrices, DB_matrices, Strand_matrices



def run_alt_only(alt_bam: str, barcode_file: str, alt_chrM_ref: str, baseq_thld_hi: int):
    """
    Build SG/DB/Strand matrices from the ALT BAM, then call raw genotypes.

    Parameters:
    - alt_bam: Path to the ALT‐shifted BAM (indexed, with “BC” tags).
    - barcode_file: One‐barcode‐per‐line text file.
    - alt_chrM_ref: Two‐column TSV of “pos” and “base” for each mito position.
    - baseq_thld_hi: Base‐quality threshold (e.g. 30).
    """
    # 1) Read barcodes into a list
    with open(barcode_file, "r") as fh:
        barcodes = [line.strip() for line in fh]

    # 2) Build read‐pair and molecule dictionaries from ALT BAM
    read_pair_dict_alt, molecule_dict_alt = build_molecule_dict(alt_bam, "BC")

    # 3) Load ALT reference into a DataFrame
    alt_ref_df = pd.read_table(alt_chrM_ref, names=["pos", "base"])

    # 4) Generate SG/DB/Strand matrices
    sg_dict, db_dict, strand_dict = generate_alt_matrices(
        molecule_dict=molecule_dict_alt,
        read_pair_dict=read_pair_dict_alt,
        barcodes=barcodes,
        alt_ref_df=alt_ref_df,
        baseq_thld_hi=baseq_thld_hi
    )

    # 5) Write out SG/DB/Strand CSVs under “alt_matrices/”
    out_dir = "alt_matrices"
    os.makedirs(out_dir, exist_ok=True)
    for mol_id in sg_dict:
        sg_df = pd.DataFrame(sg_dict[mol_id], columns=["A", "C", "G", "T", "N"])
        sg_df.index += 1
        sg_df.to_csv(os.path.join(out_dir, f"{mol_id}_SG.csv"), index_label="Position")

        db_df = pd.DataFrame(db_dict[mol_id], columns=["A", "C", "G", "T", "N"])
        db_df.index += 1
        db_df.to_csv(os.path.join(out_dir, f"{mol_id}_DB.csv"), index_label="Position")

        strand_df = pd.DataFrame(strand_dict[mol_id], columns=["Forward", "Reverse"])
        strand_df.index += 1
        strand_df.to_csv(os.path.join(out_dir, f"{mol_id}_Strand.csv"), index_label="Position")

    # 6) Call raw genotype‐calling on ALT molecules
    raw_out = "alt_raw_genotypes"
    os.makedirs(raw_out, exist_ok=True)
    barcode_set = set(barcodes)

    generate_genotype_matrices(
        molecule_dict=molecule_dict_alt,
        read_pair_dict=read_pair_dict_alt,
        bcs=barcode_set,
        mito_ref=alt_ref_df,
        BaseQ_thld_hi=baseq_thld_hi,
        out_genotypeTotal_file=os.path.join(raw_out, "ALT.RawGenotypes.Total"),
        out_genotypeVerySensitive_file=os.path.join(raw_out, "ALT.RawGenotypes.VerySensitive"),
        out_genotypeSensitive_file=os.path.join(raw_out, "ALT.RawGenotypes.Sensitive"),
        out_genotypeSpecific_file=os.path.join(raw_out, "ALT.RawGenotypes.Specific"),
        out_totalCts_file=os.path.join(raw_out, "ALT.QualifiedTotalCts")
    )
    logger.info("Finished writing raw‐genotype tables for ALT molecules.")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build SG/DB/Strand matrices from an ALT-shifted BAM."
    )
    parser.add_argument(
        "--alt_bam", type=str, required=True, help="Path to the ALT-shifted BAM."
    )
    parser.add_argument(
        "--barcode", type=str, required=True, help="Path to the barcode file."
    )
    parser.add_argument(
        "--alt_chrM_ref", type=str, required=True, help="Path to the ALT chrM reference."
    )
    parser.add_argument(
        "--BaseQ_thld_hi",
        type=int,
        default=30,
        help="Base quality threshold (default: 30)."
    )
    args = parser.parse_args()

    run_alt_only(
        alt_bam=args.alt_bam,
        barcode_file=args.barcode,
        alt_chrM_ref=args.alt_chrM_ref,
        baseq_thld_hi=args.BaseQ_thld_hi
    )
