#!/usr/bin/env python3
"""
Extract CDS from a GTF annotation + genome FASTA and translate to proteins.

Pure Python: no gffread / AGAT dependency.

Usage:
  python3 extract_proteins.py --gtf annotation.gtf --genome genome.fasta \
      --out-cds sparrow_cds.fna --out-prot sparrow_prot_raw.faa

Defaults (if an argument is omitted):
  --gtf       <script dir>/work/input/annotation.gtf
  --genome    <script dir>/work/input/genome.fasta
  --out-cds   <script dir>/work/sparrow_cds.fna
  --out-prot  <script dir>/work/sparrow_prot_raw.faa
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path

# ---- Default paths (relative to this script) ----
BASE = Path(__file__).resolve().parent / "work"

# ---- Standard genetic code ----
CODON_TABLE = {
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L", "TTA": "L", "TTG": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TTC": "F", "TTT": "F", "TAC": "Y", "TAT": "Y",
    "TGC": "C", "TGT": "C", "TGA": "*", "TAA": "*", "TAG": "*",
    "TGG": "W", "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
}

COMPLEMENT = str.maketrans("ATGCatgc", "TACGtacg")


def parse_gtf_cds(gtf_file):
    """Parse the GTF and return {transcript_id: [(chrom, start, end, strand), ...]}"""
    print("[1/5] Parsing GTF for CDS features...", file=sys.stderr)
    cds_dict = defaultdict(list)
    total_lines = 0
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue
            total_lines += 1
            chrom, _, _, start, end, _, strand, _, attrs = cols

            # Extract transcript_id
            tid = None
            for token in attrs.split(";"):
                token = token.strip()
                if not token:
                    continue
                # GTF style: key "value" or key value
                if token.lower().startswith("transcript_id"):
                    if " " in token:
                        tid = token.split(" ", 1)[1].strip().strip('"')
                    elif "=" in token:
                        tid = token.split("=", 1)[1].strip().strip('"')
                    break

            if tid:
                cds_dict[tid].append((chrom, int(start) - 1, int(end), strand))

    print(f"   Found {total_lines} CDS lines across {len(cds_dict)} transcripts", file=sys.stderr)
    return cds_dict


def load_genome(genome_file):
    """Load the genome FASTA -> {contig_name: sequence}"""
    print("[2/5] Loading genome...", file=sys.stderr)
    genome = {}
    current_chrom = None
    current_seq = []

    with open(genome_file) as f:
        for line in f:
            if line.startswith(">"):
                if current_chrom is not None:
                    genome[current_chrom] = "".join(current_seq)
                current_chrom = line[1:].split()[0].strip()
                current_seq = []
            else:
                current_seq.append(line.strip())

    if current_chrom is not None:
        genome[current_chrom] = "".join(current_seq)

    total_bp = sum(len(v) for v in genome.values())
    print(f"   Loaded {len(genome)} contigs, {total_bp:,} bp", file=sys.stderr)
    return genome


def extract_cds(cds_dict, genome, out_cds):
    """Concatenate CDS exons and write the nucleotide sequences"""
    print("[3/5] Extracting CDS sequences...", file=sys.stderr)
    cds_count = 0

    with open(out_cds, "w") as out:
        for tid, exons in cds_dict.items():
            # Sort exons by strand direction
            strand = exons[0][3]
            if strand == "+":
                exons.sort(key=lambda x: x[1])
            else:
                exons.sort(key=lambda x: x[1], reverse=True)

            seq_parts = []
            for chrom, start, end, s in exons:
                if start >= end:
                    continue
                # Match the contig name
                chrom_match = chrom
                if chrom_match not in genome:
                    # Some GTFs use chr1 while the genome uses 1, etc.
                    for alt in [chrom, chrom.split(".")[0],
                                chrom.replace("chr", ""),
                                "chr" + chrom]:
                        if alt in genome:
                            chrom_match = alt
                            break
                    else:
                        continue

                if end <= len(genome[chrom_match]):
                    chunk = genome[chrom_match][start:end]
                    if s == "-":
                        chunk = chunk.translate(COMPLEMENT)[::-1]
                    seq_parts.append(chunk)

            seq = "".join(seq_parts)
            if len(seq) >= 3:
                out.write(f">{tid}\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i+60] + "\n")
                cds_count += 1

    print(f"   Wrote {cds_count} CDS sequences", file=sys.stderr)
    return cds_count


def translate_cds(in_cds, out_prot):
    """Translate CDS nucleotides -> protein sequences"""
    print("[4/5] Translating CDS -> protein...", file=sys.stderr)
    prot_count = 0

    header = None
    cds_seq = ""

    with open(in_cds) as fin, open(out_prot, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                # Translate the previous record
                if header and len(cds_seq) >= 3:
                    prot = _translate(cds_seq)
                    _write_fasta(fout, header, prot)
                    prot_count += 1
                header = line[1:].strip()
                cds_seq = ""
            else:
                cds_seq += line.strip()

        # Last record
        if header and len(cds_seq) >= 3:
            prot = _translate(cds_seq)
            _write_fasta(fout, header, prot)
            prot_count += 1

    print(f"   Translated {prot_count} proteins", file=sys.stderr)
    return prot_count


def _translate(dna_seq):
    """Translate DNA -> protein in codon triplets"""
    dna = dna_seq.upper()
    aa = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3]
        aa.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa)


def _write_fasta(f, header, seq):
    f.write(f">{header}\n")
    for i in range(0, len(seq), 60):
        f.write(seq[i:i+60] + "\n")


def show_sample(out_prot, n=3):
    """Print the headers of the first n proteins"""
    print(f"\n[5/5] Sample headers:", file=sys.stderr)
    count = 0
    with open(out_prot) as f:
        for line in f:
            if line.startswith(">"):
                print(f"   {line[1:].rstrip()}", file=sys.stderr)
                count += 1
                if count >= n:
                    break


def main():
    parser = argparse.ArgumentParser(
        description="Extract CDS from GTF + genome and translate to proteins.")
    parser.add_argument("--gtf", default=str(BASE / "input" / "annotation.gtf"),
                        help="input GTF annotation")
    parser.add_argument("--genome", default=str(BASE / "input" / "genome.fasta"),
                        help="input genome FASTA")
    parser.add_argument("--out-cds", default=str(BASE / "sparrow_cds.fna"),
                        help="output CDS FASTA")
    parser.add_argument("--out-prot", default=str(BASE / "sparrow_prot_raw.faa"),
                        help="output protein FASTA")
    args = parser.parse_args()

    # 1. Parse the GTF
    cds_dict = parse_gtf_cds(args.gtf)
    if not cds_dict:
        print("ERROR: No CDS features found in GTF!", file=sys.stderr)
        sys.exit(1)

    # 2. Load the genome
    genome = load_genome(args.genome)

    # 3. Extract CDS
    n_cds = extract_cds(cds_dict, genome, args.out_cds)
    if n_cds == 0:
        print("ERROR: No CDS extracted! Check contig name matching.", file=sys.stderr)
        sys.exit(1)

    # 4. Translate
    translate_cds(args.out_cds, args.out_prot)

    # 5. Preview
    show_sample(args.out_prot)

    print(f"\nDone: CDS={args.out_cds} | Protein={args.out_prot}", file=sys.stderr)


if __name__ == "__main__":
    main()
