#!/usr/bin/env python3
"""
Filter a protein FASTA to keep only the longest isoform per gene.

Usage:
  python3 filter_longest_isoform.py <input.faa> <output.faa>

The output header of each retained sequence is:
  ><original header> gene_id=<gene_id> len=<sequence length>
"""
import re
import sys
from collections import defaultdict


def read_fasta(path):
    """Read FASTA -> {header: sequence}"""
    seqs = {}
    header = None
    parts = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(parts)
                header = line[1:].strip()
                parts = []
            else:
                parts.append(line.strip())
        if header:
            seqs[header] = "".join(parts)
    return seqs


def get_gene_id(header):
    """Extract a gene_id from the FASTA header (several formats supported)"""
    h = header.strip()

    # Format 1: gene_id "XXX" (NCBI GTF style)
    m = re.search(r'gene_id\s+"([^"]+)"', h)
    if m:
        return m.group(1)

    # Format 2: gene_id=XXX, gene=XXX, gene XXX, Parent=XXX
    for pat in [r'gene_id[=:\s]+(\S+)', r'gene[=:\s]+(\S+)', r'Parent=([^;\s]+)']:
        m = re.search(pat, h, re.IGNORECASE)
        if m:
            return m.group(1).rstrip(";")

    # Format 3: NCBI accession only (LOC / XM_ / XP_), no explicit gene_id.
    # In this case the transcript_id in the header is used as the gene key.
    return h.split()[0] if " " in h else h


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <input.faa> <output.faa>", file=sys.stderr)
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2]

    print(f"Reading: {infile}", file=sys.stderr)
    proteins = read_fasta(infile)
    total = len(proteins)
    print(f"  Total sequences: {total}", file=sys.stderr)

    # Group by gene_id and keep the longest sequence
    gene_best = {}  # gene_id -> (header, length, sequence)
    for header, seq in proteins.items():
        gid = get_gene_id(header)
        seq_len = len(seq)
        if gid not in gene_best or seq_len > gene_best[gid][1]:
            gene_best[gid] = (header, seq_len, seq)

    # Write output
    with open(outfile, "w") as out:
        for gid, (hdr, slen, seq) in sorted(gene_best.items()):
            out.write(f">{hdr} gene_id={gid} len={slen}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")

    removed = total - len(gene_best)
    print(f"  Genes retained: {len(gene_best)}", file=sys.stderr)
    print(f"  Isoforms removed: {removed}", file=sys.stderr)
    print(f"Written: {outfile}", file=sys.stderr)

    # Sample check of header format
    print("\nSample headers:", file=sys.stderr)
    for i, (gid, (hdr, slen, seq)) in enumerate(sorted(gene_best.items())):
        if i >= 3:
            break
        print(f"  Gene: {gid} | Len: {slen} | Header: {hdr[:120]}", file=sys.stderr)


if __name__ == "__main__":
    main()
