#!/usr/bin/env python3
"""
Parse OrthoFinder Orthogroups.tsv -> gene-symbol-annotated ortholog maps.

Usage:
  python3 parse_orthogroups.py <Orthogroups.tsv> <human_proteome.faa> \
      <annotation.gtf> <output_dir>

Outputs (in <output_dir>):
  sparrow_human_one2one.tsv        Sparrow_Gene, Human_Symbol, Human_Ensembl, Orthogroup
  sparrow_human_all_orthologs.tsv  all ortholog pairs with a relation label
  ortholog_summary.txt             summary statistics
"""
import os
import re
import sys
from collections import defaultdict


def build_human_ensp_to_symbol(human_faa):
    """Build ENSP -> gene symbol from an Ensembl human protein FASTA"""
    print(f"[1/3] Building human ENSP->symbol map from: {human_faa}", file=sys.stderr)
    ensp2symbol = {}
    total = 0

    with open(human_faa) as f:
        for line in f:
            if line.startswith(">"):
                total += 1
                # >ENSP00000493376.1 pep ... gene:ENSG00000279973.2 ... gene_symbol:XXX ...
                ensp_m = re.match(r'>(\S+)', line)
                sym_m = re.search(r'gene_symbol:(\S+)', line)

                if ensp_m and sym_m:
                    ensp_id = ensp_m.group(1)              # e.g. ENSP00000493376.1
                    ensp_clean = re.sub(r'\.\d+$', '', ensp_id)  # ENSP00000493376
                    if ensp_clean not in ensp2symbol:
                        ensp2symbol[ensp_clean] = sym_m.group(1)

    print(f"  {total} proteins -> {len(ensp2symbol)} ENSP->symbol mappings", file=sys.stderr)
    return ensp2symbol


def build_sparrow_transcript_to_gene(gtf_file):
    """Build transcript_id -> gene_id from the GTF"""
    print("[2/3] Building sparrow transcript->gene map from GTF", file=sys.stderr)
    tx2gene = {}
    tx_count = 0

    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue

            attrs = cols[8]
            tid = None
            gid = None

            # Extract transcript_id and gene_id
            m = re.search(r'transcript_id\s+"([^"]+)"', attrs)
            if m:
                tid = m.group(1)
            m = re.search(r'gene_id\s+"([^"]+)"', attrs)
            if m:
                gid = m.group(1)

            if tid and gid and tid not in tx2gene:
                tx2gene[tid] = gid
                tx_count += 1

    gene_count = len(set(tx2gene.values()))
    print(f"  {tx_count} transcripts -> {gene_count} unique genes", file=sys.stderr)
    return tx2gene


def resolve_sparrow(seq_id, tx2gene):
    """Map a Sparrow sequence ID to a gene name.

    Priority: gene_id in the header, then transcript_id -> gene_id via the GTF,
    then the XP_ accession with the version suffix removed, then the raw ID.
    """
    # seq_id may look like:
    #   XM_039733698.1 gene_id=LOC120513861 len=xxx  (from filter_longest_isoform.py)
    #   XM_039733698.1                               (plain transcript ID)
    #   XP_039589632.1                               (protein ID)

    # Try the gene_id= field in the header first
    m = re.search(r'gene_id=(\S+)', seq_id)
    if m:
        return m.group(1).split()[0]  # drop a trailing len=xxx

    # Otherwise transcript_id -> gene_id (GTF mapping)
    clean = seq_id.split()[0]
    if clean in tx2gene:
        return tx2gene[clean]

    # XP_ protein ID: drop the version suffix and use it as the gene name
    if clean.startswith("XP_"):
        return re.sub(r'\.\d+$', '', clean)

    # Fallback
    return re.sub(r'\.\d+$', '', clean)


def resolve_human(seq_id, ensp2symbol):
    """Map an Ensembl human protein ID (ENSP) to a gene symbol"""
    # seq_id may look like:
    #   ENSP00000493376.1 gene_id=ENSG00000279973 len=xxx
    #   ENSP00000493376.1

    # If the header carries gene_id=, keep that value (may already be a symbol)
    m = re.search(r'gene_id=(\S+)', seq_id)
    if m:
        return m.group(1).split()[0]

    # ENSP -> symbol
    clean = seq_id.split()[0]
    clean = re.sub(r'\.\d+$', '', clean)

    if clean in ensp2symbol:
        return ensp2symbol[clean]

    # Fallback: return the raw ENSP accession
    return clean


def main():
    if len(sys.argv) < 5:
        print(f"Usage: {sys.argv[0]} <Orthogroups.tsv> <human_proteome.faa> <annotation.gtf> <output_dir>",
              file=sys.stderr)
        sys.exit(1)

    tsv_file = sys.argv[1]
    human_faa = sys.argv[2]
    gtf_file = sys.argv[3]
    outdir = sys.argv[4]

    os.makedirs(outdir, exist_ok=True)

    # ---- Build the two mapping tables ----
    ensp2symbol = build_human_ensp_to_symbol(human_faa)
    tx2gene = build_sparrow_transcript_to_gene(gtf_file)

    # ---- Parse Orthogroups.tsv ----
    # Column order: Orthogroup | Homo_sapiens | Passer_montanus
    print(f"\n[3/3] Parsing Orthogroups: {tsv_file}", file=sys.stderr)

    mappings = []
    sparrow_only = []
    human_only = []

    with open(tsv_file) as f:
        header = f.readline().strip().split("\t")
        print(f"  Columns: {header}", file=sys.stderr)
        # col[0]=Orthogroup, col[1]=Homo_sapiens, col[2]=Passer_montanus

        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 3:
                continue

            og = cols[0]
            hs_raw = [x.strip() for x in cols[1].split(",") if x.strip()]  # human
            sp_raw = [x.strip() for x in cols[2].split(",") if x.strip()]  # sparrow

            for s in sp_raw:
                if not hs_raw:
                    sparrow_only.append((og, s))
                for h in hs_raw:
                    mappings.append((og, s, h))

            if not sp_raw and hs_raw:
                for h in hs_raw:
                    human_only.append((og, h))

    # ---- Build indexes ----
    sp2hu = defaultdict(set)
    hu2sp = defaultdict(set)
    for og, s, h in mappings:
        sg = resolve_sparrow(s, tx2gene)
        hg = resolve_human(h, ensp2symbol)
        sp2hu[sg].add(hg)
        hu2sp[hg].add(sg)

    one2one = sum(1 for v in sp2hu.values() if len(v) == 1)

    # ---- All orthologs table ----
    out_all = f"{outdir}/sparrow_human_all_orthologs.tsv"
    with open(out_all, "w") as f:
        f.write("Orthogroup\tSparrow_Gene\tHuman_Symbol\tHuman_Ensembl\tRelation\n")
        for og, s, h in mappings:
            sg = resolve_sparrow(s, tx2gene)
            hs = resolve_human(h, ensp2symbol)
            h_ensp = re.sub(r'\.\d+$', '', h.split()[0])  # protein ID kept for reference

            ns, nh = len(sp2hu[sg]), len(hu2sp[hs])
            if ns == 1 and nh == 1:
                rel = "1-to-1"
            elif ns == 1 and nh > 1:
                rel = "1-to-many"
            elif ns > 1 and nh == 1:
                rel = "many-to-1"
            else:
                rel = "many-to-many"
            f.write(f"{og}\t{sg}\t{hs}\t{h_ensp}\t{rel}\n")
    print(f"\n  Wrote: {out_all}", file=sys.stderr)

    # ---- 1-to-1 table ----
    out_o2o = f"{outdir}/sparrow_human_one2one.tsv"
    with open(out_o2o, "w") as f:
        f.write("Sparrow_Gene\tHuman_Symbol\tHuman_Ensembl\tOrthogroup\n")
        seen = set()
        for og, s, h in mappings:
            sg = resolve_sparrow(s, tx2gene)
            hs = resolve_human(h, ensp2symbol)
            h_ensp = re.sub(r'\.\d+$', '', h.split()[0])
            pair = (sg, hs)

            if len(sp2hu[sg]) == 1 and len(hu2sp[hs]) == 1 and pair not in seen:
                seen.add(pair)
                f.write(f"{sg}\t{hs}\t{h_ensp}\t{og}\n")
    print(f"  Wrote: {out_o2o}", file=sys.stderr)

    # ---- Summary ----
    total_ogs = len(set(m[0] for m in mappings))
    summary = f"""\
{'=' * 55}
OrthoFinder: Passer montanus vs Homo sapiens
{'=' * 55}
Shared orthogroups:         {total_ogs}
Sparrow-specific groups:    {len(sparrow_only)}
Human-specific groups:      {len(human_only)}
Total ortholog pairs:       {len(mappings)}
Sparrow genes mapped:       {len(sp2hu)}
Human genes mapped:         {len(hu2sp)}
1-to-1 orthologs:           {one2one}
ENSP->Symbol in human FASTA: {len(ensp2symbol)}
Transcript->Gene in GTF:     {len(tx2gene)}
"""

    out_sum = f"{outdir}/ortholog_summary.txt"
    with open(out_sum, "w") as f:
        f.write(summary)

    print(summary, file=sys.stderr)
    print(f"  Wrote: {out_sum}", file=sys.stderr)

    # ---- Preview ----
    print("\n--- Sample 1-to-1 rows ---", file=sys.stderr)
    with open(out_o2o) as f:
        for i, line in enumerate(f):
            if i <= 5:
                print(f"  {line.rstrip()}", file=sys.stderr)


if __name__ == "__main__":
    main()
