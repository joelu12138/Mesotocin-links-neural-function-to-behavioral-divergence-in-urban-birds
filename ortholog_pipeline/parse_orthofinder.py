#!/usr/bin/env python3
"""
Alternative parser: parse an OrthoFinder Orthologues TSV into mapping tables.

Usage:
  python3 parse_orthofinder.py <orthologues.tsv> <output_dir>

Outputs (in <output_dir>):
  sparrow_human_one2one.tsv      1-to-1 orthologs
  sparrow_human_all_orthologs.tsv all ortholog pairs with a relation label
  ortholog_summary.txt           summary statistics

This parser keeps raw accessions (no gene-symbol translation). Use
parse_orthogroups.py when gene symbols are needed.
"""
import re
import sys
from collections import defaultdict


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <orthologues.tsv> <output_dir>", file=sys.stderr)
        sys.exit(1)

    tsv_file = sys.argv[1]
    outdir = sys.argv[2]

    print(f"Parsing: {tsv_file}", file=sys.stderr)

    mappings = []       # (orthogroup, sparrow_gene, human_gene)
    sparrow_only = []

    with open(tsv_file) as f:
        f.readline()  # skip header
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 3:
                continue
            og = cols[0]
            sp = [x.strip() for x in cols[1].split(",") if x.strip()]
            hs = [x.strip() for x in cols[2].split(",") if x.strip()]

            for s in sp:
                if not hs:
                    sparrow_only.append((og, s))
                for h in hs:
                    mappings.append((og, s, h))

    # Build indexes
    sp2hu = defaultdict(set)
    hu2sp = defaultdict(set)
    for og, s, h in mappings:
        sp2hu[s].add(h)
        hu2sp[h].add(s)

    one2one = sum(1 for v in sp2hu.values() if len(v) == 1)

    # ---- All orthologs table ----
    out_all = f"{outdir}/sparrow_human_all_orthologs.tsv"
    with open(out_all, "w") as f:
        f.write("Orthogroup\tSparrow_Gene\tHuman_Ensembl_Gene\tRelation\n")
        for og, s, h in mappings:
            hc = re.sub(r"\.\d+$", "", h)
            ns, nh = len(sp2hu[s]), len(hu2sp[h])
            if ns == 1 and nh == 1:
                rel = "1-to-1"
            elif ns == 1 and nh > 1:
                rel = "1-to-many"
            elif ns > 1 and nh == 1:
                rel = "many-to-1"
            else:
                rel = "many-to-many"
            f.write(f"{og}\t{s}\t{hc}\t{rel}\n")
    print(f"  Wrote: {out_all}", file=sys.stderr)

    # ---- 1-to-1 table ----
    out_o2o = f"{outdir}/sparrow_human_one2one.tsv"
    with open(out_o2o, "w") as f:
        f.write("Sparrow_Gene\tHuman_Ensembl_Gene\tOrthogroup\n")
        seen = set()
        for og, s, h in mappings:
            hc = re.sub(r"\.\d+$", "", h)
            if len(sp2hu[s]) == 1 and len(hu2sp[h]) == 1 and (s, hc) not in seen:
                seen.add((s, hc))
                f.write(f"{s}\t{hc}\t{og}\n")
    print(f"  Wrote: {out_o2o}", file=sys.stderr)

    # ---- Summary ----
    summary = f"""\
{'=' * 50}
OrthoFinder Results: Passer montanus vs Homo sapiens
{'=' * 50}
Shared orthogroups:            {len(set(m[0] for m in mappings))}
Sparrow-only groups:           {len(sparrow_only)}
Total ortholog pairs:          {len(mappings)}
Sparrow genes mapped:          {len(sp2hu)}
Human genes mapped:            {len(hu2sp)}
1-to-1 orthologs:              {one2one}
1-to-many (sparrow view):      {len(sp2hu) - one2one}
"""
    out_sum = f"{outdir}/ortholog_summary.txt"
    with open(out_sum, "w") as f:
        f.write(summary)

    print(summary, file=sys.stderr)
    print(f"  Wrote: {out_sum}", file=sys.stderr)


if __name__ == "__main__":
    main()
