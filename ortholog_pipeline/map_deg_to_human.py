# -*- coding: utf-8 -*-
# Map a tree sparrow DEG table (gene ID in the first column) to human ortholog SYMBOLs.
# Usage: from the repository root, run  python ortholog_pipeline/map_deg_to_human.py
import csv
import sys
from collections import defaultdict

DEG_IN = 'data/DEG_results.csv'          # input DEG table
OUT    = 'data/DEG_results_human.csv'    # output
ALLMAP = 'ortholog_pipeline/sparrow_human_all_orthologs.tsv'

# 1. Read the ortholog map: sparrow gene -> set of human symbols
sym = defaultdict(set)
with open(ALLMAP) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        sym[r['Sparrow_Gene']].add(r['Human_Symbol'])

# 2. Read the DEG table (first column is the gene ID, header cell is empty)
with open(DEG_IN) as f:
    rd = csv.reader(f)
    header = next(rd)
    header[0] = 'gene'                       # fix the rowname column name
    out_rows = [['gene', 'Human_Symbol', 'map_type'] + header[1:]]
    for row in rd:
        g = row[0]
        if not g:
            continue
        ids = sorted(sym.get(g, set()))
        symbols = [s for s in ids if not s.startswith('ENSP')]
        if symbols:
            hsym = ';'.join(symbols)
            mtype = '1:1' if len(symbols) == 1 else 'many'
        elif ids:                            # only ENSP IDs, no official symbol
            hsym = ';'.join(ids)
            mtype = 'many'
        elif not g.startswith('LOC'):        # already a gene symbol, no ortholog record
            hsym, mtype = g, 'self'
        else:
            hsym, mtype = '', 'none'
        out_rows.append([g, hsym, mtype] + row[1:])

# 3. Write output
with open(OUT, 'w', newline='') as f:
    csv.writer(f).writerows(out_rows)

from collections import Counter
c = Counter(r[2] for r in out_rows[1:])
print(f'Wrote {OUT}: {len(out_rows)-1} genes')
print(f'1:1={c["1:1"]}  many={c["many"]}  self={c["self"]}  none={c["none"]}')
