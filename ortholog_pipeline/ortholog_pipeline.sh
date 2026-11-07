#!/bin/bash
# =============================================================================
# Ortholog identification between tree sparrow and human (OrthoFinder + DIAMOND)
#
# Required inputs (place in ${WORKDIR}/input/):
#   genome.fasta    tree sparrow genome assembly
#   annotation.gtf  tree sparrow annotation (CDS features with transcript_id/gene_id)
#   Homo_sapiens.GRCh38.pep.all.fa  Ensembl human proteome (Phase 2 prints the URL)
#
# Helper scripts (same directory as this script):
#   extract_proteins.py        CDS extraction + translation (pure Python)
#   filter_longest_isoform.py  keep the longest isoform per gene
#   parse_orthogroups.py       parse OrthoFinder Orthogroups.tsv -> symbol-annotated maps
#   parse_orthofinder.py       alternative parser for Orthologues/<species>.tsv
#
# Environment variables (all optional):
#   WORKDIR     working directory (default: <script dir>/work)
#   CONDA_ENV   conda env with orthofinder/diamond/seqkit (default: orthofinder_env)
#   THREADS     number of threads (default: all available cores)
#
# Usage: bash ortholog_pipeline.sh
# =============================================================================
set -euo pipefail

TOOLDIR="$(cd "$(dirname "$0")" && pwd)"
WORKDIR="${WORKDIR:-${TOOLDIR}/work}"
CONDA_ENV="${CONDA_ENV:-orthofinder_env}"
GENOME="${GENOME:-${WORKDIR}/input/genome.fasta}"
GTF="${GTF:-${WORKDIR}/input/annotation.gtf}"
THREADS="${THREADS:-$(nproc)}"
ENSEMBL_RELEASE="113"

mkdir -p "${WORKDIR}"/{logs,results,tmp}
export TMPDIR="${WORKDIR}/tmp"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

# ---- Helper: locate a Python script ----
find_py() {
    local name="$1"
    for d in "${WORKDIR}" "${TOOLDIR}"; do
        if [[ -f "${d}/${name}" ]]; then
            echo "${d}/${name}"
            return
        fi
    done
    echo "ERROR: ${name} not found in ${WORKDIR}/ or ${TOOLDIR}/" >&2
    exit 1
}

# ============================================================================
# Phase 0: Environment check
# ============================================================================
echo "===== Phase 0: Environment ====="
echo "Threads: ${THREADS}"
echo "RAM:     $(free -h | awk '/Mem/{print $2}')"
echo "Disk:    $(df -h ${WORKDIR} | tail -1 | awk '{print $4}')"
echo ""

for tool in orthofinder diamond python3; do
    command -v ${tool} &>/dev/null || { echo "MISSING: ${tool}"; exit 1; }
done
for f in "${GENOME}" "${GTF}"; do
    [[ -f "${f}" ]] || { echo "MISSING: ${f}"; exit 1; }
done
echo "GTF feature types:"
cut -f3 "${GTF}" | sort | uniq -c | sort -rn | head -10

# Locate all Python scripts
EXTRACT_PY=$(find_py "extract_proteins.py")
FILTER_PY=$(find_py "filter_longest_isoform.py")
PARSE_PY=$(find_py "parse_orthogroups.py")
echo ""
echo "Python scripts:"
echo "  extract:  ${EXTRACT_PY}"
echo "  filter:   ${FILTER_PY}"
echo "  parse:    ${PARSE_PY}"

# ============================================================================
# Phase 1: Extract tree sparrow proteins
# ============================================================================
echo ""
echo "===== Phase 1: Extracting sparrow proteins ====="

echo "[1.1] Extracting CDS and translating via Python..."
python3 "${EXTRACT_PY}" \
    --gtf "${GTF}" \
    --genome "${GENOME}" \
    --out-cds "${WORKDIR}/sparrow_cds.fna" \
    --out-prot "${WORKDIR}/sparrow_prot_raw.faa" \
    2>&1 | tee "${WORKDIR}/logs/extract_proteins.log"

RAW=$(grep -c '^>' "${WORKDIR}/sparrow_prot_raw.faa" 2>/dev/null || echo 0)
echo "  Raw proteins: ${RAW}"

if [[ ${RAW} -lt 1000 ]]; then
    echo "ERROR: Too few proteins extracted (${RAW})."
    echo "Check: ${WORKDIR}/logs/extract_proteins.log"
    exit 1
fi

echo ""
echo "[1.2] Sample headers:"
head -6 "${WORKDIR}/sparrow_prot_raw.faa"

echo ""
echo "[1.3] Filtering longest isoform per gene..."
python3 "${FILTER_PY}" \
    "${WORKDIR}/sparrow_prot_raw.faa" \
    "${WORKDIR}/sparrow_prot_longest.faa"

SPARROW=$(grep -c '^>' "${WORKDIR}/sparrow_prot_longest.faa")
echo "  Sparrow genes retained: ${SPARROW}"

echo ""
echo "[1.4] Protein quality check:"
echo "  Amino acid composition (top 10):"
grep -v '^>' "${WORKDIR}/sparrow_prot_longest.faa" | grep -o . | sort | uniq -c | sort -rn | head -10
echo "  Internal stops: $(grep -v '^>' "${WORKDIR}/sparrow_prot_longest.faa" | awk '{sub(/\*$/,""); if (/\*/) n++} END{print n+0}')"

# ============================================================================
# Phase 2: Human reference proteome
# ============================================================================
echo ""
echo "===== Phase 2: Human proteome ====="

# Look for the human proteome file by priority
HUMAN=""
for f in \
    "${WORKDIR}/Homo_sapiens.GRCh38.pep.all.fa" \
    "${WORKDIR}/input/Homo_sapiens.GRCh38.pep.all.fa"; do
    if [[ -f "${f}" ]]; then
        HUMAN="${f}"
        break
    fi
done

if [[ -z "${HUMAN}" ]]; then
    echo "ERROR: Human proteome file not found!"
    echo ""
    echo "Please download and place it as:"
    echo "  ${WORKDIR}/Homo_sapiens.GRCh38.pep.all.fa"
    echo ""
    echo "Download URL (on a machine with internet):"
    echo "  wget https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"
    echo "  gunzip Homo_sapiens.GRCh38.pep.all.fa.gz"
    echo ""
    echo "Then re-run this script."
    exit 1
fi

echo "Using: ${HUMAN}"

echo "[2.1] Filtering longest isoform..."
python3 "${FILTER_PY}" \
    "${HUMAN}" \
    "${WORKDIR}/Homo_sapiens.pep.longest.faa"

HUMAN_N=$(grep -c '^>' "${WORKDIR}/Homo_sapiens.pep.longest.faa")
echo "  Human genes retained: ${HUMAN_N}"

# ============================================================================
# Phase 3: Run OrthoFinder
# ============================================================================
echo ""
echo "===== Phase 3: OrthoFinder ====="

# Prepare input directory (one FASTA per species)
ORTHO_IN="${WORKDIR}/orthofinder_input"
rm -rf "${ORTHO_IN}"
mkdir -p "${ORTHO_IN}"
cp "${WORKDIR}/sparrow_prot_longest.faa" "${ORTHO_IN}/Passer_montanus.faa"
cp "${WORKDIR}/Homo_sapiens.pep.longest.faa" "${ORTHO_IN}/Homo_sapiens.faa"

echo "Input files:"
ls -lh "${ORTHO_IN}/"

echo ""
echo "Starting OrthoFinder (DIAMOND --ultra-sensitive, ${THREADS} threads)..."
echo "  Estimated runtime: 30-90 minutes"
echo "  Monitor progress:  tail -f ${WORKDIR}/logs/orthofinder_run.log"
echo ""

cd "${WORKDIR}"

orthofinder \
    -f "${ORTHO_IN}" \
    -t "${THREADS}" \
    -a "${THREADS}" \
    -S diamond \
    -og \
    -o "${WORKDIR}/results" \
    2>&1 | tee "${WORKDIR}/logs/orthofinder_run.log"

# Check the output directory
RESULTS_DIR=$(ls -d "${WORKDIR}/results/Results_"* 2>/dev/null | tail -1)

if [[ -z "${RESULTS_DIR}" ]]; then
    echo ""
    echo "ERROR: OrthoFinder did not produce a Results directory."
    echo "Check the log: ${WORKDIR}/logs/orthofinder_run.log"
    exit 1
fi

echo ""
echo "OrthoFinder finished! Results: ${RESULTS_DIR}"

# ============================================================================
# Phase 4: Parse results
# ============================================================================
echo ""
echo "===== Phase 4: Parse orthologs ====="

# parse_orthogroups.py takes Orthogroups.tsv and needs the human proteome (for
# ENSP -> gene symbol) and the GTF (for transcript_id -> gene_id).
OG_FILE=$(find "${RESULTS_DIR}/Orthogroups" -name "Orthogroups.tsv" 2>/dev/null | head -1)

if [[ -z "${OG_FILE}" ]]; then
    echo "ERROR: No Orthogroups.tsv found!"
    echo "Contents of ${RESULTS_DIR}/Orthogroups/:"
    find "${RESULTS_DIR}/Orthogroups" -type f
    exit 1
fi

echo "Parsing: ${OG_FILE}"

python3 "${PARSE_PY}" "${OG_FILE}" "${HUMAN}" "${GTF}" "${WORKDIR}/results"

# ============================================================================
# Done
# ============================================================================
echo ""
echo "================================================"
cat "${WORKDIR}/results/ortholog_summary.txt"
echo "================================================"
echo ""
echo "===== PIPELINE COMPLETE ====="
echo ""
echo "Output files:"
echo "  1-to-1:  ${WORKDIR}/results/sparrow_human_one2one.tsv"
echo "  All:     ${WORKDIR}/results/sparrow_human_all_orthologs.tsv"
echo "  Summary: ${WORKDIR}/results/ortholog_summary.txt"
