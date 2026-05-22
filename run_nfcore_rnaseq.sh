#!/usr/bin/env bash
set -euo pipefail

# Usage:
# MOUSE_FASTA=/absolute/path/to/GRCm39.primary_assembly.genome.fa \
# MOUSE_GFF3=/absolute/path/to/gencode.vM38.basic.annotation.gff3 \
# NFCORE_WORK_DIR=/absolute/path/to/nextflow_work \
# ./run_nfcore_rnaseq.sh

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
files_dir="${project_dir}/files"

mouse_fasta="${MOUSE_FASTA}"
mouse_gff3="${MOUSE_GFF3}"
work_dir="${NFCORE_WORK_DIR:-${project_dir}/work}"

samplesheet="${project_dir}/samplesheet.csv"
outdir="${project_dir}/res"

sars2_fasta="${files_dir}/SARS2sequence.fasta"
narita1_fasta="${files_dir}/Narita1sequence.fasta"
sars2_gff3="${files_dir}/SARS2sequence.gff3"
narita1_gff3="${files_dir}/Narita1sequence.gff3"

combined_fasta="${files_dir}/mouse_virus.fa"
combined_gff3="${files_dir}/mouse_virus.gff3"
combined_gtf="${files_dir}/mouse_virus.gtf"
fixed_gtf="${files_dir}/mouse_virus.fixed.gtf"

cat "${mouse_fasta}" \
    "${sars2_fasta}" \
    "${narita1_fasta}" \
    > "${combined_fasta}"

{
    grep "^#" "${mouse_gff3}"
    grep -v "^#" "${mouse_gff3}"
    grep -v "^#" "${sars2_gff3}"
    grep -v "^#" "${narita1_gff3}"
} > "${combined_gff3}"

gffread "${combined_gff3}" -T -o "${combined_gtf}"

awk -F'\t' '
BEGIN { OFS = "\t" }
{
    attributes = $9
    if (attributes ~ /transcript_id "/ && attributes !~ /gene_id "/) {
        transcript_id = attributes
        sub(/.*transcript_id "/, "", transcript_id)
        sub(/".*/, "", transcript_id)
        sub(/;[[:space:]]*$/, "", attributes)
        $9 = attributes "; gene_id \"" transcript_id "\";"
    }
    print
}
' "${combined_gtf}" > "${fixed_gtf}"

nextflow run nf-core/rnaseq \
    --input "${samplesheet}" \
    --outdir "${outdir}" \
    --fasta "${combined_fasta}" \
    --gtf "${fixed_gtf}" \
    -r 3.23.0 \
    -profile docker \
    -work-dir "${work_dir}"
