#!/usr/bin/env bash
# Build one append-only, client-facing P42 package including the requested data.
set -euo pipefail

PROJECT=/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis
SCRIPT_DIR="${PROJECT}/scripts/20261024"
MANIFEST="${PROJECT}/analysis/20261024_customer_request_response_v1/data_delivery_availability_manifest.tsv"
DELIVERY="${PROJECT}/custom_research_report_20261024"
AUDIT_DIR="${PROJECT}/analysis/20261024_client_package_self_audit"
LOG_DIR="${PROJECT}/logs/20261024_client_package"
RUN_ID=$(date +%Y%m%dT%H%M%S%z)
LOG_FILE="${LOG_DIR}/client_package_${RUN_ID}.log"
WORK_DIR="${LOG_DIR}/work_${RUN_ID}"
PY=/usr/bin/python3
SAMTOOLS=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/samtools
BCFTOOLS=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/bcftools
TABIX=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/tabix

mkdir -p "${LOG_DIR}" "${WORK_DIR}"
exec >"${LOG_FILE}" 2>&1

die() { echo "[$(date --iso-8601=seconds)] ERROR: $*" >&2; exit 1; }
copy_data() {
    local source=$1 destination=$2
    [[ -s "${source}" ]] || die "Source file is missing or empty: ${source}"
    mkdir -p "$(dirname "${destination}")"
    rsync -a --partial --append-verify --info=progress2 "${source}" "${destination}"
}
verify_identical() {
    local source=$1 destination=$2 label=$3
    local output="${WORK_DIR}/verify_${label}.txt"
    rsync -acn --out-format='%i %n' "${source}" "${destination}" >"${output}"
    [[ ! -s "${output}" ]] || die "Checksum verification differs for ${label}; see ${output}"
}
sanitize_bam_header() {
    local source=$1 header=$2 removed_file=$3
    "${SAMTOOLS}" view --no-PG -H "${source}" | awk -v removed_file="${removed_file}" '
        /^@PG/ { removed++; next }
        (index($0, "/home/") || index($0, "/mnt/") || index($0, "/Work_bio/") || index($0, "file:///")) { removed++; next }
        { print }
        END { print removed + 0 > removed_file }
    ' >"${header}"
}
sanitize_bam_source() {
    local source=$1 target=$2 sample_id=$3 header removed temporary source_count target_count
    header="${WORK_DIR}/${sample_id}.bam.header.sam"
    removed="${WORK_DIR}/${sample_id}.bam_header_removed.txt"
    temporary="${WORK_DIR}/${sample_id}.client.bam"
    sanitize_bam_header "${source}" "${header}" "${removed}"
    "${SAMTOOLS}" reheader -P "${header}" "${source}" >"${temporary}"
    "${SAMTOOLS}" index -@ 4 "${temporary}"
    "${SAMTOOLS}" quickcheck -v "${temporary}"
    source_count=$("${SAMTOOLS}" view -c "${source}")
    target_count=$("${SAMTOOLS}" view -c "${temporary}")
    [[ "${source_count}" == "${target_count}" ]] || die "Alignment-count mismatch for sanitized BAM ${sample_id}"
    mv "${temporary}" "${target}"
    mv "${temporary}.bai" "${target}.bai"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\talignment_records_preserved_header_paths_removed\n' \
        "${sample_id}" "${target#${DELIVERY}/}" "$(stat -c %s "${source}")" "$(stat -c %s "${target}")" \
        "${source_count}" "${target_count}" "$(<"${removed}")" >>"${WORK_DIR}/bam_header_sanitization_audit.tsv"
}
sanitize_header() {
    local source=$1 header=$2 removed_file=$3
    "${BCFTOOLS}" view --no-version -h "${source}" | awk -v removed_file="${removed_file}" '
        /^##reference=/ { normalized++; print "##reference=GRCh38"; next }
        /^##/ && (index($0, "/home/") || index($0, "/mnt/") || index($0, "/Work_bio/") || index($0, "file:///") || $0 ~ /^##cmdline=/ || $0 ~ /^##bcftools_/) { removed++; next }
        { print }
        END { print removed + normalized > removed_file }
    ' >"${header}"
}
variant_body_hash() {
    "${BCFTOOLS}" view --no-version -H "$1" | sha256sum | awk '{print $1}'
}
audit_sanitized_vcf() {
    local source=$1 target=$2 sample_id=$3 file_type=$4 removed_file=$5
    local source_hash target_hash removed_count
    source_hash=$(variant_body_hash "${source}")
    target_hash=$(variant_body_hash "${target}")
    [[ "${source_hash}" == "${target_hash}" ]] || die "Variant-body hash mismatch for ${sample_id} ${file_type}"
    removed_count=$(<"${removed_file}")
    printf '%s\t%s\t%s\t%s\t%s\t%s\tvariant_records_preserved_header_paths_removed\n' \
        "${sample_id}" "${file_type}" "${target#${DELIVERY}/}" "${source_hash}" "${target_hash}" "${removed_count}" \
        >>"${WORK_DIR}/vcf_header_sanitization_audit.tsv"
}
sanitize_vcf_source() {
    local source=$1 target=$2 sample_id=$3 file_type=$4 stem header removed temporary
    stem="${sample_id}_${file_type}"
    header="${WORK_DIR}/${stem}.header.txt"
    removed="${WORK_DIR}/${stem}.removed.txt"
    temporary="${WORK_DIR}/${stem}.vcf.gz"
    sanitize_header "${source}" "${header}" "${removed}"
    "${BCFTOOLS}" reheader -h "${header}" -o "${temporary}" "${source}"
    "${TABIX}" -f -p vcf "${temporary}"
    mv "${temporary}" "${target}"
    mv "${temporary}.tbi" "${target}.tbi"
    audit_sanitized_vcf "${source}" "${target}" "${sample_id}" "${file_type}" "${removed}"
}
sanitize_bcf_source() {
    local source=$1 target=$2 sample_id=$3 file_type=$4 stem header removed temporary_bcf temporary_vcf
    stem="${sample_id}_${file_type}"
    header="${WORK_DIR}/${stem}.header.txt"
    removed="${WORK_DIR}/${stem}.removed.txt"
    temporary_bcf="${WORK_DIR}/${stem}.bcf"
    temporary_vcf="${WORK_DIR}/${stem}.vcf.gz"
    sanitize_header "${source}" "${header}" "${removed}"
    "${BCFTOOLS}" reheader -h "${header}" -o "${temporary_bcf}" "${source}"
    "${BCFTOOLS}" view --no-version -Oz -o "${temporary_vcf}" "${temporary_bcf}"
    "${TABIX}" -f -p vcf "${temporary_vcf}"
    mv "${temporary_vcf}" "${target}"
    mv "${temporary_vcf}.tbi" "${target}.tbi"
    audit_sanitized_vcf "${source}" "${target}" "${sample_id}" "${file_type}" "${removed}"
}

[[ -f "${MANIFEST}" ]] || die "Missing verified source manifest: ${MANIFEST}"
[[ ! -e "${AUDIT_DIR}" ]] || die "Refusing to overwrite an existing audit directory: ${AUDIT_DIR}"
if [[ -e "${DELIVERY}/${DELIVERY##*/}.md" || -e "${DELIVERY}/P42_Vedita_18_WGS_1024.md" ]]; then
    die "Refusing to overwrite a completed report in ${DELIVERY}"
fi
if [[ -e "${DELIVERY}" ]]; then
    echo "[$(date --iso-8601=seconds)] Resuming only an unfinished package without a root report: ${DELIVERY}"
else
    mkdir -p "${DELIVERY}"
fi
mkdir -p "${DELIVERY}/figures" "${DELIVERY}/results" "${DELIVERY}/methods" \
    "${DELIVERY}/genome_wide_data/bam" \
    "${DELIVERY}/genome_wide_data/manta/filtered_sv" \
    "${DELIVERY}/genome_wide_data/manta/candidate_sv" \
    "${DELIVERY}/genome_wide_data/manta/candidate_small_indel" \
    "${DELIVERY}/genome_wide_data/delly/discovery" \
    "${DELIVERY}/genome_wide_data/delly/region_genotyped"

echo "[$(date --iso-8601=seconds)] Starting transfer into ${DELIVERY##*/}"
printf 'sample_id\tstage\tsource_type\tclient_relative_path\tstatus\n' >"${WORK_DIR}/delly_conversion_audit.tsv"
printf 'sample_id\tfile_type\tclient_relative_path\tsource_variant_body_sha256\tclient_variant_body_sha256\tsanitized_header_lines\tstatus\n' >"${WORK_DIR}/vcf_header_sanitization_audit.tsv"
printf 'sample_id\tclient_relative_path\tsource_bytes\tclient_bytes\tsource_alignment_records\tclient_alignment_records\tsanitized_header_lines\tstatus\n' >"${WORK_DIR}/bam_header_sanitization_audit.tsv"

while IFS=$'\t' read -r sample_id file_type source_path file_exists size_bytes index_path index_exists note; do
    [[ "${sample_id}" != "sample_id" ]] || continue
    [[ "${file_exists}" == "1" ]] || continue
    case "${file_type}" in
        BAM)
            target="${DELIVERY}/genome_wide_data/bam/$(basename "${source_path}")"
            [[ ! -e "${target}" && ! -e "${target}.bai" ]] || die "Refusing to overwrite sanitized BAM target: ${target}"
            sanitize_bam_source "${source_path}" "${target}" "${sample_id}"
            ;;
        Manta_filtered_SV)
            target="${DELIVERY}/genome_wide_data/manta/filtered_sv/${sample_id}.$(basename "${source_path}")"
            [[ ! -e "${target}" && ! -e "${target}.tbi" ]] || die "Refusing to overwrite sanitized Manta target: ${target}"
            sanitize_vcf_source "${source_path}" "${target}" "${sample_id}" "Manta_filtered_SV"
            ;;
        Manta_candidate_SV)
            target="${DELIVERY}/genome_wide_data/manta/candidate_sv/${sample_id}.$(basename "${source_path}")"
            [[ ! -e "${target}" && ! -e "${target}.tbi" ]] || die "Refusing to overwrite sanitized Manta target: ${target}"
            sanitize_vcf_source "${source_path}" "${target}" "${sample_id}" "Manta_candidate_SV"
            ;;
        Manta_candidate_small_indel)
            target="${DELIVERY}/genome_wide_data/manta/candidate_small_indel/${sample_id}.$(basename "${source_path}")"
            [[ ! -e "${target}" && ! -e "${target}.tbi" ]] || die "Refusing to overwrite sanitized Manta target: ${target}"
            sanitize_vcf_source "${source_path}" "${target}" "${sample_id}" "Manta_candidate_small_indel"
            ;;
        DELLY_discovery_BCF|DELLY_region_genotyped_BCF)
            if [[ "${file_type}" == "DELLY_discovery_BCF" ]]; then
                category=discovery
            else
                category=region_genotyped
            fi
            target="${DELIVERY}/genome_wide_data/delly/${category}/${sample_id}.${category}.vcf.gz"
            if [[ -e "${target}" || -e "${target}.tbi" ]]; then
                die "Refusing to overwrite existing DELLY conversion target: ${target}"
            fi
            sanitize_bcf_source "${source_path}" "${target}" "${sample_id}" "${file_type}"
            printf '%s\t%s\t%s\t%s\tconverted_and_indexed\n' "${sample_id}" "${category}" "${file_type}" "genome_wide_data/delly/${category}/$(basename "${target}")" >>"${WORK_DIR}/delly_conversion_audit.tsv"
            ;;
    esac
done <"${MANIFEST}"

expected_conversion_rows=36
actual_conversion_rows=$(( $(wc -l <"${WORK_DIR}/delly_conversion_audit.tsv") - 1 ))
[[ "${actual_conversion_rows}" -eq "${expected_conversion_rows}" ]] || die "DELLY conversion audit has ${actual_conversion_rows}, expected ${expected_conversion_rows} rows"
expected_bam_audit_rows=18
actual_bam_audit_rows=$(( $(wc -l <"${WORK_DIR}/bam_header_sanitization_audit.tsv") - 1 ))
[[ "${actual_bam_audit_rows}" -eq "${expected_bam_audit_rows}" ]] || die "BAM header audit has ${actual_bam_audit_rows}, expected ${expected_bam_audit_rows} rows"
expected_header_audit_rows=90
actual_header_audit_rows=$(( $(wc -l <"${WORK_DIR}/vcf_header_sanitization_audit.tsv") - 1 ))
[[ "${actual_header_audit_rows}" -eq "${expected_header_audit_rows}" ]] || die "VCF header audit has ${actual_header_audit_rows}, expected ${expected_header_audit_rows} rows"
"${PY}" "${SCRIPT_DIR}/10_build_unified_client_package.py" --delivery "${DELIVERY}" --conversion-audit "${WORK_DIR}/delly_conversion_audit.tsv" --header-sanitization-audit "${WORK_DIR}/vcf_header_sanitization_audit.tsv" --bam-sanitization-audit "${WORK_DIR}/bam_header_sanitization_audit.tsv"
"${PY}" "${SCRIPT_DIR}/11_validate_unified_client_package.py" --delivery "${DELIVERY}" --audit-dir "${AUDIT_DIR}" --samtools "${SAMTOOLS}" --bcftools "${BCFTOOLS}"
echo "[$(date --iso-8601=seconds)] Succeeded: unified client package built and self-audited"
