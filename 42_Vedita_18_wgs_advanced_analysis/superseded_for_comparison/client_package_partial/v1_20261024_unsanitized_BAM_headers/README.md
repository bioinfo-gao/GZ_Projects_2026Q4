# Archived incomplete client-package build: unsanitized BAM headers

## Status

This is an internal, non-deliverable partial build archived on 2026-10-24 before any client-facing report or VCF was created. It must never be supplied to a client.

## Why it was archived

The partial build copied source BAM files byte-for-byte. Inspection showed that their `@PG` header records contain internal execution paths and raw FASTQ paths. A client-ready package must not expose those paths. The partial package was therefore moved intact from `custom_research_report_20261024/` into this archive; no file was deleted and no source BAM was modified.

## What is preserved

| Item | Purpose |
| :--- | :---: |
| `custom_research_report_20261024_partial/` | The incomplete 140-GiB transfer state, including copied/partial BAM files; retained only for comparison and recovery audit |
| `initial_partial_copy.log` | First transfer log before VCF header-path screening was added |
| `first_resumed_copy.log` | First resumed-transfer log after VCF screening was added but before BAM header screening |
| `run_unified_client_package_corrective_snapshot.sh` | Corrective script snapshot retained with the archive; it documents the replacement build approach but is not the exact unversioned in-memory script that produced the initial partial copy |

## Corrected build

The replacement build creates a fresh `custom_research_report_20261024/`. It removes `@PG` and any remaining path-bearing BAM header lines using `samtools reheader -P`, rebuilds the BAM index, verifies BAM readability, and checks source/client alignment-record counts. It similarly strips local execution paths from Manta/DELLY VCF headers while preserving and hashing their variant bodies. The corrected package is the only candidate for client review.
