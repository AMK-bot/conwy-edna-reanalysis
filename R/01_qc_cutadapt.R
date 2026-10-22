#!/usr/bin/env Rscript
# 01_qc_cutadapt.R  --  Phase 2: read QC + MiFish-U primer removal.
# Plan §3 Phase 2.  Input: data/ena/raw/*.fastq.gz  ->  Output: data/qc/trimmed/*.fastq.gz
# Orchestrates CLI tools (run inside the env: `micromamba run -n edna ...`).
#
# Steps:
#   1. FastQC + MultiQC on raw reads (read the per-base quality curves; high
#      duplication is normal for amplicons).
#   2. cutadapt as LINKED adapters, --discard-untrimmed, IUPAC degeneracy direct:
#        F  GTCGGTAAAACTCGTGCCAGC
#        R  CATAGTGGGGTATCTAATCCCAGTTTG
#   3. Hand truncLen/maxEE choices to 02_dada2.R based on the quality curves.
#
# Reading: cutadapt linked-adapters https://cutadapt.readthedocs.io/  · FastQC/MultiQC docs.
# TODO: implement (cutadapt/fastqc/multiqc are not yet in the env — add when Phase 2 starts).
stop("01_qc_cutadapt.R is a Phase-2 stub — not yet implemented.")
