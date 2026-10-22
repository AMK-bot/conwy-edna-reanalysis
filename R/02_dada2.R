#!/usr/bin/env Rscript
# 02_dada2.R  --  Phase 2: DADA2 denoising -> ASV table.
# Plan §3 Phase 2.  Input: data/qc/trimmed/  ->  Output: data/asv/seqtab.rds, repseqs.fasta
#
# Steps (DADA2): filterAndTrim (truncLen/maxEE from QC curves; keep >=12 bp merge
#   overlap) -> learnErrors PER SEQUENCING RUN -> dada(pool="pseudo") -> mergePairs
#   -> makeSequenceTable -> removeBimeraDenovo(method="consensus").
# Validation targets (authors): ~2726 12S ASVs pre-curation; ASVs >20% longer than
#   172 bp discarded. Note version drift (authors used DADA2 v1.14.1; env has 1.34.0).
#
# Reading: DADA2 tutorial https://benjjneb.github.io/dada2/tutorial.html · Callahan 2017 (ASVs).
stop("02_dada2.R is a Phase-2 stub — not yet implemented.")
