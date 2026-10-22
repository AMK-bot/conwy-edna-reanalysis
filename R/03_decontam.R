#!/usr/bin/env Rscript
# 03_decontam.R  --  Phase 2: decontamination + tag-jump control (where eDNA rigor lives).
# Plan §3 Phase 2.  Input: data/asv/seqtab.rds + controls  ->  Output: data/asv/seqtab_clean.rds
#
# Steps: decontam (prevalence method) using the 66 field + 68 lab negative controls
#   (BLANK*/EBlank/ENC/PC families in PRJEB48362); an abundance threshold for
#   index-hopping WITH a sensitivity analysis; metabaR optional for tag-jump modeling.
# State every threshold and sensitivity-test it (rigor checklist, plan §5).
#
# Reading: decontam vignette https://benjjneb.github.io/decontam/ · MIEM reporting guidelines.
stop("03_decontam.R is a Phase-2 stub — not yet implemented.")
