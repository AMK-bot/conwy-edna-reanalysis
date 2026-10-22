#!/usr/bin/env Rscript
# =============================================================================
# R/10_sim_recovery.R  --  Phase-3 SIMULATION RECOVERY STUDY for the occupancy model
# -----------------------------------------------------------------------------
# QUESTION (eDNAplan.md Sec 2.2, the governing caveat): across recovered taxa and
#   642 shared libraries, the reprocessed/source-study pooled detection-cell ratio
#   is ~84%. Miss indicators are associated with read depth, source-study relative
#   abundance, and the upstream-to-downstream biological-signal gradient. No
#   independent biomass was measured and no occupancy-bias direction is assumed.
#
#   This study asks, by simulating from occumb's OWN generative model with a
#   KNOWN psi and intentionally injecting a depth x ecological-gradient confound:
#     Does occumb's read-depth detection safeguard (theta ~ logdepth) recover the
#     true psi under the confound, or does a residual bias remain at low psi?
#
# DESIGN (deliberately SMALL so the whole study runs in a few minutes on 8 GB):
#   * occumb's generative model (Fukaya et al. 2022), exactly as occumb() fits it:
#       z[i,j]    ~ Bernoulli(psi[i,j])           site occupancy (the estimand)
#       u[i,j,k]  ~ Bernoulli(z * theta[i,j,k])   sequence capture (detection)
#       r[i,j,k]  ~ Gamma(phi[i], 1)              relative sequence dominance
#       y[.,j,k]  ~ Multinomial(N[j,k], pi)       read counts; N = per-rep DEPTH
#   * I = 6 species with KNOWN psi spanning rare->common (0.15 .. 0.85 at mean x).
#   * J = 16 site x timepoint units, K = 3 replicates.
#   * THE CONFOUND, injected on purpose:
#       - site covariate x  (a latent ecological gradient) drives psi: logit psi = a_i + x
#       - sequencing depth N tracks x too:  log10 N = base + gN * x   (deep downstream)
#       - CAPTURE prob theta drops at low depth: logit theta = c_i + sth * logdepth
#     so low-x cells are simultaneously low-psi AND low-depth AND under-detected --
#     creating a deliberately difficult known-truth scenario rather than a fitted
#     biological model of the Conwy.
#   * The simulation also gives rare species smaller sequence-dominance phi -> fewer
#     reads -> possible zero-read non-detections even when captured. This is a
#     synthetic residual channel the depth safeguard cannot directly identify.
#
# TWO fits per simulated dataset (both keep psi ~ x, the occupancy covariate):
#   (A) SAFEGUARD : formula_theta_shared = ~ logdepth   (depth-aware detection)
#   (B) NAIVE     : formula_theta        = ~ 1          (depth-blind detection)
#
# REPORTS, per true-psi bin x model: bias (mean psi_hat - psi_true), RMSE, and
#   95% CrI coverage, but only as an internal diagnostic. A convergence receipt
#   gates the entire numerical result. If any fit-level max psi Rhat exceeds 1.05,
#   tables are marked WITHHELD_UNCONVERGED and the public figure shows only the
#   failed gate; direction and magnitude must not be reported as findings.
#
# OUTPUTS (idempotent -- overwritten each run):
#   data/derived/sim_recovery_results.tsv
#   data/metadata/layerA_recovery_status.tsv
#   figures/phase3_sim_recovery.png
#
# HARD RULES honoured: runs SYNCHRONOUSLY, single R process, small + short MCMC.
# =============================================================================

suppressMessages(library(occumb))
set.seed(20260705)

# ---- study size knobs (kept small; shrink further if runtime creeps up) -----
I      <- 6L          # species
J      <- 16L         # site x timepoint units
K      <- 3L          # replicates per unit
NSIM   <- 10L         # simulated datasets (each fit twice: safeguard + naive)
N_CHAINS <- 2L; N_ADAPT <- 300L; N_BURN <- 1800L; N_ITER <- 3800L; N_THIN <- 2L
PRIOR_ULIM <- 6       # regularise among-species SD (default 1e4 is too diffuse
PRIOR_PREC <- 0.1     #   for 6 species and lets MCMC wander -> poor mixing)

# ---- fixed true generative parameters ---------------------------------------
psi_level <- seq(0.15, 0.85, length.out = I)   # KNOWN psi at mean gradient (x=0)
a_i       <- qlogis(psi_level)                 # species psi intercepts
b_psi     <- 1.0                               # latent ecological-gradient effect on psi
x_design  <- as.numeric(scale(seq(-1.5, 1.5, length.out = J)))  # ecological gradient
# --- the confound + deficit (tuned to ~86% overall, ~62% upstream sensitivity)
base_ld   <- 2.95     # log10 mean reads at mean gradient
g_depth   <- 0.7      # depth tracks the gradient: deep downstream, shallow upstream
s_theta   <- 1.4      # capture prob falls with depth  (secondary, depth-driven channel)
ct_i      <- qlogis(0.45) + 0.2 * (psi_level - mean(psi_level))  # capture intercept
phi_i     <- exp(-0.4 + 1.3 * (psi_level - mean(psi_level)))     # synthetic sequence-
                                                                  # dominance channel

# ---- simulate ONE dataset from occumb's generative model --------------------
simulate_one <- function() {
  psi_true <- matrix(0, I, J)
  for (j in seq_len(J)) psi_true[, j] <- plogis(a_i + b_psi * x_design[j])
  z_true <- matrix(rbinom(I * J, 1, psi_true), I, J)          # occupancy (per unit)

  logN <- base_ld + g_depth * x_design + rnorm(J * K, 0, 0.2) # depth ~ gradient
  N    <- pmax(round(10^matrix(logN, J, K)), 8L)

  y  <- array(0L, dim = c(I, J, K))
  ld <- matrix(0, J, K)
  for (j in seq_len(J)) for (k in seq_len(K)) {
    logdz <- log10(N[j, k]) - base_ld                        # centred log depth
    theta <- plogis(ct_i + s_theta * logdz)                  # depth-driven capture
    u     <- rbinom(I, 1, z_true[, j] * theta)               # captured sequences
    r     <- rgamma(I, shape = phi_i, rate = 1) * u           # relative dominance
    if (sum(r) > 0) y[, j, k] <- as.integer(rmultinom(1, N[j, k], r / sum(r)))
    ld[j, k] <- log10(N[j, k])                               # library depth covariate
  }
  ld_std <- (ld - mean(ld)) / sd(ld)                          # z-scored log depth
  det <- apply(y > 0, c(1, 2), any)                           # observed occupancy
  sens <- mean(det[z_true == 1])                              # realised detection|present
  dat  <- occumbData(y = y, site_cov = list(x = x_design),
                     repl_cov = list(logdepth = ld_std))
  list(dat = dat, psi_true = psi_true, sensitivity = sens)
}

fit_quiet <- function(theta_safeguard, data) {
  args <- list(formula_psi = ~ x, data = data,
               prior_ulim = PRIOR_ULIM, prior_prec = PRIOR_PREC,
               n.chains = N_CHAINS, n.adapt = N_ADAPT, n.burnin = N_BURN,
               n.thin = N_THIN, n.iter = N_ITER, verbose = FALSE)
  if (theta_safeguard) args$formula_theta_shared <- ~ logdepth  # (A) depth-aware
  suppressWarnings(do.call(occumb, args))                       # (B) depth-blind default
}

# ---- run the study ----------------------------------------------------------
message(sprintf("[10_sim] %d sims x 2 fits | I=%d J=%d K=%d | MCMC %d chains x %d iter",
                NSIM, I, J, K, N_CHAINS, N_ITER))
t_start <- Sys.time()
cells   <- list()      # per-(sim,model,species,unit) records
sens_v  <- numeric(0)  # realised detection sensitivity per sim
rhat_v  <- numeric(0)  # max psi Rhat per fit (convergence audit)

for (s in seq_len(NSIM)) {
  set.seed(1000 + s)
  sim <- simulate_one()
  sens_v <- c(sens_v, sim$sensitivity)

  for (mod in c("safeguard", "naive")) {
    fit <- fit_quiet(theta_safeguard = (mod == "safeguard"), data = sim$dat)
    ps  <- get_post_samples(fit, "psi")                # draws x I x J
    pm  <- apply(ps, c(2, 3), mean)
    lo  <- apply(ps, c(2, 3), quantile, probs = 0.025)
    hi  <- apply(ps, c(2, 3), quantile, probs = 0.975)
    rhat_v <- c(rhat_v, max(get_post_summary(fit, "psi")[, "Rhat"], na.rm = TRUE))
    for (i in seq_len(I)) for (j in seq_len(J)) {
      pt <- sim$psi_true[i, j]
      cells[[length(cells) + 1]] <- data.frame(
        sim = s, model = mod, species = i, unit = j,
        psi_true = pt, psi_hat = pm[i, j],
        lo95 = lo[i, j], hi95 = hi[i, j],
        covered = (pt >= lo[i, j] & pt <= hi[i, j]))
    }
  }
  message(sprintf("  sim %2d/%d done  (detection|present = %.2f)", s, NSIM, sim$sensitivity))
}

rec <- do.call(rbind, cells)

# ---- aggregate: per true-psi bin x model ------------------------------------
brks <- c(0, 0.15, 0.35, 0.55, 0.75, 1.0)
labs <- c("[0,0.15]", "(0.15,0.35]", "(0.35,0.55]", "(0.55,0.75]", "(0.75,1]")
rec$psi_bin <- factor(labs[findInterval(rec$psi_true, brks, rightmost.closed = TRUE)],
                      levels = labs)

agg_list <- list()
for (m in c("safeguard", "naive")) for (b in labs) {
  sub <- rec[rec$model == m & rec$psi_bin == b, ]
  if (nrow(sub) == 0) next
  agg_list[[length(agg_list) + 1]] <- data.frame(
    model         = m,
    psi_bin       = b,
    psi_true_mean = mean(sub$psi_true),
    n_cells       = nrow(sub),
    bias          = mean(sub$psi_hat - sub$psi_true),
    rmse          = sqrt(mean((sub$psi_hat - sub$psi_true)^2)),
    coverage95    = mean(sub$covered))
}
agg <- do.call(rbind, agg_list)
agg <- agg[order(match(agg$psi_bin, labs), agg$model), ]

# ---- hard convergence gate --------------------------------------------------
rhat_max <- max(rhat_v)
rhat_median <- median(rhat_v)
rhat_threshold <- 1.05
convergence_pass <- is.finite(rhat_max) && rhat_max <= rhat_threshold
claim_status <- if (convergence_pass) "NUMERIC_DIAGNOSTIC_AVAILABLE" else "WITHHELD_UNCONVERGED"
agg$claim_status <- claim_status

status <- data.frame(
  analysis = "Layer A smoke-scale known-truth recovery study",
  fit_scale = sprintf("%d simulations x 2 fits; %d chains x %d iterations", NSIM, N_CHAINS, N_ITER),
  max_psi_Rhat = rhat_max,
  median_fit_max_psi_Rhat = rhat_median,
  convergence_threshold = rhat_threshold,
  convergence_pass = convergence_pass,
  claim_status = claim_status,
  allowed_claim = if (convergence_pass)
    "Numerical recovery diagnostics may be discussed with smoke-scale limitations"
  else
    "No recovery magnitude or direction; current rare-species occupancy claims remain demoted",
  stringsAsFactors = FALSE)

# ---- write TSV --------------------------------------------------------------
dir.create("data/derived", showWarnings = FALSE, recursive = TRUE)
dir.create("data/metadata", showWarnings = FALSE, recursive = TRUE)
dir.create("figures",      showWarnings = FALSE, recursive = TRUE)
out_tsv <- "data/derived/sim_recovery_results.tsv"
write.table(agg, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message("[10_sim] wrote ", out_tsv)
status_tsv <- "data/metadata/layerA_recovery_status.tsv"
write.table(status, status_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message("[10_sim] wrote ", status_tsv)

# A failed recovery gate must suppress every cached Layer-A numerical figure,
# not only the recovery plot below. The helper is also directly runnable when a
# completed status receipt already exists.
source("scripts/apply_layerA_claim_gate.R")
if (!convergence_pass) write_layerA_withheld_figures(status)

# ---- figure: bias / coverage / safeguard shift, vs true psi -----------------
sg <- agg[agg$model == "safeguard", ]
nv <- agg[agg$model == "naive", ]
sg <- sg[order(sg$psi_true_mean), ]; nv <- nv[match(sg$psi_bin, nv$psi_bin), ]
col_sg <- "#1b6ca8"; col_nv <- "#c8102e"

png("figures/phase3_sim_recovery.png", width = 1500, height = 520, res = 130)
if (!convergence_pass) {
  par(mar = c(1, 1, 3, 1))
  plot.new()
  title("Layer A numerical recovery result withheld", cex.main = 1.8)
  text(0.5, 0.60, sprintf("Convergence gate failed: max psi Rhat = %.2f (required <= %.2f)",
                          rhat_max, rhat_threshold), cex = 1.4)
  text(0.5, 0.46, "No recovery magnitude or direction is an allowed finding.", cex = 1.4)
  text(0.5, 0.31, "Headline conclusion: current rare-species occupancy estimates remain exploratory.",
       cex = 1.2, col = "#a33d46")
} else {
  op <- par(mfrow = c(1, 3), mar = c(4.4, 4.4, 3.2, 1.0), oma = c(0, 0, 2.2, 0),
            cex.lab = 1.05)

  ## panel 1: bias vs true psi
  yl <- range(c(sg$bias, nv$bias, 0)) + c(-0.02, 0.02)
  plot(sg$psi_true_mean, sg$bias, type = "b", pch = 19, col = col_sg, lwd = 2,
       ylim = yl, xlim = c(0, 1), xlab = "true psi (bin mean)",
       ylab = "bias  (mean psi_hat - psi_true)", main = "A. Recovery bias")
  abline(h = 0, lty = 2, col = "grey50")
  lines(nv$psi_true_mean, nv$bias, type = "b", pch = 17, col = col_nv, lwd = 2)
  legend("topright", c("safeguard (theta ~ logdepth)", "naive (depth-blind)"),
         col = c(col_sg, col_nv), pch = c(19, 17), lwd = 2, bty = "n", cex = 0.85)

  ## panel 2: coverage vs true psi
  plot(sg$psi_true_mean, sg$coverage95, type = "b", pch = 19, col = col_sg, lwd = 2,
       ylim = c(0, 1), xlim = c(0, 1), xlab = "true psi (bin mean)",
       ylab = "95% CrI coverage", main = "B. Calibration")
  abline(h = 0.95, lty = 2, col = "grey50")
  lines(nv$psi_true_mean, nv$coverage95, type = "b", pch = 17, col = col_nv, lwd = 2)

  ## panel 3: safeguard's psi shift (A - B)
  delta <- sg$bias - nv$bias
  plot(sg$psi_true_mean, delta, type = "b", pch = 18, col = "#5b3a8c", lwd = 2,
       xlim = c(0, 1), xlab = "true psi (bin mean)",
       ylab = "psi shift  (safeguard - naive)", main = "C. Where the safeguard acts")
  abline(h = 0, lty = 2, col = "grey50")

  mtext(sprintf("occumb psi recovery under a depth x ecological-gradient confound  |  %d sims, I=%d J=%d K=%d  |  realised detection|present = %.0f%%",
                NSIM, I, J, K, 100 * mean(sens_v)), outer = TRUE, cex = 0.8)
  par(op)
}
dev.off()
message("[10_sim] wrote figures/phase3_sim_recovery.png")

# ---- self-check + one-line verdict ------------------------------------------
stopifnot(
  "results table is empty"          = nrow(agg) > 0,
  "coverage outside [0,1]"          = all(agg$coverage95 >= 0 & agg$coverage95 <= 1),
  "non-finite bias/rmse"            = all(is.finite(agg$bias)) && all(is.finite(agg$rmse)),
  "rmse negative"                   = all(agg$rmse >= 0)
)

low <- labs[1]                                            # lowest true-psi bin
gv  <- function(m, col, b = low) agg[[col]][agg$model == m & agg$psi_bin == b]
bA  <- gv("safeguard", "bias"); bB <- gv("naive", "bias")
cA  <- gv("safeguard", "coverage95"); cB <- gv("naive", "coverage95")
rA  <- gv("safeguard", "rmse");        rB <- gv("naive", "rmse")
shift_low <- bA - bB                                      # safeguard's low-psi shift
runtime   <- as.numeric(Sys.time() - t_start, units = "mins")

# data-driven, honest classification (no assumed direction)
dir_naive <- if (bB < -0.03) "DOWNWARD (as feared)" else if (bB > 0.03) "UPWARD (over-estimate)" else "negligible"
helps     <- abs(bA) < abs(bB) - 1e-3   # does the safeguard bring low-psi bias closer to 0?
sg_effect <- if (helps) "REDUCES" else "does NOT reduce (it enlarges)"

if (convergence_pass) {
  cat("\n================= SIMULATION RECOVERY: results (per true-psi bin) =================\n")
  print(format(agg, digits = 3), row.names = FALSE)
  cat(sprintf("\nconvergence (psi Rhat over all %d fits): max %.2f | median-of-fit-maxes %.2f\n",
              length(rhat_v), max(rhat_v), median(rhat_v)))
  cat(sprintf("realised detection|present = %.1f%%\n", 100 * mean(sens_v)))
  cat(sprintf(paste0(
    "VERDICT: naive low-psi bias is %s (%+0.3f); the depth safeguard %s it ",
    "(safeguard bias %+0.3f, shift %+0.3f vs naive). 95%% CrI coverage at low psi: ",
    "safeguard %.2f vs naive %.2f. RESIDUAL at the lowest psi: |bias|=%.3f, RMSE=%.3f ",
    "(safeguard). These numerical diagnostics are smoke-scale and are not a biological ",
    "occupancy result.\n"),
    dir_naive, bB, sg_effect, bA, shift_low, cA, cB, abs(bA), rA))
  cat(sprintf("total runtime: %.1f min\n", runtime))
  cat("==================================================================================\n")
} else {
  cat("\n================= SIMULATION RECOVERY: CLAIM WITHHELD =================\n")
  cat(sprintf("status: %s\n", claim_status))
  cat(sprintf("convergence gate: max psi Rhat %.2f; required <= %.2f\n",
              rhat_max, rhat_threshold))
  cat("No recovery magnitude or direction is allowed from these unconverged fits.\n")
  cat("Cached numerical Layer-A displays were replaced by withheld-status cards.\n")
  cat(sprintf("total runtime: %.1f min\n", runtime))
  cat("=======================================================================\n")
}
