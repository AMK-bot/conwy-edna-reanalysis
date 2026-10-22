# shiny/app.R  --  River Conwy fish eDNA: occupancy & reclassification explorer (Phase-4 deliverable).
# Run:  micromamba run -n edna Rscript -e "shiny::runApp('shiny', launch.browser = FALSE, port = 8080)"
# Deploy: shinyapps.io (needs the owner's account — the one user-input step).
#
# Reads local Phase-3 summary tables. Layer-A tables are currently smoke-scale diagnostics,
# not production estimates, and are git-ignored until their provenance is finalized. The app
# therefore labels them explicitly and never upgrades itself based on a guessed/nonexistent file.
# Occurrence framing throughout (reads != abundance). Deps kept to shiny + ggplot2 (both in the env).

library(shiny)
library(ggplot2)

# --- locate data/derived relative to this app, wherever it's launched from ---
find_derived <- function() {
  for (p in c("data/derived", "../data/derived",
              file.path(dirname(getwd()), "data/derived"))) if (dir.exists(p)) return(p)
  "data/derived"
}
DERIVED <- find_derived()
APP_ISSUES <- character()
rd <- function(f, required = character(), base_dir = DERIVED) { # read a TSV or return NULL
  p <- file.path(base_dir, f)
  if (!file.exists(p)) {
    APP_ISSUES <<- c(APP_ISSUES, paste0(f, " is missing"))
    return(NULL)
  }
  out <- tryCatch(
    read.delim(p, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) {
      APP_ISSUES <<- c(APP_ISSUES, paste0(f, " could not be read: ", conditionMessage(e)))
      NULL
    })
  if (is.null(out)) return(NULL)
  missing_cols <- setdiff(required, names(out))
  if (length(missing_cols)) {
    APP_ISSUES <<- c(
      APP_ISSUES,
      paste0(f, " has the wrong schema; missing: ", paste(missing_cols, collapse = ", ")))
    return(NULL)
  }
  out
}

# Layer A has no `occumb_covariate_estimates.tsv` output. The analysis-label summary that
# actually exists is the intercept-model table; covariate smoke outputs have different schemas
# and are loaded separately below. Production output must carry explicit provenance before this
# label changes (see docs/phase3-claim-status.md).
OCC_FILE <- "occumb_intercept_estimates.tsv"
OCC <- rd(OCC_FILE, c(
  "species", "psi_mean", "psi_l95", "psi_u95", "psi_Rhat", "psi_neff",
  "detection_perrep_mean", "detection_perrep_l95", "detection_perrep_u95",
  "false_neg_perrep_mean", "false_neg_3rep_mean"))
MODEL_LABEL <- "intercept-only smoke fit (2 chains x 2,000 iterations; diagnostic only)"

ALPHA <- rd("layerE_alpha_site_time.tsv", c("section", "hill_q1"))
INDIC <- rd("layerE_indicator_species.tsv", c("species", "IndVal_stat"))
DRIFT <- rd("reclassB_drift_catalogue.tsv", c(
  "db_species", "Genus", "audit_category", "audit_hypothesis", "evidence_flag",
  "n_asv", "pct_reads_full", "audit_note"))
RATES <- rd("reclassB_assignment_rates.tsv")
CVR <- rd("occumb_covariate_corrected_vs_raw.tsv", c(
  "species", "raw_naive_occ", "corrected_psi", "corrected_l95", "corrected_u95"))
DOPT <- rd("occumb_covariate_design_opt.tsv", c("K", "N", "Utility", "source"))
SIM <- rd("sim_recovery_results.tsv", c(
  "model", "psi_bin", "psi_true_mean", "bias", "rmse", "coverage95", "claim_status"))
RECOVERY_STATUS <- rd(
  "layerA_recovery_status.tsv",
  c("max_psi_Rhat", "convergence_threshold", "convergence_pass", "claim_status", "allowed_claim"),
  file.path(dirname(DERIVED), "metadata"))

pretty_sp <- function(x) gsub("_", " ", x)

ui <- navbarPage(
  title = "River Conwy fish eDNA — occupancy & reclassification explorer",
  header = tagList(
    tags$div(style = "padding:8px 14px;background:#fff3cd;color:#664d03;font-size:92%;",
      tags$b("Layer A status: exploratory diagnostic, not a final result. "),
      "Production MCMC and the planned depth-threshold/model sensitivity sweep have not run. ",
      "Do not cite the rare-species occupancy values or their correction direction."),
    tags$div(style = "padding:6px 14px;color:#555;font-size:90%;",
      "Reanalysis of Perry et al. 2024 (Nat. Commun. 15:4372). Detections, not abundances. ",
      tags$b("Layer A table shown: "), MODEL_LABEL,
      if (length(APP_ISSUES)) tags$span(
        style = "color:#9c2f2f;", paste0(" | unavailable inputs: ", length(APP_ISSUES))))),

  # ---- Tab 1: occupancy & detection ----
  tabPanel("Occupancy & detection",
    if (is.null(OCC)) h4("occumb estimates not found — run R/06_occupancy_occumb.R first.") else
    sidebarLayout(
      sidebarPanel(width = 3,
        selectInput("sp", "Species", choices = setNames(OCC$species, pretty_sp(OCC$species))),
        htmlOutput("sp_stats"),
        tags$hr(),
        tags$small(tags$b("Diagnostic definitions only. "), tags$b("psi"),
                   " = modelled occupancy probability. ",
                   tags$b("Detection/rep"), " = probability a single water-sample replicate detects the ",
                   "species where present. ", tags$b("False-negative rate"), " = chance of missing it ",
                   "despite presence (per replicate, and over the nominal 3). ",
                   "The values below come from a short smoke fit; rare-taxon values are not inferential.")),
      mainPanel(width = 9,
        fluidRow(
          column(7, h4("Smoke-fit occupancy diagnostic"), plotOutput("psi_forest", height = 460)),
          column(5, h4("Smoke-fit detection diagnostic"), plotOutput("det_scatter", height = 460))))
    )),

  # ---- Tab 2: failed recovery gate + provisional design diagnostic ----
  tabPanel("Detection correction & design",
    if (is.null(CVR)) h4("Corrected-vs-raw output not found — run R/06_occupancy_occumb.R (covariate model).") else
    fluidPage(
      h4("Why the smoke-fit correction is not a final occupancy estimate"),
      tags$div(style = "padding:10px 12px;background:#f8d7da;color:#842029;margin-bottom:10px;",
        tags$b("Recovery gate failed: "),
        "the smoke-scale known-truth MCMC did not converge, so neither its magnitude nor direction is ",
        "reported. Eighteen observed libraries are zero within the 14-label harmonised Layer-A input ",
        "and are dropped by the community model but counted as non-detections by the naïve comparator; ",
        "those libraries contain reads under raw taxonomy labels."),
      p("Blue points show the covariate smoke fit only to document why the claim was demoted; grey points are naïve ",
        "detection. Treat all rare-taxon smoke-fit adjustments as uninterpretable diagnostics, not discoveries. ",
        "Longer production chains can address Monte Carlo convergence but cannot, by themselves, solve this ",
        "structural recoverability or comparison-mask problem."),
      plotOutput("corrected_plot", height = 470),
      if (!is.null(RECOVERY_STATUS)) tagList(
        h4("Recovery-study gate receipt"),
        tableOutput("recovery_status_tbl"),
        tags$small("Numerical recovery rows are suppressed whenever convergence fails.")),
      tags$hr(),
      h4("Survey-design diagnostic (provisional)"),
      p("The current smoke-fit utility calculation favours more water-sample replicates over deeper ",
        "sequencing. This is a hypothesis for production/sensitivity testing, not yet a recommendation."),
      if (!is.null(DOPT)) plotOutput("design_plot", height = 300)
    )),

  # ---- Tab 3: community baseline (Layer E) ----
  tabPanel("Community (Layer E)",
    if (is.null(ALPHA)) h4("Layer E outputs not found — run R/04_baseline_E.R first.") else
    fluidRow(
      column(6, h4("Site-time alpha diversity (Hill q1) by river section"),
             plotOutput("alpha_box", height = 420),
             tags$small("Available replicate libraries are combined within site-time units; effort remains 1-3 libraries per unit.")),
      column(6, h4("Exploratory indicators by section"),
             tableOutput("indic_tbl"),
             tags$small("Restricted-series screen on 14 sites; no indicator reaches p <= 0.05.")))),

  # ---- Tab 3: reclassification (Layer B) ----
  tabPanel("Reclassification (Layer B)",
    if (is.null(DRIFT)) h4("Layer B outputs not found — run R/07_reclass_B.R first.") else
    fluidPage(
      h4("Reference-label assignment audit"),
      p("Three comparisons spanning two current reference resources disagree on the species label of the ",
        "dominant short ", tags$i("Salmo"), " ASV. The stable inference is genus; none of the species ",
        "labels is treated as a proven local identity."),
      tableOutput("drift_tbl"),
      if (!is.null(RATES)) tagList(tags$hr(), h4("Assignment rate by rank"), tableOutput("rates_tbl"))
    ))
)

server <- function(input, output, session) {
  sel <- reactive({
    req(!is.null(OCC), input$sp %in% OCC$species)
    OCC[OCC$species == input$sp, , drop = FALSE]
  })

  output$sp_stats <- renderUI({
    s <- sel(); if (nrow(s) == 0) return(NULL)
    fmt <- function(m, l, u) sprintf("%.2f  [%.2f-%.2f]", m, l, u)
    HTML(paste0(
      "<h4 style='margin-top:8px'>", pretty_sp(s$species), "</h4>",
      "<p><b>Occupancy psi:</b> ", fmt(s$psi_mean, s$psi_l95, s$psi_u95), "<br>",
      "<b>Detection / replicate:</b> ", fmt(s$detection_perrep_mean, s$detection_perrep_l95, s$detection_perrep_u95), "<br>",
      "<b>False-negative / replicate:</b> ", sprintf("%.2f", s$false_neg_perrep_mean), "<br>",
      "<b>False-negative over 3 reps:</b> ", sprintf("%.3f", s$false_neg_3rep_mean), "</p>",
      "<p style='color:#842029;font-size:85%'><b>Smoke-fit diagnostic only.</b> Do not report this ",
      "rare-species point estimate as corrected occupancy.</p>",
      paste0("<p style='color:#888;font-size:85%'>convergence Rhat ",
             sprintf("%.3f", s$psi_Rhat), ", n.eff ", s$psi_neff, "</p>")))
  })

  output$psi_forest <- renderPlot({
    req(!is.null(OCC), input$sp %in% OCC$species)
    d <- OCC[order(OCC$psi_mean), ]; d$species <- factor(pretty_sp(d$species), levels = pretty_sp(d$species))
    d$hl <- OCC$species[order(OCC$psi_mean)] == input$sp
    ggplot(d, aes(psi_mean, species, colour = hl)) +
      geom_errorbar(aes(xmin = psi_l95, xmax = psi_u95), orientation = "y", width = 0.3) +
      geom_point(size = 2.6) +
      scale_colour_manual(values = c(`FALSE` = "#4C6EF5", `TRUE` = "#E8590C"), guide = "none") +
      scale_x_continuous(limits = c(0, 1)) +
      labs(x = "Smoke-fit psi (95% CrI; diagnostic only)", y = NULL,
           caption = "Rare-species magnitude and direction are not reportable; see recovery-study audit.") +
      theme_minimal(base_size = 13)
  })

  output$det_scatter <- renderPlot({
    req(!is.null(OCC), input$sp %in% OCC$species)
    d <- OCC; d$hl <- d$species == input$sp
    ggplot(d, aes(psi_mean, detection_perrep_mean, colour = hl)) +
      geom_point(size = 3) +
      scale_colour_manual(values = c(`FALSE` = "#4C6EF5", `TRUE` = "#E8590C"), guide = "none") +
      scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) +
      labs(x = "Smoke-fit psi", y = "Smoke-fit detection / replicate",
           subtitle = "Exploratory model output; not evidence that a reported absence is a false negative") +
      theme_minimal(base_size = 13)
  })

  output$alpha_box <- renderPlot({
    req(!is.null(ALPHA))
    d <- ALPHA; d$section <- factor(d$section, levels = c("upper", "middle", "lower"))
    ggplot(d, aes(section, hill_q1, fill = section)) +
      geom_boxplot(alpha = 0.85, outlier.size = 0.6) +
      scale_fill_manual(values = c(upper = "#63A6FF", middle = "#4C6EF5", lower = "#3B2F9E"), guide = "none") +
      labs(x = NULL, y = "Hill q1 (effective # species)") + theme_minimal(base_size = 13)
  })

  output$indic_tbl <- renderTable({
    if (is.null(INDIC)) return(NULL)
    d <- INDIC[order(-INDIC$IndVal_stat), ]; d$species <- pretty_sp(d$species); d
  }, digits = 3)

  output$drift_tbl <- renderTable({
    req(!is.null(DRIFT))
    cols <- intersect(c("db_species", "audit_hypothesis", "audit_category", "evidence_flag",
                        "n_asv", "pct_reads_full", "Genus"), names(DRIFT))
    d <- DRIFT[order(-DRIFT$pct_reads_full), cols]
    names(d)[names(d) == "pct_reads_full"] <- "% reads"
    head(d, 40)
  }, digits = 2)

  output$rates_tbl <- renderTable(RATES, digits = 1)

  output$corrected_plot <- renderPlot({
    req(!is.null(CVR))
    dd <- CVR[order(CVR$corrected_psi), ]
    dd$species <- factor(pretty_sp(dd$species), levels = pretty_sp(dd$species))
    ggplot(dd) +
      geom_segment(aes(y = species, yend = species, x = raw_naive_occ, xend = corrected_psi), colour = "#CED4DA") +
      geom_errorbar(aes(y = species, xmin = corrected_l95, xmax = corrected_u95),
                    orientation = "y", width = 0.25, colour = "#4C6EF5") +
      geom_point(aes(raw_naive_occ, species), colour = "#868E96", size = 2.4) +
      geom_point(aes(corrected_psi, species), colour = "#3B2F9E", size = 2.9) +
      scale_x_continuous(limits = c(0, 1)) +
      labs(x = "diagnostic only  (grey = naïve detection; blue = confounded smoke-fit psi)", y = NULL,
           caption = "Recovery MCMC failed convergence; blue magnitude and direction are not results.") +
      theme_minimal(base_size = 13)
  })

  output$recovery_status_tbl <- renderTable({
    req(!is.null(RECOVERY_STATUS))
    RECOVERY_STATUS[, c("max_psi_Rhat", "convergence_threshold", "convergence_pass",
                        "claim_status", "allowed_claim")]
  }, digits = 3)

  output$design_plot <- renderPlot({
    req(!is.null(DOPT))
    b <- DOPT[grepl("budget", DOPT$source), ]
    req(nrow(b) > 0)
    ggplot(b, aes(factor(K), Utility)) +
      geom_col(fill = "#4C6EF5", alpha = 0.85) +
      labs(x = "replicates K (at fixed budget)", y = "smoke-fit expected utility",
           caption = "Provisional diagnostic; production and sensitivity runs are still required.") +
      theme_minimal(base_size = 13)
  })
}

# Headless, non-launching schema/reactive smoke check:
#   EDNA_APP_SCHEMA_CHECK=1 micromamba run -n edna Rscript shiny/app.R
if (identical(Sys.getenv("EDNA_APP_SCHEMA_CHECK"), "1")) {
  layer_a <- list(OCC = OCC, CVR = CVR, DOPT = DOPT, SIM = SIM,
                  RECOVERY_STATUS = RECOVERY_STATUS)
  app_data <- c(layer_a, list(DRIFT = DRIFT, RATES = RATES))
  missing_app_data <- names(app_data)[vapply(app_data, is.null, logical(1))]
  if (length(missing_app_data)) {
    stop("Shiny schema check failed: ", paste(missing_app_data, collapse = ", "),
         ". Details: ", paste(APP_ISSUES, collapse = " | "))
  }
  stopifnot(nrow(OCC) > 0, nrow(CVR) > 0, nrow(DOPT) > 0, nrow(SIM) > 0,
            nrow(RECOVERY_STATUS) == 1, nrow(DRIFT) > 0, nrow(RATES) > 0,
            all(OCC$psi_mean >= 0 & OCC$psi_mean <= 1),
            all(CVR$corrected_psi >= 0 & CVR$corrected_psi <= 1))
  ui_html <- htmltools::renderTags(ui)$html
  stopifnot(grepl("Layer A status: exploratory diagnostic", ui_html, fixed = TRUE),
            grepl("Recovery gate failed", ui_html, fixed = TRUE))
  shiny::testServer(server, {
    session$setInputs(sp = OCC$species[[1]])
    stopifnot(nrow(sel()) == 1L)
  })
  cat("PASS: Shiny data schemas, ranges, UI construction, and species reactive are valid.\n")
  invisible(TRUE)
} else {
  shinyApp(ui, server)
}
