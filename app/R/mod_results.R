# mod_results.R : base-case results and export ----------------------------------

mod_results_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("boxes")),
    bslib::layout_columns(
      col_widths = c(7, 5),
      bslib::card(bslib::card_header("Results by arm"), shiny::uiOutput(ns("table"))),
      bslib::card(bslib::card_header("Interpretation"), shiny::uiOutput(ns("text")))
    ),
    bslib::layout_columns(
      col_widths = c(6, 6),
      bslib::card(bslib::card_header("State occupancy: comparator"), shiny::plotOutput(ns("occ_comp"), height = 330)),
      bslib::card(bslib::card_header("State occupancy: intervention"), shiny::plotOutput(ns("occ_int"), height = 330))
    ),
    bslib::card(
      bslib::card_header("Model trace (first 26 cycles)"),
      shiny::radioButtons(ns("trace_arm"), NULL, c("Comparator" = "comp", "Intervention" = "int"), inline = TRUE),
      shiny::uiOutput(ns("trace"))
    )
  )
}

mod_results_server <- function(id, model) {
  shiny::moduleServer(id, function(input, output, session) {
    output$boxes <- shiny::renderUI({
      m <- model()
      r <- m$res
      cur <- m$spec$settings$currency
      bslib::layout_columns(
        col_widths = c(3, 3, 3, 3),
        bslib::value_box("Incremental cost", paste(cur, fmt_num(r$inc_cost, 0))),
        bslib::value_box("Incremental QALYs", fmt_num(r$inc_qaly, 3)),
        bslib::value_box("ICER", if (is.finite(r$icer)) paste(cur, fmt_num(r$icer, 0)) else "n/a",
                         shiny::p(class = "small", r$icer_label)),
        bslib::value_box("Net monetary benefit",
                         paste(cur, fmt_num(r$nmb, 0)),
                         shiny::p(class = "small", paste0("at ", fmt_num(m$spec$settings$wtp, 0), " per QALY")))
      )
    })

    output$table <- shiny::renderUI({
      m <- model()
      html_table(results_display(m$res, m$spec$arms, m$spec$settings$currency))
    })

    output$text <- shiny::renderUI({
      m <- model()
      r <- m$res
      st <- m$spec$settings
      v <- m$validation
      n_fail <- sum(v$status == "fail")
      n_warn <- sum(v$status == "warn")
      ce <- is.finite(r$icer) && r$inc_qaly > 0 && r$icer <= st$wtp || grepl("Dominant", r$icer_label)
      shiny::tagList(
        shiny::p(sprintf("%s gains %s discounted QALYs at an extra cost of %s %s versus %s.",
                         m$spec$arms$int, fmt_num(r$inc_qaly, 3), st$currency,
                         fmt_num(r$inc_cost, 0), m$spec$arms$comp)),
        shiny::p(sprintf("At a threshold of %s %s per QALY the intervention %s cost-effective.",
                         st$currency, fmt_num(st$wtp, 0), if (ce) "is" else "is not")),
        if (is.finite(r$share_ly_after_fu)) shiny::p(sprintf(
          "%s of the incremental life-years accrue after the end of trial follow-up (year %s), so they rest on extrapolation.",
          fmt_pct(r$share_ly_after_fu), round(st$fu_years, 2))),
        shiny::p(class = if (n_fail) "text-danger" else if (n_warn) "text-warning" else "text-success",
                 sprintf("Validation: %d failed, %d warnings. See the Validation tab.", n_fail, n_warn))
      )
    })

    output$occ_comp <- shiny::renderPlot({
      m <- model()
      plot_occupancy(m$curves, "comp", m$spec$arms$comp) + ggplot2::coord_cartesian(xlim = c(0, 15))
    })
    output$occ_int <- shiny::renderPlot({
      m <- model()
      plot_occupancy(m$curves, "int", m$spec$arms$int) + ggplot2::coord_cartesian(xlim = c(0, 15))
    })

    output$trace <- shiny::renderUI({
      m <- model()
      tr <- utils::head(m$res[[input$trace_arm]]$trace, 26)
      html_table(tr[, c("cycle", "t_start", "pf", "pd", "dead", "on_tx_prop", "ly", "qaly",
                        "drug", "pf_state", "pd_state", "subsequent", "terminal")], digits = 4)
    })
  })
}

# ---- export -------------------------------------------------------------------

mod_export_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::card(
    bslib::card_header("Export"),
    shiny::p("One Excel workbook with inputs, results, full traces for both arms, curves,",
             "goodness of fit, validation checks, saved scenarios, DSA and PSA runs."),
    shiny::downloadButton(ns("xlsx"), "Download Excel workbook", class = "btn-primary"),
    shiny::p(class = "small text-muted mt-3",
             "Every number can be reproduced outside Shiny: see run_model() in R/fct_spec.R and the README.")
  )
}

mod_export_server <- function(id, model, fits, scenarios, dsa, psa,
                              psa_bands = \() NULL, psa_compare = \() NULL,
                              waning_explorer = \() NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    output$xlsx <- shiny::downloadHandler(
      filename = \() paste0("psm_results_", format(Sys.time(), "%Y%m%d_%H%M"), ".xlsx"),
      content = \(file) {
        m <- model()
        sp <- m$spec
        sp$settings$lifetable <- NULL
        cv <- m$curves
        sheets <- list(
          README = data.frame(note = c(
            "Partitioned survival model export.",
            "Inputs: model specification as used. Results: discounted unless stated.",
            "Trace sheets: per-cycle occupancy, life-years, QALYs and costs (undiscounted per cycle; discount factors given).",
            "Curves: final model survival after effects, waning, general-population floor and constraints.",
            "PSA ICER is the ratio of mean incremental cost to mean incremental QALYs.",
            "PSA_runs: one row per run; invalid runs keep valid = FALSE and a reason (excluded from summaries, not fixed).")),
          Inputs = flatten_spec(sp),
          Results = m$res$table,
          Trace_comparator = m$res$comp$trace,
          Trace_intervention = m$res$int$trace,
          Curves = data.frame(time = cv$t,
                              comp_OS = cv$comp$OS, comp_PFS = cv$comp$PFS, comp_ToT = cv$comp$ToT,
                              int_OS = cv$int$OS, int_PFS = cv$int$PFS, int_ToT = cv$int$ToT,
                              gen_pop = cv$gp %||% NA),
          GOF_OS = fits()$OS$gof,
          GOF_PFS = fits()$PFS$gof,
          Validation = m$validation
        )
        if (!is.null(m$spec$settings$lifetable)) sheets$Life_table <- m$spec$settings$lifetable
        sc <- scenarios()
        if (!is.null(sc) && nrow(sc)) sheets$Scenarios <- sc
        d <- dsa()
        if (!is.null(d)) sheets$DSA <- as.data.frame(d)
        p <- psa()
        if (!is.null(p)) {
          st <- attr(p, "settings")
          sheets$PSA_settings <- data.frame(
            setting = c("runs", "seed", "survival uncertainty", "paired by patient id",
                        "curves sampled", "cost CV", "OS-PFS log-HR correlation", "workers"),
            value = c(st$n, st$seed, st$surv_method, st$paired,
                      paste(st$sampled_curves, collapse = "; "), st$cost_cv, st$hr_corr, st$cores))
          attr(p, "bands") <- NULL
          attr(p, "settings") <- NULL
          sheets$PSA_runs <- p
        }
        b <- psa_bands()
        if (!is.null(b)) sheets$PSA_survival_bands <- b
        cmp <- psa_compare()
        if (!is.null(cmp)) sheets$PSA_vs_deterministic <- cmp
        wx <- waning_explorer()
        if (!is.null(wx)) {
          wx$matches_current_inputs <- identical(attr(wx, "spec"), explorer_key(m$spec))
          attr(wx, "spec") <- NULL
          sheets$Waning_explorer <- wx
        }
        wd <- waning_diagnostics(m$curves, m$spec)
        sheets$Waning <- data.frame(
          item = c("structure", "anchor", "start (model year)", "end (model year)", "target HR",
                   "frailty variance", "applied to OS", "applied to PFS", "justification"),
          value = c(m$spec$waning$type, m$spec$waning$anchor, wd$start, wd$end, m$spec$waning$target_hr,
                    m$spec$waning$frailty, m$spec$waning$apply_os, m$spec$waning$apply_pfs,
                    m$spec$waning$rationale))
        writexl::write_xlsx(sheets, file)
      }
    )
  })
}
