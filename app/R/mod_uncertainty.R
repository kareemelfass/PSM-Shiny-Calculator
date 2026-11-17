# mod_uncertainty.R : scenarios, one-way DSA, PSA --------------------------------

mod_scenarios_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    intro_card(
      purpose = "Test how much the conclusion depends on your assumptions, one change at a time.",
      do = "Scenarios: change any input on another tab, give the scenario a name, and click Save current result. The first one you save is the base case. DSA: click Run DSA to move each input across a plausible range and rank them in a tornado chart.",
      check = "Which inputs move the result most, and whether any scenario flips the decision at your willingness-to-pay threshold.",
      next_step = "the PSA tab."
    ),
    bslib::card(
      bslib::card_header("Scenario analysis"),
      bslib::layout_columns(
        col_widths = c(6, 3, 3),
        shiny::textInput(ns("name"), NULL, placeholder = "Scenario name, e.g. 'Log-normal OS, waning yr 5'"),
        shiny::actionButton(ns("save"), "Save current result", class = "btn-primary"),
        shiny::actionButton(ns("clear"), "Clear all")
      ),
      shiny::helpText("Change any input on the other tabs, then save. The first saved scenario is the base case.",
                      "A scenario matters when it changes the conclusion at the threshold (Doc 13)."),
      shiny::uiOutput(ns("table"))
    ),
    bslib::card(
      bslib::card_header("One-way deterministic sensitivity analysis"),
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        shiny::div(
          shiny::sliderInput(ns("cost_pct"), help_label("Cost range (+/- %)", "Every cost input is moved up and down by this percentage."), 5, 50, 20),
          shiny::sliderInput(ns("wane_delta"), help_label("Waning timing range (+/- years)", "Years added to and taken from the waning start time. Only matters if waning is switched on."), 0, 5, 2, step = 0.5)
        ),
        shiny::radioButtons(ns("outcome"), help_label("Outcome", "ICER: the cost per QALY gained. NMB: net monetary benefit, which is QALYs times willingness to pay, minus cost. NMB is easier to read when the ICER jumps around."), c("ICER" = "icer", "Net monetary benefit" = "nmb"), inline = TRUE),
        shiny::actionButton(ns("run_dsa"), "Run DSA", class = "btn-primary")
      ),
      shiny::helpText("HRs use their 95% interval; utilities use +/- 1.96 SE; discount rates 0% to twice the base.",
                      "If waning is on, its start (and duration, if gradual) move by the timing range.",
                      "Survival distributions and waning structures: test as scenarios or in the waning explorer."),
      shiny::plotOutput(ns("tornado"), height = 420),
      shiny::uiOutput(ns("dsa_table"))
    )
  )
}

mod_scenarios_server <- function(id, model) {
  shiny::moduleServer(id, function(input, output, session) {
    store <- shiny::reactiveVal(NULL)
    shiny::observeEvent(input$save, {
      m <- model()
      r <- m$res
      nm <- if (nzchar(trimws(input$name))) trimws(input$name) else paste("Scenario", NROW(store()) + 1)
      row <- data.frame(scenario = nm, inc_cost = r$inc_cost, inc_qaly = r$inc_qaly,
                        icer = r$icer, label = r$icer_label, nmb = r$nmb,
                        os_comp = m$spec$endpoints$OS$comp_model, os_approach = m$spec$endpoints$OS$approach,
                        waning = if (m$spec$waning$type == "none") "none" else paste0(m$spec$waning$type, " @", m$spec$waning$start, " (", m$spec$waning$anchor, ")"),
                        wtp = m$spec$settings$wtp)
      store(rbind(store(), row))
      shiny::updateTextInput(session, "name", value = "")
    })
    shiny::observeEvent(input$clear, store(NULL))

    output$table <- shiny::renderUI({
      s <- store()
      if (is.null(s)) return(shiny::p(class = "text-muted", "No scenarios saved yet."))
      base <- s[1, ]
      ce <- \(r) (r$inc_qaly > 0 & is.finite(r$icer) & r$icer <= r$wtp) | grepl("Dominant", r$label)
      show <- data.frame(
        Scenario = s$scenario,
        `Inc. cost` = fmt_num(s$inc_cost, 0), `Inc. QALYs` = fmt_num(s$inc_qaly, 3),
        ICER = ifelse(is.finite(s$icer), fmt_num(s$icer, 0), "n/a"),
        `Change vs base` = ifelse(is.finite(s$icer) & is.finite(base$icer),
                                  fmt_pct(s$icer / base$icer - 1), "n/a"),
        `Cost-effective?` = ifelse(ce(s), "Yes", "No"),
        status = ifelse(ce(s) == ce(base), "pass", "warn"),
        check.names = FALSE)
      show$status[1] <- "info"
      names(show)[names(show) == "status"] <- "Conclusion vs base"
      shiny::tagList(
        html_table(show, status_col = "Conclusion vs base"),
        shiny::helpText("warn = this scenario changes the cost-effectiveness conclusion. With very small",
                        "or negative incremental QALYs, read costs and QALYs rather than the ICER change.")
      )
    })

    dsa <- shiny::eventReactive(input$run_dsa, {
      m <- model()
      shiny::withProgress(message = "Running one-way sensitivity analysis",
                          run_dsa(m$pieces, m$spec, input$cost_pct / 100, input$wane_delta))
    })

    output$tornado <- shiny::renderPlot({
      d <- dsa()
      base <- attr(d, "base")
      plot_tornado(d, input$outcome, if (input$outcome == "icer") base$icer else base$nmb)
    })
    output$dsa_table <- shiny::renderUI({
      d <- dsa()
      html_table(data.frame(Parameter = d$parameter, Low = signif(d$low_value, 4), High = signif(d$high_value, 4),
                            `ICER low` = fmt_num(d$icer_low, 0), `ICER high` = fmt_num(d$icer_high, 0),
                            `NMB low` = fmt_num(d$nmb_low, 0), `NMB high` = fmt_num(d$nmb_high, 0),
                            check.names = FALSE))
    })

    list(scenarios = store, dsa = shiny::reactive(if (input$run_dsa > 0) dsa() else NULL))
  })
}

mod_psa_ui <- function(id) {
  ns <- shiny::NS(id)
  max_cores <- max(1, parallel::detectCores(logical = FALSE) %||% 1)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 330,
      shiny::radioButtons(ns("surv"), help_label("Survival curve uncertainty", "Parameter covariance: draws new curve parameters from the uncertainty of the fitted models. Fast. Bootstrap: resamples the patients and refits the models in every run. Slower, and it keeps each patient's OS, PFS and ToT rows together when the data have an id column. None: the curves stay fixed, so only the other inputs vary."), c(
        "Parameter covariance (multivariate normal)" = "mvn",
        "Bootstrap patients and refit" = "boot",
        "None: curves fixed (V1)" = "none")),
      shiny::uiOutput(ns("surv_note")),
      shiny::numericInput(ns("n"), help_label("Number of runs", "How many random draws are made. More runs give a steadier answer but take longer. In the browser, start with 100 to 200."), 1000, min = 50, max = 10000, step = 100),
      shiny::numericInput(ns("seed"), help_label("Random seed", "Fixes the random numbers, so the same settings always give the same results. Change it to check that the answer is stable."), 2026),
      shiny::numericInput(ns("cores"), help_label(paste0("Parallel workers (1 to ", max_cores, ")"), "How many processors work at the same time. This speeds up big runs on your own computer. The browser version runs on one."), 1,
                          min = 1, max = max_cores, step = 1),
      shiny::sliderInput(ns("cost_cv"), help_label("Cost uncertainty (coefficient of variation, %)", "How much costs vary between draws, as a percentage of their value. A larger value gives a wider spread."), 5, 50, 20),
      shiny::sliderInput(ns("corr"), help_label("Correlation of OS and PFS log-HRs", "How strongly the uncertainty in the OS and PFS hazard ratios moves together. 0 means they vary independently."), -0.9, 0.9, 0, step = 0.1),
      shiny::checkboxInput(ns("wane_on"), "Sample the waning start time", FALSE),
      shiny::conditionalPanel(
        sprintf("input['%s']", ns("wane_on")),
        shiny::sliderInput(ns("wane_rng"), "Waning start, uniform range (same anchor as the Treatment effects tab)",
                           0, 15, c(3, 7), step = 0.5)
      ),
      shiny::actionButton(ns("run"), "Run PSA", class = "btn-primary"),
      shiny::helpText("Also sampled: entered log-HRs (normal), utilities (beta), costs (gamma).",
                      "The app is busy while the PSA runs.")
    ),
    intro_card(
      purpose = "Measure how uncertain the result is when all the uncertain inputs vary together.",
      do = "Choose how to include survival-curve uncertainty, set the number of runs, then click Run PSA. Start with a small number of runs as a test, then increase it.",
      check = "The spread on the cost-effectiveness plane and the acceptability curve. Also that the convergence plot has flattened, which means more runs would not change the answer, and that the invalid-runs table at the bottom is understood.",
      next_step = "the Export tab."
    ),
    shiny::uiOutput(ns("summary")),
    bslib::layout_columns(
      col_widths = c(6, 6),
      bslib::card(shiny::plotOutput(ns("plane"), height = 380)),
      bslib::card(shiny::plotOutput(ns("ceac"), height = 380))
    ),
    bslib::card(
      bslib::card_header("Survival uncertainty"),
      shiny::numericInput(ns("band_x"), "Plot to year", 15, min = 1, width = "160px"),
      shiny::plotOutput(ns("bands"), height = 380)
    ),
    bslib::layout_columns(
      col_widths = c(6, 6),
      bslib::card(
        bslib::card_header("Deterministic vs probabilistic"),
        shiny::uiOutput(ns("compare")),
        shiny::helpText("Means differ when the model is non-linear in its inputs (expected with",
                        "survival curves). A large gap is worth explaining, not hiding.")
      ),
      bslib::card(
        bslib::card_header("Convergence"),
        shiny::radioButtons(ns("conv_out"), NULL, c("ICER" = "icer", "Net monetary benefit" = "nmb"), inline = TRUE),
        shiny::plotOutput(ns("conv"), height = 300)
      )
    ),
    shiny::uiOutput(ns("invalid"))
  )
}

mod_psa_server <- function(id, model, fits, ipd, code_dir) {
  shiny::moduleServer(id, function(input, output, session) {
    output$surv_note <- shiny::renderUI({
      m <- model()
      tasks <- surv_tasks(m$spec)
      what <- if (length(tasks)) paste(vapply(tasks, \(x) paste(x$ep, paste(m$spec$arms[x$arms], collapse = " + ")), ""), collapse = "; ") else "none (all curves come from HRs)"
      d <- ipd()
      paired <- ipd_is_paired(d[d$arm %in% unlist(m$spec$arms), ])
      shiny::div(class = "small text-muted mb-2",
        shiny::p(shiny::strong("Fitted curves sampled: "), what),
        if (input$surv == "mvn") shiny::p("Fast. Draws each model's parameters from their estimated covariance; arms of a joint model move together. OS, PFS and ToT models are sampled independently.") else
        if (input$surv == "boot") shiny::p(paste("Resamples patients and refits the chosen models each run (slower; consider 500 runs and parallel workers).",
          if (paired) "Patient ids found: OS, PFS and ToT stay linked." else "No patient ids: endpoints are resampled separately.")) else
        shiny::p("Survival curves are the same in every run (V1 behaviour). Understates uncertainty.")
      )
    })

    runs <- shiny::eventReactive(input$run, {
      m <- model()
      shiny::withProgress(message = "Running PSA", value = 0, {
        run_psa(m$pieces, m$spec, n = input$n, seed = input$seed, cost_cv = input$cost_cv / 100,
                hr_corr = input$corr, surv_method = input$surv, fits = fits(), ipd = ipd(),
                cores = min(max(1, input$cores %||% 1), parallel::detectCores() %||% 1),
                code_dir = code_dir,
                wane_start_range = if (isTRUE(input$wane_on) && m$spec$waning$type != "none") input$wane_rng else NULL,
                progress = \(p) shiny::setProgress(p, detail = sprintf("%d of %d runs", round(p * input$n), input$n)))
      })
    })
    wtp <- shiny::reactive(model()$spec$settings$wtp)
    summ <- shiny::reactive(psa_summary(runs(), wtp()))

    output$summary <- shiny::renderUI({
      s <- summ()
      cur <- model()$spec$settings$currency
      method_lab <- c(none = "curves fixed", mvn = "parameter covariance", boot = "bootstrap")[[s$settings$surv_method]]
      shiny::tagList(
        bslib::layout_columns(
          col_widths = c(3, 3, 3, 3),
          bslib::value_box("Mean incremental cost", paste(cur, fmt_num(s$mean_d_cost, 0))),
          bslib::value_box("Mean incremental QALYs", fmt_num(s$mean_d_qaly, 3)),
          bslib::value_box("Probabilistic ICER", paste(cur, fmt_num(s$icer, 0)),
                           shiny::p(class = "small", "Ratio of means (Doc 13 Case 17)")),
          bslib::value_box("Probability cost-effective", fmt_pct(s$p_ce),
                           shiny::p(class = "small", paste("at", fmt_num(wtp(), 0))))
        ),
        shiny::p(class = if (s$n_invalid) "text-danger" else "text-success",
                 sprintf("%d of %d runs valid (survival never rises, occupancy within 0-1 and sums to 1). Survival uncertainty: %s.",
                         s$n_valid, s$n, method_lab)),
        if (!is.null(s$settings$wane_start_range)) shiny::p(class = "small text-muted",
          sprintf("Waning start sampled uniformly between %s and %s years (structure and duration fixed).",
                  s$settings$wane_start_range[1], s$settings$wane_start_range[2])),
        if (s$n_constraint) shiny::p(class = "text-warning",
          sprintf("PFS exceeded OS before the explicit constraint in %d runs (%s).", s$n_constraint, fmt_pct(s$n_constraint / s$n)))
      )
    })
    enough <- \() shiny::validate(shiny::need(sum(runs()$valid) > 1, "Not enough valid runs: see the reasons below."))
    output$plane <- shiny::renderPlot({ enough(); plot_ce_plane(runs(), wtp()) })
    output$ceac <- shiny::renderPlot({ enough(); plot_ceac(ceac(runs(), 3 * max(wtp(), 1)), wtp()) })

    bands <- shiny::reactive(psa_bands(runs(), model()$spec$arms))
    output$bands <- shiny::renderPlot({
      shiny::validate(shiny::need(sum(runs()$valid) > 1, "Not enough valid runs."))
      plot_psa_bands(bands(), model()$curves, model()$spec$arms, input$band_x %||% 15)
    })

    output$compare <- shiny::renderUI({
      enough()
      cmp <- psa_vs_deterministic(runs(), model()$res)
      lab <- c(cost_comp = "Cost, comparator", cost_int = "Cost, intervention",
               qaly_comp = "QALYs, comparator", qaly_int = "QALYs, intervention",
               ly_comp = "Life-years, comparator", ly_int = "Life-years, intervention",
               d_cost = "Incremental cost", d_qaly = "Incremental QALYs")
      is_cost <- grepl("cost", cmp$measure)
      f <- \(x) ifelse(is_cost, fmt_num(x, 0), fmt_num(x, 3))
      html_table(data.frame(Measure = lab[cmp$measure], Deterministic = f(cmp$deterministic),
                            `PSA mean` = f(cmp$probabilistic), Difference = fmt_pct(cmp$rel_diff),
                            check.names = FALSE))
    })

    output$conv <- shiny::renderPlot({
      shiny::validate(shiny::need(sum(runs()$valid) > 1, "Not enough valid runs."))
      plot_convergence(psa_convergence(runs(), wtp()), input$conv_out)
    })

    output$invalid <- shiny::renderUI({
      s <- summ()
      if (is.null(s$reasons)) return(NULL)
      bslib::card(
        bslib::card_header("Invalid runs (excluded, not fixed)"),
        html_table(data.frame(Reason = names(s$reasons), Runs = as.integer(s$reasons))),
        shiny::helpText("Find the cause before relying on the PSA. Run numbers and reasons are in the Excel export (Doc 13 Case 16).")
      )
    })

    list(runs = shiny::reactive(if (input$run > 0) runs() else NULL),
         summary = shiny::reactive(if (input$run > 0) summ() else NULL),
         bands = shiny::reactive(if (input$run > 0) bands() else NULL),
         compare = shiny::reactive(if (input$run > 0) psa_vs_deterministic(runs(), model()$res) else NULL))
  })
}

# ---- validation tab ---------------------------------------------------------------

mod_validation_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    intro_card(
      purpose = "Automated checks that show what you can trust and what still needs evidence.",
      do = "Read the counts at the top, then go through the table. Look first at anything flagged as a problem.",
      check = "For each warning, ask whether it could change the decision. If it could, fix it or test it as a scenario."
    ),
    shiny::uiOutput(ns("counts")),
    bslib::card(bslib::card_header("Automated checks"), shiny::uiOutput(ns("table"))),
    bslib::card(
      bslib::card_header("How to read this"),
      shiny::p("Computational checks test that the code does what it should. Structural checks test that",
               "states, costs and pathways fit together. Input checks flag assumptions that need evidence.",
               "Nothing is fixed silently: when a constraint binds, it is reported here with where it binds.",
               "Judge each warning by whether it could change the decision (Doc 13, materiality).")
    )
  )
}

mod_validation_server <- function(id, model, psa_summary, psa_compare = \() NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    checks <- shiny::reactive({
      v <- model()$validation
      s <- psa_summary()
      if (!is.null(s)) {
        v <- rbind(v, check_row("Computational", "PSA runs valid",
          if (s$n_invalid == 0) "pass" else "fail",
          sprintf("%d of %d runs valid. Probabilistic ICER = ratio of means.%s", s$n_valid, s$n,
                  if (s$n_invalid) paste0(" Reasons: ", paste(names(s$reasons), s$reasons, sep = " x", collapse = "; "), ".") else ""),
          "Doc 13 Cases 16, 17"))
        st <- s$settings
        v <- rbind(v, check_row("Input", "Survival uncertainty in the PSA",
          if (st$surv_method == "none" && length(st$sampled_curves)) "warn" else "info",
          switch(st$surv_method,
            none = "Fitted survival curves are fixed across runs: the PSA understates uncertainty.",
            mvn = paste("Parameters sampled from their covariance for:", paste(st$sampled_curves, collapse = "; "),
                        ". OS, PFS and ToT models sampled independently."),
            boot = paste0("Patients resampled and models refitted for: ", paste(st$sampled_curves, collapse = "; "),
                          if (isTRUE(st$paired)) ". Paired by patient id (endpoints stay linked)." else ". No patient ids: endpoints resampled separately.")),
          "Doc 13 Case 16"))
        if (!is.null(st$wane_start_range)) {
          v <- rbind(v, check_row("Input", "Waning timing in the PSA", "info",
            sprintf("Waning start sampled uniformly from %s to %s years; structure fixed (test structures in the waning explorer).",
                    st$wane_start_range[1], st$wane_start_range[2]), "Doc 12"))
        }
        cmp <- psa_compare()
        if (!is.null(cmp)) {
          r <- cmp$rel_diff[cmp$measure == "d_qaly"]
          v <- rbind(v, check_row("Computational", "PSA mean vs deterministic (incremental QALYs)",
            if (is.finite(r) && abs(r) > 0.1) "warn" else "info",
            paste0("PSA mean differs from the deterministic value by ", fmt_pct(r),
                   ". Some difference is expected (non-linearity); a large one needs explaining."),
            "Doc 13"))
        }
      }
      v
    })
    output$counts <- shiny::renderUI({
      v <- checks()
      bslib::layout_columns(
        col_widths = c(3, 3, 3, 3),
        bslib::value_box("Passed", sum(v$status == "pass"), theme = "success"),
        bslib::value_box("Warnings", sum(v$status == "warn"), theme = "warning"),
        bslib::value_box("Failed", sum(v$status == "fail"), theme = if (any(v$status == "fail")) "danger" else "secondary"),
        bslib::value_box("Information", sum(v$status == "info"), theme = "info")
      )
    })
    output$table <- shiny::renderUI({
      v <- checks()
      ord <- order(match(v$status, c("fail", "warn", "info", "pass")))
      v <- v[ord, ]
      html_table(data.frame(Status = v$status, Domain = v$domain, Check = v$check,
                            Detail = v$detail, Source = v$source), status_col = "Status")
    })
  })
}
