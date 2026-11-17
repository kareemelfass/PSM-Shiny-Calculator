# mod_effects.R : treatment-effect duration (waning) and crossover (Docs 10, 12) ----
# V3: structures labelled with their appraisal role, timing anchored to treatment
# start or treatment stop, target HR, individual-level (frailty) waning,
# justification, structure comparison, absolute gap, and a timing explorer.

mod_effects_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 360,
      shiny::h5("Treatment-effect duration"),
      shiny::radioButtons(ns("type"), NULL, c(
        "Sustained effect (no waning)" = "none",
        "Equal hazards from a time point" = "equal",
        "Gradual equalisation (HR moves to target)" = "gradual",
        "Stop and drop (stress test only)" = "stop_drop"), selected = "none"),
      shiny::uiOutput(ns("role")),
      shiny::conditionalPanel(
        sprintf("input['%s'] != 'none'", ns("type")),
        shiny::radioButtons(ns("anchor"), "Count waning time from", c(
          "Treatment start (model start)" = "start",
          "Treatment stop (maximum duration, Time on treatment tab)" = "tx_stop")),
        shiny::numericInput(ns("start"), "Waning starts after (years)", 5, min = 0, step = 0.5),
        shiny::conditionalPanel(
          sprintf("input['%s'] == 'gradual'", ns("type")),
          shiny::numericInput(ns("end"), "Fully waned after (years)", 7, min = 0, step = 0.5)
        ),
        shiny::uiOutput(ns("abs_time")),
        shiny::checkboxInput(ns("apply_os"), "Apply to OS", TRUE),
        shiny::checkboxInput(ns("apply_pfs"), "Apply to PFS", TRUE),
        bslib::accordion(
          open = FALSE,
          bslib::accordion_panel(
            "Advanced",
            shiny::numericInput(ns("target_hr"), "HR after waning (1 = no relative effect left)",
                                1, min = 0.05, step = 0.05),
            shiny::numericInput(ns("frailty"), "Frailty variance (0 = wane population hazards)",
                                0, min = 0, step = 0.1),
            shiny::helpText(
              "Frailty > 0 wanes the individual-level effect. Surviving intervention patients are",
              "frailer, so their population hazard ends up above the comparator's: the marginal HR",
              "goes above 1. Forcing the marginal HR to 1 overestimates benefit (Jennings 2024).",
              "Gamma frailty; use as a sensitivity analysis (e.g. 0.25 to 1).")
          )
        )
      ),
      shiny::textAreaInput(ns("rationale"), "Evidence and justification", "", rows = 4,
        placeholder = "e.g. Treatment stops at 2 years; no long-term data beyond year 4; mechanism does not support lifetime effect; scenarios tested at 3 and 5 years after stopping."),
      shiny::hr(),
      shiny::h5("Treatment crossover"),
      shiny::radioButtons(ns("basis"), "Comparator OS represents", c(
        "Observed ITT (patients could switch to the intervention)" = "itt",
        "Adjusted to a no-switching pathway (RPSFTM, IPCW, two-stage)" = "adjusted")),
      shiny::checkboxInput(ns("subseq_int"),
                           "Comparator subsequent-therapy cost includes the intervention", TRUE),
      shiny::helpText("Survival and costs must describe the same treatment pathway (Doc 10).")
    ),
    bslib::card(bslib::card_header("Timeline"), shiny::plotOutput(ns("timeline"), height = 170)),
    bslib::layout_columns(
      col_widths = c(6, 6),
      bslib::card(bslib::card_header("Intervention OS under each assumption"),
                  shiny::plotOutput(ns("structures"), height = 360)),
      bslib::card(bslib::card_header("Absolute OS advantage over the comparator"),
                  shiny::plotOutput(ns("gap"), height = 360))
    ),
    bslib::layout_columns(
      col_widths = c(6, 6),
      bslib::card(bslib::card_header("Implied hazard ratio (current choice)"),
                  shiny::plotOutput(ns("hr"), height = 330)),
      bslib::card(bslib::card_header("OS hazards (current choice)"),
                  shiny::plotOutput(ns("haz"), height = 330))
    ),
    bslib::card(
      bslib::card_header("Waning explorer: timings x structures"),
      bslib::layout_columns(
        col_widths = c(4, 3, 3, 2),
        shiny::textInput(ns("ex_starts"), "Start times to test (same anchor as above)", "2, 3, 5, 7, 10"),
        shiny::numericInput(ns("ex_width"), "Gradual: years to full waning", 2, min = 0.25, step = 0.5),
        shiny::radioButtons(ns("ex_out"), "Show", c("ICER" = "icer", "Inc. QALYs" = "inc_qaly",
                                                   "NMB" = "nmb"), inline = TRUE),
        shiny::actionButton(ns("ex_run"), "Run explorer", class = "btn-primary mt-4")
      ),
      shiny::helpText("Reviewers test timings and convergence structures rather than accept indefinite",
                      "efficacy (NICE TA798, TA1099). Uses all current settings, including target HR and frailty."),
      shiny::plotOutput(ns("ex_plot"), height = 360),
      shiny::uiOutput(ns("ex_table"))
    ),
    shiny::uiOutput(ns("note"))
  )
}

# What the explorer depends on: the whole specification except free-text notes
explorer_key <- function(spec) {
  spec$waning$rationale <- NULL
  spec
}

mod_effects_server <- function(id, model, arms, fits) {
  shiny::moduleServer(id, function(input, output, session) {
    out <- shiny::reactive({
      list(
        waning = list(type = input$type, anchor = input$anchor %||% "start",
                      start = input$start %||% Inf,
                      end = if (input$type == "gradual") input$end %||% Inf else input$start %||% Inf,
                      target_hr = if (is.finite(input$target_hr %||% NA) && input$target_hr > 0) input$target_hr else 1,
                      frailty = max(input$frailty %||% 0, 0),
                      apply_os = isTRUE(input$apply_os), apply_pfs = isTRUE(input$apply_pfs),
                      rationale = input$rationale %||% ""),
        crossover = list(comp_os_basis = input$basis, subseq_includes_int = isTRUE(input$subseq_int))
      )
    })

    output$role <- shiny::renderUI({
      ws <- waning_structures()
      role <- ws$role[ws$type == input$type]
      cls <- switch(input$type, none = "text-primary", stop_drop = "text-danger", "text-success")
      shiny::p(class = paste("small", cls), role)
    })

    output$abs_time <- shiny::renderUI({
      m <- model()
      wt <- m$curves$waning
      msg <- if (input$type == "gradual") sprintf("Waning runs from year %s to year %s of the model.", round(wt$start, 2), round(wt$end, 2)) else
        sprintf("Waning at year %s of the model.", round(wt$start, 2))
      if (isTRUE(wt$anchor_missing_cap)) msg <- paste(msg, "No maximum treatment duration is set, so this counts from treatment start.")
      shiny::p(class = "small text-muted", msg)
    })

    xmax <- shiny::reactive({
      wt <- model()$curves$waning
      max(10, min(30, (if (is.finite(wt$end)) wt$end else 5) + 8))
    })

    output$timeline <- shiny::renderPlot({
      m <- model()
      wt <- m$curves$waning
      plot_waning_timeline(m$spec$tot$int$cap %||% NA, m$spec$settings$fu_years %||% NA,
                           wt$start, wt$end, wt$type, xmax())
    })

    structures <- shiny::reactive(compare_structures(model()$pieces, model()$spec, "OS"))
    current_label <- shiny::reactive(waning_structures()$label[waning_structures()$type == input$type])
    output$structures <- shiny::renderPlot(plot_structures_compare(structures(), current_label(), xmax()))
    output$gap <- shiny::renderPlot({
      wt <- model()$curves$waning
      lm <- unique(round(c(if (is.finite(wt$start)) wt$start, 5, 10), 1))
      plot_gap_compare(structures(), xmax(), lm[lm <= xmax()])
    })
    output$hr <- shiny::renderPlot(plot_implied_hr(model()$curves, xmax()))
    output$haz <- shiny::renderPlot(plot_model_hazards(model()$curves, arms(), "OS") +
                                      ggplot2::coord_cartesian(xlim = c(0, xmax())))

    explorer <- shiny::eventReactive(input$ex_run, {
      st <- suppressWarnings(as.numeric(trimws(strsplit(input$ex_starts, ",")[[1]])))
      st <- sort(unique(st[is.finite(st) & st >= 0]))
      shiny::validate(shiny::need(length(st) > 0, "Enter start times, e.g. 2, 3, 5"))
      m <- model()
      g <- shiny::withProgress(message = "Running waning explorer",
                               waning_grid(m$pieces, m$spec, st, input$ex_width %||% 2))
      attr(g, "spec") <- explorer_key(m$spec)
      g
    })
    # Results go stale when any model input changes after the run
    stale <- shiny::reactive(!identical(attr(explorer(), "spec"), explorer_key(model()$spec)))
    output$ex_plot <- shiny::renderPlot({
      plot_waning_explorer(explorer(), model()$spec$settings$wtp, input$ex_out, model()$spec$settings$currency)
    })
    output$ex_table <- shiny::renderUI({
      g <- explorer()
      lm <- grep("^int_os_yr", names(g), value = TRUE)
      show <- data.frame(
        Structure = g$structure,
        Starts = ifelse(is.na(g$start_year), "", fmt_num(g$start_year, 1)),
        `Fully waned` = ifelse(is.na(g$end_year), "", fmt_num(g$end_year, 1)),
        `Inc. QALYs` = fmt_num(g$inc_qaly, 3), `Inc. cost` = fmt_num(g$inc_cost, 0),
        ICER = ifelse(is.finite(g$icer), fmt_num(g$icer, 0), "n/a"),
        check.names = FALSE)
      for (k in lm) show[[paste("Int. OS yr", sub("int_os_yr", "", k))]] <- fmt_pct(g[[k]])
      show$status <- ifelse(g$cost_effective == g$cost_effective[g$type == "none"], "pass", "warn")
      show$status[g$type == "none"] <- "info"
      names(show)[names(show) == "status"] <- "Conclusion vs sustained"
      shiny::tagList(
        if (stale()) shiny::div(class = "alert alert-warning",
          "Model inputs changed after this run. Click 'Run explorer' again before using these results."),
        html_table(show, status_col = "Conclusion vs sustained"),
                     shiny::helpText("warn = this timing/structure changes the cost-effectiveness conclusion",
                                     "at the current threshold. Included in the Excel export."))
    })

    output$note <- shiny::renderUI({
      if (input$type == "gradual" && is.finite(input$end %||% NA) && input$end <= input$start) {
        shiny::div(class = "alert alert-warning", "End is not after start: treated as equal hazards from the start time.")
      }
    })

    list(spec = out, explorer = shiny::reactive(if (input$ex_run > 0) explorer() else NULL))
  })
}
