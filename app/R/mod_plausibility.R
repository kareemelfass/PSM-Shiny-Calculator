# mod_plausibility.R : landmarks, benchmarks, hazards (Doc 7) --------------------

mod_plausibility_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 330,
      shiny::textInput(ns("years"), help_label("Landmark years", "Time points, in years, at which survival is reported. Separate them with commas."), "1, 2, 5, 10, 20"),
      shiny::h6("Clinical benchmark ranges for OS (%), one per landmark"),
      shiny::textInput(ns("bench_comp"), "Comparator", "", placeholder = "e.g. 60-75, 35-50, 8-15, 1-4, 0-1"),
      shiny::textInput(ns("bench_int"), "Intervention", "", placeholder = "same format"),
      shiny::helpText("Ranges from clinical experts, registries or longer follow-up of related trials.",
                      "Record their source. Leave blank to skip the check."),
      shiny::numericInput(ns("xmax"), "Plot to year", 30, min = 1)
    ),
    intro_card(
      purpose = "Check that the long-term predictions are believable, using clinical knowledge and not only fit statistics.",
      do = "Enter the landmark years you want to see. If you have them, enter expected survival ranges from clinical experts, registries or longer follow-up of related trials. Then review the curves, the hazard chart and the landmark table.",
      check = "Survival at each landmark falls inside the expected range, and the modelled risk of death does not drop below that of the general population.",
      next_step = "the Treatment effects tab."
    ),
    bslib::layout_columns(
      col_widths = c(6, 6),
      bslib::card(bslib::card_header("Model survival (after effects, waning and floor)"),
                  shiny::plotOutput(ns("curves"), height = 380)),
      bslib::card(bslib::card_header("Modelled OS hazard vs general population"),
                  shiny::plotOutput(ns("hazard"), height = 380))
    ),
    bslib::card(bslib::card_header("Landmark survival of the selected curves"),
                shiny::uiOutput(ns("landmarks"))),
    bslib::card(
      bslib::card_header("Candidate models for comparator OS: fit vs long-term prediction"),
      shiny::uiOutput(ns("candidates")),
      shiny::helpText("Raw fitted curves (no effects, waning or floor). Mean = area under the curve over",
                      "the horizon, undiscounted. Similar AIC can hide very different futures (Doc 7).")
    )
  )
}

mod_plausibility_server <- function(id, curves, fits_os, arms, settings) {
  shiny::moduleServer(id, function(input, output, session) {
    years <- shiny::reactive({
      y <- suppressWarnings(as.numeric(trimws(strsplit(input$years, ",")[[1]])))
      y <- y[is.finite(y) & y > 0]
      shiny::validate(shiny::need(length(y) > 0, "Enter landmark years, e.g. 1, 2, 5, 10"))
      y
    })

    output$curves <- shiny::renderPlot({
      plot_model_curves(curves(), arms(), NULL, c("OS", "PFS"), TRUE, input$xmax)
    })
    output$hazard <- shiny::renderPlot({
      plot_model_hazards(curves(), arms(), "OS") +
        ggplot2::coord_cartesian(xlim = c(0, input$xmax))
    })

    output$landmarks <- shiny::renderUI({
      y <- years()
      tb <- landmark_table(curves(), arms(), y)
      flag <- rep("", nrow(tb))
      for (a in c("comp", "int")) {
        b <- parse_bench(input[[paste0("bench_", a)]], length(y))
        row <- which(tb$arm == arms()[[a]] & tb$endpoint == "OS")
        if (!is.null(b) && length(row)) {
          v <- unlist(tb[row, paste0("yr_", y)])
          out <- which(v < b[, 1] | v > b[, 2])
          flag[row] <- if (length(out)) paste("Outside benchmark at year", paste(y[out], collapse = ", ")) else
            "Within benchmark"
        }
      }
      show <- tb
      for (cn in paste0("yr_", y)) show[[cn]] <- fmt_pct(tb[[cn]])
      names(show) <- c("Arm", "Endpoint", paste("Year", y))
      show$Benchmark <- flag
      html_table(show)
    })

    output$candidates <- shiny::renderUI({
      fe <- fits_os()
      y <- years()
      st <- settings()
      tt <- make_grid(st$horizon, 1 / 52)
      rows <- lapply(fe$indep$comp, \(f) {
        if (!isTRUE(f$ok)) return(NULL)
        s <- f$S(tt, fe$arms$comp)
        v <- stats::approx(tt, s, xout = y, rule = 2)$y
        data.frame(Model = f$label, AIC = round(f$aic, 1), BIC = round(f$bic, 1),
                   t(setNames(fmt_pct(v), paste("Year", y))),
                   `Mean (yr)` = round(sum(diff(tt) * (s[-1] + s[-length(s)]) / 2), 2),
                   check.names = FALSE)
      })
      tb <- do.call(rbind, rows)
      b <- parse_bench(input$bench_comp, length(y))
      if (!is.null(b)) {
        tb$Benchmark <- vapply(fe$indep$comp[vapply(fe$indep$comp, \(f) isTRUE(f$ok), logical(1))], \(f) {
          v <- stats::approx(tt, f$S(tt, fe$arms$comp), xout = y, rule = 2)$y
          if (any(v < b[, 1] | v > b[, 2])) "Outside" else "Within"
        }, character(1))
      }
      html_table(tb[order(tb$AIC), ], digits = 1)
    })
  })
}
