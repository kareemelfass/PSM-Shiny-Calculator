# mod_endpoint.R : fit and choose curves for one endpoint (OS or PFS) -----------

mod_endpoint_ui <- function(id, endpoint) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 330,
      shiny::h5(paste(endpoint, "model choice")),
      shiny::radioButtons(ns("approach"), help_label("Treatment-effect approach", "Independent fits: each arm gets its own curve. Most flexible, but needs data on both arms. Joint fit: one model for both arms with arm as a covariate. Fewer parameters, and it assumes the arms differ in a simple, consistent way. Comparator curve + hazard ratio: fit only the comparator, then apply a hazard ratio to get the intervention. Use this when you have no patient data for the intervention."), c(
        "Independent fits per arm" = "indep",
        "Joint fit (arm as covariate)" = "joint",
        "Comparator curve + hazard ratio" = "hr")),
      shiny::selectInput(ns("comp_model"), help_label("Comparator model (joint: both arms)", "The statistical distribution used to extend survival over time, such as exponential, Weibull or log-normal. Different choices can give very different long-term results, so compare them on the Fit statistics and Curves tabs."), NULL),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'indep'", ns("approach")),
        shiny::selectInput(ns("int_model"), "Intervention model", NULL)
      ),
      shiny::checkboxInput(ns("hybrid"), "Hybrid: KM to cut-point, then model", FALSE),
      shiny::helpText("Uses the observed curve up to the cut-point, then the fitted model after it."),
      shiny::conditionalPanel(
        sprintf("input['%s']", ns("hybrid")),
        shiny::numericInput(ns("cut"), help_label("Cut-point (years)", "The time after which the fitted model takes over from the observed curve. Choose a time where the observed data are still reliable, with enough patients still at risk."), 1, min = 0.05, step = 0.25),
        shiny::uiOutput(ns("cut_info"))
      ),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'hr'", ns("approach")),
        shiny::hr(),
        shiny::radioButtons(ns("hr_src"), help_label("HR source", "Enter: type a hazard ratio you already have, such as from a published trial or a network meta-analysis. Cox model: estimate it from the data you loaded."),
                            c("Enter (e.g. NMA, published)" = "manual", "Cox model on these data" = "cox")),
        shiny::conditionalPanel(
          sprintf("input['%s'] == 'manual'", ns("hr_src")),
          shiny::numericInput(ns("hr"), help_label("HR point estimate", "Below 1 means the intervention lowers the risk of the event compared with the comparator."), 0.7, min = 0.01, step = 0.01),
          shiny::numericInput(ns("hr_lo"), help_label("95% interval, lower", "The 95% interval shows how uncertain the HR is. It is used in the sensitivity analyses (DSA and PSA)."), 0.55, min = 0.01, step = 0.01),
          shiny::numericInput(ns("hr_hi"), "95% interval, upper", 0.89, min = 0.01, step = 0.01),
          shiny::radioButtons(ns("hr_dir"), help_label("The HR is reported as", "Pick the direction the source uses. If it reports comparator vs intervention, the app inverts it for you."), c(
            "Intervention vs comparator" = "int_vs_comp",
            "Comparator vs intervention (will be inverted)" = "comp_vs_int"))
        ),
        shiny::uiOutput(ns("hr_info"))
      )
    ),
    intro_card(
      purpose = paste0("Choose how the ", endpoint, " curve is built for each arm, and test whether it can be trusted beyond the trial data."),
      do = "In the left panel choose the approach and a model for each arm. Then open Fit statistics to compare models, Curves to see them against the observed data, and Hazards and PH diagnostics to check the assumptions.",
      check = "The chosen curve follows the observed data AND its long-term tail is clinically believable. The model with the best statistical fit is not always the best extrapolation.",
      next_step = if (identical(endpoint, "OS")) "the Progression-free survival tab." else "the Clinical plausibility tab."
    ),
    bslib::navset_card_tab(
      bslib::nav_panel("Fit statistics",
        shiny::uiOutput(ns("engine_note")),
        shiny::uiOutput(ns("gof")),
        shiny::helpText("AIC and BIC rank models within one dataset only. The best statistical fit",
                        "is not necessarily the best extrapolation (Doc 6). BIC uses n = patients.")),
      bslib::nav_panel("Curves",
        bslib::layout_columns(col_widths = c(3, 9),
          shiny::div(
            shiny::radioButtons(ns("plot_arm"), "Arm", c("Comparator" = "comp", "Intervention" = "int")),
            shiny::checkboxGroupInput(ns("plot_models"), "Models", NULL),
            shiny::numericInput(ns("plot_h"), "Plot to year", 20, min = 1)
          ),
          shiny::plotOutput(ns("fit_plot"), height = 430))),
      bslib::nav_panel("Hazards", shiny::plotOutput(ns("haz_plot"), height = 430)),
      bslib::nav_panel("PH diagnostics",
        shiny::uiOutput(ns("ph_text")),
        bslib::layout_columns(col_widths = c(6, 6),
          shiny::plotOutput(ns("lch_plot"), height = 360),
          shiny::plotOutput(ns("zph_plot"), height = 360)))
    )
  )
}

mod_endpoint_server <- function(id, endpoint, data, settings) {
  shiny::moduleServer(id, function(input, output, session) {
    engine <- fit_engine()
    cat <- available_models(engine)

    fits <- shiny::reactive({
      d <- data$ipd()
      a <- data$arms()
      shiny::validate(shiny::need(any(d$endpoint == endpoint & d$arm == a$comp),
                                  paste("No", endpoint, "data for the comparator arm.")))
      shiny::withProgress(message = paste("Fitting", endpoint, "models"), value = 0, {
        n_steps <- nrow(cat) * 3
        i <- 0
        fit_endpoint(d, endpoint, a, engine, progress = \(msg) {
          i <<- i + 1
          shiny::incProgress(1 / n_steps, detail = msg)
        })
      })
    })

    # Keep choices valid for the data available
    shiny::observe({
      fe <- fits()
      ch <- setNames(cat$id, cat$label)
      shiny::updateSelectInput(session, "comp_model", choices = ch,
                               selected = keep_or(shiny::isolate(input$comp_model), cat$id, "weibull"))
      shiny::updateSelectInput(session, "int_model", choices = ch,
                               selected = keep_or(shiny::isolate(input$int_model), cat$id, "weibull"))
      shiny::updateCheckboxGroupInput(session, "plot_models", choices = ch,
                                      selected = intersect(c("exp", "weibull", "lnorm", "llogis", "gompertz"), cat$id))
      if (!fe$have[["int"]]) {
        shiny::updateRadioButtons(session, "approach", selected = "hr")
        shiny::updateRadioButtons(session, "hr_src", selected = "manual")
      }
    })

    cox <- shiny::reactive({
      fe <- fits()
      if (all(fe$have)) cox_summary(fe) else NULL
    })

    hr_spec <- shiny::reactive({
      if (identical(input$hr_src, "cox")) {
        cx <- cox()
        shiny::validate(shiny::need(!is.null(cx), "A Cox HR needs data in both arms."))
        z <- stats::qnorm(0.975)
        list(hr = cx$hr, se = cx$se_log, lower = exp(log(cx$hr) - z * cx$se_log),
             upper = exp(log(cx$hr) + z * cx$se_log), source = "Cox model on trial data")
      } else {
        shiny::req(input$hr > 0)
        h <- hr_to_log(input$hr, input$hr_lo %||% NA, input$hr_hi %||% NA, input$hr_dir)
        c(h, source = "Entered")
      }
    })

    spec <- shiny::reactive({
      fe <- fits()
      shiny::req(input$comp_model, input$approach)
      if (input$approach != "hr") {
        shiny::validate(shiny::need(fe$have[["int"]],
          paste("No intervention", endpoint, "data: use the 'Comparator curve + hazard ratio' approach.")))
      }
      list(approach = input$approach, comp_model = input$comp_model,
           int_model = input$int_model %||% input$comp_model,
           hybrid = isTRUE(input$hybrid), cut = input$cut %||% 1,
           hr = if (input$approach == "hr") hr_spec() else list(hr = 1, se = NA))
    })

    output$engine_note <- shiny::renderUI({
      if (engine == "flexsurv") return(NULL)
      shiny::div(class = "alert alert-warning",
        "The flexsurv package is not installed, so only 4 models (survreg) are available.",
        "Install flexsurv for Gompertz, gamma, generalised gamma and Royston-Parmar splines.")
    })

    output$gof <- shiny::renderUI({
      g <- fits()$gof
      shiny::validate(shiny::need(nrow(g) > 0, "No models fitted."))
      show <- data.frame(Approach = g$approach, Arm = g$arm, Model = g$label, k = g$npar,
                         AIC = round(g$aic, 1), `dAIC` = round(g$d_aic, 1),
                         BIC = round(g$bic, 1), `dBIC` = round(g$d_bic, 1),
                         Note = g$note, check.names = FALSE)
      html_table(show, digits = 1)
    })

    output$cut_info <- shiny::renderUI({
      fe <- fits()
      txt <- vapply(names(fe$km), \(a) {
        km <- fe$km[[a]]
        sprintf("%s: %d of %d at risk", fe$arms[[a]], km_at_risk(km, input$cut), km$n)
      }, character(1))
      shiny::helpText(paste(txt, collapse = "; "),
                      ". No universal rule: justify with risk set, events and hazard stability (Doc 6).")
    })

    output$hr_info <- shiny::renderUI({
      h <- hr_spec()
      cx <- cox()
      shiny::div(class = "small mt-2",
        shiny::p(sprintf("Applied: HR = %.3f (%.3f to %.3f), intervention vs comparator.",
                         h$hr, h$lower %||% NA, h$upper %||% NA)),
        shiny::p(sprintf("S_int(t) = S_comp(t)^%.3f. Survival is never multiplied by the HR (Doc 3).", h$hr)),
        if (!is.null(cx)) shiny::p(sprintf("For reference, trial Cox HR = %.3f; PH test p = %.3f.",
                                           cx$hr, cx$p_ph))
      )
    })

    output$fit_plot <- shiny::renderPlot({
      fe <- fits()
      shiny::validate(shiny::need(fe$have[[input$plot_arm]], "No data for this arm."),
                      shiny::need(length(input$plot_models) > 0, "Select at least one model."))
      plot_km_fits(fe, input$plot_arm, input$plot_models, input$plot_h %||% 20)
    })

    output$haz_plot <- shiny::renderPlot({
      fe <- fits()
      shiny::validate(shiny::need(fe$have[[input$plot_arm]], "No data for this arm."),
                      shiny::need(length(input$plot_models) > 0, "Select models on the Curves tab."))
      plot_fit_hazards(fe, input$plot_arm, input$plot_models, input$plot_h %||% 20)
    })

    output$ph_text <- shiny::renderUI({
      cx <- cox()
      shiny::validate(shiny::need(!is.null(cx), "PH diagnostics need data in both arms."))
      shiny::p(sprintf(
        "Cox HR (intervention vs comparator) = %.3f, 95%% CI %.3f to %.3f. Schoenfeld test p = %.3f. %s",
        cx$hr, exp(log(cx$hr) - 1.96 * cx$se_log), exp(log(cx$hr) + 1.96 * cx$se_log), cx$p_ph,
        "A non-significant test does not establish that PH holds (Doc 3)."))
    })
    output$lch_plot <- shiny::renderPlot({
      fe <- fits()
      shiny::req(length(fe$km) > 0)
      plot_logcumhaz(fe)
    })
    output$zph_plot <- shiny::renderPlot({
      cx <- cox()
      shiny::req(!is.null(cx), !is.null(cx$zph))
      plot_schoenfeld(cx$zph)
    })

    list(fits = fits, spec = spec)
  })
}
