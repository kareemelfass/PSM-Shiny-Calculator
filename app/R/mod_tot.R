# mod_tot.R : time on treatment (Doc 11) ----------------------------------------

mod_tot_ui <- function(id) {
  ns <- shiny::NS(id)
  arm_controls <- function(a, title, cap_default, beyond_default) {
    bslib::card(
      bslib::card_header(shiny::textOutput(ns(paste0("title_", a)), inline = TRUE)),
      shiny::radioButtons(ns(paste0("method_", a)), help_label("Time-on-treatment curve", "Equal to PFS: patients stay on treatment for as long as they are progression-free. Simple, and useful when there are no time-on-treatment data. Fitted: a curve fitted to your time-on-treatment (ToT) data. Needs ToT rows in the data file."),
                          c("Equal to PFS" = "pfs", "Fitted to ToT data" = "fit")),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'fit'", ns(paste0("method_", a))),
        shiny::selectInput(ns(paste0("model_", a)), "ToT model", NULL)
      ),
      shiny::numericInput(ns(paste0("cap_", a)), help_label("Maximum treatment duration (years); blank = none", "A stopping rule, such as a 2-year course. Leave blank if treatment continues until it stops working."),
                          cap_default, min = 0, step = 0.25),
      shiny::checkboxInput(ns(paste0("beyond_", a)), "Treatment beyond progression allowed", beyond_default),
      shiny::helpText("If ticked, patients can stay on treatment after their disease progresses.")
    )
  }
  shiny::tagList(
    intro_card(
      purpose = "Set how long patients stay on treatment. Drug costs follow time on treatment, not progression-free survival.",
      do = "For each arm choose whether treatment time equals PFS or is fitted to your time-on-treatment data. Add a maximum duration if there is a stopping rule, and say whether treatment can continue after progression.",
      check = "The plot looks sensible: time on treatment never goes above overall survival, and a stopping rule shows as a cut-off.",
      next_step = "the Costs and utilities tab."
    ),
    bslib::layout_columns(
      col_widths = c(6, 6),
      arm_controls("comp", "Comparator", NA, FALSE),
      arm_controls("int", "Intervention", 2, TRUE)
    ),
    bslib::layout_columns(
      col_widths = c(6, 6),
      bslib::card(bslib::card_header("ToT, PFS and OS by arm"), shiny::plotOutput(ns("plot"), height = 380)),
      bslib::card(bslib::card_header("ToT model fit"), shiny::uiOutput(ns("gof")),
                  shiny::helpText("Drug costs follow time on treatment, not PFS. ToT is always capped",
                                  "at OS; it can exceed PFS only if treatment beyond progression is allowed."))
    )
  )
}

mod_tot_server <- function(id, data, curves) {
  shiny::moduleServer(id, function(input, output, session) {
    engine <- fit_engine()
    cat <- available_models(engine)

    output$title_comp <- shiny::renderText(paste("Comparator:", data$arms()$comp))
    output$title_int <- shiny::renderText(paste("Intervention:", data$arms()$int))

    fits <- shiny::reactive({
      d <- data$ipd()
      a <- data$arms()
      if (!any(d$endpoint == "ToT" & d$arm %in% unlist(a))) return(NULL)
      shiny::withProgress(message = "Fitting ToT models",
                          fit_endpoint(d, "ToT", a, engine))
    })

    shiny::observe({
      fe <- fits()
      ch <- setNames(cat$id, cat$label)
      for (a in c("comp", "int")) {
        has <- !is.null(fe) && isTRUE(fe$have[[a]])
        shiny::updateSelectInput(session, paste0("model_", a), choices = ch,
                                 selected = keep_or(shiny::isolate(input[[paste0("model_", a)]]), cat$id, "weibull"))
        if (!has) shiny::updateRadioButtons(session, paste0("method_", a), selected = "pfs")
      }
    })

    spec <- shiny::reactive({
      fe <- fits()
      one <- function(a) {
        method <- input[[paste0("method_", a)]] %||% "pfs"
        if (method == "fit") {
          shiny::validate(shiny::need(!is.null(fe) && isTRUE(fe$have[[a]]),
                                      "No ToT data for this arm: use 'Equal to PFS'."))
        }
        list(method = method, model = input[[paste0("model_", a)]] %||% "weibull",
             cap = input[[paste0("cap_", a)]] %||% NA,
             beyond_prog = isTRUE(input[[paste0("beyond_", a)]]))
      }
      list(comp = one("comp"), int = one("int"))
    })

    output$plot <- shiny::renderPlot({
      cv <- curves()
      a <- data$arms()
      t <- cv$t
      d <- do.call(rbind, lapply(c("comp", "int"), \(k) rbind(
        data.frame(time = t, s = cv[[k]]$OS, curve = "OS", arm = a[[k]]),
        data.frame(time = t, s = cv[[k]]$PFS, curve = "PFS", arm = a[[k]]),
        data.frame(time = t, s = cv[[k]]$ToT, curve = "ToT", arm = a[[k]]))))
      p <- ggplot2::ggplot(d, ggplot2::aes(time, s, colour = curve)) +
        ggplot2::geom_line(linewidth = 0.8) + ggplot2::facet_wrap(~arm)
      fe <- fits()
      if (!is.null(fe)) {
        kms <- do.call(rbind, lapply(names(fe$km), \(k) km_df(fe$km[[k]], a[[k]])))
        p <- p + ggplot2::geom_step(data = kms, ggplot2::aes(time, surv), inherit.aes = FALSE,
                                    colour = "grey40", alpha = 0.6)
      }
      p + ggplot2::coord_cartesian(xlim = c(0, 8)) + pct_axis() + theme_psm() +
        ggplot2::scale_colour_manual(values = c(OS = "#3E5C8A", PFS = "#5DA9A4", ToT = "#D9822B")) +
        ggplot2::labs(x = "Years", y = "Proportion", colour = NULL,
                      subtitle = "Grey steps: ToT Kaplan-Meier (if data supplied).")
    })

    output$gof <- shiny::renderUI({
      fe <- fits()
      if (is.null(fe)) return(shiny::p(class = "text-muted", "No ToT data supplied."))
      g <- fe$gof[fe$gof$approach == "Independent", ]
      html_table(data.frame(Arm = g$arm, Model = g$label, AIC = round(g$aic, 1),
                            dAIC = round(g$d_aic, 1), BIC = round(g$bic, 1)), digits = 1)
    })

    list(fits = fits, spec = spec)
  })
}
