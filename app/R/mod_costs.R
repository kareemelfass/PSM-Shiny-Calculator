# mod_costs.R : costs and utilities ------------------------------------------

econ_fields <- function() {
  data.frame(
    id = c("upfront", "drug_cost", "admin_cost", "interval_days", "rdi", "c_pf", "c_pd",
           "c_subseq", "prog_share", "c_terminal", "ae_cost", "ae_qaly"),
    label = c("Upfront one-off cost (e.g. cell therapy)", "Drug cost per administration",
              "Administration cost per administration", "Dosing interval (days)",
              "Relative dose intensity (0-1)", "Progression-free care, per year",
              "Progressed-disease care, per year", "Subsequent therapy, per progression",
              "Share of PFS events that are progression (0-1)", "Terminal care, per death",
              "Adverse events, one-off cost", "Adverse events, one-off QALY loss")
  )
}

mod_costs_ui <- function(id) {
  ns <- shiny::NS(id)
  f <- econ_fields()
  arm_card <- function(a) {
    d <- default_econ(a)
    bslib::card(
      bslib::card_header(shiny::textOutput(ns(paste0("title_", a)), inline = TRUE)),
      lapply(seq_len(nrow(f)), \(i) shiny::numericInput(
        ns(paste0(a, "_", f$id[i])), f$label[i], d[[f$id[i]]], min = 0))
    )
  }
  shiny::tagList(
    bslib::layout_columns(
      col_widths = c(4, 4, 4),
      bslib::card(
        bslib::card_header("Utilities"),
        shiny::numericInput(ns("u_pf"), "Progression-free", 0.75, min = 0, max = 1, step = 0.01),
        shiny::numericInput(ns("u_pf_se"), "SE (for PSA and DSA)", 0.03, min = 0, step = 0.005),
        shiny::numericInput(ns("u_pd"), "Progressed disease", 0.60, min = 0, max = 1, step = 0.01),
        shiny::numericInput(ns("u_pd_se"), "SE (for PSA and DSA)", 0.04, min = 0, step = 0.005),
        shiny::helpText("Drug costs use time on treatment x (365.25 / interval) x relative dose intensity.",
                        "Subsequent-therapy cost applies at each new progression and should match",
                        "the crossover setting on the Treatment effects tab (Doc 10).")
      ),
      arm_card("comp"),
      arm_card("int")
    )
  )
}

mod_costs_server <- function(id, arms) {
  shiny::moduleServer(id, function(input, output, session) {
    output$title_comp <- shiny::renderText(paste("Comparator:", arms()$comp))
    output$title_int <- shiny::renderText(paste("Intervention:", arms()$int))
    shiny::reactive({
      f <- econ_fields()$id
      get_arm <- \(a) setNames(lapply(f, \(x) input[[paste0(a, "_", x)]] %||% 0), f)
      e <- list(comp = get_arm("comp"), int = get_arm("int"),
                utility = list(pf = input$u_pf, pd = input$u_pd,
                               pf_se = input$u_pf_se %||% 0, pd_se = input$u_pd_se %||% 0))
      shiny::validate(
        shiny::need(e$comp$interval_days > 0 && e$int$interval_days > 0, "Dosing interval must be above 0."),
        shiny::need(is.finite(e$utility$pf) && is.finite(e$utility$pd), "Enter both utilities.")
      )
      e
    })
  })
}
