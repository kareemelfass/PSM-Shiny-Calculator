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

# Plain-language explanation of each cost input, shown in an info popup
econ_help <- function() {
  list(
    upfront = "A one-off cost paid once at the start, for example a cell therapy or a procedure.",
    drug_cost = "The price of the drug each time it is given.",
    admin_cost = "The cost of giving each dose, for example infusion time.",
    interval_days = "Days between doses. Used with time on treatment to count doses per year.",
    rdi = "The share of the planned dose that patients actually receive. 1 means every dose in full.",
    c_pf = "Yearly cost of routine care while the disease has not progressed.",
    c_pd = "Yearly cost of routine care after the disease has progressed.",
    c_subseq = "Cost of the next treatments given when the disease progresses.",
    prog_share = "Of all progression-free survival events, the share that are progression. The rest are deaths.",
    c_terminal = "One-off cost of end-of-life care, applied at each death.",
    ae_cost = "One-off cost of treating adverse events.",
    ae_qaly = "One-off loss of quality of life from adverse events, in QALYs."
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
        ns(paste0(a, "_", f$id[i])), help_label(f$label[i], econ_help()[[f$id[i]]]),
        d[[f$id[i]]], min = 0))
    )
  }
  shiny::tagList(
    intro_card(
      purpose = "Enter what the treatment and care cost, and the quality of life in each health state. With survival, these give the costs and QALYs.",
      do = "Enter a utility for the progression-free and the progressed state, between 0 and 1. Then, for each arm, enter the drug, administration, care, subsequent-therapy, terminal-care and adverse-event inputs. Use the upfront cost for one-off treatments.",
      check = "Every default has been replaced by your own evidence, and the subsequent-therapy cost matches the crossover choice on the Treatment effects tab.",
      next_step = "the Results menu, starting with Base case."
    ),
    bslib::layout_columns(
      col_widths = c(4, 4, 4),
      bslib::card(
        bslib::card_header("Utilities"),
        shiny::numericInput(ns("u_pf"), help_label("Progression-free", "Quality of life while the disease has not progressed. 1 is full health and 0 is dead."), 0.75, min = 0, max = 1, step = 0.01),
        shiny::numericInput(ns("u_pf_se"), help_label("SE (for PSA and DSA)", "Standard error: how uncertain the utility is. A larger value gives a wider spread in the sensitivity analyses."), 0.03, min = 0, step = 0.005),
        shiny::numericInput(ns("u_pd"), help_label("Progressed disease", "Quality of life after the disease has progressed. 1 is full health and 0 is dead."), 0.60, min = 0, max = 1, step = 0.01),
        shiny::numericInput(ns("u_pd_se"), help_label("SE (for PSA and DSA)", "Standard error: how uncertain the utility is. A larger value gives a wider spread in the sensitivity analyses."), 0.04, min = 0, step = 0.005),
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
