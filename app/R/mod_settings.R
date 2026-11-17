# mod_settings.R : model settings ----------------------------------------------

mod_settings_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
  intro_card(
    purpose = "Set the basic rules of the model: how long it runs, how future costs and benefits are discounted, and who the patients are.",
    do = "Check each box. The defaults are a sensible starting point, but they are illustrative. Change the time horizon, discount rates, willingness to pay, and patient age and sex to fit your appraisal.",
    check = "The time horizon is long enough to cover the patients' lifetime, and the discount rates and willingness to pay match what your decision-maker requires.",
    next_step = "the Data tab."
  ),
  bslib::layout_columns(
    col_widths = c(4, 4, 4),
    bslib::card(
      bslib::card_header("Time"),
      shiny::numericInput(ns("horizon"), help_label("Time horizon (years)", "How far into the future the model runs. For a life-limiting disease it should be long enough for nearly all patients to have died. A horizon that is too short cuts off lifetime costs and benefits."), 30, min = 1, max = 100),
      shiny::selectInput(ns("cycle"), help_label("Cycle length", "The size of each model step. Weekly is more precise but has more steps. Monthly is coarser."),
                         c("Weekly" = "week", "Monthly" = "month"), "week"),
      shiny::checkboxInput(ns("hcc"), "Half-cycle correction (trapezoid)", TRUE),
      shiny::numericInput(ns("fu"), help_label("Trial follow-up (years); blank = longest OS time in data", "How long patients were followed in the trial. The app uses it to tell observed time from extrapolated time, for example how much of the result comes from beyond the data. Leave blank to use the longest OS time in your data."),
                          NA, min = 0)
    ),
    bslib::card(
      bslib::card_header("Economics"),
      shiny::numericInput(ns("dr_cost"), help_label("Discount rate, costs (%)", "Yearly rate at which future costs count for less than costs today. Use the rate your decision-maker requires (3.5% is the NICE reference case)."), 3.5, min = 0, step = 0.5),
      shiny::numericInput(ns("dr_out"), help_label("Discount rate, outcomes (%)", "Yearly rate at which future health gains (QALYs) count for less than gains today."), 3.5, min = 0, step = 0.5),
      shiny::textInput(ns("currency"), help_label("Currency label", "Only changes the text shown on screen. It does not convert any amounts."), "CU"),
      shiny::numericInput(ns("wtp"), help_label("Willingness to pay per QALY", "The most the decision-maker will pay for one extra QALY. The ICER is compared with this value to judge whether the intervention is good value."), 30000, min = 0, step = 1000),
      shiny::helpText("All default inputs are illustrative. Replace them with your own evidence.")
    ),
    bslib::card(
      bslib::card_header("Population and mortality"),
      shiny::numericInput(ns("age"), help_label("Starting age (years)", "Average age when patients start. With the life table it sets the lowest death rate the model allows, that of the general population of the same age."), 60, min = 18, max = 99),
      shiny::sliderInput(ns("male"), help_label("Male (%)", "Share of male patients. Used to blend the male and female columns of the life table."), 0, 100, 50),
      shiny::checkboxInput(ns("floor_on"), "Apply general-population mortality floor", TRUE),
      shiny::helpText("Stops modelled death rates from falling below those of the general population of the same age and sex."),
      shiny::radioButtons(ns("lt_src"), help_label("Life table", "Background death rates of the general population by age and sex. The built-in table is illustrative only. For real use, upload a national life table with columns age, qx_male and qx_female, where qx is the chance of dying within one year."),
                          c("Illustrative (replace for real use)" = "example", "Upload" = "upload")),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'upload'", ns("lt_src")),
        shiny::fileInput(ns("lt_file"), "CSV: age, qx_male, qx_female", accept = ".csv")
      ),
      shiny::checkboxInput(ns("pfs_le_os"), "Constrain PFS to OS (reported in Validation)", TRUE),
      shiny::helpText("PFS cannot logically be higher than OS. If the curves cross, this forces PFS down to OS and reports where it happened.")
    )
  )
  )
}

mod_settings_server <- function(id, example_dir) {
  shiny::moduleServer(id, function(input, output, session) {
    lifetable <- shiny::reactive({
      if (input$lt_src == "upload") {
        shiny::req(input$lt_file)
        read_lifetable(input$lt_file$datapath)
      } else {
        read_lifetable(file.path(example_dir, "illustrative_lifetable.csv"))
      }
    })
    shiny::reactive({
      shiny::req(input$horizon > 0, input$age)
      list(
        horizon = input$horizon,
        cycle = if (input$cycle == "week") 7 / 365.25 else 1 / 12,
        hcc = input$hcc,
        dr_cost = (input$dr_cost %||% 0) / 100,
        dr_out = (input$dr_out %||% 0) / 100,
        currency = input$currency,
        wtp = input$wtp %||% 0,
        age = input$age,
        male = input$male / 100,
        floor_on = input$floor_on,
        lifetable = if (isTRUE(input$floor_on)) lifetable() else NULL,
        pfs_le_os = input$pfs_le_os,
        fu_override = if (is.finite(input$fu %||% NA)) input$fu else NA
      )
    })
  })
}
