# mod_settings.R : model settings ----------------------------------------------

mod_settings_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_columns(
    col_widths = c(4, 4, 4),
    bslib::card(
      bslib::card_header("Time"),
      shiny::numericInput(ns("horizon"), "Time horizon (years)", 30, min = 1, max = 100),
      shiny::selectInput(ns("cycle"), "Cycle length",
                         c("Weekly" = "week", "Monthly" = "month"), "week"),
      shiny::checkboxInput(ns("hcc"), "Half-cycle correction (trapezoid)", TRUE),
      shiny::numericInput(ns("fu"), "Trial follow-up (years); blank = longest OS time in data",
                          NA, min = 0)
    ),
    bslib::card(
      bslib::card_header("Economics"),
      shiny::numericInput(ns("dr_cost"), "Discount rate, costs (%)", 3.5, min = 0, step = 0.5),
      shiny::numericInput(ns("dr_out"), "Discount rate, outcomes (%)", 3.5, min = 0, step = 0.5),
      shiny::textInput(ns("currency"), "Currency label", "CU"),
      shiny::numericInput(ns("wtp"), "Willingness to pay per QALY", 30000, min = 0, step = 1000),
      shiny::helpText("All default inputs are illustrative. Replace them with your own evidence.")
    ),
    bslib::card(
      bslib::card_header("Population and mortality"),
      shiny::numericInput(ns("age"), "Starting age (years)", 60, min = 18, max = 99),
      shiny::sliderInput(ns("male"), "Male (%)", 0, 100, 50),
      shiny::checkboxInput(ns("floor_on"), "Apply general-population mortality floor", TRUE),
      shiny::radioButtons(ns("lt_src"), "Life table",
                          c("Illustrative (replace for real use)" = "example", "Upload" = "upload")),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'upload'", ns("lt_src")),
        shiny::fileInput(ns("lt_file"), "CSV: age, qx_male, qx_female", accept = ".csv")
      ),
      shiny::checkboxInput(ns("pfs_le_os"), "Constrain PFS to OS (reported in Validation)", TRUE)
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
