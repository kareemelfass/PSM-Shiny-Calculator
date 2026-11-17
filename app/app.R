# Partitioned Survival Model (PSM) app, V3 ---------------------------------------
# V2 added survival-curve uncertainty to the PSA. V3 reworks treatment-effect waning (Doc 12).
# Three-state PSM (progression-free, progressed, dead) for HTA submissions.
# Run: shiny::runApp("path/to/PSM_Shiny_V3"), or open app.R in RStudio and click Run App.

library(shiny)
library(bslib)
library(survival)

# Find this app's folder however the file is run (runApp, source(), or line by line)
find_app_dir <- function() {
  has_code <- \(d) !is.null(d) && length(d) == 1 && file.exists(file.path(d, "R", "fct_engine.R"))
  candidates <- list(
    tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = \(e) NULL),     # source("app.R")
    getwd(),                                                                    # runApp() sets this
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      tryCatch(dirname(rstudioapi::getSourceEditorContext()$path), error = \(e) NULL)  # RStudio, line by line
    }
  )
  for (d in candidates) if (has_code(d)) return(normalizePath(d))
  stop("Cannot find the app's R/ folder. Run shiny::runApp(\"path/to/PSM_Shiny_V3\"), ",
       "or setwd() to the PSM_Shiny_V3 folder first.", call. = FALSE)
}
app_dir <- find_app_dir()

# Load the app's functions explicitly, so they exist even without Shiny's autoload
for (f in sort(list.files(file.path(app_dir, "R"), pattern = "\\.[Rr]$", full.names = TRUE))) {
  source(f, local = FALSE)
}

example_dir <- file.path(app_dir, "example")

ui <- page_navbar(
  title = "Partitioned survival model V3",
  theme = bs_theme(version = 5, preset = "shiny", primary = "#3E5C8A"),
  fillable = FALSE,
  nav_panel("Settings", mod_settings_ui("settings")),
  nav_panel("Data", mod_data_ui("data")),
  nav_menu("Survival",
    nav_panel("Overall survival", mod_endpoint_ui("os", "OS")),
    nav_panel("Progression-free survival", mod_endpoint_ui("pfs", "PFS")),
    nav_panel("Clinical plausibility", mod_plausibility_ui("plaus"))
  ),
  nav_panel("Treatment effects", mod_effects_ui("effects")),
  nav_panel("Time on treatment", mod_tot_ui("tot")),
  nav_panel("Costs and utilities", mod_costs_ui("costs")),
  nav_menu("Results",
    nav_panel("Base case", mod_results_ui("results")),
    nav_panel("Scenarios and DSA", mod_scenarios_ui("scen")),
    nav_panel("PSA", mod_psa_ui("psa")),
    nav_panel("Export", mod_export_ui("export"))
  ),
  nav_panel("Validation", mod_validation_ui("valid")),
  nav_spacer(),
  nav_item(tags$span(class = "small text-muted",
                     paste("Fitting engine:", fit_engine())))
)

server <- function(input, output, session) {
  settings <- mod_settings_server("settings", example_dir)
  data <- mod_data_server("data", example_dir)
  os <- mod_endpoint_server("os", "OS", data, settings)
  pfs <- mod_endpoint_server("pfs", "PFS", data, settings)

  # Curves are defined before the modules that both use them and feed them;
  # reactives are lazy, so `effects` and `tot` exist by the time this runs.
  fits <- reactive(list(OS = os$fits(), PFS = pfs$fits(), ToT = tot$fits()))

  spec <- reactive({
    st <- settings()
    st$fu_years <- if (is.finite(st$fu_override)) st$fu_override else data$fu()
    ef <- effects$spec()
    list(
      arms = data$arms(),
      settings = st,
      endpoints = list(OS = os$spec(), PFS = pfs$spec()),
      waning = ef$waning,
      tot = tot$spec(),
      econ = econ(),
      crossover = ef$crossover
    )
  })

  pieces <- reactive({
    sp <- spec()
    t <- make_grid(sp$settings$horizon, sp$settings$cycle)
    tryCatch(predict_pieces(sp, fits(), t),
             error = \(e) validate(need(FALSE, conditionMessage(e))))
  })

  model <- reactive({
    sp <- spec()
    pc <- pieces()
    cv <- assemble_curves(pc, sp)
    res <- run_engine(cv, sp)
    list(spec = sp, pieces = pc, curves = cv, res = res,
         validation = validate_model(sp, pc, cv, res, fits()))
  })
  curves <- reactive(model()$curves)

  effects <- mod_effects_server("effects", model, data$arms, fits)
  tot <- mod_tot_server("tot", data, curves)
  econ <- mod_costs_server("costs", data$arms)
  mod_plausibility_server("plaus", curves, os$fits, data$arms, settings)
  mod_results_server("results", model)
  sc <- mod_scenarios_server("scen", model)
  ps <- mod_psa_server("psa", model, fits, data$ipd, file.path(app_dir, "R"))
  mod_validation_server("valid", model, ps$summary, ps$compare)
  mod_export_server("export", model, fits, sc$scenarios, sc$dsa, ps$runs, ps$bands, ps$compare,
                    effects$explorer)
}

shinyApp(ui, server)
