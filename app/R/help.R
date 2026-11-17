# help.R : in-app guidance (intro cards, info popups, Start here tab, glossary) ----
# Pure interface helpers. They add explanation to the screen and never touch the model.

# A small "i" icon. Clicking it opens a popup with the explanation (works on phones,
# where hover does not exist). Put it in the label of an input, NOT inside a checkbox
# or radio-choice label, because clicking a checkbox label would toggle the checkbox.
help_info <- function(text, title = NULL) {
  bslib::popover(
    shiny::tags$span(
      class = "text-primary ms-1", style = "cursor: pointer;",
      tabindex = "0", role = "button", `aria-label` = "More information", "\u24D8"
    ),
    text,
    title = title,
    placement = "right"
  )
}

# An input label with an info popup next to it
help_label <- function(label, text, title = NULL) {
  shiny::tagList(label, help_info(text, title))
}

# A box at the top of a tab: what it is for, what to do, what to check before moving on.
# Uses plain HTML <details>, so users can collapse it without any JavaScript.
intro_card <- function(purpose, do, check, next_step = NULL) {
  row <- function(head, body) shiny::tags$p(class = "mb-1", shiny::tags$strong(head), " ", body)
  bslib::card(
    class = "mb-3 border-start border-4 border-primary bg-light",
    shiny::tags$details(
      open = NA,
      shiny::tags$summary(class = "fw-semibold", "About this tab (click to hide or show)"),
      shiny::div(
        class = "mt-2 small",
        row("What it is for:", purpose),
        row("What to do here:", do),
        row("Check before moving on:", check),
        if (!is.null(next_step)) row("Next:", next_step)
      )
    )
  )
}

# Glossary: one place for the terms that appear across the app ----------------------
glossary_terms <- function() {
  list(
    c("Partitioned survival model (PSM)",
      "A model that splits each patient's time into three states: progression-free, progressed, and dead. The share of patients in each state over time comes from survival curves, not from transition probabilities."),
    c("Kaplan-Meier (KM) curve",
      "The observed survival curve from trial data. Models are fitted to it and judged against it."),
    c("Censored",
      "A patient who had not had the event when follow-up ended. Their true time is longer than what was observed. In the data file, event = 0."),
    c("OS, PFS, ToT",
      "Overall survival (time to death), progression-free survival (time to progression or death), and time on treatment."),
    c("Extrapolation",
      "Extending a fitted curve past the end of the trial data, up to the model time horizon. Most of the uncertainty in a PSM sits here."),
    c("AIC and BIC",
      "Scores that rank how well different survival models fit the same data. Lower is better. They only describe fit to the observed data, not whether the long-term prediction is believable."),
    c("Hazard and hazard ratio (HR)",
      "The hazard is the instantaneous risk of the event at a given time. An HR below 1 means the intervention lowers that risk compared with the comparator."),
    c("Proportional hazards (PH)",
      "The assumption that the HR stays constant over time. The PH diagnostics tab helps to judge whether it holds."),
    c("Treatment-effect waning",
      "The benefit of the intervention fading over time. The app offers several ways to model it, because the choice can change the result a lot."),
    c("Frailty",
      "Unmeasured differences between patients. When the intervention helps frailer patients survive, the people left in that arm are frailer on average, which is why the population-level effect can fade faster than the individual-level effect."),
    c("Landmark survival",
      "The share of patients alive at a chosen time point, for example 5 years. Easy to compare with clinical expectations."),
    c("QALY",
      "Quality-adjusted life year: one year of life weighted by quality of life (utility), where 1 is full health and 0 is dead."),
    c("ICER",
      "Incremental cost-effectiveness ratio: extra cost divided by extra QALYs of the intervention versus the comparator."),
    c("Willingness to pay (WTP) and NMB",
      "WTP is the most the decision-maker will pay for one QALY. Net monetary benefit = QALYs x WTP minus cost. Positive NMB means the intervention is good value at that WTP."),
    c("Half-cycle correction",
      "Assumes events happen in the middle of each model cycle, on average, which reduces bias from using discrete cycles."),
    c("Relative dose intensity (RDI)",
      "The share of the planned dose that patients actually receive, from 0 to 1."),
    c("DSA",
      "Deterministic sensitivity analysis: change one input at a time across a plausible range and see how the result moves. Shown as a tornado chart."),
    c("PSA",
      "Probabilistic sensitivity analysis: change all uncertain inputs at once, many times, by random draws, to see the spread of results."),
    c("Cost-effectiveness plane and CEAC",
      "The plane plots each PSA run as extra QALYs against extra cost. The CEAC (acceptability curve) shows the chance that the intervention is cost-effective at each WTP value."),
    c("Bootstrap",
      "Re-sampling the patients many times and refitting the models each time, to capture uncertainty in the survival curves."),
    c("Hybrid curve",
      "Uses the observed KM curve up to a cut-point and the fitted model after it."),
    c("Guyot reconstruction",
      "A method that rebuilds approximate patient-level data from a published KM figure plus its numbers at risk. Useful when you have no patient data.")
  )
}

glossary_ui <- function() {
  terms <- glossary_terms()
  bslib::accordion(
    open = FALSE,
    bslib::accordion_panel(
      "Glossary of terms used in this app",
      shiny::tags$dl(
        lapply(terms, function(t) shiny::tagList(shiny::tags$dt(t[1]), shiny::tags$dd(t[2])))
      )
    )
  )
}

# The "Start here" tab ---------------------------------------------------------------
# `go` buttons are wired up in app.R: each one jumps to a tab through bslib::nav_select().
start_step <- function(num, tab, title, text, input_id, go_label = NULL) {
  shiny::tags$tr(
    shiny::tags$td(class = "align-top fw-bold text-primary", num),
    shiny::tags$td(class = "align-top", shiny::tags$strong(title), shiny::tags$br(),
                   shiny::tags$span(class = "small", text)),
    shiny::tags$td(class = "align-top text-end",
                   shiny::actionButton(input_id, go_label %||% paste("Go to", tab),
                                       class = "btn-sm btn-outline-primary"))
  )
}

start_here_ui <- function() {
  shiny::tagList(
    bslib::card(
      bslib::card_header("Welcome"),
      shiny::p(shiny::HTML(paste0(
        "This calculator estimates the costs, health gains and cost-effectiveness of an intervention ",
        "against a comparator, using a three-state partitioned survival model: ",
        "<em>progression-free</em>, <em>progressed</em> and <em>dead</em>."))),
      shiny::p(shiny::HTML(paste0(
        "<strong>Fastest way to see it work:</strong> leave everything as it is, then open ",
        "<strong>Results &gt; Base case</strong>. The built-in example is already filled in."))),
      shiny::p(class = "text-danger",
               shiny::strong("All example data, the life table, and the default costs and utilities are synthetic or illustrative."),
               "Replace them with your own evidence before using any result.")
    ),
    bslib::card(
      bslib::card_header("Work through the tabs in this order"),
      shiny::tags$table(
        class = "table table-sm align-middle mb-0",
        shiny::tags$tbody(
          start_step("1", "Settings", "Set the time horizon, discounting and population",
                     "How long the model runs, how often it ticks (weekly or monthly), discount rates, willingness to pay per QALY, and the age and sex of patients.",
                     "go_settings"),
          start_step("2", "Data", "Load patient data and choose the two arms",
                     "Use the synthetic example or upload your own CSV. Pick the comparator and the intervention. No patient data? Rebuild it from a published curve at the bottom of this tab.",
                     "go_data"),
          start_step("3", "Overall survival", "Choose survival models for OS, then PFS",
                     "Pick how each arm's curve is fitted, compare models, and check the curves against the observed data. Repeat for progression-free survival.",
                     "go_os", "Go to Overall survival"),
          start_step("4", "Clinical plausibility", "Check the long-term predictions make sense",
                     "Look at survival at landmark years, compare with clinical expectations, and see if the modelled risk of death ever drops below the general population.",
                     "go_plaus"),
          start_step("5", "Treatment effects", "Decide how long the benefit lasts",
                     "Choose whether the treatment effect lasts for life or fades (wanes), and say how crossover to the intervention is handled.",
                     "go_effects"),
          start_step("6", "Time on treatment", "Set how long patients stay on treatment",
                     "Drug costs follow time on treatment, not progression-free survival. Set any maximum duration here.",
                     "go_tot"),
          start_step("7", "Costs and utilities", "Enter costs and quality of life",
                     "Drug, administration, care and adverse-event costs for each arm, and utilities for the progression-free and progressed states.",
                     "go_costs"),
          start_step("8", "Results", "Read the base case, then test it",
                     "Base case gives costs, QALYs and the ICER. Scenarios and DSA test one change at a time. PSA tests everything at once. Export saves it all to Excel.",
                     "go_results", "Go to Base case"),
          start_step("9", "Validation", "Check the warnings before you trust a result",
                     "Automated checks on the code, the structure, and the assumptions that still need evidence.",
                     "go_valid")
        )
      )
    ),
    bslib::layout_columns(
      col_widths = c(6, 6),
      bslib::card(
        bslib::card_header("Tips for using it in the browser"),
        shiny::tags$ul(
          class = "small mb-0",
          shiny::tags$li("The first load is slow because R and its packages download into your browser. After that it is faster."),
          shiny::tags$li("Your data stays on your own device. Nothing is sent to a server."),
          shiny::tags$li("Refreshing the page resets all inputs. Download the Excel export before you leave."),
          shiny::tags$li("The PSA can be slow here, because a browser cannot run parallel workers. Try 100 to 200 runs first."),
          shiny::tags$li(shiny::HTML("Look for the <span class='text-primary'>\u24D8</span> icon next to an input. Click it for a plain-language explanation."))
        )
      ),
      bslib::card(
        bslib::card_header("Your CSV file"),
        shiny::tags$ul(
          class = "small mb-0",
          shiny::tags$li("Columns: ", shiny::code("arm"), ", ", shiny::code("endpoint"), ", ",
                         shiny::code("time"), ", ", shiny::code("event"),
                         ", and optionally ", shiny::code("id"), "."),
          shiny::tags$li(shiny::code("endpoint"), " must be OS, PFS or ToT."),
          shiny::tags$li(shiny::code("event"), ": 1 = event happened, 0 = censored."),
          shiny::tags$li("Tell the app the time unit (days, weeks, months or years) on the Data tab."),
          shiny::tags$li("Not sure of the format? Use ", shiny::strong("Download example CSV"), " on the Data tab.")
        )
      )
    ),
    glossary_ui()
  )
}
