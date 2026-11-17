# mod_data.R : patient data, arms, and Guyot reconstruction ---------------------

mod_data_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 320,
      shiny::radioButtons(ns("src"), "Data source",
                          c("Synthetic example" = "example", "Upload CSV" = "upload")),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'upload'", ns("src")),
        shiny::fileInput(ns("file"), "CSV columns: arm, endpoint, time, event (+ optional id)", accept = ".csv"),
        shiny::selectInput(ns("unit"), help_label("Time unit in file", "The unit your time column uses. The app converts everything to years."),
                           c("days", "weeks", "months", "years"), "months")
      ),
      shiny::selectInput(ns("comp"), help_label("Comparator arm", "The arm you compare against, for example standard care."), NULL),
      shiny::selectInput(ns("int"), help_label("Intervention arm", "The new treatment being assessed. Choose (no data: use HR) if you have no patient data for it and will apply a hazard ratio instead, for example from a network meta-analysis."), NULL),
      shiny::helpText("endpoint must be OS, PFS or ToT. event: 1 = event, 0 = censored.",
                      "Optional 'id' column: links each patient's OS, PFS and ToT rows (used by the PSA bootstrap).",
                      "The intervention arm may have no data if you will apply an HR (e.g. from an NMA)."),
      shiny::downloadButton(ns("dl_example"), "Download example CSV", class = "btn-sm")
    ),
    intro_card(
      purpose = "Load the patient-level survival data, and tell the app which arm is the comparator and which is the intervention.",
      do = "Use the synthetic example, or choose Upload CSV and set the time unit. Then pick the comparator and intervention arms in the left panel. If you only have a published curve and no patient data, use the reconstruction box at the bottom of this tab.",
      check = "The arm names, patient counts and event counts look right, and the Kaplan-Meier curves look like the trial you expect.",
      next_step = "Overall survival, then Progression-free survival."
    ),
    bslib::layout_columns(
      col_widths = c(5, 7),
      bslib::card(bslib::card_header("Data summary"), shiny::uiOutput(ns("summary"))),
      bslib::card(
        bslib::card_header("Kaplan-Meier"),
        shiny::radioButtons(ns("km_ep"), NULL, c("OS", "PFS", "ToT"), inline = TRUE),
        shiny::plotOutput(ns("km_plot"), height = 330)
      )
    ),
    bslib::card(
      bslib::card_header("Reconstruct patient data from a published KM curve (Guyot 2012)"),
      shiny::helpText("Use this when you only have a published figure. Steps: (1) Digitise the curve with a plot-digitising tool and save time and survival, with survival in percent as in the example file. (2) Save the numbers-at-risk table as time and n. (3) Upload both files, or tick the example box to see how it works. (4) Set the arm name, endpoint and time unit. (5) Click Reconstruct and check the plot. (6) Click Add to dataset."),
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        shiny::div(
          shiny::fileInput(ns("g_digi"), "Digitised curve CSV (time, survival)", accept = ".csv"),
          shiny::fileInput(ns("g_risk"), "Numbers at risk CSV (time, n)", accept = ".csv"),
          shiny::checkboxInput(ns("g_example"), "Use the example digitised files", FALSE)
        ),
        shiny::div(
          shiny::textInput(ns("g_arm"), "Arm name", "Reconstructed arm"),
          shiny::selectInput(ns("g_ep"), "Endpoint", c("OS", "PFS", "ToT")),
          shiny::selectInput(ns("g_unit"), "Time unit of the files",
                             c("days", "weeks", "months", "years"), "months"),
          shiny::numericInput(ns("g_events"), "Total events reported (optional)", NA, min = 0)
        ),
        shiny::div(
          shiny::actionButton(ns("g_run"), "Reconstruct", class = "btn-primary"),
          shiny::actionButton(ns("g_add"), "Add to dataset"),
          shiny::uiOutput(ns("g_msg"))
        )
      ),
      shiny::plotOutput(ns("g_plot"), height = 300)
    )
  )
}

mod_data_server <- function(id, example_dir) {
  shiny::moduleServer(id, function(input, output, session) {
    example_csv <- file.path(example_dir, "example_ipd_months.csv")
    added <- shiny::reactiveVal(NULL)

    base <- shiny::reactive({
      if (input$src == "upload") {
        shiny::req(input$file)
        tryCatch(read_ipd(input$file$datapath, input$unit),
                 error = \(e) shiny::validate(shiny::need(FALSE, conditionMessage(e))))
      } else {
        read_ipd(example_csv, "months")
      }
    })
    ipd <- shiny::reactive(rbind(base(), added()))

    shiny::observe({
      a <- unique(ipd()$arm)
      cur_comp <- shiny::isolate(input$comp)
      cur_int <- shiny::isolate(input$int)
      comp_sel <- if (isTRUE(cur_comp %in% a)) cur_comp else if ("Standard care" %in% a) "Standard care" else a[1]
      int_sel <- if (isTRUE(cur_int %in% c(a, "(no data: use HR)"))) cur_int else setdiff(a, comp_sel)[1]
      shiny::updateSelectInput(session, "comp", choices = a, selected = comp_sel)
      shiny::updateSelectInput(session, "int", choices = c(a, "(no data: use HR)"),
                               selected = int_sel %||% "(no data: use HR)")
    })

    arms <- shiny::reactive({
      shiny::req(input$comp, input$int)
      shiny::validate(shiny::need(input$comp != input$int,
                                  "Choose different comparator and intervention arms."))
      int <- if (input$int == "(no data: use HR)") "Intervention" else input$int
      list(comp = input$comp, int = int)
    })

    fu <- shiny::reactive({
      d <- ipd()
      d <- d[d$endpoint == "OS" & d$arm %in% unlist(arms()), ]
      if (nrow(d)) max(d$time) else NA
    })

    output$summary <- shiny::renderUI({
      d <- ipd()
      g <- split(d, list(d$arm, d$endpoint), drop = TRUE)
      tb <- do.call(rbind, lapply(g, \(x) data.frame(
        Arm = x$arm[1], Endpoint = x$endpoint[1], Patients = nrow(x), Events = sum(x$event),
        `Max follow-up (yr)` = round(max(x$time), 2), check.names = FALSE)))
      tb <- tb[order(tb$Arm, tb$Endpoint), ]
      d_arms <- d[d$arm %in% unlist(arms()), ]
      shiny::tagList(
        html_table(tb),
        shiny::p(class = "small text-muted",
          if (ipd_is_paired(d_arms)) "Patient ids found: the PSA bootstrap resamples whole patients (OS, PFS and ToT together)." else
            "No complete patient ids: the PSA bootstrap resamples OS, PFS and ToT separately, so their correlation is not kept. Add an 'id' column to link them.")
      )
    })

    output$km_plot <- shiny::renderPlot({
      d <- ipd()
      d <- d[d$endpoint == input$km_ep, ]
      shiny::validate(shiny::need(nrow(d) > 0, paste("No", input$km_ep, "data.")))
      kms <- do.call(rbind, lapply(split(d, d$arm), \(x) {
        km_df(survival::survfit(survival::Surv(time, event) ~ 1, data = x), x$arm[1])
      }))
      ggplot2::ggplot(kms, ggplot2::aes(time, surv, colour = arm)) +
        ggplot2::geom_step(linewidth = 0.8) + pct_axis() + theme_psm() +
        ggplot2::labs(x = "Years", y = input$km_ep, colour = NULL)
    })

    output$dl_example <- shiny::downloadHandler(
      filename = "example_ipd_months.csv",
      content = \(file) file.copy(example_csv, file)
    )

    # ---- Guyot ----
    recon <- shiny::eventReactive(input$g_run, {
      if (isTRUE(input$g_example)) {
        digi <- utils::read.csv(file.path(example_dir, "example_digitised_os_standard_care.csv"))
        risk <- utils::read.csv(file.path(example_dir, "example_nrisk_os_standard_care.csv"))
      } else {
        shiny::validate(shiny::need(input$g_digi, "Upload the digitised curve."),
                        shiny::need(input$g_risk, "Upload the numbers-at-risk table."))
        digi <- utils::read.csv(input$g_digi$datapath)
        risk <- utils::read.csv(input$g_risk$datapath)
      }
      rec <- tryCatch(guyot_reconstruct(digi, risk, input$g_events %||% NA),
                      error = \(e) shiny::validate(shiny::need(FALSE, conditionMessage(e))))
      list(rec = rec, digi = prepare_digitised(digi))
    })

    output$g_msg <- shiny::renderUI({
      r <- recon()
      w <- attr(r$rec, "warnings")
      fixed <- attr(r$rec, "n_monotone_fixed") %||% 0
      shiny::div(class = "mt-2",
        shiny::p(sprintf("Reconstructed %d patients, %d events.", nrow(r$rec), sum(r$rec$event))),
        if (fixed > 0) shiny::p(class = "text-warning",
                                sprintf("%d digitised points rose and were made non-increasing.", fixed)),
        if (length(w)) shiny::p(class = "text-warning", paste(w, collapse = "; "))
      )
    })

    output$g_plot <- shiny::renderPlot({
      r <- recon()
      km <- survival::survfit(survival::Surv(time, event) ~ 1, data = r$rec)
      ggplot2::ggplot() +
        ggplot2::geom_point(data = r$digi, ggplot2::aes(time, surv), colour = "grey40", size = 1.4) +
        ggplot2::geom_step(data = km_df(km, "rec"), ggplot2::aes(time, surv),
                           colour = "#D9822B", linewidth = 0.9) +
        pct_axis() + theme_psm() +
        ggplot2::labs(x = paste0("Time (", input$g_unit, ")"), y = "Survival",
                      title = "Digitised points (grey) vs KM of reconstructed data (orange)")
    })

    shiny::observeEvent(input$g_add, {
      r <- recon()
      new <- data.frame(arm = input$g_arm, endpoint = input$g_ep,
                        time = r$rec$time * time_factor(input$g_unit), event = r$rec$event,
                        id = NA_character_)   # reconstructed curves cannot be linked across endpoints
      added(rbind(added(), new))
      shiny::showNotification(sprintf("Added %d records for %s, %s.", nrow(new),
                                      input$g_arm, input$g_ep), type = "message")
    })

    list(ipd = ipd, arms = arms, fu = fu)
  })
}
