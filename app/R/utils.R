# utils.R : small shared helpers ---------------------------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

# Time units -> years
time_factor <- function(unit) {
  switch(unit,
    days   = 1 / 365.25,
    weeks  = 7 / 365.25,
    months = 1 / 12,
    years  = 1,
    stop("Unknown time unit: ", unit)
  )
}

# Model time grid: 0, cycle, 2*cycle, ..., horizon
make_grid <- function(horizon, cycle) {
  k <- round(horizon / cycle)
  seq(0, k) * cycle
}

# Per-cycle cumulative hazard from a survival vector on the grid.
# Not clipped at 0: a rising curve gives a negative value, which validation reports.
cycle_hazard <- function(s) {
  ls <- log(pmax(s, 1e-300))
  pmin(-diff(ls), 700)
}

# Survival vector rebuilt from per-cycle cumulative hazards
from_hazard <- function(dh) exp(-c(0, cumsum(dh)))

# Mid-point of each cycle
cycle_mid <- function(t) (t[-length(t)] + t[-1]) / 2

fmt_num <- function(x, digits = 0) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits, big.mark = ","))
}

fmt_pct <- function(x, digits = 1) {
  ifelse(is.na(x), "NA", paste0(formatC(100 * x, format = "f", digits = digits), "%"))
}

# Bootstrap-styled HTML table; optional status column becomes a coloured badge
html_table <- function(df, status_col = NULL, digits = NULL) {
  if (!nrow(df)) return(shiny::p(class = "text-muted", "Nothing to show yet."))
  badge <- function(s) {
    cls <- switch(s, pass = "bg-success", warn = "bg-warning text-dark",
                  fail = "bg-danger", info = "bg-info text-dark", "bg-secondary")
    shiny::span(class = paste("badge", cls), s)
  }
  head_row <- shiny::tags$tr(lapply(names(df), shiny::tags$th))
  rows <- lapply(seq_len(nrow(df)), \(i) {
    shiny::tags$tr(lapply(names(df), \(nm) {
      v <- df[[nm]][i]
      if (!is.null(status_col) && nm == status_col) return(shiny::tags$td(badge(v)))
      if (is.integer(v)) v <- format(v) else if (is.numeric(v)) v <- fmt_num(v, digits %||% 3)
      shiny::tags$td(v)
    }))
  })
  shiny::div(class = "table-responsive",
    shiny::tags$table(class = "table table-sm table-striped align-middle",
      shiny::tags$thead(head_row), shiny::tags$tbody(rows)))
}

# Flatten a nested list to key/value rows (for export)
flatten_spec <- function(x, prefix = "") {
  out <- list()
  for (nm in names(x)) {
    key <- if (nzchar(prefix)) paste(prefix, nm, sep = ".") else nm
    v <- x[[nm]]
    if (is.list(v) && !is.data.frame(v)) {
      out <- c(out, list(flatten_spec(v, key)))
    } else if (is.data.frame(v)) {
      out <- c(out, list(data.frame(key = key, value = paste0("<table, ", nrow(v), " rows>"))))
    } else {
      out <- c(out, list(data.frame(key = key, value = paste(format(v), collapse = "; "))))
    }
  }
  do.call(rbind, out)
}

# Keep the current selection if still valid, otherwise use the default
keep_or <- function(x, choices, default) {
  if (length(x) == 1 && !is.na(x) && x %in% choices) x else default
}
