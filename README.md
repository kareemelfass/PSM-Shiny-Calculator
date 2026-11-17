# PSM Shiny Calculator

Three-state partitioned survival model (progression-free, progressed, dead) for HTA work, as an interactive calculator.

**Live app:** https://kareemelfass.github.io/PSM-Shiny-Calculator/

The app runs entirely in your browser (R compiled to WebAssembly through [Shinylive](https://posit-dev.github.io/r-shinylive/)). Nothing is uploaded to a server. First load takes a while because R and its packages are downloaded.

All example data, the life table and default costs and utilities are synthetic or illustrative.

## Built-in guidance

- **Start here tab:** the workflow in 9 numbered steps, with a button that jumps to each tab, plus tips and a glossary.
- **"About this tab" box** at the top of every tab: what it is for, what to do, what to check, and what comes next. Click to hide it.
- **Info icons** next to inputs: click for a plain-language explanation.
- All of this lives in `app/R/help.R`. It is interface text only and never changes the model.

## Repo layout

- `app/` : the Shiny app (see `app/README.md` for the model documentation)
- `.github/workflows/deploy.yml` : builds the static site and publishes it to GitHub Pages on every push to `main`

## Run locally with full R speed

```r
install.packages(c("shiny", "bslib", "survival", "flexsurv", "ggplot2", "writexl"))
shiny::runApp("app")
```

## Browser limits

- Parallel workers do not exist in the browser, so the PSA bootstrap runs on one core and is slow. Use fewer runs, or run locally.
- If `flexsurv` is not available in the browser build, the app falls back to the `survreg` engine (fewer distributions). The engine in use is shown top right.
