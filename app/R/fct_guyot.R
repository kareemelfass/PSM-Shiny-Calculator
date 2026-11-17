# fct_guyot.R : reconstruct IPD from a digitised KM curve ----------------------
# Guyot P, Ades AE, Ouwens MJ, Welton NJ. BMC Med Res Methodol 2012;12:9.
# Follows the published algorithm; loops are capped so bad input cannot hang the app.

# Counts of x in the right-closed bins (b[1], b[2]], (b[2], b[3]], ...
bin_counts <- function(x, breaks) {
  tabulate(findInterval(x, breaks, left.open = TRUE), nbins = length(breaks) - 1)
}

# Tidy digitised points: sort, survival as a proportion, start at (0, 1),
# non-increasing, one point per time (lowest survival at a vertical drop)
prepare_digitised <- function(digi) {
  d <- data.frame(time = as.numeric(digi[[1]]), surv = as.numeric(digi[[2]]))
  d <- d[is.finite(d$time) & is.finite(d$surv) & d$time >= 0, ]
  if (max(d$surv) > 1.5) d$surv <- d$surv / 100
  d <- d[order(d$time, -d$surv), ]
  d$surv <- pmin(pmax(d$surv, 0), 1)
  n_fixed <- sum(diff(d$surv) > 0)
  d$surv <- cummin(d$surv)
  d <- stats::aggregate(surv ~ time, data = d, FUN = min)
  if (d$time[1] > 0) d <- rbind(data.frame(time = 0, surv = 1), d)
  d$surv[d$time == 0] <- 1
  attr(d, "n_monotone_fixed") <- n_fixed
  d
}

guyot_reconstruct <- function(digi, nrisk, tot_events = NA, max_iter = 1000) {
  d <- prepare_digitised(digi)
  t_S <- d$time
  S <- d$surv
  r <- data.frame(time = as.numeric(nrisk[[1]]), n = as.numeric(nrisk[[2]]))
  r <- r[order(r$time), ]
  r <- r[r$time <= max(t_S), ]
  if (r$time[1] != 0) stop("The numbers-at-risk table must start at time 0.")

  t_risk <- r$time
  n_risk <- round(r$n)
  n_int <- length(n_risk)
  lower <- vapply(t_risk, \(x) min(which(t_S >= x)), integer(1))
  upper <- c(lower[-1] - 1L, length(t_S))
  if (any(upper < lower)) {
    stop("Each numbers-at-risk interval needs at least one digitised point. ",
         "Check intervals starting at: ", paste(t_risk[upper < lower], collapse = ", "))
  }
  n_t <- upper[n_int]

  n_censor <- rep(0, n_int)
  n_hat <- rep(n_risk[1] + 1, n_t + 1)
  cen <- rep(0, n_t)
  d_ev <- rep(0, n_t)
  km_hat <- rep(1, n_t)
  last_i <- rep(1, n_int)
  warn <- character()

  spread_censor <- function(i, ncen, from, to) {
    if (ncen <= 0) return(rep(0, to - from + 1))
    ct <- t_S[from] + seq_len(ncen) * (t_S[to + 1] - t_S[from]) / (ncen + 1)
    bin_counts(ct, t_S[from:(to + 1)])
  }

  if (n_int > 1) {
    for (i in 1:(n_int - 1)) {
      n_censor[i] <- round(n_risk[i] * S[lower[i + 1]] / S[lower[i]] - n_risk[i + 1])
      iter <- 0
      last <- last_i[i]
      repeat {
        iter <- iter + 1
        if (n_censor[i] <= 0) {
          cen[lower[i]:upper[i]] <- 0
          n_censor[i] <- 0
        } else {
          cen[lower[i]:upper[i]] <- spread_censor(i, n_censor[i], lower[i], upper[i])
        }
        n_hat[lower[i]] <- n_risk[i]
        last <- last_i[i]
        for (k in lower[i]:upper[i]) {
          if (i == 1 && k == lower[i]) {
            d_ev[k] <- 0
            km_hat[k] <- 1
          } else {
            d_ev[k] <- if (n_hat[k] > 0) round(n_hat[k] * (1 - S[k] / km_hat[last])) else 0
            km_hat[k] <- km_hat[last] * (1 - d_ev[k] / max(n_hat[k], 1))
          }
          n_hat[k + 1] <- n_hat[k] - d_ev[k] - cen[k]
          if (d_ev[k] != 0) last <- k
        }
        n_censor[i] <- n_censor[i] + (n_hat[lower[i + 1]] - n_risk[i + 1])
        done <- !(n_hat[lower[i + 1]] > n_risk[i + 1] ||
                  (n_hat[lower[i + 1]] < n_risk[i + 1] && n_censor[i] > 0))
        if (done) break
        if (iter >= max_iter) {
          warn <- c(warn, paste("Interval", i, "did not fully converge"))
          break
        }
      }
      if (n_hat[lower[i + 1]] < n_risk[i + 1]) n_risk[i + 1] <- n_hat[lower[i + 1]]
      last_i[i + 1] <- last
    }
  }

  # Last interval: censoring at the average earlier rate
  li <- n_int
  n_censor[li] <- if (n_int > 1) {
    min(round(sum(n_censor[1:(n_int - 1)]) * (t_S[upper[li]] - t_S[lower[li]]) /
                (t_S[upper[li - 1]] - t_S[lower[1]])), n_risk[li])
  } else 0

  run_last <- function() {
    if (upper[li] > lower[li]) {
      idx <- lower[li]:(upper[li] - 1)
      cen[idx] <<- if (n_censor[li] > 0) {
        ct <- t_S[lower[li]] + seq_len(n_censor[li]) *
          (t_S[upper[li]] - t_S[lower[li]]) / (n_censor[li] + 1)
        bin_counts(ct, t_S[lower[li]:upper[li]])
      } else 0
    }
    n_hat[lower[li]] <<- n_risk[li]
    last <- last_i[li]
    for (k in lower[li]:upper[li]) {
      d_ev[k] <<- if (km_hat[last] != 0 && n_hat[k] > 0) {
        round(n_hat[k] * (1 - S[k] / km_hat[last]))
      } else 0
      km_hat[k] <<- km_hat[last] * (1 - d_ev[k] / max(n_hat[k], 1))
      n_hat[k + 1] <<- n_hat[k] - d_ev[k] - cen[k]
      if (n_hat[k + 1] < 0) {
        n_hat[k + 1] <<- 0
        cen[k] <<- n_hat[k] - d_ev[k]
      }
      if (d_ev[k] != 0) last <- k
    }
  }
  if (n_censor[li] < 0) n_censor[li] <- 0
  run_last()

  # Optional: adjust last-interval censoring so total events match the report
  if (is.finite(tot_events)) {
    sum_before <- if (n_int > 1) sum(d_ev[1:upper[n_int - 1]]) else 0
    if (n_int > 1 && sum_before >= tot_events) {
      # As published: no further events or censoring in the last interval
      d_ev[lower[li]:upper[li]] <- 0
      if (upper[li] > lower[li]) cen[lower[li]:(upper[li] - 1)] <- 0
      warn <- c(warn, "Reported total events already reached before the last interval")
    } else {
      iter <- 0
      sumd <- sum(d_ev[1:upper[li]])
      while ((sumd > tot_events || (sumd < tot_events && n_censor[li] > 0)) && iter < max_iter) {
        iter <- iter + 1
        n_censor[li] <- max(n_censor[li] + (sumd - tot_events), 0)
        run_last()
        sumd <- sum(d_ev[1:upper[li]])
      }
      if (iter >= max_iter) warn <- c(warn, "Total-events adjustment did not converge")
    }
  }

  # Write individual records: events at digitised times, censoring at midpoints
  ev_t <- rep(t_S[seq_len(n_t)], times = pmax(d_ev, 0))
  cen_n <- pmax(cen[seq_len(n_t - 1)], 0)
  cen_t <- rep((t_S[seq_len(n_t - 1)] + t_S[2:n_t]) / 2, times = cen_n)
  n_left <- n_risk[1] - length(ev_t) - length(cen_t)
  ipd <- data.frame(
    time = c(ev_t, cen_t, rep(t_S[n_t], max(n_left, 0))),
    event = c(rep(1L, length(ev_t)), rep(0L, length(cen_t) + max(n_left, 0)))
  )
  if (n_left < 0) warn <- c(warn, "More reconstructed records than patients at risk at time 0")
  attr(ipd, "warnings") <- warn
  attr(ipd, "n_monotone_fixed") <- attr(d, "n_monotone_fixed")
  ipd[order(ipd$time), ]
}
