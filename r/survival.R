#!/usr/bin/env Rscript
# Prognostic modelling that the Python side does not cover.
# Basic KM and single Cox stay in omics_stats (lifelines); this adds
# penalised model building and the validation figures reviewers ask for.
#
# params:
#   method        lasso_cox | timeroc | nomogram | calibration | dca
#   data_path     one row per subject (or `data` records)
#   time, event   column names; event must be 0/1 (1 = event occurred)
#   predictors    candidate variable columns
#   alpha         glmnet mixing: 1 lasso (default), 0 ridge
#   nfolds        cross-validation folds (default 10)
#   lambda        "1se" (default, sparser) or "min"
#   times         evaluation horizons for timeroc/calibration/dca
#   risk_column   pre-computed risk score column, when scoring an existing model

OMICS_R_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(OMICS_R_DIR, "lib", "common.R"))

prepare <- function(params, require_predictors = TRUE) {
  df <- omics_read_table(params, "data_path", "data")
  time_col <- as.character(params$time %||% "time")[1]
  event_col <- as.character(params$event %||% "event")[1]
  for (col in c(time_col, event_col)) {
    if (!col %in% colnames(df)) {
      stop(sprintf("Column '%s' not found. Available: %s", col, paste(colnames(df), collapse = ", ")), call. = FALSE)
    }
  }

  predictors <- as.character(params$predictors %||% character())
  if (require_predictors) {
    if (!length(predictors)) stop("predictors is required.", call. = FALSE)
    missing <- setdiff(predictors, colnames(df))
    if (length(missing)) {
      stop(sprintf("predictors not found in data: %s", paste(missing, collapse = ", ")), call. = FALSE)
    }
  }

  keep <- c(time_col, event_col, predictors, as.character(params$risk_column %||% character()))
  keep <- intersect(unique(keep), colnames(df))
  df <- df[, keep, drop = FALSE]
  df[[time_col]] <- as.numeric(df[[time_col]])
  df[[event_col]] <- as.numeric(df[[event_col]])
  before <- nrow(df)
  df <- df[stats::complete.cases(df), , drop = FALSE]
  dropped <- before - nrow(df)

  if (!all(df[[event_col]] %in% c(0, 1))) {
    stop("event must be coded 0/1 (1 = event occurred).", call. = FALSE)
  }
  if (any(df[[time_col]] <= 0, na.rm = TRUE)) {
    stop("time must be positive; drop or correct non-positive follow-up times.", call. = FALSE)
  }

  events <- sum(df[[event_col]] == 1)
  if (events < 5L) {
    stop(sprintf("Only %d event(s) observed. Prognostic modelling is not interpretable at this event count.", events), call. = FALSE)
  }

  list(df = df, time = time_col, event = event_col, predictors = predictors,
       events = events, n = nrow(df), dropped = dropped)
}

# Events per variable: the standard guard against overfit prognostic models.
epv_warning <- function(events, n_vars) {
  if (n_vars < 1) return(NULL)
  epv <- events / n_vars
  if (epv < 10) {
    sprintf("Events per variable is %.1f (%d events / %d variables). Below 10 the model is prone to overfitting; report this and validate externally.", epv, events, n_vars)
  } else NULL
}

run_lasso_cox <- function(params) {
  omics_require(c("glmnet", "survival"))
  io <- prepare(params)
  x <- stats::model.matrix(
    stats::as.formula(paste("~", paste(io$predictors, collapse = " + "))),
    data = io$df
  )[, -1, drop = FALSE]
  if (ncol(x) < 2L) stop("LASSO needs at least 2 predictor columns after encoding.", call. = FALSE)
  y <- survival::Surv(io$df[[io$time]], io$df[[io$event]])

  alpha <- as.numeric(params$alpha %||% 1)
  nfolds <- as.integer(params$nfolds %||% 10)
  seed <- as.integer(params$seed %||% 42)
  set.seed(seed)
  fit <- glmnet::cv.glmnet(x, y, family = "cox", alpha = alpha, nfolds = nfolds)

  choice <- as.character(params$lambda %||% "1se")[1]
  lambda <- if (identical(choice, "min")) fit$lambda.min else fit$lambda.1se
  coefs <- as.matrix(stats::coef(fit, s = lambda))
  selected <- rownames(coefs)[coefs[, 1] != 0]
  if (!length(selected)) {
    stop(sprintf("No variables retained at lambda.%s. Try lambda='min' or reduce the candidate set.", choice), call. = FALSE)
  }

  # Risk score from the penalised model, then an unpenalised Cox on that single
  # score so the reported HR is interpretable.
  risk <- as.numeric(stats::predict(fit, newx = x, s = lambda, type = "link"))
  io$df$risk_score <- risk
  io$df$risk_group <- ifelse(risk > stats::median(risk), "High", "Low")
  cox <- survival::coxph(
    survival::Surv(io$df[[io$time]], io$df[[io$event]]) ~ risk_score,
    data = io$df
  )
  cindex <- unname(summary(cox)$concordance[1])

  name <- as.character(params$output_name %||% "lasso_cox")[1]
  coef_table <- data.frame(
    variable = selected,
    coefficient = coefs[selected, 1],
    hazard_ratio = exp(coefs[selected, 1]),
    stringsAsFactors = FALSE
  )
  coef_path <- omics_save_table(coef_table, "survival", paste0(name, "_coefficients"))
  score_path <- omics_save_table(
    data.frame(sample = rownames(io$df), risk_score = risk, risk_group = io$df$risk_group, stringsAsFactors = FALSE),
    "survival", paste0(name, "_risk_scores")
  )

  cv_figure <- omics_save_base_plot(
    function() { graphics::plot(fit); graphics::title(main = "Cross-validated partial likelihood deviance", line = 2.5) },
    "survival", paste0(name, "_cv"), width = 6.5, height = 5
  )

  km_figure <- NULL
  if (requireNamespace("survminer", quietly = TRUE)) {
    surv_df <- io$df
    surv_df$.time <- surv_df[[io$time]]
    surv_df$.event <- surv_df[[io$event]]
    sfit <- survival::survfit(survival::Surv(.time, .event) ~ risk_group, data = surv_df)
    gg <- survminer::ggsurvplot(
      sfit, data = surv_df, pval = TRUE, risk.table = TRUE, conf.int = TRUE,
      palette = c("#b5483a", "#2f5f8f"), legend.title = "Risk",
      xlab = "Time", ylab = "Survival probability"
    )
    km_figure <- omics_save_base_plot(
      function() print(gg), "survival", paste0(name, "_km"), width = 7, height = 6.5
    )
  }

  list(
    method = "lasso_cox",
    alpha = alpha,
    lambda_rule = choice,
    lambda = lambda,
    seed = seed,
    n = io$n,
    events = io$events,
    rows_dropped_missing = io$dropped,
    candidates = length(io$predictors),
    selected_variables = selected,
    coefficients = coef_table,
    concordance_index = cindex,
    risk_group_cut = "median risk score",
    coefficient_table = coef_path,
    risk_score_table = score_path,
    cv_figure = cv_figure,
    km_figure = km_figure,
    warnings = c(
      epv_warning(io$events, length(selected)),
      "The C-index above is computed on the same data used to fit the model and is optimistic. Report an external or held-out validation C-index before claiming prognostic value."
    )
  )
}

run_timeroc <- function(params) {
  omics_require(c("timeROC"))
  risk_column <- as.character(params$risk_column %||% "risk_score")[1]
  params$predictors <- params$predictors %||% risk_column
  io <- prepare(params, require_predictors = FALSE)
  if (!risk_column %in% colnames(io$df)) {
    stop(sprintf("risk_column '%s' not found. Run lasso_cox first or supply a score column.", risk_column), call. = FALSE)
  }

  times <- as.numeric(params$times %||% stats::quantile(io$df[[io$time]], c(0.25, 0.5, 0.75), names = FALSE))
  roc <- timeROC::timeROC(
    T = io$df[[io$time]], delta = io$df[[io$event]], marker = io$df[[risk_column]],
    cause = 1, times = times, iid = TRUE
  )
  auc <- as.numeric(roc$AUC)
  ci <- tryCatch(stats::confint(roc)$CI_AUC / 100, error = function(e) NULL)

  name <- as.character(params$output_name %||% "timeroc")[1]
  figure <- omics_save_base_plot(
    function() {
      colours <- c("#2f5f8f", "#b5483a", "#5f8f2f", "#8f6f2f")
      for (i in seq_along(times)) {
        graphics::plot(roc, time = times[i], col = colours[(i - 1) %% length(colours) + 1], add = i > 1, title = FALSE)
      }
      graphics::legend(
        "bottomright", bty = "n",
        legend = sprintf("t = %.3g (AUC %.3f)", times, auc),
        col = colours[(seq_along(times) - 1) %% length(colours) + 1], lty = 1, lwd = 2
      )
      graphics::title(main = "Time-dependent ROC")
    },
    "survival", name, width = 6, height = 6
  )

  list(
    method = "timeroc",
    marker = risk_column,
    n = io$n,
    events = io$events,
    times = times,
    auc = auc,
    auc_ci = if (!is.null(ci)) as.data.frame(ci) else NULL,
    figure = figure,
    notes = "AUC computed on the data supplied; if this is the training set the value is optimistic."
  )
}

run_nomogram <- function(params) {
  omics_require(c("rms", "survival"))
  io <- prepare(params)
  df <- io$df
  df$.time <- df[[io$time]]
  df$.event <- df[[io$event]]
  dd <- rms::datadist(df)
  old <- options(datadist = dd)
  on.exit(options(old), add = TRUE)

  formula <- stats::as.formula(paste("survival::Surv(.time, .event) ~", paste(io$predictors, collapse = " + ")))
  fit <- rms::cph(formula, data = df, x = TRUE, y = TRUE, surv = TRUE)
  times <- as.numeric(params$times %||% stats::quantile(df$.time, c(0.5, 0.75), names = FALSE))
  surv_fun <- rms::Survival(fit)
  funs <- lapply(times, function(t) function(lp) surv_fun(t, lp))
  names(funs) <- sprintf("%.3g-time survival", times)

  nom <- rms::nomogram(fit, fun = funs, lp = FALSE, funlabel = names(funs))
  name <- as.character(params$output_name %||% "nomogram")[1]
  figure <- omics_save_base_plot(
    function() graphics::plot(nom, xfrac = 0.35),
    "survival", name, width = as.numeric(params$width %||% 9), height = as.numeric(params$height %||% 7)
  )

  list(
    method = "nomogram",
    n = io$n,
    events = io$events,
    predictors = io$predictors,
    times = times,
    concordance_index = unname(fit$stats["Dxy"] / 2 + 0.5),
    figure = figure,
    warnings = c(
      epv_warning(io$events, length(io$predictors)),
      "A nomogram displays the fitted model; it is not itself validation. Pair it with a calibration curve and external validation."
    )
  )
}

run_calibration <- function(params) {
  omics_require(c("rms", "survival"))
  io <- prepare(params)
  df <- io$df
  df$.time <- df[[io$time]]
  df$.event <- df[[io$event]]
  dd <- rms::datadist(df)
  old <- options(datadist = dd)
  on.exit(options(old), add = TRUE)

  horizon <- as.numeric((params$times %||% stats::median(df$.time))[1])
  formula <- stats::as.formula(paste("survival::Surv(.time, .event) ~", paste(io$predictors, collapse = " + ")))
  units <- as.integer(params$groups %||% max(20, floor(io$n / 5)))
  fit <- rms::cph(formula, data = df, x = TRUE, y = TRUE, surv = TRUE, time.inc = horizon)
  cal <- rms::calibrate(fit, cmethod = "KM", method = "boot", u = horizon, m = units, B = as.integer(params$bootstrap %||% 200))

  name <- as.character(params$output_name %||% "calibration")[1]
  figure <- omics_save_base_plot(
    function() {
      graphics::plot(cal, xlab = sprintf("Predicted survival at t = %.3g", horizon),
                     ylab = "Observed survival", subtitles = FALSE)
      graphics::abline(0, 1, lty = 3, col = "#666666")
    },
    "survival", name, width = 6, height = 6
  )

  list(
    method = "calibration",
    n = io$n,
    events = io$events,
    horizon = horizon,
    per_group = units,
    bootstrap = as.integer(params$bootstrap %||% 200),
    figure = figure,
    notes = "Bootstrap-corrected calibration is internal validation only; external validation remains necessary."
  )
}

run_dca <- function(params) {
  omics_require(c("survival"))
  io <- prepare(params)
  df <- io$df
  horizon <- as.numeric((params$times %||% stats::median(df[[io$time]]))[1])
  formula <- stats::as.formula(paste("survival::Surv(", io$time, ",", io$event, ") ~", paste(io$predictors, collapse = " + ")))
  fit <- survival::coxph(formula, data = df)
  surv <- summary(survival::survfit(fit, newdata = df), times = horizon)
  predicted_risk <- 1 - as.numeric(surv$surv)

  thresholds <- as.numeric(params$thresholds %||% seq(0.01, 0.60, by = 0.01))
  km_overall <- summary(survival::survfit(stats::as.formula(paste("survival::Surv(", io$time, ",", io$event, ") ~ 1")), data = df), times = horizon)
  event_rate <- 1 - as.numeric(km_overall$surv)

  net_benefit <- vapply(thresholds, function(pt) {
    flagged <- predicted_risk >= pt
    if (!any(flagged)) return(0)
    sub <- df[flagged, , drop = FALSE]
    km <- summary(survival::survfit(stats::as.formula(paste("survival::Surv(", io$time, ",", io$event, ") ~ 1")), data = sub), times = horizon)
    rate <- 1 - as.numeric(km$surv)
    if (!length(rate) || is.na(rate)) return(NA_real_)
    tp <- rate * mean(flagged)
    fp <- (1 - rate) * mean(flagged)
    tp - fp * (pt / (1 - pt))
  }, numeric(1))

  treat_all <- vapply(thresholds, function(pt) event_rate - (1 - event_rate) * (pt / (1 - pt)), numeric(1))
  curve <- data.frame(threshold = thresholds, model = net_benefit, treat_all = treat_all, treat_none = 0)

  omics_require("ggplot2")
  long <- data.frame(
    threshold = rep(curve$threshold, 3),
    net_benefit = c(curve$model, curve$treat_all, curve$treat_none),
    strategy = rep(c("Model", "Treat all", "Treat none"), each = nrow(curve)),
    stringsAsFactors = FALSE
  )
  ymin <- min(-0.02, stats::quantile(long$net_benefit, 0.02, na.rm = TRUE))
  plot <- ggplot2::ggplot(long, ggplot2::aes(x = threshold, y = net_benefit, colour = strategy)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = c(Model = "#b5483a", `Treat all` = "#2f5f8f", `Treat none` = "#888888")) +
    ggplot2::coord_cartesian(ylim = c(ymin, max(long$net_benefit, na.rm = TRUE) * 1.05)) +
    ggplot2::labs(x = "Threshold probability", y = "Net benefit", colour = NULL,
                  title = sprintf("Decision curve at t = %.3g", horizon)) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  name <- as.character(params$output_name %||% "dca")[1]
  figure <- omics_save_plot(plot, "survival", name, width = 6.5, height = 5)
  table_path <- omics_save_table(curve, "survival", paste0(name, "_net_benefit"))

  list(
    method = "dca",
    n = io$n,
    events = io$events,
    horizon = horizon,
    event_rate_at_horizon = event_rate,
    net_benefit_table = table_path,
    figure = figure,
    notes = "The model is clinically useful only over the threshold range where its net benefit exceeds both treat-all and treat-none."
  )
}

handler <- function(params) {
  method <- tolower(as.character(params$method %||% "lasso_cox")[1])
  switch(
    method,
    lasso_cox = run_lasso_cox(params),
    timeroc = run_timeroc(params),
    nomogram = run_nomogram(params),
    calibration = run_calibration(params),
    dca = run_dca(params),
    stop(sprintf("Unknown method '%s'. Use lasso_cox, timeroc, nomogram, calibration, or dca.", method), call. = FALSE)
  )
}

omics_main(handler)
