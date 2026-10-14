#' Refit nuisance models and estimators in each bootstrap resample.
#' Resample within the original folds, so duplicate participants stay in one fold.
#' A supplied fixed_propensity is retained in every resample.
bootstrap_tree <- function(obs_data, trt_predictors, outcome_predictors, B = 100, V = 5,
                           K = 4, n_timepoints = length(outcome_columns), seed = 1,
                           estimators = c("aipw", "tmle", "aipw_cf", "tmle_cf"),
                           tol = 1e-8, maxit = 1000, outcome_columns,
                           outcome_model = "tree", regime = NULL, folds = NULL,
                           tmle_scale = "original", tree_control = list(), fixed_propensity = NULL) {
  n <- nrow(obs_data)
  out <- setNames(lapply(estimators, function(est) matrix(NA_real_, B, n_timepoints)), estimators)
  use_full <- any(estimators %in% c("aipw", "tmle"))
  use_cf <- any(estimators %in% c("aipw_cf", "tmle_cf", "tmle_cf_pooled"))
  original_folds <- if (is.null(folds)) kfold_indices_cs(n, V, seed) else folds
  fold_ids <- sort(unique(original_folds))
  progress <- bootstrap_progress_logger(outcome_model, regime, n, K, V, seed, B, estimators)
  progress$write("started", 0L, 0L, out)
  completed <- b <- 0L
  finished <- FALSE
  on.exit(if (!finished) progress$write("interrupted", completed, b, out), add = TRUE)
  for (b in seq_len(B)) {
    report <- b == 1L || b %% progress$every == 0L || b == B
    if (report) progress$write("running", completed, b, out)
    set.seed(seed + b)
    sampled_rows <- unlist(lapply(fold_ids, function(v)
      sample(which(original_folds == v), sum(original_folds == v), replace = TRUE)))
    folds_b <- original_folds[sampled_rows]
    dat_b <- obs_data[sampled_rows, , drop = FALSE]
    rownames(dat_b) <- NULL
    if (nuisance_support_complete(dat_b, trt_predictors, outcome_predictors, K, V, seed + b,
                                  folds = folds_b, cf = use_cf, outcome_columns = outcome_columns)) {
      nuis <- nuis_CF <- NULL
      if (use_full) nuis <- tryCatch(
        fit_nuisance_ind(dat_b, trt_predictors, outcome_predictors, K, n_timepoints, outcome_model,
                         varc = FALSE, outcome_columns = outcome_columns, tree_control = tree_control,
                         fixed_propensity = fixed_propensity),
        error = function(e) NULL)
      if (use_cf) nuis_CF <- tryCatch(
        fit_nuisance_ind_CF(dat_b, V, trt_predictors, outcome_predictors, K, n_timepoints, seed + b,
                            outcome_model, folds = folds_b, varc = FALSE,
                            outcome_columns = outcome_columns, tree_control = tree_control,
                            fixed_propensity = fixed_propensity),
        error = function(e) NULL)
      for (est in estimators) {
        if ((est %in% c("aipw", "tmle") && is.null(nuis)) ||
            (est %in% c("aipw_cf", "tmle_cf", "tmle_cf_pooled") && is.null(nuis_CF))) next
        fit <- tryCatch(switch(est,
                               aipw = AIPW(dat_b, nuis, K, n_timepoints, outcome_columns),
                               aipw_cf = AIPW_CF(dat_b, nuis_CF, K, n_timepoints, outcome_columns),
                               tmle = TMLE(dat_b, nuis, K, n_timepoints, tol, maxit, outcome_columns, tmle_scale = tmle_scale),
                               tmle_cf = TMLE_CF(dat_b, nuis_CF, K, n_timepoints, tol, maxit, outcome_columns, tmle_scale = tmle_scale),
                               tmle_cf_pooled = TMLE_CF_pooled(dat_b, nuis_CF, K, n_timepoints, tol, maxit, outcome_columns,
                                                               tmle_scale = tmle_scale)
        ), error = function(e) NULL)
        if (!is.null(fit) && all(is.finite(fit$psi)) &&
            (is.null(fit$target_convergence) || fit$target_convergence == 0L)) out[[est]][b, ] <- fit$psi
      }
    }
    completed <- b
    if (report) progress$write("completed_resample", completed, b, out)
  }
  progress$write("finished", completed, B, out)
  finished <- TRUE
  list(progress_log = progress$path, estimates = out,
       valid_B = vapply(out, function(x) sum(complete.cases(x)), integer(1)))
}



#' Record bootstrap progress and successful resample counts in a tab-separated log.
bootstrap_progress_logger <- function(outcome_model, regime, n, K, V, seed, B,
                                      estimators) {
  every <- suppressWarnings(as.integer(Sys.getenv("BOOTSTRAP_PROGRESS_EVERY", "1")))
  if (length(every) != 1L || is.na(every) || every < 1L)
    stop("BOOTSTRAP_PROGRESS_EVERY must be a positive integer")
  directory <- Sys.getenv("BOOTSTRAP_LOG_DIR",
                          file.path(Sys.getenv("RESULT_DIR", "results_DOOR18"), "bootstrap_progress"))
  regime_label <- if (is.null(regime)) "custom" else as.character(regime)
  safe <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)
  prefix <- sprintf("%s_R%s_n%d_K%d_seed%d_pid%d_", safe(outcome_model),
                    safe(regime_label), n, K, seed, Sys.getpid())
  path <- NA_character_
  enabled <- TRUE
  start <- proc.time()[["elapsed"]]
  tryCatch({
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    path <- tempfile(prefix, tmpdir = normalizePath(directory, mustWork = TRUE),
                     fileext = ".tsv")
  }, error = function(e) {
    enabled <<- FALSE
    warning("Bootstrap progress logging unavailable: ", conditionMessage(e))
  })
  write <- function(status, completed, current, estimates) {
    if (!enabled) return(invisible(NULL))
    counts <- vapply(estimates, function(x) sum(complete.cases(x)), integer(1))
    row <- data.frame(timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3%z"),
                      status = status, learner = outcome_model, regime = regime_label,
                      n = n, K = K, V = V, seed = seed, probability_floor = probability_clip(),
                      pid = Sys.getpid(), completed = completed, B = B, current = current,
                      elapsed_seconds = round(proc.time()[["elapsed"]] - start, 3),
                      stringsAsFactors = FALSE)
    for (est in estimators) row[[paste0("valid_", est)]] <- counts[[est]]
    tryCatch(write.table(row, path, sep = "\t", row.names = FALSE,
                         col.names = !file.exists(path), append = file.exists(path),
                         quote = FALSE), error = function(e) {
                           enabled <<- FALSE
                           warning("Bootstrap progress logging stopped: ", conditionMessage(e))
                         })
    invisible(NULL)
  }
  list(write = write, every = every, path = path)
}
