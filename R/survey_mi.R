#' Simple survey MI analysis: one imputation per bootstrap replicate
#'
#' Draws multinomial bootstrap replicate weights `Wt_1, ..., Wt_B`. For each
#' replicate `b`, imputes missing values once using `Wt_b`, fits the substantive
#' weighted model, and stores one coefficient vector. Returns all `B` estimates
#' plus an empirical summary (mean and percentile 95% interval).
#'
#' @param design A [survey::svydesign()] object.
#' @param analysis Two-sided formula for the substantive model (e.g. `y ~ x1 + x2`).
#' @param formulas,variables,methods Imputation specification; same rules as
#'   [survey_mi] (provide `formulas` or `variables` + `methods`).
#' @param analysis_family `"gaussian"`, `"binomial"`, or `"poisson"` for the
#'   substantive model fit.
#' @param B Number of multinomial bootstrap replicate columns from
#'   [multinom_weights()]. Defaults to `100`.
#' @param n_rounds Number of chained-equation rounds for each imputation.
#' @param seed Optional seed for reproducible bootstrap weights and imputation draws.
#' @param fitters Optional fitter overrides; same as [survey_mi].
#'
#' @return A list with:
#'   \describe{
#'     \item{coefs}{matrix (`B x p`) of coefficient estimates by bootstrap replicate.}
#'     \item{coef_names}{coefficient names.}
#'     \item{summary}{data frame with `term`, `estimate`, `conf.low`, `conf.high`.}
#'     \item{B}{number of bootstrap replicates used.}
#'   }
#' @export
survey_mi <- function(design,
                      analysis,
                      formulas = NULL,
                      variables = NULL,
                      methods = NULL,
                      analysis_family = c("gaussian", "binomial", "poisson"),
                      B = 100L,
                      n_rounds = 10L,
                      seed = NULL,
                      fitters = NULL) {
  stopifnot(inherits(design, "survey.design"))
  stopifnot(inherits(analysis, "formula"))
  analysis_family <- match.arg(analysis_family)
  B <- as.integer(B)
  stopifnot(B >= 1L, n_rounds >= 1L)

  if (!is.null(formulas) && !is.null(variables)) {
    stop("Specify either `formulas` or (`variables` and `methods`), not both.")
  }

  df <- .design_variables_df(design)

  if (!is.null(seed)) {
    set.seed(seed)
  }
  wt_df <- multinom_weights(design, B = B)
  wcols <- grep("^Wt_", names(wt_df), value = TRUE)
  if (!length(wcols)) {
    stop("multinom_weights() must return columns `Wt_1`, ..., `Wt_B`.")
  }
  mat_w <- as.matrix(wt_df[, wcols, drop = FALSE])
  stopifnot(ncol(mat_w) == B)

  if (!is.null(formulas)) {
    parsed <- .formulas_resolve(formulas, df)
    df_fit <- df
  } else if (!is.null(variables)) {
    if (is.null(methods)) {
      stop("With `variables`, supply named `methods` for incomplete outcomes.")
    }
    vars <- unique(as.character(variables))
    unk <- setdiff(vars, names(df))
    if (length(unk)) {
      stop("Unknown `variables`: ", paste(unk, collapse = ", "))
    }
    df_fit <- df[, vars, drop = FALSE]
    parsed <- .variables_methods_resolve(variables, methods, df_fit)
  } else {
    stop("Supply `formulas` or (`variables` and `methods`).")
  }

  av <- all.vars(stats::terms(analysis))
  if (!all(av %in% names(df_fit))) {
    stop(
      "All variables in `analysis` must appear in the imputation working data ",
      "(use `variables` or `formulas` on the appropriate columns)."
    )
  }

  fitters_map <- .resolve_fitters(fitters)
  tgt_set <- parsed$tgt_set
  if (!length(tgt_set)) {
    stop("No incomplete variables to impute; use `survey::svyglm` directly.")
  }

  preds_pool <- names(df_fit)[!names(df_fit) %in% tgt_set]
  frm_by_target <- parsed$frm_by_target[tgt_set]
  method_map <- parsed$method_map
  mc <- vapply(df_fit[tgt_set], function(x) sum(is.na(x)), integer(1))
  vars_ord <- tgt_set[order(mc[tgt_set] / nrow(df_fit), tgt_set)]

  coefs <- NULL
  cnm <- NULL

  for (b in seq_len(B)) {
    if (!is.null(seed)) {
      set.seed(.survey_mi_subseed(seed, b, 1L))
    }
    wb <- mat_w[, b, drop = TRUE]
    imp_b <- .survey_mi_single_draw(
      df_fit,
      wb,
      frm_by_target,
      method_map,
      tgt_set,
      n_rounds,
      fitters_map,
      vars_ord,
      preds_pool
    )
    fm <- .fit_weighted_analysis(analysis, imp_b, wb, analysis_family)

    if (is.null(coefs)) {
      cnm <- names(fm$coef)
      coefs <- matrix(
        NA_real_,
        nrow = B,
        ncol = length(cnm),
        dimnames = list(paste0("b", seq_len(B)), cnm)
      )
    }
    coefs[b, ] <- unname(as.numeric(fm$coef))
  }

  estimate <- colMeans(coefs, na.rm = TRUE)
  conf_low <- apply(coefs, 2L, stats::quantile, probs = 0.025, na.rm = TRUE, names = FALSE)
  conf_high <- apply(coefs, 2L, stats::quantile, probs = 0.975, na.rm = TRUE, names = FALSE)
  summary_df <- data.frame(
    term = colnames(coefs),
    estimate = unname(estimate),
    conf.low = unname(conf_low),
    conf.high = unname(conf_high),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      coefs = coefs,
      coef_names = cnm,
      summary = summary_df,
      B = B,
      analysis = analysis,
      analysis_family = analysis_family
    ),
    class = c("survey_mi", "list")
  )
}


.survey_mi_subseed <- function(seed, b, a) {
  s <- as.numeric(seed) + b * 10007L + a * 30011L
  as.integer(abs(s) %% 2147483629L) + 1L
}


.fit_weighted_analysis <- function(analysis, data, weights, family) {
  need_vars <- all.vars(stats::terms(analysis))

  ok <- stats::complete.cases(data[, need_vars, drop = FALSE])
  xf <- data[ok, , drop = FALSE]
  ww <- as.numeric(weights)[ok]

  mf <- stats::model.frame(analysis, data = xf, na.action = stats::na.fail)
  tt <- stats::terms(analysis)
  tt <- stats::terms(tt, data = mf)

  x <- stats::model.matrix(tt, mf)
  y <- stats::model.response(mf, type = "any")

  if (nrow(x) != length(ww)) {
    stop("internal mismatch: model matrix nrow vs weights.")
  }

  nm <- colnames(x)

  if (family == "binomial") {
    if (is.factor(y)) {
      y <- as.integer(y) - 1L
    } else if (is.logical(y)) {
      y <- as.integer(y)
    } else {
      y <- as.numeric(y)
    }
  } else if (family == "poisson") {
    y <- as.numeric(y)
  } else {
    y <- as.numeric(y)
  }

  cpp <- function(name) get0(name, inherits = TRUE, mode = "function")

  if (family == "gaussian") {
    wg <- cpp("weighted_gaussian_cpp")
    if (!is.null(wg)) {
      rr <- wg(x, y, ww)
      b <- .coef_nm(rr$coef, nm)
      V <- rr$vcov
      dimnames(V) <- list(nm, nm)
      return(list(coef = b, vcov = V))
    }
    ff <- stats::lm.wfit(x = x, y = y, w = ww)
    ne <- sum(ww > 0)
    p <- ncol(x)
    rss <- sum(ww * (y - drop(x %*% ff$coefficients))^2)
    sg2 <- if (ne > p) rss / (ne - p) else max(rss, 1e-8)
    xtwx <- crossprod(x, x * ww)
    vc <- tryCatch(sg2 * solve(xtwx), error = function(e) sg2 * MASS::ginv(xtwx))
    dimnames(vc) <- list(nm, nm)
    list(coef = .coef_nm(ff$coefficients, nm), vcov = vc)
  } else if (family == "binomial") {
    wg <- cpp("weighted_binomial_cpp")
    if (!is.null(wg)) {
      rr <- wg(x, y, ww)
      b <- .coef_nm(rr$coef, nm)
      V <- rr$vcov
      dimnames(V) <- list(nm, nm)
      return(list(coef = b, vcov = V))
    }
    z <- stats::glm.fit(x = x, y = y, weights = ww, family = stats::binomial())
    mu <- pmin(pmax(z$fitted.values, 1e-6), 1 - 1e-6)
    xtwx <- crossprod(x, x * (ww * mu * (1 - mu)))
    vc <- tryCatch(solve(xtwx), error = function(e) MASS::ginv(xtwx))
    dimnames(vc) <- list(nm, nm)
    list(coef = .coef_nm(z$coefficients, nm), vcov = vc)
  } else if (family == "poisson") {
    wg <- cpp("weighted_poisson_cpp")
    if (!is.null(wg)) {
      rr <- wg(x, y, ww)
      b <- .coef_nm(rr$coef, nm)
      V <- rr$vcov
      dimnames(V) <- list(nm, nm)
      return(list(coef = b, vcov = V))
    }
    z <- stats::glm.fit(x = x, y = y, weights = ww, family = stats::poisson())
    mu <- pmax(z$fitted.values, 1e-6)
    xtwx <- crossprod(x, x * (ww * mu))
    vc <- tryCatch(solve(xtwx), error = function(e) MASS::ginv(xtwx))
    dimnames(vc) <- list(nm, nm)
    list(coef = .coef_nm(z$coefficients, nm), vcov = vc)
  }

  stop("unknown `family`")
}


.coef_nm <- function(coef, nm) {
  coef <- unname(as.numeric(coef))
  stats::setNames(coef, nm)
}

#' Sequential weighted imputation for complex surveys
#'
#' Sequential chained equations (least to most missing) with multinomial PSU
#' bootstrap replicate weights from [multinom_weights()] whenever `design`
#' inherits `"survey.design"`.
#'
#' @param design Either a survey design (`inherits(, "survey.design")`) —
#'   typically from [survey::svydesign()] — **or** a `data.frame` (legacy usage).
#'
#' @param formulas Explicit conditional specifications (optional): a named list
#'   grouping models by family (`gaussian`, `binomial`, `poisson`), each entry one
#'   formula or list of formulas, LHS a single bare name. Omit this if you supply
#'   **`variables` + `methods`** instead (PROC MI–style VAR + method).
#'
#'   Do not pass **`formulas`** together with **`variables`** — choose one style.
#'
#' @param variables **PROC MI–style VAR list.** Character vector of analysis
#'   variables that enter the imputation (complete and incomplete). For each
#'   incomplete `y` among these names, a model is built automatically as
#'   `y ~ <all other names in variables>` (additive main effects only).
#'   **Only these columns** are passed to the imputation engine; all other
#'   columns in `design`/`data` are ignored for fitting. **Returned imputations
#'   are data frames with exactly these columns** (same row order as the survey
#'   design). Add any auxiliary predictors (e.g. stratification controls) to
#'   `variables` if you want them in the model.
#'
#' @param ... Currently unused (reserved).
#'
#' @param B Number of multinomial bootstrap replicate weight columns from
#'   [multinom_weights()] (**survey designs only**). Defaults to **`100`** when
#'   **`NULL`** (with **`n_imputations`** defaulting to **`5`**).
#'
#' @param weights Analysis weights (`data.frame` path only).
#' @param methods **With `variables` (PROC MI–style):** named list `list(gaussian =
#'   ..., binomial = ..., poisson = ...)` with entries as character vectors.
#'   List **only incomplete variables** that need imputation, **each exactly once**
#'   across the three families. Omit a family (or set it to `character(0)` or
#'   `NULL`) when unused.
#'
#'   **`data.frame`, `formulas = NULL`, `variables = NULL`:** legacy behaviour
#'   over all missing columns in the full frame (see `.resolve_methods()`).
#' @param n_rounds Rounds `c` in Raghunathan et al..
#' @param n_imputations Number of stochastic MI draws **per bootstrap replicate**
#'   (**survey:** denoted **`A`**; total completed datasets = **`B * A`**). Defaults
#'   to **`5`** when **`NULL`** on survey inputs. For a plain **`data.frame`**,
#'   this is the total number of imputations (same default **`5`** when **`NULL`**).
#' @param seed Optional; passed to \code{\link{set.seed}} prior to RNG use.
#' @param fitters Optional overrides for `gaussian`, `binomial`, `poisson` fitters.
#' @param parallel If **`TRUE`** (survey designs only, with **`block_info`**),
#'   run bootstrap blocks on parallel workers via **[parallel::mclapply]** (Unix/macOS
#'   fork only). Windows falls back to sequential with a message.
#' @param n_cores Maximum worker processes when **`parallel`** is **`TRUE`**;
#'   default **`detectCores()-1`** (minimum **`1`**).
#'
#' @return List of completed `data.frame`s: length **`n_imputations`** (plain
#'   **`data.frame`**), or **`B * A`** for survey designs ( **`A`** =
#'   **`n_imputations`** per **`Wt_b`**). Attributes **`B`**, **`A`**, and
#'   **`bootstrap_rep`** map each list element to its replicate index **`b`**.
#' @keywords internal
survey_mi_engine <- function(design,
                      formulas = NULL,
                      ...,
                      weights = NULL,
                      variables = NULL,
                      methods = NULL,
                      B = NULL,
                      n_rounds = 10,
                      n_imputations = NULL,
                      seed = NULL,
                      fitters = NULL,
                      parallel = FALSE,
                      n_cores = NULL) {
  stopifnot(n_rounds >= 1L)
  if (!inherits(design, "survey.design")) {
    if (is.null(n_imputations)) {
      n_imputations <- 5L
    }
    stopifnot(n_imputations >= 1L)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!is.null(formulas) && !is.null(variables)) {
    stop("Specify either `formulas` or (`variables` and `methods`), not both.")
  }

  if (inherits(design, "survey.design")) {
    if (!is.null(weights)) {
      stop("`weights` is not used for survey inputs; multinom replicate weights fill that role.")
    }
    if (is.null(B) && is.null(n_imputations)) {
      B <- 100L
      A <- 5L
    } else if (is.null(B)) {
      B <- 100L
      A <- as.integer(n_imputations)
    } else if (is.null(n_imputations)) {
      B <- as.integer(B)
      A <- 5L
    } else {
      B <- as.integer(B)
      A <- as.integer(n_imputations)
    }
    stopifnot(B >= 1L, A >= 1L)

    df <- .design_variables_df(design)

    if (!is.null(seed)) {
      set.seed(seed)
    }
    wt_df <- multinom_weights(design, B = B)
    wcols <- grep("^Wt_", names(wt_df), value = TRUE)
    if (!length(wcols)) {
      stop("multinom_weights() must return columns `Wt_1`, …, `Wt_B`.")
    }
    mat_w <- as.matrix(wt_df[, wcols, drop = FALSE])
    if (ncol(mat_w) != B) {
      stop("Internal mismatch: multinom_weights columns vs `B`.")
    }
    n_imputations <- B * A
    which_wt <- rep(seq_len(B), each = A)

    if (!is.null(formulas)) {
      parsed <- .formulas_resolve(formulas, df)
      df_fit <- df
    } else if (!is.null(variables)) {
      if (is.null(methods)) {
        stop("With `variables`, supply named `methods` listing incomplete variables by family.")
      }

      vars <- unique(as.character(variables))
      unk <- setdiff(vars, names(df))
      if (length(unk)) {
        stop("Unknown `variables`: ", paste(unk, collapse = ", "))
      }

      df_fit <- df[, vars, drop = FALSE]

      parsed <- .variables_methods_resolve(variables, methods, df_fit)
    } else {
      stop(
        "For survey designs, supply `formulas` or (`variables` and `methods`)."
      )
    }

    imps <- survey_mi_fit(
      data = df_fit,
      weights_rep = mat_w,
      weight_col_idx = which_wt,
      frm_by_target = parsed$frm_by_target,
      method_map = parsed$method_map,
      tgt_set = parsed$tgt_set,
      n_rounds = n_rounds,
      n_imputations = n_imputations,
      fitters = fitters,
      parallel = parallel,
      n_cores = n_cores,
      block_info = list(B = B, A = A)
    )
    attr(imps, "B") <- B
    attr(imps, "A") <- A
    attr(imps, "bootstrap_rep") <- rep(seq_len(B), each = A)
    imps
  } else if (is.data.frame(design)) {
    if (!is.null(formulas)) {
      parsed <- .formulas_resolve(formulas, design)

      nr <- nrow(design)
      w <- if (is.null(weights)) rep(1, nr) else weights
      stopifnot(length(w) == nr)
      if (any(!is.finite(w)) || any(w < 0)) {
        stop("`weights` must be finite and non-negative.")
      }
      survey_mi_fit(
        data = design,
        weights_rep = matrix(as.numeric(w), nrow = nr, ncol = 1L),
        weight_col_idx = rep(1L, n_imputations),
        frm_by_target = parsed$frm_by_target,
        method_map = parsed$method_map,
        tgt_set = parsed$tgt_set,
        n_rounds = n_rounds,
        n_imputations = n_imputations,
        fitters = fitters
      )
    } else if (!is.null(variables)) {
      if (is.null(methods)) {
        stop("With `variables`, supply named `methods` listing incomplete variables by family.")
      }

      vars <- unique(as.character(variables))
      unk <- setdiff(vars, names(design))
      if (length(unk)) {
        stop("Unknown `variables`: ", paste(unk, collapse = ", "))
      }

      nr <- nrow(design)
      w <- if (is.null(weights)) rep(1, nr) else weights
      stopifnot(length(w) == nr)
      if (any(!is.finite(w)) || any(w < 0)) {
        stop("`weights` must be finite and non-negative.")
      }

      df_fit <- design[, vars, drop = FALSE]

      parsed <- .variables_methods_resolve(variables, methods, df_fit)
      survey_mi_fit(
        data = df_fit,
        weights_rep = matrix(as.numeric(w), nrow = nr, ncol = 1L),
        weight_col_idx = rep(1L, n_imputations),
        frm_by_target = parsed$frm_by_target,
        method_map = parsed$method_map,
        tgt_set = parsed$tgt_set,
        n_rounds = n_rounds,
        n_imputations = n_imputations,
        fitters = fitters
      )
    } else {
      survey_mi_legacy_df(
        data = design,
        weights = weights,
        methods = methods,
        n_rounds = n_rounds,
        n_imputations = n_imputations,
        fitters = fitters
      )
    }
  } else {
    stop("`design` must inherit \"survey.design\" or be a data.frame.")
  }
}


survey_mi_legacy_df <- function(data,
                               weights,
                               methods,
                               n_rounds,
                               n_imputations,
                               fitters) {
  stopifnot(is.data.frame(data))
  n <- nrow(data)
  w <- if (is.null(weights)) rep(1, n) else weights
  stopifnot(length(w) == n)
  if (any(!is.finite(w)) || any(w < 0)) {
    stop("`weights` must be finite and non-negative.")
  }

  mc <- vapply(data, function(x) sum(is.na(x)), integer(1))
  vmiss <- names(mc[mc > 0])
  if (!length(vmiss)) {
    return(replicate(n_imputations, data, simplify = FALSE))
  }

  vmiss <- vmiss[order(mc[vmiss] / n, vmiss)]

  preds <- names(mc[mc == 0L])

  frms <- vector("list", length(vmiss))
  names(frms) <- vmiss
  for (tg in vmiss) {
    rhs <- preds[preds %in% names(data)]
    frms[[tg]] <- if (length(rhs)) stats::reformulate(rhs, response = tg) else
      stats::as.formula(paste(tg, "~ 1"))
  }

  mth <- .resolve_methods(data, vmiss, methods)

  survey_mi_fit(
    data = data,
    weights_rep = matrix(as.numeric(w), nrow = n, ncol = 1L),
    weight_col_idx = rep(1L, n_imputations),
    frm_by_target = frms,
    method_map = mth,
    tgt_set = vmiss,
    n_rounds = n_rounds,
    n_imputations = n_imputations,
    fitters = fitters
  )
}


.design_variables_df <- function(design) {
  v <- design$variables
  if (!is.data.frame(v)) {
    return(as.data.frame(v))
  }
  v
}


.lhs_symbol <- function(frm) {
  if (!inherits(frm, "formula") || length(frm) != 3L) {
    stop("Each conditional model must be a two-sided formula `y ~ ....`.")
  }
  lhs <- frm[[2L]]
  if (!is.symbol(lhs)) {
    stop(
      "Left-hand sides must be a single variable symbol; ",
      deparse(lhs, width.cutoff = 500), " not supported yet."
    )
  }
  as.character(lhs)
}


.restrict_formula_rhs <- function(frm_template, tgt, avail, anchor) {
  tt <- stats::terms(frm_template, data = anchor)
  lbl <- attr(tt, "term.labels")
  rsp <- tgt

  if (!length(lbl)) {
    f <- stats::as.formula(paste(rsp, "~ 1"))
    environment(f) <- environment(stats::formula(frm_template))
    return(f)
  }

  ok <- vapply(lbl, function(tm) {
    av <- tryCatch(
      all.vars(as.formula(paste("~", tm))),
      error = function(...) NA_character_,
      silent = TRUE
    )
    if (!length(av) || anyNA(av)) {
      FALSE
    } else {
      all(av %in% avail)
    }
  }, logical(1))

  lbl2 <- lbl[ok]
  f <- if (length(lbl2)) stats::reformulate(lbl2, response = rsp) else
    stats::as.formula(paste(rsp, "~ 1"))

  environment(f) <- environment(stats::formula(frm_template))
  f
}


.formulas_resolve <- function(formulas, data) {
  fams <- c("gaussian", "binomial", "poisson")

  stopifnot(is.list(formulas))
  if (!length(formulas)) {
    stop("`formulas` is empty.")
  }
  bn <- names(formulas)
  if (anyDuplicated(bn) || !all(bn %in% fams)) {
    stop("`formulas` must have unique names in {gaussian, binomial, poisson}.")
  }

  frm_out <- list()
  meth_by_tgt <- character(0)

  for (fam in names(formulas)) {
    elm <- formulas[[fam]]
    if (inherits(elm, "formula")) {
      elm <- list(elm)
    } else if (!is.list(elm)) {
      stop("`formulas$", fam, "` must be one formula or a list of formulas.")
    }
    for (f_i in elm) {
      tg <- .lhs_symbol(f_i)
      if (!tg %in% names(data)) {
        stop("`", tg, "` (LHS formula) does not occur in design data.")
      }
      if (tg %in% names(frm_out)) {
        stop("Duplicate outcome `", tg, "` in `formulas`.")
      }

      frm_out[[tg]] <- stats::formula(stats::terms(f_i, data = data))

      meth_by_tgt[tg] <- fam
    }
  }

  mc_all <- vapply(data, function(x) sum(is.na(x)), integer(1))
  tgt_nm <- names(frm_out)
  need_na <- names(mc_all[mc_all > 0L])

  miss_spec <- setdiff(need_na, tgt_nm)
  if (length(miss_spec)) {
    stop(
      "Columns with NA must each appear as an LHS exactly once:\nmissing: ",
      paste(miss_spec, collapse = ", ")
    )
  }

  surplus <- tgt_nm[!(tgt_nm %in% need_na)]
  if (length(surplus)) {
    warning(
      "LHS has no NA in data (models still skipped if nowhere needed):\n",
      paste(surplus, collapse = ", "),
      call. = FALSE
    )
  }

  tg_order <- tgt_nm[
    order(mc_all[tgt_nm] / nrow(data), tgt_nm)
  ]

  frm_ord <- frm_out[tg_order]
  mlist <- stats::setNames(vector("list", length(tg_order)), tg_order)
  for (nm in tg_order) {
    mlist[[nm]] <- meth_by_tgt[[nm]]
  }

  list(
    frm_by_target = frm_ord,
    method_map = mlist,
    tgt_set = tg_order
  )
}


.flatten_family_methods <- function(methods) {
  fams <- c("gaussian", "binomial", "poisson")
  if (!is.list(methods) || !length(methods)) {
    stop("`methods` must be a non-empty named list.")
  }
  nm <- names(methods)
  if (anyDuplicated(nm) || !all(nm %in% fams)) {
    stop("`methods` names must be unique and from {gaussian, binomial, poisson}.")
  }
  mp <- character(0)
  for (fam in nm) {
    v <- methods[[fam]]
    if (is.null(v)) next
    v <- unique(as.character(unlist(v, use.names = FALSE)))
    v <- v[!is.na(v) & nzchar(v)]
    dup <- intersect(v, names(mp))
    if (length(dup)) {
      stop(
        "Variable(s) appear in more than one `methods` family: ",
        paste(dup, collapse = ", ")
      )
    }
    mp[v] <- fam
  }
  mp
}


.variables_methods_resolve <- function(variables, methods, data) {
  vars <- unique(as.character(variables))
  if (!length(vars)) {
    stop("`variables` is empty.")
  }
  unk <- setdiff(vars, names(data))
  if (length(unk)) {
    stop("Unknown `variables`: ", paste(unk, collapse = ", "))
  }

  mmap <- .flatten_family_methods(methods)
  extra <- setdiff(names(mmap), vars)
  if (length(extra)) {
    stop(
      "These names are in `methods` but not in `variables`: ",
      paste(extra, collapse = ", ")
    )
  }
  imputed <- vars[vapply(vars, function(v) anyNA(data[[v]]), logical(1))]

  if (!length(imputed)) {
    warning("No variable in `variables` has missing values.", call. = FALSE)
    return(list(
      frm_by_target = list(),
      method_map = list(),
      tgt_set = character(0)
    ))
  }

  mc_all <- vapply(data, function(x) sum(is.na(x)), integer(1))
  missed <- setdiff(imputed, names(mmap))
  if (length(missed)) {
    stop(
      "Assign a method (gaussian/binomial/poisson) for incomplete variable(s): ",
      paste(missed, collapse = ", ")
    )
  }

  frm_out <- list()
  for (tg in imputed) {
    rhs <- setdiff(vars, tg)
    frm_out[[tg]] <- if (length(rhs)) {
      stats::reformulate(rhs, response = tg)
    } else {
      stats::as.formula(paste(tg, "~ 1"))
    }
  }

  tg_order <- imputed[
    order(mc_all[imputed] / nrow(data), imputed)
  ]

  frm_ord <- frm_out[tg_order]
  mlist <- stats::setNames(vector("list", length(tg_order)), tg_order)
  for (nm in tg_order) {
    mlist[[nm]] <- mmap[[nm]]
  }

  list(
    frm_by_target = frm_ord,
    method_map = mlist,
    tgt_set = tg_order
  )
}


.survey_mi_single_draw <- function(data,
                                   wcol,
                                   frm_by_target,
                                   method_map,
                                   tgt_set,
                                   n_rounds,
                                   fitters_map,
                                   vars_ord,
                                   preds_pool) {
  y_imp <- data
  for (tg in tgt_set) {
    y_imp[[tg]] <- .initialize_missing(y_imp[[tg]])
  }

  for (rd in seq_len(n_rounds)) {
    for (j in seq_along(vars_ord)) {
      tg <- vars_ord[j]

      if (!any(is.na(data[[tg]]))) {
        next
      }

      if (rd == 1L) {
        avail <- unique(c(preds_pool, vars_ord[seq_len(j - 1L)]))
      } else {
        avail <- unique(c(setdiff(vars_ord, tg), preds_pool))
      }

      frm_r <- .restrict_formula_rhs(frm_by_target[[tg]], tg, avail = avail, anchor = data)

      tt <- stats::terms(frm_r, data = y_imp)
      mf_fit <- stats::model.frame(tt, data = y_imp, na.action = stats::na.pass)

      tgt_raw_vec <- stats::model.response(mf_fit, type = "any")
      if (is.null(tgt_raw_vec)) {
        next
      }
      cc_fit <- stats::complete.cases(mf_fit) & !is.na(tgt_raw_vec)

      if (!any(cc_fit) || sum(cc_fit) < 3L) {
        next
      }

      mf_cc <- mf_fit[cc_fit, , drop = FALSE]
      x_fit <- stats::model.matrix(tt, data = mf_cc)
      target_raw <- data[[tg]]
      y_fit <- stats::model.response(mf_cc, type = "any")
      meth <- method_map[[tg]]
      if (meth == "binomial") {
        if (is.factor(y_fit)) {
          y_fit <- as.integer(y_fit) - 1L
        } else if (is.logical(y_fit)) {
          y_fit <- as.integer(y_fit)
        } else {
          y_fit <- as.numeric(y_fit)
        }
      } else {
        y_fit <- as.numeric(y_fit)
      }
      w_fit <- as.numeric(wcol[cc_fit])

      fh <- fitters_map[[meth]]
      fh_fit <- fh(x_fit, y_fit, w_fit)

      is_mis <- is.na(target_raw)
      mf_miss <- stats::model.frame(
        tt,
        data = y_imp[is_mis, , drop = FALSE],
        na.action = stats::na.pass
      )
      ok_miss <- stats::complete.cases(mf_miss)
      if (!any(ok_miss)) {
        next
      }

      mf_ok <- mf_miss[ok_miss, , drop = FALSE]

      mf_pred_full <- mf_ok
      mf_pred_full[[tg]] <- rep.int(0, nrow(mf_pred_full))

      x_new <- stats::model.matrix(stats::formula(tt), data = mf_pred_full)
      cn <- colnames(x_fit)
      if (!all(cn %in% colnames(x_new))) {
        next
      }
      x_new <- x_new[, cn, drop = FALSE]

      dd <- .draw_imputations(fh_fit, meth, x_new, target_raw)
      ix <- seq_len(nrow(y_imp))[is_mis][ok_miss]
      y_imp[[tg]][ix] <- dd
    }
  }

  y_imp
}


survey_mi_fit <- function(data,
                          weights_rep,
                          weight_col_idx,
                          frm_by_target,
                          method_map,
                          tgt_set,
                          n_rounds,
                          n_imputations,
                          fitters,
                          parallel = FALSE,
                          n_cores = NULL,
                          block_info = NULL) {
  stopifnot(is.matrix(weights_rep), nrow(weights_rep) == nrow(data))
  nw <- ncol(weights_rep)
  stopifnot(
    length(weight_col_idx) == n_imputations,
    all(as.integer(weight_col_idx) >= 1L &
        as.integer(weight_col_idx) <= nw)
  )
  weight_col_idx <- as.integer(weight_col_idx)

  preds_pool <- names(data)[!names(data) %in% tgt_set]

  fitters_map <- .resolve_fitters(fitters)

  frm_by_target <- frm_by_target[tgt_set]

  mc <- vapply(data[tgt_set], function(x) sum(is.na(x)), integer(1))

  vars_ord <- tgt_set[
    order(mc[tgt_set] / nrow(data), tgt_set)
  ]

  ## quick exit
  if (!length(vars_ord)) {
    return(replicate(n_imputations, data, simplify = FALSE))
  }

  use_parallel <- isTRUE(parallel) &&
    !is.null(block_info) &&
    isTRUE(block_info$B > 1L) &&
    .survey_mi_fork_parallel_ok()

  if (isTRUE(parallel) && !is.null(block_info) && block_info$B > 1L &&
      !.survey_mi_fork_parallel_ok()) {
    warning(
      "`parallel = TRUE` requires Unix/macOS forking; running sequentially.",
      call. = FALSE
    )
  }

  if (use_parallel) {
    B_blk <- as.integer(block_info$B)
    A_blk <- as.integer(block_info$A)
    stopifnot(B_blk * A_blk == n_imputations)
    ncr <- .survey_mi_n_cores(n_cores)
    blocks <- parallel::mclapply(
      seq_len(B_blk),
      function(b) {
        tryCatch(
          .survey_mi_fit_one_b(
            b,
            data,
            weights_rep,
            frm_by_target,
            method_map,
            tgt_set,
            n_rounds,
            fitters_map,
            vars_ord,
            preds_pool,
            A_blk
          ),
          error = function(e) {
            list(.worker_error = conditionMessage(e), .b = b)
          }
        )
      },
      mc.cores = ncr,
      mc.set.seed = FALSE,
      mc.preschedule = FALSE
    )
    inv <- .survey_mi_parallel_blocks_invalid_fit(blocks, A_blk)
    if (any(inv)) {
      msg <- .survey_mi_parallel_failure_msg(blocks, inv)
      warning(
        "Parallel `survey_mi` workers failed; re-running sequentially.\n",
        msg,
        call. = FALSE
      )
    } else {
      out <- vector("list", n_imputations)
      idx <- 0L
      for (bb in seq_len(B_blk)) {
        for (aa in seq_len(A_blk)) {
          idx <- idx + 1L
          out[[idx]] <- blocks[[bb]][[aa]]
        }
      }
      return(out)
    }
  }

  out <- vector("list", n_imputations)

  for (m in seq_len(n_imputations)) {
    wcol <- weights_rep[, weight_col_idx[m], drop = TRUE]

    out[[m]] <- .survey_mi_single_draw(
      data,
      wcol,
      frm_by_target,
      method_map,
      tgt_set,
      n_rounds,
      fitters_map,
      vars_ord,
      preds_pool
    )
  }

  out
}


.survey_mi_fit_one_b <- function(b,
                                 data,
                                 weights_rep,
                                 frm_by_target,
                                 method_map,
                                 tgt_set,
                                 n_rounds,
                                 fitters_map,
                                 vars_ord,
                                 preds_pool,
                                 A) {
  wcol <- weights_rep[, b, drop = TRUE]
  out_a <- vector("list", A)
  for (a in seq_len(A)) {
    out_a[[a]] <- .survey_mi_single_draw(
      data,
      wcol,
      frm_by_target,
      method_map,
      tgt_set,
      n_rounds,
      fitters_map,
      vars_ord,
      preds_pool
    )
  }
  out_a
}


.resolve_methods <- function(data, vars, methods = NULL) {
  out <- stats::setNames(vector("list", length(vars)), vars)
  supported <- c("gaussian", "binomial", "poisson")

  ## variable -> family char
  m_by_v <- stats::setNames(character(0), character(0))
  if (!is.null(methods)) {
    if (is.list(methods) && all(names(methods) %in% supported)) {
      for (mtd in names(methods)) {
        vs <- methods[[mtd]]
        if (inherits(vs, "formula")) {
          stop("`methods` list entries must name variables, not formulas (use `formulas=`).")
        }
        vs <- unique(as.character(vs))
        vs <- vs[!is.na(vs) & nzchar(vs)]
        overlap <- intersect(vs, names(m_by_v))
        if (length(overlap)) {
          stop("Duplicates across method groups: ",
               paste(overlap, collapse = ", "))
        }
        m_by_v[vs] <- mtd
      }
      badnms <- setdiff(names(m_by_v), names(data))
      if (length(badnms)) {
        stop("Unknown variables named in methods: ",
             paste(badnms, collapse = ", "))
      }
    } else if (is.character(methods) && !is.null(names(methods))) {
      bm <- setdiff(unique(unname(methods)), supported)
      if (length(bm)) {
        stop("Unsupported families: ",
             paste(bm, collapse = ", "))
      }
      ux <- setdiff(names(methods), names(data))
      if (length(ux)) stop("methods names not in data: ", paste(ux, collapse = ", "))
      m_by_v <- methods
    } else {
      stop("`methods`: use list(...) by family names or legacy `c(var='gaussian',)`.")
    }
  }

  for (v in vars) {
    if (v %in% names(m_by_v)) {
      mtd <- unname(m_by_v[[v]])
    } else {
      x <- data[[v]]
      ux <- unique(stats::na.omit(x))
      if (
        is.logical(x) ||
          (is.factor(x) && nlevels(x) == 2L) ||
          (
            is.numeric(x) &&
              length(ux) <= 2L &&
              all(ux %in% c(0, 1))
          )
      ) {
        mtd <- "binomial"
      } else if (is.integer(x) && length(ux) && all(ux >= 0L)) {
        mtd <- "poisson"
      } else {
        mtd <- "gaussian"
      }
    }
    if (!mtd %in% supported) stop("Unsupported `", mtd, "` for ", v)
    out[[v]] <- mtd
  }
  out
}


.resolve_fitters <- function(fitters = NULL) {
  native <- list(
    gaussian = get0("weighted_gaussian_cpp", mode = "function"),
    binomial = get0("weighted_binomial_cpp", mode = "function"),
    poisson = get0("weighted_poisson_cpp", mode = "function")
  )

  defs <- list(
    gaussian = if (is.function(native$gaussian)) {
      function(x, y, w) native$gaussian(x, y, w)
    } else {
      function(x, y, w) {
        ff <- stats::lm.wfit(x = x, y = y, w = w)
        ne <- sum(w > 0)
        p <- ncol(x)
        rss <- sum(w * (y - drop(x %*% ff$coefficients))^2)
        sg2 <- if (ne > p) rss / (ne - p) else max(rss, 1e-8)
        xtwx <- crossprod(x, x * w)
        vc <- tryCatch(sg2 * solve(xtwx), error = function(e) MASS::ginv(xtwx) * sg2)
        list(coef = ff$coefficients, sigma2 = sg2, vcov = vc)
      }
    },
    binomial = if (is.function(native$binomial)) {
      function(x, y, w) native$binomial(x, y, w)
    } else {
      function(x, y, w) {
        z <- stats::glm.fit(x = x, y = y, weights = w, family = stats::binomial())
        mu <- pmin(pmax(z$fitted.values, 1e-6), 1 - 1e-6)
        xtwx <- crossprod(x, x * (w * mu * (1 - mu)))
        vc <- tryCatch(solve(xtwx), error = function(e) MASS::ginv(xtwx))
        list(coef = z$coefficients, vcov = vc)
      }
    },
    poisson = if (is.function(native$poisson)) {
      function(x, y, w) native$poisson(x, y, w)
    } else {
      function(x, y, w) {
        z <- stats::glm.fit(x = x, y = y, weights = w, family = stats::poisson())
        mu <- pmax(z$fitted.values, 1e-6)
        xtwx <- crossprod(x, x * (w * mu))
        vc <- tryCatch(solve(xtwx), error = function(e) MASS::ginv(xtwx))
        list(coef = z$coefficients, vcov = vc)
      }
    }
  )

  if (is.null(fitters)) return(defs)

  for (kk in names(fitters)) {
    if (!kk %in% names(defs)) {
      stop("Unknown fitters slot ", kk)
    }
    defs[[kk]] <- fitters[[kk]]
  }
  defs
}


.draw_imputations <- function(fit, method, x_new, original_target) {
  beta <- fit$coef
  if (anyNA(beta)) {
    beta[is.na(beta)] <- 0
  }

  ## Do not double-perturb coefficients:
  ## bootstrap replicate weights already inject parameter uncertainty.
  bdraw <- as.numeric(beta)

  eta <- drop(x_new %*% bdraw)

  if (method == "gaussian") {
    s2 <- fit$sigma2
    if (!is.finite(s2) || s2 <= 0L) {
      s2 <- 1e-8
    }
    eta + stats::rnorm(length(eta), sd = sqrt(s2))
  } else if (method == "binomial") {
    p <- 1 / (1 + exp(-eta))
    p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
    yb <- stats::rbinom(length(p), 1L, p)
    if (is.factor(original_target)) {
      lv <- levels(original_target)
      return(factor(lv[yb + 1L], levels = lv))
    }
    if (is.logical(original_target)) {
      return(as.logical(yb))
    }
    yb
  } else {
    mu <- pmax(exp(eta), 1e-8)
    stats::rpois(length(mu), lambda = mu)
  }
}


.initialize_missing <- function(x) {
  if (!any(ii <- is.na(x))) return(x)

  oo <- stats::na.omit(x)

  x[ii] <- sample(oo, sum(ii), TRUE)
  x
}
