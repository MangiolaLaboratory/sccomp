# Smooth-term utilities for sccomp formulas (brms-style splines)
#
# This file implements the R-side machinery that turns smooth terms
# (`s()`, `t2()`) in a sccomp formula into design-matrix columns that the
# existing Stan model can consume without modification.
#
# The decomposition is the standard mgcv / brms one:
#   * `Xf` — unpenalised null-space columns of the smooth (linear part).
#     These are appended to the fixed-effect design matrix `X` and
#     estimated as ordinary `beta` coefficients.
#   * `Xr` — penalised "wiggly" basis columns. These occupy a sccomp
#     random-effect slot with `n_factors = 1` and `n_groups = ncol(Xr)`,
#     so the slot's per-cell-group SD plays the role of brms's `sds_*`.
#
# Built via:
#   1. `mgcv::smoothCon(s(...), data, absorb.cons = TRUE,
#                       diagonal.penalty = TRUE)[[1]]` — basis + penalty.
#   2. `mgcv::smooth2random(sm, vnames = "", type = 2)` — re-parameterise
#      into `Xf` + standardised `rand$Xr` (so the wiggly coefficients have
#      a `N(0, sds * I)` prior).
#
# At prediction time, `predict_smooth_at_newdata()` evaluates the basis at
# new covariate values via `mgcv::PredictMat()` and re-applies the same
# `trans.D` (diagonal rescaling) and `trans.U` (orthonormal rotation)
# captured at fit time, so the new columns line up with the columns the
# model was trained on.

#' Detect smooth specials in a formula
#'
#' @param fm A one-sided formula.
#' @return TRUE if `fm` contains any `s(...)` or `t2(...)` term.
#' @keywords internal
#' @noRd
has_smooth_terms <- function(fm) {
  if (is.null(fm)) return(FALSE)
  trm <- stats::terms(fm, specials = c("s", "t2"))
  length(unlist(attr(trm, "specials"))) > 0
}


#' Strip smooth specials from a formula
#'
#' Returns a new formula with `s()` / `t2()` term labels removed. All other
#' terms, including random-effect clauses (`(... | g)`), are preserved. The
#' intercept setting (`+ 0` / `- 1`) is preserved. If nothing remains, returns
#' `~ 1` (or `~ 0`).
#'
#' @param fm A one-sided formula.
#' @return A one-sided formula without smooth term labels.
#' @keywords internal
#' @noRd
strip_smooth_terms <- function(fm) {
  if (is.null(fm)) return(fm)
  trm <- stats::terms(fm, specials = c("s", "t2"))
  smooth_idx <- unlist(attr(trm, "specials"))
  
  vars <- attr(trm, "variables")
  smooth_labels <- if (length(smooth_idx) > 0) {
    vapply(
      smooth_idx + 1L,
      function(i) deparse(vars[[i]], width.cutoff = 500L),
      character(1)
    )
  } else {
    character(0)
  }
  
  all_labels <- attr(trm, "term.labels")
  keep_labels <- setdiff(all_labels, smooth_labels)
  keep_labels <- ifelse(grepl("\\|", keep_labels), paste0("(", keep_labels, ")"), keep_labels)
  
  has_intercept <- attr(trm, "intercept") == 1L
  
  if (length(keep_labels) == 0) {
    return(if (has_intercept) ~ 1 else ~ 0)
  }
  
  rhs <- paste(keep_labels, collapse = " + ")
  if (!has_intercept) rhs <- paste0(rhs, " + 0")
  stats::as.formula(paste("~", rhs), env = environment(fm))
}


#' Strip random-effect clauses from a formula
#'
#' Returns a formula with `(... | g)` term labels removed, preserving all
#' ordinary fixed-effect terms and the intercept setting. This is used just
#' before calling `model.matrix()`, which cannot evaluate sccomp/brms-style
#' random-effect clauses.
#'
#' @param fm A one-sided formula.
#' @return A one-sided formula without random-effect term labels.
#' @keywords internal
#' @noRd
strip_random_effect_terms <- function(fm) {
  if (is.null(fm)) return(fm)
  trm <- stats::terms(fm)
  
  keep_labels <- attr(trm, "term.labels")
  keep_labels <- keep_labels[!grepl("\\|", keep_labels)]
  
  has_intercept <- attr(trm, "intercept") == 1L
  
  if (length(keep_labels) == 0) {
    return(if (has_intercept) ~ 1 else ~ 0)
  }
  
  rhs <- paste(keep_labels, collapse = " + ")
  if (!has_intercept) rhs <- paste0(rhs, " + 0")
  stats::as.formula(paste("~", rhs), env = environment(fm))
}


#' Parse smooth terms from a sccomp formula
#'
#' For each `s(...)` / `t2(...)` term, builds an mgcv smooth via
#' `smoothCon(..., absorb.cons = TRUE, diagonal.penalty = TRUE)` and the
#' brms-style re-parameterisation via `smooth2random(..., type = 2)`. Single-
#' penalty smooths (e.g. plain `s(x)`) yield one penalised block; multi-
#' penalty smooths (e.g. `t2(x, y)`, `s(x, z, bs = "fs")`) yield 2+ blocks,
#' one per penalty. Each block is mapped to its own sccomp random-effect
#' slot, with its own per-cell-group smoothing SD — matching the way brms
#' generates one `sds_*` parameter per penalty.
#'
#' @param fm A one-sided formula, possibly containing `s()` / `t2()` terms.
#' @param data A data frame containing every variable referenced inside the
#'   smooth terms (in the same row order as the design matrices that will be
#'   built downstream).
#' @return A list with components:
#'   * `parametric_formula` — `fm` with smooth specials stripped.
#'   * `smooth_labels`      — character vector of original smooth labels
#'     (one per `s()` / `t2()` term).
#'   * `Xf_list`            — list of `N × Ks_k` matrices (unpenalised cols),
#'     one per smooth term.
#'   * `Xr_list`            — **flat** list of `N × k_b` matrices, one per
#'     penalty block. A single-penalty smooth contributes 1 entry; a multi-
#'     penalty smooth contributes 2+. Each consumes one RE slot.
#'   * `Xr_to_smooth`       — integer vector parallel to `Xr_list` mapping
#'     each block to its parent smooth index.
#'   * `Xr_slot_labels`     — character vector parallel to `Xr_list` giving
#'     each block its sccomp RE-slot label: just the smooth's label for
#'     single-penalty smooths, or `<label>__b<b>` for multi-penalty smooths
#'     (so per-block `gfi` rownames stay unique).
#'   * `smooth_specs`       — list of mgcv `smoothCon` objects (for predict).
#'   * `smooth_re_objs`     — list of mgcv `smooth2random` objects (rotations).
#'   When `fm` has no smooth specials, lists are length-0 and
#'   `parametric_formula` is `fm` unchanged.
#'
#' @keywords internal
#' @noRd
parse_formula_smooths <- function(fm, data) {
  empty <- list(
    parametric_formula = fm,
    smooth_labels      = character(0),
    Xf_list            = list(),
    Xr_list            = list(),
    Xr_to_smooth       = integer(0),
    Xr_slot_labels     = character(0),
    smooth_specs       = list(),
    smooth_re_objs     = list()
  )
  if (!has_smooth_terms(fm)) return(empty)
  
  check_and_install_packages("mgcv")
  
  trm <- stats::terms(fm, specials = c("s", "t2"))
  vars <- attr(trm, "variables")
  # +1 to skip the `list` call head when indexing into `vars`.
  smooth_idx <- unlist(attr(trm, "specials")) + 1L
  
  smooth_calls  <- lapply(smooth_idx, function(i) vars[[i]])
  smooth_labels <- vapply(
    smooth_calls,
    function(e) deparse(e, width.cutoff = 500L),
    character(1)
  )
  smooth_labels <- unname(smooth_labels)
  
  Xf_list        <- vector("list", length(smooth_calls))
  smooth_specs   <- vector("list", length(smooth_calls))
  smooth_re_objs <- vector("list", length(smooth_calls))
  Xr_list        <- list()   # flat across smooths × blocks
  Xr_to_smooth   <- integer(0)
  Xr_slot_labels <- character(0)
  
  for (k in seq_along(smooth_calls)) {
    label <- smooth_labels[k]
    
    # Evaluate the s()/t2() call against the data frame. `s` and `t2` are
    # bound in a small lookup env so the call resolves even when the
    # caller's formula environment doesn't have mgcv attached.
    fm_env  <- environment(fm)
    if (is.null(fm_env)) fm_env <- baseenv()
    eval_env <- new.env(parent = fm_env)
    eval_env$s  <- mgcv::s
    eval_env$t2 <- mgcv::t2
    sm_spec <- eval(smooth_calls[[k]], envir = as.list(data), enclos = eval_env)
    
    # absorb.cons = TRUE removes the constant function from the basis so the
    # smooth doesn't fight the intercept (matches brms).
    # diagonal.penalty = TRUE re-parameterises so the penalty matrix is
    # diagonal, which makes the random-effect prior N(0, sds^2 I) clean.
    sm <- tryCatch(
      mgcv::smoothCon(
        sm_spec,
        data            = as.data.frame(data),
        absorb.cons     = TRUE,
        diagonal.penalty = TRUE
      )[[1]],
      error = function(e) stop(sprintf(
        "sccomp says: failed to build smooth basis for `%s`: %s",
        label, conditionMessage(e)
      ))
    )
    
    re <- mgcv::smooth2random(sm, vnames = "", type = 2)
    
    Xf <- if (is.null(re$Xf)) matrix(0, nrow = nrow(data), ncol = 0) else re$Xf
    attr(Xf, "s.label") <- NULL
    if (ncol(Xf) > 0) {
      colnames(Xf) <- sprintf("%s__lin%d", label, seq_len(ncol(Xf)))
    }
    
    # Pull all penalty blocks (`re$rand$Xr`, `re$rand$Xr.0`, `re$rand$Xr.1`,
    # ...). Each one becomes its own sccomp RE slot.
    if (is.null(re$rand) || length(re$rand) == 0L) {
      stop(sprintf(
        "sccomp says: smooth `%s` produced no penalised basis (unsupported basis type?).",
        label
      ))
    }
    block_names    <- names(re$rand)
    n_blocks_smooth <- length(block_names)
    
    for (b in seq_along(block_names)) {
      Xr_b <- re$rand[[block_names[b]]]
      attr(Xr_b, "s.label") <- NULL
      # Single-block smooths keep the historical `<label>` slot label (and the
      # matching `<label>___basis<jj>` column naming); multi-block smooths
      # disambiguate with `<label>__b<b>`.
      slot_label <- if (n_blocks_smooth == 1L) label
                    else sprintf("%s__b%d", label, b)
      colnames(Xr_b) <- sprintf("%s___basis%02d", slot_label, seq_len(ncol(Xr_b)))
      Xr_list[[length(Xr_list) + 1L]] <- Xr_b
      Xr_to_smooth                    <- c(Xr_to_smooth, k)
      Xr_slot_labels                  <- c(Xr_slot_labels, slot_label)
    }
    
    Xf_list[[k]]        <- Xf
    smooth_specs[[k]]   <- sm
    smooth_re_objs[[k]] <- re
  }
  
  list(
    parametric_formula = strip_smooth_terms(fm),
    smooth_labels      = smooth_labels,
    Xf_list            = Xf_list,
    Xr_list            = Xr_list,
    Xr_to_smooth       = Xr_to_smooth,
    Xr_slot_labels     = Xr_slot_labels,
    smooth_specs       = smooth_specs,
    smooth_re_objs     = smooth_re_objs
  )
}


#' Evaluate the per-level base basis (`Xb`) of a factor-smooth at new data
#'
#' `mgcv::smooth2random.fs.interaction()` ignores `object$X` and rebuilds the
#' basis from `object$Xb` (the per-level base basis) whenever `Xb` is non-NULL.
#' To predict at new data we must therefore recompute `Xb` ourselves, which is
#' the same logic `mgcv:::Predict.matrix.fs.interaction()` runs internally
#' before its level-by-level expansion. We replicate just that initial step.
#'
#' @keywords internal
#' @noRd
compute_fs_Xb_at_newdata <- function(sm, newdata) {
  newdata <- as.data.frame(newdata)
  fac_dropped <- newdata
  fac_dropped[[sm$fterm]] <- NULL
  
  obj <- sm
  class(obj)         <- sm$base$bs
  obj$rank           <- sm$base$rank
  obj$null.space.dim <- sm$base$null.space.dim
  obj$bs.dim         <- sm$base$bs.dim
  obj$term           <- sm$base$term
  mgcv::Predict.matrix(obj, fac_dropped) %*% sm$P
}


#' Evaluate a smooth basis at new covariate values
#'
#' Reproduces the (Xf, Xr_blocks) decomposition obtained at fit time, but on
#' `newdata`. The approach is *clone-and-rerun*: we feed `mgcv::smooth2random`
#' a shallow clone of the smooth object whose basis fields point at the new
#' data, then re-invoke the same `smooth2random(..., type = 2)` decomposition
#' that ran at fit time. Because `smooth2random` is purely a function of the
#' smooth's basis matrices and penalty / level structure, this reproduces the
#' fit-time `re$rand$Xr*` and `re$Xf` column-for-column when fed the original
#' training rows.
#'
#' Supported smooth families:
#'   * `tp.smooth`, `cr.smooth`, single-penalty smooths: clone-and-rerun on
#'     `$X = PredictMat(sm, newdata)`.
#'   * `fs.interaction` (`bs = "fs"`): factor smooths store a per-level base
#'     basis in `$Xb` that `smooth2random` reads in preference to `$X`. We
#'     recompute `Xb` at new rows via `compute_fs_Xb_at_newdata()` and also
#'     refresh `$fac` so the level masks `fac == flev[j]` line up.
#'
#' Tensor-product `t2.smooth` is not yet supported at prediction time because
#' `mgcv::PredictMat` applies the t2 constraint differently from `smoothCon`,
#' producing a basis that does not match the fit-time decomposition.
#'
#' @param sm      An mgcv smooth object (the result of `smoothCon()[[1]]`).
#' @param re_obj  The matching `smooth2random()` result (unused, kept for
#'   API parity with prior versions).
#' @param newdata A data frame containing every variable the smooth uses.
#' @return A list with `Xf` (`N_new × Ks`) and `Xr_blocks`, a list of
#'   `N_new × k_b` matrices (one per penalty block, in `re$rand` order).
#'
#' @keywords internal
#' @noRd
predict_smooth_at_newdata <- function(sm, re_obj, newdata) {
  check_and_install_packages("mgcv")
  if (inherits(sm, "t2.smooth")) {
    stop(
      "sccomp says: prediction with t2() smooths is not yet supported ",
      "(constraint absorption mismatch between PredictMat() and smoothCon()). ",
      "Use s(...) or s(..., bs = \"fs\") for now."
    )
  }
  
  newdata <- as.data.frame(newdata)
  
  sm_new <- sm
  if (inherits(sm, "fs.interaction")) {
    # fs smooths: smooth2random uses `$Xb` + `$base$S` rather than `$X` + `$S`,
    # so we refresh `$Xb` (per-level base basis at new rows) and `$fac`.
    fac_new <- newdata[[sm$fterm]]
    if (is.null(fac_new)) {
      stop(sprintf(
        "sccomp says: factor `%s` is missing from new data (required by smooth `%s`).",
        sm$fterm, sm$label
      ))
    }
    sm_new$Xb  <- compute_fs_Xb_at_newdata(sm, newdata)
    sm_new$fac <- factor(fac_new, levels = levels(sm$fac))
  } else {
    Xp <- mgcv::PredictMat(sm, data = newdata)
    attr(Xp, "s.label") <- NULL
    sm_new$X <- Xp
  }
  
  re_new <- mgcv::smooth2random(sm_new, vnames = "", type = 2)
  
  Xf_new <- if (is.null(re_new$Xf)) {
    matrix(0, nrow = nrow(newdata), ncol = 0)
  } else {
    as.matrix(re_new$Xf)
  }
  attr(Xf_new, "s.label") <- NULL
  
  Xr_blocks <- if (is.null(re_new$rand)) {
    list()
  } else {
    lapply(re_new$rand, function(m) {
      m <- as.matrix(m)
      attr(m, "s.label") <- NULL
      m
    })
  }
  
  list(Xf = Xf_new, Xr_blocks = Xr_blocks)
}


#' Wrap an Xr matrix into a sccomp random-effect slot
#'
#' Each smooth occupies exactly one slot with `n_factors = 1` and
#' `n_groups = ncol(Xr)`: the slot's per-cell-group SD then plays the role of
#' brms's smooth-wise `sds_*`. The slot structure (X / X_unseen / ncol / gfi /
#' n_groups / n_factors) mirrors what `prepare_re_slot()` returns for an
#' ordinary random-effect clause.
#'
#' @param Xr     A numeric matrix (N rows × K wiggly basis columns).
#' @param label  Smooth label, used as the slot's single "factor" name.
#' @return A list with the same fields as `empty_re_slot` in
#'   `data_spread_to_model_input()`.
#'
#' @keywords internal
#' @noRd
build_smooth_slot <- function(Xr, label) {
  stopifnot(is.matrix(Xr) || is.data.frame(Xr))
  Xr <- as.matrix(Xr)
  K  <- ncol(Xr)
  
  if (K == 0L) {
    return(list(
      X         = Xr,
      X_unseen  = Xr[, integer(0), drop = FALSE],
      ncol      = 0L,
      gfi       = matrix(integer(0), nrow = 0L, ncol = 0L),
      n_groups  = 0L,
      n_factors = 0L
    ))
  }
  
  # One factor (= the smooth), K "groups" (= the basis columns).
  # `gfi[1, j] = j` means basis-column j is at position j in X_random_effect_k.
  gfi <- matrix(
    seq_len(K),
    nrow = 1L, ncol = K,
    dimnames = list(label, colnames(Xr))
  )
  
  list(
    X         = Xr,
    X_unseen  = Xr[, integer(0), drop = FALSE],
    ncol      = K,
    gfi       = gfi,
    n_groups  = K,
    n_factors = 1L
  )
}


#' Build replicate-time smooth design pieces from fit-time metadata
#'
#' For every `s()` / `t2()` term recorded in `smooth_results`, evaluates the
#' basis at `new_data_tail` via [predict_smooth_at_newdata()] and packages the
#' two design contributions that `prepare_replicate_data()` needs:
#'
#' * the unpenalised (null-space) columns are renamed `<label>__lin*` and
#'   appended to `parametric_X` (so they line up positionally with the
#'   fit-time `beta` coefficients);
#' * each penalised block becomes one random-effect slot, using the same
#'   `<label>___basis*` (single-penalty) or `<label>__b<b>___basis*`
#'   (multi-penalty) column naming as fit time, so matching against the
#'   original `X_random_effect_*` matrices still succeeds.
#'
#' When `smooth_results` is empty (no smooth terms in the original formula)
#' the function returns `parametric_X` unchanged and an empty slot list.
#'
#' @param parametric_X Numeric matrix; the parametric (non-smooth) design
#'   matrix evaluated on the replicate rows. New unpenalised smooth columns
#'   are appended on the right.
#' @param new_data_tail Data frame of the replicate rows (must contain every
#'   variable referenced by any smooth term).
#' @param smooth_results A list with `smooth_specs`, `smooth_re_objs`, and
#'   `smooth_labels` (parallel lists captured at fit time), or `NULL` /
#'   length-0 to signal "no smooths".
#' @return A list with:
#'   * `new_X` — `parametric_X` with the unpenalised smooth columns appended.
#'   * `smooth_replicate_slots` — flat list of slot lists (one per penalised
#'     block, in fit-time order), each with `X`, `X_unseen`, `which`.
#'
#' @keywords internal
#' @noRd
build_smooth_replicate_design <- function(parametric_X,
                                          new_data_tail,
                                          smooth_results) {
  smooth_specs   <- smooth_results$smooth_specs
  smooth_re_objs <- smooth_results$smooth_re_objs
  smooth_labels  <- smooth_results$smooth_labels
  
  if (length(smooth_specs) == 0L) {
    return(list(new_X = parametric_X, smooth_replicate_slots = list()))
  }
  
  # One basis evaluation per smooth term, carrying the label alongside the
  # (Xf, Xr_blocks) pair so the downstream maps don't need parallel vectors.
  smooth_evals <- purrr::pmap(
    list(spec = smooth_specs, re = smooth_re_objs, label = smooth_labels),
    function(spec, re, label) {
      pred <- predict_smooth_at_newdata(spec, re, new_data_tail)
      pred$label <- label
      pred
    }
  )
  
  # Append unpenalised (Xf) columns in fit-time order.
  new_X <- purrr::reduce(smooth_evals, .init = parametric_X, function(X, e) {
    if (ncol(e$Xf) == 0L) return(X)
    colnames(e$Xf) <- sprintf("%s__lin%d", e$label, seq_len(ncol(e$Xf)))
    cbind(X, e$Xf)
  })
  
  # One slot per penalised block, flattened across smooths.
  smooth_replicate_slots <- smooth_evals |>
    purrr::map(function(e) {
      n_blocks <- length(e$Xr_blocks)
      purrr::map2(e$Xr_blocks, seq_along(e$Xr_blocks), function(Xr_b, b) {
        colnames(Xr_b) <- if (n_blocks == 1L) {
          sprintf("%s___basis%02d", e$label, seq_len(ncol(Xr_b)))
        } else {
          sprintf("%s__b%d___basis%02d", e$label, b, seq_len(ncol(Xr_b)))
        }
        list(
          X        = Xr_b,
          X_unseen = Xr_b[, integer(0), drop = FALSE],
          which    = seq_len(ncol(Xr_b)) |> as.array()
        )
      })
    }) |>
    unlist(recursive = FALSE)
  
  list(new_X = new_X, smooth_replicate_slots = smooth_replicate_slots)
}

#' Smooth metadata attached at fit time (R-only, not Stan data)
#'
#' @param x A fitted `sccomp_tbl` or its `model_input` list.
#' @keywords internal
#' @noRd
get_smooth_results <- function(x) {
  mi <- if (inherits(x, "sccomp_tbl")) attr(x, "model_input") else x
  if (is.null(mi)) return(NULL)
  attr(mi, "smooth_results")
}