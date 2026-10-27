test_that("sccomp_test_smooth() errors on a fit without smooth terms", {
  skip_if_not_installed("mgcv")
  skip_cmdstan()
  data("counts_obj")

  fit <- sccomp_estimate(
    counts_obj,
    formula_composition = ~ type,
    formula_variability = ~ 1,
    sample              = "sample",
    cell_group          = "cell_group",
    abundance           = "count",
    cores               = 1,
    inference_method    = "pathfinder",
    max_sampling_iterations = 200,
    verbose = FALSE
  )

  expect_error(
    sccomp_test_smooth(fit, from = 1, to = 2),
    "no smooth"
  )
})

test_that("sccomp_test_smooth() endpoints contrast returns the standard columns", {
  skip_if_not_installed("mgcv")
  skip_cmdstan()
  data("counts_obj")

  samples <- levels(counts_obj$sample)
  covs <- tibble::tibble(
    sample     = samples,
    pseudotime = seq(0, 6, length.out = length(samples))
  )
  counts_pt <- counts_obj |> dplyr::left_join(covs, by = "sample")

  fit <- sccomp_estimate(
    counts_pt,
    formula_composition = ~ s(pseudotime, k = 4),
    formula_variability = ~ 1,
    sample              = "sample",
    cell_group          = "cell_group",
    abundance           = "count",
    cores               = 1,
    inference_method    = "pathfinder",
    max_sampling_iterations = 200,
    verbose = FALSE
  )

  res <- fit |>
    sccomp_test_smooth(from = 1, to = 5, number_of_draws = 100)

  n_groups <- counts_pt |> dplyr::distinct(cell_group) |> nrow()
  expect_equal(nrow(res), n_groups)

  expect_true(all(
    c("cell_group", "smooth", "from", "to",
      "c_lower", "c_effect", "c_upper", "c_pH0", "c_FDR") %in% colnames(res)
  ))

  # auto-detected smooth variable, correct labels
  expect_equal(unique(res$smooth), "pseudotime")
  expect_equal(unique(res$from), "1")
  expect_equal(unique(res$to), "5")

  # probabilities and FDR are in [0, 1]
  expect_true(all(res$c_pH0 >= 0 & res$c_pH0 <= 1))
  expect_true(all(res$c_FDR >= 0 & res$c_FDR <= 1))
  # CI brackets the point effect
  expect_true(all(res$c_lower <= res$c_effect & res$c_effect <= res$c_upper))
})

test_that("sccomp_test_smooth() average mode requires length-2 intervals", {
  skip_if_not_installed("mgcv")
  skip_cmdstan()
  data("counts_obj")

  samples <- levels(counts_obj$sample)
  covs <- tibble::tibble(
    sample     = samples,
    pseudotime = seq(0, 6, length.out = length(samples))
  )
  counts_pt <- counts_obj |> dplyr::left_join(covs, by = "sample")

  fit <- sccomp_estimate(
    counts_pt,
    formula_composition = ~ s(pseudotime, k = 4),
    formula_variability = ~ 1,
    sample              = "sample",
    cell_group          = "cell_group",
    abundance           = "count",
    cores               = 1,
    inference_method    = "pathfinder",
    max_sampling_iterations = 200,
    verbose = FALSE
  )

  # wrong shape for average
  expect_error(
    sccomp_test_smooth(fit, from = 1, to = 5, comparison = "average"),
    "length-2"
  )

  res <- fit |>
    sccomp_test_smooth(
      from = c(0, 2), to = c(4, 6),
      comparison = "average", resolution = 6, number_of_draws = 100
    )

  n_groups <- counts_pt |> dplyr::distinct(cell_group) |> nrow()
  expect_equal(nrow(res), n_groups)
  expect_equal(unique(res$from), "[0, 2]")
  expect_equal(unique(res$to), "[4, 6]")
})

test_that("sccomp_test_smooth() selects a curve via `at` for factor smooths", {
  skip_if_not_installed("mgcv")
  skip_cmdstan()
  data("counts_obj")

  samples <- levels(counts_obj$sample)
  covs <- tibble::tibble(
    sample     = samples,
    pseudotime = seq(0, 6, length.out = length(samples)),
    tissue     = factor(rep(c("blood", "lymph", "tumor"),
                            length.out = length(samples)))
  )
  counts_pt <- counts_obj |> dplyr::left_join(covs, by = "sample")

  fit <- sccomp_estimate(
    counts_pt,
    formula_composition = ~ s(pseudotime, tissue, bs = "fs", k = 5),
    formula_variability = ~ 1,
    sample              = "sample",
    cell_group          = "cell_group",
    abundance           = "count",
    cores               = 1,
    inference_method    = "pathfinder",
    max_sampling_iterations = 200,
    verbose = FALSE
  )

  res <- fit |>
    sccomp_test_smooth(
      smooth = "pseudotime", from = 1, to = 5,
      at = list(tissue = "tumor"), number_of_draws = 100
    )

  n_groups <- counts_pt |> dplyr::distinct(cell_group) |> nrow()
  expect_equal(nrow(res), n_groups)
  expect_equal(unique(res$smooth), "pseudotime")
})
