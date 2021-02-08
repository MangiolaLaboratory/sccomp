
#' Formula parser
#'
#' @param fm A formula
#'
#' @return A character vector
#'
#'
parse_formula <- function(fm) {
	if (attr(terms(fm), "response") == 1)
		stop("The formula must be of the kind \"~ covariates\" ")
	else
		as.character(attr(terms(fm), "variables"))[-1]
}

#' Get matrix from tibble
#'
#' @import dplyr
#' @importFrom tidyr gather
#' @importFrom magrittr set_rownames
#'
#' @param tbl A tibble
#' @param rownames A character string of the rownames
#'
#' @return A matrix
as_matrix <- function(tbl, rownames = NULL) {
	tbl %>%

		ifelse_pipe(
			tbl %>%
				ifelse_pipe(!is.null(rownames),		~ .x %>% dplyr::select(-contains(rownames))) %>%
				summarise_all(class) %>%
				gather(variable, class) %>%
				pull(class) %>%
				unique() %>%
				`%in%`(c("numeric", "integer")) %>% `!`() %>% any(),
			~ {
				warning("to_matrix says: there are NON-numerical columns, the matrix will NOT be numerical")
				.x
			}
		) %>%
		as.data.frame() %>%

		# Deal with rownames column if present
		ifelse_pipe(!is.null(rownames),
								~ .x  %>%
									set_rownames(tbl %>% pull(!!rownames)) %>%
									select(-!!rownames)) %>%

		# Convert to matrix
		as.matrix()
}

#' vb_iterative
#'
#' @description Runs iteratively variational bayes until it suceeds
#'
#' @importFrom rstan vb
#'
#' @param model A Stan model
#' @param output_samples An integer of how many samples from posteriors
#' @param iter An integer of how many max iterations
#' @param tol_rel_obj A real
#' @param additional_parameters_to_save A character vector
#' @param ... List of paramaters for vb function of Stan
#'
#' @return A Stan fit object
#'
vb_iterative = function(model,
												output_samples,
												iter,
												tol_rel_obj,
												additional_parameters_to_save = c(),
												data,
												...) {
	res = NULL
	i = 0
	while (res %>% is.null | i > 5) {
		res = tryCatch({
			my_res = vb(
				model,
				data = data,
				output_samples = output_samples,
				iter = iter,
				tol_rel_obj = tol_rel_obj,
				#seed = 654321,
				#pars=c("counts_rng", "exposure_rate", "alpha_sub_1", additional_parameters_to_save),
				...
			)
			boolFalse <- T
			return(my_res)
		},
		error = function(e) {
			i = i + 1
			writeLines(sprintf("Further attempt with Variational Bayes: %s", e))
			return(NULL)
		},
		finally = {
		})
	}

	return(res)
}

#' function to pass initialisation values
#'
#' @return A list
inits_fx =
	function () {
		pars =
			res_discovery %>%
			filter(`.variable` != "counts_rng") %>%
			distinct(`.variable`) %>%
			pull(1)

		foreach(
			par = pars,
			.final = function(x)
				setNames(x, pars)
		) %do% {
			res_discovery %>%
				filter(`.variable` == par) %>%
				mutate(init = rnorm(n(), mean, sd)) %>%
				mutate(init = 0) %>%
				select(`.variable`, S, G, init) %>%
				pull(init)
		}
	}

#' Produce generated quantities plots with marked outliers
#'
#' @importFrom purrr pmap
#' @importFrom purrr map_int
#' @importFrom purrr when
#' @import ggplot2
#'
#' @param .x A tibble
#' @param symbol A symbol object
#' @param .count A symbol object
#' @param .sample A symbol object
#' @param covariate A character string
#'
#' @return A ggplot
produce_plots = function(.x,
												 symbol,
												 .count,
												 .sample,
												 covariate) {
	# Set plot theme
	my_theme =
		theme_bw() +
		theme(
			panel.border = element_blank(),
			axis.line = element_line(),
			panel.grid.major = element_line(size = 0.2),
			panel.grid.minor = element_line(size = 0.1),
			text = element_text(size = 12),
			#aspect.ratio = 1,
			axis.text.x = element_text(
				angle = 90,
				hjust = 1,
				vjust = 0.5
			),
			strip.background = element_blank(),
			axis.title.x  = element_text(margin = margin(
				t = 10,
				r = 10,
				b = 10,
				l = 10
			)),
			axis.title.y  = element_text(margin = margin(
				t = 10,
				r = 10,
				b = 10,
				l = 10
			))
		)

	max_y  = .x %>% summarise(a = max(!!as.symbol(.count)), b = max(.upper_2)) %>% as.numeric %>% max

	{
		ggplot(data = .x, aes(
			y = !!as.symbol(.count),
			x = !!as.symbol(.sample)
		))
		# +
		# 	geom_errorbar(
		# 		aes(ymin = `.lower_1`,
		# 				ymax = `.upper_1`),
		# 		width = 0,
		# 		linetype = "dashed",
		# 		color = "#D3D3D3"
		# 	)
	} %>%
		when(
			".lower_2" %in% colnames(.x) ~ (.) +
				geom_errorbar(aes(
					ymin = `.lower_2`,
					ymax = `.upper_2`,
					color = `deleterious outliers`
				),
				width = 0),
			~ (.)

		) %>%
		ifelse_pipe(
			covariate %>% is.null %>% `!`,
			~ .x + geom_point(aes(
				size = `exposure rate`, fill = !!as.symbol(covariate)
			), shape = 21),
			~ .x + geom_point(
				aes(size = `exposure rate`),
				shape = 21,
				fill = "black"
			)
		) +
		scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
		coord_cartesian(ylim = c(NA, max_y)) +
		my_theme +
		ggtitle(symbol)
}

#' fit_to_counts_rng
#'
#' @importFrom tidyr separate
#' @importFrom tidyr nest
#' @importFrom rstan summary
#'
#' @param fit A fit object
#' @param adj_prob_theshold fit real
#'
fit_to_counts_rng = function(fit, adj_prob_theshold){

	writeLines(sprintf("executing %s", "fit_to_counts_rng"))

	fit %>%
		rstan::summary("counts_rng",
									 prob = c(adj_prob_theshold, 1 - adj_prob_theshold)) %$%
		summary %>%
		as_tibble(rownames = ".variable") %>%
		separate(.variable,
						 c(".variable", "S", "G"),
						 sep = "[\\[,\\]]",
						 extra = "drop") %>%
		mutate(S = S %>% as.integer, G = G %>% as.integer) %>%
		select(-one_of(c("n_eff", "Rhat", "khat"))) %>%
		rename(`.lower` = (.) %>% ncol - 1,
					 `.upper` = (.) %>% ncol)
}

#' draws_to_tibble_x_y
#'
#' @importFrom tidyr pivot_longer
draws_to_tibble_x_y = function(fit, par, x, y) {

	par_names = names(fit) %>% grep(sprintf("%s", par), ., value = T)

	fit %>%
		rstan::extract(par_names, permuted=F) %>%
		as.data.frame %>%
		as_tibble() %>%
		mutate(.iteration = 1:n()) %>%
		pivot_longer(
		  names_to = c("dummy", ".chain", ".variable", x, y),
		  cols = contains(par),
		  names_sep = "\\.|\\[|,|\\]|:",
		  names_ptypes = list(
		    ".variable" = character()),
		  values_to = ".value"
		) %>%
	  mutate(
	    !!as.symbol(x) := as.integer(!!as.symbol(x)),
	    !!as.symbol(y) := as.integer(!!as.symbol(y))
	    ) %>%
		select(-dummy) %>%
		arrange(.variable, !!as.symbol(x), !!as.symbol(y), .chain) %>%
		group_by(.variable, !!as.symbol(x), !!as.symbol(y)) %>%
		mutate(.draw = 1:n()) %>%
		ungroup() %>%
		select(!!as.symbol(x), !!as.symbol(y), .chain, .iteration, .draw ,.variable ,     .value) %>%
	  filter(.variable == par)

}

draws_to_tibble_x = function(fit, par, x) {

	par_names = names(fit) %>% grep(sprintf("%s", par), ., value = T)

	fit %>%
		rstan::extract(par_names, permuted=F) %>%
		as.data.frame %>%
		as_tibble() %>%
		mutate(.iteration = 1:n()) %>%
		pivot_longer(names_to = c("dummy", ".chain", ".variable", x),  cols = contains(par), names_sep = "\\.|\\[|,|\\]|:", names_ptypes = list(".chain" = integer(), ".variable" = character(), "A" = integer(), "C" = integer()), values_to = ".value") %>%
		select(-dummy) %>%
		arrange(.variable, !!as.symbol(x), .chain) %>%
		group_by(.variable, !!as.symbol(x)) %>%
		mutate(.draw = 1:n()) %>%
		ungroup() %>%
		select(!!as.symbol(x), .chain, .iteration, .draw ,.variable ,     .value)

}

summary_to_tibble = function(fit, par, x, y = NULL) {

	par_names = names(fit) %>% grep(sprintf("%s", par), ., value = T)

	fit %>%
		rstan::summary(par_names) %$%
		summary %>%
		as_tibble(rownames = ".variable") %>%
		when(
			is.null(y) ~ (.) %>% tidyr::extract(col = .variable, into = c(".variable", x), "(.+)\\[(.+)\\]", convert = T),
			~ (.) %>% tidyr::extract(col = .variable, into = c(".variable", x, y), "(.+)\\[(.+),(.+)\\]", convert = T)
		) %>%
	  filter(.variable == "beta")

}

#' @importFrom tibble enframe
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
generate_quantities = function(fit, N, M, exposure){


  rstan::gqs(
    stanmodels$generated_quantities,
    #rstan::stan_model("inst/stan/generated_quantities.stan"),
    draws =  as.matrix(fit)
  ) %>%

    rstan::extract("counts") %$% counts %>%
    as.data.frame() %>%
    as_tibble(rownames = "draw") %>%
    gather(N_M, generated_quantity, -draw) %>%

    nest(data = -N_M) %>%
    mutate(quantiles = map(
      data,
      ~ quantile(
        .x$generated_quantity,
        probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
      ) %>%
        enframe() %>%
        spread(name, value)
    )) %>%

    separate(N_M, c("N", "M")) %>%
    mutate(N = as.integer(N), M = as.integer(M)) %>%

    select(-data) %>%
    unnest(quantiles)


}

#' do_inference
#'
#' @description This function calls the stan model.
#'
#' @importFrom tibble tibble
#' @import rstan
#' @import dplyr
#' @importFrom tidyr spread
#' @import tidybayes
#' @importFrom foreach foreach
#' @importFrom foreach %do%
#' @importFrom magrittr %$%
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom purrr map2
#' @importFrom purrr map_int
#' @importFrom tidybulk scale_abundance
#' @importFrom tidybayes gather_draws
#'
#' @param my_df A tibble including a cell_type name column | sample name column | read counts column | covariates column
#' @param formula A formula
#' @param .sample A column name as symbol
#' @param .cell_type A column name as symbol
#' @param .count A column name as symbol
#' @param .significance A column name as symbol
#' @param .do_check A column name as symbol
#' @param approximate_posterior_inference A boolean
#' @param approximate_posterior_analysis A boolean
#' @param C An integer
#' @param X A tibble
#' @param lambda_mu_mu A real
#' @param cores An integer
#' @param exposure_rate_multiplier A real
#' @param intercept_shift_scale A real
#' @param additional_parameters_to_save A character vector
#' @param adj_prob_theshold A real
#' @param how_many_posterior_draws A real number of posterior draws needed
#' @param to_exclude A boolean
#' @param truncation_compensation A real
#' @param save_generated_quantities A boolean
#' @param inits_fx A function
#' @param prior_from_discovery A tibble
#' @param pass_fit A fit
#' @param tol_rel_obj A real
#' @param write_on_disk A boolean
#' @param seed an integer
#'
#' @return A tibble with additional columns
#'
do_inference = function(.data,
												formula,
												approximate_posterior_inference = F,
												approximate_posterior_analysis = F,
												C,
												X,
												cores,
												additional_parameters_to_save,
												to_exclude = tibble(N = integer(), M = integer()),
												truncation_compensation = 1,
												save_generated_quantities = F,
												inits_fx = "random",
												prior_from_discovery = tibble(`.variable` = character(),
																											mean = numeric(),
																											sd = numeric()),
												pass_fit = F,
												tol_rel_obj = 0.01,
												write_on_disk = F,
												seed,
												output_samples= output_samples) {




	# # if analysis approximated
	# # If posterior analysis is approximated I just need enough
	# how_many_posterior_draws_practical = ifelse(approximate_posterior_analysis, 1000, how_many_posterior_draws)
	# additional_parameters_to_save = additional_parameters_to_save %>% c("lambda_log_param", "sigma_raw") %>% unique

  how_namy_to_exclude = to_exclude %>% nrow

  fit =
    vb_iterative(
      stanmodels$glm_dirichlet_multinomial,
      #pcc_seq_model, #
      output_samples = output_samples,
      iter = 5000,
      tol_rel_obj = 0.01,
      data = list(
        N = nrow(.data),
        M = ncol(.data)-1,
        y = .data %>% dplyr::select(-N),
        X = X
      )
    )

	# fit =
	#   sampling(
	#     stanmodels$glm_dirichlet_multinomial,
	#     #stan_model("glm_dirichlet_multinomial.stan"),
	#     data = list(
	#       N = nrow(.data),
	#       M = ncol(.data)-1,
	#       y = .data %>% dplyr::select(-N),
	#       X = X
	#     ),
	#     cores = 4
	#     #, iter = 5000, warmup = 300
	#   )

	# # Plot results
	# fit %>%
	#   draws_to_tibble_x_y("beta", "C", "M")
	#
	#   tidybayes::gather_draws(beta[C, M]) %>%
	#   median_qi() %>%
	#   filter(C==2) %>%
	#   bind_cols(cell_type = colnames(.data)[-1]) %>%
	#   ggplot(aes(forcats::fct_reorder(cell_type, .value), .value)) +
	#   geom_point() +
	#   geom_errorbar(aes(ymin = .lower, ymax =.upper)) +
	#   geom_hline(yintercept = 0) +
	#   theme_bw() +
	#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  fit

}




do_inference_imputation = function(.data,
                        formula,
                        approximate_posterior_inference = F,
                        approximate_posterior_analysis = F,
                        C,
                        X,
                        cores,
                        additional_parameters_to_save,
                        to_include = tibble(N = integer(), M = integer()),
                        truncation_compensation = 1,
                        save_generated_quantities = F,
                        inits_fx = "random",
                        prior_from_discovery = tibble(`.variable` = character(),
                                                      mean = numeric(),
                                                      sd = numeric()),
                        pass_fit = F,
                        tol_rel_obj = 0.01,
                        write_on_disk = F,
                        seed,
                        precision = precision,
                        exposure = exposure) {




  # # if analysis approximated
  # # If posterior analysis is approximated I just need enough
  # how_many_posterior_draws_practical = ifelse(approximate_posterior_analysis, 1000, how_many_posterior_draws)
  # additional_parameters_to_save = additional_parameters_to_save %>% c("lambda_log_param", "sigma_raw") %>% unique





  # Correct for 0 prop ##############################
  ###################################################

  #https://www.rdocumentation.org/packages/DirichletReg/versions/0.3-0/topics/DR_data
  fix_zeros = function(proportions){
    ( proportions*(nrow(proportions)-1) + (1/ncol(proportions)) ) / nrow(proportions)
  }
  # proportions[proportions==0] = min(proportions[proportions>0])
  # proportions = proportions/apply(proportions,1,sum)

  ###################################################
  ###################################################

  # Convert to log ratios
  .data_clr =
    .data %>%
    tidybulk::as_matrix(rownames = N) %>%
    apply(1, function(x) x/sum(x)) %>%
    t() %>%
    fix_zeros() %>%
    as.data.frame() %>%
    as_tibble() %>%
    rowid_to_column("N") %>%
    gather(M, proportions, -N) %>%
    mutate(M = as.integer(M)) %>%
    group_by(N) %>%
    mutate(centered_log_ratio = proportions %>% boot::logit() %>% scale(scale = F) %>% as.numeric) %>%
    ungroup(N) %>%
    select(-proportions) %>%
    spread(M, centered_log_ratio)


  how_namy_to_include = to_include %>% nrow
  I  = precision %>% nrow

  fit =
    vb_iterative(
      #stanmodels$glm_imputation,
      stan_model("inst/stan/glm_imputation.stan"),
      output_samples = 2000,
      iter = 5000,
      tol_rel_obj = 0.01,
      data = list(
        N = nrow(.data_clr),
        M = ncol(.data_clr)-1,
        y = .data_clr %>% dplyr::select(-N),
        X = X
      )
    )

  # fit =
  #   sampling(
  #     stanmodels$glm_dirichlet_multinomial,
  #     #stan_model("glm_dirichlet_multinomial.stan"),
  #     data = list(
  #       N = nrow(.data),
  #       M = ncol(.data)-1,
  #       y = .data %>% dplyr::select(-N),
  #       X = X
  #     ),
  #     cores = 4
  #     #, iter = 5000, warmup = 300
  #   )

  # # Plot results
  # fit %>%
  #   draws_to_tibble_x_y("beta", "C", "M")
  #
  #   tidybayes::gather_draws(beta[C, M]) %>%
  #   median_qi() %>%
  #   filter(C==2) %>%
  #   bind_cols(cell_type = colnames(.data)[-1]) %>%
  #   ggplot(aes(forcats::fct_reorder(cell_type, .value), .value)) +
  #   geom_point() +
  #   geom_errorbar(aes(ymin = .lower, ymax =.upper)) +
  #   geom_hline(yintercept = 0) +
  #   theme_bw() +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  fit




}
