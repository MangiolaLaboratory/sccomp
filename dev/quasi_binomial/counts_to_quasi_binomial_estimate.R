

counts_to_quasi_binomial = function(.data, .sample, .cell_group, .count, .factor_of_interest){

  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .factor_of_interest = enquo(.factor_of_interest)
  .count = enquo(.count)

  .data %>%
    with_groups(!!.sample, ~ mutate(.x, n = sum(!!.count))) %>%
    mutate(proportion := !!.count/n) %>%
    rename(factor_of_interest := !!.factor_of_interest) %>%
    nest(data = -!!.cell_group) %>%
    mutate(estimate = map(
      data,
      ~ {
        fit = glm(proportion~factor_of_interest,family=quasibinomial(), weights=n, data=.x)

        fit %>%
        anova(test="F") %>%
        broom::tidy() %>%
        suppressWarnings() %>%
        filter(term !="NULL") %>%
          mutate(estimate = fit$coefficients[2])
      }
    )) %>%
    unnest(estimate)


}
