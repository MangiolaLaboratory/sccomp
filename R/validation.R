check_columns_exist = function(.data, columns){

  if((!columns %in% (.data %>% colnames)) %>% any)
    stop(
      sprintf(
        "The columns %s are not present in your tibble",
        paste(columns[(!columns %in% (.data %>% colnames))], collapse=" ")
      )
    )
}

check_if_count_integer = function(.data, .count){
  .count = enquo(.count)

  # Check if the counts column is an integer
  if (.data %>% select(!!.count) %>% sapply(class) != "integer")
    stop(
      sprintf(
        "The column %s must be of class integer. You can do as mutate(`%s` = `%s` %%>%% as.integer)",
        quo_name(.count),
        quo_name(.count),
        quo_name(.count)
      )
    )
}

#' Check if NA
#'
#' @importFrom tidyr drop_na
#' @importFrom dplyr enquo
#'
#' @param .data A tibble including a gene name column | sample name column | read counts column | covariates column
#' @param .sample A column name
#' @param .cell_type A column name
#' @param .count A column name
#' @param .significance A column name
#' @param .do_check A column name
#' @param formula_columns A symbol vector
#'
check_if_any_NA = function(.data, columns){


  if(
    .data %>%
    drop_na(columns) %>%
    nrow %>% `<`
    (
      .data %>% nrow
    )
  )
    stop(sprintf("There are NA values in you tibble for any of the column %s", paste(columns, collapse=", ")))
}

check_if_within_posterior = function(.data, my_df, .do_check, .count){

  writeLines(sprintf("executing %s", "check_if_within_posterior"))

  .data %>%
    left_join(my_df, by = c("S", "G")) %>%
    filter((!!.do_check)) %>% # Filter only DE genes
    rowwise() %>%
    mutate(`ppc` = !!.count %>% between(`.lower`, `.upper`)) %>%
    mutate(`is higher than mean` = (!`ppc`) &
             (!!.count > mean)) %>%
    ungroup
}
