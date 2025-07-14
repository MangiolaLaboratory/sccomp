#' @keywords internal
#' @noRd
#'
.subset <- 		function(.data,	 .column)	{
  # Make col names
  .column <- enquo(.column)

  # Check if column present
  if(quo_names(.column) %in% colnames(.data) %>% all %>% `!`)
    stop("sccomp says: some of the .column specified do not exist in the input data frame.")

  .data %>%

    # Selecting the right columns
    select(	!!.column,	all_of(get_specific_annotation_columns(.data, !!.column)	)) %>%
    distinct()

}

#' @importFrom dplyr distinct_at
#' @importFrom rlang enquo
#' @importFrom magrittr equals
#' @importFrom dplyr select
#'
#' @keywords internal
#' @noRd
#'
#'
get_specific_annotation_columns <- function(.data, .col){


  # Comply with CRAN NOTES
  . <- NULL

  # Make col names
  .col <- enquo(.col)

  # x-annotation df
  n_x <- .data %>% distinct(across(all_of(quo_names(.col)))) %>% nrow

  # element wise columns
  .data %>%
    select(-!!.col) %>%
    colnames %>%
    map(
      ~ {
        condition <- .data %>%
          distinct(across(all_of(c(quo_names(.col), .x)))) %>%
          nrow %>%
          equals(n_x)
        
        if (condition) {
          .x
        } else {
          NULL
        }
      }
    ) %>%

    # Drop NULL
    {	(.)[lengths((.)) != 0]	} %>%
    unlist

}

#' Convert array of quosure (e.g. c(col_a, col_b)) into character vector
#'
#' @keywords internal
#'
#' @importFrom rlang quo_name
#' @importFrom rlang quo_squash
#'
#' @keywords internal
#' @noRd
#'
#' @param v A array of quosures (e.g. c(col_a, col_b))
#'
#' @return A character vector
quo_names <- function(v) {

  v <- quo_name(quo_squash(v))
  gsub('^c\\(|`|\\)$', '', v) %>%
    strsplit(', ') %>%
    unlist
}




