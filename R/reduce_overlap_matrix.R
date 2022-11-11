#'Reduces overlap matrix using repeated column and row names
#'
#'\code{reduce_overlap_matrix} allows generic reduction of a symmetric overlap
#'matrix by using repeated column and row names as group labels. New overlap
#'values i,j are the average of overlap values from all original elements m,n
#'where m is in group i and n is in group j.
#'
#'@param overlap_mat matrix of overlap coefficients between organisms with
#'  column and row names corresponding to group assignment.
#'@param grouping vector containing the group assignment of each column and row.
#'
#'@return reduced overlap matrix with column and row names set to group labels.
#'@importFrom assertthat assert_that
#'@importFrom tibble column_to_rownames
#'@importFrom tidyr pivot_wider
#'@importFrom dplyr group_by summarise mutate case_when arrange
#'@importFrom magrittr %>% 
#'@importFrom rlang .data
reduce_overlap_matrix <- function(overlap_mat, grouping) {

  #check consistency of overlap_mat
  assertthat::assert_that(nrow(overlap_mat) == length(grouping) &
                            ncol(overlap_mat) == length(grouping),
                    msg = "Dimensions of overlap matrix should match grouping length")
  
  #use group_by reduction to take mean between groups
  reduced_mat <- expand.grid(grouping, grouping,
                             stringsAsFactors = FALSE) %>%
    dplyr::mutate(Overlap = as.vector(overlap_mat)) %>%
    dplyr::group_by(.data$Var1, .data$Var2) %>%
    dplyr::summarise(Overlap = mean(.data$Overlap),
                     .groups = "drop") %>%
    dplyr::mutate(
      Overlap = dplyr::case_when(.data$Var1 == .data$Var2 ~ 1.,
                                 .data$Var1 != .data$Var2 ~ .data$Overlap)) %>%
    dplyr::arrange(.data$Var1,.data$Var2) %>%
    tidyr::pivot_wider(
      id_cols = .data$Var1,
      names_from = .data$Var2,
      values_from = .data$Overlap,
      values_fill = 0.
    ) %>%
    tibble::column_to_rownames("Var1") %>%
    as.matrix()
  
  return(reduced_mat)
}
