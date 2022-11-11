


#' Make heatmap by taxonomic group
#'
#' @param score_file_path path to csv with organism scores
#' @param reduced_matrix_path path to csv with reduced matrix
#' @param output_dir Path to the output directory, defaults to the current
#'   working directory.
#' @param output_type Extension associated with format to save the plot in. Must
#'   be one of "eps", "ps", "tex" (pictex), "pdf", "jpeg", "tiff", "png", "bmp",
#'   "svg" or "wmf".
#' @param min_pos_score minimim pos_score to filter groups by
#' @param score_suffix suffix to remove from the score_file_path to leave the
#'   dataset name, defaults to "_organism_scores.csv"
#' @param ... additional arguments passed to \code{\link[ggplot2]{ggsave}}
#'
#' @return returns a ggplot2 object
#'
#' @author Sarah C. Jenson
#'
#' @export
#'
#' @importFrom readr read_csv cols col_character col_double
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr select mutate filter distinct inner_join full_join
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom tidyr pivot_longer
#' @import ggplot2
group_heatmap <- function(score_file_path,
                          reduced_matrix_path,
                          output_dir = ".",
                          output_type = "jpeg",
                          min_pos_score = 1,
                          score_suffix = "_organism_scores.csv",
                          ...) {

  #checking inputs
  assertthat::assert_that(file.exists(score_file_path), msg = "score_file_path not found.")

  assertthat::assert_that(dir.exists(output_dir), msg = "output_dir not found.")

  assertthat::assert_that(any(output_type %in% c("eps", "ps", "tex",
                                                 "pdf", "jpeg", "tiff",
                                                 "png", "bmp", "svg", "wmf")),
                          msg = "invalid output_type, must be one of: eps, ps, tex,
                                                 pdf, jpeg, tiff,
                                                 png, bmp, svg, or wmf")

  assertthat::assert_that(is.numeric(min_pos_score) & (min_pos_score >= 0),
                          msg = "min_pos_score must be numeric and >= 0")

  assertthat::assert_that(grepl(score_suffix, score_file_path, fixed = TRUE),
                          msg = "score_suffix not found in score_file_path.")

  #getting dataset name
  dataset_name <- gsub(score_suffix, "", basename(score_file_path), fixed = TRUE)

  print(paste0("Starting ", dataset_name))

  #Loading inputs
  scores <- readr::read_csv(score_file_path, col_types = readr::cols(
    kegg_id = readr::col_character(),
    organism_id = readr::col_double(),
    kegg_org_code = readr::col_character(),
    organism_name = readr::col_character(),
    species_name = readr::col_character(),
    genus_name = readr::col_character(),
    peptide_count = readr::col_double(),
    protein_count = readr::col_double(),
    strong_peptide_count = readr::col_double(),
    taxonomic_group = readr::col_double(),
    organism_hits = readr::col_double(),
    group_hits = readr::col_double(),
    orgscores = readr::col_double(),
    pos_scores = readr::col_double()
  ))

  reduced_matrix <- readr::read_csv(reduced_matrix_path,
                                    col_types = readr::cols(
                                      .default = readr::col_double())) %>%
    tibble::column_to_rownames(var ="rowname")



  #filtering scores and reduced_matrix
  filt_scores <- scores %>%
    dplyr::filter(.data$pos_scores > min_pos_score)

  #making sure not empty
  if(length(unique(filt_scores$taxonomic_group)) <= 1){
    print(paste0("Only ",
                 length(unique(filt_scores$taxonomic_group)),
                 " groups passed filtering. Returning NA" ))
    return(NA)
  }

  #adding top group info and sorting by taxononmic group
  filt_scores <- filt_scores %>%
    dplyr::mutate(is_top_group = dplyr::if_else(.data$pos_scores ==
                                                  max(.data$pos_scores),
                                                TRUE, FALSE)) %>%
    dplyr::arrange(.data$taxonomic_group)

 #getting taxonomic groups
  tax_groups <- filt_scores %>%
    dplyr::pull(.data$taxonomic_group) %>%
    unique() %>%
    as.character()


  #verifying groups exist in matrix.  Need to check why this error keeps showing up!!!
  assertthat::assert_that(all(tax_groups %in% colnames(reduced_matrix)) &
                            all(tax_groups %in% rownames(reduced_matrix)),
                          msg = "All tax_groups not found in reduced_matrix")


  filt_group_matrix <- reduced_matrix[tax_groups,tax_groups]

  #making long version of score matrix
  long_overlap <- filt_group_matrix %>%
    tibble::rownames_to_column(var = "group_y") %>%
    tidyr::pivot_longer(cols = -.data$group_y,
                        names_to = "group_x",
                        values_to = "overlap") %>%
    dplyr::mutate(group_y = as.character(.data$group_y),
                  group_x = as.character(.data$group_x))

  #adding scores
  long_overlap <- filt_scores %>%
    dplyr::select(.data$taxonomic_group,
                  .data$is_top_group,
                  .data$pos_scores) %>%
    dplyr::mutate(taxonomic_group = as.character(.data$taxonomic_group)) %>%
    dplyr::mutate(pos_scores = format(.data$pos_scores, digits=NULL)) %>%
    dplyr::rename(group_x = .data$taxonomic_group) %>%
    dplyr::mutate(group_y = .data$group_x) %>%
    dplyr::full_join(long_overlap, by = c("group_y", "group_x"))


  #making groups ordered factors
  long_overlap <- long_overlap %>%
    dplyr::mutate(group_y = factor(.data$group_y, levels = tax_groups, ordered = TRUE),
           group_x = factor(.data$group_x, levels = tax_groups, ordered = TRUE))

  #making graph
  graph <- ggplot2::ggplot(long_overlap, aes(x = .data$group_x,
                                             y = .data$group_y,
                                             fill = .data$overlap)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(aes(label = .data$pos_scores,
                           color = .data$is_top_group),
                       na.rm = TRUE) +
    ggplot2::scale_color_manual(breaks = c(TRUE, FALSE), values = c("red", "black")) +
    ggplot2::guides(color = "none") +
    ggplot2::theme_bw() +
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::scale_y_discrete(expand = c(0,0)) +
    ggplot2::scale_fill_distiller(palette = "Blues", direction = -1) +
    ggplot2::ggtitle(dataset_name,
                     subtitle = paste0("Group level with minimum pos_score: ",
                                       min_pos_score)) +
    ggplot2::xlab("Group Number") + ggplot2::ylab("Group Number")


  #saving graph to output_dir
  #making filename for output
  output_filename <- paste0(dataset_name,
                            "_group_heatmap_min_p_",
                            min_pos_score,
                            ".",
                            output_type)

  ggplot2::ggsave(filename = output_filename,
                  path = normalizePath(output_dir),
                  plot = graph,
                  ...)

  return(graph)

}






