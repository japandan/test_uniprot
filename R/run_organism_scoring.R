#' Preprocess, filter, and calculate organism scores
#'
#' @param filt_counts_path path to filtered counts file produced by
#'   \code{\link{filter_counts}}
#' @param fasta_for_filter path to fasta file to be used for filtering. All tags
#'   that are found within this file are removed before the organisms are
#'   scored. If NULL then no tags are removed.
#' @param organism_metadata_df dataframe containing metadata for every organism
#'   in the database, must include a column for taxonomic groups as calculated
#'   by \code{\link[CandidateSearchDatabase]{leader_similarity_cluster}}
#' @param scoring_overlap_mat overlap matrix to use for scoring
#' @param min_total_tags Minimum number of total tags in the filt_counts_path
#'   file required to proceed with scoring. Defaults to 10.
#' @param min_total_strong_tags Minimum total number of strong tags required to
#'   continue scoring after strong tags are calculated. Defaults to 10.
#' @param output_dir directory to write the selected counts and organism scores.
#'   Defaults to the working directory.
#' @param export_matrix If TRUE then the reduced matrix will be written as a csv
#'   in the output_dir directory. Defaults to FALSE
#' @param filt_count_suffix Suffix to remove from the `filt_counts_path`
#'   filename such that the dataset name remains. The dataset name is used to
#'   name the result files. It defaults to "_filt_sample_counts.csv" which is
#'   the default for the output of \code{\link{filter_counts}}
#' @inheritParams calculate_taxonomic_scores
#' @inheritParams assign_tags
#' @inheritParams preprocess_counts
#'
#' @return If any organisms passed filtering the function will return the path
#'   to the result csv. It will return NA if no organisms passed filtering.
#'
#' @author Sarah C. Jenson \email{sarah.jenson@@pnnl.gov}
#'
#' @export
#'
#'@importFrom assertthat assert_that has_extension
#'@importFrom dplyr select mutate filter left_join inner_join
#'@importFrom readr read_csv write_csv
#'@importFrom tibble column_to_rownames rownames_to_column
#'@importFrom magrittr %>%
#'@importFrom rlang .data
run_organism_scoring <- function(filt_counts_path,
                                fasta_for_filter,
                                organism_metadata_df,
                                scoring_overlap_mat,
                                Pthrsh = 0.4,
                                min_total_tags = 10,
                                min_total_strong_tags = 10,
                                min_organism_strong_tags = 1,
                                min_group_hits = 2,
                                output_dir = ".",
                                filt_count_suffix = "_filt_sample_counts.csv",
                                export_matrix = FALSE)
{

  #checking inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  #checking that files exist
  assertthat::assert_that(file.exists(filt_counts_path),
                          msg = "filt_counts_path not found.")
  if(!is.null(fasta_for_filter)){
    assertthat::assert_that(file.exists(fasta_for_filter),
                            msg = "fasta_for_filter not found.")
    #checking extension of fasta_for_filter
    assertthat::has_extension(fasta_for_filter, "fasta")
  }



  assertthat::assert_that(dir.exists(output_dir), msg = "output_dir not found.")


  #verifying all organisms in scoring_overlap_mat exist in organism_metadata_df
  assertthat::assert_that(setequal(colnames(scoring_overlap_mat),
                                   organism_metadata_df$kegg_id),
      msg = "column names of scoring_overlap_mat must contain the same set of kegg identifiers as organism_metadata_df$kegg_id")

  assertthat::assert_that(setequal(rownames(scoring_overlap_mat),
                                   organism_metadata_df$kegg_id),
                          msg = "row names of scoring_overlap_mat must contain the same set of kegg identifiers as organism_metadata_df$kegg_id")


  # Defining tag counts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  #defining the dataset name by removing the filt_count_suffix which has a
  #default set to match the output of filter_counts()
  dataset_name <- gsub(filt_count_suffix, "", basename(filt_counts_path))

  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print(paste0("Starting Dataset: ", dataset_name))




  #reading filtered sample counts
  sample_counts <- readr::read_csv(filt_counts_path,
                                  col_types = readr::cols(
                                     .default = readr::col_double(),
                                     peptide_tag = readr::col_character()
                                   )) %>%
                   tibble::column_to_rownames("peptide_tag")

  print(paste0("Starting with ", nrow(sample_counts), " tags."))

  #testing for enough tags
  if(nrow(sample_counts) < min_total_tags)
  {
    print(paste0("To few tags to pass min_total_tags of ",
                 min_total_tags))
    return(NULL)
  }



  #filtering out tags contained in fasta_for_filter if it is not NA

  if (!is.null(fasta_for_filter))
  {
    sample_counts <- sample_counts %>%
      filter_by_fasta(fasta_file = fasta_for_filter)

    if(is.null(sample_counts) | nrow(sample_counts) < min_total_tags) {
      print(paste0("To few tags to pass min_total_tags of ",
                   min_total_tags))
      return(NULL)
    }
  }


  #calling preprocess_counts to filter for strong tags at the taxonomic group level
  selected_counts <- sample_counts %>%
    preprocess_counts(organism_metadata_df,
                      min_organism_strong_tags = min_organism_strong_tags)

  #checking that enough tags remain
  if(any(is.null(selected_counts),
         (nrow(selected_counts) < min_total_strong_tags),
         (ncol(selected_counts) == 0))){
    print(paste0("To few strong tags or organisms to pass minimum_total_strong_tags of ",
                min_total_strong_tags,
                " and min_organism_strong_tags of ",
                min_organism_strong_tags))
    return(NULL)
  } else
  {
    #exporting selected counts
    selected_counts_path <- file.path(
      normalizePath(output_dir, winslash = "/"),
      paste0(dataset_name, "_selected_counts.csv"),
      fsep = "/")

    selected_counts %>% tibble::rownames_to_column(var = "tag") %>%
      readr::write_csv(file = selected_counts_path)
  }

  # Calculating taxonomic scores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  #making export matrix path if export_matrix == TRUE
  if(export_matrix == TRUE){
    export_matrix_path <- file.path(
      normalizePath(output_dir, winslash = "/"),
      paste0(dataset_name, "_reduced_matrix.csv"),
      fsep = "/")
  } else {
    export_matrix_path <- NULL
  }

  #making sure that all KEGG Ids match up correctly between data structures
  selected_organisms <- colnames(selected_counts)

  overlap_mat_for_scoring <- scoring_overlap_mat[selected_organisms, selected_organisms]

  subset_organism_metadata <- tibble::tibble(kegg_id = selected_organisms) %>%
    dplyr::left_join(
      organism_metadata_df,
      by="kegg_id"
    )

  #calling function to calculate taxonomic scores
  print(paste0("run_organism_scoring is saving scores to: ",export_matrix_path))
  organism_scores <- calculate_taxonomic_scores(
    count_mat = as.matrix(selected_counts),
    overlap_mat = as.matrix(overlap_mat_for_scoring),
    organism_df = subset_organism_metadata,
    Pthrsh = Pthrsh,
    min_group_hits = min_group_hits,
    export_matrix_path = export_matrix_path
  )


  #writing results to csv if any organisms passed filtering
  if (nrow(organism_scores) > 0) {
    #exporting organism scores
    org_score_path <- file.path(
      normalizePath(output_dir, winslash = "/"),
      paste0(dataset_name, "_organism_scores.csv"),
      fsep = "/"
    )

    readr::write_csv(organism_scores, org_score_path)

    print(paste0("Dataset: ", dataset_name, " complete."))

    return(org_score_path)
  } else {
    print(paste0("Dataset: ", dataset_name, " had no groups pass filtering."))
    return(NULL)
  }

}
