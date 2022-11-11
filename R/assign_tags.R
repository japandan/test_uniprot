#' Produces tag assignments from the total peptide counts and the cosine similarity matrix
#'
#' @param P cosine similarity matrix, has same dimensions as count_mat
#' @param pepcount total peptide counts, `rowSums(count_mat)` in
#'   `calculate_taxonomic_scores()`. Should be a vector of the same length as
#'   number of tags / number of rows in `count_mat`
#' @param Pthrsh minimum threshold for keeping the tag assignment, default 0.4.
#'
#' @author Kristin H Jarman
#'
assign_tags <- function(P, pepcount, Pthrsh = 0.4) {
  # Take total peptide counts and cosine similarity matrix and produce tag assignments

  # each value in P is the cosine similarity value between a tag and an organism
  #Npep is number of rows in count_mat and number of tags
  Npep <- length(pepcount)

  print(paste0("assign_tags() has ", length(pepcount), " input tags."))

  #make length pre-assigned to speed up loop
  imax <- vector()

  #pepcount is one element per tag, the value in each element corriponds to the
  #number of peptides the tag matches too across all organisms because it is
  #calculated using the row sums of count_df

  #ptmpvec is one element per organism

  #P is similarity mat and is the same dimensions as count mat, rows = tags,
  #columns = organisms
  for (nn in 1:Npep) {
    #for each tag, take the tag's row in similarity matrix and divide by how
    #many peptides the tag matches to across all organisms
    ptmpvec <- P[nn, ]/pepcount[[nn]]
    #if a tag matches to tons of peptides then ptmpvec will be smaller and Pmax
    #less likely to be bigger than Pthrsh
    Pmax <- max(ptmpvec)
    if (pepcount[[nn]] > 0 & Pmax > Pthrsh) {
      # if the tag has any hits anywhere and Pmax is bigger than the threshold
      # which(ptmpvec == Pmax) gives the positional index of Pmax within ptmpvec
      # this positional index also corresponds to an organism
      #if Pmax occurs in multiple spots then tags are assigned to multiple organisms
      imax <- c(imax, which(ptmpvec == Pmax))
    }
  }



  return(imax)
}
