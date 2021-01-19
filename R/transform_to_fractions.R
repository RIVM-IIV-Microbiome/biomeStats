#' Transform to Fractions
#'
#' @details Prepares the \code{\link{phyloseq-class}} object to
#'          taxa abundance table for testing.
#'
#' @param x \code{\link{phyloseq-class}} object
#' @param det.thres Remove taxa with less than this abundance in total dataset
#' @param prev.thres Remove taxa with less than this prevalence in total dataset
#'
#' @examples
#' library(biomestats)
#' data("ili_data")
#' ps <- ili_data
#' gen_tab <- transform_to_fractions(ps,
#'                                   det.thres= 0.0001,
#'                                   prev.thres=5/100)
#'
#' @importFrom tidyr pivot_wider
#' @importFrom tibble rownames_to_column column_to_rownames
#'
#' @return A tibble
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' Ferreira, J.A. and Fuentes, S., 2020. Some comments on certain statistical
#' aspects of the study of the microbiome.
#' \emph{Briefings in bioinformatics} 21(4), pp.1487-1494.
#' \url{https://doi.org/10.1093/bib/bbz077}
#'
#' @export

transform_to_fractions <- function(x,
                                   det.thres = 0.001,
                                   prev.thres = 5 / 100) {

  taxa_abund <- taxa_abund_tib <- sample_ldf_sub <- NULL
  Shannon <- value.reabund <- value <- variable <- SD <- NULL
  prop.positive <- CV <- entropy <- average <- NULL
  taxa_wide <- core_taxa <- taxa_ldf <- taxa_wide_core <- NULL
  ID <- sum.i <- NULL

  # x <- ps

  taxa_abund <- t(as.data.frame(microbiome::abundances(x)))
  taxa_abund_tib <- dplyr::as_tibble(taxa_abund, rownames = "ID")

  # names(asv_tab)[1] <- "ID"
  # get column of bacterial abundances
  # message("calculating ....")
  Shannon <- function(p) {
    entropy <- ifelse(p > 0, -p * log(p), 0)
    return(entropy)
  }

  taxa_ldf <- reshape2::melt(taxa_abund_tib)
  taxa_ldf <- taxa_ldf %>%
    group_by(ID) %>%
    mutate(
      average = mean(value, na.rm = TRUE),
      SD = sd(value, na.rm = TRUE),
      prop.positive = mean(value > 0, na.rm = TRUE),
      sum.i = sum(value, na.rm = TRUE),
      value.reabund = value / sum.i,
      entropy = sum(Shannon(value.reabund)),
      CV = SD / average
    ) %>%
    ungroup()

  # head(subset(asv_tab2, ID=="052_056_S56"))

  # other <- c("average","SD","prop.positive","entropy","CV")
  sample_ldf_sub <- taxa_ldf %>%
    select(ID, average, SD, prop.positive, entropy, CV) %>%
    distinct(ID, .keep_all = T)

  taxa_wide <- taxa_ldf %>%
    tidyr::pivot_wider(
      id_cols = ID,
      names_from = variable,
      values_from = value.reabund
    ) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("ID")

  # Filter our data
  core_taxa <- core_members(t(taxa_wide),
               detection = det.thres,
               prevalence = prev.thres
  )
  taxa_wide_core <- taxa_wide[, core_taxa]

  # Join the selected taxa with properties of sample.
  taxa_wide_core <- taxa_wide_core %>%
    tibble::rownames_to_column("ID") %>%
    left_join(sample_ldf_sub, by = "ID")


  # message("done ....")

  return(taxa_wide_core)
}
