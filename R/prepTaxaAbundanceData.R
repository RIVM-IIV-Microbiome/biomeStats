#' Prepares Taxa Abundance Data for Association analysis
#'
#' @name prepTaxaAbundanceData
#'
#' @details Prepares the \code{\link{phyloseq-class}} object to
#'          taxa abundance table for check associations
#'
#' @param x \code{\link{phyloseq-class}} object
#'
#' @param detection_threshold Remove taxa with less than this abundance in total dataset
#'
#' @param prevalence_threshold Remove taxa with less than this prevalence in total dataset
#'
#' @return A tibble
#' @examples
#' library(biomeStats)
#' data("FuentesIliGutData")
#' ps <- FuentesIliGutData
#' gen_tab <- prepTaxaAbundanceData(ps,
#'                         detection_threshold = 0.01,
#'                         prevalence_threshold = 10/100)
#' head(gen_tab)
#'
#' @author Sudarshan A. Shetty
#'
#' @references
#' Ferreira JA, Fuentes S. (2020). Some comments on certain statistical
#' aspects of the study of the microbiome.
#' \emph{Briefings in bioinformatics} 21(4), pp.1487-1494.
#' \url{https://doi.org/10.1093/bib/bbz077}
#'
#' @importFrom microbiome core_members abundances
#' @importFrom tidyr pivot_wider
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr as_tibble group_by mutate ungroup select distinct left_join
#' @importFrom reshape2 melt
#' @importFrom stats sd
#' @export
prepTaxaAbundanceData <- function(x,
                                   detection_threshold = 0.001,
                                   prevalence_threshold = 5 / 100) {

  taxa_abund <- taxa_abund_tib <- sample_ldf_sub <- NULL
  Shannon <- value.reabund <- value <- variable <- SD <- NULL
  prop.positive <- CV <- entropy <- average <- NULL
  taxa_wide <- sample_ids <- core_taxa <- taxa_ldf <- taxa_wide_core <- NULL
  ID <- sum.i <- NULL

  taxa_abund <- t(as.data.frame(microbiome::abundances(x)))
  taxa_abund_tib <- dplyr::as_tibble(taxa_abund, rownames = "sample_ids")

  # names(asv_tab)[1] <- "ID"
  # get column of bacterial abundances
  # message("calculating ....")
  Shannon <- function(p) {
    entropy <- ifelse(p > 0, -p * log(p), 0)
    return(entropy)
  }

  taxa_ldf <- reshape2::melt(taxa_abund_tib)
  taxa_ldf <- taxa_ldf %>%
    dplyr::group_by(sample_ids) %>%
    dplyr::mutate(
      average = mean(value, na.rm = TRUE),
      SD = sd(value, na.rm = TRUE),
      prop.positive = mean(value > 0, na.rm = TRUE),
      sum.i = sum(value, na.rm = TRUE),
      value.reabund = value / sum.i,
      entropy = sum(Shannon(value.reabund)),
      CV = SD / average
    ) %>%
    dplyr::ungroup()

  # head(subset(asv_tab2, ID=="052_056_S56"))

  # other <- c("average","SD","prop.positive","entropy","CV")
  sample_ldf_sub <- taxa_ldf %>%
    dplyr::select(sample_ids, average, SD, prop.positive, entropy, CV) %>%
    dplyr::distinct(sample_ids, .keep_all = T)

  taxa_wide <- taxa_ldf %>%
    tidyr::pivot_wider(
      id_cols = sample_ids,
      names_from = variable,
      values_from = value.reabund
    ) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("sample_ids")

  # Filter our data
  core_taxa <- core_members(t(taxa_wide),
               detection = detection_threshold,
               prevalence = prevalence_threshold
  )
  taxa_wide_core <- taxa_wide[, core_taxa]

  # Join the selected taxa with properties of sample.
  taxa_wide_core <- taxa_wide_core %>%
    tibble::rownames_to_column("sample_ids") %>%
    dplyr::left_join(sample_ldf_sub, by = "sample_ids")


  # message("done ....")

  return(taxa_wide_core)
}
