#' Do Wilcox Test Taxa and Paire-wise on Groups
#'
#' @name doWilcoxTest
#'
#' @details Performs one and two sample Wilcoxon tests
#'          wrapper using \code{\link[rstatix]{wilcox_test}},
#'          \code{\link[rstatix]{adjust_pvalue}} and
#'          \code{\link[rstatix]{add_significance}} for an
#'          input \code{phyloseq} object. Useful for comparing
#'          large number of taxa between groups.
#'
#' @param x A phyloseq object
#'
#' @param sample_group Variable in \code{sample_data} to
#'                     test for differences in taxa abundances
#'
#' @param paired Logical TRUE or FALSE passed to
#'               \code{\link[rstatix]{wilcox_test}}
#'
#' @param adj_method One of "holm", "hochberg", "hommel",
#'                   "bonferroni", "BH", "BY", "fdr", "none". Default="BH".
#'                   For more information see
#'                   \code{\link[rstatix]{adjust_pvalue}}.
#'
#' @param ... Additional parameters that are passed on to
#'            \code{\link[rstatix]{wilcox_test}}
#'
#' @return A tibble.
#'
#'
#' @author Sudarshan A. Shetty
#'
#' @examples
#' library(biomeStats)
#' library(phyloseq)
#' library(microbiome)
#' library(dplyr)
#' data("FuentesIliGutData")
#' ps <- subset_samples(FuentesIliGutData, ILI %in% c("C", "L1"))
#' ps.rel <- microbiome::transform(ps, "compositional")
#' ps.rel <- core(ps.rel, detection=0.0001, prevalence=50/100)
#'
#' wilcox_results <- doWilcoxTest(ps.rel,
#'                                sample_group = "ILI",
#'                                paired = FALSE,
#'                                adj_method = "BH",
#'                                alternative = "greater") %>%
#'                   dplyr::filter(p.adj <= 0.05)
#' print(wilcox_results)
#'
#' @importFrom microbiome abundances
#' @importFrom phyloseq sample_data
#' @importFrom reshape2 melt
#' @importFrom dplyr group_by
#' @importFrom magrittr %>%
#' @importFrom rstatix wilcox_test adjust_pvalue add_significance
#'
#' @seealso \code{\link[rstatix]{wilcox_test}}
#' @export

doWilcoxTest <- function(x,
                         sample_group = NULL,
                         paired=FALSE,
                         adj_method = "BH",
                         ...){

  feat_tib <- Comparison <- variable <- wilcox_out <- NULL

  feat_tib <- microbiome::abundances(x) %>%
    as.data.frame() %>%
    t()


  feat_tib <- cbind(feat_tib,sample_data(x)[,sample_group])

  colnames(feat_tib)[ncol(feat_tib)] <- "Comparison"

  wilcox_out <- feat_tib %>%
    reshape2::melt() %>%
    dplyr::group_by(variable) %>%
    rstatix::wilcox_test(value ~ Comparison,
                         paired = paired,
                         ...) %>%
    rstatix::adjust_pvalue(method = adj_method) %>%
    rstatix::add_significance("p.adj")

  return(wilcox_out)

}

