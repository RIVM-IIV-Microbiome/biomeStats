#' Do Kruskal-Wallis Test Taxa and Multiple Groups
#'
#' @name doKruskalWallisTest
#'
#' @details Performs Kruskal-Wallis Test one or multiple levels and
#'          is a wrapper using \code{\link[rstatix]{kruskal_test}},
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
#' @param adj_method One of "holm", "hochberg", "hommel",
#'                   "bonferroni", "BH", "BY", "fdr", "none". Default="BH".
#'                   For more information see
#'                   \code{\link[rstatix]{adjust_pvalue}}
#'
#' @param effect_size Logical. Default is TRUE. Calculates effect size for
#'                    Kruskal-Wallis test using
#'                    \code{\link[rstatix]{kruskal_effsize}}
#'
#' @param ... Additional parameters that are passed on to
#'            \code{\link[rstatix]{kruskal_effsize}}
#'
#' @author Sudarshan A. Shetty
#'
#' @return A tibble.
#'
#' @examples
#' library(biomeStats)
#' library(phyloseq)
#' library(microbiome)
#' library(dplyr)
#' data("FuentesIliGutData")
#' ps.rel <- microbiome::transform(FuentesIliGutData, "compositional")
#' ps.rel <- core(ps.rel, detection=0.001, prevalence=50/100)
#'
#' kw_results <- doKruskalWallisTest(ps.rel,
#'                                   sample_group = "ILI",
#'                                   adj_method = "BH",
#'                                   effect_size = TRUE)
#'
#' print(kw_results)
#' @seealso \code{\link[rstatix]{kruskal_test}}
#'          \code{\link[rstatix]{kruskal_effsize}}
#' @importFrom microbiome abundances
#' @importFrom phyloseq sample_data
#' @importFrom reshape2 melt
#' @importFrom dplyr group_by rename
#' @importFrom magrittr %>%
#' @importFrom rstatix kruskal_test adjust_pvalue add_significance kruskal_effsize
#' @export
#'

doKruskalWallisTest <- function(x,
                         sample_group = NULL,
                         adj_method = "BH",
                         effect_size = TRUE,
                         ...){

  feat_tib <- Comparison <- variable <- wilcox_out <- NULL

  feat_tib <- microbiome::abundances(x) %>%
    as.data.frame() %>%
    t()


  feat_tib <- cbind(feat_tib,sample_data(x)[,sample_group])

  colnames(feat_tib)[ncol(feat_tib)] <- "Comparison"

  feat_melt <- feat_tib %>%
    reshape2::melt()

  kw_out <- feat_melt %>%
    dplyr::group_by(variable) %>%
    rstatix::kruskal_test(value ~ Comparison) %>%
    rstatix::adjust_pvalue(method = adj_method) %>%
    rstatix::add_significance("p.adj") %>%
    dplyr::mutate_if(is.factor, as.character)

  if(effect_size){
    ks_efsize <- feat_melt %>%
      dplyr::group_by(variable) %>%
      rstatix::kruskal_effsize(value ~ Comparison, ...) %>%
      dplyr::rename(effsize_method = "method") %>%
      dplyr::mutate_if(is.factor, as.character)


    kw_out <- kw_out %>%
      left_join (ks_efsize)

    return(kw_out)

  } else {

    return(kw_out)

    }



}


