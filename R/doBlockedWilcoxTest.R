#' Do a Blocked Wilcox Test Taxa on Groups
#'
#' @name doWilcoxTest
#'
#' @details Performs one and two sample blocked Wilcoxon tests
#'          wrapper using \code{\link[coin]{wilcox_test}}.
#'          Useful for comparing large number of taxa between groups while
#'          stratifying by confounder like sex, age, etc.
#'
#' @param x A phyloseq or data.frame. If data.frame, it must have features and
#'           sample_group
#'
#' @param features Features to test
#'
#' @param sample_group Variable in \code{sample_data} to
#'                     test for differences in taxa abundances. Must be factor.
#'
#' @param adj_method One of "holm", "hochberg", "hommel",
#'                   "bonferroni", "BH", "BY", "fdr", "none". Default="BH".
#'                   For more information see
#'                   \code{p.adjust}.
#'
#' @param block_fct Column variable to stratify/block in testing.
#'
#' @param ... Additional parameters that are passed on to
#'            \code{\link[coin]{wilcox_test}}
#'
#' @return A tibble.
#'
#'
#' @author Sudarshan A. Shetty
#'
#' @examples
#' library(biomeStats)
#' library(biomeUtils)
#' library(phyloseq)
#' library(microbiome)
#' library(dplyr)
#' data("FuentesIliGutData")
#' feat.to.test <- core_members(getProportions(FuentesIliGutData), 0.01, 25/100)
#' ps <- getProportions(FuentesIliGutData) |>
#'   filterSampleData(ILI !="L2") |>
#'   filterTaxaByNames(ids = feat.to.test, keep = TRUE)
#'
#' wilcox_results <- doBlockedWilcoxTest(x = ps,
#'                                       features = feat.to.test,
#'                                       sample_group = "ILI",
#'                                       block_fct = "sex",
#'                                       adj_method = "BH")
#' print(wilcox_results)
#'
#' @importFrom microbiome abundances meta
#' @importFrom coin wilcox_test statistic pvalue
#' @importFrom stats p.adjust
#'
#' @seealso \code{\link[rstatix]{wilcox_test}}
#' @export

doBlockedWilcoxTest <- function(x,
                                features=NULL,
                                sample_group = NULL,
                                block_fct=NULL,
                                adj_method = "BH",
                                ...){

  tst.df <- test.tib <- feature <- test_group <- block.var <- NULL
  block_variable<-block <-p.adj <- NULL

  if(is.null(features)){
    stop("Provide features to test 'features' cannot be NULL")
  }
  if(is.null(sample_group)){
    stop("Provide sample_group to test 'sample_group' cannot be NULL")
  }
  if(is.null(block_fct)){
    stop("Provide block_fct to stratify 'block_fct' cannot be NULL")
  }

  if (!is(x, "phyloseq")){
    test.tib <- x
  } else if (is(x, "phyloseq")) {
    test.tib <- t(abundances(x))

    if(!is.null(block_fct)){
      sam.tib <- meta(x)[,c(sample_group,block_fct)]
    }
    # else {
    #   sam.tib <- meta(x)[,c(sample_group)]
    # }
    test.tib <- cbind(sam.tib,test.tib[,features])
  }

  all.feat.tib <- NULL
  for (f in features) {

    tst.df <- cbind(feature  = test.tib[,f],
                    group = as.factor(test.tib[,sample_group]),
                    block.var = as.factor(test.tib[,block_fct])) |>
      as.data.frame()
    tst.df$group <- as.factor(tst.df$group)
    tst.df$block.var <- as.factor(tst.df$block.var)
    wilcox.tib <- coin::wilcox_test(feature ~ group | block.var,
                                    data = tst.df,...)

    res.tib <- tibble::tibble(feature= f,
                              test_group  = sample_group,
                              block_variable = block_fct,
                              Z = statistic(wilcox.tib),
                              pval = pvalue(wilcox.tib))

    all.feat.tib <- dplyr::bind_rows(all.feat.tib, res.tib)

  }
  all.feat.tib <- all.feat.tib |>
    dplyr::mutate(p.adj = p.adjust(all.feat.tib$pval, method = adj_method)) |>
    dplyr::mutate(sig = ifelse(p.adj < 0.05, "Y", "N"))
  return(all.feat.tib)
}

