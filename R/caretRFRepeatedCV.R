#' Random Forest with Repeated Cross Validation
#'
#' @name caretRFRepeatedCV
#'
#' @details Performs a Random forest with repeated Cross-validation. The function
#'          as the name suggests uses caret under the hood. Make sure to cite the
#'          original \code{caret} package with correct version number.
#'
#' @param x \code{Phyloseq} object with \code{sample_data} and \code{otu_table}
#'
#' @param features Features to use for RF model.
#'
#' @param sample_group Variable in \code{sample_data} to
#'                     test for differences in taxa abundances
#'
#' @param folds Repeated hold-out to apply. Default is '5'
#'
#' @param repeats Number of repeats. Default is '100'
#'
#' @param seed Variable in \code{sample_data} to
#'                     test for differences in taxa abundances
#'
#' @param metric See \code{\link[caret]{train}}. Default is 'Accuracy'
#'
#' @param print_final Logical. Default is TRUE. Prints the final model output
#'
#' @param ... Additional arguments passed on to \code{\link[caret]{train}}.
#'
#' @return A list is returned of class train. See \code{\link[caret]{train}}
#'
#' @author Sudarshan A. Shetty
#'
#' @examples
#'
#' library(biomeStats)
#' library(microbiome)
#' library(biomeUtils)
#' ps <- FuentesIliGutData |>
#'              microbiome::transform("compositional") |>
#'              biomeUtils::mutateTaxaTable(FeatureID = taxa_names(FuentesIliGutData)) |>
#'              biomeUtils::filterSampleData(ILI != "L2")
#' # select features reduced for speed in example
#' features.to.use <- core_members(ps, 0.01, 25/100)
#'
#' # for example reduce folds and repeats
#' rf.fit <- caretRFRepeatedCV(ps,
#'                             sample_group = "ILI",
#'                             features = features.to.use,
#'                             folds = 3,
#'                             repeats = 5,
#'                             set.seed = 1819,
#'                             metric = "Accuracy",
#'                             print_final=TRUE)
#' print(rf.fit)
#' #Check splits for each fold
#' #table(rf.fit$pred$Resample, rf.fit$pred$obs)
#'
#' @importFrom microbiome abundances meta
#' @import caret
#'
#' @seealso \code{\link[caret]{train}}
#' @export
#'
NULL
caretRFRepeatedCV <- function(x,
                              sample_group = NULL,
                              features = NULL,
                              folds = 5,
                              repeats = 100,
                              seed = 1819,
                              metric = "Accuracy",
                              print_final = TRUE,
                              ...){

  if (!is(x, "phyloseq")){
    stop("Input must be a phyloseq object with `otu_table` and `sample_data`")
  }

  if(is.null(features)){
    stop("Provide features to test 'features' cannot be NULL")
  }
  if(is.null(sample_group)){
    stop("Provide sample_group to test 'sample_group' cannot be NULL")
  }

  prop <- as.data.frame(t(abundances(x)))
  prop <- prop[,features]
  prop$Group <- meta(x)[[sample_group]]

  set.seed(seed)

  #x folds repeat y times
  control <- caret::trainControl(method='repeatedcv',
                                 number=folds,
                                 repeats=repeats,
                                 classProbs=T,
                                 savePredictions = T)

  # Randomly variable selected is mtry
  mtry <- sqrt(ncol(prop))
  cat("Tuning grid... ")
  tunegrid <- expand.grid(.mtry=mtry)
  cat("mtry: ", tunegrid[[1]])

  rf_default <- caret::train(Group ~.,
                             data= prop,
                             method='rf',
                             metric=metric,
                             tuneGrid=tunegrid,
                             trControl=control,
                             ...)

  if(print_final && metric=="Accuracy"){
    print(rf_default$finalModel)
    cat("Accuracy:", round(rf_default$results[[2]],2),
        "(s.d:", round(rf_default$results[[4]],2), ")")
  }
  return(rf_default)
}
