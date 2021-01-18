#' @title Illustrate.ONE.association.BY.STRATUM
#' @details Tests based on the whole data set blocked by stratum.
#' This version uses the asymptotic distribution
#' @param columns dataset variable to specifically test
#' @param data.for.testing dataset to use for testing
#' @param TYPES variable properties, ordinal, categorical, continous, etc.
#' @param variables.of.interest variables to test
#' @export

test.ONE.association.0 <- function(columns, data.for.testing, TYPES, variables.of.interest) {
  # 	print(columns)
  # 	columns <- c(1,193); B <- 100000; Spearman.rather.than.AD <- TRUE; TYPES[columns]
  p.value <- test <- sample.characteristics <- Sign <- NA
  stratum <- NULL
  response.by.treatment <- na.omit(subset(data.for.testing,
    select = c(variables.of.interest[columns], "stratum")
  ))
  # 	head(response.by.treatment); str(response.by.treatment)

  checks.by.stratum <- t(sapply(response.by.treatment$stratum, check.strata, response.by.treatment))
  good.strata <- response.by.treatment$stratum[checks.by.stratum > 1]

  tested <- FALSE

  if (length(good.strata) > 1) {
    response.by.treatment <- subset(response.by.treatment, is.element(stratum, good.strata))
    names(response.by.treatment) <- c("treatment", "response", "stratum")
    rownames(response.by.treatment) <- NULL

    type.1 <- TYPES[columns[1]]
    type.2 <- TYPES[columns[2]]
    if ((type.1 == "binary" & type.2 != "binary") |
      (type.1 == "categorical" & (type.2 != "categorical" & type.2 != "binary")) |
      (type.1 == "ordinal" & (type.2 != "ordinal" & type.2 != "categorical" & type.2 != "binary"))) {
      response.by.treatment <- response.by.treatment[, c(2, 1, 3)]

      type.2 <- TYPES[columns[1]]
      type.1 <- TYPES[columns[2]]
    }
    names(response.by.treatment) <- c("response", "treatment", "stratum")
    head(response.by.treatment)
    response.by.treatment$stratum <- as.factor(response.by.treatment$stratum)

    if ((is.element(type.2, c("ordinal")) &
      is.element(type.1, c("dirac.and.continuous", "continuous", "ordinal"))) |
      (is.element(type.2, c("binary")) &
        is.element(type.1, c("ordinal")))) {
      independence.test <- independence_test(response ~ treatment | stratum,
        data = response.by.treatment,
        distribution = "asymptotic", teststat = "scalar"
      )
      p.value <- as.numeric(pvalue(independence.test))
      Sign <- sign(statistic(spearman_test(response ~ treatment | stratum,
        data = response.by.treatment,
        distribution = "asymptotic"
      )))
      list.of.samples <- split(response.by.treatment$response, f = response.by.treatment$treatment)
      test <- "Sum statistic"
      k <- length(list.of.samples)
      if (length(names(list.of.samples)) <= 5) {
        sample.characteristics <- NULL
        for (j in (1:k)) {
          sample.characteristics <- c(sample.characteristics, length(list.of.samples[[j]]))
        }
        sample.characteristics <- paste(sample.characteristics, collapse = "/")
      } else {
        sample.characteristics <- "too many"
      }
      tested <- TRUE
    }

    if (type.1 == "continuous" & type.2 == "continuous") {
      Spearman.test <- spearman_test(response ~ treatment | stratum,
        data = response.by.treatment,
        distribution = "asymptotic"
      )
      p.value <- as.numeric(pvalue(Spearman.test))
      Sign <- sign(statistic(Spearman.test))
      test <- "Spearman"
      sample.characteristics <- NA
      tested <- TRUE
    }

    if (is.element(type.2, c("binary", "categorical")) &
      is.element(type.1, c("dirac.and.continuous", "continuous")) & (tested == FALSE)) {
      KW.test <- kruskal_test(response ~ as.factor(treatment) | stratum,
        data = response.by.treatment,
        distribution = "asymptotic"
      )
      p.value <- as.numeric(pvalue(KW.test))
      if (type.2 == "binary") {
        Sign <- sign(cor(
          y = response.by.treatment$response,
          x = response.by.treatment$treatment
        ))
      } else {
        Sign <- 0
      }
      list.of.samples <- split(response.by.treatment$response, f = response.by.treatment$treatment)
      test <- "Kruskal-Wallis"
      k <- length(list.of.samples)
      sample.characteristics <- NULL
      for (j in (1:k)) {
        sample.characteristics <- c(sample.characteristics, length(list.of.samples[[j]]))
      }
      sample.characteristics <- paste(sample.characteristics, collapse = "/")
      tested <- TRUE
    }

    if (is.element(type.1, c("binary", "categorical", "ordinal")) &
      is.element(type.2, c("binary", "categorical", "ordinal")) & (tested == FALSE)) {
      CHM.test <- cmh_test(as.factor(response) ~ as.factor(treatment) | stratum,
        data = response.by.treatment, distribution = "asymptotic"
      )
      p.value <- as.numeric(pvalue(CHM.test))
      if ((type.1 == "binary" & type.2 != "categorical") | (type.2 == "binary" & type.1 != "categorical")) {
        Sign <- sign(cor(
          y = response.by.treatment$response,
          x = response.by.treatment$treatment
        ))
      } else {
        Sign <- 0
      }
      test <- "Cochran-Mantel-Haenszel"
      list.of.samples <- split(response.by.treatment$response, f = response.by.treatment$treatment)
      k <- length(list.of.samples)
      sample.characteristics <- NULL
      for (j in (1:k)) {
        sample.characteristics <- c(sample.characteristics, length(list.of.samples[[j]]))
      }
      sample.characteristics <- paste(sample.characteristics, collapse = "/")
      tested <- TRUE
    }
  }

  if (tested == FALSE) {
    print(paste("Did not test ", paste(columns, collapse = ","), "!", sep = ""))
  }

  return(c(
    names(data.for.testing)[columns], as.character(p.value), Sign, test,
    sample.characteristics
  ))
}
#