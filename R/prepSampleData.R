#' @title Prepare Metadata
#'
#' @name prepSampleData
#'
#' @details Prepares Metadata for Testing.
#'
#' @param x \code{\link{phyloseq-class}} object
#'
#' @param sample_ids columns with sample ids matching rownames of otu_table(ps)
#'
#' @param reference_group e.g. "no_ili" this will be coded as zero.
#'                        The ref group witin compare vars of interest.
#'
#' @param compare The column name with case-control description
#'
#' @param confounder_variables Confounders to create stratum e.g. c("gender", "age_yrs_oct2014")
#'
#' @param variable_interest variable of interest to test
#'
#' @examples
#' # comord.vars <- c("BMI_2014","respiratory_disease_2014", "Roken_2014")
#' # med.vars <- c("ace_inhibitors_2014","arbs_2014","flu_shot_2014_2015")
#' # patho.vars <- c("influenza_2014_any","haemophilus_2014_any","coronavirus_2014")
#' # all.vars <- c(comord.vars,med.vars,patho.vars)
#' # phenotypic.data.1 <- prepSampleData(ps.genus.ctrl.ili,
#' # sample_ids = "seq_sam_name",
#' # ref.group="no_ili",
#' # compare= "condition_status",
#' # confounder_variables = c("gender", "age_yrs_oct2014"),
#' # variable_interest=all.vars)
#'
#' @author Sudarshan A. Shetty
#' @importFrom microbiome meta abundances
#' @importFrom tidyr pivot_wider
#' @export

prepSampleData <- function(x,
                          sample_ids = NULL,
                          reference_group = NULL,
                          compare = NULL,
                          confounder_variables = NULL,
                          variable_interest = NULL) {

  if (is.null(reference_group) | is.null(confounder_variables) | is.null(compare)) {
    stop("Please supply main.group &  variable_interest")
  }

  phenotypic.data <- age <- NULL

  # x <- ps

  if (!is.null(variable_interest)) {
    var_ast_select <- c(sample_ids, compare, confounder_variables, variable_interest)
  } else {
    var_ast_select <- c(sample_ids, compare, confounder_variables)
  }
  # x <-
  phenotypic.data <- microbiome::meta(x)[, var_ast_select]

  rownames(phenotypic.data) <- NULL
  # names(phenotypic.data)[1] <- "ID"
  #names(phenotypic.data)[2] <- "ILI"
  phenotypic.data[, compare] <- ifelse(phenotypic.data[, compare] == reference_group, 0, 1)

  # phenotypic.data <- base::transform(phenotypic.data,ILI=ifelse(ILI==ref.group,0,1))

  # covert all to numeric
  # phenotypic.data <- phenotypic.data %>%
  # mutate_if(is.character, as.numeric)
  # phenotypic.data <- phenotypic.data %>% replace(., is.na(.), -1)

  # names(phenotypic.data)[4] <- "age"
  # phenotypic.data <- phenotypic.data[,-5] # remove age_group

  # range(phenotypic.data$age)

  .check_confounders(confounder_variables)

  age.breaks <- seq(min(phenotypic.data$age), max(phenotypic.data$age), 5)
  phenotypic.data <- base::transform(phenotypic.data,
                                     age.group = cut(age,
                                                     breaks = age.breaks,
                                                     right = TRUE,
                                                     include.lowest = TRUE))

  phenotypic.data$age.group <- as.character(phenotypic.data$age.group)

  # Create strata of "confounders":

  name.of.stratification <- confounder_variables
  aux.data.set <- subset(phenotypic.data,
                         select = confounder_variables)

  aux.data.set <- base::transform(aux.data.set,
                                  stratum = rep(NA, nrow(aux.data.set)))

  for (i in (1:nrow(aux.data.set))) {

    aux.data.set[i, ]$stratum <- paste(aux.data.set[i, -length(aux.data.set)], collapse = "/")

  }
  aux.data.set$stratum <- as.character(aux.data.set$stratum)
  # str(aux.data.set); head(aux.data.set)
  #
  phenotypic.data <- base::transform(phenotypic.data, stratum = aux.data.set$stratum)
  phenotypic.data$stratum <- as.character(phenotypic.data$stratum)
  # table(phenotypic.data$stratum)

  return(phenotypic.data)
}



.check_confounders <- function(confounder) {
  if (isFALSE(any(confounder %in% "age"))) {
    stop("Provide confounder_variables must have age as one of the confounders")
  }
}
