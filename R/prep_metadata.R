#' @title Prepare Metadata
#' @details Prepares Metadata for Testing.
#' @param x \code{\link{phyloseq-class}} object
#' @param sam.ids columns with sample ids matching rownames of otu_table(ps)
#' @param ref.group e.g. "no_ili" this will be coded as zero.
#' The ref group witin compare vars of interest
#' @param compare The column name with case-control description
#' @param confounder.vars confounders to create stratum e.g. c("gender", "age_yrs_oct2014")
#' @param var.interest variable of interest to test
#' @examples
#' # comord.vars <- c("BMI_2014","respiratory_disease_2014", "Roken_2014")
#' # med.vars <- c("ace_inhibitors_2014","arbs_2014","flu_shot_2014_2015")
#' # patho.vars <- c("influenza_2014_any","haemophilus_2014_any","coronavirus_2014")
#' # all.vars <- c(comord.vars,med.vars,patho.vars)
#' # phenotypic.data.1 <- prep_metadata_2(ps.genus.ctrl.ili,
#' # sam.ids = "seq_sam_name",
#' # ref.group="no_ili",
#' # compare= "condition_status",
#' # confounder.vars = c("gender", "age_yrs_oct2014"),
#' # var.interest=all.vars)
#'
#' @export

prep_metadata <- function(x,
                          sam.ids = NULL,
                          ref.group= NULL,
                          compare= NULL,
                          confounder.vars = NULL,
                          var.interest=NULL) {

  if(is.null(ref.group) | is.null(confounder.vars) | is.null(compare)){
    stop("Please supply main.group &  var.interest")
  }

  phenotypic.data <- ILI <- age <- NULL

  if(!is.null(var.interest)){
    var_ast_select <- c(sam.ids, compare, confounder.vars, var.interest)
  } else {
    var_ast_select <- c(sam.ids, compare, confounder.vars)
  }
  #x <-
  phenotypic.data <- microbiome::meta(x)[,var_ast_select]
  rownames(phenotypic.data) <- NULL
  names(phenotypic.data)[1] <- "ID"
  names(phenotypic.data)[2] <- "ILI"
  phenotypic.data <- base::transform(phenotypic.data,ILI=ifelse(ILI==ref.group,0,1))

  # covert all to numeric
  phenotypic.data <- phenotypic.data %>%
    mutate_if(is.character, as.numeric)
  #phenotypic.data <- phenotypic.data %>% replace(., is.na(.), -1)

  names(phenotypic.data)[4] <- "age"
  #phenotypic.data <- phenotypic.data[,-5] # remove age_group

  #range(phenotypic.data$age)
  age.breaks <- seq(60,90,5)
  phenotypic.data <- base::transform(phenotypic.data,
                                     age.group=cut(age,breaks=age.breaks,
                                                   right=TRUE,
                                                   include.lowest=TRUE))

  phenotypic.data$age.group <- as.character(phenotypic.data$age.group)
  #table(phenotypic.data$age.group)
  # Create strata of "confounders":
  confounders <- c("gender","age.group")
  name.of.stratification <- "Sex and age-group"
  #
  aux.data.set <- subset(phenotypic.data,select=confounders)
  aux.data.set <- base::transform(aux.data.set,stratum=rep(NA,nrow(aux.data.set)))
  for(i in (1:nrow(aux.data.set))){
    aux.data.set[i,]$stratum <- paste(aux.data.set[i,-length(aux.data.set)],collapse="/")
  }
  aux.data.set$stratum <- as.character(aux.data.set$stratum)
  # str(aux.data.set); head(aux.data.set)
  #
  phenotypic.data <- base::transform(phenotypic.data,stratum=aux.data.set$stratum)
  phenotypic.data$stratum <- as.character(phenotypic.data$stratum)
  #table(phenotypic.data$stratum)

  return(phenotypic.data)
}
