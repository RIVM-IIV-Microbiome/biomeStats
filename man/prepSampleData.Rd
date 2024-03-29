% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepSampleData.R
\name{prepSampleData}
\alias{prepSampleData}
\title{Prepare Metadata}
\usage{
prepSampleData(
  x,
  sample_ids = NULL,
  reference_group = NULL,
  compare = NULL,
  confounder_variables = NULL,
  variable_interest = NULL
)
}
\arguments{
\item{x}{\code{\link{phyloseq-class}} object}

\item{sample_ids}{columns with sample ids matching rownames of otu_table(ps)}

\item{reference_group}{e.g. "no_ili" this will be coded as zero.
The ref group witin compare vars of interest.}

\item{compare}{The column name with case-control description}

\item{confounder_variables}{Confounders to create stratum e.g. c("gender", "age_yrs_oct2014")}

\item{variable_interest}{variable of interest to test}
}
\description{
Prepare Metadata
}
\details{
Prepares Metadata for Testing.
}
\examples{
# comord.vars <- c("BMI_2014","respiratory_disease_2014", "Roken_2014")
# med.vars <- c("ace_inhibitors_2014","arbs_2014","flu_shot_2014_2015")
# patho.vars <- c("influenza_2014_any","haemophilus_2014_any","coronavirus_2014")
# all.vars <- c(comord.vars,med.vars,patho.vars)
# phenotypic.data.1 <- prepSampleData(ps.genus.ctrl.ili,
# sample_ids = "seq_sam_name",
# ref.group="no_ili",
# compare= "condition_status",
# confounder_variables = c("gender", "age_yrs_oct2014"),
# variable_interest=all.vars)

}
\author{
Sudarshan A. Shetty
}
