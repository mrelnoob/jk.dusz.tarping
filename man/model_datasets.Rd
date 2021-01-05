% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_data_preparation.R
\name{model_datasets}
\alias{model_datasets}
\title{Datasets for Models Building}
\usage{
model_datasets(
  response.var = c("efficiency", "edges", "overlaps", "latest_condition", "fixation")
)
}
\arguments{
\item{response.var}{A character string specifying which modelling dataset should be produced (either:
efficiency", "edges", "overlaps", "latest_condition", or "fixation"):
\itemize{
\item "efficiency" will produce the dataset for the 3 efficiency evaluation variables (namely \emph{
eff_eradication}, \emph{eff_expansion} and \emph{eff_vigour});
\item "edges" will produce the dataset having for response variable the tarping operations that observed
regrowth at the edge of the tarped area;
\item "overlaps" produce the dataset having for response variable the tarping operations that observed
regrowth at at strip overlaps;
\item "latest_condition" will produce the dataset having for response variable the condition (good, bad, etc.)
of the fabric at the date of the latest observation;
\item "fixation" will produce the dataset having for response variable the tarping operations that reported
fixation problems during the operation.
}}
}
\value{
A tibble.
}
\description{
The \code{model_datasets} function creates the respective reduced datasets that should be used
to model each response variables. For instance, if \code{response.var = "fixation"}, the function will
produce a dataset containing \emph{pb_fixation} as the \strong{response} variable and a subset of
variables to be used as \strong{predictors/explanatory variables} (removing all the variables that
should not be used to model \emph{pb_fixation}'s variations).
}
\examples{
\dontrun{
eff_model <- model_datasets(response.var = "efficiency")
}
}