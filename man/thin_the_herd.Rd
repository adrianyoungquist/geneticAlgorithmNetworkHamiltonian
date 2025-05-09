% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evolution.R
\name{thin_the_herd}
\alias{thin_the_herd}
\title{Select the Top Survivors}
\usage{
thin_the_herd(sims, fibril_type, top_frac = 0.5, max_survivors = 20)
}
\arguments{
\item{sims}{List of simulation replicates.}

\item{fibril_type}{Type of fibril topology (e.g., "2-ribbon").}

\item{top_frac}{Top fraction to keep.}

\item{max_survivors}{Maximum number of survivors.}
}
\value{
A list of best parameters and their yields.
}
\description{
Filters top-performing parameter sets based on fibril yield.
}
\examples{
# thin_the_herd(sim_results, fibril_type = "2-ribbon")
}
