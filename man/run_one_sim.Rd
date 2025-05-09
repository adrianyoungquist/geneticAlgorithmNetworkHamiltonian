% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{run_one_sim}
\alias{run_one_sim}
\title{Simulate a Single Parameter Set}
\usage{
run_one_sim(
  param,
  suff_stats,
  node_count,
  burn_time = 2^20,
  simulation_function = default_ergm_simulation,
  sim_args = list(burn_time = burn_time, node_count = node_count),
  disable_parallel = TRUE
)
}
\arguments{
\item{param}{Numeric parameter vector.}

\item{suff_stats}{ERGM sufficient statistics formula (string).}

\item{node_count}{Number of nodes in the network.}

\item{burn_time}{MCMC burn-in time.}

\item{simulation_function}{Function for simulation}

\item{sim_args}{Arguments for simulation}

\item{disable_parallel}{(boolean) If true, set cores to 0 (useful when already instide parallelized loop)}
}
\value{
A network object.
}
\description{
Simulates a network model given a parameter vector using ERGM.
}
\examples{
run_one_sim(c(100, -30), "edges+kstar(2)", 48)
}
