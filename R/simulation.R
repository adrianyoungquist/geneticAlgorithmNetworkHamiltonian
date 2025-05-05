# ------------------------------------------------------------------------
# simulation.R
#
# Purpose: Functions to simulate networks from parameter lists and measure fibril formation.
#
# Author: Gianmarc Grazioli and contributors
# License: GPL-3
#
# This file is part of the fibga R package.
# See GPL-3 License at <https://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------

#' @title Simulate Networks from Parameter Lists
#' @description Simulates networks in parallel from a list of parameter vectors and calculates their fibril yields.
#' @name simulation
NULL


#' Simulate a Single Parameter Set
#'
#' @description
#' Simulates a network model given a parameter vector using ERGM.
#'
#' @param param Numeric parameter vector.
#' @param suff_stats ERGM sufficient statistics formula (string).
#' @param node_count Number of nodes in the network.
#' @param burn_time MCMC burn-in time.
#' @param simulation_function Function for simulation
#' @param sim_args Arguments for simulation
#' @param disable_parallel (boolean) If true, set cores to 0 (useful when already instide parallelized loop)
#' @return A network object.
#' @examples
#' run_one_sim(c(100, -30), "edges+kstar(2)", 48)
#' @importFrom ergm control.simulate.formula
#' @importFrom network network.initialize
#' @importFrom stats as.formula simulate
#' @export
run_one_sim <- function(param, suff_stats, node_count, burn_time = 2^20,
                        simulation_function = default_ergm_simulation, sim_args = list(burn_time = burn_time, node_count = node_count),
                        disable_parallel = TRUE) {
  log_message(level = "INFO", "Run one sim for ", param)

  if (disable_parallel) {
    sim_args$np <- 1
  }
  simulation_function(param, suff_stats, sim_args)
}


#' Simulate Networks from a List of Parameters
#'
#' @description
#' Given a list of parameter sets, simulate networks and evaluate fibril yields based on the chosen topology.
#'
#' @param param_list List of parameter vectors.
#' @param suff_stats ERGM sufficient statistics formula as a string.
#' @param fib_type Type of fibril topology to measure (e.g., "2-ribbon").
#' @param node_count Number of nodes in the network.
#' @param np Number of parallel cores.
#' @param burn_time MCMC burn-in time.
#' @param simulation_function Function for simulation (default default_ergm_simulation)
#' @param sim_args Arguments for simulation_function
#' @return A list containing parameters, the statistics formula, and their corresponding fibril yields.
#' @examples
#' simulate_parallel_with_yield(list(c(100, -30), c(101, -25)),
#'     "edges+kstar(2)", "2-ribbon", 48, np = 2)
#' @importFrom parallel makeCluster parLapply stopCluster
#' @export
simulate_parallel_with_yield <- function(param_list, suff_stats, fib_type, node_count,
                                         np = 16, burn_time = sum(2^(1:20)),
                                         simulation_function = default_ergm_simulation, sim_args = list(burn_time = burn_time, node_count = node_count)) {
  log_message(level = "INFO", "simulate_parallel_with_yield")

  # Helper function to calculate fibril yield from network
  get_yield <- function(net_in) {
    counts <- assay(net_in)
    sum(counts[[fib_type]]) / nrow(counts)
  }

  # Helper function to simulate one parameter set
  run_simulation <- function(param, fib, disable_parallel = TRUE) {

    sim_result <- run_one_sim(param, suff_stats, node_count, burn_time, simulation_function, sim_args = sim_args, disable_parallel = disable_parallel)
    counts <- assay(sim_result)
    yield <- sum(counts[[fib]]) / nrow(counts)
    log_message(level = "INFO", "yield: ", yield)
    list(sim_result = sim_result, yield = yield, param = param)
  }

  func_args <-  list(sim_args = sim_args, fib_type = fib_type, suff_stats = suff_stats, node_count = node_count, burn_time = burn_time, simulation_function = simulation_function)
  if (np > 1) {
    cl <- parallel::makeCluster(np)
    on.exit(parallel::stopCluster(cl))
    log_message(level = "INFO", "created cluster")
    results <- parallel::parLapply(cl, param_list, function(param, args = func_args) {
      log_message(level = "INFO", "starting run_one_sim for ", param)
      sim_result <- run_one_sim(param, args$suff_stats, args$node_count, args$burn_time, args$simulation_function, sim_args = args$sim_args, disable_parallel = TRUE)
      counts <- assay(sim_result)
      yield <- sum(counts[[args$fib_type]]) / nrow(counts)
      log_message(level = "INFO", "yield: ", yield)
      list(sim_result = sim_result, yield = yield, param = param)
    })
  } else {
    results <- lapply(param_list, function(x) {run_simulation(x, fib_type, disable_parallel = FALSE)})
  }

  list(
    params = lapply(results, `[[`, "param"),
    suff_stats = suff_stats,
    yields = sapply(results, `[[`, "yield")
  )
}


#' Simulate All Parameter Sets
#'
#' @description
#' Simulates multiple parameter sets in parallel.
#'
#' @param param_list List of parameter vectors.
#' @param suff_stats ERGM sufficient statistics formula.
#' @param node_count Number of nodes.
#' @param burn_time MCMC burn-in time.
#' @param np Number of parallel cores.
#' @param simulation_function (default default_ergm_simulation)
#' @param sim_args Arguments for simulation
#' @return List of simulation results.
#' @examples
#' simulate_parallel(list(c(1, -1), c(2, -2)), "edges+kstar(2)", 48, np = 2)
#' @importFrom parallel makeCluster parLapply stopCluster
#' @export
simulate_parallel <- function(param_list, suff_stats, node_count, burn_time = 2^20, np = 16,
                              simulation_function = default_ergm_simulation, sim_args = list(burn_time=burn_time, node_count=node_count)) {

  func_args <-  list(sim_args = sim_args, suff_stats = suff_stats, node_count = node_count, burn_time = burn_time, simulation_function = simulation_function)
  if (np > 1) {
    cl <- parallel::makeCluster(np)
    on.exit(parallel::stopCluster(cl))
    log_message(level = "INFO", "created cluster")
    results <- parallel::parLapply(cl, param_list, function(param, args = func_args) {
      log_message(level = "INFO", "starting run_one_sim for ", param)
      sim_result <- run_one_sim(param, args$suff_stats, args$node_count, args$burn_time, args$simulation_function, sim_args = args$sim_args, disable_parallel = TRUE)
      list(sim_result = sim_result, params = param, suff_stats = args$suff_stats)
    })
  } else {
    results <- lapply(param_list, function(param) {
      list(sim = run_one_sim(param, suff_stats, node_count, burn_time, simulation_function, sim_args = sim_args, disable_parallel = FALSE),
           params = param,
           suff_stats = suff_stats)
    })
  }
  return(results)
}

#' @title Default ERGM Simulation Function
#' @description A basic ERGM simulation using ergm model
#' @param theta Parameter vector (coefficients).
#' @param suff_stats ERGM statistics formula string.
#' @param sim_args List of extra simulation arguments (burn time, etc.)
#' @return A simulated network.
#' @importFrom ergm control.simulate.formula
#' @importFrom stats as.formula simulate
#' @export
default_ergm_simulation <- function(theta, suff_stats, sim_args) {
  nw <- network::network.initialize(sim_args$node_count, directed = FALSE)

  formula <- stats::as.formula(paste("nw", suff_stats, sep = "~"))
  stats::simulate(formula, coef = theta,
                 control = ergm::control.simulate.formula(MCMC.burn = sim_args$burn_time,
                                                          MCMC.interval = 1),
                 constraint = ~bd(maxout = 12),
                 output = "network",
                 basis = nw)
}



