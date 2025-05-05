# ------------------------------------------------------------------------
# full_evolution_script.R
#
# Purpose: # This main function carries out automated discovery of network Hamiltonian models
#  that can self-assemble into the 2-ribbon amyloid fibril topological structure.
#  The code implements the genetic algorithm described in the article:
# 
#  "A Genetic Algorithm for Automated Parameterization of Network Hamiltonian 
#                     Models of Amyloid Fibril Formation"
# by: Gianmarc Grazioli, Andy Tao, Inika Bhatia, and Patrick Regan
#
# currently under review at Digital Discovery - Royal Society of Chemistry
# preprint: https://chemrxiv.org/engage/chemrxiv/article-details/64a3d048ba3e99daef7c671a
#
# The example used in this code is very similar to the code used to generate
# the figure in the aforementioned article showing the evolution of different
# generations of models as the genetic algorithm converges on a region of parameter
# space that produces maximal fibril yield for 2-ribbon type amyloid fibril structures.
# To give an idea of run time, this main function runs in under 10 minutes on a
# 2019 MacBook Pro with a 2.4 GHz 8-Core Intel Core i9 processor and 64 GB of RAM.
#
# Author: Gianmarc Grazioli and contributors
# License: GPL-3
#
# This file is part of the fibga R package.
# See GPL-3 License at <https://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------

#' @title Run Full Evolution Pipeline
#' @description Executes the full pipeline: sample initial parameter sets, simulate, evolve generations, and evaluate final fibril yield.
#' @details
#' This pipeline implements the genetic algorithm described in:
#' "A Genetic Algorithm for Automated Parameterization of Network Hamiltonian Models of Amyloid Fibril Formation" (Gianmarc Grazioli et al., 2023).
#' It samples initial points around a known center, selects fibril-forming models, and evolves them.
#' @name full_evolution_script
NULL

#' Run the Full Evolution Pipeline
#'
#' @param init_param Initial parameter vector to sample around.
#' @param spheroid_pt_count Number of initial points to sample.
#' @param radius Radius of the initial hyperspheroid.
#' @param node_count Number of nodes in networks.
#' @param suff_stats ERGM sufficient statistics formula (string).
#' @param fib_type Fibril topology to optimize.
#' @param num_lin_kids (default 3)
#' @param num_mutants_init (default 2)
#' @param num_initial_survivors how many values to have in gen_zero (default 8)
#' @param num_generations Number of generations to evolve.
#' @param burn_time MCMC burn-in time for simulations during genetic algorithm
#' @param spheroid_sim_burn_time burn-in time for initial spheroid simulations
#' @param final_sim_burn_time burn-in time for reps of final simulation
#' @param np Number of parallel cores.
#' @param scale_fac (default 0.15)
#' @param reps (default 16)
#' @param noise_var_init (default 0.10)
#' @param top_frac (default 0.25)
#' @param max_survivors (default 4)
#' @param fix_edge (default TRUE)
#' @param fix_edge_value (default init_param\[1\])
#' @param init_yield_cut (default 0.0)
#' @param smart_var (default TRUE)
#' @param smart_pts (default TRUE)
#' @param use_line_density (default TRUE)
#' @param min_line_density (default 3.3)
#' @param print_params (default FALSE)
#' @param pair_max (default 1000)
#' @param child_max (default 20)
#' @param final_reps (default 32)
#' @param simulation_function Function for simulation (default is ERGM simulation)
#' @param sim_args A list with extra arguments for the simulation
#' @param spheroid_sim_args Arguments for initial spheroid simulationns
#' @param final_sim_args Arguments for final simulations
#' @return A list containing evolved generations and final simulation results.
#' @examples
#' \dontrun{
#' run_full_evolution(init_param = c(100, -30, 10, 25),
#'                    spheroid_pt_count = 10,
#'                    node_count = 48,
#'                    suff_stats = "edges+kstar(2)+nsp(1:2)",
#'                    generations = 2)
#' }
#' @seealso \code{\link{evolve}}: Description of parameters
#' @importFrom graphics plot points
#' @importFrom stats as.formula
#' @importFrom utils head
#' @export
run_full_evolution <- function(init_param, spheroid_pt_count = 50, radius = 2, node_count = 48,
                               suff_stats = "edges+kstar(2)+nsp(1:2)", fib_type = "2-ribbon",
                               num_lin_kids = 3, num_mutants_init = 2, num_initial_survivors = 8,
                               num_generations = 3, scale_fac = 0.15,
                               top_frac = 0.25, noise_var_init = 0.10, reps = 1,
                               burn_time = 2^20, spheroid_sim_burn_time = 2^21 - 1, final_sim_burn_time = 2^21 -1, np = 16,
                               max_survivors = 4, fix_edge = TRUE, fix_edge_value = init_param[1],
                               init_yield_cut = 0.0, smart_var = TRUE, smart_pts = TRUE,
                               use_line_density = TRUE, min_line_density = 3.3,
                               print_params = FALSE, pair_max = 1000, child_max = 20, final_reps = 32,
                               simulation_function = default_ergm_simulation, sim_args = list(burn_time = burn_time, node_count=node_count),
                               spheroid_sim_args = NULL, final_sim_args = NULL) {
  if (is.null(spheroid_sim_args)) {
    spheroid_sim_args <- sim_args
  }
  if (is.null(final_sim_args)) {
    final_sim_args <- sim_args
  }

  begin <- proc.time()[[3]]

  # Sample initial points
  gen_zero_matrix <- sample_hyperspheroid(length(init_param), spheroid_pt_count, radius, init_param,
                                          scale_fac = scale_fac, fix_edge = TRUE, fix_edge_value = init_param[1])
  gen_zero_params <- split(gen_zero_matrix, seq_len(nrow(gen_zero_matrix)))

  log_message(level = "INFO", "Simulating initial hyperspheroid points...\n")
  spheroid_sims <- simulate_parallel_with_yield(gen_zero_params, suff_stats, fib_type, node_count, np = np, burn_time = spheroid_sim_burn_time,
                                                simulation_function = simulation_function, sim_args = spheroid_sim_args)
  log_message(level = "TRACE", "Spheroid sims: ", spheroid_sims)
  # Select initial survivors
  gen_zero <- spheroid_sims$params[spheroid_sims$yields > 0]
  gen_zero <- utils::head(gen_zero, num_initial_survivors)

  log_message(level = "TRACE", "gen zero:", gen_zero)

  if (length(gen_zero) == 0) {
    stop("No valid starting points found from hyperspheroid sampling.")
  }

  # Evolve generations
  log_message(level = "INFO", "Starting evolution...\n")
  evolved_generations <- evolve(gen_zero, num_generations, num_lin_kids, num_mutants_init, fib_type, suff_stats, node_count,
                                reps = reps, noise_var_init = noise_var_init, top_frac = top_frac,
                                max_survivors = max_survivors, fix_edge = fix_edge, fix_edge_value = fix_edge_value, init_yield_cut = init_yield_cut,
                                smart_var = smart_var, smart_pts = smart_pts, use_line_density = use_line_density,
                                min_line_density = min_line_density, print_params = print_params,
                                burn_time = burn_time, np = np, pair_max = pair_max, child_max = child_max,
                                simulation_function = simulation_function, sim_args = sim_args)
  log_message(level = "DEBUG", "Evolved generations: ", evolved_generations)

  log_message(level = "INFO", "Running replications of best model")
  # Final evaluation
  reps_of_best <- list(
    replicate(final_reps, init_param, simplify = FALSE),
    replicate(final_reps, spheroid_sims$params[[which.max(spheroid_sims$yields)]], simplify = FALSE)
  )
  for (gen in seq_along(evolved_generations)) {
    reps_of_best[[2 + gen]] <- replicate(final_reps, evolved_generations[[gen]]$params[[1]], simplify = FALSE)
  }

  final_sims <- lapply(reps_of_best, function(rep_set) {
    simulate_parallel_with_yield(rep_set, suff_stats, fib_type, node_count, np = np, burn_time = final_sim_burn_time,
                                 simulation_function = simulation_function, sim_args = final_sim_args)
  })

  mean_fibril_fractions <- sapply(final_sims, function(fs) mean(fs$yields))

  # Plot mean fibril yields across evolution
 # plot(mean_fibril_fractions, type = 'l', main = "Mean Fibril Yields Across Evolution", ylab = "Mean Yield")
 # points(mean_fibril_fractions)

  end <- proc.time()[[3]]
  time_sec <- end - begin
  log_message(level = "INFO", sprintf("Total time: %d minutes and %d seconds.\n", floor(time_sec / 60), round(time_sec %% 60)))

  list(
    generations = evolved_generations,
    final_sims = final_sims,
    mean_fibril_fractions = mean_fibril_fractions
  )
}
