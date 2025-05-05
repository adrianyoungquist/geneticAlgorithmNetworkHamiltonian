# ------------------------------------------------------------------------
# evolution.R
#
# Purpose: Functions for running the genetic algorithm to evolve network models
# toward forming fibril structures.
#
# Author: Gianmarc Grazioli and contributors
# License: GPL-3
#
# This file is part of the fibga R package.
# See GPL-3 License at <https://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------

#' @title Evolution Functions for Genetic Algorithm
#' @description Functions to simulate, evaluate, and evolve parameter sets across multiple generations to optimize fibril yield in networks.
#' @name evolution
NULL

#' Select the Top Survivors
#'
#' @description
#' Filters top-performing parameter sets based on fibril yield.
#'
#' @param sims List of simulation replicates.
#' @param fibril_type Type of fibril topology (e.g., "2-ribbon").
#' @param top_frac Top fraction to keep.
#' @param max_survivors Maximum number of survivors.
#' @return A list of best parameters and their yields.
#' @examples
#' # thin_the_herd(sim_results, fibril_type = "2-ribbon")
#' @export
thin_the_herd <- function(sims, fibril_type, top_frac = 0.5, max_survivors = 20) {
  #print("thin_the_herd arguments")
  #args <- as.list(match.call())
  #print(args)

  fib_fracs_with_params <- lapply(seq_along(sims[[1]]), function(par_idx) {
    fib_frac_tot <- sum(sapply(sims, function(rep) {
      nets <- rep[[par_idx]]$sim
      fib_frac_data <- get(fibril_type, assay(nets))
      sum(fib_frac_data) / length(fib_frac_data)
    }))
    log_message("INFO", "fib_frac tot:", fib_frac_tot)
    list(fibFrac = fib_frac_tot / length(sims), params = sims[[1]][[par_idx]]$params)
  })

  top_n <- min(max(2, floor(top_frac * length(fib_fracs_with_params))), max_survivors)

  fib_fracs <- sapply(fib_fracs_with_params, function(x) x$fibFrac)
  params <- lapply(fib_fracs_with_params, function(x) x$params)

  fittest_indices <- order(fib_fracs, decreasing = TRUE)[1:top_n]

  list(params = params[fittest_indices],
       suff_stats = sims[[1]][[1]]$suff_stats,
       yields = sort(fib_fracs, decreasing = TRUE)[1:top_n])
}

#' Evolve Generations
#'
#' @description
#' Breeds, mutates, simulates, and selects survivors across multiple generations.
#'
#' @param orig_parents Initial parent parameter sets.
#' @param num_generations Number of generations.
#' @param num_lin_kids Number of interpolated children per pair.
#' @param num_mutants_init Initial number of mutants per child.
#' @param fibril_type Fibril type to optimize (e.g., "2-ribbon").
#' @param suff_stats ERGM sufficient statistics formula.
#' @param node_count Number of nodes.
#' @param top_frac Fraction of top survivors kept each generation.
#' @param noise_var_init Initial mutation noise variance.
#' @param reps Number of replicates per simulation.
#' @param burn_time MCMC burn-in time.
#' @param np Number of parallel cores.
#' @param max_survivors Maximum survivors kept each generation.
#' @param fix_edge Logical; fix first parameter value.
#' @param fix_edge_value Fixed value for the first parameter.
#' @param init_yield_cut Initial minimum yield cutoff.
#' @param smart_var Logical; adaptively increase mutation variance if stagnation occurs.
#' @param fudge Allowed drop in yield to still accept a generation.
#' @param smart_pts Logical; increase number of mutants if needed.
#' @param use_line_density Logical; control number of interpolated points by distance.
#' @param min_line_density Minimum density of interpolated points.
#' @param print_params Save parameter vectors to file.
#' @param pair_max Maximum allowed number of parent pairs.
#' @param child_max Maximum total number of offspring.
#' @param simulation_function Function for simulation; default is basic ergm simulation
#' @param sim_args Arguments passed to simulation function
#' @return List of generations, each containing survivors and yields.
#' @examples
#' # evolve(list(c(1, -1, 0.5)), 3, 3, 2, "2-ribbon", "edges+kstar(2)", 48)
#' @export
evolve <- function(orig_parents, num_generations, num_lin_kids, num_mutants_init,
                   fibril_type, suff_stats, node_count,
                   top_frac = 0.25, noise_var_init = 0.1, reps = 1,
                   burn_time = 2^20, np = 16, max_survivors = 20,
                   fix_edge = FALSE, fix_edge_value = 100,
                   init_yield_cut = 0.0, smart_var = TRUE,
                   fudge = 0.0, smart_pts = TRUE,
                   use_line_density = TRUE, min_line_density = 0.5,
                   print_params = FALSE, pair_max = 10000, child_max = 100000,
                   simulation_function = default_ergm_simulation, sim_args = list(burn_time = burn_time, node_count=node_count)
) {

  generations <- list()
  parents <- orig_parents
  yield_cut <- init_yield_cut
  noise_var <- noise_var_init
  num_mutants <- num_mutants_init

  if (print_params) write("", file = "all_params.csv")

  for (g in seq_len(num_generations)) {
    log_message("INFO", "Generation", g, "is born!\n")

    new_generation <- breed_generation(parents, num_lin_kids, num_mutants,
                                       noise_var, use_line_density, min_line_density,
                                       pair_max, child_max)

    if (fix_edge) {
      new_generation <- lapply(new_generation, function(p) {
        p[1] <- fix_edge_value
        p
      })
    }

    if (print_params) {
      param_lines <- sapply(new_generation, function(p) paste(c(p, g), collapse = ", "))
      write(param_lines, file = "all_params.csv", append = TRUE)
    }

    log_message("INFO", "Simulating", length(new_generation) * reps, "parameters...\n")
    sim_reps <- replicate(reps, simulate_parallel(new_generation, suff_stats, node_count, burn_time, np, simulation_function = simulation_function, sim_args = sim_args), simplify = FALSE)
    survivors <- thin_the_herd(sim_reps, fibril_type, top_frac, max_survivors)

    log_message("INFO", length(survivors$params), "survivors remaining\n")
    log_message("INFO", "Best yields =\n", survivors$yields)

    generations[[g]] <- survivors

    if (g > 1) {
      if (survivors$yields[1] >= generations[[g-1]]$yields[1] * (1.0 - fudge) && survivors$yields[1] >= yield_cut) {
        parents <- survivors$params
        yield_cut <- survivors$yields[1]
        noise_var <- noise_var_init
        num_mutants <- num_mutants_init
      } else {
        log_message("INFO", "Offspring did not surpass parents, rerunning parents.")
        if (smart_var) {
          noise_var <- noise_var * 1.05
          log_message("INFO", "Increasing noise variance to", noise_var)
          if (smart_pts) {
            num_mutants <- num_mutants + 1
            log_message("INFO", "Increasing number of mutants to", num_mutants)
          }
        }
      }
    } else {
      if (survivors$yields[1] >= yield_cut) {
        parents <- survivors$params
      } else {
        log_message("INFO", "First generation did not surpass initial yield cut, rerunning originals.")
      }
    }
  }

  return(generations)
}
