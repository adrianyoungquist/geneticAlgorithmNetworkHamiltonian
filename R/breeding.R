# ------------------------------------------------------------------------
# breeding.R
#
# Purpose: Functions for breeding parent vectors and generating new generations
# for the genetic algorithm in network fibril evolution.
#
# Author: Gianmarc Grazioli and contributors
# License: GPL-3
#
# This file is part of the fibga R package.
# See GPL-3 License at <https://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------

#' @title Breeding Functions for Parameter Evolution
#' @description Functions to create offspring parameter vectors by interpolation and mutation, supporting the genetic algorithm.
#' @name breeding
NULL

#' Breed Two Parents
#'
#' @description
#' Creates evenly spaced interpolated points between two parent vectors based on the distance between them.
#'
#' @param vec1 First parent vector.
#' @param vec2 Second parent vector.
#' @param num_points_init Initial number of children points to generate.
#' @param use_line_density Logical; adjust number of points based on distance.
#' @param min_line_density Minimum allowed point density along the line.
#' @return A list of offspring vectors.
#' @examples
#' parent1 <- c(0, 0)
#' parent2 <- c(1, 1)
#' children <- breed_two_parents(parent1, parent2, 5)
#' print(children)
#' @export
breed_two_parents <- function(vec1, vec2, num_points_init, use_line_density = TRUE, min_line_density = 0.5) {
  basis <- vec2 - vec1
  distance <- sqrt(sum(basis^2))

  # Adjust number of points based on the distance if desired
  num_points <- if (use_line_density) {
    line_density <- num_points_init / distance
    if (line_density < min_line_density) round(min_line_density * distance) else num_points_init
  } else {
    num_points_init
  }

  steps <- seq(0, 1, length.out = num_points)
  children <- lapply(steps, function(s) vec1 + s * basis)

  if (length(children) > 2) children <- children[-c(1, length(children))]
  return(children)
}

#' Add Noise to a Point
#'
#' @description
#' Adds independent Gaussian noise to each coordinate of a vector.
#'
#' @param vec A numeric vector.
#' @param noise_var Standard deviation of noise.
#' @return A noisy vector.
#' @examples
#' add_noise_to_point(c(1, 2, 3), noise_var = 0.1)
#' @importFrom stats rnorm
#' @export
add_noise_to_point <- function(vec, noise_var = 0.05) {
  vec + rnorm(length(vec), mean = 0, sd = noise_var)
}

#' Breed a New Generation
#'
#' @description
#' Breeds a new generation by interpolating between all parent pairs and mutating children and parents.
#'
#' @param parents List of parent vectors.
#' @param num_lin_kids Number of children between each parent pair.
#' @param num_mutants Number of mutants per child.
#' @param noise_var Variance of the Gaussian noise for mutation.
#' @param use_line_density Logical; adapt number of points to distance.
#' @param min_line_density Minimum line density.
#' @param pair_max Maximum number of parent pairs.
#' @param child_max Maximum number of total offspring.
#' @return A list of offspring vectors.
#' @examples
#' parents <- list(c(0, 0), c(1, 1), c(2, 2))
#' new_gen <- breed_generation(parents, num_lin_kids = 5, num_mutants = 2)
#' length(new_gen)
#'
#' @importFrom utils combn
#' @export
breed_generation <- function(parents, num_lin_kids, num_mutants, noise_var = 0.1,
                             use_line_density = TRUE, min_line_density = 0.5,
                             pair_max = 10000, child_max = 100000) {
  if (num_lin_kids < 3) warning("If num_lin_kids < 3, no points between will be created.")

  pairs <- utils::combn(seq_along(parents), 2, simplify = FALSE)
  log_message("INFO", "pairs:")
  log_message("INFO", pairs)
  if (length(pairs) > pair_max) {
    warning("Max pairs exceeded, sampling subset.")
    pairs <- sample(pairs, pair_max)
  }

  children <- do.call(c, lapply(pairs, function(pair) {
    par1 <- parents[[pair[1]]]
    par2 <- parents[[pair[2]]]
    dist <- sqrt(sum((par1 - par2)^2))
    if (dist * min_line_density > 1) {
      breed_two_parents(par1, par2, num_lin_kids, use_line_density, min_line_density)
    } else {
      NULL
    }
  }))

  if (length(children) < 2) {
    warning("No breeding pairs formed; returning parents.")
    return(parents)
  }

  if (num_mutants > 0) {
    mutant_children <- lapply(rep(children, each = num_mutants), add_noise_to_point, noise_var = noise_var)
    if (length(mutant_children) > child_max) mutant_children <- sample(mutant_children, child_max)
    parent_mutants <- lapply(rep(parents, each = num_mutants), add_noise_to_point, noise_var = noise_var)
    children <- c(mutant_children, parent_mutants, parents)
  }

  return(children)
}
