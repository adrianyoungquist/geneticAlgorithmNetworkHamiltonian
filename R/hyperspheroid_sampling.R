# ------------------------------------------------------------------------
# hyperspheroid_sampling.R
#
# Purpose: Functions to generate random samples around a center point inside a scaled hyperspheroid.
#
# Author: Gianmarc Grazioli and contributors
# License: GPL-3
#
# This file is part of the fibga R package.
# See GPL-3 License at <https://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------

#' @title Sample Points from a Hyperspheroid
#' @description Generates random points inside a distorted hypersphere to initialize parameter sets for exploration.
#' @details
#' Sampling points from a hyperspheroid allows exploration around a center parameter set,
#' while preserving anisotropy based on the center's dimension scaling.
#' @name sample_hyperspheroid
NULL

#' Normalize a Vector
#'
#' @description
#' Normalize a numeric vector to have unit L2 norm.
#'
#' @param vec A numeric vector.
#' @return Normalized vector.
#' @examples
#' normalize(c(3, 4))
#' @export
normalize <- function(vec) {
  vec / sqrt(sum(vec^2))
}

#' Generate Random Points in a Hypersphere
#'
#' @description
#' Create random points inside a hypersphere of given radius and center.
#'
#' @param d Number of dimensions.
#' @param n Number of points.
#' @param radius Radius of the hypersphere.
#' @param center Center coordinates.
#' @return A matrix where each column is a sampled point.
#' @examples
#' get_hypersphere(3, 10, 1, c(0, 0, 0))
#' @importFrom stats rnorm
#' @export
get_hypersphere <- function(d, n, radius, center) {
  coords <- replicate(n, {
    raw_pts <- stats::rnorm(d + 1)
    radius * normalize(raw_pts)[1:d]
  }, simplify = "matrix")
  coords <- t(coords)
  sweep(coords, 2, center, FUN = "+")
}

#' Sample a Hyperspheroid
#'
#' @description
#' Samples random points from a distorted hypersphere, scaled based on the center's magnitude per dimension.
#'
#' @param d Number of dimensions.
#' @param n Number of points.
#' @param radius Sphere radius.
#' @param center Center point.
#' @param scale_fac Scaling factor applied to each dimension.
#' @param fix_edge Logical; fix first coordinate to a set value.
#' @param fix_edge_value Value to fix first coordinate to if `fix_edge = TRUE`.
#' @return A matrix of sampled points.
#' @examples
#' set.seed(42)
#' points <- sample_hyperspheroid(3, 10, radius = 2, center = c(5, 10, 15))
#' plot(points[1, ], points[2, ], main = "Hyperspheroid Sample (First 2 dimensions)")
#' @export
sample_hyperspheroid <- function(d, n, radius, center,
                                 scale_fac = 0.15, fix_edge = FALSE, fix_edge_value = 100) {
  sphere <- get_hypersphere(d, n, radius, center)
  sphere <- sweep(sphere, 2, center, `-`)

  scaling_factors <- center * scale_fac
  scaling_factors[scaling_factors == 0] <- 1e-8  # Prevent divide-by-zero when scaling

  sphere <- sweep(sphere, 2, apply(sphere, 2, max), FUN = "/")
  sphere <- sweep(sphere, 2, scaling_factors, FUN = "*")
  sphere <- sweep(sphere, 2, center, `+`)

  if (fix_edge) {
    sphere <- cbind(rep(fix_edge_value, nrow(sphere)), sphere[, -1, drop = FALSE])
  }

  rbind(sphere, center)  # Always include center point at the end
}
