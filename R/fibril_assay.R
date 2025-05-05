# ------------------------------------------------------------------------
# fibril_assay.R
#
# Purpose: Classify nodes in networks into different fibril topologies based on ORCA orbit counts.
#
# Author: Gianmarc Grazioli and contributors
# License: GPL-3
#
# This file is part of the fibga R package.
# See GPL-3 License at <https://www.gnu.org/licenses/>.
#
# Note:
# This fibril assay method biases the genetic algorithm toward forming longer fibril segments
# by only counting interior nodes. It may slightly undercount nodes on fibril ends.
# Based on ORCA orbit counts; see:
# - https://www.jstatsoft.org/article/view/v071i10/0
# - ergm.graphlets package: http://CRAN.R-project.org/package=ergm.graphlets
# ------------------------------------------------------------------------

#' @title Fibril Assay Functions
#' @description Classify nodes of networks into fibril topologies using ORCA orbit counts.
#' @details
#' This function calculates the yield of different fibril topologies.
#' It counts only nodes that are part of the interior of a fibrillar section of a graph,
#' leading to slight undercounting at fibril ends but encouraging longer fibrils.
#'
#' For background on ORCA orbit counting, see:
#' \url{https://www.jstatsoft.org/article/view/v071i10/0}.
#' Also see ergm.graphlets package:
#' \url{http://CRAN.R-project.org/package=ergm.graphlets}.
#' @name fibril_assay
NULL

#' Classify Fibril Structures
#'
#' @description
#' Takes a network and classifies each node into various fibril topology types based on graphlet orbit counts.
#'
#' @param net A network object.
#' @return A data frame where each row is a node and each column indicates membership in a fibril topology.
#' @examples
#' set.seed(0)
#' net <- network::network.initialize(24, directed = FALSE)
#' net <- default_ergm_simulation(c(100, -50), "edges+kstar(2)",
#'              list(burn_time = 2^20, node_count = 24))
#' assay(net)
#' @importFrom network as.edgelist
#' @importFrom orca count5
#' @export
assay <- function(net) {
  edgelist <- network::as.edgelist(net)
  od <- as.data.frame(orca::count5(edgelist))

  # Classify node into different fibril types based on ORCA orbit counts
  classify_node <- function(row) {
    c(
      "1-ribbon" = as.integer(row["O15"] == 2 & row["O16"] == 2 & row["O8"] == 0 & row["O7"] == 0),
      "2-ribbon" = as.integer(row["O15"] == 12 & row["O16"] == 12 & row["O37"] == 4),
      "1,2 2-ribbon" = as.integer((row["O59"] == 1 & row["O60"] == 1 & row["O13"] == 1) |
                                    (row["O45"] == 2 & row["O46"] == 2 & row["O47"] == 2 &
                                      row["O48"] == 4 & row["O61"] == 1 & row["O60"] == 2 & row["O59"] == 2)),
      "double 1,2 2-ribbon" = as.integer((row["O67"] == 4 & row["O66"] == 4 & row["O65"] == 2) |
                                           (row["O67"] == 0 & row["O66"] == 2 & row["O65"] == 1)),
      "3-prism" = as.integer(row["O53"] == 4 & row["O52"] == 2 & row["O51"] == 4 &
                               row["O38"] == 4 & row["O37"] == 8 & row["O36"] == 4 &
                               row["O35"] == 4 & row["O33"] == 1)
    )
  }
  membership <- t(apply(od, 1, classify_node))
  as.data.frame(membership)
}
