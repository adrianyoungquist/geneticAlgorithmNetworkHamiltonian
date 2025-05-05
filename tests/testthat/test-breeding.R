test_that("breed_two_parents creates expected number of children", {
  parents <- list(c(0, 0), c(1, 1))
  children <- breed_two_parents(parents[[1]], parents[[2]], num_points_init = 5)
  expect_type(children, "list")
  expect_true(length(children) >= 1)
})

test_that("add_noise_to_point adds noise correctly", {
  vec <- c(1, 2, 3)
  noisy_vec <- add_noise_to_point(vec, noise_var = 0.1)
  expect_equal(length(vec), length(noisy_vec))
  expect_false(all(vec == noisy_vec))  # Should almost never be exactly equal
})

test_that("breed_generation produces non-empty list", {
  parents <- list(c(0, 0), c(5, 5))
  kids <- breed_generation(parents, num_lin_kids = 5, num_mutants = 2)
  log_message(level = "INFO", "kids:")
  log_message(level = "INFO", kids)
  expect_type(kids, "list")
  expect_equal(length(kids), 12)
})
