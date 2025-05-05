
test_that("get_hypersphere returns valid sphere", {
  center <- c(5, 10, 15)
  coords <- get_hypersphere(3, 10, 2, center)
  expect_true(all(apply((t(coords) - center)^2, 2, function(x) {sum(x)}) < 4))
})

test_that("sample_hyperspheroid returns expected shape", {
  center <- c(5, 10, 15)
  result <- sample_hyperspheroid(d = 3, n = 10, radius = 2, center = center)
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 3)
})
