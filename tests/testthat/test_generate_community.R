context("generate_community")

test_that("generate_community runs without error", {
  expect_silent(generate_community(100))
  expect_silent(generate_community(1))
  expect_silent(generate_community())
})

test_that("relative abundances sum to 100", {
  comm <- generate_community()
  expect_equal(sum(comm$rDNA_relabund), 100)
  expect_equal(sum(comm$rRNA_relabund), 100)
})

test_that("expected number of OTUs are generated", {
  n <- sample(1:99999, 1)
  expect_equal(nrow(generate_community(n)), n)
})
