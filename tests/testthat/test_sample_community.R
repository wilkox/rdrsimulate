context("sample_community")
comm <- generate_community()

test_that("sample_community runs without error", {
  expect_silent(sample_community(comm, 1000, 1000))
  expect_silent(sample_community(comm, 1, 1))
})

test_that("the sampled community has the expected number of rDNA and rRNA sequences", {
  nDNA <- sample(1:999, 1)
  nRNA <- sample(1:999, 1)
  sampled_comm <- sample_community(comm, nDNA, nRNA)
  expect_equal(sum(sampled_comm$rDNA_abund), nDNA)
  expect_equal(sum(sampled_comm$rRNA_abund), nRNA)
})

test_that("number of input and output OTUs are equal", {
  nDNA <- sample(1:999, 1)
  nRNA <- sample(1:999, 1)
  sampled_comm <- sample_community(comm, nRNA, nRNA)
  expect_equal(nrow(comm), nrow(sampled_comm))
})

