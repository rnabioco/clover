test_that("read_bcerror works", {
  bcerr_path <- clover_example("yeast/grande.bcerr.tsv.gz")
  bcerr <- read_bcerror(bcerr_path)

  expect_equal(nrow(bcerr), 37902)
  expect_equal(ncol(bcerr), 15)
})
