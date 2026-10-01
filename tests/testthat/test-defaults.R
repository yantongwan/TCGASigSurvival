test_that("default Treg signature is available", {
  sig <- default_treg_signature()
  expect_true("FOXP3" %in% sig)
  expect_true("IL2RA" %in% sig)
  expect_gte(length(sig), 6)
})
