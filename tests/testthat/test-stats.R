test_that("kl_divergence", {
  # 1. Identical Distributions
  # When two distributions are identical, the KL divergence should be 0
  p1 <- c(0.25, 0.25, 0.25, 0.25) # Uniform distribution
  q1 <- c(0.25, 0.25, 0.25, 0.25) # Same uniform distribution

  # Expected: zero
  expect_equal(kl_divergence(p1, q1), 0)
  expect_equal(kl_divergence(q1, p1), 0)

  expect_equal(symmetric_kl(p1, q1), 0)

  # 2. Slightly Different Distributions
  p2 <- c(0.3, 0.2, 0.2, 0.3)
  q2 <- c(0.25, 0.25, 0.25, 0.25)

  # Expected: small but non-zero
  res <- symmetric_kl(p2, q2)
  expect_gt(res, 0)
  expect_lt(res, 1)

  # 3. Very Different Distributions - with smoothing to avoid infinite values
  p3 <- c(0.9, 0.1, 0.0001, 0.0001) # Adding small values to avoid zeros
  q3 <- c(0.0001, 0.0001, 0.1, 0.9)

  # Expected: large value
  res <- symmetric_kl(p3, q3)
  expect_gt(res, 0)
  expect_gt(log10(res), 1)

  # 4. One-hot vs Uniform
  p4 <- c(1, 0, 0, 0) # One-hot vector (first element is 1, rest are 0)
  q4 <- c(0.25, 0.25, 0.25, 0.25) # Uniform distribution

  symmetric_kl(p4, q4) # Expected: ~log(4) ≈ 1.386

  # 5. Dealing with Zeros - using Laplace smoothing
  p5 <- c(0.5, 0.5, 0.0, 0.0)
  q5 <- c(0.4, 0.3, 0.3, 0.0)

  # Add pseudocount
  p5_smooth <- p5 + 0.0001
  q5_smooth <- q5 + 0.0001
  # Renormalize
  p5_smooth <- p5_smooth / sum(p5_smooth)
  q5_smooth <- q5_smooth / sum(q5_smooth)

  symmetric_kl(p5_smooth, q5_smooth)
})
