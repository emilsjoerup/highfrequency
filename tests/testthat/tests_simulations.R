
library(highfrequency)
library(testthat)
library(data.table)
test_that("Heston", {
  set.seed(1)
  volatilityModel <- list(modelType = "Heston", sigma = matrix(c(0.2, 0.04, 0.03,
                                                                   0.04, 0.4, -0.4,
                                                                   0.03, -0.4, 0.8), ncol = 3, nrow = 3))
  driftModel <- list(modelType = "Vasicek")
  spec <- hfSimSpec(volatilityModel = volatilityModel, driftModel = driftModel, nObs = 100, nDays = 2, nSeries = 3)
  spec
  sim <- hfSim(spec)
  sim
  cols <- colnames(sim$prices)[-1]
  x <- copy(sim$prices)[,(cols) := lapply(.SD, exp), .SDcols = cols][]
  
  rc <- rCov(x, makeReturns = TRUE)
})


test_that("LiLinton", {
  set.seed(1)
  volatilityModel <- list(modelType = "LiLinton", sigma = c(0.2))
  nObs <- 391
  nDays <- 6
  nSeries <- 1
  spec <- hfSimSpec(volatilityModel = volatilityModel, nObs = nObs, nDays = nDays, nSeries = nSeries)
  sim <- hfSim(spec)
  cols <- colnames(sim$prices)[-1]
  x <- copy(sim$prices)[,(cols) := lapply(.SD, exp), .SDcols = cols][]
  rv <- rCov(x, makeReturns = TRUE)
  
  expect_equal(mean(rv[, RV]) , 0.199388218)
  expect_equal(prod(dim(x) - c(0,1)), prod(c(nObs, nDays, nSeries)))
})

