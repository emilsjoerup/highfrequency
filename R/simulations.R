#' simulate hf prices
#' @param spec object of type spec created by \code{\link{hfSimSpec}}
#' @return an object of class \code{hfSim} which is a list containing the simulated data and any applicable additional information, e.g. underlying variance process, bursts in the drift or volatility.
#' @details The class \code{hfSim} has the following methods implemented:
#' \enumerate{
#' \item{\code{print}} prints the simulated prices and information about the spec and simulation
#' }
#' 
#' @author Emil Sjoerup
#' @export
hfSim <- function(spec){
  UseMethod('hfSim', spec)
}
#' @export
hfSim.highfrequencySimSpec <- function(spec){
  tic <- Sys.time()
  ## Check that the spec is correctly specified
  ## We replace any missing values with reasonable choices and error out on outright wrong settings, e.g. negative variance.
  spec <- hfSimSpec(volatilityModel = spec$volatilityModel, driftModel = spec$driftModel, jumpModel = spec$jumpModel, diurnalModel = spec$diurnalModel,
                    burstModel = spec$burstModel, noiseModel = spec$noiseModel, timeSettings = spec$timeSettings, nSeries = spec$nSeries, nDays = spec$nDays, 
                    nObs = spec$nObs, discretize = spec$discretize, returnType = spec$returnType)
  
  # Unpacking the spec
  includeJumps <- spec$jumpModel$includeJumps
  driftModel <- spec$driftModel
  volatilityModel <- spec$volatilityModel
  diurnalModel <- spec$diurnalModel
  burstModel <- spec$burstModel
  jumpModel <- spec$jumpModel
  noiseModel <- spec$noiseModel
  timeSettings <- spec$timeSettings
  nObs <- spec$nObs
  nSeries <- spec$nSeries
  nDays <- spec$nDays
  volatilityModel$sigma <- as.matrix(volatilityModel$sigma)
  
  prices <- NULL

  if (spec$timeSettings$sampling != "equidistant") {
    timestamps <- rep(0:(nDays-1), each = nObs) * 86400
  } else {
    timestamps <- rep(0:(nDays-1), each = nObs) * 86400 + seq(spec$timeSettings$tradingStart, spec$timeSettings$tradingEnd, length.out = nObs)
    dt <- matrix(rep(1/nObs, nDays * nObs * nSeries), ncol= nSeries)
  }

  # Here, we lower the volatility of the volatility process to account for the jump variation
  if (jumpModel$modelType != "none" && includeJumps) {
    jumpModel$jumpVolatility <- diag(volatilityModel$sigma * jumpModel$jumpComponent) # We create a variable to contain the jump volatility which will be used in the RNG
    diag(volatilityModel$sigma) <- diag(volatilityModel$sigma * (1-jumpModel$jumpComponent))
  }

  # Adjust the volatility such that the average stays the same.
  if (burstModel$volModel$modelType == "constantBurst") {
    diag(volatilityModel$sigma) <- diag(volatilityModel$sigma) * (1/(1+((burstModel$volModel$burstInterval[2] - burstModel$volModel$burstInterval[1]) * burstModel$volModel$burstMultiplier)))
  }

  ### Returns coming from mu(t)
  driftReturns <- switch(driftModel$modelType,
                         none = 0,
                         constant = driftModel$drift * dt,
                         Vasicek = vasicekModel(driftModel, nObs, nSeries, nDays, dt)
                         )
  ### Returns coming from sigma(t)
  volatilityReturns <- switch(volatilityModel$modelType,
                              constant = constantVolatilitySim(volatilityModel, nDays, nSeries, nObs, dt),
                              Heston = hestonModel(volatilityModel, nObs, nSeries, nDays, dt),
                              HuangTauchen = huangTauchen(volatilityModel, nObs, nSeries, nDays, dt),
                              LiLinton = liLinton(volatilityModel, nObs, nSeries, nDays, dt)
                              )

  ### Returns from drift bursts
  driftBursts <- switch(burstModel$driftModel$modelType,
                        none = 0,
                        singularityBurst = singularityDriftBurst(burstModel$driftModel, nDays, nSeries, nObs, dt))

  ### What we need to multiply to the volatility returns
  volBursts <- switch(burstModel$volModel$modelType,
                      none = 1,
                      constantBurst = FoFVolatilitySim(burstModel$volModel, nDays, nSeries, nObs, dt),
                      singularityBurst = singularityVolBurst(burstModel$volModel, nDays, nSeries, nObs, dt))

  ### Returns coming from J(t)
  jumps <- switch (jumpModel$modelType,
    none = NULL,
    PA = preAnnouncedJumpSim(jumpModel, nDays, nSeries, nObs)
  )

  diurnality <- switch (diurnalModel$modelType,
                        none = 1,
                        revJ = reverseJDiurnality(diurnalModel, nDays, nSeries, nObs, dt)
  )
  
  
  noise <- switch(noiseModel$noiseType,
                  none = 0,
                  additiveGaussian = NULL,
                  ratio = NULL,
                  ARMA = ARIMAnoise(noiseModel, nDays, nSeries, nObs))
  
  if(volatilityModel$modelType != "LiLinton"){ ## The todorovTauchen model returns both jumps and prices!
    #Construct our returns that comes from volatility and diurnality of this
    returns <- (driftReturns + driftBursts) + volatilityReturns$returns * volBursts * diurnality
    
    # If we need to include jumps, then we do it here
    if(!is.null(jumps)){
      jumps$jumpIndices <- jumps$jumpIndices + 0:(nDays-1) * nObs
      for (j in 1:nSeries){
        returns[jumps$jumpIndices[, j], j] <- returns[jumps$jumpIndices[, j], j] + jumps$jumps[,j]
      }
    }
  } else {
    
    ## Unpack Li Linton values.
    jumps <- volatilityReturns$jumps
    jumpIntensities <- volatilityReturns$jumpIntensities
    volJumps <- volatilityReturns$volJumps
    prices <- volatilityReturns$prices
    sigma <- volatilityReturns$sigma
  }
  
  ## we have the output list
  out <- list()

  ## Preparing the objects to return
  if(spec$returnType == "DT"){
    if(spec$discretize){
      ## Discretized prices
      if(is.null(prices)){
        out$prices <- data.table(DT = as.POSIXct(timestamps, origin = as.POSIXct(spec$timeSettings$origin)), log(round(100 * exp(colCumsum(returns) + noise))/ 100), key = "DT")[]
      } else {
        out$prices <- data.table(DT = as.POSIXct(timestamps, origin = as.POSIXct(spec$timeSettings$origin)), log(round(100 * exp(prices + noise))/ 100), key = "DT")[]
      }
    } else {
      if(is.null(prices)){
        out$prices <- data.table(DT = as.POSIXct(timestamps, origin = as.POSIXct(spec$timeSettings$origin)), colCumsum(returns) + noise, key = "DT")[]
      } else {
        out$prices <- data.table(DT = as.POSIXct(timestamps, origin = as.POSIXct(spec$timeSettings$origin)), prices + noise, key = "DT")[]
      }
      
    }
    setnames(out$prices, c("DT", paste0("PRICE", 1:nSeries)))
    
    if(!is.null(volatilityReturns$sigma)){
      out$sigma <- data.table(DT = as.POSIXct(timestamps, origin = as.POSIXct(spec$timeSettings$origin)), volatilityReturns$sigma * volBursts * diurnality, key = "DT")[]
      setnames(out$sigma, c("DT", paste0("SIGMA", 1:nSeries)))
    }
    if(!is.null(volatilityReturns$volatilityFactor)){
      out$volatilityFactor <- data.table(DT = as.POSIXct(timestamps, origin = as.POSIXct(spec$timeSettings$origin)),  volatilityReturns$volatilityFactor, key = "DT")[]
      setnames(out$volatilityFactor, c("DT", paste0("FACTOR", 1:nSeries)))
    }
    
    if(!is.null(volatilityReturns$volatilityFactor2)){
      out$volatilityFactor2 <- data.table(DT = as.POSIXct(timestamps, origin = as.POSIXct(spec$timeSettings$origin)), VOLFACTOR2 = volatilityReturns$volatilityFactor2, key = "DT")[]
      setnames(out$volatilityFactor2, c("DT", paste0("FACTOR", 1:nSeries)))
    }
    
    if(any(diurnality != 1)){
      out$diurnality <- data.table(DT = as.POSIXct(timestamps, origin = as.POSIXct(spec$timeSettings$origin)), volatilityReturns$diurnality, key = "DT")[]
      setnames(out$diurnality, c("DT", paste0("DIURNALITY")))
    }
    
    if(driftModel$modelType != "none"){
      out$drift <- data.table(DT = as.POSIXct(timestamps, origin = as.POSIXct(spec$timeSettings$origin)), DRIFT = driftReturns + driftBursts, key = "DT")[]
      setnames(out$drift, c("DT", paste0("DRIFT", 1:nSeries)))
    }
    
    if(!is.null(jumps)){
      out$jumps <- data.table(DT = as.POSIXct(timestamps, origin = as.POSIXct(spec$timeSettings$origin)), JUMPS = jumps, key = "DT")[]
      setnames(out$jumps, c("DT", paste0("JUMPS", 1:nSeries)))
    }
    
    
  } else {
    if(is.null(prices)){
      out$prices <- xts(colCumsum(returns), as.POSIXct(timestamps, origin = spec$timeSettings$origin))
    } else {
      out$prices <- xts(prices, as.POSIXct(timestamps, origin = spec$timeSettings$origin))
    }
    
    
    if(spec$discretize){
      out$prices <- log(round(100 * exp(out$prices)) / 100)
    }
    colnames(out$prices) <- paste0("PRICE", 1:nSeries)
    if(!is.null(volatilityReturns$sigma)){
      out$sigma <- volatilityReturns$sigma * volBursts * diurnality
      colnames(out$sigma) <- paste0("SIGMA", 1:nSeries)
    }
    
    if(!is.null(volatilityReturns$volatilityFactor)){
      out$volatilityFactor <- volatilityReturns$volatilityFactor
      colnames(out$volatilityFactor) <- paste0("FACTOR", 1:nSeries)
    }
    if(!is.null(volatilityReturns$volatilityFactor2)){
      out$volatilityFactor2 <- volatilityReturns$volatilityFactor2
      colnames(out$volatilityFactor2) <- paste0("FACTOR", 1:nSeries)
    }
    if(any(diurnality != 1)){
      out$diurnality <- diurnality
      colnames(out$diurnality) <- paste0("DIURNALITY")
    }
    if(driftModel$modelType != "none"){
      out$drift <- driftReturns + driftBursts
      colnames(out$drift) <- paste0("DRIFT", 1:nSeries)
    }

    if(!is.null(jumps)){
      out$jumps <- jumps
      colnames(out$jumps) <- paste0("JUMP", 1:nSeries)
    }
    
    
  }
  
  
  class(out) <- "hfSim"
  attr(out, "spec") <- spec
  toc <- Sys.time()
  out$timer <- toc - tic
  return(out)
}




#' @importFrom mvtnorm rmvnorm
#' @keywords internal
constantVolatilitySim <- function(model, nDays, nSeries, nObs, dt){
  # This is simply a brownian motion with covariance.
  returns <- rmvnorm(nDays * nObs, mean = rep(0, nrow(model$sigma)), model$sigma) * sqrt(dt)
  return(list("returns" = returns))
}

#' @importFrom data.table between
#' @keywords internal
FoFVolatilitySim <- function(model, nDays, nSeries, nObs, dt){
  burstIndices <- round(nObs * model$burstInterval)
  #returns <- rmvnorm(nDays * nObs, mean = rep(0, nrow(model$sigma)), model$sigma)
  #returns[between(1:nObs %% nObs, burstIndices[1], burstIndices[2]),] <- returns[between(1:nObs %% nObs, burstIndices[1], burstIndices[2]),] * sqrt(model$burstModel$burstMultiplier)
  #returns <- returns * sqrt(dt)
  returns <- matrix(rep(rep(1, nSeries), nDays * nObs), ncol = nSeries)
  returns[between(1:nObs %% nObs, burstIndices[1], burstIndices[2]),] <- sqrt(model$burstMultiplier)
  return(returns)

}

#' @importFrom data.table between
#' @keywords internal
singularityVolBurst <- function(model, nDays, nSeries, nObs, dt){
  timestamps <- round(colCumsum(dt) %% 1, 5) # round for numerical stability
  nonBurstIndices <- which(!between(timestamps, model$burstInterval[1], model$burstInterval[2]), arr.ind = TRUE)
  pivot <- mean(model$burstInterval)
  # Volatility burst.
  volBurst <- matrix(model$b * 1/abs(pivot - timestamps)^model$beta, ncol = nSeries)
  volBurst[which(volBurst > 100)] <- mean(volBurst[which(volBurst > 100) + c(-1,1)])
  volBurst[nonBurstIndices] <- 1
  return(volBurst)

}



#' @importFrom data.table between
#' @keywords internal
singularityDriftBurst <- function(model, nDays, nSeries, nObs, dt){
  timestamps <- round(colCumsum(dt) %% 1, 5)  # round for numerical stability
  nonBurstIndices <- which(!between(timestamps, model$burstInterval[1], model$burstInterval[2]), arr.ind = TRUE)
  pivot <- mean(model$burstInterval)
  #pivot <- timestamps[round(mean(round(nObs * model$burstInterval)))] ## Find the pivot and make sure it's actually an obseration, otherwise we may get NaNs or INFS
  # Drift burst
  driftDB <- matrix(model$a * (sign(timestamps - pivot)/(abs(pivot - timestamps)^model$alpha)), ncol = nSeries)
  #driftDB[abs(driftDB) > 10000] <- 0
  driftDB[nonBurstIndices] <- 0
  driftDB[is.nan(driftDB)] <- 0
  return(driftDB * dt)

}



############# Simulate jumps
#' @keywords internal
preAnnouncedJumpSim <- function(model, nDays, nSeries, nObs){
  jumps <- matrix(rnorm(nDays * nSeries, sd = model$jumpVolatility), ncol = nSeries, byrow = TRUE)
  jumpIndices <- round(matrix(sample((model$jumpTime[1] * nObs):(model$jumpTime[2] * nObs), nSeries * nDays, replace = TRUE), nrow = nDays, ncol = nSeries))
  # make sure the jumps don't happen during trading and not after (i.e. we try to put it in to indices that dont exits)
  jumpIndices <- matrix(apply(jumpIndices, 2, FUN = function(x) pmin.int(x, nObs)), ncol = nSeries)
  jumpIndices <- matrix(apply(jumpIndices, 2, FUN = function(x) pmax.int(x, 1)), ncol = nSeries)

  out <- list("jumps" = jumps, "jumpIndices" = jumpIndices)
  return(out)
}


#' @keywords internal
reverseJDiurnality <- function(model, nDays, nSeries, nObs, dt){
  times <- colCumsum(dt) %% 1
  diurnality <- model$C + model$A * exp(-model$a * times) + model$B * exp(-model$b * (1-times))
  return(diurnality)
}


#' @importFrom stats arima.sim
#' @keywords  internal
ARIMAnoise <- function(model, nDays, nSeries, nObs){
  noise <- matrix(NA, nrow = nDays * nObs, ncol = nSeries)
  for (j in 1:nSeries) {
    noise[,j] <- arima.sim(n = nDays * nObs, model = model$model, rand.gen = model$rand.gen, n.start = model$n.start, sd = sqrt(model$variance[j]))
  }
  return(noise)
  
}

## Signal to noise ratio
#' @keywords internal
ratioNoise <- function(model, nDays, nSeries, nObs, sigma){
  noise <- matrix(NA, nrow = nDays * nObs, ncol = nSeries)
  for (j in 1:nSeries) {
    noise[,j] <- rnorm(nDays * nObs, sd = sigma[,j])
  }
  return(noise)
}


## Additive Gaussian noise
#' @keywords internal
additiveGaussianNoise <- function(model, nDays, nSeries, nObs){
  noise <- matrix(NA, nrow = nDays * nObs, ncol = nSeries)
  for (j in 1:nSeries) {
    noise[,j] <- rnorm(nDays * nObs, sd = sqrt(model$variance[j]))
  }
  return(noise)
  
}


#' @export
print.hfSim <- function(x, ...){
  spec <- attr(x, "spec")
  print(spec)
  
  cat("#-----------------------------------------------#\n")
  cat("#----------- High frequency simulation ---------#\n")
  cat("#-----------------------------------------------#\n")
  
  prettyPrintXtsAndData.table(x$prices)
  cat("\nElapsed time\n")
  print(x$timer)
  
  
}

#' @keywords internal
prettyPrintXtsAndData.table <- function(x, n = 5){
  if("xts" %in% class(x)){
    if(nrow(x) < 10){
      print(x)
    } else {
      x <- as.data.table(rbind(head(x, n+1), tail(x,n +1)))
      setnames(x, "index", " ")
      print(x, row.names = FALSE, topn = n)
    }
  } else {
    print(x)
  }
  
}
