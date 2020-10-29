#' create high frequency simulation spec
#' @param volatilityModel List containing the volatility model specification. Possible inputs are: 
#' \enumerate{
#'    \item \code{modelType} A string denoting which model type to use, available models and descriptions thereof are available through the \code{\link{listAvailableVolatilityModels}} function
#'    \item \code{sigma} A symmetric matrix denoting the variance-covariance matrix. Not used when \code{modelType} is "huangTauchen"
#'    \item \code{volOfVol} numeric of length 1 or \code{nSeries} denoting the volatility of the volatility process for one-factor models
#'    \item \code{meanReversion} numeric of or length 1 or \code{nSeries} denoting the strength of mean reversion for the volatility process in one-factor models
#'    \item \code{rho} numeric of length 1 or \code{nSeries} denoting the correlation between \eqn{w_{t}} and \eqn{b_{t}} in the Heston model. It is analogous for \eqn{w_{1}} and \eqn{w_{2}} in the Li Linton model
#'    \item \code{alpha1} numeric of length 1 or \code{nSeries} denoting the strength of mean reversion in the persistent volatility component in the Huang & Tauchen model
#'    \item \code{alpha2} numeric of length 1 or \code{nSeries} denoting the strength of mean reversion in the transient volatility component in the Huang & Tauchen model
#'    \item \code{beta0} numeric of length 1 or \code{nSeries} denoting the level of the total volatility component in the Huang & Tauchen model
#'    \item \code{beta1} numeric of length 1 or \code{nSeries} denoting the level of the persistent volatility component in the Huang & Tauchen model
#'    \item \code{beta2} numeric of length 1 or \code{nSeries} denoting the level of the transient volatility component in the Huang & Tauchen model
#'    \item \code{phi} numeric of length 1 or \code{nSeries} denoting the strength of the feedback in the transient volatility factor in the Huang & Tauchen model
#'    \item \code{rho1} numeric of length 1 or \code{nSeries} denoting the leverage correlations between \eqn{w_{1,t}} and \eqn{w_{2,t}} in the Huang & Tauchen model
#'    \item \code{volJumps} numeric of length 1 or \code{nSeries} denoting the intensity of the jumps in the volatility factor in the Li & Linton model
#'    \item \code{lambda} numeric of length 1 or \code{nSeries} denoting the intensity of the poisson process that governs jumps in the price and volatility in the Li & Linton model
#'    \item \code{mu1} numeric of length 1 or \code{nSeries} denoting the long term mean of the price in the Li & Linton model
#'    \item \code{kappa1} numeric of length 1 or \code{nSeries} denoting the mean reversion of the price process in the Li & Linton model
#' }
#' @param driftModel List containing the drift model specification. Possible inputs are:
#'  \enumerate{
#'    \item \code{modelType} A string denoting which model type to use, available models and descriptions thereof are available through the \code{\link{listAvailableDriftModels}} function
#'    \item \code{drift} numeric of length 1 or \code{nSeries} denoting the daily drift in case modelType = "constant", and the long term mean of drift in the Vasicek model. Default is 0.
#'    \item \code{meanReversion} numeric of length 1 or \code{nSeries} denoting the mean reversion of drift. Default is 2.
#'    \item \code{driftVol} numeric of length 1 or \code{nSeries} denoting the volatility of drift in case modelType = "vasicek".
#' }
#' @param jumpModel List containing the jump model specification. Possible inputs are:
#' \enumerate{
#'    \item \code{modelType} A string denoting which model type to use, available models and descriptions thereof are available through the \code{\link{listAvailableJumpModels}} function
#'    \item \code{jumpComponent} numeric of length 1 or \code{nSeries} between 0 and 1 denoting how much of the total sigma of the model that should come from the jump variation.
#' }
#' @param diurnalModel List containing the diurnal model specification. Possible inputs are:
#' \enumerate{
#'    \item \code{modelType} A string denoting which model type to use, available models and descriptions thereof are available through the \code{\link{listAvailableDiurnalModels}} function.
#'    \item \code{C} numeric of length 1 denoting the level of the diurnal pattern. (Usually calibrated to make sure that when integrating over the diurnal factor, it gives 1)
#'    \item \code{A} numeric of length 1 denoting the A value of the diurnal pattern. This is used to control the level of the diurnal factor during the beginning of the day. Default is 0.75
#'    \item \code{B} numeric of length 1 denoting the B value of the diurnal pattern. This is used to control the level of the diurnal factor during the end of the day. Default is 0.25
#'    \item \code{a} numeric of length 1 denoting the a value of the diurnal pattern. This is used to control the steepness of the drop during the beginning of the day. Default is 10
#'    \item \code{b} numeric of length 1 denoting the b value of the diurnal pattern. This is used to control the steepness of the increase during the end of the day. Default is 10
#' }
#' @param burstModel List containing the burst model specifications. Possible inputs are: 
#' \enumerate{
#'    \item driftModel list with entries 
#'    \enumerate{
#'          \item \code{modelType} A string denoting which model type to use, available models and descriptions thereof are available through the \code{\link{listAvailableBurstModels}} function
#'          \item \code{alpha} numeric of length one denoting the alpha exponent in the drift burst. Default is 0.65 
#'          \item \code{a} numeric of length one denoting the a scalar in the drift burst this can be used to control the severity of the drift burst. Default is 1 
#'          \item \code{burstInterval} numeric of length two denoting the interval when the drift burst is to take place. Default is \code{c(15/32, 17/32)} 
#'    }
#'    \item volModel list with entries 
#'    \enumerate{
#'          \item \code{modelType} A string denoting which model type to use, available models and descriptions thereof are available through the \code{\link{listAvailableBurstModels}} function
#'          \item \code{beta} (Used only when modelType is "singularityBurst") A numeric of length one denoting the beta exponent in the volatility burst. Default is 0.4 
#'          \item \code{b} (Used only when modelType is "singularityBurst") A numeric of length one denoting the b scalar in the volatility burst this can be used to control the severity of the volatility burst. Default is 1 
#'          \item \code{burstMultiplier} (Used only when modelType is "constantBurst") A numeric of length one denoting the severity of the volatility burst. Default is 3 denoting a temporary three-fold increase of volatility
#'          \item \code{burstInterval} A numeric of length two when the volatility burst is to take place. Default is \code{c(15/32, 17/32)} 
#'    }
#' }
#' @param noiseModel List containing the noise model specification. Possible inputs are: 
#' \enumerate{
#'    \item \code{modelType} A string denoting which model type to use, available models and descriptions thereof are available through the \code{\link{listAvailableNoiseModels}} function.
#'    \item \code{signalToNoise} A numeric of either length 1, nrow(sigma) or nSeries denoting the signal to noise ratio of the noise.
#'    \item \code{variance} A numeric of either length 1, nrow(sigma) or nSeries denoting the variance of the noise when the variance is constant.
#' }
#' @param timeSettings List containing the settings for the time in the simulation. Possible inputs are:
#' \enumerate{
#'    \item \code{tradingStart} A numeric of length one denoting in seconds after midnight the start of trading. Default is 34200 which corresponds to 09:30
#'    \item \code{tradingEnd} A numeric of length one denoting in seconds after midnight the start of trading. Default is 57600 which corresponds to 16:00
#'    \item \code{origin} A character which can be coerced to a date, denoting the first day of the simulated series. Default is "1970-01-01"
#'    \item \code{sampling} A character denoting which sampling scheme to use. Currently only "equidistant" is implemented.
#' }
#' @param nSeries Integer-valued numeric of length one denoting how many series to simulate. This is overwritten in case nrow of the variance covariance matrix is greater than one. Then, nSeries will be set to match \code{nrow(sigma)}
#' @param nDays Integer-valued numeric of length one denoting how many days to simulate the data over.
#' @param nObs Integer-valued numeric of length one denoting how many days to simulate the data over.
#' @param discretize Logical denoting whether to discretize the prices. The prices are discretized using log(round(100 * exp(prices)) / 100). i.e. to nearest cent.
#' @param returnType character determining the type of object to use for the storage of returns, prices, jumps etc. Defaults to "DT", which corresponds to \code{data.table} any other input will return \code{xts} objects.
#' @returns an object of type "highfrequencySimSpec"
#' @author Emil Sjoerup
#' @export
#' @usage 
#' hfSimSpec(volatilityModel = list(modelType = "constant", sigma = 0.2),
#'           driftModel = list(modelType = "none", drift = 0),
#'           jumpModel = list(modelType = "none"),
#'           diurnalModel = list(modelType = "none"),
#'           burstModel = list(driftModel = list(modelType = "none"), 
#'                             volModel = list(modelType = "none")),
#'           noiseModel = list(noiseType = "none"),
#'           timeSettings = list(tradingStart = 34200, tradingEnd = 57600, 
#'                               origin = "1970-01-01", sampling = "equidistant"),
#'           nSeries = 1,
#'           nDays = 1,
#'           nObs = 23401,
#'           discretize = FALSE,
#'           returnType = "DT")
hfSimSpec <- function(volatilityModel = list(modelType = "constant", sigma = 0.2),
                      driftModel = list(modelType = "none", drift = 0),
                      jumpModel = list(modelType = "none"),
                      diurnalModel = list(modelType = "none"),
                      burstModel = list(driftModel = list(modelType = "none"), 
                                        volModel = list(modelType = "none")),
                      noiseModel = list(noiseType = "none"),
                      timeSettings = list(tradingStart = 34200, tradingEnd = 57600, 
                                          origin = "1970-01-01", sampling = "equidistant"),
                      nSeries = 1, nDays = 1, nObs = 23401, discretize = FALSE, returnType = "DT"){
  ## The code for ensuring the models are valid is inspired by code from Alexios Ghalanos' rugarch.

  ######## Volatility model checking starts ######## Jump model checking ends ########
  vm <- match(names(volatilityModel), c("modelType","meanReversion", "sigma","volOfVol", "rho", "alpha1", "alpha2", "beta0",
                                        "beta1", "beta2", "phi", "rho1", "rho2", "rho", "kappa1", "delta", "lambda", "mu1"))
  if(any(is.na(vm))){
    idx <- which(is.na(vm))
    enx <- NULL
    for (i in 1:length(idx)) enx <- c(enx, names(volatilityModel)[idx[i]])
    warning(paste(c("unidentified option(s) in volatilityModel:\n", enx), sep="", collapse= " "), call. = FALSE, domain=NULL)
  }

  # If we don't get the input, we set it to the defaults
  if(is.null(volatilityModel$modelType)) volatilityModel$modelType <- "constant"
  if(volatilityModel$modelType == "LiLinton"){
    if(is.null(volatilityModel$sigma)) volatilityModel$sigma <- matrix(rep(0.2, nSeries), ncol = nSeries)
  } else {
    if(is.null(volatilityModel$sigma)) volatilityModel$sigma <- diag(0.2, nrow = nSeries, ncol = nSeries)
  }
  browser()
  if(nSeries != nrow(volatilityModel$sigma)){
    if(nrow(volatilityModel$sigma) != 1){
      warning("nSeries does not match the number of assets mandated in the sigma matrix. Setting nSeries to match nrow(sigma)")
      nSeries <- nrow(volatilityModel$sigma)
    } else { ## We have 1 value for sigma, so we make it into a nSeries by nSeries diagonal matrix
      if(volatilityModel$modelType == "LiLinton"){ ## Here we just repeat the vector until it fits in the length.
        volatilityModel$sigma <- matrix(rep(volatilityModel$sigma, nSeries)[1:nSeries], ncol = nSeries)
      } else {
        volatilityModel$sigma <- diag(as.numeric(volatilityModel$sigma), ncol = nSeries, nrow = nSeries)
      }
    }
  }
  
  volatilityModel$sigma <- matrix(volatilityModel$sigma, ncol = nSeries)


  validVolatilityModelTypes <- listAvailableVolatilityModels()[,1]
  if(!is.character(volatilityModel$modelType)){
    stop("volatility model type must be a character.\n", call.=FALSE)
  }
  if(!(volatilityModel$modelType %in% validVolatilityModelTypes)){
    stop("volatility model type specified does not appear in valid model types. See listAvalableVolatilityModels() for valid types.\n", call.=FALSE)
  }
  if(!is.numeric(volatilityModel$sigma) | any(diag(volatilityModel$sigma) < 0) | (!isSymmetric(volatilityModel$sigma) & (nrow(volatilityModel$sigma) > 1))){
    stop("sigma must be symmetric matrix with non-negative diagonal")
  }
  
  if(nrow(volatilityModel$sigma) > 1){
    ev <- eigen(volatilityModel$sigma, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("sigma is numerically not positive semidefinite")
    }
  }



  if(volatilityModel$modelType == "Heston"){
    # Defaults based on the slides of Bezirgen Veliyev for the 2018 high frequency econometrics course at Aarhus University:
    if(is.null(volatilityModel$meanReversion)) volatilityModel$meanReversion <- rep(5/250, nSeries)
    if(is.null(volatilityModel$sigma)) volatilityModel$sigma <- diag(0.2/250, nSeries)
    if(is.null(volatilityModel$volOfVol)) volatilityModel$volOfVol <- rep(0.5/250, nSeries)
    if(is.null(volatilityModel$rho)) volatilityModel$rho <- rep(-0.5, nSeries)
  }


  if(volatilityModel$modelType == "HuangTauchen"){
    if(is.null(volatilityModel$alpha1)) volatilityModel$alpha1 <- rep(-0.00137, nSeries)
    if(is.null(volatilityModel$alpha2)) volatilityModel$alpha2 <- rep(-1.386, nSeries)
    if(is.null(volatilityModel$beta0)) volatilityModel$beta0 <- rep(-1.2, nSeries)
    if(is.null(volatilityModel$beta1)) volatilityModel$beta1 <- rep(0.04, nSeries)
    if(is.null(volatilityModel$beta2)) volatilityModel$beta2 <- rep(1.5, nSeries)
    if(is.null(volatilityModel$phi)) volatilityModel$phi <- rep(0.25, nSeries)
    if(is.null(volatilityModel$rho1)) volatilityModel$rho1 <- rep(-0.3, nSeries)
    if(is.null(volatilityModel$rho2)) volatilityModel$rho2 <- rep(-0.3, nSeries)
  }

  if(volatilityModel$modelType == "LiLinton"){
    if(is.null(volatilityModel$meanReversion)) volatilityModel$meanReversion <- rep(5/252, nSeries)
    if(is.null(volatilityModel$mu1)) volatilityModel$mu1 <- rep(3.6, nSeries)
    if(is.null(volatilityModel$kappa1)) volatilityModel$kappa1 <- rep(0.5, nSeries)
    if(is.null(volatilityModel$sigma)) volatilityModel$sigma <- rep(0.04/252, nSeries)
    if(is.null(volatilityModel$volOfVol)) volatilityModel$volOfVol <- rep(0.05/252, nSeries)
    if(is.null(volatilityModel$rho)) volatilityModel$rho <- rep(-0.5, nSeries)
    if(is.null(volatilityModel$lambda)) volatilityModel$lambda <- rep(3, nSeries)
    if(is.null(volatilityModel$volJumps)) volatilityModel$volJumps <- volatilityModel$volOfVol
    jumpModel$modelType <- "none" ## We don't allow for separate jumps.
    burstModel$driftModel$modelType <- NULL # We don't allow for drift bursts
    burstModel$volModel$modelType <- NULL # We don't allow for vol bursts
    
    
    if(nrow(volatilityModel$sigma) > 1){
      stop("sigma must be a numeric with 1 row. The LiLinton model is not implemented to allow for correlation between assets.")
    }
  }


  ######## Drift model checking starts ######## Jump model checking ends ########
  dm <- match(names(driftModel), c("modelType", "drift", "meanReversion", "driftVol"))
  if(any(is.na(dm))){
    idx <- which(is.na(dm))
    enx <- NULL
    for (i in 1:length(idx)) enx <- c(enx, names(driftModel)[idx[i]])
    warning(paste(c("unidentified option(s) in driftModel:\n", enx), sep="", collapse= " "), call. = FALSE, domain=NULL)
  }

  # If we don't get the input, we set it to the defaults
  if(is.null(driftModel$modelType)) driftModel$modelType <- "none"
  if(is.null(driftModel$drift)) driftModel$drift <- 0
  # drift must be same length as sigma has columns
  if(length(driftModel$drift) != ncol(volatilityModel$sigma)) driftModel$drift <- rep(driftModel$drift, ncol(volatilityModel$sigma))[1:ncol(volatilityModel$sigma)]
  if(!is.character(driftModel$modelType)){
    stop("drift model type must be a character.\n", call.=FALSE)
  }
  validDriftModelTypes <- listAvailableDriftModels()[,1]
  if(!(driftModel$modelType %in% validDriftModelTypes)){
    stop("drift model type specified does not appear in valid model types. See listAvailableDriftModels() for valid types.\n", call.=FALSE)
  }

  ## Check that vasicek model is correctly specified
  if(driftModel$modelType == "Vasicek"){

    if(is.null(driftModel$drift)) driftModel$drift = rep(0, ncol(volatilityModel$sigma))
    if(is.null(driftModel$meanReversion)) driftModel$meanReversion = rep(2, ncol(volatilityModel$sigma))
    if(is.null(driftModel$driftVol)) driftModel$driftVol = rep(0.0391, ncol(volatilityModel$sigma))

    if(length(driftModel$drift) != ncol(volatilityModel$sigma)) driftModel$drift <- rep(driftModel$drift, ncol(volatilityModel$sigma))[1:ncol(volatilityModel$sigma)]
    if(length(driftModel$meanReversion) != ncol(volatilityModel$sigma)) driftModel$meanReversion <- rep(driftModel$meanReversion, ncol(volatilityModel$sigma))[1:ncol(volatilityModel$sigma)]
    if(length(driftModel$driftVol) != ncol(volatilityModel$sigma)) driftModel$driftVol <- rep(driftModel$driftVol, ncol(volatilityModel$sigma))[1:ncol(volatilityModel$sigma)]
  }

  if(!is.numeric(driftModel$drift) | (length(driftModel$drift) != 1 & length(driftModel$drift) != nSeries) & length(driftModel$drift) != ncol(volatilityModel$sigma)){
    stop("drift must be a numeric with length equal to 1 or equal to nSeries, or ncol(sigma) ")
  }
  ######## Drift model checking ends


  ######## Jump model checking starts ######## Jump model checking ends ########
  jm <- match(names(jumpModel), c("modelType", "jumpComponent", "jumpTime", "includeJumps"))
  if(any(is.na(jm))){
    idx <- which(is.na(jm))
    enx <- NULL
    for (i in 1:length(idx)) enx <- c(enx, names(jumpModel)[idx[i]])
    warning(paste(c("unidentified option(s) in jumpModel:\n", enx), sep="", collapse= " "), call. = FALSE, domain=NULL)
  }


  # We set the includeJumps tag to TRUE if we should actually include jumps
  # And we set it FALSE if we should not.
  if(jumpModel$modelType != "none"){
    jumpModel$includeJumps <- TRUE
  } else if(jumpModel$modelType == "none"){
    jumpModel$includeJumps <- FALSE
  }


  if(is.null(jumpModel$modelType)) jumpModel$modelType <- "none"
  if(is.null(jumpModel$jumpComponent)) jumpModel$jumpComponent <- 1/5
  if(is.null(jumpModel$jumpTime)) jumpModel$jumpTime <- c(15/32, 17/32)

  # Check whether the jump component makes sense.
  if(jumpModel$includeJumps && jumpModel$jumpComponent >= 1){
    stop("Jump component equal to or greater than 1, this is not allowed. Jump component should be between 0 and 1", call. = FALSE, domain=NULL)
  }

  if(jumpModel$includeJumps && jumpModel$jumpComponent <= 0){
    stop("Jump component equal to or less than 0, this is not allowed. Jump component should be between 0 and 1", call. = FALSE, domain=NULL)
  }

  if(jumpModel$includeJumps && length(jumpModel$jumpTime) != 2){
    stop("jumpTime must be of length 2", call.=FALSE, domain = NULL)
  }

  if(jumpModel$includeJumps && jumpModel$jumpTime[1] > jumpModel$jumpTime[2]){
    stop("First entry of jumpTime must be lower than second entry of jumpTime", call.=FALSE, domain=NULL)
  }

  if(jumpModel$includeJumps && (min(jumpModel$jumpTime)<0 | max(jumpModel$jumpTime)>1)){
    stop("jumpTime must be between 0 and 1", call. = FALSE, domain=NULL)
  }

  


  ######### Burst model checking starts ########
  bm <- match(names(burstModel), c("driftModel", "volModel"))
  if(any(is.na(bm))){
    idx <- which(is.na(bm))
    enx <- NULL
    for (i in 1:length(idx)) enx <- c(enx, names(burstModel)[idx[i]])
    warning(paste(c("unidentified option(s) in burstModel:\n", enx), sep="", collapse= " "), call. = FALSE, domain=NULL)
  }

  if(is.null(burstModel$driftModel$modelType)) burstModel$driftModel$modelType = "none"
  if(is.null(burstModel$volModel$modelType)) burstModel$volModel$modelType = "none"

  validBurstModels <- listAvailableBurstModels()[,1]
  if(!(burstModel$driftModel$modelType %in% validBurstModels[c(1,3)])){
    stop("burst model for the drift component's model type specified does not appear in valid model types. See listAvailableBurstModels() for valid types.\n", call.=FALSE)
  }
  if(!(burstModel$volModel$modelType %in% validBurstModels)){
    stop("burst model for the volatility component's model type specified does not appear in valid model types. See listAvailableBurstModels() for valid types.\n", call.=FALSE)
  }

  if(burstModel$driftModel$modelType == "singularityBurst" && is.null(burstModel$driftModel)) burstModel$driftModel <- list(alpha = 0.6, a = 0.002,
                                                                                                                 burstInterval = c(15/32, 17/32), modelType = "singularityBurst")
  if(burstModel$driftModel$modelType == "singularityBurst" && !is.null(burstModel$driftModel)){
    if(is.null(burstModel$driftModel$alpha)) burstModel$driftModel$alpha <- 0.6
    if(length(burstModel$driftModel$alpha) != 1){
      stop("the alpha scalar in the driftModel of the burstModel must be of length 1", call.=FALSE, domain = NULL)
    }
    if(is.null(burstModel$driftModel$a)) burstModel$driftModel$a <- 1
    if(length(burstModel$driftModel$a) != 1){
      stop("the a scalar in the driftModel of the burstModel must be of length 1", call.=FALSE, domain = NULL)
    }
    if(is.null(burstModel$driftModel$burstInterval)) burstModel$driftModel$burstInterval <- c(15/32, 17/32)
    if(length(burstModel$driftModel$burstInterval) != 2){
      stop("burstInterval must be of length 2", call.=FALSE, domain = NULL)
    }
  }



  # We have received a burst model, which we must check
  if(burstModel$volModel$modelType == "constantBurst"){
    if(is.null(burstModel$volModel$burstMultiplier)) burstModel$volModel$burstMultiplier <- 3
    if(length(burstModel$volModel$burstMultiplier) != 1){
      stop("the burstMultiplier scalar in the volModel of the burstModel must be of length 1", call.=FALSE, domain = NULL)
    }
    if(is.null(burstModel$volModel$burstInterval)) burstModel$volModel$burstInterval <- c(15/32, 17/32)
    if(length(burstModel$volModel$burstInterval) != 2){
      stop("burstInterval must be of length 2", call.=FALSE, domain = NULL)
    }
  }


  if(burstModel$volModel$modelType == "singularityBurst"){
    if(is.null(burstModel$volModel$b)) burstModel$volModel$b <- 1
    if(length(burstModel$volModel$b) != 1){
      stop("the b scalar in the volModel of the burstModel must be of length 1", call.=FALSE, domain = NULL)
    }
    if(is.null(burstModel$volModel$beta)) burstModel$volModel$beta <- 0.4
    if(length(burstModel$volModel$beta) != 1){
      stop("the beta scalar in the volModel of the burstModel must be of length 1", call.=FALSE, domain = NULL)
    }
    if(is.null(burstModel$volModel$burstInterval)) burstModel$volModel$burstInterval <- c(15/32, 17/32)
    if(length(burstModel$volModel$burstInterval) != 2){
      stop("burstInterval must be of length 2", call.=FALSE, domain = NULL)
    }
  }


  if(burstModel$driftModel$modelType == "singularityBurst" && burstModel$volModel$modelType == "singularityBurst"){
    if(burstModel$volModel$beta < 0 | burstModel$volModel$beta > 1/2 | burstModel$driftModel$alpha - burstModel$volModel$beta > 1/2){
      warning("The singularity burst model set up does not ensure absence of arbitrage. The conditions required are: 0 < beta < 1/2 and alpha - beta < 1/2.")
    }
    if(!all(burstModel$driftModel$burstInterval == burstModel$volModel$burstInterval)){
      warning("The singularity burst interval of the drift does not match the singularity burst of the volatility.
               Typically the drift and volatility co-explode.")
    }
  }



  ######## diurnal model checking starts ########
  dm <- match(names(diurnalModel), c("modelType","C", "A", "B", "a", "b"))
  if(any(is.na(dm))){
    idx <- which(is.na(dm))
    enx <- NULL
    for (i in 1:length(idx)) enx <- c(enx, names(diurnalModel)[idx[i]])
    warning(paste(c("unidentified option(s) in diurnalModel:\n", enx), sep="", collapse= " "), call. = FALSE, domain=NULL)
  }

  if(is.null(diurnalModel$modelType)) diurnalModel$modelType = "none"

  validDiurnalModels <- listAvailableDiurnalModels()[,1]
  if(!(diurnalModel$modelType %in% validDiurnalModels)){
    stop("diurnal model type specified does not appear in valid model types. See listAvailableDiurnalModels() for valid types.\n", call.=FALSE)
  }

  # We have received a diurnal model, which we must check
  if(diurnalModel$modelType == "revJ"){
    if(is.null(diurnalModel$C)) diurnalModel$C <- 0.88929198
    if(is.null(diurnalModel$A)) diurnalModel$A <- 0.75
    if(is.null(diurnalModel$B)) diurnalModel$B <- 0.25
    if(is.null(diurnalModel$a)) diurnalModel$a <- 10
    if(is.null(diurnalModel$b)) diurnalModel$b <- 10
  }

  
  ######## noise model checking starts ########
  nm <- match(names(noiseModel), c("noiseType", "signalToNoiseRatio", "variance", "model", "rand.gen", "n.start", "extra"))
  if(any(is.na(nm))){
    idx <- which(is.na(nm))
    enx <- NULL
    for (i in 1:length(nm)) enx <- c(enx, names(noiseModel)[idx[i]])
    warning(paste(c("unidentified option(s) in noiseModel:\n", enx), sep="", collapse= " "), call. = FALSE, domain=NULL)
  }
  
  if(is.null(noiseModel)){
    noiseModel$noiseType <- "none"
  }
  
  if(noiseModel$noiseType == "additiveGaussian"){
    if(is.null(noiseModel$variance)) noiseModel$variance <- 0.02/252
    noiseModel$variance <- rep(noiseModel$variance, length(nSeries))[1:nSeries] ## Make sure we have the right amount of variances in the noise model - one for each series
  }
  
  if(noiseModel$noiseType == "ratio"){
    if(is.null(noiseModel$signalToNoiseRatio)) noiseModel$signalToNoiseRatio <- 2
  }
  
  if(noiseModel$noiseType == "ARMA"){
    if(is.null(noiseModel$model)) noiseModel$model <- list(ar = 0.7)
    if(is.null(noiseModel$variance)) noiseModel$variance <- 0.3 * 1e-4
    noiseModel$variance <- rep(noiseModel$variance, nSeries)[1:nSeries] ## Make sure we have the right amount of variances in the noise model
    if(is.null(noiseModel$rand.gen)) noiseModel$rand.gen <- rnorm
    if(is.null(noiseModel$n.start)) noiseModel$n.start <- NA
    if(is.null(noiseModel$extra)) noiseModel$extra <- list(NULL)
  }
  

  ######## time settings checking starts ########
  ts <- match(names(timeSettings), c("tradingStart", "tradingEnd", "origin", "sampling"))
  if(any(is.na(ts))){
    idx <- which(is.na(ts))
    enx <- NULL
    for (i in 1:length(idx)) enx <- c(enx, names(timeSettings)[idx[i]])
    warning(paste(c("unidentified option(s) in timeSettings:\n", enx), sep="", collapse= " "), call. = FALSE, domain=NULL)
  }

  if(is.null(timeSettings$tradingStart)) timeSettings$tradingStart <- 34200
  if(length(timeSettings$tradingStart) != 1 | !is.numeric(timeSettings$tradingStart)){
    stop("the tradingStart argument must be a numeric of length 1", call.=FALSE, domain = NULL)
  }

  if(is.null(timeSettings$tradingEnd)) timeSettings$tradingEnd <- 57600
  if(length(timeSettings$tradingEnd) != 1 | !is.numeric(timeSettings$tradingEnd)){
    stop("the tradingEnd argument must be a numeric of length 1", call.=FALSE, domain = NULL)
  }

  if(is.null(timeSettings$origin)) timeSettings$origin <- "1970-01-01"
  if(length(timeSettings$origin) != 1){
    stop("the origin argument must be of length 1", call.=FALSE, domain = NULL)
  }
  timeSettings$origin <- try(as.Date(timeSettings$origin), silent = TRUE)
  if (inherits(timeSettings$origin, "try-error")) {
    stop("origin must be coercible to a date")
  }
  timeSettings$origin <- as.character(timeSettings$origin)


  if(is.null(timeSettings$sampling)) timeSettings$sampling <- "equidistant"


  if(timeSettings$tradingEnd < timeSettings$tradingStart){
    stop("Start of trading must be before end of trading")
  }

  if(!is.numeric(nSeries) | nSeries <= 0 | nSeries %% 1 != 0){
    stop("nSeries must be a positive integer", call. = FALSE, domain = NULL)
  }

  if(!is.numeric(nDays) | nDays <= 0 | nDays %% 1 != 0){
    stop("nDays must be a positive integer", call. = FALSE, domain = NULL)
  }

  if(!is.numeric(nObs) | nObs <= 0 | nObs %% 1 != 0){
    stop("nObs must be a positive integer", call. = FALSE, domain = NULL)
  }

  if(!is.logical(discretize)){
    stop("discretize must a logical", call. = FALSE, domain = NULL)
  }
  
  if(!returnType %in% c("DT", "xts")) stop('returnType must be either \"DT\" or \"xts\" corresponding to data.table and xts objects respectively.')
  simSpec <- list(volatilityModel = volatilityModel, driftModel = driftModel, jumpModel = jumpModel, burstModel = burstModel, diurnalModel = diurnalModel,
                  noiseModel = noiseModel, nSeries = nSeries, nDays = nDays, nObs = nObs, discretize = discretize, timeSettings = timeSettings, returnType = returnType)

  # Should probably be an S4 class, but it's bad practice to mix and match S3 and S4 classes AFAIK
  class(simSpec) <- "highfrequencySimSpec"

  return(simSpec)

}


#' List the available drift models for simulations
#' @export
#' @author Emil Sjoerup
#' @return This function returns the available drift models in a matrix
listAvailableDriftModels <- function(){
  models <- matrix(
    c("none", "no drift",
      "constant", "constant drift",
      "Vasicek", "Vasicek model"
      ), ncol = 2, byrow=TRUE
  )

  colnames(models) <- c("Abbreviation", "Description")
  return(models)
}


#' List the available burst models for simulations
#' @export
#' @author Emil Sjoerup
#' @return This function returns the available burst models in a matrix
listAvailableBurstModels <- function(){
  models <- matrix(
    c("none", "no burst",
      "constantBurst", "Piecewise constant volatility with a constant burst of volatility in a pre-defined interval (FoF) (Only available for volatility)",
      "singularityBurst", "Burst that approaches a singularity as in CRO(2018)"), ncol = 2, byrow=TRUE
  )
  colnames(models) <- c("Abbreviation", "Description")
  return(models)
}


#' List the available diurnal models for simulations
#' @export
#' @author Emil Sjoerup
#' @return This function returns the available diurnal models in a matrix
listAvailableDiurnalModels <- function(){
  models <- matrix(
    c(
      "none", "no diurnality",
      "revJ", "reverse J shaped diurnality (FORM: C + A * exp(-a * t) + B * exp(-b * (1 - t)) )"
    ), ncol = 2, byrow = TRUE
  )
}

#' List the available noise models for simulations
#' @export
#' @author Emil Sjoerup
#' @return This function returns the available noise models in a matrix
listAvailableNoiseModels <- function(){
  models <- matrix(
    c(
      "none", "no noise",
      "additiveGaussian", "additive Gaussian noise term",
      "ratio", "additive Gaussian noise with constant signal to noise ratio",
      "ARMA", "AR(I)MA noise simulated by the arima.sim function"
    ), ncol = 2, byrow = TRUE
  )
  colnames(models) <- c("Abbreviation", "Description")
  return(models)
}


#' List the available volatility models for simulations
#' @export
#' @author Emil Sjoerup
#' @return This function returns the available volatility models in a matrix
listAvailableVolatilityModels <- function(){
  models <- matrix(
    c("constant", "constant volatility",
      "Heston", "heston stochastic volatility model",
      "HuangTauchen", "two factor stochastic volatility model of Huang and Tauchen (2005)",
      "LiLinton", "One-factor with simultaneous jumps in both price and volatility Li and Linton (2020)"),
    ncol = 2, byrow=TRUE
  )

  colnames(models) <- c("Abbreviation", "Description")

  return(models)
}

#' List the available jump models for simulations
#' @export
#' @author Emil Sjoerup
#' @return This function returns the available Jump models in a matrix
listAvailableJumpModels <- function(){
  models <- matrix(
    c("none", "No jumps",
      "PA", "Single daily pre announced jump"), ncol = 2, byrow=TRUE
  )
  colnames(models) <- c("Abbreviation", "Description")
  return(models)
}


#' @export
print.highfrequencySimSpec <- function(x, ...){
  cat(  "#-----------------------------------------------#")
  cat("\n#-------- High frequency simulation spec -------#")
  cat("\n#-----------------------------------------------#")
  cat("\nDrift model type:", x$driftModel$modelType)
  # if(x$driftModel$modelType != "none"){
  #   drift <- x$driftModel$drift
  #   names(drift) <- paste0("PRICE", 1:x$nSeries)
  #   if(x$driftModel$modelType != "Vasicek"){
  #     cat("#- Drift components: -#\n")
  #   } else {
  #     cat("#- Mean drift level: -#\n")
  #   }
  #   print(drift)
  #   if(x$driftModel$modelType == "Vasicek"){
  #     MR <- x$driftModel$meanReversion
  #     names(MR) <- paste0("PRICE", 1:x$nSeries)
  #     DV <- x$driftModel$driftVol
  #     names(DV) <- paste0("PRICE", 1:x$nSeries)
  #     cat("#- Mean reversion of drift: -#\n")
  #     print(MR)
  #     cat("#- Volatility of drift: -# \n")
  #     print(DV)
  #   }
  #   
  #   
  # }
  cat("\nVolatility model type:", x$volatilityModel$modelType)
  # cat("#- Covariance matrix: -#\n")
  # sig <- x$volatilityModel$sigma
  # if(nrow(sig) > 1){
  #   colnames(sig) <- rownames(sig) <- paste0("PRICE", 1:x$nSeries)
  # } else {
  #   colnames(sig) <- paste0("PRICE", 1:x$nSeries)
  #   rownames(sig) <- " "
  # }
  # 
  # print(sig)
  # cat("#-----------------------------------------------#\n")
  if(x$jumpModel$modelType != "none"){
    cat("\nJump model type:", x$jumpModel$modelType)
  } else {
    cat("\nJumps: No jumps included")
  }
  if(x$burstModel$driftModel$modelType != "none"){
    cat("\nDrift burst model type: Singularity")
  }
  if(x$burstModel$volModel$modelType != "none"){
    if(x$burstModel$volModel$modelType == "singularityBurst") cat("\nVolatility burst model type: Singularity")
    if(x$burstModel$volModel$modelType == "constant") cat("\nVolatility burst model type: Constant burst")
  }
  
  if(x$diurnalModel$modelType == "revJ"){
    cat("\nDiurnal pattern: Reverse J")
  } else {
    cat("\nDiurnality: no diurnal pattern")
  }
  
  if(x$noiseModel$noiseType != "none"){
    cat("\nNoise model:", x$noiseModel$noiseType)
  } else {
    cat("\nNoise: no noise included")
  }
  
  if(x$discretize){
    cat("\nPrices are discretized to nearest cent")
  } else {
    cat("\nPrices are not discretized")
  }
  
  
  cat("\nTrading opens: ", format(as.POSIXct(x$timeSettings$tradingStart, origin = as.POSIXct(x$timeSettings$origin)), "%H:%M:%S"))
  cat("\nTrading closes:", format(as.POSIXct(x$timeSettings$tradingEnd, origin = as.POSIXct(x$timeSettings$origin)), "%H:%M:%S"))
  cat("\nNumber of observations to simulate", x$nSeries * x$nDays * x$nObs, "\n")
  invisible(x)
}
