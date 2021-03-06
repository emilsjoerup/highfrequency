#' @title Get high frequency data from Alpha Vantage
#' 
#' @description Function to retrieve high frequency data from Alpha Vantage - wrapper around quantmod's getSymbols.av function
#' 
#' @param symbols character vector with the symbols to import.
#' @param interval the sampling interval of the data retrieved. Should be one of one of "1min", "5min", "15min", "30min", or "60min"
#' @param outputType string either \code{"xts"} or \code{"DT"} to denote the type of output wanted. \code{"xts"} will yield an xts object, \code{"DT"} will yield a data.table object.
#' @param apiKey string with the api key provided by Alpha Vantage. 
#' @param doSleep logical when the length of symbols > 5 the function will sleep for 12 seconds by default. 
#' 
#' @details The \code{doSleep} argument is set to true as default because Alpha Vantage has a limit of five calls per minute. 
#' The function does not try to extract when the last API call was made which means that if
#' you made successive calls to get 3 symbols in rapid succession, the function may not retrieve all the data.
#' 
#' @return An object of type \code{xts} or \code{data.table} in case the length of symbols is 1. If the length of symbols > 1 the \code{xts} and 
#' \code{data.table} objects will be put into a list.
#' 
#' @author Emil Sjoerup (wrapper only) Paul Teetor (for quantmod's getSymbols.av)
#' 
#' @seealso The getSymbols.av function in the quantmod package
#' @examples
#' \dontrun{
#' # Get data for SPY at an interval of 1 minute in the standard xts format.
#' data <- getAlphaVantageData(symbols = "SPY", apiKey = "yourKey", interval = "1min")
#' 
#' # Get data for 3M and Goldman Sachs  at a 5 minute interval in the data.table format. 
#' # The data.tables will be put in a list.
#' data <- getAlphaVantageData(symbols = c("MMM", "GS"), interval = "5min", 
#'                   outputType = "DT", apiKey = 'yourKey')
#' 
#' # Get data for JPM and Citicorp at a 15 minute interval in the xts format. 
#' # The xts objects will be put in a list.
#' data <- getAlphaVantageData(symbols = c("JPM", "C"), interval = "15min", 
#'                   outputType = "xts", apiKey = "yourKey")
#' }
#' 
#' @importFrom data.table as.data.table setnames
#' @keywords data
#' @export
getAlphaVantageData <- function(symbols = NULL, interval = "5min", outputType = "xts", apiKey = NULL, doSleep = TRUE) {
  # Check for API key set in quantmod if it is not present, we kindly remind the user to get it.
  if (is.null(apiKey) && is.null(quantmod::getDefaults(quantmod::getSymbols.av, "api.key")$api.key)) {
    stop("Your AlphaVantage API key is not set in the quantmod package and no API key was provided. 
          An API key is required but can be obtained for free.
          Please see the getSymbols.av function in the quantmod package for details on how to set a default API key for the current session.
          Alternatively, use the apiKey argument to provide the key.")
  }
  if (is.null(apiKey)) {
    apiKey <- quantmod::getDefaults(quantmod::getSymbols.av, "api.key")$api.key #quantmod handles the extra quotes.
  }
  
  
  # Check if we have more than one symbol
  isMultiSymbols <- length(symbols) > 1
  
  # We have more than one symbol
  if (isMultiSymbols) {
    data <- vector("list", length = length(symbols))
    names(data) <- symbols
    for (symbol in symbols) {
      # load the data for the symbol into singleSymbol.
      singleSymbol <- quantmod::getSymbols.av(symbol, periodicity = "intraday", interval = interval,
                                   api.key = apiKey, auto.assign = FALSE, output.size = "full")
      
      
      if (outputType == "DT") { #Here the user wants a data.table out. Thus we set the class DT and we change the names.
        singleSymbol <- as.data.table(singleSymbol)
        setnames(singleSymbol, colnames(singleSymbol), new = c("DT", "OPEN", "HIGH", "LOW", "CLOSE", "SIZE"))
      } else{
        # Change the names to OHLC-V 
        colnames(singleSymbol) <-  c("OPEN", "HIGH", "LOW", "CLOSE", "SIZE")
      }
      data[[symbol]] <- singleSymbol
      
      if (length(symbols) > 5 && doSleep) {
        print("Sleeping for 12 seconds because the length of symbols > 5.")
        Sys.sleep(12)
      }
    }
  }
  # We have one symbols
  else {
    # load the data for the symbol into singleSymbol.
    data <- quantmod::getSymbols.av(symbols, periodicity = "intraday", interval = interval,
                         api.key = apiKey, auto.assign = FALSE, output.size = "full")  
    if (outputType == "DT") { #Here the user wants a data table out. Thus we set the class DT and we change the names.
      data <- as.data.table(data)
      setnames(data, colnames(data), new = c("DT", "OPEN", "HIGH", "LOW", "CLOSE", "SIZE"))
    } else {
      # Change the names to OHLC-V
      colnames(data) <- c("OPEN", "HIGH", "LOW", "CLOSE", "SIZE")
    }
  }
  
  return(data)
  
}


