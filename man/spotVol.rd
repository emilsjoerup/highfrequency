\name{spotvol}
\alias{spotvol}
\title{Spot volatility estimation}
\usage{
spotvol(data, method = "detper", ..., on = "minutes", k = 5,
  marketopen = "09:30:00", marketclose = "16:00:00", tz = "GMT")
}
\arguments{
\item{data}{Either an \code{xts} object, containing price data, or a
\code{matrix} containing returns. For price data, irregularly spaced
observations are allowed. They will be aggregated to the level specified by
parameters \code{on} and \code{k}. For return data, the observations are
assumed to be equispaced, with the time between them specified by \code{on}
and \code{k}. Return data should be in matrix form, where each row
corresponds to a day, and each column to an intraday period. The output
will be in the same form as the input (\code{xts} or
\code{matrix}/\code{numeric}).} 

\item{method}{specifies which method will be used to estimate the spot
volatility. Options include \code{"detper"} and \code{"stochper"}.
See 'Details'.}

\item{on}{string indicating the time scale in which \code{k} is expressed.
Possible values are: \code{"secs", "seconds", "mins", "minutes", "hours"}.}

\item{k}{positive integer, indicating the number of periods to aggregate
over. E.g. to aggregate an \code{xts} object to the 5 minute frequency, set
\code{k = 5} and \code{on = "minutes"}.}

\item{marketopen}{the market opening time. This should be in the time zone
specified by \code{tz}. By default, \code{marketopen = "09:30:00"}.}

\item{marketclose}{the market closing time. This should be in the time zone
specified by \code{tz}. By default, \code{marketclose = "16:00:00"}.}

\item{tz}{string specifying the time zone to which the times in \code{data}
and/or \code{marketopen}/ \code{marketclose} belong. Default = \code{"GMT"}.}

\item{...}{method-specific parameters (see 'Details').}
}
\value{
A \code{spotvol} object, which is a list containing one or more of the
following outputs, depending on the method used:

\code{spot}

An \code{xts} or \code{matrix} object (depending on the input) containing
spot volatility estimates \eqn{\sigma_{t,i}}, reported for each interval
\eqn{i} between \code{marketopen} and \code{marketclose} for every day
\eqn{t} in \code{data}. The length of the intervals is specifiedby \code{k}
and \code{on}. Methods that provide this output: All.

\code{daily}

An \code{xts} or \code{numeric} object (depending on the input) containing
estimates of the daily volatility levels for each day \eqn{t} in \code{data},
if the used method decomposed spot volatility into a daily and an intraday
component. Methods that provide this output: \code{"detper"}.

\code{periodic}

An \code{xts} or \code{numeric} object (depending on the input) containing
estimates of the intraday periodicity factor for each day interval \eqn{i}
between \code{marketopen} and \code{marketclose}, if the spot volatility was
decomposed into a daily and an intraday component. If the output is in
\code{xts} format, this periodicity factor will be dated to the first day of
the input data, but it is identical for each day in the sample. Methods that
provide this output: \code{"detper"}.

\code{par}

A named list containing parameter estimates, for methods that estimate one
or more parameters. Methods that provide this output:
\code{"stochper", "kernel"}.

\code{cp}

A vector containing the change points in the volatility, i.e. the observation
indices after which the volatility level changed, according to the applied
tests. The vector starts with a 0. Methods that provide this output:
\code{"piecewise"}.

\code{ugarchfit}

A \code{ugarchfit} object, as used by the \code{rugarch}
package, containing all output from fitting the GARCH model to the data.
Methods that provide this output: \code{"garch"}.
}
\description{
The \code{spotvol} function offers several methods to estimate spot
volatility and its intraday seasonality, using high-frequency data. It
returns an object of class \code{spotvol}, which can contain various outputs,
depending on the method used. See 'Details' for a description of each method.
In any case, the output will contain the spot volatility estimates.

The input can consist of price data or return data, either tick by tick or
sampled at set intervals. The data will be converted to equispaced
high-frequency returns \eqn{r_{t,i}} (read: the \eqn{i}th return on day
\eqn{t}).
}
\details{
The following estimation methods can be specified in \code{method}:

\strong{Deterministic periodicity method (\code{"detper"})}

Parameters:
\tabular{ll}{
\code{dailyvol} \tab A string specifying the estimation method for the daily
component \eqn{s_t}. Possible values are \code{"bipower", "rv", "medrv"}.
Default = \code{"bipower"}. \cr
\code{periodicvol} \tab A string specifying the estimation method for the
component of intraday volatility, that depends in a deterministic way on the
intraday time at which the return is observed. Possible values are
\code{"TML", "SD", "WSD", "OLS"}. See Boudt et al. (2011) for details.
Default = \code{"TML"}.\cr
\code{P1} \tab A positive integer corresponding to the number of cosinus
terms used in the flexible Fourier specification of the periodicity function,
see Andersen et al. (1997) for details. Default = 5. \cr
\code{P2} \tab Same as \code{P1}, but for the sinus terms. Default = 5.\cr
\code{dummies} \tab Boolean: in case it is \code{TRUE}, the parametric
estimator of periodic standard deviation specifies the periodicity function
as the sum of dummy variables corresponding to each intraday period. If it
is \code{FALSE}, the parametric estimator uses the flexible Fourier
specification. Default = \code{FALSE}.
}
Outputs (see 'Value' for a full description of each component):
\itemize{
\item{\code{spot}}
\item{\code{daily}}
\item{\code{periodic}}
}
The spot volatility is decomposed into a deterministic periodic factor
\eqn{f_{i}} (identical for every day in the sample) and a daily factor
\eqn{s_{t}} (identical for all observations within a day). Both components
are then estimated separately. For more details, see Taylor and Xu (1997)
and Andersen and Bollerslev (1997). The jump robust versions by Boudt et al.
(2011) have also been implemented.

\strong{Stochastic periodicity method (\code{"stochper"})}

Parameters:
\tabular{ll}{
\code{P1} \tab A positive integer corresponding to the number of cosinus
terms used in the flexible Fourier specification of the periodicity function.
Default = 5. \cr
\code{P2} \tab Same as \code{P1}, but for the sinus terms. Default = 5.\cr
\code{init} \tab A named list of initial values to be used in the
optimization routine (\code{"BFGS"} in \code{optim}). Default =
\code{list(sigma = 0.03, sigma_mu = 0.005, sigma_h = 0.005, sigma_k = 0.05,
phi = 0.2, rho = 0.98, mu = c(2, -0.5), delta_c = rep(0, max(1,P1)),
delta_s = rep(0, max(1,P2)))}. See Beltratti & Morana (2001) for a definition
of each parameter. \code{init} can contain any number of these parameters.
For parameters not specified in \code{init}, the default initial value will
be used.\cr
\code{control} \tab A list of options to be passed down to
\code{optim}.
}
Outputs (see 'Value' for a full description of each component):
\itemize{
\item{\code{spot}}
\item{\code{par}}
}
This method by Beltratti and Morana (2001) assumes the periodicity factor to
be stochastic. The spot volatility estimation is split into four components:
a random walk, an autoregressive process, a stochastic cyclical process and
a deterministic cyclical process. The model is estimated using a
quasi-maximum likelihood method based on the Kalman Filter. The package
\code{FKF} is used to apply the Kalman filter. In addition to
the spot volatility estimates, all parameter estimates are returned.

\strong{Nonparametric filtering (\code{"kernel"})}

Parameters:
\tabular{ll}{
\code{type} \tab String specifying the type of kernel to be used. Options
include \code{"gaussian", "epanechnikov", "beta"}. Default =
\code{"gaussian"}.\cr
\code{h} \tab Scalar or vector specifying bandwidth(s) to be used in kernel.
If \code{h} is a scalar, it will be assumed equal throughout the sample. If
it is a vector, it should contain bandwidths for each day. If left empty,
it will be estimated. Default = \code{NULL}. \cr
\code{est} \tab String specifiying the bandwidth estimation method. Possible
values include \code{"cv", "quarticity"}. Method \code{"cv"} equals
cross-validation, which chooses the bandwidth that minimizes the Integrated
Square Error. \code{"quarticity"} multiplies the simple plug-in estimator
by a factor based on the daily quarticity of the returns. \code{est} is
obsolete if \code{h} has already been specified by the user. Default =
\code{"cv"}.\cr
\code{lower} \tab Lower bound to be used in bandwidth optimization routine,
when using cross-validation method. Default is \eqn{0.1n^{-0.2}}. \cr
\code{upper} \tab Upper bound to be used in bandwidth optimization routine,
when using cross-validation method. Default is \eqn{n^{-0.2}}. \cr
}
Outputs (see 'Value' for a full description of each component):
\itemize{
\item{\code{spot}}
\item{\code{par}}
}
This method by Kristensen (2010) filters the spot volatility in a
nonparametric way by applying kernel weights to the standard realized
volatility estimator. Different kernels and bandwidths can
be used to focus on specific characteristics of the volatility process.

Estimation results heavily depend on the bandwidth parameter \eqn{h}, so it
is important that this parameter is well chosen. However, it is difficult to
come up with a method that determines the optimal bandwidth for any kind of
data or kernel that can be used. Although some estimation methods are
provided, it is advised that you specify \eqn{h} yourself, or make sure that
the estimation results are appropiate.

One way to estimate \eqn{h}, is by using cross-validation. For each day in
the sample, \eqn{h} is chosen as to minimize the Integrated Square Error,
which is a function of \eqn{h}. However, this function often has multiple
local minima, or no minima at all (\eqn{h -> \infty}). To ensure a reasonable
optimum is reached, strict boundaries have to be imposed on \eqn{h}. These
can be specified by \code{lower} and \code{upper}, which by default are
\eqn{0.1n^{-0.2}} and \eqn{n^{-0.2}} respectively, where \eqn{n} is the
number of observations in a day.

When using the method \code{"kernel"}, in addition to the spot volatility
estimates, all used values of the bandwidth \eqn{h} are returned.

\strong{Piecewise constant volatility (\code{"piecewise"})}

Parameters:
\tabular{ll}{
\code{type} \tab String specifying the type of test to be used. Options
include \code{"MDa", "MDb", "DM"}. See Fried (2012) for details. Default =
\code{"MDa"}.\cr
\code{m} \tab Number of observations to include in reference window.
Default = \code{40}. \cr
\code{n} \tab Number of observations to include in test window.
Default = \code{20}. \cr
\code{alpha} \tab Significance level to be used in tests. Note that the test
will be executed many times (roughly equal to the total number of
observations), so it is advised to use a small value for \code{alpha}, to
avoid a lot of false positives. Default = \code{0.005}. \cr
\code{volest} \tab String specifying the realized volatility estimator to be
used in local windows. Possible values are \code{"bipower", "rv", "medrv"}.
Default = \code{"bipower"}. \cr
\code{online} \tab Boolean indicating whether estimations at a certain point
\eqn{t} should be done online (using only information available at
\eqn{t-1}), or ex post (using all observations between two change points).
Default = \code{TRUE}.  \cr
}
Outputs (see 'Value' for a full description of each component):
\itemize{
\item{\code{spot}}
\item{\code{cp}}
}

This nonparametric method by Fried (2012) assumes the volatility to be
piecewise constant over local windows. Robust two-sample tests are applied to
detect changes in variability between subsequent windows. The spot volatility
can then be estimated by evaluating regular realized volatility estimators
within each local window.

Along with the spot volatility estimates, this method will return the
detected change points in the volatility level. When plotting a
\code{spotvol} object containing \code{cp}, these change points will be
visualized.

\strong{GARCH models with intraday seasonality  (\code{"garch"})}

Parameters:
\tabular{ll}{
\code{model} \tab String specifying the type of test to be used. Options
include \code{"sGARCH", "eGARCH"}. See \code{ugarchspec} in the
\code{rugarch} package. Default = \code{"eGARCH"}. \cr
\code{garchorder} \tab Numeric value of length 2, containing the order of
the GARCH model to be estimated. Default = \code{c(1,1)}. \cr
\code{dist} \tab String specifying the distribution to be assumed on the
innovations. See \code{distribution.model} in \code{ugarchspec} for
possible options. Default = \code{"norm"}. \cr
\code{solver.control} \tab List containing solver options.
See \code{ugarchfit} for possible values. Default = \code{list()}. \cr
\code{P1} \tab A positive integer corresponding to the number of cosinus
terms used in the flexible Fourier specification of the periodicity function.
Default = 5. \cr
\code{P2} \tab Same as \code{P1}, but for the sinus terms. Default = 5.\cr
}
Outputs (see 'Value' for a full description of each component):
\itemize{
\item{\code{spot}}
\item{\code{ugarchfit}}
}
This method generates the external regressors needed to model the intraday
seasonality with a Flexible Fourier form. The \code{rugarch} package
is then employed to estimate the specified GARCH(1,1) model.

Along with the spot volatility estimates, this method will return the
\code{ugarchfit} object used by the \code{rugarch} package.
}
\examples{
# Load sample data

data(sample_real5minprices)
data(sample_returns_5min)

# Default method, deterministic periodicity
\dontshow{par.def <- par(no.readonly = TRUE)}
vol1 <- spotvol(sample_real5minprices)
plot(vol1)
\dontshow{par(par.def)}
# Compare to stochastic periodicity
\donttest{
init = list(sigma = 0.03, sigma_mu = 0.005, sigma_h = 0.007,
            sigma_k = 0.06, phi = 0.194, rho = 0.986, mu = c(1.87,-0.42),
            delta_c = c(0.25, -0.05, -0.2, 0.13, 0.02), delta_s = c(-1.2,
            0.11, 0.26, -0.03, 0.08))
# next method will take around 110 iterations
vol2 <- spotvol(sample_real5minprices, method = "stochper", init = init)
plot(as.numeric(vol1$spot[1:780]), type="l")
lines(as.numeric(vol2$spot[1:780]), col="red")
legend("topright", c("detper", "stochper"), col = c("black", "red"), lty=1)}

# Various kernel estimates
\donttest{
h1 = bw.nrd0((1:nrow(sample_returns_5min))*(5*60))
vol3 <- spotvol(sample_returns_5min, method = "kernel", h = h1)
vol4 <- spotvol(sample_returns_5min, method = "kernel", est = "quarticity")
vol5 <- spotvol(sample_returns_5min, method = "kernel", est = "cv")
plot(vol3, length = 2880)
lines(as.numeric(t(vol4$spot))[1:2880], col="red")
lines(as.numeric(t(vol5$spot))[1:2880], col="blue")
legend("topright", c("h = simple estimate", "h = quarticity corrected",
       "h = crossvalidated"), col = c("black", "red", "blue"), lty=1)}
\dontshow{par(par.def)}
# Piecewise constant volatility, using an example from Fried (2012)
\donttest{
simdata <- matrix(sqrt(5/3)*rt(3000, df = 5), ncol = 500, byrow = TRUE)
simdata <- c(1, 1, 1.5, 1.5, 2, 1)*simdata
# the volatility of the simulated now changes at 1000, 2000 and 2500
vol6 <- spotvol(simdata, method = "piecewise", m = 200, n  = 100,
                online = FALSE)
plot(vol6)}

# Compare regular GARCH(1,1) model to eGARCH, both with external regressors
\donttest{
vol7 <- spotvol(sample_returns_5min, method = "garch", model = "sGARCH")
vol8 <- spotvol(sample_returns_5min, method = "garch", model = "eGARCH")
plot(as.numeric(t(vol7$spot)), type = "l")
lines(as.numeric(t(vol8$spot)), col = "red")
legend("topleft", c("GARCH", "eGARCH"), col = c("black", "red"), lty=1)}
}
\references{
Andersen, T. G. and T. Bollerslev (1997). Intraday periodicity and volatility
persistence in financial markets. Journal of Empirical Finance 4, 115-158.

Beltratti, A. and C. Morana (2001). Deterministic and stochastic methods for estimation
of intraday seasonal components with high frequency data. Economic Notes 30, 205-234.

Boudt K., Croux C. and Laurent S. (2011). Robust estimation of intraweek periodicity
in volatility and jump detection. Journal of Empirical Finance 18, 353-367.

Fried, Roland (2012). On the online estimation of local constant volatilities.
Computational Statistics and Data Analysis 56, 3080-3090.

Kristensen, Dennis (2010). Nonparametric filtering of the realized spot volatility:
A kernel-based approach. Econometric Theory 26, 60-93.

Taylor, S. J. and X. Xu (1997). The incremental volatility information in one million
foreign exchange quotations. Journal of Empirical Finance 4, 317-340.
}
