#include <RcppArmadillo.h>
#include "internals.h"
using namespace arma;
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]] // So we can use lambda functions

constexpr double u_0 = std::log(1.5);
// s_exp function for the Huang Tauchecn simulation.
inline arma::rowvec s_exp(const arma::rowvec& x){
  arma::rowvec y = arma::zeros<rowvec>(x.n_elem);
  for(arma::uword i = 0; i < y.n_elem; i++){
    if( (x[i] <= u_0) ){
      y[i] = exp(x[i]);
    } else {
      y[i] = exp(u_0)/sqrt(u_0) * sqrt(u_0 - pow(u_0, 2.0) + x[i]);
    }
  }


  return y;
}


// 
// inline arma::rowvec vecpow(arma::rowvec x, const arma::rowvec power){
//   assert (x.n_elem == power.n_elem);
// 
//   for(arma::uword i = 0; i < x.n_elem; i++){
//     if(power[i] == 0.0){
//       x[i] = 1.0;
//     } else {
//       x[i] = pow(x[i], power[i]);
//     }
//   }
// 
//   return(x);
// 
// }



// [[Rcpp::export]]
arma::mat vasicekModel(List model, int nObs, int nSeries, int nDays, const arma::mat& dt){

  const arma::rowvec a = model["meanReversion"];
  const arma::rowvec b = model["drift"];
  arma::rowvec sigma = model["driftVol"];

  sigma = sqrt(sigma);

  arma::mat drift = arma::zeros<mat>((nObs * nDays), nSeries);

  drift.row(0) = b; //+ 1e-8; // So we don't have problems when drift = 0;

  for(int t = 1; t < nObs * nDays; t++){
    drift.row(t) = drift.row(t-1) + a % (b - drift.row(t-1)) % dt.row(t-1) + sigma % randn(1, nSeries) % sqrt(dt.row(t-1));
  }
  return drift;

}



//' Heston model:
//' dXt = mu_t * dt + sigma_t * dW_t
//' dSigma_t^2 = kappa * (theta - sigma_t^2) * dt + xi * sigma_t *dB_t
//' Where Wt and Bt are (often negatively) correlated brownian motions with correlation rho
//' We will not simulate mu_t here, that is done in R code elsewhere.
//' kappa is the mean reversion factor.
//' theta is the long term mean of volatility
//' xi is the vol of vol parameter.
//' @keywords internal
// [[Rcpp::export]]
List hestonModel(List model, int nObs, int nSeries, int nDays, const arma::mat& dt){
  arma::mat sigma = model["sigma"];
  const arma::rowvec kappa = model["meanReversion"];
  const arma::rowvec theta = sigma.diag().t();
  const arma::rowvec xi = model["volOfVol"];
  const arma::rowvec rho = model["rho"];

  // Map the variance covariance matrix into a correlation matrix to use for RNG
  sigma = sqrt(diagmat(sigma.diag())).i() * sigma * sqrt(diagmat(sigma.diag())).i();

  arma::mat wt = arma::mvnrnd(arma::zeros(nSeries), sigma, nObs * nDays).t() % sqrt(dt);
  arma::mat bt = arma::randn((nObs * nDays), nSeries) % sqrt(dt);
  // Create the bt matrix as correlated with wt using lambda functions.
  bt = bt.each_row(([&rho](const arma::rowvec& B){return B % sqrt(1.0 - square(rho)); })) + wt.each_row(([&rho](const arma::rowvec& W){ return W % rho;}));



  // Heston model in terms of dSigma^2 = ... This caused problems with NaNs (dunno why)
  // sigma^2 container
  // arma::mat sigma2 = mat(bt);
  //
  // for(int t = 1; t < nObs * nDays; t++){
  //   sigma2.row(t) = sigma2.row(t-1) + kappa % (theta - sigma2.row(t-1)) % dt.row(t) + xi % sqrt(sigma2.row(t-1)) % bt.row(t);
  // }


  // Heston model in terms of dSigma = ... Seems to not produce NaNs
  arma::mat sigmat = mat(bt);

  // Initialize sigma_t: as a draw from its unconditional distribution.
  for(int i = 0; i<nSeries; i++){
    sigmat(0,i) = as_scalar(randg(1, distr_param(2.0 * kappa(i) * sqrt(theta(i)) * pow(xi(i), -2.0), 1.0/(2.0 * kappa(i) * pow(xi(i), -2.0)))));
  }

  for(int t = 1; t < nObs * nDays; t++){
    sigmat.row(t) = sigmat.row(t-1) + kappa % (theta - square(sigmat.row(t-1))) % dt.row(t) + xi % sigmat.row(t-1) % bt.row(t);
  }

  arma::mat returns = sigmat % wt;
  
  // Create empty list and populate it. Probably use Rcpp::Create::List later
  List out;
  out["sigma"] = sigmat;
  out["returns"] = returns;
  return out;

}



//' Huang and Tauchen 2005 model
//' (we don't include drift here) (sigma_ut is calculated elsewhere too!) [[Drift and diurnality is done elsewhere]]
//' sigma_ut * nu_t (rho_1 * dWt_1 + rho_2 * dWt_2 + sqrt(1 - rho_1 ^2 - rho^2) * dWt_3)
//' nu_t^2 = s-exp(beta_0 + beta_1 * nu_t_1^2 + beta_2 * nu_t_2^2)
//' dnu_1_t^2 = alpha_1  * nu_1_t^2 * dt + dWt_1
//' dnu_2_t^2 = alpha_2 * nu_2_t^2 * dt + (1 + phi nu_2_t^2) * dW2_t
//' sigma_ut = C + A*exp(-a*t) + B*exp(-b*(1-t))
//' @keywords internal
// [[Rcpp::export]]
List huangTauchen(List model, int nObs, int nSeries, int nDays, const arma::mat& dt){

  // unpack values
  const arma::rowvec alpha1 = model["alpha1"];
  const arma::rowvec alpha2 = model["alpha2"];
  const arma::rowvec beta0 =  model["beta0"];
  const arma::rowvec beta1 = model["beta1"];
  const arma::rowvec beta2 = model["beta2"];
  const arma::rowvec phi = model["phi"];
  const arma::rowvec rho1 = model["rho1"];
  const arma::rowvec rho2 = model["rho2"];

  // Initialize containers
  const arma::mat dW1t = arma::randn(nObs * nDays, nSeries) % sqrt(dt);
  const arma::mat dW2t = arma::randn(nObs * nDays, nSeries) % sqrt(dt);
  const arma::mat dW3t = arma::randn(nObs * nDays, nSeries) % sqrt(dt);

  arma::mat volFactor1 = arma::zeros<mat>((nObs * nDays), nSeries);
  arma::mat volFactor2 = arma::zeros<mat>((nObs * nDays), nSeries);
  arma::mat nu = arma::zeros<mat>((nObs * nDays), nSeries);
  //volFactor1.row(0) = arma::randn(1, nSeries) % sqrt(-1.0/(2 * alpha1)) + dW1t.row(0);

  nu.row(0) = s_exp(beta0 + beta1 % volFactor1.row(0) + beta2 % volFactor2.row(0));
  for(int t = 1; t < nObs * nDays; t++){

    volFactor1.row(t) = volFactor1.row(t-1) + alpha1 % volFactor1.row(t-1) % dt.row(t) + dW1t.row(t);

    volFactor2.row(t) = volFactor2.row(t-1) + alpha2 % volFactor2.row(t-1) % dt.row(t) + (1.0 + phi % volFactor2.row(t-1)) % dW2t.row(t);

    nu.row(t) = s_exp(beta0 + beta1 % volFactor1.row(t) + beta2 % volFactor2.row(t));

  }


  arma::mat returns = sqrt(nu) % (dW1t.each_row([&rho1](const arma::rowvec& w1){return rho1 % w1;}) + dW2t.each_row([&rho2](const arma::rowvec& w2){return rho2 % w2;}) +
    dW3t.each_row([&rho1, &rho2](const arma::rowvec& w3){ sqrt(1.0 - rho1 - rho2) % w3;}));
  
  // Create empty list and populate it. Probably use Rcpp::Create::List later
  List out;
  out["returns"] = returns;
  out["sigma"] = sqrt(nu);
  out["volatilityFactor"] = volFactor1;
  out["volatilityFactor2"] = volFactor2;
  return out;
}


// [[Rcpp::export]]
List liLinton(List model, const int nObs, const int nSeries, const int nDays, const arma::mat& dt){
  // Parameters:
  // kappa_1, mu_1, kappa_2, mu_2, eta_1, rho, lambda, delta
  const arma::rowvec kappa1 = model["kappa1"];
  const arma::rowvec mu1 = model["mu1"];
  const arma::rowvec kappa2 =  model["meanReversion"];
  const arma::rowvec mu2 = model["sigma"];
  const arma::rowvec eta = model["volOfVol"];
  const arma::rowvec rho = model["rho"];
  const arma::rowvec lambda = model["lambda"]; 
  const arma::rowvec delta = model["volJumps"]; 
  
  
  arma::mat Nt = arma::mat(nObs * nDays, nSeries);
  arma::mat volJumps = arma::mat(nObs * nDays, nSeries);
  arma::mat jumps = arma::mat(nObs * nDays, nSeries);
  // Calculate volatililty jumps and price jumps
  for(int j = 0; j < nSeries; j++){
    Nt.col(j) = arma::colvec(Rcpp::rpois(nObs * nDays, 1.0/lambda(j))) % dt.col(j);
    volJumps.col(j) = arma::colvec(Rcpp::rexp(nObs * nDays, 1.0/delta(j))) % Nt.col(j);
    jumps.col(j) = arma::randn<colvec>(nObs * nDays) % Nt.col(j) * sqrt(mu2(j)/10.0);
  }
  
  // We have corr(w1t, w2t) = rho -> simulate independently and transform w2t to be correlated with w1t.
  arma::mat w1t = arma::randn(nObs * nDays, nSeries) % sqrt(dt);
  arma::mat w2t = arma::randn((nObs * nDays), nSeries) % sqrt(dt);
  // Create the w2t matrix as correlated with w1t using lambda functions.
  w2t = w2t.each_row(([&rho](const arma::rowvec& w2){return w2 % sqrt(1.0 - square(rho)); })) + w1t.each_row(([&rho](const arma::rowvec& w1){return w1 % rho;}));
  
  arma::mat prices = arma::mat(nObs * nDays, nSeries);
  // arma::mat sigma = arma::mat(nObs * nDays, nSeries);
  arma::mat sigma2 = arma::mat(nObs * nDays, nSeries);
  // sigma.row(0) = sqrt(mu2 + volJumps.row(0));
  sigma2.row(0) = mu2 + volJumps.row(0);
  prices.row(0) = mu1 + jumps.row(0);
  for(int t = 1; t < nObs * nDays; t++){
    // sigma.row(t) = sigma.row(t-1) + kappa2 % (sqrt(mu2) - square(sigma.row(t-1))) % dt.row(t) + eta % sigma.row(t) % w2t.row(t) + volJumps.row(t);
    sigma2.row(t) = sigma2.row(t-1) + kappa2 % (mu2 - sigma2.row(t-1)) % dt.row(t) + eta % sqrt(sigma2.row(t-1)) % w2t.row(t) + volJumps.row(t);
    
    // prices.row(t) = prices.row(t-1) + kappa1 % (mu1 - prices.row(t-1)) % dt.row(t) + sigma.row(t) % w1t.row(t) + jumps.row(t);
    prices.row(t) = prices.row(t-1) + kappa1 % (mu1 - prices.row(t-1)) % dt.row(t) + sqrt(sigma2.row(t)) % w1t.row(t) + jumps.row(t);
  }
  
  
  
  // Create empty list and populate it. Probably use Rcpp::Create::List later
  List out;
  out["jumps"] = jumps;
  out["volJumps"] = volJumps;
  out["prices"] = prices;
  arma::mat returns = diff(prices);
  returns.insert_rows(0, arma::zeros<rowvec>(nSeries));
  out["returns"] = returns;
  // out["sigma"] = sigma;
  out["sigma"] = sqrt(sigma2);
  out["jumpIntensities"] = Nt;
  return out;
}


