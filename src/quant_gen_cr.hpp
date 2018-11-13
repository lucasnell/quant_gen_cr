#ifndef __QUANTGENCR_H_
#define __QUANTGENCR_H_

#include <RcppArmadillo.h>

using namespace Rcpp;



struct TrialInfo {

    uint32_t rep;
    arma::mat traits;
    double fitness;
    double selection;
    arma::mat abund_t;
    arma::cube traits_t;

};


#endif
