/*
 --------------------------------------------------
 --------------------------------------------------

 Quantitative genetics for consumer-resource model

 --------------------------------------------------
 --------------------------------------------------
 */

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include <progress.hpp>
#include <random>
#include "pcg.hpp"
#include "quant_gen_cr.hpp"


#ifdef _OPENMP
#include <omp.h>  // omp
#endif

using namespace Rcpp;





//' Search for unique species in a matrix of species trait values.
//'
//' @noRd
//'
arma::uvec unq_spp(const arma::mat& traits,
                   const arma::vec& abund,
                   const double& precision,
                   const double& density_threshold) {

    double precision_ = precision * precision;

    uint32_t n = traits.n_rows;
    arma::uvec is_unq(n, arma::fill::ones);
    for (uint32_t i = 1; i < n; i++) {
        for (uint32_t j = 0; j < i; j++) {
            if (is_unq(j) == 0) continue; // don't want to keep looking at non-unique spp
            arma::rowvec diff_;
            diff_ = traits.row(i) - traits.row(j);
            diff_ %= diff_;
            if (arma::mean(diff_) < precision_) {
                is_unq(i) = 0;
                break;
            }
        }
    }

    arma::uvec unqs = arma::find(is_unq == 1 && abund > density_threshold);

    return unqs;
}


//' Strength of selection for traits `V`, ones for resources
//'
//'
//' @noRd
//'
inline void set_SS_V_(arma::mat& SS_V,
                      const arma::mat& V,
                      const arma::mat& U,
                      const arma::vec& N,
                      const arma::vec& P,
                      const arma::mat& C,
                      const double& r,
                      const double& f,
                      const double& b,
                      const uint32_t& p,
                      const uint32_t& q) {

    arma::mat tmp1(p, q);
    for (uint32_t i = 0; i < q; i++) tmp1.col(i) = P;

    SS_V = (-2 * r * f * arma::accu(N) * V * C + 2 * r * b * V %
        (arma::exp(-(V % V) * (U.t() % U.t())) * ((U % U) % tmp1)));

    return;

}

//' Sets matrix of traits for each resource (<# resources> x <# resources' traits>)
//'
//' @noRd
//'
inline void set_Vt_(arma::mat& Vt,
                    const arma::mat& V,
                    const arma::mat& U,
                    const arma::vec& N,
                    const arma::vec& P,
                    const arma::mat& C,
                    const double& r,
                    const double& f,
                    const double& b,
                    const double& sig2N,
                    const uint32_t& p,
                    const uint32_t& q) {

    set_SS_V_(Vt, V, U, N, P, C, r, f, b, p, q);
    Vt *= sig2N;
    Vt += V;

    return;
}





//' Strength of selection for traits `U`, ones for consumers
//'
//'
//' @noRd
//'
inline void set_SS_U_(arma::mat& SS_U,
                      const arma::mat& V,
                      const arma::mat& U,
                      const arma::vec& N,
                      const arma::mat& D,
                      const double& b,
                      const double& c,
                      const double& g,
                      const uint32_t& n,
                      const uint32_t& p,
                      const uint32_t& q) {

    arma::mat tmp1(n, q);
    arma::mat tmp2(p, q), tmp3;
    for (uint32_t i = 0; i < q; i++) tmp1.col(i) = N;
    tmp3 = arma::diagvec(U * D * U.t());
    tmp3 %= tmp3;
    tmp3 = 1 / tmp3;
    for (uint32_t i = 0; i < q; i++) tmp2.col(i) = tmp3;

    SS_U = (-2 * b * c * U) %
        (arma::exp(-(U % U) * (V.t() % V.t())) * (V % V % tmp1)) +
        ((2 * g * tmp2) % (U * D));

    return;
}

//' Sets traits for each consumer (<# consumers> x <# consumers' traits>)
//'
//'
//' @noRd
//'
inline void set_Ut_(arma::mat& Ut,
                    const arma::mat& V,
                    const arma::mat& U,
                    const arma::vec& N,
                    const arma::mat& D,
                    const double& b,
                    const double& c,
                    const double& g,
                    const double& sig2P,
                    const uint32_t& n,
                    const uint32_t& p,
                    const uint32_t& q) {

    set_SS_U_(Ut, V, U, N, D, b, c, g, n, p, q);
    Ut *= sig2P;
    Ut += U;

    return;
}



//' One trial of iterations with no mutations to new species.
//'
//'
//'
//' @noRd
//'
void one_trial(TrialInfo& out_resource,
               TrialInfo& out_consumer,
               const uint32_t& rep,
               const arma::vec& N0,
               const arma::vec& P0,
               const arma::mat& V0,
               const arma::mat& U0,
               const double& r,
               const double& a,
               const double& f,
               const double& b,
               const double& c,
               const double& m,
               const double& g,
               const double& etaN,
               const double& etaP,
               const double& sig2N,
               const double& sig2P,
               const double& delta,
               const double& delta2,
               const double& start_t,
               const double& max_t,
               pcg32& eng,
               const double& density_threshold,
               const double& precision,
               const uint32_t& return_every) {

    out_resource.rep = rep;
    out_consumer.rep = rep;

    arma::vec N = N0;
    arma::vec P = P0;
    arma::mat V = V0;
    arma::mat U = U0;

    std::lognormal_distribution<double> distr(0.0, delta);

    uint32_t n = N.n_elem;  // # resources
    uint32_t p = P.n_elem;  // # consumers
    uint32_t q = V.n_cols; // # traits

    if (return_every > 0) {

        // Total number of iteration where you'll capture output:
        double n_out_dbl = std::ceil(static_cast<double>(max_t) /
                                      static_cast<double>(return_every));
        uint32_t n_out(n_out_dbl);

        out_resource.abund_t.set_size(n, n_out);
        out_resource.traits_t.set_size(n, q, n_out);

        out_consumer.abund_t.set_size(p, n_out);
        out_consumer.traits_t.set_size(p, q, n_out);

    }


    // Matrices for non-additive effects of traits
    arma::mat C(q, q);
    C.fill(etaN);
    C.diag().fill(1);
    arma::mat D(q, q);
    D.fill(etaP);
    D.diag().fill(1);

    arma::vec A, M, Nt, Pt;
    arma::mat B, Vt, Ut;

    for (uint32_t t = 0; t < start_t; t++) {

        // DD for resources (dimensions: <# resources> x 1):
        A = a + f * arma::diagvec(V * C * V.t());
        // attack rates of consumer-->resources (<# resources> x <# consumers>):
        B = b * arma::exp(-(V % V) * (U.t() % U.t()));
        // consumer mortality rate  (<# consumers> x 1):
        M = m + g / arma::diagvec(U * D * U.t());
        // New resources' abundances (<# resources> x 1)
        Nt = N % arma::exp(r * (1 - A * arma::accu(N) - B * P));
        // New consumers' abundances (<# consumers> x 1)
        Pt = P % arma::exp(c * B.t() * N - M);

        // Traits for each resource (<# resources> x <# resources' traits>):
        set_Vt_(Vt, V, U, N, P, C, r, f, b, sig2N, p, q);
        // Traits for each consumer (<# consumers> x <# consumers' traits>)
        set_Ut_(Ut, V, U, N, D, b, c, g, sig2P, n, p, q);

        // Setting new values:
        N = Nt;
        P = Pt;
        V = arma::abs(Vt);
        U = arma::abs(Ut);

    }

    // perturbation
    for (uint32_t j = 0; j < q; j++) {
        for (uint32_t i = 0; i < n; i++) {
            V(i, j) *= distr(eng);
            // V(i, j) *= R::rlnorm(0, delta);
        }
        for (uint32_t i = 0; i < q; i++) {
            U(i, j) *= distr(eng);
            // U(i, j) *= R::rlnorm(0, delta);
        }
    }
    // Not sure why, but it now gets changed to a different SD:
    distr.param(std::lognormal_distribution<double>::param_type(0.0, delta2));

    for (uint32_t t = 0; t < max_t; t++) {

        if (t == (max_t / 2)) {
            for (uint32_t j = 0; j < q; j++) {
                for (uint32_t i = 0; i < n; i++) {
                    V(i, j) *= distr(eng);
                    // V(i, j) *= R::rlnorm(0, delta2);
                }
                for (uint32_t i = 0; i < q; i++) {
                    U(i, j) *= distr(eng);
                    // U(i, j) *= R::rlnorm(0, delta2);
                }
            }
        }

        A = a + f * arma::diagvec(V * C * V.t());
        B = b * arma::exp(-(V % V) * (U.t() % U.t()));
        M = m + g / arma::diagvec(U * D * U.t());
        Nt = N % arma::exp(r * (1 - A * arma::accu(N) - B * P));
        Pt = P % arma::exp(c * B.t() * N - M);

        set_Vt_(Vt, V, U, N, P, C, r, f, b, sig2N, p, q);
        set_Ut_(Ut, V, U, N, D, b, c, g, sig2P, n, p, q);


        // Setting new values:
        N = Nt;
        P = Pt;
        V = arma::abs(Vt);
        U = arma::abs(Ut);

        // catch low values of U
        U(arma::find(U < 1e-4)).fill(1e-4);

        if (return_every > 0) {
            if (t % return_every == 0) {
                uint32_t t_ = t / return_every;
                out_resource.abund_t.col(t_) = N;
                out_resource.traits_t.slice(t_) = V;
                out_consumer.abund_t.col(t_) = P;
                out_consumer.traits_t.slice(t_) = U;
            }
        }

    }

    // Calculated final fitnesses and selection pressure to see if we're at equilibrium
    arma::vec WN = arma::exp(r * (1 - A * arma::accu(N) - B * P));
    arma::vec WP = arma::exp(c * B.t() * N - M);
    arma::mat SV, SU;
    set_SS_V_(SV, V, U, N, P, C, r, f, b, p, q);
    set_SS_U_(SU, V, U, N, D, b, c, g, n, p, q);

    double fitness_N = arma::prod(WN(arma::find(N > density_threshold)));
    double fitness_P = arma::prod(WP(arma::find(P > density_threshold)));

    double selection_V = arma::accu(SV.rows(arma::find(N > density_threshold)));
    double selection_U = arma::accu(SU.rows(arma::find(P > density_threshold)));

    /*
     truncate any values of V or U that are over 3 to 3 so that unique species
     can be identified
    */
    V(arma::find(V > 3)).fill(3);
    U(arma::find(U > 3)).fill(3);

    // Not sure why the next line is done for consumers but not resources:
    U.rows(arma::find(P < density_threshold)).fill(100);

    // find unique species.
    arma::uvec unq_res = unq_spp(V, N, precision, density_threshold);
    arma::uvec unq_cons = unq_spp(U, P, precision, density_threshold);

    uint32_t n_res = unq_res.n_elem;
    uint32_t n_cons = unq_cons.n_elem;

    /*
     Fill in output:
     */
    out_resource.fitness = fitness_N;
    out_consumer.fitness = fitness_P;
    out_resource.selection = selection_V;
    out_consumer.selection = selection_U;
    if (n_res > 0) {
        out_resource.traits = V.rows(unq_res);
    } else {
        out_resource.traits = arma::mat(1, V.n_cols);
        out_resource.traits.fill(arma::datum::nan);
    }
    if (n_cons > 0) {
        out_consumer.traits = U.rows(unq_cons);
    } else {
        out_consumer.traits = arma::mat(1, U.n_cols);
        out_consumer.traits.fill(arma::datum::nan);
    }

    return;
}



//' Multiple trials for a given set of parameter values.
//'
//' @inheritParams quantgen_trials
//'
//' @noRd
//'
//[[Rcpp::export]]
List quantgen_trials_(const uint32_t& n_trials,
                      const arma::vec& N0,
                      const arma::vec& P0,
                      const arma::mat& V0,
                      const arma::mat& U0,
                      const double& r,
                      const double& a,
                      const double& f,
                      const double& b,
                      const double& c,
                      const double& m,
                      const double& g,
                      const double& etaN,
                      const double& etaP,
                      const double& sig2N,
                      const double& sig2P,
                      const double& delta,
                      const double& delta2,
                      const double& start_t,
                      const double& max_t,
                      const double& density_threshold,
                      const double& precision,
                      const uint32_t& n_cores,
                      const bool& show_progress,
                      const uint32_t& return_every) {

    std::vector<std::vector<TrialInfo>> infos(2);
    infos[0] = std::vector<TrialInfo>(n_trials);
    infos[1] = std::vector<TrialInfo>(n_trials);

    if (V0.n_cols != U0.n_cols) stop("V0 and U0 must have same # columns");

    // Generate seeds for random number generators (1 RNG per core)
    const std::vector<std::vector<uint64_t>> seeds = mc_seeds(n_cores);

    Progress prog_bar(n_trials, show_progress);

    #ifdef _OPENMP
    #pragma omp parallel default(shared) num_threads(n_cores) if (n_cores > 1)
    {
    #endif

    std::vector<uint64_t> active_seeds;

    // Write the active seed per core or just write one of the seeds.
    #ifdef _OPENMP
    uint32_t active_thread = omp_get_thread_num();
    active_seeds = seeds[active_thread];
    #else
    active_seeds = seeds[0];
    #endif

    pcg32 eng = seeded_pcg(active_seeds);

    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (uint32_t i = 0; i < n_trials; i++) {
        if (!Progress::check_abort()) {
            TrialInfo& res_(infos[0][i]);
            TrialInfo& cons_(infos[1][i]);
            one_trial(res_, cons_, i, N0, P0, V0, U0,
                      r, a, f, b, c, m, g, etaN, etaP, sig2N, sig2P, delta, delta2,
                      start_t, max_t, eng, density_threshold, precision,
                      return_every);
            prog_bar.increment();
        }
    }

    #ifdef _OPENMP
    }
    #endif


    /*
     -------------
     Now assembling output:
     -------------
     */

    std::vector<arma::field<arma::mat>> res_cons_mats(2);
    res_cons_mats[0] = arma::field<arma::mat>(n_trials);
    res_cons_mats[1] = arma::field<arma::mat>(n_trials);

    for (uint32_t j = 0; j < res_cons_mats.size(); j++) {
        for (uint32_t i = 0; i < n_trials; i++) {
            // Objects to interact with:
            const TrialInfo& info(infos[j][i]);
            arma::mat& rc_mat(res_cons_mats[j][i]);
            // Adjusting matrix:
            rc_mat.set_size(info.traits.n_rows, info.traits.n_cols + 2);
            rc_mat.col(0).fill(info.fitness);
            rc_mat.col(1).fill(info.selection);
            rc_mat.cols(2, rc_mat.n_cols - 1) = info.traits;
        }
    }

    List output;
    if (return_every > 0) {
        arma::field<arma::mat> res_abund(n_trials);
        arma::field<arma::cube> res_traits(n_trials);
        arma::field<arma::mat> cons_abund(n_trials);
        arma::field<arma::cube> cons_traits(n_trials);
        for (uint32_t i = 0; i < n_trials; i++) {
            res_abund[i] = infos[0][i].abund_t;
            res_traits[i] = infos[0][i].traits_t;
            cons_abund[i] = infos[1][i].abund_t;
            cons_traits[i] = infos[1][i].traits_t;
        }
        output = List::create(
            _["resource"] = List::create(
                _["final"] = res_cons_mats[0],
                _["past_traits"] = res_abund,
                _["past_abund"] = res_traits
                ),
            _["consumer"] = List::create(
                _["final"] = res_cons_mats[1],
                _["past_traits"] = cons_abund,
                _["past_abund"] = cons_traits
                )
        );
    } else {
        output = List::create(
            _["resource"] = res_cons_mats[0],
            _["consumer"] = res_cons_mats[1]
        );
    }

    return output;

}
