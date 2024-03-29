# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

Get_Power_UMP <- function(tau_pool, eg_values, CriVals) {
    .Call(`_MORST_Get_Power_UMP`, tau_pool, eg_values, CriVals)
}

Power_fixed_tau <- function(tau_pool, eg_values, tau0, CV_tau0) {
    .Call(`_MORST_Power_fixed_tau`, tau_pool, eg_values, tau0, CV_tau0)
}

davies_adaptive <- function(Q, w, is_small = FALSE) {
    .Call(`_MORST_davies_adaptive`, Q, w, is_small)
}

Get_CriVals <- function(tau_pool, eg_values, alpha) {
    .Call(`_MORST_Get_CriVals`, tau_pool, eg_values, alpha)
}

Get_CriVals <- function(tau_pool, eg_values, alpha) {
    .Call(`_MORST_Get_CriVals`, tau_pool, eg_values, alpha)
}

Davies <- function(q, lambda, lim = 10000L, acc = 0.0001) {
    .Call(`_MORST_Davies`, q, lambda, lim, acc)
}

which_abs_diff_min <- function(Power_UMP, power_x) {
    .Call(`_MORST_which_abs_diff_min`, Power_UMP, power_x)
}

Minimax_tau <- function(eg_values, alpha, n_points = 50L) {
    .Call(`_MORST_Minimax_tau`, eg_values, alpha, n_points)
}

Saddle <- function(q, egvalues) {
    .Call(`_MORST_Saddle`, q, egvalues)
}

tau_ump <- function(Eg_values, target_power, siglevel, acc_power = 0.01, acc_siglevel = 1e-06) {
    .Call(`_MORST_tau_ump`, Eg_values, target_power, siglevel, acc_power, acc_siglevel)
}

