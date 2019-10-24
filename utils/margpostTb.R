function (data, p = 1, deter = NULL, sdeter = "const", isbreak.nsea = c(1, 
    0), Tb.fix = NULL, sbreak = T, normc = T, isbreak.sea = c(1, 
    0), coding = c("mean", "effect"), tau = 3) 
{
    bayes.par = dget("bayes.par.R")
    library(Bolstad)
    N = length(data)
    stopifnot(is.ts(data))
    ii = 1
    lpdensity = hpd95 = NULL
    (Tb = if (!is.null(Tb.fix)) {
        Tb.fix
    }
    else {
        seq(tau + 1 + p, N - tau - 1)
    })
    if (sbreak) {
        for (tt in Tb) {
            (out = bayes.par(data = data, p = p, deter = deter, 
                Tb = tt, sdeter = sdeter, coding = coding[1], 
                isbreak.nsea = isbreak.nsea, isbreak.sea = isbreak.sea, 
                tau = tau))
            (W = out@W)
            (lhs = out@lhs)
            bic = out@bic
            aic = out@aic
            hq = out@hq
            fpe = out@fpe
            d = ncol(W)
            (k = N - p - 1 - d)
            nu = k + 1
            (ols = out@est.classic)
            (QF = (out@sigma2) * nu)
            (Det = det(t(W) %*% W))
            if (Det == 0) {
                Det = 1e-04
                warning("Singularity encountered. Determinant equal to zero.\n")
            }
            (log_kernel = -0.5 * log(Det) - 0.5 * k * log(QF))
            lpdensity = c(lpdensity, log_kernel)
            print(ii)
            ii = ii + 1
        }
        if (normc && is.null(Tb.fix)) {
            norm = sintegral(Tb, exp(lpdensity - max(lpdensity)), 
                n.pts = length(Tb) - 10)
            (pdensity = exp(lpdensity - max(lpdensity))/norm$value)
            (marg = mean(pdensity))
            (cdf = sintegral(Tb, pdensity)$cdf)
            (q0025 = cdf$x[min(which(cdf$y >= 0.025))])
            (q0975 = cdf$x[min(which(cdf$y >= 0.975))])
            (hpd95 = c(q0025, q0975))
        }
        if (!normc && is.null(Tb.fix)) {
            (pdensity = exp(lpdensity))
            (marg = mean(pdensity))
            message("No HPD region available if not normalized/proper density is used.\n")
            hpd95 = c(NA, NA)
        }
        if (!normc && !is.null(Tb.fix)) {
            marg = pdensity = exp(lpdensity)
        }
        if (normc && !is.null(Tb.fix)) {
            stop("Only one density point available!\n ")
        }
    }
    else {
        (out = bayes.par(data = data, p = p, deter = deter, Tb = NULL, 
            sdeter = sdeter, coding = coding[1], isbreak.nsea = isbreak.nsea, 
            isbreak.sea = isbreak.sea, tau = tau))
        (W = out@W)
        (lhs = out@lhs)
        bic = out@bic
        aic = out@aic
        hq = out@hq
        fpe = out@fpe
        d = ncol(W)
        (k = N - p - 1 - d)
        (nu = k + 1)
        (ols = out@est.classic)
        (QF = (out@sigma2) * nu)
        (Det = det(t(W) %*% W))
        if (Det == 0) {
            Det = 1e-04
            warning("Singularity encountered. Determinant equal to zero.\n")
        }
        (marg = (1/sqrt(Det)) * (QF)^(-k/2))
        if (marg == 0 || marg == Inf) {
            marg = as.numeric(exp(-0.5 * out@bic))
            message("Laplace approximation to the marginal likelihood was used!\n")
        }
        lpdensity = pdensity = hpd95 = Tb = NULL
    }
    return(list(Tb = Tb, hpd95 = hpd95, lpdensity = lpdensity, 
        pdensity = pdensity, marg = marg, info = list(bic = bic, 
            aic = aic, hq = hq, fpe = fpe)))
}
