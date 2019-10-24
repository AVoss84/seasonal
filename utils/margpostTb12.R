function (data, p = 1, deter = NULL, sdeter = "const", isbreak.nsea = c(1, 
    0), Tb.fix = NULL, sbreak = T, normc = T, isbreak.sea = c(1, 
    0), coding = c("mean", "effect"), tau = 11) 
{
    bayes.par12 = dget("bayes.par12.R")
    N = length(data)
    stopifnot(is.ts(data))
    norm = NULL
    lpdensity = pdensity = hpd95 = NULL
    (Tb = if (!is.null(Tb.fix)) {
        Tb.fix
    }
    else {
        seq(tau + 1 + p, N - tau - 1)
    })
    if (sbreak) {
        crit = NULL
        ii = 1
        for (tt in Tb) {
            (out = bayes.par12(data = data, p = p, deter = deter, 
                Tb = Tb[ii], sdeter = sdeter, coding = coding[1], 
                isbreak.nsea = isbreak.nsea, isbreak.sea = isbreak.sea, 
                tau = tau))
            (W = out@W)
            (lhs = out@lhs)
            crit = rbind(crit, c(bic = out@bic, aic = out@aic, 
                hq = out@hq, fpe = out@fpe))
            d = ncol(W)
            (k = N - p - d)
            nu = k
            (ols = out@est.classic)
            QF = (out@sigma2)
            (Det = det(t(W) %*% W))
            if (Det == 0) {
                Det = 0.001
                warning("Singularity encountered. Determinant equal to zero.\n")
            }
            log_kernel = as.numeric(-0.5 * out@bic)
            lpdensity = c(lpdensity, log_kernel)
            print(tt)
            ii = ii + 1
        }
        if (normc && is.null(Tb.fix)) {
            (pdensity = exp(lpdensity - max(lpdensity)))
            norm = sum(pdensity)
            (pdensity = pdensity/norm)
            (cdf = cumsum(pdensity))
            (q0025 = Tb[min(which(cdf >= 0.025))])
            (q0975 = Tb[min(which(cdf >= 0.975))])
            (marg = mean(pdensity))
            (hpd95 = c(q0025, q0975))
        }
        if (!normc && is.null(Tb.fix)) {
            (pdensity = exp(lpdensity))
            (marg = mean(pdensity))
            message("No HPD region available if not normalized/proper density is used.\n")
            hpd95 = c(NA, NA)
        }
        if (!normc && !is.null(Tb.fix)) {
            marg = exp(lpdensity)
        }
        if (normc && !is.null(Tb.fix)) {
            stop("Only one density point available!\n ")
        }
    }
    else {
        (out = bayes.par12(data = data, p = p, deter = deter, 
            Tb = NULL, sdeter = sdeter, coding = coding[1], isbreak.nsea = isbreak.nsea, 
            isbreak.sea = isbreak.sea, tau = tau))
        (W = out@W)
        (lhs = out@lhs)
        crit = c(bic = out@bic, aic = out@aic, hq = out@hq, fpe = out@fpe)
        d = ncol(W)
        (k = N - p - d)
        (nu = k)
        (ols = out@est.classic)
        QF = (out@sigma2)
        (Det = det(t(W) %*% W))
        if (Det == 0) {
            Det = 0.001
            warning("Singularity encountered. Determinant equal to zero.\n")
        }
        (marg = (1/sqrt(Det)) * (QF)^(-k/2))
        marg = as.numeric(exp(-0.5 * out@bic))
        message("Laplace approximation to the marginal likelihood was used!\n")
        (lpdensity = -0.5 * log(Det) - 0.5 * k * log(QF))
        pdensity = hpd95 = Tb = NULL
    }
    return(list(ols = ols, norm = norm, Tb = Tb, hpd95 = hpd95, 
        crit = crit, lpdensity = lpdensity, pdensity = pdensity, 
        marg = marg, info = list(QF = QF)))
}
