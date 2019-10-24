function (data, p = NULL, deter = NULL, sdeter = NULL, coding = c("mean", 
    "effect"), Tb = NULL, class.est = T, tau = 3, isbreak.nsea = c(1, 
    0), isbreak.sea = c(1, 0), trim = F) 
{
    dter.comp = dget("dter.comp.R")
    strbr.comp = dget("strbr.comp.R")
    setClass("bayes.par", representation(Tb = "ANY", dcomps = "ANY", 
        lhs = "ANY", p = "numeric", aic = "numeric", bic = "numeric", 
        hq = "numeric", fpe = "numeric", est.classic = "ANY", 
        sigma2 = "numeric", W = "ANY", LHS = "ANY"))
    setMethod("show", "bayes.par", function(object) {
        if (!is.null(object@Tb)) {
            cat("\n", rep("-", 18), "\n", sprintf(" Gaussian quarterly PAR(%i) model:\n", 
                object@p), rep("-", 18), "\n", sprintf(" Break date T_{b} = %i\n", 
                object@Tb), rep("-", 18), "\n\n")
        }
        else {
            cat("\n", rep("-", 18), "\n", rep("-", 18), "\n", 
                sprintf(" Gaussian quarterly PAR(%i) model:\n", 
                  object@p), rep("-", 18), "\n", rep("-", 18), 
                "\n\n")
        }
        cat(" Classical least squares estimates:\n")
        print(object@est.classic)
        cat("\n Information criteria:\n")
        print(list(aic = object@aic, bic = object@bic, hq = object@hq, 
            fpe = object@fpe))
    })
    if (p < 1) {
        stop("Autoregressive lag order must be at least one!!\n")
    }
    y = data
    stopifnot(is.ts(y))
    S = 4
    N = length(y) + S
    lags = p
    coding = coding[1]
    nseason = deter
    season = sdeter
    (Q1 = which(cycle(y) == 1))
    (Q2 = which(cycle(y) == 2))
    (Q3 = which(cycle(y) == 3))
    (Q4 = which(cycle(y) == 4))
    (le = c(length(Q1), length(Q2), length(Q3), length(Q4)))
    (ranks = rank(le, ties.method = "average"))
    (lasts = c(Q1[length(Q1)], Q2[length(Q2)], Q3[length(Q3)], 
        Q4[length(Q4)]))
    if (sum(le)%%4 != 0 && trim) {
        warning("Unequal number of observations per season. Some quarterly subseries were trimmed.\n")
        (Q1 = Q1[1:le[which.min(le)]])
        (Q2 = Q2[1:le[which.min(le)]])
        (Q3 = Q3[1:le[which.min(le)]])
        (Q4 = Q4[1:le[which.min(le)]])
        (x = matrix(0, nrow = max(lasts[which(min(raenge) == 
            ranks)]), ncol = 4))
        x[Q1, 1] = y[Q1]
        x[Q2, 2] = y[Q2]
        x[Q3, 3] = y[Q3]
        x[Q4, 4] = y[Q4]
        (cycle.trim = cycle(y)[-c((length(y) - (length(y) - nrow(x)) + 
            1):length(y))])
        (time.trim = trunc(time(y)[-c((length(y) - (length(y) - 
            nrow(x)) + 1):length(y))]))
        (phase.new = cycle.trim[length(cycle.trim)])
        (year.new = time.trim[length(time.trim)])
        (y = ts(y[-c((length(y) - (length(y) - nrow(x)) + 1):length(y))], 
            fre = 4, start = start(y), end = c(year.new, phase.new)))
    }
    else {
        (x = matrix(0, nrow = length(y), ncol = 4))
        x[Q1, 1] = y[Q1]
        x[Q2, 2] = y[Q2]
        x[Q3, 3] = y[Q3]
        x[Q4, 4] = y[Q4]
        year.new = end(y)[1]
        phase.new = end(y)[2]
    }
    (LHS = ts(x, fre = S, start = start(y), end = c(year.new, 
        phase.new)))
    (year.new.start = 1 + trunc(((start(y)[2] + lags) - 1)/S))
    (phase.new.start = cycle(y)[lags + 1])
    (rhs = ts(embed(cbind(y, x), dim = lags + 1), fre = S, start = c((year = 1 + 
        trunc(lags/S)), (st = 1 + trunc((lags)%%S)))))
    (lhs = rhs[, 1])
    (low = seq(1, (lags + 1) * (S + 1), by = S + 1))
    (rhs = rhs[, -c(low[1]:(S + 1), low[-1])])
    (x2 = ts(matrix(0, nrow = nrow(rhs), ncol = lags * S), fre = S, 
        start = c(year.new.start, phase.new.start)))
    if (lags > 0) {
        ii = 1
        (cols = rep(seq(1, S), times = 3 * (lags + 1)))
        while (ii <= lags) {
            (lags.blocks = seq(S * (ii - 1) + 1, S * ii, by = 1))
            x2[, cols[lags.blocks + ii] + S * (ii - 1)] = rhs[, 
                lags.blocks]
            ii = ii + 1
        }
    }
    (parX <- x2)
    assign("N", N, envir = environment(dter.comp))
    assign("lags", lags, envir = environment(dter.comp))
    assign("lags", lags, envir = environment(strbr.comp))
    if (!is.null(Tb)) {
        (W = strbr.comp(nseason = nseason, season = season, parX = x2, 
            y = y, Tb = Tb, tau = tau, isbreak.nsea = isbreak.nsea, 
            isbreak.sea = isbreak.sea))
    }
    else {
        if (!is.null(nseason) || !is.null(season)) {
            (W = dter.comp(nseason = nseason, season = season, 
                parX = x2, y = y, coding = coding))
        }
        else {
            W = x2
        }
    }
    if (class.est) {
        require(MASS)
        (H1 = ginv(t(W) %*% W))
        (Bhat = H1 %*% t(W) %*% lhs)
        rownames(Bhat) = colnames(W)
        (Pw = W %*% H1 %*% t(W))
        (Mw = diag(nrow(W)) - Pw)
        (RSS = as.numeric(t(lhs) %*% Mw %*% lhs))
        (npar = nrow(Bhat))
        (sig2 = RSS/(N - npar - S - 1))
        (aic = round(log(sig2) + 2 * (npar/(N - S)), digits = 10))
        (bic = round(log(sig2) + (npar * log(N - S))/(N - S), 
            digits = 10))
        (fpe = round(((N + S + npar)/(N + S - npar)) * sig2, 
            digits = 10))
        (hq = round(log(sig2) + (2 * (npar) * log(log(N - S)))/(N - 
            S), digits = 10))
    }
    else {
        sig2 = aic = bic = fpe = hq = 0
        Bhat = list(Bhat = 0)
    }
    dcomps = list(sdeter = sdeter, isbreak.nsea = isbreak.nsea, 
        isbreak.sea = isbreak.sea)
    return(invisible(new("bayes.par", dcomps = dcomps, p = lags, 
        Tb = Tb, aic = aic, bic = bic, hq = hq, fpe = fpe, lhs = lhs, 
        LHS = LHS[-c(1:lags), ], W = W, est.classic = Bhat, sigma2 = sig2)))
}
