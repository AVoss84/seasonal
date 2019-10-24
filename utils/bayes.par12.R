function (data, p = NULL, deter = c("const", "trend", "both"), 
    sdeter = c("const", "both"), coding = c("mean", "effect"), 
    Tb = NULL, tau = 3, isbreak.nsea = c(0, 1), isbreak.sea = c(1, 
        0), class.est = T) 
{
    dter.compS = dget("dter.compS.R")
    strbr.compS = dget("strbr.compS.R")
    setClass("bayes.par12", representation(lhs = "ANY", LHS = "ANY", 
        p = "numeric", aic = "numeric", bic = "numeric", hq = "numeric", 
        fpe = "numeric", est.classic = "ANY", sigma2 = "numeric", 
        W = "ANY", Tb = "ANY"))
    setMethod("show", "bayes.par12", function(object) {
        if (!is.null(object@Tb)) {
            cat("\n", rep("-", 17), "\n", sprintf(" Monthly Gaussian PAR(%i) model:\n", 
                object@p), rep("-", 17), "\n", sprintf(" Break date T_{b} = %i\n", 
                object@Tb), rep("-", 17), "\n\n")
        }
        else {
            cat("\n", rep("-", 17), "\n", rep("-", 17), "\n", 
                sprintf(" Monthly Gaussian PAR(%i) model:\n", 
                  object@p), rep("-", 17), "\n", rep("-", 17), 
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
    nseason = deter[1]
    season = sdeter[1]
    S = 12
    N = length(y) + S
    lags = p
    coding = coding[1]
    (Q1 = which(cycle(y) == 1))
    (Q2 = which(cycle(y) == 2))
    (Q3 = which(cycle(y) == 3))
    (Q4 = which(cycle(y) == 4))
    (Q5 = which(cycle(y) == 5))
    (Q6 = which(cycle(y) == 6))
    (Q7 = which(cycle(y) == 7))
    (Q8 = which(cycle(y) == 8))
    (Q9 = which(cycle(y) == 9))
    (Q10 = which(cycle(y) == 10))
    (Q11 = which(cycle(y) == 11))
    (Q12 = which(cycle(y) == 12))
    (le = c(length(Q1), length(Q2), length(Q3), length(Q4), length(Q5), 
        length(Q6), length(Q7), length(Q8), length(Q9), length(Q10), 
        length(Q11), length(Q12)))
    (year.new = end(y)[1])
    (phase.new = end(y)[2])
    (year.new.start = 1 + trunc(((start(y)[2] + lags) - 1)/S))
    (phase.new.start = cycle(y)[lags + 1])
    x = matrix(0, nrow = length(y), ncol = S)
    x[Q1, 1] = y[Q1]
    x[Q2, 2] = y[Q2]
    x[Q3, 3] = y[Q3]
    x[Q4, 4] = y[Q4]
    x[Q5, 5] = y[Q5]
    x[Q6, 6] = y[Q6]
    x[Q7, 7] = y[Q7]
    x[Q8, 8] = y[Q8]
    x[Q9, 9] = y[Q9]
    x[Q10, 10] = y[Q10]
    x[Q11, 11] = y[Q11]
    x[Q12, 12] = y[Q12]
    (LHS = ts(x, fre = S, start = start(y), end = c(year.new, 
        phase.new)))
    (rhs = ts(embed(cbind(y, x), dim = lags + 1), fre = S, start = c(year.new.start, 
        phase.new.start)))
    (lhs = rhs[, 1])
    (low = seq(1, (lags + 1) * (S + 1), by = S + 1))
    (rhs = rhs[, -c(low[1]:(S + 1), low[-1])])
    stopifnot(dim(rhs)[2]%%S == 0)
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
    assign("N", N, envir = environment(dter.compS))
    assign("lags", lags, envir = environment(dter.compS))
    if (!is.null(Tb)) {
        (W = strbr.compS(nseason = nseason, season = season, 
            parX = parX, lags = lags, y = y, Tb = Tb, isbreak.nsea = isbreak.nsea, 
            isbreak.sea = isbreak.sea, tau = tau))
    }
    else {
        if (!is.null(nseason) || !is.null(season)) {
            (W = dter.compS(nseason = nseason, season = season, 
                parX = x2, y = y, S = S, coding = coding))
        }
        else {
            W = x2
            (nm = c(paste(rep(c(paste("ys", 1:S, sep = "")), 
                times = lags), ".", rep(1:lags, each = S), sep = "")))
            colnames(W) = nm
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
        Bhat = list(coef = 0)
    }
    return(invisible(new("bayes.par12", p = lags, Tb = Tb, aic = aic, 
        bic = bic, fpe = fpe, hq = hq, lhs = lhs, LHS = LHS[-c(1:lags), 
            ], W = W, est.classic = Bhat, sigma2 = sig2)))
}
