function (y, p, n.ahead, M, burnin, restr2stat = F, sbreak = T, 
    sbreak.fix = F, s2.hyp = c(nu, lambda), in.samp.fc = T, le.W = 100, 
    low_up = c(4, 4), deter = NULL, sdeter = "const", isbreak.nsea = c(1, 
        0), isbreak.sea = c(1, 0), tau = 3) 
{
    library(MASS)
    bayes.par = dget("bayes.par.R")
    strbr.comp = dget("strbr.comp.R")
    mats01 = dget("mats01.R")
    nseason = deter
    season = sdeter
    coding = c("mean", "effect")
    S = 4
    if (frequency(y) != S) {
        stop("Time series has wrong seasonal frequency!\n")
    }
    if (in.samp.fc) {
        (y.true.T = ts(y[c((length(y) - n.ahead + 1):length(y))], 
            end = end(y), fre = S))
        (y = ts(y[-c((length(y) - n.ahead + 1):length(y))], start = start(y), 
            fre = S))
    }
    else {
        y.true.T = NULL
    }
    (n = length(y))
    (n.new = n - p)
    (z = rep(mean(y), times = n.ahead))
    (y.aug = ts(c(y, z), fre = S, start = start(y)))
    (k = c(1, floor((n - 1) * 0.5), n.new))
    (date = if (sbreak) {
        k[2]
    }
    else {
        NULL
    })
    if (sbreak.fix && !sbreak) {
        date = floor(runif(1, min = k[1] + tau + p, max = k[3] - 
            1))
        message("Draw break date from uniform instrumental density instead of posterior.\n")
    }
    if (sbreak.fix && sbreak) {
        message("Option is not supported. Set \"sbreak <- FALSE\" if uniform density should be used.\nMultinomial posterior density is used instead.\n")
    }
    MCsim = M + burnin
    (out1 = bayes.par(data = y, p = p, Tb = date, deter = deter, 
        sdeter = sdeter, isbreak.nsea = isbreak.nsea, isbreak.sea = isbreak.sea, 
        tau = tau))
    (lhs = out1@lhs)
    (X.til = out1@W)
    (X = X.til[, 1:(p * S)])
    out2 = bayes.par(data = y.aug, class.est = T, p = p, Tb = date, 
        deter = deter, sdeter = sdeter, isbreak.nsea = isbreak.nsea, 
        isbreak.sea = isbreak.sea, tau = tau)
    (W.star = out2@W)
    (W.star = ts(matrix(W.star[(nrow(W.star) - n.ahead + 1):nrow(W.star), 
        ], nrow = n.ahead), fre = S, end = end(W.star)))
    s2.eps = Tb = prod.unit = numeric(MCsim)
    if (sbreak) {
        Tb[1] = date
    }
    else {
        Tb = NULL
    }
    b = matrix(, ncol = ncol(X.til), nrow = MCsim)
    W.post = matrix(, ncol = n.ahead, nrow = MCsim)
    (W.post[1, ] = z)
    (s2.eps[1] = out2@sigma2)
    (b[1, ] = out2@est.classic)
    (OLS = out1@est.classic)
    (prod.unit[1] = abs(prod(b[1, 1:(p * S)])))
    margloglike = NULL
    (W_supp = seq(min(y) - low_up[1] * sd(y), max(y) + low_up[2] * 
        sd(y), length.out = le.W))
    print(W_supp)
    marg.Ws = matrix(0, nrow = n.ahead, ncol = length(W_supp), 
        dimnames = list(paste("W", 1:n.ahead, sep = ""), paste(1:length(W_supp), 
            sep = "")))
    if ((MCsim)%%le.W) {
        stop(" 'MCsim' must be a multiple of length(W_supp).\n")
    }
    else {
        cat("There are", MCsim/le.W, "draws per point in 'W_supp' used.\n")
    }
    c1 = 100
    c2 = 100
    (b0 = c(rep(1, times = p * S), rep(0, ncol(X.til) - p * S)))
    V = diag(c(rep(c1, times = p * S), rep(c2, ncol(X.til) - 
        p * S)))
    a = s2.hyp[1]
    lambda = s2.hyp[2]
    loglike = function(X.til, lhs, s2.eps) {
        (beta = ginv(t(X.til) %*% X.til) %*% t(X.til) %*% lhs)
        (err = lhs - X.til %*% beta)
        nNew = length(lhs)
        -(nNew/2) * log(2 * pi) - (nNew/2) * log(s2.eps) - (1/(2 * 
            s2.eps)) * t(err) %*% err
    }
    ii = 2
    hh = gg = 1
    for (ii in 2:MCsim) {
        cat("Draw nr.", ii, "\n")
        if (sbreak) {
            assign("lags", p, envir = environment(strbr.comp))
            (kBeg = k[1] + tau + p)
            (kEnd = k[3] - 1)
            (prob = numeric(kEnd - kBeg + 1))
            (j = kBeg)
            while (j <= kEnd) {
                (X.til = strbr.comp(nseason = nseason, season = season, 
                  parX = X, y = y, Tb = j, tau = tau, isbreak.nsea = isbreak.nsea, 
                  isbreak.sea = isbreak.sea))
                (beta = ginv(t(X.til) %*% X.til) %*% t(X.til) %*% 
                  lhs)
                (err = lhs[kBeg:kEnd] - X.til[kBeg:kEnd, ] %*% 
                  beta)
                (prob[j - p - tau] = -length(prob) * log(sqrt(s2.eps[ii - 
                  1])) - sum(err^2)/(2 * s2.eps[ii - 1]))
                j = j + 1
                print(j)
            }
            prob = exp(prob - max(prob))
            prob = prob/sum(prob)
            if (length(prob) > 0) {
                mn = rmultinom(n = 1, size = length(prob), prob = prob)
                postki = which.max(mn)
                Tb[ii] = postki + k[1] + p + tau - 1
            }
            cat("\nDrawn break date:", Tb[ii], "\n\n")
        }
        (y.aug = ts(c(y, W.post[ii - 1, ]), fre = S, start = start(y)))
        (date = if (sbreak) {
            Tb[ii]
        }
        else {
            NULL
        })
        if (sbreak.fix && !sbreak) {
            date = floor(runif(1, min = k[1] + tau + p, max = k[3] - 
                1))
            print(date)
        }
        out2 = bayes.par(data = y.aug, class.est = FALSE, p = p, 
            Tb = date, deter = deter, sdeter = sdeter, isbreak.nsea = isbreak.nsea, 
            isbreak.sea = isbreak.sea, tau = tau)
        (W.star = out2@W)
        (fc.rows = (nrow(W.star) - n.ahead + 1):nrow(W.star))
        (X.til = W.star[-fc.rows, ])
        (X = X.til[, 1:(p * S)])
        (W.star = matrix(W.star[fc.rows, ], nrow = n.ahead, ncol = ncol(W.star), 
            dimnames = list(paste("k=", 1:n.ahead, sep = ""), 
                NULL)))
        (R = t(X.til) %*% X.til + t(as.matrix(W.star)) %*% as.matrix(W.star) + 
            solve(V))
        (b.mean = ginv(R) %*% (t(X.til) %*% lhs + solve(V) %*% 
            b0 + t(as.matrix(W.star)) %*% as.matrix(W.post[ii - 
            1, ])))
        (b.sigma = s2.eps[ii - 1] * ginv(R))
        (U = chol(b.sigma, pivot = F))
        if (restr2stat) {
            repeat {
                (b.new = b.mean + t(U) %*% as.matrix(rnorm(length(b.mean), 
                  0, 1)))
                b[ii, ] = b.new
                mat = mats01(b[ii, 1:(p * S)], p = p, S = S)
                prod.unit[ii] = abs(mat$det01 - 1)
                if (ifelse(restr2stat, yes = !(abs(mat$det01) < 
                  0.05), no = TRUE)) {
                  break
                }
                else {
                  cat("\nPosterior draws of beta discarded.\n")
                }
            }
        }
        else {
            (b[ii, ] = b.mean + t(U) %*% as.matrix(rnorm(length(b.mean), 
                0, 1)))
            prod.unit[ii] = abs(prod(b[ii, 1:(p * S)]))
        }
        (a.star = n - p + n.ahead + length(b[ii, ] + a))
        (b.star = t(as.matrix(lhs) - X.til %*% b[ii, ]) %*% (as.matrix(lhs) - 
            X.til %*% b[ii, ]) + t(as.matrix(W.post[ii - 1, ]) - 
            W.star %*% b[ii, ]) %*% (as.matrix(W.post[ii - 1, 
            ]) - W.star %*% b[ii, ]) + t(b[ii, ] - b0) %*% solve(V) %*% 
            (b[ii, ] - b0) + lambda)
        (s2.eps[ii] = 1/rgamma(1, shape = a.star, rate = b.star))
        step = kk = lags = 1
        (W.mean = t(W.star[step, ]) %*% b[ii, ])
        (W.post[ii, step] = rnorm(1, mean = W.mean, sd = sqrt(s2.eps[ii])))
        (cycle.orig = as.numeric(cycle(y)))
        (cycle.aug = as.numeric(cycle(y.aug)))
        if (n.ahead > 1) {
            repeat {
                if (kk <= (n.ahead - 1) && lags <= p) {
                  (W.star[(kk + 1), cycle.aug[length(cycle.orig) + 
                    (kk + 1)] + (lags - 1) * S] = W.post[ii, 
                    step])
                  kk = kk + 1
                  lags = lags + 1
                }
                else {
                  break
                }
            }
            marg.Ws[step, hh] = (marg.Ws[step, hh] * (gg - 1) + 
                dnorm(W_supp[hh], mean = t(W.star[step, ]) %*% 
                  b[ii, ], sd = sqrt(s2.eps[ii])))/gg
        }
        step = 2
        while (step <= n.ahead) {
            kk = step
            (W.mean = t(W.star[step, ]) %*% b[ii, ])
            (W.post[ii, step] = rnorm(1, mean = W.mean, sd = sqrt(s2.eps[ii])))
            lags = 1
            if (n.ahead > 1) {
                repeat {
                  if (kk <= (n.ahead - 1) && lags <= p) {
                    (W.star[(kk + 1), cycle.aug[length(cycle.orig) + 
                      (kk + 1)] + (lags - 1) * S] = W.post[ii, 
                      step])
                    kk = kk + 1
                    lags = lags + 1
                  }
                  else {
                    break
                  }
                }
            }
            marg.Ws[step, hh] = (marg.Ws[step, hh] * (gg - 1) + 
                dnorm(W_supp[hh], mean = t(W.star[step, ]) %*% 
                  b[ii, ], sd = sqrt(s2.eps[ii])))/gg
            step = step + 1
        }
        if (gg < (MCsim/le.W)) {
            gg = gg + 1
        }
        else {
            gg = 1
            if (hh < le.W) {
                hh = hh + 1
            }
        }
        (marglog = loglike(X.til = X.til, lhs = lhs, s2.eps = s2.eps[ii]))
        margloglike = c(margloglike, marglog)
    }
    (marglike = exp(margloglike - max(margloglike)))
    (MLH = 1/mean(1/marglike))
    require(Bolstad)
    erw = vari = hpd95 = NULL
    for (ii in 1:nrow(marg.Ws)) {
        (con = sintegral(x = W_supp, fx = marg.Ws[ii, ])$value)
        (erw = c(erw, sintegral(x = W_supp, fx = W_supp * marg.Ws[ii, 
            ]/con)$value))
        (vari = c(vari, sintegral(x = W_supp, fx = ((W_supp - 
            erw[ii])^2) * marg.Ws[ii, ]/con)$value))
        marg.Ws[ii, ] = marg.Ws[ii, ]/con
        (cdf = sintegral(W_supp, marg.Ws[ii, ])$cdf)
        (q0025 = cdf$x[min(which(cdf$y >= 0.025))])
        (q0975 = cdf$x[min(which(cdf$y >= 0.975))])
        (hpd95 = rbind(hpd95, c(q0025, q0975)))
    }
    if (in.samp.fc) {
        (BIAS = erw - y.true.T)
    }
    else {
        BIAS = NULL
    }
    return(list(control = c(M, burnin, n.ahead), OLS = OLS, s2.eps = s2.eps, 
        Tb = Tb, prod.unit = prod.unit, b = b, W.post = W.post, 
        ytrunc = y, y.true.T = y.true.T, n.ahead = n.ahead, MLH = MLH, 
        marg.Ws = marg.Ws, W_supp = W_supp, hpd95 = hpd95, pred.mean = list(Mean = erw, 
            Variance = vari, y.true.T = y.true.T, BIAS = BIAS)))
}
