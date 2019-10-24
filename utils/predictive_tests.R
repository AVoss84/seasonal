function (g1, g2, g3 = NULL, bprior = c(alpha = 1, beta = 1), 
    level = 0.05, p0 = 0.5, ...) 
{
    signtest = dget("signtest.R")
    alpha = bprior[1]
    beta = bprior[2]
    x = seq(0, 1, by = 0.001)
    ev = NULL
    if (is.null(g3)) {
        g1 = as.vector(g1)
        g2 = as.vector(g2)
        pv = matrix(, nrow = 1, ncol = 4, dimnames = list(c("g1-g2:"), 
            c("Sign", "Wilcoxon", "HPD(low)", "HPD(up)")))
        (loss_diffs = g1 - g2)
        (ind = loss_diffs > 0)
        (S2 = sum(ind))
        TT = length(ind)
        stest = signtest(g1, g2, alternative = "two.sided")
        pv[1, 1] = stest$rval$p.val
        pv[1, 2] = wilcox.test(g1, g2, alternative = "two.sided", 
            ...)$p.value
        (bayes = structure(list(alpha = alpha, beta = beta, prior.mean = alpha/(alpha + 
            beta), prior.var = alpha * beta/(((alpha + beta)^2) * 
            (alpha + beta + 1)), fp = dbeta(x, shape1 = alpha + 
            S2, shape2 = beta + TT - S2), post.mean = (alpha + 
            S2)/(alpha + beta + TT), post.var = (alpha + S2) * 
            (beta + TT - S2)/(((alpha + beta + TT)^2) * (alpha + 
            beta + TT + 1)), hpd = qbeta(c(level/2, 1 - level/2), 
            shape1 = alpha + S2, shape2 = beta + TT - S2))))
        pv[1, 3:4] = bayes$hpd
        (post.mean = as.numeric(bayes$post.mean))
        (post.var = as.numeric(bayes$post.var))
        (num = dbinom(x = S2, size = TT, prob = 0.5, log = F))
        mx = choose(TT, S2) * as.double(beta(alpha + S2, beta + 
            TT - S2))/as.double(beta(alpha, beta))
        (probH0 = (num * p0/(num * p0 + (1 - p0) * mx)))
        (x.FBST = seq(0, 1, by = 0.002))
        (y.FBST = dbeta(x.FBST, alpha + S2, beta + TT - S2))
        (g.star = y.FBST[sum(x.FBST <= 0.5)])
        (T.set = x.FBST[y.FBST > g.star])
        (ev.over = pbeta(max(T.set), alpha + S2, beta + TT - 
            S2) - pbeta(min(T.set), alpha + S2, beta + TT - S2))
        (ev = 1 - ev.over)
    }
    else {
        g1 = as.vector(g1)
        g2 = as.vector(g2)
        g3 = as.vector(g3)
        post.mean = post.var = probH0 = NULL
        pv = matrix(, nrow = 3, ncol = 4, dimnames = list(c("g1-g2:", 
            "g1-g3:", "g2-g3:"), c("Sign", "Wilcoxon", "HPD(low)", 
            "HPD(up)")))
        (loss_diffs = g1 - g2)
        (ind = loss_diffs > 0)
        (S2 = sum(ind))
        TT = length(ind)
        stest = signtest(g1, g2, alternative = "two.sided")
        pv[1, 1] = stest$rval$p.val
        pv[1, 2] = wilcox.test(g1, g2, alternative = "two.sided", 
            ...)$p.value
        (bayes = structure(list(alpha = alpha, beta = beta, prior.mean = alpha/(alpha + 
            beta), prior.var = alpha * beta/(((alpha + beta)^2) * 
            (alpha + beta + 1)), fp = dbeta(x, shape1 = alpha + 
            S2, shape2 = beta + TT - S2), post.mean = (alpha + 
            S2)/(alpha + beta + TT), post.var = (alpha + S2) * 
            (beta + TT - S2)/(((alpha + beta + TT)^2) * (alpha + 
            beta + TT + 1)), hpd = qbeta(c(level/2, 1 - level/2), 
            shape1 = alpha + S2, shape2 = beta + TT - S2))))
        pv[1, 3:4] = bayes$hpd
        (post.mean = c(post.mean, as.numeric(bayes$post.mean)))
        (post.var = c(post.var, as.numeric(bayes$post.var)))
        (num = dbinom(x = S2, size = TT, prob = 0.5, log = F))
        mx = choose(TT, S2) * as.double(beta(alpha + S2, beta + 
            TT - S2))/as.double(beta(alpha, beta))
        (probH0 = c(probH0, (num * p0/(num * p0 + (1 - p0) * 
            mx))))
        (x.FBST = seq(0, 1, by = 0.002))
        (y.FBST = dbeta(x.FBST, alpha + S2, beta + TT - S2))
        (g.star = y.FBST[sum(x.FBST <= 0.5)])
        (T.set = x.FBST[y.FBST > g.star])
        (ev.over = pbeta(max(T.set), alpha + S2, beta + TT - 
            S2) - pbeta(min(T.set), alpha + S2, beta + TT - S2))
        (ev = c(ev, 1 - ev.over))
        (loss_diffs = g1 - g3)
        (ind = loss_diffs > 0)
        (S2 = sum(ind))
        stest = signtest(g1, g3, alternative = "two.sided")
        pv[2, 1] = stest$rval$p.val
        pv[2, 2] = wilcox.test(g1, g3, alternative = "two.sided", 
            ...)$p.value
        (bayes = structure(list(alpha = alpha, beta = beta, prior.mean = alpha/(alpha + 
            beta), prior.var = alpha * beta/(((alpha + beta)^2) * 
            (alpha + beta + 1)), fp = dbeta(x, shape1 = alpha + 
            S2, shape2 = beta + TT - S2), post.mean = (alpha + 
            S2)/(alpha + beta + TT), post.var = (alpha + S2) * 
            (beta + TT - S2)/(((alpha + beta + TT)^2) * (alpha + 
            beta + TT + 1)), hpd = qbeta(c(level/2, 1 - level/2), 
            shape1 = alpha + S2, shape2 = beta + TT - S2))))
        pv[2, 3:4] = bayes$hpd
        (post.mean = c(post.mean, as.numeric(bayes$post.mean)))
        (post.var = c(post.var, as.numeric(bayes$post.var)))
        (num = dbinom(x = S2, size = TT, prob = 0.5, log = F))
        mx = choose(TT, S2) * as.double(beta(alpha + S2, beta + 
            TT - S2))/as.double(beta(alpha, beta))
        (probH0 = c(probH0, (num * p0/(num * p0 + (1 - p0) * 
            mx))))
        (x.FBST = seq(0, 1, by = 0.002))
        (y.FBST = dbeta(x.FBST, alpha + S2, beta + TT - S2))
        (g.star = y.FBST[sum(x.FBST <= 0.5)])
        (T.set = x.FBST[y.FBST > g.star])
        (ev.over = pbeta(max(T.set), alpha + S2, beta + TT - 
            S2) - pbeta(min(T.set), alpha + S2, beta + TT - S2))
        (ev = c(ev, 1 - ev.over))
        (loss_diffs = g2 - g3)
        (ind = loss_diffs > 0)
        (S2 = sum(ind))
        stest = signtest(g2, g3, alternative = "two.sided")
        pv[3, 1] = stest$rval$p.val
        pv[3, 2] = wilcox.test(g2, g3, alternative = "two.sided", 
            ...)$p.value
        (bayes = structure(list(alpha = alpha, beta = beta, prior.mean = alpha/(alpha + 
            beta), prior.var = alpha * beta/(((alpha + beta)^2) * 
            (alpha + beta + 1)), fp = dbeta(x, shape1 = alpha + 
            S2, shape2 = beta + TT - S2), post.mean = (alpha + 
            S2)/(alpha + beta + TT), post.var = (alpha + S2) * 
            (beta + TT - S2)/(((alpha + beta + TT)^2) * (alpha + 
            beta + TT + 1)), hpd = qbeta(c(level/2, 1 - level/2), 
            shape1 = alpha + S2, shape2 = beta + TT - S2))))
        pv[3, 3:4] = bayes$hpd
        (post.mean = c(post.mean, as.numeric(bayes$post.mean)))
        (post.var = c(post.var, as.numeric(bayes$post.var)))
        (num = dbinom(x = S2, size = TT, prob = 0.5, log = F))
        mx = choose(TT, S2) * as.double(beta(alpha + S2, beta + 
            TT - S2))/as.double(beta(alpha, beta))
        (probH0 = c(probH0, (num * p0/(num * p0 + (1 - p0) * 
            mx))))
        (x.FBST = seq(0, 1, by = 0.002))
        (y.FBST = dbeta(x.FBST, alpha + S2, beta + TT - S2))
        (g.star = y.FBST[sum(x.FBST <= 0.5)])
        (T.set = x.FBST[y.FBST > g.star])
        (ev.over = pbeta(max(T.set), alpha + S2, beta + TT - 
            S2) - pbeta(min(T.set), alpha + S2, beta + TT - S2))
        (ev = c(ev, 1 - ev.over))
    }
    return(list(pvalues = pv, post.means = post.mean, post.var = post.var, 
        probH0 = probH0, evidence_value_H0 = ev))
}
