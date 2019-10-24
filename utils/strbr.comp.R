function (nseason = NULL, season = NULL, parX, y, Tb, isbreak.nsea = c(1, 
    1), isbreak.sea = c(1, 1), tau = 3) 
{
    design.break = dget("design.break.R")
    design = dget("design.R")
    assign("lags", lags, envir = environment(design))
    S = 4
    N = length(y) + S
    coding = c("mean")
    if (!is.null(nseason)) {
        stopifnot(any(nseason == c("const", "trend", "both")))
    }
    if (!is.null(season)) {
        stopifnot(any(season == c("const", "trend", "both")))
    }
    if (Tb - (lags + tau) <= 0) {
        warning("Break date must be at least T_{b} = lags + tau + 1 due to potential singularity problems.\n")
    }
    if (Tb > N - tau) {
        stop("Break date is greater than N-tau (=upper bound)!\n")
    }
    if (!is.null(nseason) && is.null(season)) {
        if ((nseason[1] == "const" && sum(isbreak.nsea) == 2)) {
            stop("Option not allowed!\n")
        }
        if (nseason[1] == "const" && all(isbreak.nsea == c(1, 
            0))) {
            const.break = array(0, dim = c(N - (lags + S), 2))
            const.break[1:(Tb - lags), 1] = rep(1, Tb - lags)
            const.break[(Tb - lags + 1):(N - (lags + S)), 2] = rep(1, 
                N - S - Tb)
            (W = ts.union(rhs = parX, const.br = const.break))
            (nm = c(paste(rep(c(paste("ys", 1:S, sep = "")), 
                times = lags), ".", rep(1:lags, each = S), sep = ""), 
                "dreg1", "dreg2"))
            colnames(W) = nm
        }
        if (nseason[1] == "both" && all(isbreak.nsea == c(0, 
            1))) {
            both.break = array(0, dim = c(N - (lags + S), 3))
            both.break[, 1] = rep(1, N - S - lags)
            both.break[1:(Tb - lags), 2] = (lags + 1):Tb
            both.break[(Tb - lags + 1):(N - S - lags), 3] = (Tb + 
                1):(N - S)
            (W = ts.union(rhs = parX, both.break = both.break))
            (nm = c(paste(rep(c(paste("ys", 1:S, sep = "")), 
                times = lags), ".", rep(1:lags, each = S), sep = ""), 
                "dreg1", "treg1", "treg2"))
            colnames(W) = nm
        }
        if (nseason[1] == "both" && all(isbreak.nsea == c(1, 
            0))) {
            both.break = array(0, dim = c(N - (lags + S), 3))
            both.break[1:(Tb - lags), 1] = rep(1, Tb - lags)
            both.break[(Tb - lags + 1):(N - (lags + S)), 2] = rep(1, 
                N - S - Tb)
            both.break[, 3] = (lags + 1):(N - S)
            (W = ts.union(rhs = parX, both.break = both.break))
            (nm = c(paste(rep(c(paste("ys", 1:S, sep = "")), 
                times = lags), ".", rep(1:lags, each = S), sep = ""), 
                "dreg1", "dreg2", "treg1"))
            colnames(W) = nm
        }
        if (nseason[1] == "both" && all(isbreak.nsea == c(1, 
            1))) {
            both.break = array(0, dim = c(N - (lags + S), 4))
            both.break[1:(Tb - lags), 1] = rep(1, Tb - lags)
            both.break[(Tb - lags + 1):(N - (lags + S)), 2] = rep(1, 
                N - S - Tb)
            both.break[1:(Tb - lags), 3] = (lags + 1):Tb
            both.break[(Tb - lags + 1):(N - S - lags), 4] = (Tb + 
                1):(N - S)
            (W = ts.union(rhs = parX, both.break = both.break))
            (nm = c(paste(rep(c(paste("ys", 1:S, sep = "")), 
                times = lags), ".", rep(1:lags, each = S), sep = ""), 
                "dreg1", "dreg2", "treg1", "treg2"))
            colnames(W) = nm
        }
    }
    if (is.null(nseason) && !is.null(season)) {
        if ((season[1] == "const" && sum(isbreak.sea) == 2)) {
            stop("Option not allowed!\n")
        }
        if (season[1] == "const" && all(isbreak.sea == c(1, 0))) {
            (W = ts.union(rhs = parX, sdum.br = design.break(mat = parX, 
                y = y, Tb = Tb, lags = lags)))
            (nm = c(paste(rep(c(paste("ys", 1:S, sep = "")), 
                times = lags), ".", rep(1:lags, each = S), sep = "")))
            colnames(W) = c(nm, rep(paste(rep(paste("dreg", 1:2, 
                sep = ""), each = 4), 1:S, sep = ""), times = 1))
        }
        if (season[1] == "both" && all(isbreak.sea == c(1, 1))) {
            sdum.br = design.break(mat = parX, y = y, Tb = Tb, 
                lags = lags)
            strend.br = design.break(mat = parX, y = y, Tb = Tb, 
                lags = lags) * ((lags + 1):(N - S))
            sboth.br = cbind(sdum.br, strend.br)
            (W = ts.union(rhs = parX, sboth.br))
            (nm = c(paste(rep(c(paste("ys", 1:S, sep = "")), 
                times = lags), ".", rep(1:lags, each = S), sep = "")))
            colnames(W) = c(nm, rep(paste(rep(paste("dreg", 1:2, 
                sep = ""), each = 4), 1:S, sep = ""), times = 1), 
                rep(paste(rep(paste("treg", 1:2, sep = ""), each = 4), 
                  1:S, sep = ""), times = 1))
        }
        if (season[1] == "both" && all(isbreak.sea == c(1, 0))) {
            sdum.br = design.break(mat = parX, y = y, Tb = Tb, 
                lags = lags)
            strend = design(mat = parX, y = y, coding = coding) * 
                ((lags + 1):(N - S))
            sboth.br = cbind(sdum.br, strend)
            (W = ts.union(rhs = parX, sboth.br))
            (nm = c(paste(rep(c(paste("ys", 1:S, sep = "")), 
                times = lags), ".", rep(1:lags, each = S), sep = "")))
            colnames(W) = c(nm, rep(paste(rep(paste("dreg", 1:2, 
                sep = ""), each = 4), 1:S, sep = ""), times = 1), 
                rep(paste(rep(paste("treg", 1, sep = ""), each = 4), 
                  1:S, sep = ""), times = 1))
        }
        if (season[1] == "both" && all(isbreak.sea == c(0, 1))) {
            sdum = design(mat = parX, y = y, coding = coding)
            strend.br = design.break(mat = parX, y = y, Tb = Tb, 
                lags = lags) * ((lags + 1):(N - S))
            sboth.br = cbind(sdum, strend.br)
            (W = ts.union(rhs = parX, sboth.br))
            (nm = c(paste(rep(c(paste("ys", 1:S, sep = "")), 
                times = lags), ".", rep(1:lags, each = S), sep = "")))
            colnames(W) = c(nm, rep(paste(rep(paste("dreg", 1:1, 
                sep = ""), each = 4), 1:S, sep = ""), times = 1), 
                rep(paste(rep(paste("treg", 1:2, sep = ""), each = 4), 
                  1:S, sep = ""), times = 1))
        }
    }
    if (!is.null(nseason) && !is.null(season)) {
        if (nseason[1] == "const" && season[1] == "const") {
            stop("Option not supported!\n")
        }
        if (nseason[1] == "trend" && season[1] == "const" && 
            all(isbreak.nsea == c(0, 0)) && all(isbreak.sea == 
            c(1, 0))) {
            trend = (lags + 1):(N - S)
            sdum.br = design.break(mat = parX, y = y, Tb = Tb, 
                lags = lags)
            (W = ts.union(rhs = parX, . = cbind(sdum.br, trend)))
            (nm = c(paste(rep(c(paste("ys", 1:S, sep = "")), 
                times = lags), ".", rep(1:lags, each = S), sep = "")))
            colnames(W) = c(nm, rep(paste(rep(paste("dreg", 1:2, 
                sep = ""), each = 4), 1:S, sep = ""), times = 1), 
                rep(paste(rep(paste("treg", 1:1, sep = ""), each = 1), 
                  1:1, sep = ""), times = 1))
        }
        if (nseason[1] == "trend" && season[1] == "const" && 
            all(isbreak.nsea == c(0, 1)) && all(isbreak.sea == 
            c(1, 0))) {
            trend.break = array(0, dim = c(N - (lags + S), 2))
            trend.break[1:(Tb - lags), 1] = (lags + 1):Tb
            trend.break[(Tb - lags + 1):(N - S - lags), 2] = (Tb + 
                1):(N - S)
            sdum.br = design.break(mat = parX, y = y, Tb = Tb, 
                lags = lags)
            (W = ts.union(rhs = parX, . = cbind(sdum.br, trend.break)))
            (nm = c(paste(rep(c(paste("ys", 1:S, sep = "")), 
                times = lags), ".", rep(1:lags, each = S), sep = "")))
            colnames(W) = c(nm, rep(paste(rep(paste("dreg", 1:2, 
                sep = ""), each = 4), 1:S, sep = ""), times = 1), 
                rep(paste(rep(paste("treg", 1:2, sep = ""), each = 1), 
                  1:1, sep = ""), times = 1))
        }
        if ((nseason[1] == "trend" && season[1] == "both") || 
            (nseason[1] == "const" && season[1] == "both")) {
            stop("Option not supported!\n")
        }
    }
    return(W)
}
