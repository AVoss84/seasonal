function (y.trunc, dcomp, sdcomp, Tb.true = NULL, n.ahead, ...) 
{
    require(TSA)
    stopifnot(n.ahead%%frequency(y.trunc) == 0)
    if (!is.null(Tb.true) && Tb.true >= length(y.trunc)) {
        stop("Break point outside the truncated sample. Set Tb.true <- NULL\n")
    }
    (AO = (seq(y.trunc) == Tb.true))
    if (is.null(dcomp) && !is.null(sdcomp)) {
        switch(sdcomp, const = {
            (month = season(y.trunc))
            if (is.null(Tb.true)) {
                (determ = cbind(model.matrix(~month - 1)))
            } else {
                (determ = as.matrix(model.matrix(~month + AO - 
                  1)))
            }
        }, both = {
            (month = season(y.trunc))
            (trend = cbind(model.matrix(~month - 1)) * seq(length(y.trunc)))
            (determ = as.matrix(model.matrix(~month + trend - 
                1)))
            if (is.null(Tb.true)) {
                (determ = as.matrix(model.matrix(~month + trend - 
                  1)))
            } else {
                (determ = as.matrix(model.matrix(~month + trend + 
                  AO - 1)))
            }
        })
        (model1b = arimax(y.trunc, xreg = determ[, -1], include.mean = T, 
            ...))
        (AO.new = ((length(y.trunc) + 1):(length(y.trunc) + n.ahead) == 
            Tb.true))
        switch(sdcomp, const = {
            (month_new = season(ts(rep(1, n.ahead), start = c(end(y.trunc)[1], 
                end(y.trunc)[2] + 1), fre = frequency(y.trunc))))
            if (is.null(Tb.true)) {
                (determ.new = as.matrix(model.matrix(~month_new - 
                  1)))
            } else {
                (determ.new = as.matrix(model.matrix(~month_new + 
                  AO.new - 1)))
            }
        }, both = {
            (month_new = season(ts(rep(1, n.ahead), start = c(end(y.trunc)[1], 
                end(y.trunc)[2] + 1), fre = frequency(y.trunc))))
            (trend_new = cbind(model.matrix(~month_new - 1)) * 
                (length(y.trunc) + seq(n.ahead)))
            if (is.null(Tb.true)) {
                (determ.new = as.matrix(model.matrix(~month_new + 
                  trend_new - 1)))
            } else {
                (determ.new = as.matrix(model.matrix(~month_new + 
                  AO.new + trend_new - 1)))
            }
        })
        (pred.sarma = predict(model1b, n.ahead = n.ahead, newxreg = determ.new[, 
            -1])$pred)
    }
    if (!is.null(dcomp) && is.null(sdcomp)) {
        switch(dcomp, const = {
            if (is.null(Tb.true)) {
                determ = NULL
            } else {
                (determ = as.matrix(AO))
            }
        }, both = {
            (trend = seq(length(y.trunc)))
            if (is.null(Tb.true)) {
                (determ = as.matrix(trend))
            } else {
                (determ = cbind(AO, trend))
            }
        })
        (model1b = arimax(y.trunc, xreg = determ, include.mean = T, 
            ...))
        (AO.new = ((length(y.trunc) + 1):(length(y.trunc) + n.ahead) == 
            Tb.true))
        switch(dcomp, const = {
            if (is.null(Tb.true)) {
                determ.new = NULL
            } else {
                (determ.new = as.matrix(AO.new))
            }
        }, both = {
            (trend_new = ts(rep(1, n.ahead), start = c(end(y.trunc)[1], 
                end(y.trunc)[2] + 1), fre = frequency(y.trunc)))
            if (is.null(Tb.true)) {
                (determ.new = as.matrix(trend_new))
            } else {
                (determ.new = cbind(AO.new, trend_new))
            }
        })
        (pred.sarma = predict(model1b, n.ahead = n.ahead, newxreg = determ.new)$pred)
    }
    if (!is.null(dcomp) && !is.null(sdcomp)) {
        stop("Only sdcomp or dcomp allowed not both!\n")
    }
    return(list(model1b = model1b, pred.sarma = pred.sarma))
}
