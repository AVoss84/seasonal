function (y.trunc, Tb.true = NULL, sdcomp, n.ahead) 
{
    if (!is.null(Tb.true) && Tb.true >= length(y.trunc)) {
        stop("Break point outside the truncated sample. Set Tb.true <- NULL\n")
    }
    switch(sdcomp, const = {
        (trend = (start(y.trunc)[2]):(length(y.trunc) + ((start(y.trunc)[2]) - 
            1)))
        (st = 1 + (trend - 1)%%frequency(y.trunc))
        if (is.null(Tb.true)) {
            (ols = lm(y.trunc ~ factor(st) - 1))
        } else {
            (AO = (seq(y.trunc) == Tb.true))
            (ols = lm(y.trunc ~ factor(st) + AO - 1))
        }
        (trend.new = (trend[length(trend)] + 1):(trend[length(trend)] + 
            n.ahead))
        (st.new = 1 + (trend.new - 1)%%frequency(y.trunc))
        (new.dat = data.frame(Seas = factor(st.new)))
        (means = coef(ols)[1:12])
        (drift = means[st.new])
        (fore.dum = ts(drift, fre = frequency(y.trunc), start = c(1, 
            st.new[1])))
    }, both = {
        (trend = (start(y.trunc)[2]):(length(y.trunc) + ((start(y.trunc)[2]) - 
            1)))
        (st = 1 + (trend - 1)%%frequency(y.trunc))
        if (is.null(Tb.true)) {
            (ols = lm(y.trunc ~ trend + factor(st) - 1))
        } else {
            (AO = (seq(y.trunc) == Tb.true))
            (ols = lm(y.trunc ~ trend + factor(st) + AO - 1))
        }
        (trend.new = (trend[length(trend)] + 1):(trend[length(trend)] + 
            n.ahead))
        (st.new = 1 + (trend.new - 1)%%frequency(y.trunc))
        (new.dat = data.frame(Seas = factor(st.new), Time = trend.new))
        means = coef(ols)[-1]
        drift = means[st.new]
        (fore.dum = ts(coef(ols)[1] * trend.new + drift, fre = frequency(y.trunc), 
            start = c(1, st.new[1])))
    })
    return(list(OLS = ols, fore.dum = fore.dum))
}
