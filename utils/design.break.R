function (mat, y, Tb, lags) 
{
    S = 4
    (fullyears = trunc(length(cycle(y))/S))
    names = rep(paste("s", 1:S, sep = ""), times = 2)
    designS = array(0, dim = c(nrow(mat), S))
    for (ii in 1:nrow(mat)) {
        designS[ii, cycle(mat)[ii]] = 1
    }
    design.brk1 = design.brk2 = array(0, dim = dim(designS))
    design.brk1[1:(Tb - lags), ] = designS[1:(Tb - lags), ]
    design.brk2[(Tb + 1 - lags):nrow(design.brk2), ] = designS[(Tb - 
        lags + 1):nrow(design.brk2), ]
    (design.br = cbind(design.brk1, design.brk2))
    (design.br = ts(design.br, fre = S, start = start(mat), names = names))
    return(design.br)
}
