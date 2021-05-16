function (mat, y, coding = c("mean", "effect"), ref = 4) 
{
    S = 4
    (fullyears = trunc(length(cycle(y))/S))
    code = coding[1]
    (names = paste("s", 1:S, sep = ""))
    designS = array(0, dim = c(nrow(mat), S))
    for (ii in 1:nrow(mat)) {
        designS[ii, cycle(mat)[ii]] = 1
    }
    if (code == "effect") {
        designS[seq(ref, fullyears * S, by = S), ] = -1
        designS = designS[, -ref]
        cat(" Last season 4 was dropped!\n")
        names = names[-ref]
    }
    (designS = ts(designS, fre = S, start = start(mat), names = names))
    return(designS)
}
