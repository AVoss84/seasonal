function (b, p, S) 
{
    (b = matrix(b, byrow = T, ncol = S, nrow = p, dimnames = list(paste("p", 
        1:p, sep = ""), paste("s", 1:S, sep = ""))))
    Omega0 = diag(S)
    dimnames(Omega0) = list(paste("i", 1:S, sep = ""), paste("j", 
        1:S, sep = ""))
    Omega1 = matrix(0, ncol = S, nrow = S)
    dimnames(Omega1) = list(paste("i", 1:S, sep = ""), paste("j", 
        1:S, sep = ""))
    for (ii in 1:S) {
        for (jj in 1:S) {
            if (ii - jj <= p && ii > jj) {
                Omega0[ii, jj] = -b[ii - jj, ii] * (ii > jj)
            }
            if ((ii - jj + S) <= p) {
                Omega1[ii, jj] = b[ii - jj + S, ii]
            }
        }
    }
    return(list(betas = b, Omega0 = Omega0, Omega1 = Omega1, 
        det01 = det(Omega0 - Omega1 * 1)))
}
