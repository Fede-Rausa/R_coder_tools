
vet2onehot = function(vet, levels){
  ## converts vec to a matrix, after levels are specified
  ## it returns a matrix where each column is a dummy
    L = length(levels)
    N = length(vet)
    mat = matrix(NA, ncol=L, nrow=N)

    for (j in 1:L){
        lev = levels[j]
        col = as.numeric(vet==lev)
        mat[,j] = col
    }

    colnames(mat) = levels
    return(mat)
}


onehot2vet = function(mat, levels){
  ## it returns a vector from a one hot matrix
    L = length(levels)
    N = nrow(mat)   
    vet = numeric(N)

    for (j in 1:L){
        lev = levels[j]
        col = mat[,j]
        vet = ifelse(col==1, lev, vet)
    }

    return(vet)
}
